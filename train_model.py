#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
AutoGluon 建模脚本：对四组特征集进行交叉验证，输出预测概率和性能指标
用法：python train_model.py --base_dir /data_d/ZJJin/NSCLC
"""

import os
import sys
import argparse
import warnings
import pandas as pd
import numpy as np
import random
import time
import json
import pickle
import shutil
from collections import defaultdict
from itertools import combinations
from autogluon.tabular import TabularPredictor
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, accuracy_score, recall_score, precision_score, f1_score, confusion_matrix

warnings.filterwarnings("ignore")
BASE_SEED = 42

def create_output_dir(output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"创建目录: {output_dir}")
    return output_dir

def load_data(file_path):
    df = pd.read_csv(file_path)
    required = ['samples', 'metastasis']
    for col in required:
        if col not in df.columns:
            raise ValueError(f"文件缺少必要列: {col}")
    print(f"数据加载: {df.shape[0]} 样本, {df.shape[1]-2} 特征, 转移比例={df['metastasis'].mean():.2%}")
    df = df.set_index('samples')
    return df

def calculate_metrics(y_true, y_pred, y_pred_proba):
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    return {
        'auc': roc_auc_score(y_true, y_pred_proba),
        'accuracy': accuracy_score(y_true, y_pred),
        'sensitivity': recall_score(y_true, y_pred),
        'specificity': tn / (tn + fp) if (tn+fp)>0 else 0,
        'precision': precision_score(y_true, y_pred),
        'f1': f1_score(y_true, y_pred)
    }

def generate_combinations(n_folds=5):
    folds = list(range(n_folds))
    combos = []
    for train_idx in combinations(folds, 3):
        remaining = [f for f in folds if f not in train_idx]
        for val_idx in remaining:
            test_idx = [f for f in remaining if f != val_idx][0]
            combos.append({'train': list(train_idx), 'val': val_idx, 'test': test_idx})
    return combos

def cross_validation(df, output_dir, n_folds=5, repeats=1, feature_set_name="clinical"):
    label = 'metastasis'
    all_results = []
    model_performances = defaultdict(list)
    all_predictions = []
    all_feature_importances = []
    combos = generate_combinations(n_folds)
    print(f"\n交叉验证: {repeats}次重复 × {len(combos)}种组合, 特征集: {feature_set_name}")

    # 创建保存CSV预测的目录
    pred_csv_dir = os.path.join(output_dir, "test_predictions_csv")
    os.makedirs(pred_csv_dir, exist_ok=True)
    fi_dir = os.path.join(output_dir, "feature_importance")
    os.makedirs(fi_dir, exist_ok=True)

    for repeat in range(repeats):
        repeat_seed = BASE_SEED + repeat
        skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=repeat_seed)
        fold_dfs = [df.iloc[test_idx] for _, test_idx in skf.split(df, df[label])]

        for combo_idx, combo in enumerate(combos):
            train_data = pd.concat([fold_dfs[f] for f in combo['train']])
            val_data = fold_dfs[combo['val']]
            test_data = fold_dfs[combo['test']]
            print(f"  重复{repeat+1} 组合{combo_idx+1}: 训练{len(train_data)} 验证{len(val_data)} 测试{len(test_data)}")

            combo_dir = os.path.join(output_dir, f"repeat_{repeat+1}", f"combo_{combo_idx+1}")
            os.makedirs(combo_dir, exist_ok=True)

            predictor = TabularPredictor(label=label, path=combo_dir, eval_metric='roc_auc', verbosity=2)
            start = time.time()
            predictor.fit(train_data=train_data, tuning_data=val_data, presets='medium_quality', ag_args_fit={'num_cpus': 4})
            train_time = time.time() - start

            leaderboard = predictor.leaderboard(test_data, silent=True)
            for model_name in leaderboard['model']:
                model = predictor._trainer.load_model(model_name)
                # 预测
                y_train_pred = model.predict(train_data)
                y_train_proba = model.predict_proba(train_data).iloc[:,1] if isinstance(model.predict_proba(train_data), pd.DataFrame) else model.predict_proba(train_data)[:,1]
                y_val_pred = model.predict(val_data)
                y_val_proba = model.predict_proba(val_data).iloc[:,1] if isinstance(model.predict_proba(val_data), pd.DataFrame) else model.predict_proba(val_data)[:,1]
                y_test_pred = model.predict(test_data)
                y_test_proba = model.predict_proba(test_data).iloc[:,1] if isinstance(model.predict_proba(test_data), pd.DataFrame) else model.predict_proba(test_data)[:,1]

                train_metrics = calculate_metrics(train_data[label], y_train_pred, y_train_proba)
                val_metrics = calculate_metrics(val_data[label], y_val_pred, y_val_proba)
                test_metrics = calculate_metrics(test_data[label], y_test_pred, y_test_proba)

                result = {
                    'repeat': repeat+1, 'combo': combo_idx+1, 'seed': repeat_seed,
                    'model_name': model_name, 'feature_set': feature_set_name,
                    'train_auc': train_metrics['auc'], 'train_accuracy': train_metrics['accuracy'],
                    'train_sensitivity': train_metrics['sensitivity'], 'train_specificity': train_metrics['specificity'],
                    'val_auc': val_metrics['auc'], 'val_accuracy': val_metrics['accuracy'],
                    'val_sensitivity': val_metrics['sensitivity'], 'val_specificity': val_metrics['specificity'],
                    'test_auc': test_metrics['auc'], 'test_accuracy': test_metrics['accuracy'],
                    'test_sensitivity': test_metrics['sensitivity'], 'test_specificity': test_metrics['specificity'],
                    'train_time': train_time
                }
                all_results.append(result)
                model_performances[model_name].append(result)

                # 保存预测CSV
                pred_df = pd.DataFrame({
                    'sample_id': test_data.index.tolist(),
                    'true_label': test_data[label].values,
                    'predicted_prob': y_test_proba,
                    'model_name': model_name,
                    'feature_set': feature_set_name,
                    'repeat': repeat+1,
                    'combo': combo_idx+1
                })
                csv_file = os.path.join(pred_csv_dir, f"{feature_set_name}_repeat{repeat+1}_combo{combo_idx+1}_{model_name}.csv")
                pred_df.to_csv(csv_file, index=False)
                all_predictions.append(pred_df)

                # 特征重要性
                try:
                    importance = predictor.feature_importance(data=test_data, model=model_name, features=test_data.columns.drop(label), num_shuffle_sets=10)
                    importance = importance.reset_index().rename(columns={'index': 'feature'})
                    importance['repeat'] = repeat+1
                    importance['combo'] = combo_idx+1
                    importance['model_name'] = model_name
                    importance['feature_set'] = feature_set_name
                    fi_file = os.path.join(fi_dir, f"{feature_set_name}_repeat{repeat+1}_combo{combo_idx+1}_{model_name}_fi.csv")
                    importance.to_csv(fi_file, index=False)
                    all_feature_importances.append(importance)
                except Exception as e:
                    print(f"特征重要性计算失败: {e}")

            # 保存当前组合结果
            with open(os.path.join(combo_dir, 'combo_results.json'), 'w') as f:
                json.dump([r for r in all_results if r['repeat']==repeat+1 and r['combo']==combo_idx+1], f, indent=2)

        # 保存当前重复结果
        pd.DataFrame([r for r in all_results if r['repeat']==repeat+1]).to_csv(os.path.join(output_dir, f"repeat_{repeat+1}", 'repeat_results.csv'), index=False)

    # 汇总所有结果
    all_results_df = pd.DataFrame(all_results)
    all_results_df.to_csv(os.path.join(output_dir, 'all_cv_results.csv'), index=False)
    print(f"所有结果保存至: {os.path.join(output_dir, 'all_cv_results.csv')}")

    # 汇总预测CSV
    pd.concat(all_predictions, ignore_index=True).to_csv(os.path.join(output_dir, 'all_test_predictions.csv'), index=False)
    # 汇总特征重要性
    if all_feature_importances:
        pd.concat(all_feature_importances, ignore_index=True).to_csv(os.path.join(output_dir, 'all_feature_importances_raw.csv'), index=False)

    # 模型性能汇总表
    summary_list = []
    for model_name, res_list in model_performances.items():
        dfm = pd.DataFrame(res_list)
        summary_list.append({
            'model_name': model_name,
            'num_runs': len(res_list),
            'feature_set': feature_set_name,
            'train_auc': f"{dfm['train_auc'].mean():.4f}±{dfm['train_auc'].std():.4f}",
            'val_auc': f"{dfm['val_auc'].mean():.4f}±{dfm['val_auc'].std():.4f}",
            'test_auc': f"{dfm['test_auc'].mean():.4f}±{dfm['test_auc'].std():.4f}",
            'test_accuracy': f"{dfm['test_accuracy'].mean():.4f}±{dfm['test_accuracy'].std():.4f}",
            'test_sensitivity': f"{dfm['test_sensitivity'].mean():.4f}±{dfm['test_sensitivity'].std():.4f}",
            'test_specificity': f"{dfm['test_specificity'].mean():.4f}±{dfm['test_specificity'].std():.4f}"
        })
    summary_df = pd.DataFrame(summary_list)
    summary_df.to_csv(os.path.join(output_dir, 'model_performance_summary.csv'), index=False)
    return all_results_df, summary_df

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--base_dir', type=str, default='/data_d/ZJJin/NSCLC', help='项目根目录')
    args = parser.parse_args()

    base_dir = args.base_dir
    feature_dir = os.path.join(base_dir, 'Figure5', 'features')
    output_root = os.path.join(base_dir, 'Figure5', 'cv_results')
    os.makedirs(output_root, exist_ok=True)

    feature_sets = {
        'clinical': os.path.join(feature_dir, 'Group1_Clinical.csv'),
        'clinical_rna': os.path.join(feature_dir, 'Group2_Clinical_Edit.csv'),
        'clinical_expression': os.path.join(feature_dir, 'Group3_Clinical_Exp.csv'),
        'clinical_rna_expression': os.path.join(feature_dir, 'Group4_All_Features.csv')
    }

    N_FOLDS = 5
    REPEATS = 1

    for name, path in feature_sets.items():
        print(f"\n{'='*50}\n运行特征集: {name}\n{'='*50}")
        output_dir = os.path.join(output_root, name)
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)
        os.makedirs(output_dir, exist_ok=True)
        df = load_data(path)
        cross_validation(df, output_dir, n_folds=N_FOLDS, repeats=REPEATS, feature_set_name=name)

    print("\n所有特征集建模完成！")

if __name__ == '__main__':
    main()