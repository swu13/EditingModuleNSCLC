#!/usr/bin/env Rscript
# Figure5 绘图脚本：雷达图、脊线图、DeLong 检验
# 用法：Rscript figure5_plots.R [项目根目录]

args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) > 0) args[1] else "."

# 定义输入输出路径
cv_root <- file.path(base_dir, "Figure5", "cv_results")
plot_dir <- file.path(base_dir, "Figure5", "plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# 加载包
library(tidyverse)
library(ggradar)
library(scales)
library(ggrepel)
library(ggridges)
library(ggsci)
library(pROC)

# ========================== (b) 雷达图 ==========================
cat("绘制雷达图...\n")
# 读取四个特征集的性能汇总
feature_sets <- c("clinical", "clinical_rna", "clinical_expression", "clinical_rna_expression")
summary_list <- list()
for (fs in feature_sets) {
  file_path <- file.path(cv_root, fs, "model_performance_summary.csv")
  if (file.exists(file_path)) {
    df <- read.csv(file_path, stringsAsFactors = FALSE) %>% mutate(feature_set = fs)
    summary_list[[fs]] <- df
  } else {
    warning("未找到文件: ", file_path)
  }
}
summary_data <- bind_rows(summary_list)

# 提取测试 AUC 均值
summary_data <- summary_data %>%
  mutate(
    auc_mean = as.numeric(sub("±.*", "", test_auc)),
    auc_sd   = as.numeric(str_extract(test_auc, "(?<=±)\\d+\\.\\d+"))
  )

# 目标模型顺序
target_models <- c("ExtraTreesEntr", "CatBoost", "ExtraTreesGini", "LightGBM", 
                   "RandomForestGini", "RandomForestEntr", "LightGBMLarge", 
                   "LightGBMXT", "NeuralNetFastAI")
model_order <- c("LightGBMXT", "CatBoost", "NeuralNetFastAI", "LightGBM", 
                 "LightGBMLarge", "ExtraTreesEntr", "RandomForestGini", 
                 "ExtraTreesGini", "RandomForestEntr")

filtered_data <- summary_data %>% filter(model_name %in% target_models)

# 准备雷达图数据
radar_data <- filtered_data %>%
  filter(feature_set %in% feature_sets) %>%
  mutate(feature_set = case_when(
    feature_set == "clinical" ~ "Clinical",
    feature_set == "clinical_rna" ~ "Clinical+RNA",
    feature_set == "clinical_expression" ~ "Clinical+Expression",
    feature_set == "clinical_rna_expression" ~ "Clinical+Both"
  )) %>%
  select(model_name, feature_set, auc_mean) %>%
  mutate(auc_mapped = 0.1 + (auc_mean - 0.5) * (0.8 / 0.2)) %>%
  pivot_wider(names_from = model_name, values_from = auc_mapped) %>%
  select(feature_set, all_of(model_order))

# 自定义雷达图函数
custom_radar <- function(data, colors, title) {
  p <- ggradar(
    plot.data = data,
    font.radar = "sans",
    grid.label.size = 5,
    axis.label.size = 5,
    axis.labels = model_order,
    group.line.width = 1.8,
    group.point.size = 5,
    legend.text.size = 13,
    legend.position = "bottom",
    grid.line.width = 0.7,
    grid.min = 0.1,
    grid.mid = 0.5,
    grid.max = 0.9,
    values.radar = c("0.5", "0.6", "0.7"),
    background.circle.colour = "white",
    gridline.mid.colour = "grey60",
    gridline.max.colour = "grey60",
    gridline.min.colour = "grey60",
    plot.title = title,
    fill = TRUE,
    fill.alpha = 0.1,
    gridline.label.offset = 0.15
  ) +
    scale_color_manual(values = colors) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 15)),
      legend.title = element_blank(),
      legend.margin = margin(t = 15, b = 5),
      legend.key = element_rect(fill = NA),
      legend.key.size = unit(1.5, "lines"),
      plot.margin = margin(1, 1, 1.5, 1, "cm")
    )
  return(p)
}

radar_colors <- c("Clinical" = "#3B4992", "Clinical+RNA" = "#EE0000",
                  "Clinical+Expression" = "#008B45", "Clinical+Both" = "#631879")
p_radar <- custom_radar(radar_data, radar_colors, 
                        "Comparison of Feature Sets for Metastasis Prediction")
ggsave(file.path(plot_dir, "Figure5_radar.png"), p_radar, width = 12, height = 10, dpi = 300)

# ========================== (c1) 脊线图：RNA编辑增益 ==========================
cat("绘制脊线图...\n")
# 读取每个特征集的 all_cv_results.csv
cv_results <- list()
for (fs in feature_sets) {
  file_path <- file.path(cv_root, fs, "all_cv_results.csv")
  if (file.exists(file_path)) {
    df <- read.csv(file_path, stringsAsFactors = FALSE) %>% mutate(feature_set = fs)
    cv_results[[fs]] <- df
  }
}
perf <- bind_rows(cv_results) %>%
  filter(model_name == "NeuralNetFastAI") %>%
  select(repeat., combo, seed, feature_set, test_auc)

# 计算增益
gain_data <- perf %>%
  pivot_wider(names_from = feature_set, values_from = test_auc) %>%
  mutate(
    RNA_Gain_Clinical = clinical_rna - clinical,
    RNA_Gain_Expression = clinical_rna_expression - clinical_expression
  ) %>%
  select(repeat., combo, seed, RNA_Gain_Clinical, RNA_Gain_Expression) %>%
  pivot_longer(cols = c(RNA_Gain_Clinical, RNA_Gain_Expression),
               names_to = "comparison", values_to = "gain") %>%
  mutate(comparison = recode(comparison,
                             RNA_Gain_Clinical = "Clinical vs Clinical+RNA",
                             RNA_Gain_Expression = "Clinical+Exp vs Clinical+Both"))

avg_gain <- gain_data %>% group_by(comparison) %>% summarise(avg_gain = mean(gain))

p_ridge <- ggplot(gain_data, aes(x = gain, y = comparison, fill = comparison)) +
  geom_density_ridges(alpha = 0.85, scale = 0.95, bandwidth = 0.004,
                      color = "white", size = 0.6,
                      jittered_points = TRUE, point_shape = 21, point_size = 2,
                      position = position_points_jitter(width = 0.05, height = 0)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#999999", linewidth = 1) +
  geom_segment(data = avg_gain, aes(x = avg_gain, xend = avg_gain,
                                    y = as.numeric(comparison) - 0.2,
                                    yend = as.numeric(comparison) + 0.2),
               color = "black", linewidth = 1.2) +
  geom_label(data = avg_gain, aes(x = avg_gain, y = as.numeric(comparison) + 0.3,
                                  label = sprintf("Mean: %+.3f", avg_gain)),
             size = 4.5, fontface = "bold", fill = "white") +
  scale_fill_manual(values = c("Clinical vs Clinical+RNA" = "#5B9BD5",
                               "Clinical+Exp vs Clinical+Both" = "#ED7D31")) +
  labs(title = "Performance Impact of RNA Editing Features",
       x = "AUC Improvement", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 12)) +
  scale_x_continuous(limits = c(-0.06, 0.12), breaks = seq(-0.05, 0.10, by = 0.05))
ggsave(file.path(plot_dir, "Figure5_ridgeplot.pdf"), p_ridge, width = 10, height = 6)

# ========================== (c2) DeLong 检验 ==========================
cat("执行 DeLong 检验...\n")
# 读取四组特征集的测试集预测结果
all_preds <- list()
for (fs in feature_sets) {
  pred_file <- file.path(cv_root, fs, "all_test_predictions.csv")
  if (file.exists(pred_file)) {
    df <- read.csv(pred_file, stringsAsFactors = FALSE) %>% mutate(feature_set = fs)
    all_preds[[fs]] <- df
  }
}
pred_all <- bind_rows(all_preds)

# 只关注 NeuralNetFastAI 模型（可根据需要修改）
model_name <- "NeuralNetFastAI"
sub <- pred_all %>% filter(model_name == !!model_name)
if (nrow(sub) == 0) {
  warning("未找到模型 ", model_name, " 的预测数据，跳过 DeLong 检验。")
} else {
  roc_list <- list()
  for (fs in feature_sets) {
    tmp <- sub %>% filter(feature_set == fs)
    if (nrow(tmp) > 0) {
      roc_obj <- roc(tmp$true_label, tmp$predicted_prob, quiet = TRUE, direction = "<")
      roc_list[[fs]] <- roc_obj
    }
  }
  
  # 定义比较
  comparisons <- list(
    list(name = "Clinical vs Clinical+RNA", fs1 = "clinical", fs2 = "clinical_rna"),
    list(name = "Clinical+Exp vs Clinical+Both", fs1 = "clinical_expression", fs2 = "clinical_rna_expression")
  )
  
  delong_res <- data.frame()
  for (comp in comparisons) {
    roc1 <- roc_list[[comp$fs1]]
    roc2 <- roc_list[[comp$fs2]]
    if (!is.null(roc1) && !is.null(roc2)) {
      test <- roc.test(roc1, roc2, method = "delong")
      diff_auc <- roc2$auc - roc1$auc
      delong_res <- rbind(delong_res, data.frame(
        comparison = comp$name,
        auc_diff = diff_auc,
        p_value = test$p.value
      ))
    }
  }
  write.csv(delong_res, file.path(plot_dir, "DeLong_test_results.csv"), row.names = FALSE)
  print(delong_res)
  
  # 可选：将 p 值标注到脊线图上（已单独保存，也可在图中绘制文字）
  # 此处仅保存结果表格，如需在图中添加，可在此基础上修改 ridge_plot。
}

cat("Figure5 绘图完成！图片保存在:", plot_dir, "\n")