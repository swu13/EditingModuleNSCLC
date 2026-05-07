#!/usr/bin/env Rscript
# Figure1 分析脚本：评估RNA编辑事件对EIF2AK2表达的预测能力
# 用法：Rscript figure1.R [项目根目录]
# 如果不指定项目根目录，默认使用当前目录

# 解析命令行参数
args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) > 0) args[1] else "."

# 构造输入文件路径
exp_file      <- file.path(base_dir, "Expression", "Exp_tpm_01A_data.csv")
editing_file  <- file.path(base_dir, "Editing", "Filled_lung_0.2.txt")
annotation_file <- file.path(base_dir, "Results", "editing_events_unique_gene_annotation.txt")

# 输出目录
output_dir <- file.path(base_dir, "Figure1_output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 加载必要的包
library(dplyr)
library(tidyr)
library(ggplot2)
library(caret)
library(ggpubr)
library(gridExtra)
library(doParallel)
library(grid)
library(forcats)

# ------------------- 读取数据 -------------------
cat("读取数据...\n")
result <- read.table(annotation_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 筛选EIF2AK2相关的编辑事件
eif2ak2_events <- result$Editing[grepl("EIF2AK2", result$gene_name)]

rna_editing <- read.table(editing_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(rna_editing) <- rna_editing$Editing
rna_editing <- rna_editing[, -1, drop = FALSE]
filtered_editing <- rna_editing[rownames(rna_editing) %in% eif2ak2_events, , drop = FALSE]

exp_tumor <- read.csv(exp_file, row.names = 1, check.names = FALSE)
if (!"EIF2AK2" %in% rownames(exp_tumor)) stop("EIF2AK2 不在表达矩阵中")
eif2ak2_exp <- data.frame(t(exp_tumor["EIF2AK2", , drop = FALSE]), check.names = FALSE)

editing_mat <- as.data.frame(t(filtered_editing))
merged <- merge(editing_mat, eif2ak2_exp, by = "row.names")
colnames(merged)[1] <- "sample"
write.csv(merged, file.path(output_dir, "EIF2AK2_tpm_with_editing_events.csv"), row.names = FALSE)

# 准备建模数据
merged_data <- read.csv(file.path(output_dir, "EIF2AK2_tpm_with_editing_events.csv"))
y <- merged_data$EIF2AK2
x <- merged_data %>% select(-sample, -EIF2AK2) %>% as.data.frame()

# ------------------- 交叉验证设置 -------------------
set.seed(123)
cv_ctrl <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 100,
  savePredictions = "final",
  returnResamp = "all",
  verboseIter = FALSE,
  allowParallel = TRUE
)

# ------------------- 单变量建模函数 -------------------
train_univariate_models <- function(x, y, cv_ctrl) {
  univariate_mse_matrix <- matrix(NA, nrow = 100, ncol = ncol(x))
  colnames(univariate_mse_matrix) <- colnames(x)
  univariate_r2_matrix <- matrix(NA, nrow = 100, ncol = ncol(x))
  colnames(univariate_r2_matrix) <- colnames(x)
  univariate_models <- list()
  
  for (i in seq_along(colnames(x))) {
    feature <- colnames(x)[i]
    current_data <- data.frame(edit = x[[feature]], EIF2AK2 = y)
    uni_model <- tryCatch({
      train(
        EIF2AK2 ~ edit,
        data = current_data,
        method = "lm",
        trControl = cv_ctrl,
        preProcess = c("center", "scale")
      )
    }, error = function(e) {
      message(sprintf("特征 %s 建模失败: %s", feature, e$message))
      return(NULL)
    })
    if (!is.null(uni_model)) {
      mse_per_repeat <- uni_model$resample %>%
        mutate(Repeat = as.numeric(gsub(".*Rep", "", Resample))) %>%
        group_by(Repeat) %>%
        summarise(MSE = mean(RMSE^2)) %>%
        pull(MSE)
      r2_per_repeat <- uni_model$resample %>%
        mutate(Repeat = as.numeric(gsub(".*Rep", "", Resample))) %>%
        group_by(Repeat) %>%
        summarise(R2 = mean(Rsquared)) %>%
        pull(R2)
      
      univariate_mse_matrix[, i] <- mse_per_repeat
      univariate_r2_matrix[, i] <- r2_per_repeat
      univariate_models[[feature]] <- uni_model
    }
  }
  return(list(
    univariate_mse_matrix = univariate_mse_matrix,
    univariate_r2_matrix = univariate_r2_matrix,
    univariate_models = univariate_models
  ))
}

cat("执行单变量交叉验证...\n")
univariate_model_results <- train_univariate_models(x, y, cv_ctrl)
univariate_mse_matrix <- univariate_model_results$univariate_mse_matrix
univariate_r2_matrix <- univariate_model_results$univariate_r2_matrix
univariate_models <- univariate_model_results$univariate_models

# 筛选有效特征并找出最佳单变量
valid_columns <- colSums(!is.na(univariate_mse_matrix)) > 0
avg_mse <- colMeans(univariate_mse_matrix[, valid_columns], na.rm = TRUE)
best_uni_feature <- names(which.min(avg_mse))
best_uni_mse <- univariate_mse_matrix[, best_uni_feature]

# ------------------- 多变量建模 -------------------
cat("执行多变量交叉验证...\n")
multivariate_model <- train(
  x = x[, valid_columns],
  y = y,
  method = "lm",
  trControl = cv_ctrl,
  preProcess = c("center", "scale")
)

process_multivariate_cv_results <- function(model) {
  multi_mse <- model$resample %>%
    mutate(Repeat = as.numeric(gsub(".*Rep", "", Resample))) %>%
    group_by(Repeat) %>%
    summarise(MSE = mean(RMSE^2)) %>%
    pull(MSE)
  multi_r2 <- model$resample %>%
    mutate(Repeat = as.numeric(gsub(".*Rep", "", Resample))) %>%
    group_by(Repeat) %>%
    summarise(R2 = mean(Rsquared)) %>%
    pull(R2)
  return(list(multi_mse = multi_mse, multi_r2 = multi_r2))
}

multivariate_cv_results <- process_multivariate_cv_results(multivariate_model)
multi_mse <- multivariate_cv_results$multi_mse
multi_r2 <- multivariate_cv_results$multi_r2

# ------------------- 绘图辅助函数 -------------------
morandi_colors <- c("#7CAEF0", "#F58787", "#8EC9A7", "#D4A5C2", "#FFD700", "#A0A0A0")

plot_uni_model <- function(feature, rank, mse_values, r2_values, univariate_models, x, y) {
  model <- univariate_models[[feature]]
  df <- data.frame(
    Editing = x[[feature]],
    Expression = y,
    Predicted = predict(model$finalModel)
  )
  best_index <- which.min(mse_values)
  best_r_squared <- r2_values[best_index]
  min_mse <- min(mse_values, na.rm = TRUE)
  mse_sd <- sd(mse_values, na.rm = TRUE)
  
  ggplot(df, aes(x = Editing, y = Expression)) +
    geom_point(color = ifelse(rank == 1, morandi_colors[1], morandi_colors[3]),
               alpha = 0.7, size = 3, shape = 19) +
    geom_smooth(method = "lm", color = morandi_colors[2], se = TRUE,
                fill = "#F0F0F0", linewidth = 1.2) +
    annotate("text",
             x = min(df$Editing) + 0.05 * diff(range(df$Editing)),
             y = max(df$Expression) - 0.1 * diff(range(df$Expression)),
             label = sprintf("R² = %.3f\nMSE = %.4f ± %.4f",
                            best_r_squared, min_mse, mse_sd),
             hjust = 0, vjust = 1, size = 4.5, color = "#333333", lineheight = 0.8) +
    labs(title = paste("Top", rank, ":", feature),
         x = "RNA Editing Level", y = "EIF2AK2 Expression") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
          panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "grey30"),
          panel.border = element_rect(color = "#EAEAEA", fill = NA),
          plot.margin = margin(10, 15, 10, 15))
}

# ------------------- 生成图P1 -------------------
cat("生成图P1（最佳单变量模型）...\n")
top2_features <- names(sort(avg_mse))[1:2]
p1 <- grid.arrange(
  plot_uni_model(top2_features[1], 1,
                 univariate_mse_matrix[, top2_features[1]],
                 univariate_r2_matrix[, top2_features[1]],
                 univariate_models, x, y),
  plot_uni_model(top2_features[2], 2,
                 univariate_mse_matrix[, top2_features[2]],
                 univariate_r2_matrix[, top2_features[2]],
                 univariate_models, x, y),
  ncol = 2,
  top = textGrob("Best Univariate Models (Showing Minimum MSE from 100 CV Repeats)",
                 gp = gpar(fontsize = 14, fontface = "bold", col = "#333333"))
)
ggsave(file.path(output_dir, "Figure1_P1_univariate_models.png"), p1, width = 12, height = 6, dpi = 300)

# ------------------- 生成图P3（系数图）-------------------
cat("生成图P3（多变量系数图）...\n")
coef_list <- list()
for (i in 1:length(multivariate_model$control$index)) {
  train_idx <- multivariate_model$control$index[[i]]
  lm_model <- lm(y ~ ., data = data.frame(x[train_idx, valid_columns], y = y[train_idx]))
  coef_values <- coef(lm_model)[-1]
  ordered_coefs <- coef_values[match(colnames(x[, valid_columns]), names(coef_values))]
  coef_list[[i]] <- data.frame(
    Repeat = i,
    Feature = colnames(x[, valid_columns]),
    Coefficient = ordered_coefs
  )
}
coef_df <- bind_rows(coef_list) %>%
  group_by(Feature) %>%
  mutate(Mean_Coeff = mean(Coefficient, na.rm = TRUE),
         SD_Coeff = sd(Coefficient, na.rm = TRUE))

top10_features <- coef_df %>%
  group_by(Feature) %>%
  summarise(Abs_Mean = abs(mean(Coefficient))) %>%
  arrange(desc(Abs_Mean)) %>%
  head(10) %>%
  pull(Feature)

plot_data <- coef_df %>%
  filter(Feature %in% top10_features) %>%
  mutate(Feature = fct_reorder(Feature, abs(Mean_Coeff), .fun = mean, .desc = TRUE),
         Sign = factor(ifelse(Mean_Coeff > 0, "Positive", "Negative")))

p3 <- ggplot(plot_data, aes(x = reorder(Feature, Mean_Coeff), y = Coefficient)) +
  geom_boxplot(aes(fill = Sign), width = 0.6,
               outlier.shape = 21, outlier.size = 1.5, outlier.color = "grey30") +
  geom_hline(yintercept = 0, color = "#333333", linewidth = 0.8, linetype = "dashed") +
  geom_point(aes(y = Mean_Coeff), shape = 18, size = 4, color = "black") +
  coord_flip() +
  scale_fill_manual(values = c("Positive" = "#4E79A7", "Negative" = "#E15759"),
                    labels = c("Positive" = "Positive effect", "Negative" = "Negative effect")) +
  labs(x = NULL, y = "Regression Coefficient", fill = "Effect Direction") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.x = element_line(color = "grey90"),
        panel.grid.major.y = element_blank(),
        axis.text.y = element_text(color = "black", size = 11),
        legend.position = "bottom",
        legend.title = element_text(face = "bold")) +
  scale_y_continuous(breaks = seq(floor(min(plot_data$Coefficient)*2)/2,
                                  ceiling(max(plot_data$Coefficient)*2)/2, by = 0.5),
                     expand = expansion(mult = c(0.05, 0.05))) +
  annotate("text", x = Inf, y = Inf, label = "Black diamonds show mean coefficients",
           hjust = 1.1, vjust = 1.5, size = 3.5, color = "grey40")

ggsave(file.path(output_dir, "Figure1_P3_coefficients.png"), p3, width = 10, height = 6, dpi = 300)

# ------------------- 生成图P4（MSE比较）-------------------
cat("生成图P4（MSE比较）...\n")
performance_df <- data.frame(
  Model = factor(rep(c("Best Univariate", "Multivariate"), each = 100),
                 levels = c("Best Univariate", "Multivariate")),
  MSE = c(best_uni_mse, multi_mse)
)

p4 <- ggplot(performance_df, aes(x = Model, y = MSE, fill = Model)) +
  geom_boxplot(width = 0.6, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.4, size = 1.8, color = "grey40") +
  scale_fill_manual(values = morandi_colors) +
  stat_compare_means(method = "t.test",
                     aes(label = after_stat(sprintf("P = %.2e", p))),
                     label.x = 1.5, label.y = max(performance_df$MSE) * 1.1,
                     size = 4.5, color = "grey30", bracket.size = 0.6, tip.length = 0.02) +
  labs(title = "100 Repeats of 10-Fold Cross-Validation",
       x = NULL, y = "Mean Squared Error") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

ggsave(file.path(output_dir, "Figure1_P4_MSE_comparison.png"), p4, width = 8, height = 6, dpi = 300)

# ------------------- 生成图P5（预测vs实际）-------------------
cat("生成图P5（多变量模型预测性能）...\n")
cv_pred <- multivariate_model$pred
best_index_multi <- which.min(multi_mse)
best_r2_multi <- multi_r2[best_index_multi]
min_mse_multi <- min(multi_mse, na.rm = TRUE)
multi_mse_sd <- sd(multi_mse, na.rm = TRUE)

p5 <- ggplot(cv_pred, aes(x = obs, y = pred)) +
  geom_point(shape = 14, color = "#D8BFD8", size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", color = "#F58787", se = FALSE, linewidth = 1.2) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", linewidth = 1.2) +
  annotate("text",
           x = min(cv_pred$obs) + 0.05 * diff(range(cv_pred$obs)),
           y = max(cv_pred$pred) - 0.1 * diff(range(cv_pred$pred)),
           label = sprintf("R² = %.3f\nMSE = %.4f ± %.4f",
                           best_r2_multi, min_mse_multi, multi_mse_sd),
           hjust = 0, vjust = 1, size = 6, color = "#333333", lineheight = 0.8) +
  labs(x = "Actual EIF2AK2 Expression", y = "Predicted Expression") +
  theme_minimal(base_size = 20) +
  theme(plot.title = element_text(face = "bold", size = 26, hjust = 0.5, margin = margin(b = 15)),
        plot.subtitle = element_text(size = 18, hjust = 0.5, color = "#666666"),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank())

ggsave(file.path(output_dir, "Figure1_P5_multivariate_performance.png"), p5, width = 8, height = 7, dpi = 300)

cat("所有图片已保存至:", output_dir, "\n")