#!/usr/bin/env Rscript
# Figure4 分析脚本：RNA编辑事件间相关性、PC1提取、表达整合、中介分析及免疫亚型差异
# 用法：Rscript figure4.R [项目根目录]
# 如果不指定项目根目录，默认使用当前目录

args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) > 0) args[1] else "."

# ---------- 输入文件路径 ----------
editing_file      <- file.path(base_dir, "Editing", "Filled_lung_0.2.txt")
exp_file          <- file.path(base_dir, "Expression", "Exp_TPM_data.csv")
clinical_immune_file <- file.path(base_dir, "Expression", "mmc2.xlsx")  # 需要用户提供

# 输出目录
output_dir <- file.path(base_dir, "Figure4_output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 加载必要的包
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(corrplot)
library(mediation)
library(readxl)

cat("Figure4 分析开始\n")

# ------------------- 1. 读取RNA编辑数据并筛选目标事件 -------------------
cat("1. 读取RNA编辑数据...\n")
if (!file.exists(editing_file)) stop("找不到编辑事件文件: ", editing_file)
rna_editing_data <- read.table(editing_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

eif2ak2_edits <- c("chr2_37102092_-", "chr2_37100891_-", "chr2_37103468_-", "chr2_37103895_-", "chr2_37103922_-")
prkar2a_edits <- c("chr3_48750200_-", "chr3_48750207_-", "chr3_48750423_-")

selected_data <- rna_editing_data %>%
  filter(Editing %in% c(eif2ak2_edits, prkar2a_edits))

# 转置
transposed_data <- as.data.frame(t(selected_data[, -1]))
colnames(transposed_data) <- selected_data$Editing
transposed_data <- transposed_data %>%
  mutate(across(everything(), as.numeric))

# ------------------- 2. 计算编辑事件间的Pearson相关性（P1） -------------------
cat("2. 计算EIF2AK2与PRKAR2A编辑事件间的相关性...\n")
cor_results <- data.frame()
for (e in eif2ak2_edits) {
  for (p in prkar2a_edits) {
    if (e %in% colnames(transposed_data) & p %in% colnames(transposed_data)) {
      ct <- cor.test(transposed_data[[e]], transposed_data[[p]], method = "pearson")
      cor_results <- rbind(cor_results, data.frame(
        EIF2AK2_Editing = e,
        PRKAR2A_Editing = p,
        Correlation = ct$estimate,
        P_value = ct$p.value,
        stringsAsFactors = FALSE
      ))
    }
  }
}
write.csv(cor_results, file.path(output_dir, "editing_event_correlations.csv"), row.names = FALSE)
cat("编辑事件相关性结果已保存。\n")

# ------------------- 3. 读取表达数据，提取目标基因的肿瘤样本 -------------------
cat("3. 读取表达数据...\n")
if (!file.exists(exp_file)) stop("找不到表达文件: ", exp_file)
expression_data <- read.csv(exp_file, row.names = 1, check.names = FALSE)
keep_cols <- grep("^.*\\..*\\..*\\.01A$", colnames(expression_data))
tumor_exp <- expression_data[, keep_cols, drop = FALSE]

target_genes <- c("EIF2AK2", "PRKAR2A", "DDX3X")
selected_rows <- tumor_exp[target_genes, , drop = FALSE]
selected_rows_df <- data.frame(t(selected_rows))
cat("肿瘤样本数:", nrow(selected_rows_df), " 基因数:", ncol(selected_rows_df), "\n")

# ------------------- 4. 对EIF2AK2和PRKAR2A编辑事件做PCA提取PC1 -------------------
cat("4. 对编辑事件进行PCA，提取PC1...\n")
# EIF2AK2 editing PCA
eif2ak2_data <- rna_editing_data %>%
  filter(Editing %in% eif2ak2_edits) %>%
  tibble::column_to_rownames("Editing")
pca_eif <- prcomp(t(eif2ak2_data), scale. = FALSE)
EIF2AK2_Editing_PC1 <- data.frame(EIF2AK2_Editing_PC1 = pca_eif$x[, 1])

# PRKAR2A editing PCA
prkar2a_data <- rna_editing_data %>%
  filter(Editing %in% prkar2a_edits) %>%
  tibble::column_to_rownames("Editing")
pca_prk <- prcomp(t(prkar2a_data), scale. = FALSE)
PRKAR2A_Editing_PC1 <- data.frame(PRKAR2A_Editing_PC1 = pca_prk$x[, 1])

# 合并PC1与表达数据
pc_df <- merge(EIF2AK2_Editing_PC1, PRKAR2A_Editing_PC1, by = 0)
rownames(pc_df) <- pc_df$Row.names
pc_df <- pc_df[, -1, drop = FALSE]
merged_data <- merge(pc_df, selected_rows_df, by = 0)
colnames(merged_data)[1] <- "Sample"
write.csv(merged_data, file.path(output_dir, "merged_data_PC1_expression.csv"), row.names = FALSE)
cat("合并数据已保存。\n")

# ------------------- 5. 相关性分析（P2, P4, P5） -------------------
cat("5. 计算变量间的相关性...\n")
pairs <- list(
  c("EIF2AK2_Editing_PC1", "EIF2AK2"),
  c("PRKAR2A_Editing_PC1", "PRKAR2A"),
  c("EIF2AK2", "PRKAR2A"),
  c("PRKAR2A_Editing_PC1", "EIF2AK2_Editing_PC1"),
  c("DDX3X", "PRKAR2A_Editing_PC1"),
  c("DDX3X", "EIF2AK2_Editing_PC1")
)
cor_res <- lapply(pairs, function(pair) {
  ct <- cor.test(merged_data[[pair[1]]], merged_data[[pair[2]]])
  data.frame(Variable1 = pair[1], Variable2 = pair[2],
             Correlation = ct$estimate, Pvalue = ct$p.value)
}) %>% bind_rows()
write.csv(cor_res, file.path(output_dir, "variable_correlations.csv"), row.names = FALSE)
print(cor_res)

# ------------------- 6. 按DDX3X表达中位数分组，对PC1进行t检验 -------------------
merged_data <- merged_data %>%
  mutate(DDX3X_group = ifelse(DDX3X > median(DDX3X), "High", "Low") %>% factor())
t_res1 <- t.test(PRKAR2A_Editing_PC1 ~ DDX3X_group, data = merged_data)
t_res2 <- t.test(EIF2AK2_Editing_PC1 ~ DDX3X_group, data = merged_data)
t_result_df <- data.frame(
  Variable = c("PRKAR2A_Editing_PC1", "EIF2AK2_Editing_PC1"),
  t.statistic = c(t_res1$statistic, t_res2$statistic),
  p.value = c(t_res1$p.value, t_res2$p.value)
)
write.csv(t_result_df, file.path(output_dir, "DDX3X_group_t_test.csv"), row.names = FALSE)
cat("按DDX3X分组的t检验结果已保存。\n")

# ------------------- 7. 中介分析（Mediation） -------------------
cat("7. 执行中介分析...\n")
library(mediation)  # 确保已安装
# 中介1: DDX3X -> EIF2AK2_Editing_PC1 -> EIF2AK2
model_m1 <- lm(EIF2AK2_Editing_PC1 ~ DDX3X, data = merged_data)
model_y1 <- lm(EIF2AK2 ~ DDX3X + EIF2AK2_Editing_PC1, data = merged_data)
mediation1 <- mediate(model_m1, model_y1, treat = "DDX3X", mediator = "EIF2AK2_Editing_PC1",
                      boot = TRUE, sims = 1000)
summary(mediation1)
saveRDS(mediation1, file.path(output_dir, "mediation_EIF2AK2.rds"))

# 中介2: DDX3X -> PRKAR2A_Editing_PC1 -> PRKAR2A
model_m2 <- lm(PRKAR2A_Editing_PC1 ~ DDX3X, data = merged_data)
model_y2 <- lm(PRKAR2A ~ DDX3X + PRKAR2A_Editing_PC1, data = merged_data)
mediation2 <- mediate(model_m2, model_y2, treat = "DDX3X", mediator = "PRKAR2A_Editing_PC1",
                      boot = TRUE, sims = 1000)
summary(mediation2)
saveRDS(mediation2, file.path(output_dir, "mediation_PRKAR2A.rds"))
cat("中介分析结果已保存。\n")

# ------------------- 8. 免疫亚型分析（需要mmc2.xlsx） -------------------
cat("8. 读取免疫亚型数据...\n")
if (!file.exists(clinical_immune_file)) {
  cat("警告: 未找到免疫亚型文件", clinical_immune_file, "跳过免疫亚型分析。\n")
} else {
  type <- read_excel(clinical_immune_file)
  lung_type <- type[type$`TCGA Study` %in% c("LUAD", "LUSC"), ]
  # 将表达数据的行名（样本ID）转换为匹配格式（如 TCGA-XX-XXXX）
  sample_names <- rownames(selected_rows_df)
  sample_names <- gsub("\\.01A$", "", sample_names)
  sample_names <- gsub("\\.", "-", sample_names)
  rownames(selected_rows_df) <- sample_names
  
  merged_immune <- merge(selected_rows_df, lung_type, by.x = 0, by.y = "TCGA Participant Barcode", all = FALSE)
  if (nrow(merged_immune) == 0) {
    cat("没有样本匹配到免疫亚型数据，跳过。\n")
  } else {
    merged_immune <- merged_immune[, c("DDX3X", "EIF2AK2", "PRKAR2A", "Immune Subtype")]
    colnames(merged_immune)[4] <- "Immune_Subtype"
    merged_immune <- merged_immune[!is.na(merged_immune$Immune_Subtype), ]
    merged_immune$group <- ifelse(merged_immune$Immune_Subtype == "C2", "C2", "Other")
    
    genes <- c("DDX3X", "EIF2AK2", "PRKAR2A")
    # Wilcoxon检验
    wilcox_res <- lapply(genes, function(g) {
      test <- wilcox.test(as.formula(paste(g, "~ group")), data = merged_immune)
      data.frame(Gene = g, Method = "Wilcoxon", W_statistic = test$statistic, p_value = test$p.value)
    }) %>% bind_rows()
    write.csv(wilcox_res, file.path(output_dir, "immune_subtype_wilcox.csv"), row.names = FALSE)
    
    # t检验
    t_res_immune <- lapply(genes, function(g) {
      test <- t.test(as.formula(paste(g, "~ group")), data = merged_immune)
      data.frame(Gene = g, Method = "t-test", t_statistic = test$statistic, p_value = test$p.value)
    }) %>% bind_rows()
    write.csv(t_res_immune, file.path(output_dir, "immune_subtype_ttest.csv"), row.names = FALSE)
    
    # 绘制小提琴图
    plot_violin <- function(gene) {
      p <- ggplot(merged_immune, aes(x = group, y = .data[[gene]], fill = group)) +
        geom_violin(alpha = 0.7, trim = FALSE) +
        geom_boxplot(width = 0.1, fill = "white", alpha = 0.8) +
        geom_point(position = position_jitter(width = 0.1), size = 2, alpha = 0.6) +
        scale_fill_viridis_d(option = "D", begin = 0.3, end = 0.7) +
        stat_compare_means(method = "wilcox.test", label = "p.signif",
                           label.y = max(merged_immune[[gene]]) * 1.1) +
        labs(title = gene, x = NULL, y = "Expression") +
        theme_minimal(base_size = 14) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5),
              axis.text.x = element_text(angle = 45, hjust = 1))
      ggsave(file.path(output_dir, paste0("immune_subtype_", gene, ".png")), p, width = 6, height = 5, dpi = 300)
      print(p)
    }
    for (g in genes) plot_violin(g)
    cat("免疫亚型分析完成。\n")
  }
}

cat("Figure4 分析完成！结果保存在:", output_dir, "\n")