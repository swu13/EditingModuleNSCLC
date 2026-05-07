#!/usr/bin/env Rscript
# 全套分析脚本：RNA编辑事件PCA、与临床N分期关联、干性分数与表达相关性
# 用法：Rscript figure6.R [项目根目录]
# 如果未指定项目根目录，默认使用当前工作目录

# 解析命令行参数
args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) > 0) args[1] else "."

# 设置输出目录
output_dir <- file.path(base_dir, "Figure6_output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ==================== 加载必要的R包 ====================
library(ggplot2)
library(ggpubr)
library(dplyr)

# ==================== 1. 读取RNA编辑数据并筛选目标事件 ====================
cat("1. 读取RNA编辑数据...\n")
rna_editing_file <- file.path(base_dir, "Editing/Filled_lung_0.2.txt")
rna_editing_data <- read.table(rna_editing_file, header = TRUE, sep = "\t")
cat("编辑事件总数:", nrow(rna_editing_data), "\n")

#AHR编辑事件
target_edits <- c(
  "chr7_17344769_+",
  "chr7_17344810_+",
  "chr7_17345061_+",
  "chr7_17345154_+",
  "chr7_17345168_+",
  "chr7_17345204_+"
)

extracted_data <- rna_editing_data[rna_editing_data$Editing %in% target_edits, ]
cat("目标编辑事件数量:", nrow(extracted_data), "\n")

# 转置数据：行为样本，列为编辑事件
rownames(extracted_data) <- extracted_data$Editing
data_no_firstcol <- extracted_data[, -1]
transposed_data <- as.data.frame(t(data_no_firstcol))

# ==================== 2. 对编辑事件进行PCA ====================
cat("2. 对编辑事件进行PCA...\n")
pca_result <- prcomp(transposed_data)
transposed_data$PC1 <- pca_result$x[, 1]
cat("PC1解释方差比例:", summary(pca_result)$importance[2, 1], "\n")

# ==================== 3. 读取临床N分期数据 ====================
cat("3. 读取临床N分期数据...\n")
clinical_file <- file.path(base_dir, "Expression/combined_clinical_01A.tsv")
phenotype_data <- read.table(clinical_file, header = TRUE, sep = "\t", row.names = 1)

# 合并PC1与临床信息
merged_data <- merge(transposed_data, phenotype_data, by = 0)
merged_data <- merged_data[, c("Row.names", "PC1", "ajcc_pathologic_n.diagnoses")]

# 清洗数据：去除缺失和NX
merged_data_not_missing_n <- merged_data[!is.na(merged_data$ajcc_pathologic_n.diagnoses), ]
merged_data_filtered_n <- merged_data_not_missing_n[merged_data_not_missing_n$ajcc_pathologic_n.diagnoses != "NX", ]
cat("N分期分布:\n")
print(table(merged_data_filtered_n$ajcc_pathologic_n.diagnoses))

# 创建二分类分组（N=0 vs N>0）
merged_data_filtered_n$n_group <- ifelse(merged_data_filtered_n$ajcc_pathologic_n.diagnoses == "N0", 
                                         "N=0", 
                                         ifelse(merged_data_filtered_n$ajcc_pathologic_n.diagnoses %in% c("N1", "N2", "N3"), 
                                                "N>0", 
                                                "Other"))
clean_data <- merged_data_filtered_n[merged_data_filtered_n$n_group != "Other", ]
cat("分组样本数: N=0:", sum(clean_data$n_group == "N=0"), ", N>0:", sum(clean_data$n_group == "N>0"), "\n")

# ==================== 4. t检验和箱线图（PC1 vs 淋巴结分期） ====================
cat("4. PC1与淋巴结分期的t检验...\n")
t_test_pc1 <- t.test(PC1 ~ n_group, data = clean_data)
print(t_test_pc1)

# 自动生成标题文本
t_stat <- round(t_test_pc1$statistic, 2)
p_val <- format(t_test_pc1$p.value, scientific = TRUE, digits = 3)
subtitle_text <- paste0("Welch's t-test: t = ", t_stat, ", p = ", p_val)

# 绘制箱线图
p4 <- ggplot(clean_data, aes(x = n_group, y = PC1, fill = n_group)) +
  geom_boxplot(width = 0.6, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.5, aes(color = n_group)) +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e")) +
  scale_color_manual(values = c("#1f77b4", "#ff7f0e")) +
  labs(
    title = "PC1 Distribution by Lymph Node Status",
    subtitle = subtitle_text,
    x = "Lymph Node Involvement",
    y = "Principal Component 1 (PC1)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "none",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) +
  stat_compare_means(
    method = "t.test",
    label = "p.format",
    label.y = max(clean_data$PC1, na.rm = TRUE) * 1.05,
    show.legend = FALSE
  )

# 保存图片
ggsave(file.path(output_dir, "PC1_by_N_status.png"), p4, width = 6, height = 5, dpi = 300)
print(p4)

# ==================== 5. 干性分数（CSscore）与N分期关联 ====================
cat("5. 读取干性分数数据...\n")
#从官网自行下载整理：https://bio-bigdata.hrbmu.edu.cn/CancerStemnessOnline/
cs_file <- file.path(base_dir, "Expression/CSscore1.csv")  # 请根据实际路径修改
CSscore_data <- read.csv(cs_file, row.names = 1)

# 合并干性分数与临床N分期
merged_cs <- merge(CSscore_data, phenotype_data, by = 0)
merged_cs <- merged_cs[, c("Row.names", "score", "ajcc_pathologic_n.diagnoses")]

# 同样清洗数据
merged_cs_not_missing <- merged_cs[!is.na(merged_cs$ajcc_pathologic_n.diagnoses), ]
merged_cs_filtered <- merged_cs_not_missing[merged_cs_not_missing$ajcc_pathologic_n.diagnoses != "NX", ]
merged_cs_filtered$n_group <- ifelse(merged_cs_filtered$ajcc_pathologic_n.diagnoses == "N0", 
                                     "N=0", 
                                     ifelse(merged_cs_filtered$ajcc_pathologic_n.diagnoses %in% c("N1", "N2", "N3"), 
                                            "N>0", 
                                            "Other"))
clean_cs <- merged_cs_filtered[merged_cs_filtered$n_group != "Other", ]

cat("干性分数与N分期的t检验...\n")
t_test_cs <- t.test(score ~ n_group, data = clean_cs)
print(t_test_cs)

subtitle_cs <- paste0("Welch's t-test: t = ", round(t_test_cs$statistic, 2), 
                      ", p = ", format(t_test_cs$p.value, scientific = TRUE, digits = 3))
p9 <- ggplot(clean_cs, aes(x = n_group, y = score, fill = n_group)) +
  geom_boxplot(width = 0.6, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.5, aes(color = n_group)) +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e")) +
  scale_color_manual(values = c("#1f77b4", "#ff7f0e")) +
  labs(title = "Stemness Score Distribution by Lymph Node Status",
       subtitle = subtitle_cs,
       x = "Lymph Node Involvement",
       y = "Stemness Score") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none",
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank()) +
  stat_compare_means(method = "t.test", label = "p.format",
                     label.y = max(clean_cs$score, na.rm = TRUE) * 1.05,
                     show.legend = FALSE)
ggsave(file.path(output_dir, "Stemness_by_N_status.png"), p9, width = 6, height = 5, dpi = 300)
print(p9)

# ==================== 6. 相关性分析：干性分数 vs MFSD2A表达 ====================
cat("6. 读取表达数据并提取MFSD2A...\n")
exp_file <- file.path(base_dir, "Expression/Exp_TPM_data.csv")
expression_data <- read.csv(exp_file, row.names = 1)
# 筛选肿瘤样本（ID含.01A）
keep_cols <- grep("^.*\\..*\\..*\\.01A$", colnames(expression_data))
new_expression_data <- expression_data[, keep_cols]
if ("MFSD2A" %in% rownames(new_expression_data)) {
  mfsd2a_exp <- data.frame(t(new_expression_data["MFSD2A", , drop = FALSE]))
  mfsd2a_exp$sample <- rownames(mfsd2a_exp)
  CSscore_data$sample <- rownames(CSscore_data)
  
  # 合并干性分数与MFSD2A表达
  merged_cor <- merge(CSscore_data, mfsd2a_exp, by = "sample")
  cat("干性分数与MFSD2A表达的相关性分析...\n")
  cor_test1 <- cor.test(merged_cor$score, merged_cor$MFSD2A, method = "pearson")
  print(cor_test1)
  
  # 可选：散点图
  p_cor <- ggplot(merged_cor, aes(x = score, y = MFSD2A)) +
    geom_point(alpha = 0.6, color = "#1f77b4") +
    geom_smooth(method = "lm", se = TRUE, color = "#ff7f0e") +
    annotate("text", x = min(merged_cor$score), y = max(merged_cor$MFSD2A),
             label = paste0("Pearson r = ", round(cor_test1$estimate, 3),
                            "\np = ", format(cor_test1$p.value, scientific = TRUE, digits = 3)),
             hjust = 0, vjust = 1, size = 4) +
    labs(title = "Correlation between Stemness Score and MFSD2A Expression",
         x = "Stemness Score", y = "MFSD2A Expression (log2TPM+1)") +
    theme_minimal()
  ggsave(file.path(output_dir, "Stemness_MFSD2A_correlation.png"), p_cor, width = 6, height = 5, dpi = 300)
  print(p_cor)
} else {
  cat("警告：MFSD2A 不在表达矩阵中，跳过相关性分析。\n")
}

# ==================== 7. PC1与MFSD2A表达相关性 ====================
# --- PC1 与 MFSD2A 相关性散点图 ---
common_samples <- intersect(rownames(transposed_data), rownames(mfsd2a_exp))
if (length(common_samples) > 2) {
    pc1_vec <- transposed_data[common_samples, "PC1"]
    mfsd2a_vec <- mfsd2a_exp[common_samples, "MFSD2A"]
    cor_test2 <- cor.test(pc1_vec, mfsd2a_vec, method = "pearson")
    cat("PC1与MFSD2A表达的相关性:\n")
    print(cor_test2)
    
    # 构建数据框用于绘图
    pc1_mfsd2a_df <- data.frame(
        PC1 = pc1_vec,
        MFSD2A = mfsd2a_vec
    )
    
    # 提取相关系数和 p 值
    r_val <- round(cor_test2$estimate, 3)
    p_val <- format(cor_test2$p.value, scientific = TRUE, digits = 3)
    
    # 绘制散点图 + 回归线
    p_pc1_mfsd2a <- ggplot(pc1_mfsd2a_df, aes(x = PC1, y = MFSD2A)) +
        geom_point(alpha = 0.6, color = "#2ca02c") +
        geom_smooth(method = "lm", se = TRUE, color = "#d62728") +
        annotate("text", 
                 x = min(pc1_mfsd2a_df$PC1, na.rm = TRUE), 
                 y = max(pc1_mfsd2a_df$MFSD2A, na.rm = TRUE),
                 label = paste0("Pearson r = ", r_val, "\np = ", p_val),
                 hjust = 0, vjust = 1, size = 4) +
        labs(title = "Correlation between PC1 and MFSD2A Expression",
             x = "Principal Component 1 (PC1)",
             y = "MFSD2A Expression (log2TPM+1)") +
        theme_minimal()
    
    # 保存图片
    ggsave(file.path(output_dir, "PC1_MFSD2A_correlation.png"), 
           p_pc1_mfsd2a, width = 6, height = 5, dpi = 300)
    print(p_pc1_mfsd2a)
} else {
    cat("共同样本不足，无法绘制 PC1 vs MFSD2A 散点图。\n")
}

cat("所有分析完成！图片已保存至:", output_dir, "\n")