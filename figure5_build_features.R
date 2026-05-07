#!/usr/bin/env Rscript
# 构建四组特征集：临床、临床+RNA编辑、临床+表达、临床+两者
# 用法：Rscript build_features.R [项目根目录]

args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) > 0) args[1] else "."

# 输入文件
module_pc_file <- file.path(base_dir, "Figure2_output", "module_PC1-10.csv")
exp_file       <- file.path(base_dir, "Expression", "Exp_TPM_data.csv")
clinical_file  <- file.path(base_dir, "Expression", "combined_clinical_01A.tsv")

# 输出目录
output_dir <- file.path(base_dir, "Figure5", "features")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 加载包
library(dplyr)
library(tidyr)

# ---------- 1. 读取数据 ----------
cat("读取模块PC分数...\n")
module_PC <- read.csv(module_pc_file, row.names = 1)

cat("读取表达数据...\n")
expression_data <- read.csv(exp_file, row.names = 1, check.names = FALSE)
Exp <- as.data.frame(t(expression_data))

cat("读取临床数据...\n")
clinical <- read.csv(clinical_file, sep = "\t", row.names = 1)

# ---------- 2. 筛选肿瘤样本并提取N分期 ----------
clinical_n <- clinical[, "ajcc_pathologic_n.diagnoses", drop = FALSE]
clinical_n <- clinical_n[clinical_n$ajcc_pathologic_n.diagnoses %in% c("N0","N1","N2","N3"), , drop = FALSE]

merged_df <- merge(clinical_n, Exp, by = 0)
rownames(merged_df) <- merged_df$Row.names
merged_df$Row.names <- NULL
merged_df$n_status <- ifelse(merged_df$ajcc_pathologic_n.diagnoses == "N0", 0, 1)
merged_df$ajcc_pathologic_n.diagnoses <- NULL

# 去除零方差基因
gene_columns <- setdiff(colnames(merged_df), "n_status")
gene_variances <- apply(merged_df[, gene_columns], 2, var, na.rm = TRUE)
zero_var_genes <- names(gene_variances[gene_variances == 0])
filtered_genes <- setdiff(gene_columns, zero_var_genes)
merged_df <- merged_df[, c(filtered_genes, "n_status")]

# t检验筛选显著差异基因 (N0 vs N+)
significant_genes <- c()
for (gene in filtered_genes) {
  g0 <- merged_df[merged_df$n_status == 0, gene]
  g1 <- merged_df[merged_df$n_status == 1, gene]
  t_test <- t.test(g0, g1, var.equal = FALSE)
  if (t_test$p.value < 0.05) significant_genes <- c(significant_genes, gene)
}
filtered_expression <- Exp[, significant_genes, drop = FALSE]
cat("显著差异基因数:", ncol(filtered_expression), "\n")

# ---------- 3. 计算表达谱PCA（前30个PC） ----------
common_samples <- intersect(rownames(module_PC), rownames(filtered_expression))
Exp_01A <- as.matrix(filtered_expression[common_samples, , drop = FALSE])
pca_res <- prcomp(Exp_01A, center = TRUE, scale. = FALSE)
pc_scores <- pca_res$x[, 1:30, drop = FALSE]
colnames(pc_scores) <- paste0("Exp_PC", 1:30)

# ---------- 4. 清理临床变量 ----------
clean_clinical <- clinical %>%
  filter(!is.na(ajcc_pathologic_n.diagnoses), ajcc_pathologic_n.diagnoses != "NX") %>%
  mutate(
    gender = ifelse(gender.demographic == "male", 0, 1),
    age_years = round(age_at_diagnosis.diagnoses / 365.25, 2),
    subtype = ifelse(project_id.project == "TCGA-LUAD", 0, 1),
    t_stage = case_when(
      ajcc_pathologic_t.diagnoses %in% c("T1","T1a","T1b") ~ 1,
      ajcc_pathologic_t.diagnoses %in% c("T2","T2a","T2b") ~ 2,
      ajcc_pathologic_t.diagnoses == "T3" ~ 3,
      ajcc_pathologic_t.diagnoses == "T4" ~ 4
    ),
    metastasis = ifelse(ajcc_pathologic_n.diagnoses == "N0", 0, 1)
  ) %>%
  select(gender, age_years, subtype, t_stage, metastasis) %>%
  na.omit()

# ---------- 5. 取三者的共同样本 ----------
s1 <- rownames(module_PC)
s2 <- rownames(pc_scores)
s3 <- rownames(clean_clinical)
common <- Reduce(intersect, list(s1, s2, s3))
module_PC_common <- module_PC[common, , drop = FALSE]
exp_pc_common    <- pc_scores[common, , drop = FALSE]
clinical_common  <- clean_clinical[common, , drop = FALSE]

# ---------- 6. 筛选与转移显著相关的编辑模块特征和表达PC ----------
filter_data <- cbind(metastasis = clinical_common$metastasis, module_PC_common, exp_pc_common)

selected_edit <- c()
for (col in colnames(module_PC_common)) {
  g0 <- filter_data[filter_data$metastasis == 0, col]
  g1 <- filter_data[filter_data$metastasis == 1, col]
  if (length(g0) >= 2 & length(g1) >= 2) {
    tt <- t.test(g0, g1, var.equal = FALSE)
    if (!is.na(tt$p.value) && tt$p.value < 0.05) selected_edit <- c(selected_edit, col)
  }
}
selected_exp <- c()
for (col in colnames(exp_pc_common)) {
  g0 <- filter_data[filter_data$metastasis == 0, col]
  g1 <- filter_data[filter_data$metastasis == 1, col]
  if (length(g0) >= 2 & length(g1) >= 2) {
    tt <- t.test(g0, g1, var.equal = FALSE)
    if (!is.na(tt$p.value) && tt$p.value < 0.05) selected_exp <- c(selected_exp, col)
  }
}
cat("筛选后编辑模块特征数:", length(selected_edit), "\n")
cat("筛选后表达PC特征数:", length(selected_exp), "\n")

# ---------- 7. 构建四组特征集（含样本ID列） ----------
clinical_features <- clinical_common[, c("gender","age_years","subtype","t_stage"), drop = FALSE]
edit_features     <- module_PC_common[, selected_edit, drop = FALSE]
exp_features      <- exp_pc_common[, selected_exp, drop = FALSE]

# 添加样本ID列（用于Python模型）
add_sample_col <- function(df, samples) {
  df <- cbind(samples = samples, df)
  rownames(df) <- NULL
  return(df)
}

group1 <- add_sample_col(clinical_features, common)
group2 <- add_sample_col(cbind(clinical_features, edit_features), common)
group3 <- add_sample_col(cbind(clinical_features, exp_features), common)
group4 <- add_sample_col(cbind(clinical_features, edit_features, exp_features), common)

# 写入CSV（label列固定为"metastasis"）
write.csv(group1, file.path(output_dir, "Group1_Clinical.csv"), row.names = FALSE)
write.csv(group2, file.path(output_dir, "Group2_Clinical_Edit.csv"), row.names = FALSE)
write.csv(group3, file.path(output_dir, "Group3_Clinical_Exp.csv"), row.names = FALSE)
write.csv(group4, file.path(output_dir, "Group4_All_Features.csv"), row.names = FALSE)

cat("特征集构建完成！保存在:", output_dir, "\n")