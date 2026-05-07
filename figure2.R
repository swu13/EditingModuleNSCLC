#!/usr/bin/env Rscript
# Figure2 分析脚本：RNA编辑事件网络构建、模块识别、生存分析与临床关联
# 用法：Rscript figure2.R [项目根目录]

args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) > 0) args[1] else "."

# 输入文件路径
editing_file      <- file.path(base_dir, "Editing", "Filled_lung_0.2.txt")
gene_anno_file    <- file.path(base_dir, "Results", "editing_events_unique_gene_annotation.txt")
exp_file          <- file.path(base_dir, "Expression", "Exp_TPM_data.csv")
survival_file     <- file.path(base_dir, "Expression", "TCGA-lung.survival.tsv")
clinical_file     <- file.path(base_dir, "Expression", "combined_clinical_01A.tsv")

# 输出目录
output_dir <- file.path(base_dir, "Figure2_output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
km_dir <- file.path(output_dir, "KM")
if (!dir.exists(km_dir)) dir.create(km_dir)

# 加载包
library(igraph)
library(survival)
library(survminer)
library(ggplot2)

# ------------------- 1. 读取编辑数据并计算相关矩阵 -------------------
cat("1. 读取RNA编辑数据...\n")
rna_editing_data <- read.table(editing_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat("编辑事件数:", nrow(rna_editing_data), "样本数:", ncol(rna_editing_data)-1, "\n")

cat("计算绝对相关系数矩阵...\n")
cor_matrix <- abs(cor(t(rna_editing_data[, -1]), method = "pearson"))
correlation_threshold <- 0.2
edges <- which(cor_matrix > correlation_threshold & upper.tri(cor_matrix), arr.ind = TRUE)
edges_df <- as.data.frame(edges)
edges_df$weight <- cor_matrix[edges]
cat("初始边数量:", nrow(edges_df), "\n")

# ------------------- 2. 网络构建与PageRank筛选 -------------------
cat("构建igraph网络...\n")
net <- graph_from_data_frame(d = edges_df, directed = FALSE)
page_rank_values <- page_rank(net)$vector
filtered_nodes <- V(net)$name[page_rank_values > quantile(page_rank_values)[4]]
filtered_edges <- subset(edges_df, row %in% filtered_nodes & col %in% filtered_nodes)
write.csv(filtered_edges, file = file.path(base_dir, "Results", "filtered_edges.csv"), row.names = FALSE)
net_filtered <- graph_from_data_frame(d = filtered_edges, directed = FALSE)
cat("筛选后节点数:", vcount(net_filtered), "边数:", ecount(net_filtered), "\n")

# ------------------- 3. 注释编辑事件到基因（使用已有注释表） -------------------
cat("读取基因注释表...\n")
gene_annotations <- read.table(gene_anno_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
filtered_node_indices <- as.integer(V(net_filtered)$name)
filtered_editing_sites <- rna_editing_data$Editing[filtered_node_indices]
node_names <- filtered_editing_sites

gene_map <- setNames(gene_annotations$gene_name, gene_annotations$Editing)
result_list <- vector("list", length(node_names))
for (i in seq_along(node_names)) {
  result_list[[i]] <- ifelse(node_names[i] %in% names(gene_map), gene_map[node_names[i]], NA)
  if (i %% 100 == 0) cat("注释进度:", i, "/", length(node_names), "\n")
}
annotation_results <- data.frame(editing_site = node_names, gene_name = unlist(result_list), stringsAsFactors = FALSE)
write.table(annotation_results, file.path(output_dir, "editing_sites_gene_annotation.txt"), row.names = FALSE, sep = "\t", quote = FALSE)

split_genes <- strsplit(annotation_results$gene_name, ",")
gene_list <- data.frame(
  editing_site = rep(annotation_results$editing_site, sapply(split_genes, length)),
  gene_name = unlist(split_genes),
  stringsAsFactors = FALSE
)
unique_genes <- unique(gene_list$gene_name)
write.table(data.frame(unique_genes), file.path(output_dir, "unique_genes_list.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)

# ------------------- 4. 表达差异分析（肿瘤 vs 正常） -------------------
cat("读取表达数据...\n")
expression_data <- read.csv(exp_file, row.names = 1, check.names = FALSE)
available_genes <- unique_genes[unique_genes %in% rownames(expression_data)]
cat("可用基因数:", length(available_genes), "/", length(unique_genes), "\n")
missing_genes <- setdiff(unique_genes, available_genes)
if (length(missing_genes) > 0) {
  cat("缺失基因:\n"); print(missing_genes)
}

genes_exp <- expression_data[available_genes, , drop = FALSE]
keep <- grepl("(01A|11A)$", colnames(genes_exp))
gen_sub <- genes_exp[, keep, drop = FALSE]
group <- ifelse(grepl("01A$", colnames(gen_sub)), "T", "N")
cat("肿瘤样本:", sum(group=="T"), "正常样本:", sum(group=="N"), "\n")

t_res <- apply(gen_sub, 1, function(x) {
  tt <- t.test(x[group=="T"], x[group=="N"])
  data.frame(mean_T = tt$estimate[1], mean_N = tt$estimate[2],
             logFC = log2(tt$estimate[1]/tt$estimate[2]), p_value = tt$p.value)
})
t_df <- do.call(rbind, t_res)
sig_count <- sum(t_df$p_value < 0.05, na.rm = TRUE)
sig_prop <- sig_count / length(unique_genes)          # 分母改为所有注释到的基因总数
cat("显著差异基因比例 (p<0.05)：", sig_prop, "\n")

# ------------------- 5. 网络聚类（Walktrap）及递归划分 -------------------
cat("执行Walktrap聚类...\n")
cwt <- cluster_walktrap(net_filtered, steps = 6)
cluster_sizes <- sizes(cwt)
large_clusters <- which(cluster_sizes >= 100)

# 递归函数
process_clusters <- function(graph, cwt, current_id, level, res_df) {
  subg <- induced_subgraph(graph, which(membership(cwt) == current_id))
  sub_cwt <- cluster_walktrap(subg, steps = 6)
  sub_sizes <- sizes(sub_cwt)
  cat("Level", level, "Cluster", current_id, "子簇数:", length(sub_sizes), "\n")
  if (length(sub_sizes) > 1) {
    for (sid in which(sub_sizes < 100)) {
      nodes_sub <- V(subg)$name[which(membership(sub_cwt) == sid)]
      res_df <- rbind(res_df, data.frame(Level = level, Cluster_ID = current_id,
                                         Subcluster_ID = sid, Subcluster_sizes = sub_sizes[sid],
                                         Nodes = I(list(nodes_sub))))
    }
    for (sid in which(sub_sizes >= 100)) {
      res_df <- process_clusters(subg, sub_cwt, sid, level + 1, res_df)
    }
  }
  return(res_df)
}

# 存储小簇
small_clusters_df <- data.frame(cluster_id = integer(), size = integer(), nodes = I(list()))
for (cid in which(cluster_sizes < 100)) {
  small_clusters_df <- rbind(small_clusters_df,
                             data.frame(cluster_id = cid, size = cluster_sizes[cid],
                                        nodes = I(list(V(net_filtered)$name[membership(cwt) == cid]))))
}
# 处理大簇
result_df <- data.frame(Level = integer(), Cluster_ID = integer(), Subcluster_ID = integer(),
                        Subcluster_sizes = integer(), Nodes = I(list()))
for (cid in large_clusters) {
  result_df <- process_clusters(net_filtered, cwt, cid, 1, result_df)
}
# 合并所有大小<100的簇
merged_df <- rbind(
  data.frame(cluster_id = small_clusters_df$cluster_id, size = small_clusters_df$size, nodes = small_clusters_df$nodes),
  data.frame(cluster_id = result_df$Cluster_ID, size = result_df$Subcluster_sizes, nodes = result_df$Nodes)
)
merged_df <- merged_df[order(merged_df$cluster_id), ]
merged_df$cluster_id <- seq_len(nrow(merged_df))
cat("最终模块数（大小<100）:", nrow(merged_df), "\n")
merged_df_filtered <- merged_df[merged_df$size > 30, ]
saveRDS(merged_df_filtered, file.path(base_dir, "Results", "merged_df_filtered.rds"))
cat("用于后续分析的模块数（大小>30）:", nrow(merged_df_filtered), "\n")
print(merged_df_filtered[, c("cluster_id", "size")])

# ------------------- 6. 计算每个模块的PC（前10个主成分） -------------------
cat("计算模块PCA向量...\n")
n_samples <- ncol(rna_editing_data) - 1
module_vectors <- data.frame(matrix(ncol = nrow(merged_df_filtered) * 10, nrow = n_samples))
rownames(module_vectors) <- colnames(rna_editing_data)[-1]
for (i in 1:nrow(merged_df_filtered)) {
  node_indices <- as.numeric(merged_df_filtered$nodes[[i]])
  module_data <- rna_editing_data[node_indices, -1, drop = FALSE]
  pca <- prcomp(t(module_data), center = TRUE, scale. = FALSE)
  pc_scores <- pca$x[, 1:min(10, ncol(pca$x)), drop = FALSE]
  if (ncol(pc_scores) < 10) {
    pc_scores <- cbind(pc_scores, matrix(NA, nrow = nrow(pc_scores), ncol = 10 - ncol(pc_scores)))
  }
  idx_start <- (i-1)*10 + 1
  module_vectors[, idx_start:(idx_start+9)] <- pc_scores
  colnames(module_vectors)[idx_start:(idx_start+9)] <- paste("Module", merged_df_filtered$cluster_id[i], "PC", 1:10)
}
write.csv(module_vectors, file.path(output_dir, "module_PC1-10.csv"), row.names = TRUE)

pc1_indices <- seq(1, ncol(module_vectors), by = 10)
module_vectors_pc1 <- module_vectors[, pc1_indices, drop = FALSE]
colnames(module_vectors_pc1) <- gsub(" PC1$", "", colnames(module_vectors_pc1))
write.csv(module_vectors_pc1, file.path(output_dir, "module_vectors.csv"), row.names = TRUE)

# ------------------- 7. 生存分析（KM曲线和Cox回归） -------------------
cat("读取生存数据...\n")
survival_data <- read.table(survival_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
survival_data$sample <- gsub("-", ".", survival_data$sample)
survival_data <- survival_data[!is.na(survival_data$OS) & !is.na(survival_data$OS.time) & survival_data$OS.time > 0, ]
survival_data <- survival_data[grepl("01A$", survival_data$sample), ]

common_samps <- intersect(rownames(module_vectors_pc1), survival_data$sample)
module_surv <- module_vectors_pc1[common_samps, , drop = FALSE]
surv_data <- survival_data[match(common_samps, survival_data$sample), ]
stopifnot(all(rownames(module_surv) == surv_data$sample))

cat("绘制KM曲线和Cox回归...\n")
cox_results <- data.frame(Module = character(), HR = numeric(), Lower_CI = numeric(),
                          Upper_CI = numeric(), P_Value = numeric(), stringsAsFactors = FALSE)

for (mod in colnames(module_surv)) {
  median_val <- median(module_surv[[mod]], na.rm = TRUE)
  group <- ifelse(module_surv[[mod]] <= median_val, "Low", "High")
  surv_data$group <- factor(group, levels = c("Low", "High"))
  
  # 正确的 survfit 调用：提供 data 参数
  km_fit <- survfit(Surv(OS.time, OS) ~ group, data = surv_data)
  
  # ggsurvplot 同时需要 fit 和 data
  km_plot <- ggsurvplot(km_fit, data = surv_data, pval = TRUE,
                        title = paste("KM Curve -", mod),
                        legend.labs = c("Low", "High"), legend.title = mod)
  ggsave(file.path(km_dir, paste0("KM_", mod, ".png")), plot = km_plot$plot, width = 8, height = 6, dpi = 300)
  
  # Cox 回归
  cox_model <- coxph(Surv(OS.time, OS) ~ module_surv[[mod]], data = surv_data)
  summ <- summary(cox_model)
  cox_results <- rbind(cox_results, data.frame(Module = mod,
                                               HR = summ$coefficients[1, "exp(coef)"],
                                               Lower_CI = summ$conf.int[1, "lower .95"],
                                               Upper_CI = summ$conf.int[1, "upper .95"],
                                               P_Value = summ$coefficients[1, "Pr(>|z|)"]))
}
write.csv(cox_results, file.path(output_dir, "Cox_regression_results.csv"), row.names = FALSE)
print(cox_results)

# ------------------- 8. 模块PC1与淋巴结分期关联 -------------------
cat("读取临床N分期数据...\n")
phenotype <- read.table(clinical_file, header = TRUE, sep = "\t", row.names = 1)
merged_clin <- merge(module_vectors_pc1, phenotype, by = 0)
rownames(merged_clin) <- merged_clin$Row.names
merged_clin$Row.names <- NULL
merged_clin <- merged_clin[!is.na(merged_clin$ajcc_pathologic_n.diagnoses) & merged_clin$ajcc_pathologic_n.diagnoses != "NX", ]
group0 <- merged_clin[merged_clin$ajcc_pathologic_n.diagnoses == "N0", , drop = FALSE]
group_pos <- merged_clin[merged_clin$ajcc_pathologic_n.diagnoses %in% c("N1","N2","N3"), , drop = FALSE]
mod_cols <- grep("Module", colnames(merged_clin), value = TRUE)
t_test_res <- data.frame(Module = character(), T_Value = numeric(), P_Value = numeric())
for (mod in mod_cols) {
  tt <- t.test(group0[[mod]], group_pos[[mod]])
  t_test_res <- rbind(t_test_res, data.frame(Module = mod, T_Value = tt$statistic, P_Value = tt$p.value))
}
write.csv(t_test_res, file.path(output_dir, "Module_N_stage_t_test.csv"), row.names = FALSE)
print(t_test_res)

cat("Figure2 分析完成！结果保存在:", output_dir, "\n")