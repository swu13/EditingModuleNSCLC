#!/usr/bin/env Rscript
# Figure3 分析脚本：特定模块（默认14）的基因网络、免疫相关性（nTreg/iTreg/Tex）、淋巴结分期箱线图
# 用法：Rscript figure3.R [项目根目录] [模块ID]

args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) > 0) args[1] else "."
target_cluster <- if (length(args) > 1) as.numeric(args[2]) else 14

# 输入文件路径
editing_file          <- file.path(base_dir, "Editing", "Filled_lung_0.2.txt")
gene_anno_file        <- file.path(base_dir, "Results", "editing_events_unique_gene_annotation.txt")
merged_cluster_file   <- file.path(base_dir, "Results", "merged_df_filtered.rds")
filtered_edges_file   <- file.path(base_dir, "Results", "filtered_edges.csv")
module_vectors_file   <- file.path(base_dir, "Figure2_output", "module_vectors.csv")
immune_file           <- file.path(base_dir, "Results", "ImmuCellAI.csv")
clinical_file         <- file.path(base_dir, "Expression", "combined_clinical_01A.tsv")

output_dir <- file.path(base_dir, "Figure3_output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

library(igraph)
library(ggplot2)
library(ggpubr)

# 辅助函数：标准化名称（空格→点）
normalize_name <- function(x) gsub(" ", ".", x)

# ------------------- 1. 读取数据 -------------------
cat("1. 读取RNA编辑数据...\n")
rna_editing_data <- read.table(editing_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cat("2. 读取基因注释表...\n")
gene_annotations <- read.table(gene_anno_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cat("3. 读取模块聚类结果...\n")
merged_df_filtered <- readRDS(merged_cluster_file)

cat("4. 读取网络边表...\n")
filtered_edges <- read.csv(filtered_edges_file, stringsAsFactors = FALSE)

cat("5. 读取模块PC1矩阵...\n")
module_vectors <- read.csv(module_vectors_file, row.names = 1)

# 匹配目标模块的PC1列
pattern <- paste0("Module[ .]", target_cluster, "[ .]?PC1?")
matched_cols <- grep(pattern, colnames(module_vectors), value = TRUE, ignore.case = TRUE)
if (length(matched_cols) == 0) {
  pattern2 <- paste0("Module[ .]", target_cluster)
  potential <- grep(pattern2, colnames(module_vectors), value = TRUE, ignore.case = TRUE)
  pc1_candidates <- potential[grepl("PC1", potential, ignore.case = TRUE)]
  if (length(pc1_candidates) >= 1) matched_cols <- pc1_candidates[1]
  else stop("无法找到模块 ", target_cluster, " 的PC1列")
}
module_pc1 <- module_vectors[, matched_cols[1], drop = FALSE]
colnames(module_pc1) <- "Module_PC1"
cat("使用的PC1列名:", matched_cols[1], "\n")
module_name_for_file <- gsub(" PC1$", "", matched_cols[1])

# ------------------- 2. 提取目标模块的编辑事件 -------------------
module_info <- merged_df_filtered[merged_df_filtered$cluster_id == target_cluster, ]
node_indices <- as.numeric(unlist(module_info$nodes))
module_events <- rna_editing_data$Editing[node_indices]
cat("模块", target_cluster, "包含编辑事件数:", length(module_events), "\n")

# 保存模块表达矩阵
module_expr <- t(rna_editing_data[node_indices, -1, drop = FALSE])
colnames(module_expr) <- module_events
write.csv(module_expr, file.path(output_dir, paste0("Module_", target_cluster, "_editing_expression.csv")), row.names = TRUE)

# ------------------- 3. 注释编辑事件到基因 -------------------
gene_map <- setNames(gene_annotations$gene_name, gene_annotations$Editing)
event_gene <- sapply(module_events, function(ev) ifelse(ev %in% names(gene_map), gene_map[[ev]], NA))
event_gene_df <- data.frame(Editing = module_events, gene_name = event_gene, stringsAsFactors = FALSE)
event_gene_df <- event_gene_df[!is.na(event_gene_df$gene_name), ]
split_rows <- strsplit(event_gene_df$gene_name, ",")
event_gene_long <- data.frame(
  Editing = rep(event_gene_df$Editing, sapply(split_rows, length)),
  gene_name = unlist(split_rows), stringsAsFactors = FALSE
)
genes <- unique(event_gene_long$gene_name)
gene_event_counts <- table(event_gene_long$gene_name)
gene_event_df <- data.frame(gene = names(gene_event_counts), event_count = as.vector(gene_event_counts))
write.csv(gene_event_df, file.path(output_dir, paste0("Module_", target_cluster, "_gene_event_counts.csv")), row.names = FALSE)
cat("模块", target_cluster, "涉及的基因数:", length(genes), "\n")

# ------------------- 4. 构建基因网络 -------------------
cat("构建基因-基因网络...\n")
if (length(genes) < 2) {
  edges_list <- data.frame(start = character(), end = character(), weight = numeric())
} else {
  n_genes <- length(genes)
  edge_weights <- matrix(0, nrow = n_genes, ncol = n_genes, dimnames = list(genes, genes))
  gene_to_edits <- split(event_gene_long$Editing, event_gene_long$gene_name)
  converted_edges <- data.frame(
    edit1 = rna_editing_data$Editing[filtered_edges$row],
    edit2 = rna_editing_data$Editing[filtered_edges$col],
    weight = filtered_edges$weight
  )
  for (i in 1:(n_genes-1)) {
    g1 <- genes[i]
    edits1 <- gene_to_edits[[g1]]
    for (j in (i+1):n_genes) {
      g2 <- genes[j]
      edits2 <- gene_to_edits[[g2]]
      shared <- intersect(edits1, edits2)
      w <- length(shared)
      idx <- (converted_edges$edit1 %in% edits1 & converted_edges$edit2 %in% edits2) |
             (converted_edges$edit1 %in% edits2 & converted_edges$edit2 %in% edits1)
      w <- w + sum(converted_edges$weight[idx], na.rm = TRUE)
      edge_weights[i, j] <- w
      edge_weights[j, i] <- w
    }
  }
  edges_list <- data.frame(start = character(), end = character(), weight = numeric())
  for (i in 1:(n_genes-1)) {
    for (j in (i+1):n_genes) {
      w <- edge_weights[i, j]
      if (w > 0) edges_list <- rbind(edges_list, data.frame(start = genes[i], end = genes[j], weight = w))
    }
  }
}
write.csv(edges_list, file.path(output_dir, paste0("Module_", target_cluster, "_gene_network_edges.csv")), row.names = FALSE)
cat("基因网络边数:", nrow(edges_list), "\n")
if (nrow(edges_list) > 0) {
  gnet <- graph_from_data_frame(edges_list, directed = FALSE, vertices = data.frame(name = genes))
  V(gnet)$size <- log(gene_event_counts[V(gnet)$name] + 1) * 5
  pdf(file.path(output_dir, paste0("Module_", target_cluster, "_gene_network.pdf")), width = 12, height = 10)
  plot(gnet, layout = layout_with_fr, vertex.label.cex = 0.8, edge.width = sqrt(E(gnet)$weight),
       main = paste("Gene co-editing network - Module", target_cluster))
  dev.off()
  cat("基因网络图已保存。\n")
}

# ------------------- 5. 免疫相关性散点图（nTreg, iTreg, Tex） -------------------
if (file.exists(immune_file)) {
  cat("读取免疫浸润数据...\n")
  imm <- read.csv(immune_file, row.names = 1)
  common <- intersect(rownames(module_pc1), rownames(imm))
  if (length(common) >= 3) {
    merged_imm <- cbind(module_pc1[common, , drop = FALSE], imm[common, , drop = FALSE])
    interest_cells <- c("nTreg", "iTreg", "Tex")
    available <- intersect(interest_cells, colnames(merged_imm))
    if (length(available) > 0) {
      for (cell in available) {
        ct <- cor.test(merged_imm$Module_PC1, merged_imm[[cell]], method = "pearson")
        r <- round(ct$estimate, 3); p <- format(ct$p.value, scientific = TRUE, digits = 3)
        p_scatter <- ggplot(merged_imm, aes(x = Module_PC1, y = .data[[cell]])) +
          geom_point(alpha = 0.6, color = "#2ca02c") +
          geom_smooth(method = "lm", se = TRUE, color = "#d62728") +
          annotate("text", x = min(merged_imm$Module_PC1), y = max(merged_imm[[cell]]),
                   label = paste0("Pearson r = ", r, ", p = ", p), hjust = 0, vjust = 1, size = 4) +
          labs(title = paste("Module", target_cluster, "PC1 vs", cell),
               x = paste("Module", target_cluster, "PC1"), y = paste(cell, "Infiltration")) +
          theme_minimal()
        ggsave(file.path(output_dir, paste0("Module_", target_cluster, "_PC1_vs_", cell, ".png")),
               p_scatter, width = 6, height = 5, dpi = 300)
        cat("散点图已保存:", cell, "\n")
      }
    } else cat("未找到 nTreg/iTreg/Tex 列。\n")
  } else cat("共同样本不足。\n")
} else cat("未找到免疫文件，跳过免疫分析。\n")

# ------------------- 6. 模块PC1与淋巴结分期的关联（t检验和箱线图） -------------------
if (file.exists(clinical_file)) {
  cat("读取临床N分期数据，进行t检验和箱线图...\n")
  clin <- read.table(clinical_file, header=TRUE, sep="\t", row.names=1)
  common_samples <- intersect(rownames(module_pc1), rownames(clin))
  if (length(common_samples) > 0) {
    df <- data.frame(PC1 = module_pc1[common_samples, "Module_PC1"],
                     N = clin[common_samples, "ajcc_pathologic_n.diagnoses"])
    df <- df[!is.na(df$N) & df$N != "NX", ]
    df$group <- ifelse(df$N == "N0", "N0", ifelse(df$N %in% c("N1","N2","N3"), "N>0", NA))
    df <- df[!is.na(df$group), ]
    if (nrow(df) > 0 && length(unique(df$group)) == 2) {
      tt <- t.test(PC1 ~ group, data=df)
      t_res <- data.frame(Module = matched_cols[1], T_Value = tt$statistic, P_Value = tt$p.value,
                          Mean_N0 = tt$estimate[1], Mean_Npos = tt$estimate[2])
      write.csv(t_res, file.path(output_dir, paste0("Module_", target_cluster, "_Nstage_t_test.csv")), row.names=FALSE)
      cat("t检验结果:\n"); print(t_res)
      
      # 箱线图
      p_box <- ggplot(df, aes(x=group, y=PC1, fill=group)) +
        geom_boxplot(outlier.shape=NA, alpha=0.8) +
        geom_jitter(width=0.2, alpha=0.5, size=1.5) +
        stat_compare_means(method="t.test", label="p.format") +
        labs(title=paste("Module", target_cluster, "PC1 vs Lymph Node Status"),
             x="Lymph Node Status", y="PC1") +
        theme_minimal() + scale_fill_manual(values=c("#1f77b4","#ff7f0e"))
      ggsave(file.path(output_dir, paste0("Module_", target_cluster, "_Nstage_boxplot.png")), p_box, width=5, height=5, dpi=300)
      cat("箱线图已保存。\n")
    } else cat("分组后样本不足或仅有一组，无法进行t检验。\n")
  } else cat("没有共同样本。\n")
} else cat("未找到临床文件，跳过N分期分析。\n")

# ------------------- 7. 生存分析结果提取（修正匹配） -------------------
cox_file <- file.path(base_dir, "Figure2_output", "Cox_regression_results.csv")
if (file.exists(cox_file)) {
  cox_res <- read.csv(cox_file, stringsAsFactors = FALSE)
  cox_res$Module_norm <- normalize_name(cox_res$Module)
  target_norm <- normalize_name(matched_cols[1])
  mod_cox <- cox_res[cox_res$Module_norm == target_norm, ]
  if (nrow(mod_cox) == 0 && module_name_for_file != matched_cols[1]) {
    mod_cox <- cox_res[cox_res$Module_norm == normalize_name(module_name_for_file), ]
  }
  if (nrow(mod_cox) > 0) {
    write.csv(mod_cox[, !names(mod_cox) %in% "Module_norm"], 
              file.path(output_dir, paste0("Module_", target_cluster, "_Cox_results.csv")), row.names = FALSE)
    cat("模块", target_cluster, "Cox回归结果已保存。\n")
  } else cat("未找到模块", target_cluster, "的Cox结果。\n")
} else cat("未找到Cox文件。\n")

km_dir <- file.path(base_dir, "Figure2_output", "KM")
if (dir.exists(km_dir)) {
  km_files <- list.files(km_dir, pattern = paste0(gsub(" ", "[ .]", matched_cols[1]), ".*\\.png$"), ignore.case = TRUE)
  if (length(km_files) == 0 && module_name_for_file != matched_cols[1]) {
    km_files <- list.files(km_dir, pattern = paste0(gsub(" ", "[ .]", module_name_for_file), ".*\\.png$"), ignore.case = TRUE)
  }
  if (length(km_files) > 0) {
    file.copy(file.path(km_dir, km_files[1]), 
              file.path(output_dir, paste0("Module_", target_cluster, "_KM_curve.png")), overwrite = TRUE)
    cat("KM曲线已复制。\n")
  } else cat("未找到对应KM曲线。\n")
} else cat("未找到KM目录。\n")

cat("Figure3 分析完成！结果保存在:", output_dir, "\n")