#!/bin/bash
# 数据预处理脚本：下载 TCGA 数据，合并表达谱、临床、生存信息，处理编辑事件并注释到基因
# 用法：./preprocess.sh [项目根目录]
# 如果不指定，默认使用当前目录

set -e  # 遇到错误立即退出

# 设置项目根目录
if [ -z "$1" ]; then
    ProjectFold="$(pwd)"
else
    ProjectFold="$1"
fi

echo "项目目录: $ProjectFold"
cd "$ProjectFold" || exit 1

# 创建子目录
mkdir -p "$ProjectFold/Expression"
mkdir -p "$ProjectFold/Editing"
mkdir -p "$ProjectFold/Results"

# ==================== 1. 表达与临床数据下载 ====================
cd "$ProjectFold/Expression" || exit 1

echo "1. 下载 TCGA 数据文件..."
wget -q https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LUAD.star_tpm.tsv.gz
wget -q https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LUSC.star_tpm.tsv.gz
wget -q https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LUAD.clinical.tsv.gz
wget -q https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LUSC.clinical.tsv.gz
wget -q https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LUAD.survival.tsv.gz
wget -q https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LUSC.survival.tsv.gz
wget -q https://gdc-hub.s3.us-east-1.amazonaws.com/download/gencode.v36.annotation.gtf.gene.probemap
wget -q https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.chr_patch_hapl_scaff.annotation.gff3.gz
gunzip -f *.gz

# ==================== 2. 合并表达谱、临床、生存数据 ====================
echo "2. 合并表达谱、临床信息和生存数据..."
Rscript -e "
    library(dplyr)

    # 读取表达数据
    luad_exp <- read.delim('TCGA-LUAD.star_tpm.tsv')
    lusc_exp <- read.delim('TCGA-LUSC.star_tpm.tsv')
    merged_exp <- merge(luad_exp, lusc_exp, by = 'Ensembl_ID')

    # 读取 probemap（基因 ID -> 基因名）
    probemap <- read.delim('gencode.v36.annotation.gtf.gene.probemap', header = TRUE, stringsAsFactors = FALSE)
    colnames(probemap)[1:2] <- c('Ensembl_ID', 'gene_name')

    # 合并并去除 Ensembl_ID 列，将 gene_name 移到第一列
    exp_with_gene <- merge(merged_exp, probemap, by = 'Ensembl_ID')
    exp_with_gene <- exp_with_gene %>%
        select(-Ensembl_ID) %>%
        select(gene_name, everything())

    # 转置：行名为基因名，列为样本
    Exp <- as.matrix(exp_with_gene[, -1])
    rownames(Exp) <- make.unique(exp_with_gene$gene_name)

    # 保存完整表达矩阵
    write.csv(Exp, file = 'Exp_tpm_full.csv', row.names = TRUE)

    # 筛选肿瘤样本（样本书包含 .01A）
    tumor_cols <- grep('\\.01A\$', colnames(Exp), value = TRUE)
    exp_tumor <- Exp[, tumor_cols, drop = FALSE]
    write.csv(exp_tumor, file = 'Exp_tpm_01A_data.csv', row.names = TRUE)

    # 处理临床数据
    luad_clin <- read.csv('TCGA-LUAD.clinical.tsv', sep = '\t', quote = '', fill = TRUE, stringsAsFactors = FALSE)
    lusc_clin <- read.csv('TCGA-LUSC.clinical.tsv', sep = '\t', quote = '', fill = TRUE, stringsAsFactors = FALSE)
    target_cols <- c('sample', 'gender.demographic', 'age_at_diagnosis.diagnoses',
                     'ajcc_pathologic_t.diagnoses', 'ajcc_pathologic_n.diagnoses')
    luad_clin <- luad_clin[, target_cols]
    lusc_clin <- lusc_clin[, target_cols]
    combined_clin <- rbind(luad_clin, lusc_clin)
    combined_clin <- combined_clin[grepl('01A\$', combined_clin$sample), ]
    combined_clin$sample <- gsub('-', '.', combined_clin$sample)
    write.table(combined_clin, file = 'combined_clinical_01A.tsv', sep = '\t', row.names = FALSE, quote = FALSE)

    # 处理生存数据
    luad_surv <- read.delim('TCGA-LUAD.survival.tsv')[, 1:3]
    lusc_surv <- read.delim('TCGA-LUSC.survival.tsv')[, 1:3]
    lung_surv <- rbind(luad_surv, lusc_surv)
    write.table(lung_surv, file = 'TCGA-lung.survival.tsv', sep = '\t', row.names = FALSE, quote = FALSE)
"

# ==================== 3. 编辑事件预处理 ====================
cd "$ProjectFold/Editing" || exit 1

echo "3. 过滤和填充编辑事件数据..."
Rscript -e "
    library(dplyr)

    luad_edit <- read.table('LUAD_Editing.txt', header = TRUE, sep = '\t', stringsAsFactors = FALSE)
    lusc_edit <- read.table('LUSC_Editing.txt', header = TRUE, sep = '\t', stringsAsFactors = FALSE)
    merged_edit <- merge(luad_edit, lusc_edit, by = 'Editing', all = TRUE)

    # 只保留肿瘤样本（列名含 .01A）
    tumor_cols <- grep('\\.01A\$', colnames(merged_edit), value = TRUE)
    edit_tumor <- merged_edit[, c('Editing', tumor_cols)]

    # 计算每行有效样本比例（>0 且非 NA）
    valid_ratio <- rowSums(edit_tumor[, -1] > 0 & !is.na(edit_tumor[, -1])) / length(tumor_cols)
    selected <- edit_tumor[valid_ratio > 0.2, ]

    # 将 NA 替换为 0
    selected[, -1] <- lapply(selected[, -1], function(x) ifelse(is.na(x), 0, x))
    write.table(selected, file = 'Filled_lung_0.2.txt', sep = '\t', row.names = FALSE, quote = FALSE)
"

# ==================== 4. 注释编辑事件到基因 ====================
cd "$ProjectFold" || exit 1

echo "4. 注释编辑事件到基因名称（使用 GFF3）..."
Rscript -e "
    library(dplyr)
    library(rtracklayer)
    library(GenomicRanges)

    # 读取 GFF3 文件（基因注释）
    gff <- import('Expression/gencode.v32.chr_patch_hapl_scaff.annotation.gff3')
    gene_gr <- GRanges(
        seqnames = seqnames(gff)[gff\$type == 'gene'],
        ranges = ranges(gff)[gff\$type == 'gene'],
        strand = strand(gff)[gff\$type == 'gene'],
        gene_name = gff\$gene_name[gff\$type == 'gene']
    )

    # 读取编辑事件
    editing <- read.table('Editing/Filled_lung_0.2.txt', header = TRUE, sep = '\t', stringsAsFactors = FALSE)
    event_names <- editing\$Editing

    # 解析事件名（chr_pos_strand）
    parse_event <- function(x) {
        parts <- strsplit(x, '_', fixed = TRUE)
        chr <- sapply(parts, '[', 1)
        pos <- as.numeric(sapply(parts, '[', 2))
        strand <- sapply(parts, '[', 3)
        GRanges(seqnames = chr, ranges = IRanges(start = pos, end = pos), strand = strand)
    }
    events_gr <- parse_event(event_names)

    # 重叠查找
    hits <- findOverlaps(events_gr, gene_gr, ignore.strand = FALSE)
    gene_list <- split(gene_gr\$gene_name[subjectHits(hits)], queryHits(hits))
    mapped <- rep(NA_character_, length(event_names))
    mapped[as.integer(names(gene_list))] <- sapply(gene_list, function(gs) paste(unique(gs), collapse = ','))

    result <- data.frame(Editing = event_names, gene_name = mapped, stringsAsFactors = FALSE)
    write.table(result, file = 'Results/editing_events_unique_gene_annotation.txt',
                sep = '\t', row.names = FALSE, quote = FALSE)
"

echo "预处理完成！生成的文件位于："
echo "  - 表达谱（肿瘤样本）: $ProjectFold/Expression/Exp_tpm_01A_data.csv"
echo "  - 编辑事件（填充后）: $ProjectFold/Editing/Filled_lung_0.2.txt"
echo "  - 基因注释结果: $ProjectFold/Results/editing_events_unique_gene_annotation.txt"
echo "  - 临床数据: $ProjectFold/Expression/combined_clinical_01A.tsv"
echo "  - 生存数据: $ProjectFold/Expression/TCGA-lung.survival.tsv"