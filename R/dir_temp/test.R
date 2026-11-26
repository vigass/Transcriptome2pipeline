rm(list = ls())

path_origin <- "~/Desktop/Project/Transcriptome2pipeline/R/"
path_origin <- path.expand(path_origin)   # 展开 ~，避免部分文件函数不识别

setwd(file.path(path_origin, "dir_temp"))

library(data.table)
library(tidyverse)

# 1. 读入数据 ---------------------------------------------------------------
exp         <- fread(file.path(path_origin, "datas", "exp_count.csv"))
sample_info <- fread(file.path(path_origin, "datas", "sample_info.csv"))

# 2. 载入全部自定义函数 -----------------------------------------------------
func_dir <- file.path(path_origin, "functions")
for (f in list.files(func_dir, full.names = TRUE, pattern = "\\.R$")) {
  source(f)
}

# 3. 构建纯数值 count 矩阵（行名=基因 symbol，列名=样本 ID）------------------
count_mat <- as.matrix(exp[, -1, with = FALSE])
rownames(count_mat) <- exp$symbol
storage.mode(count_mat) <- "integer"

# 4. 标准化：log2(CPM + 1)，用于 PCA / 热图 ---------------------------------
log2cpm_mat <- convert_counts(count_mat, target = "log2cpm")

# 5. PCA --------------------------------------------------------------------
p_pca <- plot_pca(log2cpm_mat, sample_info,
                  group_col = "Groups", case = "sepsis", control = "healthy",
                  title = "PCA Plot (sepsis vs healthy)")
print(p_pca)

# 6. 差异分析（样本量 392，自动走 limma-voom）-------------------------------
deg_res <- run_deg(count_mat, sample_info,
                   group_col = "Groups", case = "sepsis", control = "healthy",
                   covariates = NULL,
                   logfc_cutoff = 1.5, pval_cutoff = 0.05,
                   pval_type = "padj", method = "auto")
deg_res$summary

# 7. 差异结果绘图：火山图 + 热图 --------------------------------------------
p_volcano <- plot_volcano(deg_res$table,
                          pval_type = "padj", pval_cutoff = 0.05,
                          logfc_cutoff = 1.5, top_n = 10,
                          title = "Volcano Plot (sepsis vs healthy)")
print(p_volcano)

ht <- plot_heatmap(log2cpm_mat, deg_res$table, sample_info,
                   group_col = "Groups", case = "sepsis", control = "healthy",
                   top_n = 25, title = "Top DEGs Heatmap (sepsis vs healthy)")
ComplexHeatmap::draw(ht)

# 8. 富集分析（GO / KEGG × ORA / GSEA）--------------------------------------
sig_genes <- get_sig_genes(deg_res$table, direction = "both")   # Up+Down 合并
universe  <- get_universe(count_mat)                             # 全部基因作背景
ranked    <- rank_genes(deg_res$table)                           # 全基因按 log2FC 降序

go_ora        <- enrich_go(sig_genes, universe, ont = "ALL")  # GO 过表达（离线，BP/CC/MF 分面）
gsea_go_res   <- gsea_go(ranked)                   # GO GSEA（离线）
kegg_ora      <- enrich_kegg(sig_genes, universe)  # KEGG 过表达（需联网）
gsea_kegg_res <- gsea_kegg(ranked)                 # KEGG GSEA（需联网）

print(plot_enrich(go_ora$object, p_type = "dot", n = 5, ont = "ALL",
                  title = "GO Over-Representation (BP / CC / MF)"))
print(plot_enrich(gsea_go_res$object, p_type = "dot", n = 5,
                  title = "GO GSEA (BP)"))
