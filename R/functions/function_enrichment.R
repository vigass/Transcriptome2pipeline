# ==============================================================================
# 富集分析统一模块：GO / KEGG × ORA / GSEA（基于 clusterProfiler）
# 统一约定：
#   - 所有函数输入的基因 / 背景均为 gene symbol
#   - ORA  输入显著基因列表 (symbol)；GSEA 输入全基因排序向量 (名字为 symbol)
#   - 均返回 list(table = 结果数据框, object = clusterProfiler 对象)
# 说明：clusterProfiler / AnnotationDbi 会 mask dplyr 的 select/rename/filter，
#       为避免污染全局命名空间，本文件一律用命名空间限定符调用，不 attach 包。
# ==============================================================================

# -----------------------------------------------------------------------------
# 0. 内部工具：symbol -> Entrez 映射
# -----------------------------------------------------------------------------
symbol_to_entrez <- function(symbols) {
  symbols <- unique(trimws(as.character(symbols)))
  symbols <- symbols[!is.na(symbols) & symbols != ""]
  if (length(symbols) == 0) stop("输入的基因列表为空！")

  map <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                               keys    = symbols,
                               columns = "ENTREZID",
                               keytype = "SYMBOL")
  entrez <- unique(map$ENTREZID[!is.na(map$ENTREZID)])
  message(sprintf(">>> [ID转换] %d 个 symbol 映射到 %d 个唯一 Entrez ID",
                  length(symbols), length(entrez)))
  entrez
}

# 将 symbol 命名的排序向量转换为 Entrez 命名（GSEA-KEGG 用）
rank_to_entrez <- function(ranked) {
  map <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                               keys    = names(ranked),
                               columns = "ENTREZID",
                               keytype = "SYMBOL")
  map <- map[!is.na(map$ENTREZID) & !is.na(map$SYMBOL), , drop = FALSE]
  map$score <- ranked[map$SYMBOL]

  # 多个 symbol 映射到同一 Entrez 时，保留 |score| 最大的一个
  map <- map[order(abs(map$score), decreasing = TRUE), , drop = FALSE]
  map <- map[!duplicated(map$ENTREZID), , drop = FALSE]

  out <- stats::setNames(map$score, map$ENTREZID)
  out <- out[order(out, decreasing = TRUE)]
  message(sprintf(">>> [ID转换] %d 个基因保留 %d 个唯一 Entrez ID 用于 GSEA",
                  length(ranked), length(out)))
  out
}

# -----------------------------------------------------------------------------
# 1. 输入准备函数（从 DEG 结果 / 表达矩阵自动生成）
# -----------------------------------------------------------------------------

#' 提取显著差异基因列表（ORA 输入）
#' @param deg_table run_deg 返回的差异分析表（需含 symbol / Change 列）
#' @param direction "both" 合并 Up+Down；"up" 仅上调；"down" 仅下调
get_sig_genes <- function(deg_table, direction = c("both", "up", "down")) {
  direction <- match.arg(direction)
  keep <- switch(direction,
                 both = c("Up", "Down"),
                 up   = "Up",
                 down = "Down")
  deg_table$symbol[deg_table$Change %in% keep]
}

#' 生成全基因排序向量（GSEA 输入）
#' @param deg_table 差异分析表（需含 symbol / log2FC 列）
#' @return 命名数值向量，按 log2FC 降序，名字为 symbol
rank_genes <- function(deg_table) {
  ranked <- deg_table[!is.na(deg_table$log2FC), ]
  ranked <- ranked[!duplicated(ranked$symbol), ]
  out <- stats::setNames(ranked$log2FC, ranked$symbol)
  out[order(out, decreasing = TRUE)]
}

#' 提取背景基因集（ORA 的 universe）
#' @param exp_mat 表达矩阵（行名为 gene symbol）
get_universe <- function(exp_mat) {
  rownames(exp_mat)
}

# -----------------------------------------------------------------------------
# 2. ORA 富集（过表达分析）
# -----------------------------------------------------------------------------

#' GO 过表达分析（ORA）
#' @param genes 显著基因列表 (symbol)
#' @param universe 背景基因集 (symbol)，建议传入表达矩阵全部基因
#' @param ont GO 子本体："BP" / "MF" / "CC" / "ALL"
enrich_go <- function(genes, universe = NULL, ont = "BP",
                      pval_cutoff = 0.05, qval_cutoff = 0.2, ...) {
  res <- clusterProfiler::enrichGO(
    gene         = genes,
    universe     = universe,
    OrgDb        = org.Hs.eg.db::org.Hs.eg.db,
    keyType      = "SYMBOL",
    ont          = ont,
    pvalueCutoff = pval_cutoff,
    qvalueCutoff = qval_cutoff,
    readable     = TRUE,
    ...
  )
  list(table = as.data.frame(res), object = res)
}

#' KEGG 过表达分析（ORA，需联网）
#' @param genes 显著基因列表 (symbol)
#' @param universe 背景基因集 (symbol)
enrich_kegg <- function(genes, universe = NULL, organism = "hsa",
                        pval_cutoff = 0.05, qval_cutoff = 0.2, ...) {
  gene_entrez     <- symbol_to_entrez(genes)
  universe_entrez <- if (!is.null(universe)) symbol_to_entrez(universe) else NULL

  res <- tryCatch(
    clusterProfiler::enrichKEGG(
      gene         = gene_entrez,
      universe     = universe_entrez,
      organism     = organism,
      keyType      = "ncbi-geneid",
      pvalueCutoff = pval_cutoff,
      qvalueCutoff = qval_cutoff,
      ...
    ),
    error = function(e) stop("KEGG 富集失败（可能网络不可达）：",
                             conditionMessage(e), call. = FALSE)
  )
  res <- clusterProfiler::setReadable(res, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID")
  list(table = as.data.frame(res), object = res)
}

# -----------------------------------------------------------------------------
# 3. GSEA 富集（基因集富集分析）
# -----------------------------------------------------------------------------

#' GO 基因集富集分析（GSEA）
#' @param ranked 全基因排序向量（名字为 symbol，按 log2FC 降序）
gsea_go <- function(ranked, ont = "BP", pval_cutoff = 0.05, ...) {
  res <- clusterProfiler::gseGO(
    geneList     = ranked,
    OrgDb        = org.Hs.eg.db::org.Hs.eg.db,
    keyType      = "SYMBOL",
    ont          = ont,
    pvalueCutoff = pval_cutoff,
    ...
  )
  list(table = as.data.frame(res), object = res)
}

#' KEGG 基因集富集分析（GSEA，需联网）
#' @param ranked 全基因排序向量（名字为 symbol，按 log2FC 降序）
gsea_kegg <- function(ranked, organism = "hsa", pval_cutoff = 0.05, ...) {
  ranked_entrez <- rank_to_entrez(ranked)
  res <- tryCatch(
    clusterProfiler::gseKEGG(
      geneList     = ranked_entrez,
      organism     = organism,
      keyType      = "ncbi-geneid",
      pvalueCutoff = pval_cutoff,
      ...
    ),
    error = function(e) stop("KEGG GSEA 失败（可能网络不可达）：",
                             conditionMessage(e), call. = FALSE)
  )
  res <- clusterProfiler::setReadable(res, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID")
  list(table = as.data.frame(res), object = res)
}

# -----------------------------------------------------------------------------
# 4. 富集结果绘图
# -----------------------------------------------------------------------------

# 内部：把 enrichResult(ORA) / gseaResult(GSEA) 统一转成气泡图数据框
enrich_bubble_df <- function(enrich_obj, n = 10, ont = c("ALL", "BP", "CC", "MF")) {
  ont <- match.arg(ont)
  d <- as.data.frame(enrich_obj)
  if (nrow(d) == 0) stop("富集结果为空，无法绘图！")

  if (all(c("GeneRatio", "Count") %in% colnames(d))) {
    # ORA：GeneRatio 形如 "37/301"，拆成数值；Count 直接可用
    ratio <- vapply(strsplit(as.character(d$GeneRatio), "/"),
                    function(x) as.numeric(x[1]) / as.numeric(x[2]), numeric(1))
    d$GeneRatio_num <- ratio
    d$Count_num     <- d$Count
  } else if (all(c("setSize", "core_enrichment") %in% colnames(d))) {
    # GSEA：Count = 核心富集基因数，GeneRatio = Count / setSize
    core <- as.character(d$core_enrichment)
    cnt  <- ifelse(is.na(d$core_enrichment) | core == "",
                   0, lengths(strsplit(core, "/")))
    d$Count_num     <- cnt
    d$GeneRatio_num <- cnt / d$setSize
  } else {
    stop("无法识别的富集结果：缺少 GeneRatio/Count 或 setSize/core_enrichment 列")
  }

  # 按 ont 参数筛选本体（仅当结果含 ONTOLOGY 列，即 enrich_go(ont = "ALL") 时）
  has_ont <- "ONTOLOGY" %in% colnames(d)
  if (has_ont && ont != "ALL") {
    d <- d[d$ONTOLOGY == ont, , drop = FALSE]
    if (nrow(d) == 0) stop(sprintf("结果中没有 %s 的通路！", ont))
  } else if (!has_ont && ont != "ALL") {
    warning("富集结果不含 ONTOLOGY 列（仅跑了单一本体），已忽略 ont 参数。")
  }

  # 取 top n：ont = ALL 时按本体分组各取 top n，否则整体取 top n
  if (has_ont && ont == "ALL") {
    d$ONTOLOGY <- factor(d$ONTOLOGY, levels = c("BP", "CC", "MF"))
    parts <- lapply(split(d, d$ONTOLOGY), function(g) {
      g <- g[order(g$p.adjust, na.last = TRUE), , drop = FALSE]
      head(g, n)
    })
    d <- do.call(rbind, parts)
    # y 轴顺序：每个本体内最显著的在最上面
    lvls <- unlist(lapply(parts, function(g) rev(g$Description)), use.names = FALSE)
    d$Description <- factor(d$Description, levels = lvls)
  } else {
    d <- d[order(d$p.adjust, na.last = TRUE), , drop = FALSE]
    d <- head(d, n)
    d$Description <- factor(d$Description, levels = rev(d$Description))
  }
  d
}

#' 气泡图（ggplot 手绘版，ORA / GSEA 通用，Lancet 配色）
#' @param enrich_obj 富集对象
#' @param n 展示通路条数
#' @param title 图标题
#' @param ont 选择展示的 GO 本体："ALL" / "BP" / "CC" / "MF"
plot_dot <- function(enrich_obj, n = 10, title = "", ont = c("ALL", "BP", "CC", "MF")) {
  d <- enrich_bubble_df(enrich_obj, n = n, ont = ont)

  # Lancet 期刊配色：lancet[1] = 蓝 #00468B，lancet[2] = 红 #ED0000
  lancet <- ggsci::pal_lancet()(9)

  p <- ggplot2::ggplot(d, ggplot2::aes(x = GeneRatio_num, y = Description)) +
    ggplot2::geom_point(ggplot2::aes(size = Count_num, color = p.adjust)) +
    ggplot2::scale_color_gradientn(colors = lancet[c(2, 1)]) +
    ggplot2::scale_size_continuous(range = c(3, 8)) +
    ggplot2::labs(x = "GeneRatio", y = NULL,
                  size = "Count", color = "p.adjust", title = title)

  # GO 的 BP/CC/MF 分面（仅当展示多个本体时）
  if ("ONTOLOGY" %in% colnames(d) && length(unique(d$ONTOLOGY)) > 1) {
    p <- p + ggplot2::facet_wrap(~ ONTOLOGY, scales = "free_y")
  }

  if (exists("my_theme", inherits = TRUE)) p <- p + my_theme
  p
}

# 内部：按 ONTOLOGY 过滤富集对象（bar / cnet 用）
subset_enrich <- function(enrich_obj, ont) {
  res <- enrich_obj@result
  if (ont != "ALL" && "ONTOLOGY" %in% colnames(res)) {
    enrich_obj@result <- res[res$ONTOLOGY == ont, , drop = FALSE]
  }
  enrich_obj
}

#' 富集结果绘图（ORA / GSEA 通用）
#' @param enrich_obj enrich_go / enrich_kegg / gsea_go / gsea_kegg 返回的 object
#' @param p_type 绘图类型："dot" 气泡图 / "bar" 柱状图 / "cnet" 通路-基因网络图
#' @param n 展示的通路条数，默认 10
#' @param title 图标题
#' @param ont 选择展示的 GO 本体："ALL" / "BP" / "CC" / "MF"（需 enrich_go(ont = "ALL") 结果）
plot_enrich <- function(enrich_obj, p_type = c("dot", "bar", "cnet"), n = 10, title = "",
                        ont = c("ALL", "BP", "CC", "MF")) {
  # barplot 是 S3 方法（barplot.enrichResult），需先加载 enrichplot 命名空间才能正确分发
  requireNamespace("enrichplot", quietly = TRUE)
  p_type <- match.arg(p_type)
  ont    <- match.arg(ont)
  p <- switch(p_type,
              dot  = plot_dot(enrich_obj, n = n, title = title, ont = ont),
              bar  = {
                b <- graphics::barplot(subset_enrich(enrich_obj, ont), showCategory = n, title = title)
                if (exists("my_theme", inherits = TRUE)) b <- b + my_theme
                b
              },
              cnet = enrichplot::cnetplot(subset_enrich(enrich_obj, ont), showCategory = n))
  p
}
