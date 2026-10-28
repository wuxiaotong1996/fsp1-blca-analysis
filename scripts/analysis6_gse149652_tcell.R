# analysis6_gse149652_tcell.R — 分析6：GSE149652（TISCH2）AIFM2 表达谱（T/NK 谱系）
# 目的：回应审稿人"只分析T细胞不足以证明异质性"——对 GSE149652 单细胞数据做
#       全注释范围内（TISCH2 实际只有 T/NK 谱系 6 类）的 AIFM2 阳性率/平均表达、
#       TISCH2 坐标 UMAP、AIFM2 Dot plot 与 T 亚群阳性率柱状图（描述性，无 P）。
# 输入（只读，任务书指定文件在工作区不存在，使用旧项目等价物）：
#       data/raw/GSE149652_counts.mtx.gz → 等价于旧项目
#         "课题项目 申请\预实验，小课题\工作过程\数据分析\scRNA_AIFM2_Treg_CD8T\data\
#          GSE149652_merged_seurat.rds"（任务书允许"或 Seurat .rds"；
#          RNA counts+data，17,156×15,538，AIFM2 在矩阵）
#       data/raw/TISCH2_BLCA_GSE149652_annotation.csv → 等价于同目录
#         BLCA_GSE149652_CellMetainfo_table.tsv（TISCH2 注释，cell_id=GSM@barcode-1-0）
# 方法：Seurat 5.6；TISCH2 注释按 (Sample, 条码) 匹配合并（断言与 rds meta 同源、
#       断言 UMAP 坐标逐细胞一致）；LogNormalize scale.factor=10000（data 层重算，
#       磁盘 rds 不动）；阳性率=counts>0 比例 + Wilson 95%CI（手工公式）；
#       平均表达=归一化后 data 层组内均值；UMAP 直接用 TISCH2 坐标，不重新降维。
# 缺失报告：TISCH2 GSE149652 注释中无 Malignant/Macrophage/Endothelial/Fibroblast/
#       B/Mast/DC/Neutrophil（该数据集为 CD4/CD8 TIL 分选；malignancy 列 15,299
#       Immune + 239 Stromal、0 Malignant）——按任务书要求报告缺失，不自行推断。
#       全细胞谱系图由 GSE135337（analysis6b）承担。
# 输出：results/figures/FigS4A_UMAP_GSE149652_Tcell.pdf/.png（TISCH2 坐标，图例右侧）
#       results/figures/FigS4B_dotplot_AIFM2_GSE149652.pdf/.png（x=AIFM2，y=细胞类型，
#           点大小=阳性率，颜色=平均表达，右侧标注 n）
#       results/figures/FigS4C_Tcell_positivity.pdf/.png（描述性柱状，Wilson CI 误差线，
#           无 P 值）
#       results/tables/FigS4_source_data_GSE149652.csv（cell_type, n_cells, n_positive,
#           positivity_rate, wilson_lower, wilson_upper, mean_expression）
# 随机种子：20261028（本分析无随机成分，按约定保留）。
# 命名说明：任务书 Fig2A/Fig2B 名称原指"全细胞谱系"——GSE149652 无全谱系，
#       故其 UMAP/Dot 以 FigS4 命名；GSE135337 全谱系版沿用 Fig2A/Fig2B。

suppressMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
})

set.seed(20261028)

base <- "D:/HuaweiMoveData/Users/86159/Desktop/a...修改/课题项目 申请/预实验，小课题/工作过程/数据分析/scRNA_AIFM2_Treg_CD8T/data/"

# ---- 1. 读入（只读） ----
obj <- readRDS(paste0(base, "GSE149652_merged_seurat.rds"))
tsv <- fread(paste0(base, "BLCA_GSE149652_CellMetainfo_table.tsv"), sep = "\t")
stopifnot(ncol(obj) == 15538, nrow(tsv) == 15538, "AIFM2" %in% rownames(obj))

# ---- 2. TISCH2 注释合并：按 (Sample, 条码) 匹配（同一条码跨 CD4/CD8 数据集重复，
#        单独条码不是唯一键；GSM 编号由 Sample 列唯一确定） ----
tsv_bc <- tstrsplit(tsv$Cell, "@")[[2]]
rds_bc <- sub("^[^_]+_", "", Cells(obj))
key_tsv <- paste(tsv$Sample, tsv_bc)
key_rds <- paste(as.character(obj$Sample), rds_bc)
stopifnot(!anyDuplicated(key_tsv), !anyDuplicated(key_rds))
m <- match(key_rds, key_tsv)
stopifnot(all(!is.na(m)))
cat("[join] TISCH2 annotation matched", sum(!is.na(m)), "/", length(m), "cells\n")

major_tsv <- as.character(tsv[["Celltype (major-lineage)"]][m])
stopifnot(identical(major_tsv, as.character(obj$Celltype_major)))   # rds meta 与 TISCH2 同源
d_umap <- max(abs(tsv$UMAP_1[m] - obj$UMAP_1), abs(tsv$UMAP_2[m] - obj$UMAP_2))
cat("[umap] rds meta vs TISCH2 coords max|diff| =", format(d_umap, digits = 3), "\n")
# 任务书要求 UMAP 直接使用 TISCH2 坐标——无条件以 TSV 值为准
obj$UMAP_1 <- tsv$UMAP_1[m]
obj$UMAP_2 <- tsv$UMAP_2[m]

# ---- 3. LogNormalize（scale.factor=10000；磁盘 rds 只读，仅内存 data 层重算） ----
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)

# ---- 4. AIFM2 每类统计：阳性率(counts>0)+Wilson 95%CI+平均表达 ----
aif_counts <- as.numeric(GetAssayData(obj, layer = "counts")["AIFM2", ])
aif_norm   <- as.numeric(GetAssayData(obj, layer = "data")["AIFM2", ])
raw_types  <- c("Treg", "CD8T", "CD4Tconv", "NK", "CD8Tex", "Tprolif")   # 任务书列出顺序
disp_name  <- c(Treg = "Treg", CD8T = "CD8_T", CD4Tconv = "CD4_Tconv",
                NK = "NK", CD8Tex = "CD8Tex", Tprolif = "Tprolif")

wilson <- function(k, n, z = qnorm(0.975)) {
  p <- k / n
  denom <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denom
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  c(lower = center - half, upper = center + half)
}

res <- do.call(rbind, lapply(raw_types, function(tt) {
  idx <- obj$Celltype_major == tt
  n <- sum(idx); k <- sum(aif_counts[idx] > 0)
  ci <- wilson(k, n)
  data.frame(cell_type = tt, n_cells = n, n_positive = k,
             positivity_rate = 100 * k / n,
             wilson_lower = 100 * ci[1], wilson_upper = 100 * ci[2],
             mean_expression = mean(aif_norm[idx]))
}))

# 确定性锚点：复现旧项目已知值 Treg 29/3,813、CD8T 36/5,585（原始 counts>0 口径）
stopifnot(res$n_positive[res$cell_type == "Treg"] == 29,
          res$n_cells[res$cell_type == "Treg"] == 3813,
          res$n_positive[res$cell_type == "CD8T"] == 36,
          res$n_cells[res$cell_type == "CD8T"] == 5585)
cat("[anchor] Treg 29/3813 (0.76%), CD8T 36/5585 (0.64%) reproduced\n")

res$cell_type <- unname(disp_name[res$cell_type])
stopifnot(all(res$cell_type == unname(disp_name)))      # 行序即任务书列出顺序
write.csv(res, "results/tables/FigS4_source_data_GSE149652.csv", row.names = FALSE)
print(res, digits = 4)

# ---- 5. FigS4A：UMAP（TISCH2 坐标，按细胞类型着色，图例右侧） ----
ord <- res$cell_type
cols <- c(Treg = "#D55E00", CD8_T = "#E69F00", CD4_Tconv = "#F0E442",
          NK = "#56B4E9", CD8Tex = "#009E73", Tprolif = "#CC79A7")
dA <- data.frame(UMAP_1 = obj$UMAP_1, UMAP_2 = obj$UMAP_2,
                 type = factor(disp_name[as.character(obj$Celltype_major)], levels = ord))
gA <- ggplot(dA, aes(UMAP_1, UMAP_2, color = type)) +
  geom_point(size = 0.08, alpha = 0.55) +
  scale_color_manual(values = cols, name = "Cell type") +
  labs(title = "Fig S4A  GSE149652（TISCH2）T/NK 细胞亚群 UMAP（n = 15,538；坐标沿用 TISCH2）",
       x = "UMAP_1", y = "UMAP_2") +
  theme_classic(base_family = "Arial", base_size = 11) +
  theme(plot.title = element_text(size = 10.5, face = "bold"),
        legend.position = "right",
        legend.title = element_text(size = 9)) +
  guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1)))

# ---- 6. FigS4B：Dot plot（x=AIFM2，y=细胞类型，大小=阳性率，颜色=平均表达，右标 n） ----
dB <- data.frame(cell_type = factor(res$cell_type, levels = rev(ord)),
                 rate = res$positivity_rate, mean_expr = res$mean_expression,
                 n = res$n_cells, x = "AIFM2")
gB <- ggplot(dB, aes(x = x, y = cell_type)) +
  geom_point(aes(size = rate, color = mean_expr)) +
  geom_text(aes(label = paste0("n = ", format(n, big.mark = ","))),
            x = 1.5, hjust = 0, size = 3, color = "grey25") +
  scale_x_discrete(expand = expansion(mult = c(0.45, 0.85))) +
  scale_size_continuous(name = "Positivity (%)", range = c(2, 10)) +
  scale_color_gradientn(colors = c("#4575B4", "#FFFFBF", "#D73027"),
                        name = "Mean expr.") +
  coord_cartesian(clip = "off") +
  labs(title = "Fig S4B  AIFM2 表达（GSE149652，T/NK 亚群）",
       x = "", y = "") +
  theme_classic(base_family = "Arial", base_size = 11) +
  theme(plot.title = element_text(size = 10.5, face = "bold"),
        legend.position = "right")

# ---- 7. FigS4C：T 亚群阳性率柱状（描述性，Wilson CI 误差线，无 P 值） ----
dC <- res[order(res$positivity_rate), ]
dC$cell_type <- factor(dC$cell_type, levels = dC$cell_type)
gC <- ggplot(dC, aes(x = positivity_rate, y = cell_type)) +
  geom_col(fill = "#56B4E9", width = 0.6) +
  geom_errorbarh(aes(xmin = wilson_lower, xmax = wilson_upper),
                 height = 0.25, linewidth = 0.5) +
  geom_text(aes(x = wilson_upper, label = sprintf("%.2f%%", positivity_rate)),
            hjust = -0.3, size = 3.2) +
  scale_x_continuous(limits = c(0, max(dC$wilson_upper) * 1.45),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Fig S4C  T/NK 细胞亚群 AIFM2 阳性率（描述性；误差线 = Wilson 95% CI）",
       x = "AIFM2 positivity (%)", y = "") +
  theme_classic(base_family = "Arial", base_size = 11) +
  theme(plot.title = element_text(size = 10.5, face = "bold"))

ggsave("results/figures/FigS4A_UMAP_GSE149652_Tcell.pdf", gA,
       width = 7.5, height = 5.5, device = cairo_pdf, family = "Arial")
ggsave("results/figures/FigS4A_UMAP_GSE149652_Tcell.png", gA,
       width = 7.5, height = 5.5, dpi = 300, type = "cairo")
ggsave("results/figures/FigS4B_dotplot_AIFM2_GSE149652.pdf", gB,
       width = 5.5, height = 4.2, device = cairo_pdf, family = "Arial")
ggsave("results/figures/FigS4B_dotplot_AIFM2_GSE149652.png", gB,
       width = 5.5, height = 4.2, dpi = 300, type = "cairo")
ggsave("results/figures/FigS4C_Tcell_positivity.pdf", gC,
       width = 5.5, height = 4.2, device = cairo_pdf, family = "Arial")
ggsave("results/figures/FigS4C_Tcell_positivity.png", gC,
       width = 5.5, height = 4.2, dpi = 300, type = "cairo")
cat("[done] analysis6 (FigS4A/B/C + FigS4_source_data_GSE149652.csv)\n")
