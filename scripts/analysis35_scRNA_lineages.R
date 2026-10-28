# analysis35_scRNA_lineages.R — 分析3+5：scRNA 全细胞谱系 UMAP + Dot plot（图2重绘）
# 数据：GSE135337（Chen et al. Nat Commun 2021，8 例未治疗 BLCA）原始计数表，
#       来自 GEO 附 GSE135337_RAW.tar（data/raw/GSE135337/，只读）。
# 谱系定义（按任务书）：肿瘤上皮（EPCAM+）、巨噬（CD68+）、内皮、成纤维、T 细胞，
#       其余细胞归 "Other"——用 Seurat AddModuleScore 打分 + argmax 归类。
# 约束：不做以单细胞为观测单位的统计检验（伪重复）→ 本分析只有可视化与描述性计数。
# 随机种子：20261028（RunUMAP 前 set.seed）。
# 输出：results/figures/Fig2A_UMAP_all_celltypes.pdf/.png（300 dpi，Arial）
#       results/figures/Fig2B_dotplot_lineage_markers.pdf/.png
#       results/tables/Fig2_lineage_counts.csv（每样本×谱系细胞数，描述性）
#       data/processed/GSE135337_seurat.rds（QC+注释后对象，供后续分析复用）

suppressMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(ggplot2)
})

set.seed(20261028)

# ---- 1. 读取 8 个样本计数表（只读） ----
files <- list.files("data/raw/GSE135337",
                    pattern = "_gene_cell_exprs_table\\.(txt|xls)$",
                    full.names = TRUE)
stopifnot(length(files) == 8)
samples <- sub(".*_(BC[0-9]+|BCN)_gene_cell_exprs_table\\..*", "\\1", basename(files))
cat("samples:", paste(samples, collapse = ", "), "\n")

seu_list <- list()
for (i in seq_along(files)) {
  s <- samples[i]
  cat("[read]", s, "...\n")
  dt <- fread(files[i], nThread = 8)
  g <- dt$Gene_ID; sym <- dt$Symbol; dt[, c("Gene_ID", "Symbol") := NULL]
  m <- as.matrix(dt)                       # 整数计数矩阵（每样本最大 ~580MB，可承受）
  sm <- Matrix::Matrix(m, sparse = TRUE)
  rm(dt, m); gc()
  # 行名用基因符号（marker 注释按符号进行；重复符号加 .1 后缀）
  rownames(sm) <- make.unique(sym)
  seu_list[[s]] <- CreateSeuratObject(counts = sm, project = s,
                                      min.cells = 3, min.features = 200)
  cat("  cells after min.features>=200:", ncol(seu_list[[s]]), "\n")
  rm(sm); gc()
}

# ---- 2. 合并（基因取并集，细胞名加样本前缀防条码冲突） ----
seu <- merge(seu_list[[1]], y = seu_list[2:8], add.cell.ids = samples)
rm(seu_list); gc()
cat("merged cells:", ncol(seu), "genes:", nrow(seu), "\n")

# ---- 3. 轻量 QC ----
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu <- subset(seu, subset = nFeature_RNA >= 200 & nFeature_RNA <= 6000 &
                percent.mt < 25)
cat("cells after QC:", ncol(seu), "\n")

# ---- 4. 标准流程 ----
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, nfeatures = 2000, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:20, verbose = FALSE)   # 已 set.seed(20261028)

# ---- 5. 谱系注释（marker 打分 + argmax） ----
marker_sets <- list(
  Tumor       = c("EPCAM", "KRT5", "KRT7", "KRT8", "KRT18", "KRT19",
                  "UPK1A", "UPK1B", "UPK2", "UPK3A", "UPK3B"),
  Macrophage  = c("CD68", "CD163", "CD14", "CSF1R", "ITGAM", "MRC1",
                  "FCGR3A", "LYZ", "CD86"),
  Endothelial = c("PECAM1", "VWF", "CDH5", "FLT1", "KDR", "ENG", "CLDN5"),
  Fibroblast  = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA",
                  "PDGFRB", "ACTA2", "FAP", "COL6A1"),
  T           = c("CD3D", "CD3E", "CD3G", "CD2", "TRAC", "TRBC1", "IL7R",
                  "CD4", "CD8A")
)
# 核心 marker 必须在矩阵中（防文件间基因集差异）
core <- c(EPCAM = "Tumor", CD68 = "Macrophage", PECAM1 = "Endothelial",
          COL1A1 = "Fibroblast", CD3D = "T")
for (g in names(core))
  stopifnot(g %in% rownames(seu), core[g] %in% names(marker_sets))

seu <- AddModuleScore(seu, features = marker_sets, name = "LineageScore")
scores <- as.matrix(seu@meta.data[, paste0("LineageScore", 1:5)])
best <- apply(scores, 1, which.max); bestv <- apply(scores, 1, max)
lab <- c("Tumor", "Macrophage", "Endothelial", "Fibroblast", "T")
lineage <- ifelse(bestv < 0, "Other", lab[best])
seu$lineage <- factor(lineage,
                      levels = c("Tumor", "Macrophage", "Endothelial",
                                 "Fibroblast", "T", "Other"))
print(table(seu$lineage))

# ---- 6. 保存对象 + 描述性计数表 ----
saveRDS(seu, "data/processed/GSE135337_seurat.rds")
cnt <- as.data.frame.matrix(table(seu$orig.ident, seu$lineage))
cnt <- cbind(sample = rownames(cnt), cnt)
write.csv(cnt, "results/tables/Fig2_lineage_counts.csv", row.names = FALSE)
cat("[done] rds + lineage counts csv\n")

# ---- 7. 图（cairo 设备 + Arial；Windows 版 R 的 windows 设备不认 Arial） ----
pal <- c(Tumor = "#E69F00", Macrophage = "#56B4E9", Endothelial = "#009E73",
         Fibroblast = "#F0E442", T = "#0072B2", Other = "#999999")

Idents(seu) <- "lineage"
gA <- DimPlot(seu, reduction = "umap", cols = pal, pt.size = 0.15) +
  ggtitle(sprintf("Fig 2A  全细胞谱系 UMAP（GSE135337，8 例未治疗膀胱癌，n = %d 细胞）",
                  ncol(seu))) +
  theme_classic(base_family = "Arial", base_size = 11) +
  theme(plot.title = element_text(size = 11, face = "bold"))

dot_feats <- c("EPCAM", "KRT19", "CD68", "CD163", "PECAM1", "VWF",
               "COL1A1", "DCN", "CD3D", "CD3E")
gB <- DotPlot(seu, features = dot_feats) +
  scale_colour_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B") +
  labs(title = "Fig 2B  谱系 marker 表达（点大小 = 表达细胞比例，颜色 = 平均表达）") +
  theme_classic(base_family = "Arial", base_size = 11) +
  theme(plot.title = element_text(size = 11, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("results/figures/Fig2A_UMAP_all_celltypes.pdf", gA,
       width = 7, height = 6, device = cairo_pdf, family = "Arial")
ggsave("results/figures/Fig2A_UMAP_all_celltypes.png", gA,
       width = 7, height = 6, dpi = 300, type = "cairo")
ggsave("results/figures/Fig2B_dotplot_lineage_markers.pdf", gB,
       width = 7, height = 4.5, device = cairo_pdf, family = "Arial")
ggsave("results/figures/Fig2B_dotplot_lineage_markers.png", gB,
       width = 7, height = 4.5, dpi = 300, type = "cairo")
cat("[done] Fig2A_UMAP_all_celltypes + Fig2B_dotplot_lineage_markers\n")
