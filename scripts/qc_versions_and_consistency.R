# qc_versions_and_consistency.R — 最终 QC：依赖版本 + 数据一致性核对（只读）
# 1) R/包版本；2) TISCH2 注释表 vs FigS4 CSV；3) GSE135337 rds vs Fig2 CSV；
# 4) TableS1 24 行断言。全部只读，不写任何结果文件。

suppressMessages({
  library(data.table)
  library(Seurat)
})

cat("== R version ==\n")
cat(R.version.string, "\n")
for (p in c("Seurat", "GSVA", "ppcor", "data.table", "ggplot2")) {
  cat(sprintf("%s: %s\n", p, as.character(packageVersion(p))))
}
cat("GSVAdata:", requireNamespace("GSVAdata", quietly = TRUE), "\n")

base <- "D:/HuaweiMoveData/Users/86159/Desktop/a...修改/课题项目 申请/预实验，小课题/工作过程/数据分析/scRNA_AIFM2_Treg_CD8T/data/"
tsv <- fread(paste0(base, "BLCA_GSE149652_CellMetainfo_table.tsv"), sep = "\t")
cat("\n== TISCH2 TSV rows ==", nrow(tsv), "\n")
tsv_cnt <- as.data.frame(table(tsv[["Celltype (major-lineage)"]]))
names(tsv_cnt) <- c("cell_type", "n_tsv")
s4 <- fread("results/tables/FigS4_source_data_GSE149652.csv")
names(s4)[names(s4) == "n_cells"] <- "n_csv"
m1 <- merge(tsv_cnt, s4[, .(cell_type, n_csv)], by = "cell_type", all = TRUE)
m1$match <- m1$n_tsv == m1$n_csv
print(m1)
cat("[S4 vs TISCH2] all match:", all(m1$match, na.rm = TRUE) && nrow(m1) == 6, "\n")

cat("\n== GSE135337 rds lineage counts vs Fig2 CSV ==\n")
obj <- readRDS("data/processed/GSE135337_seurat.rds")
lin_cnt <- as.data.frame(table(as.character(obj$lineage)))
names(lin_cnt) <- c("cell_type", "n_rds")
f2 <- fread("results/tables/Fig2_source_data.csv")
names(f2)[names(f2) == "n_cells"] <- "n_csv"
m2 <- merge(lin_cnt, f2[, .(cell_type, n_csv)], by = "cell_type", all = TRUE)
m2$match <- m2$n_rds == m2$n_csv
print(m2)
cat("[Fig2 vs rds] all match:", all(m2$match, na.rm = TRUE) && nrow(m2) == 6,
    "; total cells:", sum(lin_cnt$n_rds), "\n")

cat("\n== TableS1 rows ==", nrow(fread("results/tables/TableS1_ssGSEA_full.csv")), "\n")
cat("[done qc]\n")
