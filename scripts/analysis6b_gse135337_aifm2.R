# analysis6b_gse135337_aifm2.R — 分析6b：GSE135337 全细胞谱系 AIFM2 表达（Fig2B 主图）
# 目的：GSE149652（TISCH2）为 CD4/CD8 TIL 分选数据，无 Malignant/Macrophage 等
#       非 T 谱系（分析6已报告缺失）。审稿人"只分析T细胞不足以证明异质性"的
#       全谱系回应由 GSE135337（8 例未治疗 BLCA，Chen et al. Nat Commun 2021）
#       承担：按分析3+5 的谱系归类（AddModuleScore argmax），计算 AIFM2 各谱系
#       阳性率/Wilson 95%CI/平均表达，出 Fig2B Dot plot（x=AIFM2，大小=阳性率，
#       颜色=平均表达，右标 n）与 Fig2_source_data.csv。
# 输入（只读）：data/processed/GSE135337_seurat.rds（21,854×42,455；meta$lineage：
#       Tumor 34,034 / Macrophage 1,652 / Endothelial 659 / Fibroblast 3,387 /
#       T 759 / Other 1,964；counts/data 按 8 样本分层 counts.BCi/data.BCi）
# 方法：按细胞名对位（顺序无关）拼接 8 层 AIFM2 行 → 阳性率=counts>0 比例，
#       Wilson 95%CI 手工公式；平均表达=归一化 data 层组内均值。
# 输出：results/tables/Fig2_source_data.csv（cell_type, n_cells, n_positive,
#           positivity_rate, wilson_lower, wilson_upper, mean_expression）
#       results/figures/Fig2B_dotplot_AIFM2.pdf/.png（300 dpi，cairo+Arial）
# 注：Fig2A_UMAP_all_celltypes（同数据全谱系 UMAP）由分析3+5 已出，本脚本不覆盖。

suppressMessages({
  library(Seurat)
  library(ggplot2)
})

set.seed(20261028)

# ---- 1. 读入（只读） ----
obj <- readRDS("data/processed/GSE135337_seurat.rds")
stopifnot(ncol(obj) == 42455, "AIFM2" %in% rownames(obj))
lin <- as.character(obj$lineage)
names(lin) <- Cells(obj)
lin_ord <- c("Tumor", "Macrophage", "T", "Endothelial", "Fibroblast", "Other")
stopifnot(setequal(unique(lin), lin_ord),
          sum(lin == "Tumor") == 34034, sum(lin == "Macrophage") == 1652,
          sum(lin == "Endothelial") == 659, sum(lin == "Fibroblast") == 3387,
          sum(lin == "T") == 759, sum(lin == "Other") == 1964)
cat("[input] GSE135337 cells:", length(lin), " lineages:", paste(lin_ord, collapse = "/"), "\n")

# ---- 2. 拼接 8 个分层（按细胞名对位，顺序无关） ----
samp <- c("BC1", "BC2", "BC3", "BC4", "BC5", "BC6", "BC7", "BCN")
aif_c <- setNames(numeric(ncol(obj)), Cells(obj))
aif_d <- setNames(numeric(ncol(obj)), Cells(obj))
for (s in samp) {
  Lc <- paste0("counts.", s); Ld <- paste0("data.", s)
  cn <- colnames(GetAssayData(obj, layer = Lc))
  stopifnot(all(cn %in% Cells(obj)))
  aif_c[cn] <- as.numeric(GetAssayData(obj, layer = Lc)["AIFM2", , drop = FALSE])
  aif_d[cn] <- as.numeric(GetAssayData(obj, layer = Ld)["AIFM2", , drop = FALSE])
}
stopifnot(all(!is.na(aif_c)), all(!is.na(aif_d)))
cat("[layers] 8 layers concatenated, no missing cells\n")

# ---- 3. 每谱系统计 ----
wilson <- function(k, n, z = qnorm(0.975)) {
  p <- k / n
  denom <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denom
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  c(lower = center - half, upper = center + half)
}

res <- do.call(rbind, lapply(lin_ord, function(tt) {
  idx <- lin == tt
  n <- sum(idx); k <- sum(aif_c[idx] > 0)
  ci <- wilson(k, n)
  data.frame(cell_type = tt, n_cells = n, n_positive = k,
             positivity_rate = 100 * k / n,
             wilson_lower = 100 * ci[1], wilson_upper = 100 * ci[2],
             mean_expression = mean(aif_d[idx]))
}))
write.csv(res, "results/tables/Fig2_source_data.csv", row.names = FALSE)
print(res, digits = 4)
# 确定性锚点：n 与阳性计数（Tumor 2,901/34,034=8.524%；Macrophage 90/1,652=5.448%）
stopifnot(res$n_positive[res$cell_type == "Tumor"] == 2901,
          res$n_positive[res$cell_type == "Macrophage"] == 90)

# ---- 4. Fig2B：Dot plot（x=AIFM2，大小=阳性率，颜色=平均表达，右标 n） ----
dB <- data.frame(cell_type = factor(res$cell_type, levels = rev(lin_ord)),
                 rate = res$positivity_rate, mean_expr = res$mean_expression,
                 n = res$n_cells, x = "AIFM2")
gB <- ggplot(dB, aes(x = x, y = cell_type)) +
  geom_point(aes(size = rate, color = mean_expr)) +
  geom_text(aes(label = paste0("n = ", format(n, big.mark = ","))),
            x = 1.5, hjust = 0, size = 3, color = "grey25") +
  scale_x_discrete(expand = expansion(mult = c(0.45, 0.85))) +
  scale_size_continuous(name = "Positivity (%)", range = c(2, 12)) +
  scale_color_gradientn(colors = c("#4575B4", "#FFFFBF", "#D73027"),
                        name = "Mean expr.") +
  coord_cartesian(clip = "off") +
  labs(title = "Fig 2B  AIFM2 全细胞谱系表达（GSE135337，8 例未治疗 BLCA）",
       x = "", y = "") +
  theme_classic(base_family = "Arial", base_size = 11) +
  theme(plot.title = element_text(size = 10.5, face = "bold"),
        legend.position = "right")

ggsave("results/figures/Fig2B_dotplot_AIFM2.pdf", gB,
       width = 5.5, height = 4.2, device = cairo_pdf, family = "Arial")
ggsave("results/figures/Fig2B_dotplot_AIFM2.png", gB,
       width = 5.5, height = 4.2, dpi = 300, type = "cairo")
cat("[done] analysis6b (Fig2B_dotplot_AIFM2 + Fig2_source_data.csv)\n")
