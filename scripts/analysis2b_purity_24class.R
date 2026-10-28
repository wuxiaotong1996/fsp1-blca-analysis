# analysis2b_purity_24class.R — 纯度校正重评估：官方 ESTIMATE 包 + 24 类 ssGSEA 偏相关
# 目的：回应审稿人"纯度混杂"质疑——用官方 ESTIMATE R 包（v1.0.11，SourceForge
#       estimateproject 源码目录安装）计算 407 例肿瘤的 Stromal/Immune/ESTIMATE 评分，
#       按 TumorPurity = cos(0.6049872018 + 0.0001467884*ESTIMATEScore) 换算纯度；
#       对 FSP1 与 24 类免疫细胞（Bindea 2013）ssGSEA 评分计算原始与纯度校正偏 ρ；
#       并报告 FSP1 与 ImmuneScore/StromalScore 的相关。
# 输入说明：任务书指定 data/processed/TCGA_BLCA_expr_clean.csv——该文件当前磁盘版本
#       为 2 列最简版（sample, AIFM2，426 行含 19 正常），无法支撑 ESTIMATE/ssGSEA
#       （需要全基因矩阵），故本脚本按分析2/4同法从 raw HiSeqV2.gz 读全矩阵取 407
#       肿瘤（只读）；AIFM2 表达与 expr_clean 完全同源。若用户提供全矩阵版
#       expr_clean，替换读取段即可。
# 输出：results/tables/Analysis2b_estimate_scores_pkg.csv（包计算 4 项评分）
#       results/tables/Analysis2b_fsp1_estimate_corr.csv（FSP1×Immune/Stromal/ESTIMATE）
#       results/tables/TableS1_ssGSEA_full.csv（cell_type, raw_rho, raw_p,
#           partial_rho, partial_p, n——本任务规格版；分析4 详细版已存
#           TableS1_ssGSEA_full_detail.csv）
#       results/figures/Fig1C_partial_correlation.pdf/.png（300 dpi，Arial，cairo）
# 随机种子：20261028。交叉校验：包计算的评分 vs 分析2 手工移植（Analysis2_estimate_scores.csv）
#       须一致（同一官方算法），防实现漂移。

suppressMessages({
  library(data.table)
  library(estimate)
  library(GSVA)
  library(ppcor)
  library(ggplot2)
})

set.seed(20261028)

# ---- 1. 全基因矩阵（407 肿瘤，raw 只读） ----
X <- as.data.frame(fread("data/raw/TCGA.BLCA.HiSeqV2.gz", nThread = 8))
rn <- toupper(as.character(X[[1]])); X[[1]] <- NULL
is_tum <- vapply(colnames(X), function(s) {
  p <- strsplit(s, "-", fixed = TRUE)[[1]]
  length(p) >= 4 && substr(p[4], 1, 2) == "01"
}, logical(1))
X <- X[, is_tum]; X <- as.matrix(X); rownames(X) <- rn
if (anyDuplicated(rn)) {
  g <- factor(rn)
  X <- rowsum(X, g) / as.numeric(table(g))
}
stopifnot(ncol(X) == 407, "AIFM2" %in% rownames(X))
f <- as.numeric(X["AIFM2", ])
cat("[input] tumor", ncol(X), "genes", nrow(X), "\n")

# ---- 2. 官方 ESTIMATE 包（GCT 输入，platform="affymetrix" 出纯度） ----
gct_in <- "data/processed/TCGA_BLCA_tumor.gct"
gct_out <- "data/processed/TCGA_BLCA_estimate_out.gct"
con <- file(gct_in, "w")
writeLines(c("#1.2", sprintf("%d\t%d", nrow(X), ncol(X)),
             paste(c("Name", "Description", colnames(X)), collapse = "\t")), con)
write.table(cbind(rownames(X), rownames(X), X), con, sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
close(con)

estimateScore(gct_in, gct_out, platform = "affymetrix")
# 输出 GCT：第1列 NAME、第2列 Description、其后 407 个样本列（4 行评分）。
# 包内 readGCT 用 check.names=TRUE 读入，条码 "-" 被换成 "."（TCGA-DK-AA6S-01 →
# TCGA.DK.AA6S.01），样本顺序不变——按同样变换比对后回映。
out <- read.delim(gct_out, header = TRUE, sep = "\t", skip = 2, check.names = FALSE)
stopifnot(nrow(out) == 4, ncol(out) == 2 + 407,
          identical(colnames(out)[-(1:2)], gsub("-", ".", colnames(X), fixed = TRUE)))
score_names <- as.character(out[[1]])
est <- as.data.frame(t(as.matrix(out[, -(1:2), drop = FALSE])))
colnames(est) <- score_names
stopifnot(nrow(est) == 407,
          setequal(colnames(est), c("StromalScore", "ImmuneScore",
                                    "ESTIMATEScore", "TumorPurity")))
purity <- as.numeric(est$TumorPurity)
purity[purity < 0] <- NA
est$TumorPurity <- purity
write.csv(cbind(sample = colnames(X), est),
          "results/tables/Analysis2b_estimate_scores_pkg.csv", row.names = FALSE)

# 交叉校验：包 vs 分析2 手工移植（同一官方算法，应逐样本一致）
old <- read.csv("results/tables/Analysis2_estimate_scores.csv", stringsAsFactors = FALSE)
m_old <- match(colnames(X), old$sample)
d_ss <- max(abs(est$StromalScore - old$StromalScore[m_old]))
d_is <- max(abs(est$ImmuneScore - old$ImmuneScore[m_old]))
d_es <- max(abs(est$ESTIMATEScore - old$ESTIMATEScore[m_old]))
d_pu <- max(abs(purity - old$TumorPurity[m_old]), na.rm = TRUE)
stopifnot(d_ss < 1e-6, d_is < 1e-6, d_es < 1e-6, d_pu < 1e-6)
cat(sprintf("[check] package vs port: max|diff| Stromal=%.2e Immune=%.2e ESTIMATE=%.2e Purity=%.2e\n",
            d_ss, d_is, d_es, d_pu))
cat(sprintf("[purity] usable %d/407  mean %.3f median %.3f\n",
            sum(!is.na(purity)), mean(purity, na.rm = TRUE),
            median(purity, na.rm = TRUE)))

# ---- 3. FSP1 与 ImmuneScore / StromalScore 相关 ----
ok <- !is.na(purity)
sc_cor <- do.call(rbind, lapply(c("ImmuneScore", "StromalScore", "ESTIMATEScore"),
  function(nm) {
    y <- as.numeric(est[[nm]])
    r1 <- cor.test(f, y, method = "spearman", exact = FALSE)
    r2 <- pcor.test(f[ok], y[ok], purity[ok], method = "spearman")
    data.frame(score = nm, rho_raw = as.numeric(r1$estimate), p_raw = r1$p.value,
               rho_partial = as.numeric(r2$estimate), p_partial = r2$p.value)
  }))
write.csv(sc_cor, "results/tables/Analysis2b_fsp1_estimate_corr.csv", row.names = FALSE)
print(sc_cor, digits = 4)

# ---- 4. 24 类 ssGSEA（Bindea 2013，GSVA fastRndWalk，与分析4同式） ----
sig <- read.csv("data/processed/bindea24_signatures.csv", stringsAsFactors = FALSE)
sets <- lapply(split(sig$Gene, sig$CellType), toupper)
sets <- sets[order(names(sets))]
stopifnot(length(sets) == 24)
ss <- t(gsva(ssgseaParam(X, sets, normalize = FALSE), verbose = FALSE))

res <- do.call(rbind, lapply(names(sets), function(nm) {
  y <- as.numeric(ss[, nm])
  r1 <- cor.test(f, y, method = "spearman", exact = FALSE)
  r2 <- pcor.test(f[ok], y[ok], purity[ok], method = "spearman")
  data.frame(cell_type = nm,
             raw_rho = as.numeric(r1$estimate), raw_p = r1$p.value,
             partial_rho = as.numeric(r2$estimate), partial_p = r2$p.value,
             n = sum(ok))
}))
res <- res[order(res$partial_rho, decreasing = TRUE), ]
write.csv(res, "results/tables/TableS1_ssGSEA_full.csv", row.names = FALSE)
print(res, digits = 4)

# ---- 5. 终端报告：原始分析 CD8⁺T（22 类教科书集 ρ=0.099/P=0.046）校正后 ----
imm <- read.csv("data/processed/immune_scores.csv",
                fileEncoding = "UTF-8-BOM", check.names = FALSE)
smp <- as.character(imm[[1]]); stopifnot(all(smp %in% colnames(X)))
fi <- as.numeric(f[match(smp, colnames(X))])
y8 <- as.numeric(imm[["CD8+ T 细胞"]])
r8_raw <- cor.test(fi, y8, method = "spearman", exact = FALSE)
r8_par <- pcor.test(fi[ok], y8[ok], purity[ok], method = "spearman")
cat(sprintf("[CD8 original] rho=%.3f P=%.3f  ->  purity-adjusted rho=%.3f P=%.3f\n",
            r8_raw$estimate, r8_raw$p.value, r8_par$estimate, r8_par$p.value))
cat(sprintf("[CD8 Bindea24] raw rho=%.3f P=%.4f -> partial rho=%.3f P=%.4f\n",
            res$raw_rho[res$cell_type == "CD8 T cells"],
            res$raw_p[res$cell_type == "CD8 T cells"],
            res$partial_rho[res$cell_type == "CD8 T cells"],
            res$partial_p[res$cell_type == "CD8 T cells"]))

# ---- 6. Fig1C：横向条形图（左 raw ρ / 右 partial ρ，按 partial ρ 排序） ----
dplot <- rbind(
  data.frame(cell_type = res$cell_type, type = "Raw ρ", rho = res$raw_rho, p = res$raw_p),
  data.frame(cell_type = res$cell_type, type = "Partial ρ", rho = res$partial_rho, p = res$partial_p))
dplot$type <- factor(dplot$type, levels = c("Raw ρ", "Partial ρ"))
dplot$cell_type <- factor(dplot$cell_type, levels = rev(res$cell_type))
dplot$sig <- factor(ifelse(dplot$p < 0.05, "P < 0.05", "n.s."),
                    levels = c("P < 0.05", "n.s."))

g <- ggplot(dplot, aes(x = rho, y = cell_type, fill = sig)) +
  geom_col(width = 0.68) +
  geom_vline(xintercept = 0, color = "grey30", linewidth = 0.3) +
  scale_fill_manual(values = c("P < 0.05" = "#D7191C", "n.s." = "#B0B0B0"),
                    name = "") +
  facet_grid(. ~ type) +
  labs(title = "Fig 1C  FSP1 × 24 类免疫细胞 ssGSEA：原始 ρ vs 肿瘤纯度校正偏 ρ（n = 342）",
       x = "Spearman ρ", y = "") +
  theme_classic(base_family = "Arial", base_size = 11) +
  theme(plot.title = element_text(size = 10.5, face = "bold"),
        legend.position = "top",
        strip.background = element_rect(fill = "grey92", color = NA),
        axis.text.y = element_text(size = 8.5))

ggsave("results/figures/Fig1C_partial_correlation.pdf", g,
       width = 9, height = 6, device = cairo_pdf, family = "Arial")
ggsave("results/figures/Fig1C_partial_correlation.png", g,
       width = 9, height = 6, dpi = 300, type = "cairo")
cat("[done] analysis2b\n")
