# analysis1b_downsampling.R — 分析1 v2：TCGA-BLCA 降采样（回应审稿人"正常 n=19 太少"）
# 输入：data/processed/TCGA_BLCA_expr_clean.csv（含 AIFM2 列）、
#       data/processed/TCGA_BLCA_group.csv（tumor/normal 标签）——均只读。
# 方法：从 407 肿瘤中无放回抽 19 例，与 19 例正常行 Mann-Whitney U（exact=FALSE，
#       表达值含并列），记录 P 值与 log2FC=mean(肿瘤)−mean(正常)；重复 1000 次；
#       set.seed(20261028)。另计算全样本原始 P（作为分位数参照）。
# 输出：results/tables/downsampling_results.csv（run, p_value, log2FC，1000 行）
#       results/figures/FigS1_downsampling.pdf/.png（300 dpi，Arial，cairo；
#       左=1000 次 P 直方图标注中位 P 与 P<0.05 比例；右=log2FC 密度图标注均值/中位数）
# 注：本版覆盖分析1旧版 FigS1（旧 summary CSV 保留）；输入文件由
#     scripts/prep_inputs_clean.R 从 raw HiSeqV2.gz 派生（AIFM2 行），可被用户版本替换。

suppressMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

set.seed(20261028)

# ---- 1. 读输入（宽矩阵/长表均兼容） ----
expr <- fread("data/processed/TCGA_BLCA_expr_clean.csv")
grp0 <- fread("data/processed/TCGA_BLCA_group.csv")

smp_col <- names(expr)[toupper(names(expr)) %in% c("SAMPLE", "SAMPLEID", "SAMPLE_ID")][1]
if (is.na(smp_col)) smp_col <- names(expr)[1]            # 无命名首列按样本 ID
aifm2_col <- names(expr)[toupper(names(expr)) == "AIFM2"]
stopifnot(length(aifm2_col) == 1)
smp_expr <- as.character(expr[[smp_col]])
x <- as.numeric(expr[[aifm2_col]])

g_smp <- names(grp0)[toupper(names(grp0)) %in% c("SAMPLE", "SAMPLEID", "SAMPLE_ID")][1]
if (is.na(g_smp)) g_smp <- names(grp0)[1]
g_grp <- names(grp0)[2]
gv <- tolower(as.character(grp0[[g_grp]]))
gv <- ifelse(grepl("^tum", gv), "tumor", ifelse(grepl("^nor", gv), "normal", NA_character_))
stopifnot(!anyNA(gv))
stopifnot(setequal(as.character(grp0[[g_smp]]), smp_expr))
m <- match(as.character(grp0[[g_smp]]), smp_expr)
x <- x[m]
stopifnot(length(x) == 426, sum(gv == "tumor") == 407, sum(gv == "normal") == 19)

tum_all <- x[gv == "tumor"]; nor <- x[gv == "normal"]
P_full <- wilcox.test(tum_all, nor, exact = FALSE)$p.value
fc_full <- mean(tum_all) - mean(nor)
cat(sprintf("[full] P=%.7f  log2FC=%.6f  (tumor mean %.5f vs normal mean %.5f)\n",
            P_full, fc_full, mean(tum_all), mean(nor)))

# ---- 2. 1000 次降采样 ----
n_tum <- length(tum_all)
res <- lapply(seq_len(1000), function(run) {
  idx <- sample(n_tum, 19)               # 无放回
  tt <- tum_all[idx]
  p <- wilcox.test(tt, nor, exact = FALSE)$p.value
  fc <- mean(tt) - mean(nor)
  data.frame(run = run, p_value = p, log2FC = fc)
})
res <- rbindlist(res)
write.csv(res, "results/tables/downsampling_results.csv", row.names = FALSE)

med_p <- median(res$p_value)
ci95 <- quantile(res$p_value, c(0.025, 0.975))
prop_lt05 <- mean(res$p_value < 0.05) * 100
prop_lt10 <- mean(res$p_value < 0.10) * 100
fc_mean <- mean(res$log2FC); fc_med <- median(res$log2FC)
pctile <- 100 * mean(res$p_value <= P_full)
cat(sprintf("[downsampling] median P=%.3f (95%% %.3f-%.3f)  P<0.05: %.1f%%  P<0.10: %.1f%%\n",
            med_p, ci95[1], ci95[2], prop_lt05, prop_lt10))
cat(sprintf("[downsampling] log2FC mean=%.6f median=%.6f  range %.6f-%.6f\n",
            fc_mean, fc_med, min(res$log2FC), max(res$log2FC)))
cat(sprintf("[percentile] full P=%.4f 位于 1000 次 P 分布的 %.1f%% 分位（%d 次 ≤ 它）\n",
            P_full, pctile, sum(res$p_value <= P_full)))

# ---- 3. FigS1（左 P 直方图 / 右 log2FC 密度，cairo+Arial 300dpi） ----
gL <- ggplot(res, aes(x = p_value)) +
  geom_histogram(binwidth = 0.05, boundary = 0,
                 fill = "#56B4E9", color = "grey30", linewidth = 0.2) +
  geom_vline(xintercept = med_p, linetype = "dashed", color = "#D55E00") +
  geom_vline(xintercept = 0.05, linetype = "dotted", color = "grey40") +
  annotate("text", x = 0.55, y = Inf, vjust = 1.6, size = 3.2,
           label = sprintf("中位 P = %.3f\nP < 0.05：%.1f%%", med_p, prop_lt05)) +
  labs(title = "A  1000 次 19 vs 19 降采样 P 值", x = "Mann-Whitney P 值", y = "频数") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme_classic(base_family = "Arial", base_size = 11) +
  theme(plot.title = element_text(size = 11, face = "bold"))

gR <- ggplot(res, aes(x = log2FC)) +
  geom_density(fill = "#009E73", alpha = 0.5) +
  geom_vline(xintercept = fc_mean, linetype = "dashed", color = "#D55E00") +
  geom_vline(xintercept = fc_med, linetype = "dotted", color = "grey30") +
  geom_vline(xintercept = fc_full, linetype = "dotdash", color = "#0072B2") +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.6, size = 3.2,
           label = sprintf("均值 = %.4f\n中位数 = %.4f\n全样本 = %.4f", fc_mean, fc_med, fc_full)) +
  labs(title = "B  1000 次降采样 log2FC（肿瘤均值 − 正常均值）",
       x = "log2FC", y = "密度") +
  theme_classic(base_family = "Arial", base_size = 11) +
  theme(plot.title = element_text(size = 11, face = "bold"))

g <- gL + gR
ggsave("results/figures/FigS1_downsampling.pdf", g,
       width = 10, height = 4.2, device = cairo_pdf, family = "Arial")
ggsave("results/figures/FigS1_downsampling.png", g,
       width = 10, height = 4.2, dpi = 300, type = "cairo")
cat("[done] analysis1b (FigS1 v2 + downsampling_results.csv)\n")
