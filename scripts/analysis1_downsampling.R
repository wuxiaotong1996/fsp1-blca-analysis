# analysis1_downsampling.R — 分析1：TCGA 降采样分析
# 目的：肿瘤组 n=407 vs 正常组 n=19 的 Mann-Whitney 检验在 n 极端不均衡下
#       功效可疑；每次随机抽 19 例肿瘤 vs 19 例正常，重复 1000 次，
#       评估 P 值分布（检验"表达无差异"结论对样本量失衡的稳健性）。
# 输入：data/processed/fsp1_expression.csv（只读）
# 输出：results/figures/FigS1_downsampling.pdf / .png（300 dpi，Arial）
# 随机种子：20261028

suppressMessages({library(ggplot2)})

set.seed(20261028)

# ---- 读取（原始数据只读，不改写） ----
dat <- read.csv("data/processed/fsp1_expression.csv",
                fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
dat$expression <- as.numeric(dat$expression)
tum <- dat$expression[dat$group == "Tumor"]
nor <- dat$expression[dat$group == "Normal"]
stopifnot(length(tum) == 407, length(nor) == 19)

# 全样本观测 P 值（与已复核结果 P=0.7932264758 对照）
p_full <- wilcox.test(tum, nor, exact = FALSE)$p.value

# ---- 1000 次降采样 ----
# exact=FALSE：表达值含并列，强制正态近似（与 ties 时自动降级一致，避免 warning）
N_ITER <- 1000
ps <- vapply(seq_len(N_ITER), function(i) {
  s <- sample(tum, size = 19)          # 无放回抽 19 例肿瘤
  wilcox.test(s, nor, exact = FALSE)$p.value
}, numeric(1))

# ---- 汇总 ----
summ <- data.frame(
  n_tumor = 407L, n_normal = 19L, n_drawn = 19L, n_iter = N_ITER,
  p_full = p_full,
  p_median = median(ps),
  p_lt_005 = mean(ps < 0.05),
  p_lt_010 = mean(ps < 0.10),
  p_min = min(ps), p_max = max(ps),
  seed = 20261028L
)
write.csv(summ, "results/tables/FigS1_downsampling_summary.csv",
          row.names = FALSE)
print(summ)

# ---- 图：P 值分布 ----
dplot <- data.frame(p = ps)
g <- ggplot(dplot, aes(x = p)) +
  geom_histogram(binwidth = 0.05, boundary = 0,
                 fill = "#4477AA", color = "white") +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "#CC3311") +
  geom_vline(xintercept = p_full, linetype = "dotted", color = "#117733") +
  annotate("text", x = 0.13, y = Inf, vjust = 2, size = 3.2, family = "Arial",
           label = sprintf("P<0.05 比例 = %.1f%%\nP<0.10 比例 = %.1f%%\n中位 P = %.3f\n全样本观测 P = %.3f",
                           summ$p_lt_005 * 100, summ$p_lt_010 * 100,
                           summ$p_median, summ$p_full)) +
  labs(title = "Fig S1  TCGA 降采样分析：FSP1 表达 肿瘤 vs 正常（19 vs 19 × 1000 次）",
       subtitle = "每次从 407 例肿瘤中无放回抽取 19 例，与 19 例正常行 Mann-Whitney 检验",
       x = "Mann-Whitney P 值", y = "频数") +
  theme_classic(base_family = "Arial", base_size = 10) +
  theme(plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8.5))

# ---- 输出：PDF + PNG（300 dpi） ----
# 字体：Windows 版 R 的 windows 位图设备不认 "Arial"（会报
# "font family not found in Windows font database"），两处均强制
# cairo 设备；cairo win32 后端经 GDI 字体枚举直接解析系统 Arial。
ggsave("results/figures/FigS1_downsampling.pdf", g,
       width = 7, height = 4.5, device = cairo_pdf, family = "Arial")
ggsave("results/figures/FigS1_downsampling.png", g,
       width = 7, height = 4.5, dpi = 300, type = "cairo")

cat("[done] FigS1_downsampling.pdf/.png + summary.csv\n")
