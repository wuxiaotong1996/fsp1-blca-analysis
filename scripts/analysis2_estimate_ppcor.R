# analysis2_estimate_ppcor.R — 分析2：ESTIMATE 肿瘤纯度 + ppcor 偏相关
# 目的：① 计算 407 例肿瘤的 ESTIMATE 免疫/基质/肿瘤纯度评分；
#       ② 校正肿瘤纯度后重估 FSP1 与 22 类免疫浸润 ssGSEA 评分的相关
#          （Spearman 原始 ρ vs 偏 ρ），检验三条显著相关是否由纯度混杂驱动。
# ESTIMATE 算法：逐样本秩归一化（10000*rank/Ng）→ 加权 ES（权重=降序秩^0.25，
#   ES=ΣRES）→ ESTIMATEScore=Stromal+Immune → 纯度=cos(0.6049872018+0.0001467884*score)
#   —— 精确复刻 Yoshihara 2013 (Nat Commun 4:2612) 官方 estimate R 包 estimateScore.R
#   （包取自 SourceForge estimateproject 官方发布）。注：纯度公式官方标注为
#   Affymetrix 平台专用，RNA-seq 应用为文献通行做法，此处如实记录。
# 签名：data/processed/ESTIMATE_SI_geneset.gmt（141 基质 + 141 免疫，官方包原文件）
# 输入只读：data/raw/TCGA.BLCA.HiSeqV2.gz、data/processed/immune_scores.csv
# 输出：results/tables/Analysis2_estimate_scores.csv / Analysis2_purity_summary.csv /
#       Analysis2_ppcor_22class.csv；results/figures/FigS2_partial_corr.pdf/.png

suppressMessages({
  library(data.table)
  library(ppcor)
  library(ggplot2)
})

# ---- 1. 读 HiSeqV2（只读），肿瘤 407 例，基因符号去重取均值 ----
X <- as.data.frame(fread("data/raw/TCGA.BLCA.HiSeqV2.gz", nThread = 8))
rn <- toupper(as.character(X[[1]])); X[[1]] <- NULL
is_tum <- vapply(colnames(X), function(s) {
  p <- strsplit(s, "-", fixed = TRUE)[[1]]
  length(p) >= 4 && substr(p[4], 1, 2) == "01"
}, logical(1))
X <- X[, is_tum]
X <- as.matrix(X)
rownames(X) <- rn
if (anyDuplicated(rn)) {   # 该文件无重复符号（已查：0 重复），此分支为通用性保留
  g <- factor(rn)
  X <- rowsum(X, g) / as.numeric(table(g))
}
stopifnot(ncol(X) == 407, "AIFM2" %in% rownames(X))
fsp1 <- X["AIFM2", ]
cat("[input] tumor samples:", ncol(X), " genes:", nrow(X), "\n")

# ---- 2. ESTIMATE 评分（精确复刻 estimateScore.R） ----
gmt <- readLines("data/processed/ESTIMATE_SI_geneset.gmt")
sets <- lapply(strsplit(gmt, "\t"), function(f) f[-c(1, 2)])
names(sets) <- c("StromalSignature", "ImmuneSignature")
ov <- vapply(sets, function(s) sum(s %in% rownames(X)), integer(1))
cat("[ESTIMATE] signature overlap:", paste(names(ov), ov, sep = "=", collapse = ", "), "\n")

m <- 10000 * apply(X, 2, function(col) rank(col, ties.method = "average") / length(col))

est_es <- function(m, genes_set) {
  in_set <- rownames(m) %in% genes_set
  Nh <- sum(in_set); Ng <- nrow(m); Nm <- Ng - Nh
  ord <- apply(m, 2, order, decreasing = TRUE)
  w <- abs(apply(m, 2, sort, decreasing = TRUE))^0.25   # 降序秩^0.25
  tag01 <- matrix(in_set[ord], nrow = Ng) * 1           # 向量[矩阵下标]会丢 dim，须显式 matrix()
  sw <- colSums(tag01 * w)
  Pn <- tag01 * w / rep(sw, each = Ng)
  P0 <- (1 - tag01) / Nm
  Fn <- apply(Pn, 2, cumsum); F0 <- apply(P0, 2, cumsum)
  colSums(Fn - F0)                                      # ES = ΣRES
}
es <- sapply(sets, est_es, m = m)
est_score <- rowSums(es)
purity <- cos(0.6049872018 + 0.0001467884 * est_score)
purity[purity < 0] <- NA   # 官方"out of bounds"处理

sc <- data.frame(sample = colnames(X),
                 StromalScore = es[, 1], ImmuneScore = es[, 2],
                 ESTIMATEScore = est_score, TumorPurity = purity)
write.csv(sc, "results/tables/Analysis2_estimate_scores.csv", row.names = FALSE)

summ <- data.frame(
  n = sum(!is.na(purity)),
  mean = mean(purity, na.rm = TRUE), sd = sd(purity, na.rm = TRUE),
  median = median(purity, na.rm = TRUE),
  min = min(purity, na.rm = TRUE), max = max(purity, na.rm = TRUE))
write.csv(summ, "results/tables/Analysis2_purity_summary.csv", row.names = FALSE)
print(summ)

# ---- 3. 22 类 ssGSEA 原始 ρ vs 纯度校正偏 ρ ----
imm <- read.csv("data/processed/immune_scores.csv",
                fileEncoding = "UTF-8-BOM", check.names = FALSE)
smp <- as.character(imm[[1]])
stopifnot(all(smp %in% colnames(X)))
f <- as.numeric(fsp1[smp]); p <- as.numeric(purity[smp])
ok <- !is.na(p)

res <- do.call(rbind, lapply(2:ncol(imm), function(j) {
  y <- as.numeric(imm[[j]])
  r_raw <- cor.test(f, y, method = "spearman", exact = FALSE)
  idx <- ok
  pr <- pcor.test(f[idx], y[idx], p[idx], method = "spearman")
  data.frame(class = names(imm)[j],
             n = sum(idx),
             rho_raw = as.numeric(r_raw$estimate), p_raw = r_raw$p.value,
             rho_partial = as.numeric(pr$estimate), p_partial = pr$p.value)
}))
write.csv(res, "results/tables/Analysis2_ppcor_22class.csv", row.names = FALSE)

# 与已复核基线核对（旧项目 step2 标准式结果）
chk <- res[res$class %in% c("中性粒细胞", "浆样树突细胞", "细胞毒性细胞"),
           c("class", "rho_raw", "p_raw")]
print(chk)
stopifnot(abs(res$rho_raw[res$class == "中性粒细胞"] - (-0.141)) < 0.01,
          res$p_raw[res$class == "中性粒细胞"] < 0.01,
          abs(res$rho_raw[res$class == "浆样树突细胞"] - 0.140) < 0.01,
          res$p_raw[res$class == "浆样树突细胞"] < 0.01,
          abs(res$rho_raw[res$class == "细胞毒性细胞"] - 0.114) < 0.01)
cat("[check] 三条显著相关与基线一致\n")

# ---- 4. FigS2：原始 ρ vs 偏 ρ ----
dplot <- rbind(
  data.frame(class = res$class, type = "原始 ρ", rho = res$rho_raw, p = res$p_raw),
  data.frame(class = res$class, type = "纯度校正偏 ρ", rho = res$rho_partial, p = res$p_partial))
dplot$type <- factor(dplot$type, levels = c("原始 ρ", "纯度校正偏 ρ"))
dplot$class <- factor(dplot$class, levels = rev(res$class[order(res$rho_raw)]))

g <- ggplot(dplot, aes(x = rho, y = class, color = type, shape = type)) +
  geom_point(size = 2.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_text(data = dplot[dplot$p < 0.05, ],
            aes(label = "*"), color = "black", nudge_y = 0.25, show.legend = FALSE) +
  labs(title = "Fig S2  FSP1 表达与 22 类免疫浸润：原始 ρ vs 肿瘤纯度校正偏 ρ",
       x = "Spearman ρ", y = "", color = "", shape = "") +
  theme_classic(base_family = "Arial", base_size = 11) +
  theme(plot.title = element_text(size = 11, face = "bold"),
        legend.position = "top")

ggsave("results/figures/FigS2_partial_corr.pdf", g,
       width = 7, height = 6.5, device = cairo_pdf, family = "Arial")
ggsave("results/figures/FigS2_partial_corr.png", g,
       width = 7, height = 6.5, dpi = 300, type = "cairo")
cat("[done] analysis2\n")
