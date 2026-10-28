# analysis4_ssgsea_24class.R — 分析4：补充表 S1 完整 ssGSEA（Bindea 2013 24 类免疫细胞）
# 目的：在 407 例 TCGA-BLCA 肿瘤上计算 Bindea 2013 Table S1 全部 24 类免疫细胞
#       ssGSEA 评分，报告每类与 FSP1(AIFM2) 的原始 Spearman ρ 与肿瘤纯度校正偏 ρ。
# 签名：data/processed/bindea24_signatures.csv（24 类、512 基因；来源 Bindea 2013
#       Immunity 39:782-795 Table S1——IOBR signature_collection 转写 23 类 + TimiGP
#       转录补 NK CD56bright 10 基因，NK CD56dim 恢复 Table S1 15 基因。原表
#       pDC/Treg 各仅 1 个、Th17 仅 3 个标志物（TimiGP 构建脚本注明），如实保留并报告）
# 方法：GSVA 2.6.14 gsva(method="ssgsea", normalize=FALSE)（.fastRndWalk 标准式，
#       与旧项目 step2 基线同式；已与手动 R 移植对照：三类 ρ/P 完全一致；
#       Spearman ρ 对评分尺度线性变换不变，normalize 开关不影响 ρ/P）。
#       原始/偏 ρ 均 Spearman（exact=FALSE 确定性），ppcor::pcor.test 校正肿瘤纯度
#       （results/tables/Analysis2_estimate_scores.csv）；BH FDR 按 24 类。
# 输出：results/tables/TableS1_ssGSEA_full.csv（n_genes×2 / ρ / P / FDR ×原始与偏）
#       results/tables/TableS1_ssgsea_scores_matrix.csv（407×24 评分矩阵，供复查）
# 输入只读：data/raw/TCGA.BLCA.HiSeqV2.gz、data/processed/bindea24_signatures.csv、
#           results/tables/Analysis2_estimate_scores.csv

suppressMessages({
  library(data.table)
  library(GSVA)
  library(ppcor)
})

set.seed(20261028)

# ---- 1. 读 HiSeqV2（与分析 2 同一预处理：407 肿瘤、符号去重） ----
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
f <- as.numeric(X["AIFM2", ])
cat("[input] tumor samples:", ncol(X), " genes:", nrow(X), "\n")

# ---- 2. 24 类签名 ----
sig <- read.csv("data/processed/bindea24_signatures.csv", stringsAsFactors = FALSE)
sets <- lapply(split(sig$Gene, sig$CellType), toupper)
sets <- sets[order(names(sets))]
n_sig <- lengths(sets)
n_used <- vapply(sets, function(s) sum(s %in% rownames(X)), integer(1))
stopifnot(length(sets) == 24, sum(n_sig) == 512, all(n_used >= 1))
cat("[signature] 24 sets, total genes:", sum(n_sig),
    " matrix intersection per set: min", min(n_used), "max", max(n_used), "\n")

# ---- 3. ssGSEA（GSVA .fastRndWalk；normalize=FALSE 原始 ES，ρ 不受尺度影响） ----
param <- ssgseaParam(X, sets, normalize = FALSE)
ss <- t(gsva(param, verbose = FALSE))          # 样本 × 类别

# 确定性锚点：与手动 fastRndWalk 移植（temp_check_ssgsea.R）对照的三类 ρ
chk <- sapply(c("Neutrophils", "pDC", "Cytotoxic cells"),
              function(nm) cor(f, ss[, nm], method = "spearman"))
stopifnot(abs(chk["Neutrophils"] - (-0.0501)) < 1e-3,
          abs(chk["pDC"] - 0.1763) < 1e-3,
          abs(chk["Cytotoxic cells"] - 0.1695) < 1e-3)
cat("[check] 三类 ρ 与手动移植对照一致（-0.0501 / 0.1763 / 0.1695）\n")

# ---- 4. 原始 ρ vs 纯度校正偏 ρ（24 类） ----
pur <- read.csv("results/tables/Analysis2_estimate_scores.csv", stringsAsFactors = FALSE)
purity <- setNames(pur$TumorPurity, pur$sample)[colnames(X)]   # 与列序对齐
ok <- !is.na(purity)
cat("[purity] usable samples:", sum(ok), "/", length(ok), "\n")

res <- do.call(rbind, lapply(names(sets), function(nm) {
  y <- as.numeric(ss[, nm])
  r_raw <- cor.test(f, y, method = "spearman", exact = FALSE)
  pr <- pcor.test(f[ok], y[ok], purity[ok], method = "spearman")
  data.frame(class = nm,
             n_genes_signature = n_sig[[nm]],
             n_genes_used = n_used[[nm]],
             rho_raw = as.numeric(r_raw$estimate),
             p_raw = r_raw$p.value,
             rho_partial = as.numeric(pr$estimate),
             p_partial = pr$p.value)
}))
res$fdr_raw <- p.adjust(res$p_raw, "BH")
res$fdr_partial <- p.adjust(res$p_partial, "BH")
res <- res[order(res$p_raw), ]
write.csv(res, "results/tables/TableS1_ssGSEA_full_detail.csv", row.names = FALSE)

sc_mat <- cbind(sample = rownames(ss), as.data.frame(ss))
write.csv(sc_mat, "results/tables/TableS1_ssgsea_scores_matrix.csv", row.names = FALSE)

print(res, digits = 3)
cat("[done] analysis4\n")
