# prep_inputs_clean.R — 派生任务所需的两个输入文件（raw 只读）
# 背景：任务指定的 data/processed/TCGA_BLCA_expr_clean.csv 与 TCGA_BLCA_group.csv
#       在磁盘上不存在（已全盘检索）。其内容可由 raw 完全推导，故从
#       data/raw/TCGA.BLCA.HiSeqV2.gz 派生（AIFM2 行 = 第 11,443 行），
#       分组按样本条码第 4 段 01=肿瘤 / 11=正常。若用户另有版本，直接覆盖这两个
#       文件后重跑 analysis1b_downsampling.R 即可（该脚本对宽/长格式均兼容）。
# 输出：data/processed/TCGA_BLCA_expr_clean.csv（sample, AIFM2）
#       data/processed/TCGA_BLCA_group.csv（sample, group: tumor/normal）

suppressMessages(library(data.table))

X <- fread("data/raw/TCGA.BLCA.HiSeqV2.gz", nThread = 8)
rn <- toupper(as.character(X[[1]]))
i_aifm2 <- which(rn == "AIFM2")
stopifnot(length(i_aifm2) == 1)
v <- as.numeric(X[i_aifm2, -1])
samples <- names(X)[-1]

grp <- vapply(samples, function(s) {
  p <- strsplit(s, "-", fixed = TRUE)[[1]]
  if (length(p) >= 4 && substr(p[4], 1, 2) == "01") "tumor" else
  if (length(p) >= 4 && substr(p[4], 1, 2) == "11") "normal" else NA_character_
}, character(1))
stopifnot(!anyNA(grp), sum(grp == "tumor") == 407, sum(grp == "normal") == 19)

fwrite(data.table(sample = samples, AIFM2 = v),
       "data/processed/TCGA_BLCA_expr_clean.csv")
fwrite(data.table(sample = samples, group = grp),
       "data/processed/TCGA_BLCA_group.csv")
cat("[prep] expr_clean:", length(v), "rows (tumor", sum(grp == "tumor"),
    "normal", sum(grp == "normal"), ") written\n")
