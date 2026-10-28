# install_packages.R — FSP1×BLCA 论文补充分析依赖安装
# 运行：Rscript install_packages.R
# 镜像：CRAN=aliyun、Bioconductor=tuna（2026-10-28 实测可达；R libcurl 直连失败，
# 改用 wininet 下载栈）

options(repos = c(CRAN = "https://mirrors.aliyun.com/CRAN"))
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
options(timeout = 900)
options(download.file.method = "wininet")

# 1) CRAN 侧
install.packages(c("BiocManager", "ppcor", "survminer", "patchwork"),
                 dependencies = TRUE)

# 2) Bioconductor 侧（GSVA 与 ComplexHeatmap 依赖链较长，逐个装便于定位失败）
library(BiocManager)
BiocManager::install(c("GSVA", "ComplexHeatmap"), ask = FALSE, update = FALSE)

# 3) 验证
ip <- rownames(installed.packages())
want <- c("Seurat", "GSVA", "ppcor", "ggplot2", "ComplexHeatmap",
          "survival", "survminer", "patchwork", "dplyr", "data.table")
cat("installed:", paste(want[want %in% ip], collapse = ", "), "\n")
cat("missing  :", paste(want[!want %in% ip], collapse = ", "), "\n")
if (all(want %in% ip)) cat("[OK] all packages ready\n") else cat("[WARN] some missing\n")
