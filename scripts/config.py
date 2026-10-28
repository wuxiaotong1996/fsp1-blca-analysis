# config.py
import os

# 全局配置
GENE_SYMBOL = "AIFM2"        # FSP1 的官方基因名（不要用 FSP1，数据库不识别）
GENE_ALIAS = "FSP1"          # 蛋白名，用于图表标题
CANCER_CODE = "BLCA"         # 膀胱癌 TCGA 代码
CANCER_NAME = "Bladder Cancer"

# GEPIA2 配置
GEPIA2_URL = "http://gepia2.cancer-pku.cn/"
GEPIA2_LOG2FC_CUTOFF = 1
GEPIA2_PVALUE_CUTOFF = 0.01

# 输出路径
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(BASE_DIR, "output")
FIGURES_DIR = os.path.join(OUTPUT_DIR, "figures")
DATA_DIR = os.path.join(OUTPUT_DIR, "data")

os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(DATA_DIR, exist_ok=True)

print(f"配置就绪：基因={GENE_SYMBOL}，癌种={CANCER_CODE}")
print(f"图片输出: {FIGURES_DIR}")
print(f"数据输出: {DATA_DIR}")

# ==================== 以下为三个脚本运行所需（在您的配置基础上追加） ====================
# 注：GEPIA2 的匿名 JSON 接口已下线（2026-10 实测 404），Step1 表达+生存改用
# 与 GEPIA2 同源的 UCSC Xena TCGA 数据（GEPIA2 底层数据同为 TCGA）；GEPIA2_URL 保留备用。
from pathlib import Path

ENTREZ_ID = 84883              # AIFM2 (FSP1) 的 NCBI Entrez Gene ID

# 数据源
EXPR_URL = "https://tcga.xenahubs.net/download/TCGA.BLCA.sampleMap/HiSeqV2.gz"
CLIN_URL = "https://tcga.xenahubs.net/download/TCGA.BLCA.sampleMap/BLCA_clinicalMatrix"
CBIOPORTAL_API = "https://www.cbioportal.org/api"    # cBioPortal REST API（Step3，无需登录）
STUDY_ID = "blca_tcga_pan_can_atlas_2018"            # TCGA BLCA PanCancer Atlas
TIMER2_URL = "https://timer.cistrome.org/"           # TIMER2.0 网页版人工交叉验证用

# 路径（转为 pathlib.Path 供脚本用 / 拼接；raw/ 为下载原始文件缓存）
DATA_DIR = Path(DATA_DIR)
FIGURES_DIR = Path(FIGURES_DIR)
RAW_DIR = DATA_DIR / "raw"
RAW_DIR.mkdir(parents=True, exist_ok=True)

# 网络
REQUEST_TIMEOUT = 60
USER_AGENT = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) fsp1-blca-analysis/1.0"
RETRIES = 2

# 分析参数
EXPR_GROUP_CUT = "median"    # 生存分析分组：按肿瘤样本 FSP1 表达中位数分高/低组
SURV_END = "OS"              # 生存终点（OS / RFS，列缺失时自动回退）
SSGSEA_ALPHA = 0.25          # ssGSEA 秩权重指数（Barbie 2009 默认）
CORR_METHOD = "spearman"     # 免疫浸润相关性方法
