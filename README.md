# fsp1_blca_paper — FSP1×膀胱癌生信分析论文工作区

FSP1（铁死亡抑制蛋白1，基因名 **AIFM2**，Entrez 84883）在膀胱尿路上皮癌（BLCA）中的**表达、生存、免疫浸润、遗传变异**分析，附 scRNA 层面 Treg vs CD8⁺T 表达的补充分析。全部结果经 2026-10-15 独立复核通过；Step2 已按用户选定方案改为 GSVA fastRndWalk 标准式并重跑（2026-10-15）。

## 目录结构

```
fsp1_blca_paper\
├─ data\raw\           Xena 原始下载（HiSeqV2.gz 22.8MB、clinicalMatrix）
├─ data\processed\     五张结果表 + GEPIA2 网页版证据（PDF×4+JSON）
├─ scripts\            7 个脚本 + config.py（自定位路径，可整体搬迁重跑）
├─ results\figures\    10 张结果图（step1-3 + GEPIA2 交叉验证 + fig1 拼图）
├─ results\tables\     step2 相关表 + 三步 summary JSON
└─ supplementary\      旧 README/SUMMARY/REPORT + scRNA 补充项目（脚本+结果）
```

> 旧项目原件在 `课题项目 申请\预实验，小课题\工作过程\数据分析\fsp1_blca_analysis\`，未删除，仍是可运行的工作副本。

## 分析流程（三步，全部跑通且复核全绿）

| 步骤 | 脚本 | 数据源 | 输出 |
|---|---|---|---|
| Step1 表达+生存 | `step1_gepia2.py` | Xena TCGA BLCA（GEPIA2 同源） | `fsp1_expression.csv`、`blca_survival.csv`、图 `step1_*` |
| Step1 交叉验证 | `step1_gepia2_browser.py` | GEPIA2 网页 UI（browser-use 驱动） | `gepia2_browser\`（PDF 文本层核实，非 OCR） |
| Step2 免疫浸润 | `step2_timer2.py` | Xena 表达矩阵 + 22 类 ssGSEA（**GSVA fastRndWalk 标准式**） | `immune_scores.csv`、`step2_correlation.csv`、图 `step2_*` |
| Step3 遗传变异 | `step3_cbioportal.py` | cBioPortal REST（blca_tcga_pan_can_atlas_2018） | `fsp1_mutations.csv`、`fsp1_cna.csv`、图 `step3_*` |
| 一键重跑 | `run_all.py` | — | 全套输出（缓存命中则跳过下载） |

## 核心结果（真实数值，已独立复核）

- **表达**：肿瘤 407 vs 正常 19，中位数 8.9795/8.977，Mann-Whitney **P=0.79 无差异**；GEPIA2 网页版交叉验证 num(T)=404、num(N)=28。
- **生存**（Xena 版，中位切点）：log-rank **P=0.15，HR=0.80（95%CI 0.60–1.08）**；GEPIA2 网页版 **HR(high)=0.86、logrank p=0.3、n=201/201**——方向一致，均不显著。
- **免疫浸润**（ssGSEA，GSVA 标准式更新后）：
  - 中性粒细胞 **ρ=−0.141, P=0.004**（负相关）
  - 浆样树突细胞 **ρ=0.140, P=0.005**（正相关）
  - 细胞毒性细胞 **ρ=0.114, P=0.022**（正相关）
  - CD8⁺T ρ=0.099, P=0.046（边缘）
- **遗传变异**（cBioPortal）：突变 3/408（0.7%；S362C/H104Y/L241V）、扩增 3、深缺失 2、浅层缺失 139（34%）。
- **scRNA 补充**（GSE149652，Treg n=3,813 vs CD8⁺T n=5,585）：AIFM2 阳性率 0.76% vs 0.64%，**Wilcoxon p=0.51 不显著**——采用"差异有限"诚实表述（见 `supplementary\scRNA_AIFM2_Treg_CD8T\results\README.md`）。

## 论文写作口径（诚实原则）

1. 表达无差异、生存方向一致但不显著——**不要**写成"高表达预后差/表达上调"。
2. 三条显著免疫相关（中性粒−/浆样树突+/细胞毒+）是当前最可写的点。
3. Treg vs CD8 差异不显著——基金申请书中"Treg 0.76% 阳性"只能作为背景事实引用，不能作为"Treg 富集"证据。

## 复现方法

```powershell
cd scripts
python run_all.py          # 缓存命中跳过下载；重下需 python download_xena.py
```

已知坑（勿重踩）：Xena 老域名 `tcga.xenahubs.net` 单连接 ~18KB/s，需 8 线程 Range 并行；cBioPortal mutations 端点必须传 `sampleListId`；browser-use 须 `enable_default_extensions=False` 且关闭 thinking；GEPIA2/TIMER API 已下线，只能走网页 UI 或 Xena 同源数据。

## 论文补充分析（2026-10-28 全部完成，目标期刊 BMC Cancer / Journal of Cancer）

五项任务书分析（#3 与 #5 合并），全部只读原始数据、种子 20261028、图 300 dpi PDF+PNG 双格式 Arial、逐项日志见 `analysis_log.md`：

| # | 分析 | 脚本 | 核心结果 | 输出 |
|---|---|---|---|---|
| 1 | TCGA 降采样稳健性 | `scripts/analysis1_downsampling.R` / `analysis1b_downsampling.R`（v2 审稿回应版） | 1000×19vs19 中位 P=0.559，P<0.05 仅 3.2%，原始 P 位于第 76 百分位 → "表达无差异"对 n 失衡稳健 | FigS1_downsampling、downsampling_results.csv |
| 2 | ESTIMATE 纯度 + 偏相关（22 类） | `scripts/analysis2_estimate_ppcor.R` | 中性粒−/浆样树突+ 校正后存续；细胞毒/CD8 校正后不显著；纯度 n=342 可算 | FigS2_partial_corr、Analysis2_* |
| 3+5 | scRNA 全谱系（GSE135337，8 例） | `scripts/analysis35_scRNA_lineages.R` | 42,455 细胞 → 肿瘤 34,034/巨噬 1,652/内皮 659/成纤维 3,387/T 759/Other 1,964；仅可视化，无按细胞统计检验 | Fig2A_UMAP、Fig2B_dotplot、GSE135337_seurat.rds |
| 4 | Table S1 24 类 ssGSEA（Bindea 2013） | `scripts/analysis4_ssgsea_24class.R` | Tcm ρ=−0.238/FDR=2.7e-5 最强；NK CD56bright/Tgd/Cytotoxic/CD8/pDC 校正后存续；中性粒校正后负相关显现（两套签名体系同向） | TableS1_ssGSEA_full、TableS1_ssgsea_scores_matrix |
| 2b | 纯度重评估（官方 ESTIMATE 包 + 24 类偏相关，回应"纯度混杂"） | `scripts/analysis2b_purity_24class.R` | 包 vs 手工移植交叉校验 max|diff|≈1e-12；FSP1×Immune/Stromal 评分均 ns（ρ=0.051/0.025）；24 类正相关校正后全存续，中性粒负相关显现（ρ=−0.168/P=0.002）；原始 CD8⁺T ρ=0.099/P=0.046 → 校正 ρ=0.080/P=0.138 | TableS1_ssGSEA_full（任务书规格版）、Analysis2b_*、Fig1C_partial_correlation |
| 6 | GSE149652（TISCH2）AIFM2 表达谱 + GSE135337 全谱系（回应"只分析T细胞不足以证明异质性"） | `scripts/analysis6_gse149652_tcell.R` / `analysis6b_gse135337_aifm2.R` | TISCH2 注释仅 T/NK 6 类（**无 Malignant/Macrophage 等，已报告缺失**）；GSE149652 全亚群阳性率 <2%（Treg 0.76%/CD8⁺T 0.64%，旧值精确复现）；GSE135337 全谱系 5–9%（肿瘤上皮 8.5%/内皮 9.4%/巨噬 5.4%）→ FSP1 非 T 细胞特异、广泛多谱系表达；跨数据集不可直接互比（测序深度） | Fig2B_dotplot_AIFM2（GSE135337）、FigS4A/B/C（GSE149652）、Fig2_source_data.csv、FigS4_source_data_GSE149652.csv |

## 待办

- [ ] 用户教程文档中 step2 代码块仍是旧归一化公式，需替换为 GSVA fastRndWalk 标准式（用户侧）
- [ ] 论文骨架与图表升级（下一步：I-M-Ra-D 骨架、fig1_composite 拆分/重建、期刊选择）
- [ ] 补充分析结果写入论文正文/图表（用户侧改稿）
- [ ] scRNA 大数据文件（rds/mtx/gz，共数百 MB）未复制，仍在旧项目 `data\` 下，如需引用再搬运
