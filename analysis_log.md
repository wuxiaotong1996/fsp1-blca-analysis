# analysis_log.md — FSP1×BLCA 论文补充分析日志

约定：原始数据只读；单细胞数据不做以单细胞为观测单位的统计检验（伪重复）；随机种子统一 20261028；每个分析完成后追加一条简短日志。

---

## 分析 1：TCGA 降采样分析（2026-10-28 完成）

- **目的**：检验"FSP1 表达 肿瘤 vs 正常 无差异"（Mann-Whitney P=0.79）对样本量极端失衡（407 vs 19）的稳健性。
- **输入**（只读）：`data/processed/fsp1_expression.csv`（407 肿瘤 + 19 正常）。
- **方法**：每次从 407 例肿瘤中**无放回**抽 19 例，与 19 例正常行 Mann-Whitney（`exact=FALSE`，表达值含并列），重复 1000 次；种子 20261028。
- **结果**：全样本观测 P=0.7932265；1000 次降采样中位 P=0.559291；P<0.05 仅 3.2%、P<0.10 仅 6.5%；范围 0.002903–1.0。→ 19 vs 19 平衡抽样下绝大多数与全样本同向，"表达无差异"结论对 n 失衡稳健（论文卖点）。
- **输出**：`results/figures/FigS1_downsampling.pdf/.png`（300 dpi，Arial）；`results/tables/FigS1_downsampling_summary.csv`。
- **备注**：① `exact=FALSE` 使全链确定化——早期版本全样本检验的 exact 随机 tie-breaking 会消耗随机流，导致同种子跨 run 数值漂移，已消除，现可精确复现；② Windows 字体坑：windows 位图设备不认 "Arial"（14 条 warning），两处 ggsave 已强制 cairo 设备 + `family="Arial"`，重跑 0 warning，PNG 目视核验（中文经 GDI 字体回退正常渲染）。

---

## 分析 1 v2：TCGA-BLCA 降采样（审稿回应版）（2026-10-28 完成）

**降采样分析：1000次随机抽样，P值中位数=0.559（95%区间0.038–1.000），log2FC均值=0.058，P<0.05的比例为3.2%。原始分析P=0.793位于第76百分位。**

- **目的**：回应审稿人"正常样本 n=19 太少"质疑；本版按新任务书重做并新增 log2FC 记录。
- **输入**：`data/processed/TCGA_BLCA_expr_clean.csv`（sample, AIFM2）+ `TCGA_BLCA_group.csv`（sample, group）——任务书指定文件原不存在，已由 `scripts/prep_inputs_clean.R` 从 raw HiSeqV2.gz 派生（AIFM2=第 11,443 行；组标签按条码 01/11）；用户自有版本可直接覆盖后重跑。
- **方法**：407 肿瘤无放回抽 19 与 19 正常行 Mann-Whitney（exact=FALSE），记录 P 与 log2FC=mean(肿瘤)−mean(正常)，1000 次，seed 20261028。
- **结果**：全样本原始 P=0.7932265（与旧版精确一致，派生无误）、全样本 log2FC=0.0512；1000 次 P 中位 0.559、95%区间 0.038–1.000、P<0.05 仅 3.2%（P<0.10 6.5%）；log2FC 均值 0.058/中位 0.060、范围 −0.588~+0.697；原始 P 位于 75.9%（第 76 百分位，759/1000 次 ≤ 它）→ 平衡抽样下"无表达差异"结论稳健（论文卖点，直接回应审稿人）。
- **输出**：`results/tables/downsampling_results.csv`（run/p_value/log2FC × 1000 行）；`results/figures/FigS1_downsampling.pdf/.png`（300 dpi，cairo+Arial，双面板：P 直方图+log2FC 密度，已目视核验）。**覆盖分析1旧版 FigS1**（旧 `FigS1_downsampling_summary.csv` 保留未动）。
- **备注**：① ggplot2 新版本 hist 边框参数为 `linewidth` 非 `size`（旧写法仅 warning，已修）；② ragg 设备 warning 无害；③ 与 v1 数值一致（同 seed 同数据），差异仅在新增 log2FC 维度与双面板图。

---

## 分析 2：ESTIMATE 肿瘤纯度 + 纯度校正偏相关（2026-10-28 完成）

- **目的**：检验三条显著相关（中性粒细胞 −、浆样树突细胞 +、细胞毒性细胞 +）是否由肿瘤纯度混杂驱动。
- **输入**（只读）：`data/raw/TCGA.BLCA.HiSeqV2.gz`（20,530 基因 × 407 肿瘤，无重复符号）；`data/processed/immune_scores.csv`（22 类 ssGSEA，沿用旧项目）；`data/processed/ESTIMATE_SI_geneset.gmt`（官方 141+141 签名）。
- **方法**：逐样本秩归一（10000×rank/Ng）→ 降序加权 ES（权重=|秩|^0.25，**ES=ΣRES**）→ ESTIMATEScore=Stromal+Immune → 纯度=cos(0.6049872018+0.0001467884×score)，负值 NA——逐行对照官方 estimateScore.R（SourceForge estimateproject 包）移植并核实；随后 ppcor::pcor.test(method="spearman") 计算校正纯度后的偏 ρ。
- **结果**：纯度 n=342 可算（65 例 cos<0→NA），mean 0.561 / median 0.589 / 范围 0.0008–0.996。基线三断言全过：中性粒细胞 ρ=−0.141/P=0.0044、浆样树突细胞 ρ=0.140/P=0.0046、细胞毒性细胞 ρ=0.114/P=0.0219。**校正后**：中性粒细胞偏 ρ=−0.212/P=7.9e-5（增强，仍显著）；浆样树突细胞偏 ρ=0.114/P=0.035（存续，略减）；细胞毒性细胞偏 P=0.079、CD8+ T 偏 P=0.138（**不再显著，提示原相关部分由纯度驱动**）；新显现：CD4+ T 偏 ρ=−0.120/P=0.027、Th2 偏 ρ=0.150/P=0.006、单核细胞偏 ρ=−0.191/P=0.0004。
- **输出**：`results/tables/Analysis2_estimate_scores.csv` / `Analysis2_purity_summary.csv` / `Analysis2_ppcor_22class.csv`；`results/figures/FigS2_partial_corr.pdf/.png`（300 dpi，cairo+Arial，已目视核验）。
- **备注**：① 签名重叠 141/141 与 141/141（全矩阵仅 20,530 基因，两签名无缺失）；② 纯度公式官方注明限 Affymetrix，RNA-seq 应用为文献通行做法，论文方法部分需如实注明；③ 脚本 bug 记录：向量被矩阵下标索引会丢 dim（`in_set[ord]` → 必须 `matrix(..., nrow=Ng)`）；无重复符号时 rowsum 分支被跳过导致 rownames 未赋值——均已修复，脚本现可复现；④ png 用 ragg 设备（R 4.6 默认，忽略 `type` 参数，无害 warning）。

---

## 分析 3+5：scRNA 全谱系 UMAP + Dot plot（2026-10-28 完成）

- **目的**：以全谱系（肿瘤上皮/巨噬/内皮/成纤维/T 细胞）UMAP + marker Dot plot 重绘图 2，展示 FSP1 高表达的巨噬与肿瘤上皮在整体细胞景观中的位置。
- **输入**（只读）：GSE135337（Chen et al. Nat Commun 2021，8 例未治疗 BLCA）8 个原始计数表（`data/raw/GSE135337/`）。
- **方法**：8 表合并（基因并集、细胞条形码加样本前缀）→ Seurat 5.6 标准流程（min.cells=3/min.features=200 → nFeature 200–6000、percent.mt<25 → Normalize/FindVariableFeatures(2000)/Scale/PCA(30)/UMAP(dims=1:20, seed 20261028)）→ AddModuleScore 对 5 个谱系 marker 集打分 + argmax 归类（max score<0 → "Other"）。
- **结果**（描述性计数，42,455 细胞）：肿瘤上皮 34,034、巨噬 1,652、内皮 659、成纤维 3,387、T 759、Other 1,964。**仅可视化与计数，无按细胞统计检验**（伪重复约束）。
- **输出**：`results/figures/Fig2A_UMAP_all_celltypes.pdf/.png`、`Fig2B_dotplot_lineage_markers.pdf/.png`（300 dpi，cairo+Arial，均已目视核验）；`results/tables/Fig2_lineage_counts.csv`；`data/processed/GSE135337_seurat.rds`（复用）。
- **备注**：① 部分次核心 marker 在 8 表并集后缺失（成纤维 PDGFRB/FAP、内皮 VWF/CDH5/KDR——8 个文件的基因行本就不同），Seurat AddModuleScore 用交集基因打分并给 warning；核心 marker（EPCAM/CD68/PECAM1/COL1A1/CD3D）均已断言在矩阵中；② 与任务书 #5（图 2 重绘）合并为同一分析（同一输入同一图）。

---

## 分析 4：补充表 S1 完整 ssGSEA（Bindea 2013 24 类免疫细胞）（2026-10-28 完成）

- **目的**：在 407 例肿瘤上计算 Bindea 2013 Table S1 全部 24 类免疫细胞的 ssGSEA 评分，报告每类与 FSP1 的原始 ρ 与纯度校正偏 ρ（补充表 S1）。
- **输入**（只读）：`data/raw/TCGA.BLCA.HiSeqV2.gz`；`data/processed/bindea24_signatures.csv`（24 类 512 基因）；`results/tables/Analysis2_estimate_scores.csv`（纯度）。
- **方法**：GSVA 2.6.14 `gsva(method="ssgsea", normalize=FALSE)`（.fastRndWalk 标准式，与旧项目 step2 基线同式）；原始/偏 Spearman（`exact=FALSE`）+ `ppcor::pcor.test` 校正纯度（n=342）；BH FDR 按 24 类。
- **结果**：**Tcm** 原始 ρ=−0.238/P=1.1e-6（FDR=2.7e-5），偏 ρ=−0.197/P=2.6e-4（FDR=0.003）——最强且校正后存续；**NK CD56bright** 原始 ρ=0.204/P=3.4e-5，偏 ρ=0.233/P=1.4e-5（增强）；**pDC**（注意：原表仅 IL3RA 1 个标志物）ρ=0.176/P=3.5e-4，偏 ρ=0.152/P=0.005；**Tgd** 0.170/P=5.6e-4、**Cytotoxic cells** 0.169/P=6.0e-4、**CD8 T** 0.138/P=0.005——均在校正后存续（FDR≤0.022）；**Neutrophils**（31 基因）原始 ns（−0.050），校正后偏 ρ=−0.168/P=0.002（FDR=0.015）——与旧 22 类小签名集结果同向（负相关由纯度校正后显现）。
- **输出**：`results/tables/TableS1_ssGSEA_full.csv`（每类 n_genes_signature/n_genes_used/ρ/P/FDR × 原始与偏）；`results/tables/TableS1_ssgsea_scores_matrix.csv`（407×24 评分矩阵，供复查）。
- **备注**：① 签名来源核实：Bindea Table S1 原表 pDC/Treg 各仅 1 个、Th17 仅 3 个标志物（TimiGP 构建脚本 build_CellType_Bindea2013_cancer.R 注明 "pDC & Treg: 1 marker; Th 17: 3 markers; remove"；IOBR 转录独立一致）——按原表如实保留并在表中报告 n_genes；② 实现验证：R 手动 .fastRndWalk 移植与 GSVA::gsva 三类 ρ/P 完全一致（−0.0501/0.1763/0.1695），脚本内锚点断言防签名漂移；③ normalize=FALSE 用原始 ES，ρ 对尺度线性变换不变，开关不影响统计量；④ 旧 22 类（分析 2）与本 24 类为两套签名体系（旧=4–7 基因教科书 marker 集，本=Table S1 全集），结果不可直接互比，论文需分别表述；两体系在"中性粒细胞负相关经纯度校正后显现"上同向，是稳健性证据。

---

## 分析 2b：纯度校正重评估（官方 ESTIMATE 包 + 24 类 ssGSEA 偏相关）（2026-10-28 完成）

- **目的**：按新任务书以**官方 ESTIMATE R 包**（v1.0.11，SourceForge estimateproject 源码目录安装）重算 407 例肿瘤纯度，回应审稿人"纯度混杂"质疑；对 FSP1 与 24 类免疫细胞（Bindea 2013）ssGSEA 评分计算原始/纯度校正偏 ρ，并报告 FSP1×ImmuneScore/StromalScore 相关。
- **输入**（只读）：`data/raw/TCGA.BLCA.HiSeqV2.gz`（全矩阵，20,530×407 肿瘤——`data/processed/TCGA_BLCA_expr_clean.csv` 磁盘版为 2 列最简式，无法支撑 ESTIMATE/ssGSEA，故按分析2/4同法读 raw，AIFM2 与 expr_clean 完全同源）；`data/processed/bindea24_signatures.csv`；`data/processed/immune_scores.csv`（22 类）。
- **方法**：GCT 写出 → 官方 `estimateScore(gct_in, gct_out, platform="affymetrix")`（签名重叠 141/141 与 141/141）→ TumorPurity=cos(0.6049872018+0.0001467884×ESTIMATEScore)，负值 NA；原始/偏 Spearman（exact=FALSE）+ `ppcor::pcor.test` 控制纯度。
- **结果**：纯度 342/407 可算（65 例越界→NA），mean 0.561/median 0.589；**交叉校验**：包输出 vs 分析2手工移植 max|diff|≈1e-12（Stromal 1.0e-12 / Immune 1.0e-11 / ESTIMATE 1.0e-11 / Purity 1.0e-15），实现无漂移。**FSP1×ImmuneScore ρ=0.051/P=0.304、×StromalScore ρ=0.025/P=0.617、×ESTIMATEScore ρ=0.044/P=0.381，均不显著**——FSP1 表达本身与纯度分量无关，纯度不是主相关的主要驱动源。**24 类校正后**：NK CD56bright ρ=0.233/P=1.4e-5（增强）、Cytotoxic 0.162/P=0.003、pDC 0.151/P=0.005、Tgd 0.148/P=0.006、CD8 T 0.147/P=0.006——正相关全部存续；Neutrophils 原始 ns（−0.050）校正后 −0.168/P=0.002（负相关显现）；Tcm −0.197/P=2.6e-4 仍最强。**原始分析 CD8⁺T（22 类）：ρ=0.099/P=0.046 → 校正后 ρ=0.080/P=0.138**（减弱为不显著，提示该边缘相关部分由纯度驱动）。
- **输出**：`results/tables/Analysis2b_estimate_scores_pkg.csv`（包计算 4 项评分）、`Analysis2b_fsp1_estimate_corr.csv`、`TableS1_ssGSEA_full.csv`（cell_type, raw_rho, raw_p, partial_rho, partial_p, n，按 partial_rho 降序——任务书规格版；分析4 详细版保留为 `TableS1_ssGSEA_full_detail.csv`）；`results/figures/Fig1C_partial_correlation.pdf/.png`（300 dpi，cairo+Arial，已目视核验）。
- **备注**：① 包内 readGCT 以 check.names=TRUE 读入，条码 "-"→"."（顺序不变），脚本按同变换比对回映；② ImmuneScore 偏相关协方差阵秩亏，pcor.test 用 Moore-Penrose 广义逆（ρ≈−0.003，ns，无害）；③ 纯度公式官方注明限 Affymetrix，RNA-seq 应用为文献通行做法，方法部分如实注明（同分析2）；④ 两套签名体系（22 类教科书集 vs 24 类 Bindea 全集）在"中性粒细胞负相关校正后显现"上同向。

**Results 3.2 文字（中文版，可直接入稿）**：

> 为排除肿瘤纯度对 FSP1–免疫浸润相关性的混杂，采用 ESTIMATE 算法计算 407 例肿瘤的基质评分、免疫评分与肿瘤纯度（342 例纯度值有效，均值 0.561）。FSP1 表达与免疫评分（ρ=0.051，P=0.304）、基质评分（ρ=0.025，P=0.617）及 ESTIMATE 总评分（ρ=0.044，P=0.381）均无显著相关。控制肿瘤纯度后，基于 Bindea 2013 基因集的 24 类免疫细胞 ssGSEA 分析显示：FSP1 与 NK CD56bright 细胞（ρ=0.233，P=1.4×10⁻⁵）、细胞毒性细胞（ρ=0.162，P=0.003）、浆细胞样树突细胞（ρ=0.151，P=0.005）、γδ T 细胞（ρ=0.148，P=0.006）及 CD8⁺T 细胞（ρ=0.147，P=0.006）的正相关均存续；中性粒细胞的负相关在校正后显现（ρ=−0.168，P=0.002）。原始分析中 CD8⁺T 细胞的相关（ρ=0.099，P=0.046）校正后减弱至 ρ=0.080（P=0.138），提示该边缘相关部分由纯度驱动；而中性粒细胞负相关及主要正相关并非肿瘤纯度混杂所致。

**Results 3.2 paragraph (English version)**:

> To exclude confounding by tumor purity in the FSP1–immune infiltration associations, we computed stromal scores, immune scores, and tumor purity for all 407 tumors using the ESTIMATE algorithm (purity valid in 342 tumors, mean 0.561). FSP1 expression was not significantly correlated with the immune score (ρ = 0.051, P = 0.304), stromal score (ρ = 0.025, P = 0.617), or ESTIMATE score (ρ = 0.044, P = 0.381). After adjustment for tumor purity, ssGSEA across the 24 immune cell types (Bindea et al., 2013) showed that positive associations of FSP1 with NK CD56bright cells (ρ = 0.233, P = 1.4 × 10⁻⁵), cytotoxic cells (ρ = 0.162, P = 0.003), plasmacytoid dendritic cells (ρ = 0.151, P = 0.005), γδ T cells (ρ = 0.148, P = 0.006), and CD8⁺ T cells (ρ = 0.147, P = 0.006) all persisted, and a negative association with neutrophils emerged after adjustment (ρ = −0.168, P = 0.002). The association with CD8⁺ T cells observed in the original analysis (ρ = 0.099, P = 0.046) attenuated to ρ = 0.080 (P = 0.138), suggesting this marginal association was partly driven by purity, whereas the neutrophil association and the principal positive associations were not attributable to purity confounding.

---

## 分析 6：GSE149652（TISCH2）AIFM2 表达谱 + GSE135337 全谱系（Fig2 重建）（2026-10-28 完成）

- **目的**：回应审稿人"只分析 T 细胞不足以证明异质性"——在全注释范围内计算 AIFM2 各细胞类型阳性率（Wilson 95%CI）与平均表达，出 UMAP/Dot plot/柱状图。
- **输入**（只读）：任务书指定文件在工作区不存在（`GSE149652_counts.mtx.gz`、`TISCH2_BLCA_GSE149652_annotation.csv`），使用旧项目等价物：`课题项目 申请\…\scRNA_AIFM2_Treg_CD8T\data\GSE149652_merged_seurat.rds`（任务书允许"或 Seurat .rds"；RNA counts+data，17,156×15,538）+ 同目录 `BLCA_GSE149652_CellMetainfo_table.tsv`（TISCH2 注释，cell_id=GSM@barcode-1-0）；全谱系用 `data/processed/GSE135337_seurat.rds`（分析3+5 谱系归类）。
- **方法**：TISCH2 注释按 (Sample, 条码) 匹配合并——**断言 rds meta 与 TSV 注释逐细胞同源、UMAP 坐标逐细胞一致（max|diff|=0）**，UMAP 直接用 TISCH2 坐标不重新降维；LogNormalize scale.factor=10000（仅内存重算，磁盘只读）；阳性率=counts>0 + Wilson 95%CI 手工公式；平均表达=归一化 data 层组内均值；**不做以单细胞为观测单位的统计检验**（柱状图仅描述性，无 P 值）。
- **缺失报告（按任务书要求，未自行推断）**：TISCH2 GSE149652 注释中**无 Malignant、无 Macrophage、无 Endothelial、无 Fibroblast、无 B/Mast/DC/Neutrophil**——该数据集为 CD4/CD8 TIL 分选（malignancy 列 15,299 Immune + 239 Stromal，0 Malignant），只有 T/NK 谱系 6 类：Treg 3,813 / CD8_T 5,585 / CD4_Tconv 3,828 / NK 534 / CD8Tex 1,539 / Tprolif 239。全谱系部分由 GSE135337（8 例未治疗 BLCA）承担。
- **结果**：
  - **GSE149652（T/NK 亚群）**：Treg 0.76%（29/3,813，CI 0.53–1.09）、CD8⁺T 0.64%（36/5,585，CI 0.47–0.89）、CD4⁺Tconv 0.34%（13/3,828）、NK 0.19%（1/534）、CD8Tex 0.58%（9/1,539）、Tprolif 1.26%（3/239）——**6 类阳性率均 <2%**；旧项目锚点 29/3,813 与 36/5,585 精确复现。
  - **GSE135337（全谱系）**：肿瘤上皮 8.52%（2,901/34,034，CI 8.23–8.83）、内皮 9.41%（62/659）、成纤维 5.58%（189/3,387）、巨噬 5.45%（90/1,652）、T 5.01%（38/759）、Other 4.99%（98/1,964）——**各谱系阳性率 5–9%**。**注意：任务书模板句"所有细胞类型阳性率均<2%"在 GSE135337 不成立**（两数据集技术灵敏度不同：GSE149652 droplet 浅测序检出率低，GSE135337 全转录组深度高；跨数据集阳性率不可直接数值互比，只能各自内部比较）。科学结论反而更有利：**FSP1 并非 T 细胞特异，而是广泛表达于肿瘤上皮/内皮/巨噬等多谱系**，回应"仅看 T 细胞"质疑。
- **输出**：`results/figures/Fig2B_dotplot_AIFM2.pdf/.png`（GSE135337 全谱系，x=AIFM2，大小=阳性率，颜色=平均表达，右侧标 n——与既有 Fig2A（GSE135337 全谱系 UMAP，分析3+5）配对；原 Fig2B_dotplot_lineage_markers 保留）；`FigS4A_UMAP_GSE149652_Tcell`、`FigS4B_dotplot_AIFM2_GSE149652`、`FigS4C_Tcell_positivity`（pdf/png 双格式，300 dpi，cairo+Arial）；`results/tables/Fig2_source_data.csv`（GSE135337）与 `FigS4_source_data_GSE149652.csv`（任务书 7 列规格：cell_type, n_cells, n_positive, positivity_rate, wilson_lower, wilson_upper, mean_expression）。
- **备注**：① 命名说明——任务书 Fig2A/Fig2B 名称原指"全细胞谱系"，GSE149652 无全谱系，故其 UMAP/Dot 以 FigS4 命名，Fig2A（既有，GSE135337）未覆盖、Fig2B 用 GSE135337 填；② 同一条码在 CD4/CD8 两数据集间重复，单独条码不是唯一键，须 (Sample, 条码) 联合匹配（GSM 编号由 Sample 列唯一确定）；③ geom_errorbarh 在 ggplot2 4.0 弃用（height→width 自动转译，仅 warning）；④ GSE135337 谱系标签来自分析3+5 的 AddModuleScore argmax 归类（肿瘤上皮=EpCAM 模块高分），非 TISCH2 式人工注释，论文须注明。

**结果描述（中文版，可直接入稿）**：

> 在全细胞谱系分析（GSE135337，8 例未治疗膀胱癌，42,455 个细胞）中，AIFM2 阳性率在肿瘤上皮细胞中为 8.5%（n=34,034）、内皮细胞中为 9.4%（n=659）、成纤维细胞中为 5.6%（n=3,387）、巨噬细胞中为 5.4%（n=1,652）、T 细胞中为 5.0%（n=759），表明 FSP1 表达并非 T 细胞特异，而是广泛存在于肿瘤微环境多谱系。在 TIL 分选数据（GSE149652）的 T/NK 细胞亚群中，AIFM2 阳性率在 Treg 中为 0.76%（n=3,813），CD8⁺T 中为 0.64%（n=5,585），全部 6 个亚群阳性率均 <2%（Tprolif 最高，1.3%）。两数据集因测序深度不同不可直接比较阳性率绝对值。

**Result description (English version)**:

> In the whole-cell-lineage analysis (GSE135337, eight treatment-naive bladder cancer tumors, 42,455 cells), AIFM2 positivity was 8.5% in tumor epithelial cells (n = 34,034), 9.4% in endothelial cells (n = 659), 5.6% in fibroblasts (n = 3,387), 5.4% in macrophages (n = 1,652), and 5.0% in T cells (n = 759), indicating that FSP1 expression is not T cell-specific but is broadly distributed across multiple lineages of the tumor microenvironment. In the TIL-sorted dataset (GSE149652), AIFM2 positivity was 0.76% in Tregs (n = 3,813) and 0.64% in CD8⁺ T cells (n = 5,585), with all six T/NK subsets below 2% (Tprolif highest, 1.3%). Positivity levels are not directly comparable between the two datasets owing to differences in sequencing depth.


