# QC_report.md — FSP1×BLCA 补充分析最终质量检查

生成日期：2026-10-28 ｜ 检查范围：6 项任务书分析 + 纯度重评估（2b）+ scRNA 全谱系（6/6b）
结论：**静态核查全绿；2 项发现需处理（Fig2B 重名、Fig1A/B/D 未拆分）；数字一致性核查因 manuscript_draft.md 缺失而阻塞；7 脚本从头重跑验证完成（见 §5）。**

---

## 1. 图表完整性

### 1.1 稿件引用名 vs 实际文件

| 稿件引用 | 实际文件（results/figures/） | PDF (KB) | PNG (KB) | >10KB | 状态 |
|---|---|---|---|---|---|
| Fig1A | — | — | — | — | ❌ 缺失（fig1_composite.png 尚未拆分） |
| Fig1B | — | — | — | — | ❌ 缺失（同上） |
| Fig1C | Fig1C_partial_correlation | 50.5 | 129.6 | ✓✓ | ✅ |
| Fig1D | — | — | — | — | ❌ 缺失（同上） |
| Fig2A | Fig2A_UMAP_all_celltypes | 1369.8 | 308.2 | ✓✓ | ✅ |
| Fig2B | Fig2B_dotplot_AIFM2 | 43.6 | 98.7 | ✓✓ | ✅（⚠️ 重名，见发现 #1） |
| Fig2C | — | — | — | — | ❌ 缺失（内容即 FigS4C 的 T 亚群阳性率柱状图） |
| FigS1 | FigS1_downsampling | 46.8 | 147.1 | ✓✓ | ✅ |

其余存在且双格式齐全：FigS2_partial_corr（22 类旧版）、FigS4A_UMAP_GSE149652_Tcell（594.7/1181）、FigS4B_dotplot_AIFM2_GSE149652（43.1/96.9）、FigS4C_Tcell_positivity（41.8/78.1）——全部 >10KB。

旧项目图（仅 PNG、无 PDF，早于 300dpi 双格式约定）：fig1_composite、step1_expression、step1_survival、step2_correlation、step2_scatter、step3_alteration、step3_lollipop、step3_cna、gepia2_expression、gepia2_survival。

### 1.2 对应关系说明（与任务书命名的差异）

- 任务书 Fig2A/Fig2B 原指"全细胞谱系"；GSE149652 是 TIL 分选数据（无全谱系），故其 UMAP/Dot/柱状图命名为 **FigS4A/B/C**，全谱系版为 Fig2A（分析3+5）+ Fig2B（分析6b）。
- 若稿件正文引用了"图2C"，应改为引用 FigS4C。

---

## 2. 数据一致性

### 2.1 FigS4_source_data_GSE149652.csv vs TISCH2 注释表（15,538 行）——✅ 全一致

| 类型 | TISCH2 注释计数 | CSV 计数 | 一致 |
|---|---|---|---|
| Treg | 3,813 | 3,813 | ✓ |
| CD8T | 5,585 | 5,585 | ✓ |
| CD4Tconv | 3,828 | 3,828 | ✓ |
| NK | 534 | 534 | ✓ |
| CD8Tex | 1,539 | 1,539 | ✓ |
| Tprolif | 239 | 239 | ✓ |

（QC 脚本初检时 CD8T/CD4Tconv 因 CSV 显示名改为 CD8_T/CD4_Tconv 显示"不匹配"，数值一致，仅为命名差异。）

### 2.2 Fig2_source_data.csv vs GSE135337_seurat.rds 谱系计数——✅ 全一致（总 42,455）

Tumor 34,034 / Macrophage 1,652 / Endothelial 659 / Fibroblast 3,387 / T 759 / Other 1,964，6 类全对。

### 2.3 TableS1_ssGSEA_full.csv —— ✅ 24 行

24 类免疫细胞（表头 + 24 行，n=342）。另有分析4 详细版 `TableS1_ssGSEA_full_detail.csv`（含 n_genes_signature 等扩展列）并存，引用时二选一（见发现 #5）。

### 2.4 确定性锚点复核——✅ 全部复现

- Treg 29/3,813（0.7606%）、CD8T 36/5,585（0.6446%）——与旧项目已知值精确一致；
- GSE135337 Tumor 2,901/34,034（8.524%）、Macrophage 90/1,652（5.448%）。

---

## 3. 数字一致性（analysis_log vs 正文）—— ⛔ 阻塞

**manuscript_draft.md 在磁盘上不存在**。已全盘搜索：工作区（根目录仅 README/analysis_log/requirements）、上级目录、papers/、C/D 盘桌面、文档、下载、当前 VSCode 打开文件——均无此文件。

待文件到位后执行以下比对（锚点数字清单已备好）：
- 摘要：降采样中位 P=0.559 / P<0.05 占 3.2%；CD8⁺T 校正后 ρ=0.080/P=0.138（原 0.099/P=0.046 必须撤下）；全谱系阳性率 8.5%/5.4%（不得写"所有类型均<2%"）；
- 3.1：原始 P=0.793、第 76 百分位、log2FC 中位 0.060；
- 3.2：FSP1×Immune/Stromal ρ=0.051/0.025（P=0.304/0.617）；24 类校正后 NK CD56bright ρ=0.233/P=1.4×10⁻⁵、Cytotoxic 0.162/P=0.003、pDC 0.151/P=0.005、Tgd 0.148/P=0.006、CD8 T 0.147/P=0.006、Neutrophils −0.168/P=0.002、Tcm −0.197/P=2.6×10⁻⁴；
- 3.4：GSE135337 42,455 细胞、6 谱系计数、阳性率 5–9%；GSE149652 6 亚群均 <2%。

---

## 4. 随机种子核查

| 脚本 | set.seed(20261028) | 随机成分 | 结论 |
|---|---|---|---|
| analysis1_downsampling.R（v1，已被 1b 取代） | ✓ | 降采样 | ✅ |
| analysis1b_downsampling.R | ✓ | 降采样 | ✅ |
| analysis2_estimate_ppcor.R | 无 | **无随机成分**（ppcor/手工 ssGSEA 均确定性） | ✅ 无需种子 |
| analysis2b_purity_24class.R | ✓ | GSVA ssGSEA | ✅ |
| analysis4_ssgsea_24class.R | ✓ | GSVA ssGSEA | ✅ |
| analysis35_scRNA_lineages.R | ✓ | PCA/UMAP | ✅ |
| analysis6_gse149652_tcell.R | ✓ | 无（按约定保留） | ✅ |
| analysis6b_gse135337_aifm2.R | ✓ | 无（按约定保留） | ✅ |
| prep_inputs_clean.R | 无 | 无随机 | ✅ 无需种子 |
| step1_gepia2.py（旧项目） | default_rng(7) | 箱线图抖动（仅可视化） | ⚠️ 旧种子早于约定 |
| step2_timer2.py / step3_cbioportal.py（旧项目） | 无随机调用 | — | ✅ |

结论：全部涉及随机的补充分析脚本均设 20261028。

---

## 5. 代码可复现

### 5.1 环境（实测）

- R 4.6.1（2026-06-24 ucrt）；Seurat 5.6.0；GSVA 2.6.14；ppcor 1.1；data.table 1.18.6.1；ggplot2 4.0.3；GSVAdata **未安装**（无需——24 类签名来自本地 bindea24_signatures.csv）。
- 安装脚本：`scripts/install_packages.R`（aliyun CRAN + tuna Bioconductor 镜像，wininet 下载栈）。

### 5.2 从头重跑验证——✅ 7/7 通过（QC_rerun.log，22:10:18–22:16:04，约 5.8 分钟，exit 0，无 FAILED/Error 行）

| 脚本 | 结果 | 关键证据（日志） |
|---|---|---|
| analysis1b_downsampling.R | ✅ | [done] analysis1b（FigS1 + downsampling_results.csv 重出） |
| analysis2_estimate_ppcor.R | ✅ | [done] analysis2（输入断言 407×20,530 过） |
| analysis2b_purity_24class.R | ✅ | [done] analysis2b |
| analysis4_ssgsea_24class.R | ✅ | [done] analysis4 |
| analysis35_scRNA_lineages.R | ✅ | [done] rds + lineage counts csv；Fig2A/Fig2B_dotplot_lineage_markers 重出 |
| analysis6_gse149652_tcell.R | ✅ | join 15,538/15,538；UMAP max\|diff\|=0；锚点 Treg 29/3,813、CD8T 36/5,585 复现；FigS4A/B/C 重出 |
| analysis6b_gse135337_aifm2.R | ✅ | 42,455 细胞、8 层拼接无缺失；Fig2B_dotplot_AIFM2 + Fig2_source_data.csv 重出 |

- 跳过说明：analysis1 v1 被 1b 取代（其 FigS1 输出已被 1b 覆盖，重跑反而回退）；prep_inputs_clean.R 为数据派生工具（重跑会覆盖用户自有版输入文件，未执行）。
- 重跑后图件与 CSV 已重新生成，尺寸复核见 §1。

### 5.3 Python 旧脚本（step1–3）

- 依赖（requirements.txt）：requests≥2.28、pandas≥2.0、numpy≥1.24、matplotlib≥3.7、seaborn≥0.13、scipy≥1.10（Windows + Python 3.12 已验证，2026-10-15）。
- **未重跑**：需 Xena/cBioPortal 联网下载，且 GEPIA2/TIMER API 已下线（只能浏览器自动化）——属旧项目范围，2026-10-15 已整体复核；与本次补充分析无数据流依赖。

---

## 6. 投稿前待办清单

- [ ] **参考文献核实（用户手动，两篇已预检索命中，请对照稿件引文逐字段核对）**：
  - **Hu J 2026** → 预检索命中：Hu J, Wang Y, Du Y, et al. SESN2 promotes ferroptosis in bladder cancer by negatively regulating the NRF2/FSP1 signaling axis. *Transl Androl Urol* 2026;15(9). DOI: **10.21037/tau-2026-0472**（PMID 42866843，2026-08-28 在线）。核对点：题目/卷期/页码与稿件引用一致。
  - **Sung H 2026** → 预检索命中：Sung H, Filho AM, Laversanne M, et al. Global cancer statistics 2024: GLOBOCAN estimates of incidence and mortality worldwide for 34 cancers in 186 countries. *CA Cancer J Clin* 2026;76(4):e70090. DOI: **10.3322/caac.70090**（PMID 42417444）。核对点：文章内容为 **2024 年数据、2026 年刊出**——若稿件引作"Global cancer statistics 2026"则为题目错误。
- [ ] **Acknowledgements**：待作者填写（致谢对象、基金号归属机构）。
- [ ] **Funding**：待作者填写（基金名称+编号）。
- [ ] **Data availability**：原始数据均为公开源（Xena/TCGA、GEO GSE135337、GSE149652、cBioPortal）；脚本+结果表可打包至 GitHub 并生成 Zenodo DOI——待建仓后填 DOI。
- [ ] **期刊格式**：请提供目标期刊名 → 我对照检查（字数限制、图注格式、参考文献样式、补充材料命名规范、Fig/Table 顺序）。

---

## 7. 发现汇总（需处理项）

1. **Fig2B 重名**：`Fig2B_dotplot_AIFM2.pdf`（GSE135337，正文用）与 `Fig2B_dotplot_lineage_markers.pdf`（分析3+5 旧版）并存——投稿只交前者，后者删除或改名（如 FigS5）。
2. **Fig1A/B/D 缺失**：fig1_composite.png（旧四合一拼图）尚未拆分为 Fig1A–D；若正文按 A–D 引用，需拆分重建（README 待办已列）。
3. **Fig2C 缺失**：任务书"图2C"为可选，实际内容（T 亚群阳性率柱状图）在 FigS4C；正文引用需改。
4. **manuscript_draft.md 缺失**：§3 数字一致性核查阻塞；文件到位后立即可补。
5. **TableS1 两版本并存**：任务书规格 6 列版（`TableS1_ssGSEA_full.csv`）与分析4 详细版（`_detail.csv`）同名不同列——投稿引用需明确指向其中一个。
6. **旧项目图无 PDF**：step1–3 图仅 PNG（早于双格式约定）；若作为正文图需重新导出 PDF。
7. **旧 Python 脚本种子为 7**（早于 20261028 约定；仅影响箱线图抖动，不影响任何统计量）。
8. **FigS2（22 类旧图）定位待确认**：正文 3.2 已改用 24 类体系（Fig1C），FigS2 是否仍引用需定；若引用，两套签名体系不可互比（analysis_log 已注明）。
