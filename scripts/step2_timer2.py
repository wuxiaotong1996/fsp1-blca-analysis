# -*- coding: utf-8 -*-
"""Step 2：FSP1 与免疫浸润的关联分析（TIMER2.0 网页功能的脚本等价实现）。

TIMER2.0 匿名 API 已下线（301 → TIMER3 需注册 token），本脚本采用
**本地 ssGSEA（单样本基因集富集分析，Barbie 2009）**：以内置的免疫细胞
经典 marker 基因签名为基因集，对每个肿瘤样本计算各免疫细胞浸润评分，
再与 FSP1 表达做 Spearman 相关分析。
TIMER2.0 网页版（compbio.cn/timer3）可作为人工交叉验证。

输出：
- figures/step2_correlation.png  FSP1 与各免疫细胞评分的相关系数条形图
- figures/step2_scatter.png      FSP1 与相关最强 4 类细胞评分的散点图
- data/immune_scores.csv、data/step2_correlation.csv、data/step2_summary.json
"""
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import config as C

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei']
plt.rcParams['axes.unicode_minus'] = False

# 免疫细胞经典 marker 基因签名（教科书级 lineage marker，参考 Newman 2015 LM22、
# Bindea 2013 与免疫细胞图谱惯例；名称按官方基因符号）
IMMUNE_SIGNATURES = {
    'CD8+ T 细胞': ['CD8A', 'CD8B', 'CD3D', 'CD3E'],
    'CD4+ T 细胞': ['CD4', 'IL7R', 'CD3G'],
    '调节性 T 细胞': ['FOXP3', 'CTLA4', 'IL2RA', 'TNFRSF18'],
    'Th1 细胞': ['TBX21', 'IFNG', 'IL12RB2', 'STAT4'],
    'Th2 细胞': ['GATA3', 'IL4', 'IL5', 'IL13'],
    'Th17 细胞': ['RORC', 'IL17A', 'IL22', 'CCR6'],
    'NK 细胞': ['NCAM1', 'KLRD1', 'KLRK1', 'NCR1', 'NKG7', 'GNLY'],
    '细胞毒性细胞': ['GZMB', 'PRF1', 'GZMA', 'FASLG'],
    '耗竭 T 细胞': ['PDCD1', 'LAG3', 'HAVCR2', 'TIGIT', 'TOX', 'ENTPD1'],
    'γδ T 细胞': ['TRGC1', 'TRGC2', 'TRDC', 'TRDV2'],
    'B 细胞': ['CD19', 'MS4A1', 'CD79A', 'CD79B', 'BLK'],
    '浆细胞': ['SDC1', 'XBP1', 'MZB1', 'TNFRSF17'],
    '单核细胞': ['CD14', 'LYZ', 'S100A8', 'S100A9', 'FCN1'],
    'M1 巨噬细胞': ['NOS2', 'IL1B', 'TNF', 'IL6', 'CXCL9', 'CXCL10'],
    'M2 巨噬细胞': ['MRC1', 'CD163', 'MSR1', 'ARG1', 'IL10'],
    '中性粒细胞': ['CSF3R', 'FCGR3B', 'S100A12', 'CEACAM8'],
    '树突状细胞': ['ITGAX', 'CD1C', 'CLEC9A', 'FCER1A'],
    '浆样树突细胞': ['LILRA4', 'IL3RA', 'CLEC4C', 'TLR7'],
    '肥大细胞': ['TPSAB1', 'TPSB2', 'CPA3', 'MS4A2', 'KIT'],
    '嗜酸性粒细胞': ['SIGLEC8', 'IL5RA', 'CCL11', 'EPX'],
    '内皮细胞': ['PECAM1', 'VWF', 'CDH5', 'ENG', 'KDR'],
    '成纤维细胞': ['COL1A1', 'COL1A2', 'COL3A1', 'ACTA2', 'FAP', 'PDGFRB', 'DCN'],
    '免疫检查点': ['CD274', 'PDCD1LG2', 'IDO1', 'LILRB2'],
}


def load_expr_tumor():
    f = C.RAW_DIR / 'TCGA.BLCA.HiSeqV2.gz'
    expr = pd.read_csv(f, sep='\t', compression='gzip', index_col=0)
    expr.index = expr.index.astype(str).str.upper()
    is_tumor = expr.columns.to_series().apply(
        lambda s: s.count('-') >= 3 and s.split('-')[3][:2] == '01')
    x = expr.loc[:, is_tumor]
    # 同一基因符号多行时取均值，避免重复行
    x = x.groupby(x.index).mean()
    print(f'[数据] 肿瘤样本 n={x.shape[1]}，基因 n={x.shape[0]}')
    return x


def ssgsea_scores(expr, signatures, alpha=0.25):
    """ssGSEA（Barbie 2009）标准实现：GSVA 2.6.12 的 .fastRndWalk 公式。

    对每个样本：ES = Σ_S w·pos / Σ_S w − (n(n+1)/2 − Σ_S pos) / (n−k)
    其中 w = int(平均秩)^α（GSVA 对平均秩取整），pos = n − 降序位置 + 1。
    """
    rank_df = expr.rank(axis=0, method='average')          # 每列（样本）内基因平均秩
    r_int = rank_df.astype(int)                            # GSVA: type(R) <- "integer"（截断）
    n_genes = rank_df.shape[0]
    order_idx = np.argsort(-r_int.values, axis=0, kind='stable')   # 每样本降序（行号）
    pos = np.empty((n_genes, rank_df.shape[1]), dtype=np.int32)    # 行 -> 1-based 降序位置
    pos[order_idx, np.arange(rank_df.shape[1])] = np.arange(1, n_genes + 1)[:, None]
    scores = {}
    for name, genes in signatures.items():
        g = [g.upper() for g in genes if g.upper() in rank_df.index]
        if not g:
            print(f'[跳过] {name}：签名基因均不在表达矩阵中')
            continue
        rows = [list(rank_df.index).index(gg) for gg in g]
        w = (r_int.values[rows, :] ** alpha).astype(float)         # |r|^α
        pos_asc = n_genes - pos[rows, :] + 1
        step_in = (w * pos_asc).sum(axis=0) / w.sum(axis=0)
        step_out = (n_genes * (n_genes + 1) / 2 - pos_asc.sum(axis=0)) / (n_genes - len(g))
        scores[name] = pd.Series(step_in - step_out, index=rank_df.columns)
    return pd.DataFrame(scores)


def main():
    expr = load_expr_tumor()
    fsp1 = next((g for g in (C.GENE_SYMBOL, C.GENE_ALIAS) if g in expr.index), None)
    if fsp1 is None:
        raise SystemExit('表达矩阵中未找到 FSP1/AIFM2')
    label = fsp1 if fsp1 == C.GENE_ALIAS else f'{C.GENE_ALIAS} ({fsp1})'
    print(f'[基因] 使用行 {fsp1}（显示 {label}）')
    x = expr.loc[fsp1].astype(float)

    scores = ssgsea_scores(expr, IMMUNE_SIGNATURES, C.SSGSEA_ALPHA)
    scores.to_csv(C.DATA_DIR / 'immune_scores.csv', encoding='utf-8-sig')
    print(f'[评分] 免疫细胞类别 n={scores.shape[1]}')

    rows = []
    for col in scores.columns:
        rho, p = stats.spearmanr(x, scores[col])
        rows.append({'cell_type': col, 'rho': rho, 'p_value': p})
    corr = pd.DataFrame(rows).sort_values('rho', key=abs, ascending=False)
    corr.to_csv(C.DATA_DIR / 'step2_correlation.csv', index=False, encoding='utf-8-sig')

    # 图1：相关系数条形图（按 rho 排序，着色=方向，标注显著性）
    fig, ax = plt.subplots(figsize=(5.2, 7.2), dpi=300)
    cc = corr.sort_values('rho')
    colors = ['#C44E52' if v > 0 else '#4C72B0' for v in cc['rho']]
    ax.barh(cc['cell_type'], cc['rho'], color=colors, alpha=0.85, height=0.65)
    for i, (_, r) in enumerate(cc.iterrows()):
        star = '***' if r.p_value < 0.001 else '**' if r.p_value < 0.01 else '*' if r.p_value < 0.05 else ''
        ax.text(r.rho + (0.01 if r.rho >= 0 else -0.01), i,
                f'{r.rho:.2f}{star}', va='center',
                ha='left' if r.rho >= 0 else 'right', fontsize=7.5)
    ax.axvline(0, color='k', lw=0.8)
    ax.set_xlabel(f'{label} 与免疫浸润评分的 Spearman 相关系数')
    ax.set_title(f'{label} 与 BLCA 免疫微环境浸润的关联（ssGSEA，n={len(x)}）', fontsize=10)
    ax.set_xlim(-1.05 * max(0.05, abs(cc['rho']).max()), 1.05 * max(0.05, abs(cc['rho']).max()))
    fig.tight_layout()
    out1 = C.FIGURES_DIR / 'step2_correlation.png'
    fig.savefig(out1, dpi=300); plt.close(fig)
    print(f'[图] {out1}')

    # 图2：相关最强 4 类散点
    top4 = corr.head(4)
    fig, axes = plt.subplots(1, 4, figsize=(13, 3.1), dpi=300, sharey=False)
    for ax, (_, r) in zip(axes, top4.iterrows()):
        ax.scatter(scores[r.cell_type], x, s=12, alpha=0.5, color='#2E5E8C', linewidths=0)
        k, b = np.polyfit(scores[r.cell_type], x, 1)
        xs = np.linspace(scores[r.cell_type].min(), scores[r.cell_type].max(), 50)
        ax.plot(xs, k * xs + b, color='#C44E52', lw=1.4)
        ax.set_xlabel(r.cell_type, fontsize=9)
        ax.set_ylabel(f'{label} 表达' if ax is axes[0] else '', fontsize=9)
        ax.set_title(f'ρ = {r.rho:.2f}, P = {r.p_value:.3g}', fontsize=8.5)
    fig.suptitle(f'{label} 表达与免疫浸润评分（相关最强 4 类）', fontsize=10)
    fig.tight_layout()
    out2 = C.FIGURES_DIR / 'step2_scatter.png'
    fig.savefig(out2, dpi=300); plt.close(fig)
    print(f'[图] {out2}')

    summary = {'n_samples': int(len(x)), 'n_cell_types': int(scores.shape[1]),
               'top_correlations': corr.head(5).to_dict('records')}
    with open(C.DATA_DIR / 'step2_summary.json', 'w', encoding='utf-8') as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    print('[Step2 完成] 前5:', [(r.cell_type, round(r.rho, 3), f'{r.p_value:.2e}')
                               for _, r in corr.head(5).iterrows()])


if __name__ == '__main__':
    main()
