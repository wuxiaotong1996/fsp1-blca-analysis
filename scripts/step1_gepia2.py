# -*- coding: utf-8 -*-
"""Step 1：FSP1 在 BLCA 中的表达差异 + 生存分析（GEPIA2 网页功能的脚本等价实现）。

数据源：UCSC Xena TCGA hub（GEPIA2 底层数据同为 TCGA）。
- 表达：TCGA.BLCA.sampleMap/HiSeqV2.gz（log2(x+1) RSEM）
- 生存：BLCA_clinicalMatrix（OS 等）
输出：
- figures/step1_expression.png  肿瘤 vs 正常箱线图（Mann-Whitney U 检验）
- figures/step1_survival.png    Kaplan-Meier 曲线（中位数分组，log-rank + Cox HR）
- data/fsp1_expression.csv、data/blca_survival.csv、data/step1_summary.json
"""
import json
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize

import config as C

# 中文字体（Windows 自带）
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei']
plt.rcParams['axes.unicode_minus'] = False


# ---------- 下载 ----------
def download(url, dest, desc):
    import requests
    if dest.exists() and dest.stat().st_size > 1000:
        print(f'[缓存] {desc} 已存在：{dest}')
        return
    print(f'[下载] {desc} …')
    for attempt in range(1, C.RETRIES + 1):
        try:
            with requests.get(url, stream=True, timeout=C.REQUEST_TIMEOUT,
                              headers={'User-Agent': C.USER_AGENT}) as r:
                r.raise_for_status()
                with open(dest, 'wb') as f:
                    for chunk in r.iter_content(1 << 16):
                        f.write(chunk)
            print(f'[完成] {dest}（{dest.stat().st_size / 1e6:.1f} MB）')
            return
        except Exception as e:
            print(f'[重试 {attempt}/{C.RETRIES}] {e}')
            time.sleep(3)
    raise RuntimeError(f'下载失败：{url}')


# ---------- 表达 ----------
def load_expr():
    """返回 (Series FSP1 表达, 样本类型 Series)。"""
    f = C.RAW_DIR / 'TCGA.BLCA.HiSeqV2.gz'
    download(C.EXPR_URL, f, 'BLCA 表达矩阵 HiSeqV2')
    expr = pd.read_csv(f, sep='\t', compression='gzip', index_col=0)
    expr.index = expr.index.astype(str).str.upper()
    hits = [g for g in (C.GENE_SYMBOL, C.GENE_ALIAS) if g in expr.index]
    if not hits:
        raise SystemExit('表达矩阵中未找到 FSP1/AIFM2 行')
    gene = hits[0]
    x = expr.loc[gene]
    # 样本类型：TCGA-XX-XXXX-XX（第 4 段前两位：01=原发肿瘤，11=正常）
    btype = x.index.to_series().apply(lambda s: s.split('-')[3][:2] if s.count('-') >= 3 else '??')
    is_tumor = btype.str.startswith('01')
    is_normal = btype == '11'
    print(f'[表达] 基因行={gene}，肿瘤 n={is_tumor.sum()}，正常 n={is_normal.sum()}，'
          f'其他 n={(~is_tumor & ~is_normal).sum()}')
    return x, is_tumor, is_normal, gene


def fig_expression(x, is_tumor, is_normal, gene, label):
    tumor = x[is_tumor].astype(float).dropna()
    normal = x[is_normal].astype(float).dropna()
    u, p = stats.mannwhitneyu(tumor, normal, alternative='two-sided')
    fig, ax = plt.subplots(figsize=(4.2, 4.6), dpi=300)
    bp = ax.boxplot([normal, tumor], tick_labels=['正常组织', '肿瘤组织'],
                    patch_artist=True, widths=0.55, showfliers=False)
    for patch, color in zip(bp['boxes'], ['#4C72B0', '#C44E52']):
        patch.set_facecolor(color); patch.set_alpha(0.75)
    for i, vals in enumerate((normal, tumor)):
        jitter = np.random.default_rng(7).normal(0, 0.03, len(vals))
        ax.scatter(i + 1 + jitter, vals, s=7, alpha=0.4, color='#333333', linewidths=0)
    ymax = float(x.max())
    ax.plot([1, 1, 2, 2], [ymax, ymax + 0.25, ymax + 0.25, ymax], lw=1, color='k')
    ax.text(1.5, ymax + 0.28, f'P = {p:.3g}', ha='center', fontsize=10)
    ax.set_ylabel(f'{label} 表达量（log2 RSEM+1）')
    ax.set_title(f'{label} 在{C.CANCER_NAME}中的表达\n（TCGA，n={len(normal)}+{len(tumor)}）', fontsize=10)
    ax.set_ylim(float(x.min()) - 0.3, ymax + 0.9)
    fig.tight_layout()
    out = C.FIGURES_DIR / 'step1_expression.png'
    fig.savefig(out, dpi=300); plt.close(fig)
    df = pd.DataFrame({'sample': list(normal.index) + list(tumor.index),
                       'group': ['Normal'] * len(normal) + ['Tumor'] * len(tumor),
                       'expression': list(normal) + list(tumor)})
    df.to_csv(C.DATA_DIR / 'fsp1_expression.csv', index=False, encoding='utf-8-sig')
    print(f'[图] {out}')
    return {'n_tumor': int(len(tumor)), 'n_normal': int(len(normal)),
            'median_tumor': float(tumor.median()), 'median_normal': float(normal.median()),
            'mannwhitney_p': float(p)}


# ---------- 生存 ----------
def load_clin(expr_index):
    """临床矩阵 → 生存表（OS 优先，缺失回退 vital_status+days_to_* 构造）。"""
    f = C.RAW_DIR / 'TCGA.BLCA.clinicalMatrix'
    download(C.CLIN_URL, f, 'BLCA 临床矩阵')
    clin = pd.read_csv(f, sep='\t', index_col=0)
    clin.index = clin.index.astype(str)
    df = clin.reindex(expr_index)
    # Xena 矩阵通常直接提供 OS/OS.time；不同版本列名略有差异，做自适应
    if C.SURV_END in df.columns and f'{C.SURV_END}.time' in df.columns:
        surv = pd.DataFrame({'time': pd.to_numeric(df[f'{C.SURV_END}.time'], errors='coerce'),
                             'event': pd.to_numeric(df[C.SURV_END], errors='coerce')}, index=df.index)
    else:
        t = pd.Series(np.nan, index=df.index); e = pd.Series(np.nan, index=df.index)
        vital = df.get('vital_status', pd.Series('', index=df.index)).astype(str).str.upper()
        dead = vital.isin(['DECEASED', 'DEAD', '1'])
        t[dead] = pd.to_numeric(df['days_to_death'], errors='coerce').reindex(t.index)[dead]
        t[~dead] = pd.to_numeric(df['days_to_last_followup'], errors='coerce').reindex(t.index)[~dead]
        e[dead] = 1; e[~dead] = 0
        surv = pd.DataFrame({'time': t / 30.44, 'event': e})  # 天 → 月
    surv = surv.dropna(subset=['time'])
    surv = surv[surv['time'] > 0]
    print(f'[生存] 有效记录 n={len(surv)}，事件数={int(surv["event"].sum())}')
    return surv


def logrank(time_a, event_a, time_b, event_b):
    """两样本 log-rank 检验（手写）。"""
    t_all = np.sort(np.unique(np.concatenate([time_a, time_b])))
    o1 = e1 = 0.0
    for t in t_all:
        n1 = np.sum(time_a >= t); n2 = np.sum(time_b >= t)
        d1 = np.sum((time_a == t) & (event_a == 1)); d2 = np.sum((time_b == t) & (event_b == 1))
        d = d1 + d2
        if d == 0 or (n1 + n2) == 0:
            continue
        o1 += d1
        e1 += d * n1 / (n1 + n2)
    v1 = 0.0
    for t in t_all:
        n1 = np.sum(time_a >= t); n2 = np.sum(time_b >= t)
        d1 = np.sum((time_a == t) & (event_a == 1)); d2 = np.sum((time_b == t) & (event_b == 1))
        d = d1 + d2
        if d == 0 or n1 + n2 <= 1:
            continue
        v1 += d * (n1 + n2 - d) * n1 * n2 / ((n1 + n2) ** 2 * (n1 + n2 - 1))
    z = (o1 - e1) / np.sqrt(v1) if v1 > 0 else 0.0
    return 2 * (1 - stats.norm.cdf(abs(z))), float(o1), float(e1)


def cox_ph(time, event, group):
    """单变量 Cox 比例风险（部分似然 BFGS 拟合），返回 HR 与 95%CI。"""
    order = np.argsort(-time)  # 降序（事件时间从大到小便于风险集）
    t = time[order]; e = event[order]; g = group[order]
    mask = e == 1
    t_e, g_e = t[mask], g[mask]
    def nll(beta):
        b = beta[0]
        # 风险集：所有 t_j >= t_i 的个体
        s = 0.0
        for i in range(len(t_e)):
            at_risk = t >= t_e[i]
            s += b * g_e[i] - np.log(np.sum(np.exp(b * g[at_risk])))
        return -s
    res = minimize(nll, x0=[0.0], method='BFGS')
    beta = res.x[0]
    # Hessian（数值）
    eps = 1e-4
    h = (nll([beta + eps]) - 2 * nll([beta]) + nll([beta - eps])) / eps ** 2
    se = np.sqrt(1.0 / h) if h > 0 else np.nan
    hr = float(np.exp(beta))
    return hr, float(np.exp(beta - 1.96 * se)), float(np.exp(beta + 1.96 * se)), float(beta)


def km_curve(time, event):
    """Kaplan-Meier 估计（手写），返回 (x, y) 阶梯坐标。"""
    t = np.sort(np.unique(time))
    n = len(time); at_risk = n; surv = 1.0
    xs, ys = [0.0], [1.0]
    for ti in t:
        d = np.sum((time == ti) & (event == 1))
        if d > 0:
            surv *= (at_risk - d) / at_risk
        at_risk -= np.sum(time == ti)
        xs.append(ti); ys.append(surv)
        if at_risk == 0:
            break
    return np.array(xs), np.array(ys)


def fig_survival(surv, expr, is_tumor, label):
    common = surv.index.intersection(expr.index[is_tumor])
    e = expr.loc[common].astype(float)
    s = surv.loc[common]
    cut = e.median()
    hi = e >= cut
    t_hi, t_lo = s['time'][hi].values, s['time'][~hi].values
    ev_hi, ev_lo = s['event'][hi].values, s['event'][~hi].values
    p_lr, o1, e1 = logrank(t_hi, ev_hi, t_lo, ev_lo)
    hr, lo95, hi95, _ = cox_ph(s['time'].values, s['event'].values, hi.values.astype(float))

    fig, ax = plt.subplots(figsize=(5.4, 4.4), dpi=300)
    for tt, ee, lab, color in ((t_hi, ev_hi, f'高表达 (n={hi.sum()})', '#C44E52'),
                               (t_lo, ev_lo, f'低表达 (n={(~hi).sum()})', '#4C72B0')):
        xx, yy = km_curve(tt, ee)
        ax.step(xx, yy, where='post', color=color, lw=1.8, label=lab)
    ax.set_xlabel('总生存时间（月）')
    ax.set_ylabel('总生存概率')
    ax.set_title(f'{label} 表达与{C.CANCER_NAME}总生存（TCGA）')
    ax.text(0.03, 0.06,
            f'Log-rank P = {p_lr:.3g}\nHR = {hr:.2f}（95%CI {lo95:.2f}–{hi95:.2f}）',
            transform=ax.transAxes, fontsize=9, va='bottom',
            bbox=dict(boxstyle='round', fc='white', ec='gray', alpha=0.9))
    ax.legend(loc='upper right', fontsize=9)
    ax.set_ylim(-0.05, 1.05)
    fig.tight_layout()
    out = C.FIGURES_DIR / 'step1_survival.png'
    fig.savefig(out, dpi=300); plt.close(fig)
    s_out = s.copy(); s_out['fsp1_high'] = hi.astype(int)
    s_out.to_csv(C.DATA_DIR / 'blca_survival.csv', encoding='utf-8-sig')
    print(f'[图] {out}')
    return {'n_total': int(len(common)), 'events': int(s['event'].sum()),
            'median_cut': float(cut), 'logrank_p': float(p_lr),
            'hr': float(hr), 'hr_ci_low': float(lo95), 'hr_ci_high': float(hi95)}


def main():
    x, is_tumor, is_normal, gene = load_expr()
    label = gene if gene == C.GENE_ALIAS else f'{C.GENE_ALIAS} ({gene})'
    res = {'gene': gene, 'label': label,
           'expression': fig_expression(x, is_tumor, is_normal, gene, label)}
    surv = load_clin(x.index)
    res['survival'] = fig_survival(surv, x, is_tumor, label)
    with open(C.DATA_DIR / 'step1_summary.json', 'w', encoding='utf-8') as f:
        json.dump(res, f, ensure_ascii=False, indent=2)
    print('[Step1 完成]', json.dumps(res, ensure_ascii=False))


if __name__ == '__main__':
    main()
