# -*- coding: utf-8 -*-
"""Step 3：FSP1 在 BLCA 中的遗传变异分析（cBioPortal 官方 REST API，无需登录）。

数据：TCGA BLCA PanCancer Atlas（blca_tcga_pan_can_atlas_2018）
- 突变：mutation profile（FSP1/AIFM2，Entrez 84883）
- 拷贝数变异：GISTIC 离散值（-2 深缺失 … +2 扩增）

输出：
- figures/step3_alteration.png  变异频率汇总（突变/扩增/缺失/多事件）
- figures/step3_lollipop.png    突变位点棒棒糖图（蛋白位置 × 次数）
- figures/step3_cna.png         拷贝数状态分布
- data/fsp1_mutations.csv、data/fsp1_cna.csv、data/step3_summary.json
"""
import json
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests

import config as C

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei']
plt.rcParams['axes.unicode_minus'] = False


def api_get(path, **params):
    url = f'{C.CBIOPORTAL_API}{path}'
    for attempt in range(1, C.RETRIES + 1):
        try:
            r = requests.get(url, params=params, timeout=C.REQUEST_TIMEOUT,
                             headers={'User-Agent': C.USER_AGENT})
            r.raise_for_status()
            return r.json()
        except Exception as e:
            print(f'[重试 {attempt}/{C.RETRIES}] {e}')
            time.sleep(3)
    raise RuntimeError(f'cBioPortal 请求失败：{path}')


def api_post(path, json_body):
    """cBioPortal * /fetch 端点：POST + JSON body。"""
    url = f'{C.CBIOPORTAL_API}{path}'
    for attempt in range(1, C.RETRIES + 1):
        try:
            r = requests.post(url, json=json_body, timeout=C.REQUEST_TIMEOUT,
                              headers={'User-Agent': C.USER_AGENT})
            r.raise_for_status()
            return r.json()
        except Exception as e:
            print(f'[重试 {attempt}/{C.RETRIES}] {e}')
            time.sleep(3)
    raise RuntimeError(f'cBioPortal 请求失败：{path}')


def get_profiles():
    """定位 mutations 与 gistic（DISCRETE）分子谱 ID。"""
    data = api_get(f'/studies/{C.STUDY_ID}/molecular-profiles')
    mut = next((p['molecularProfileId'] for p in data
                if p['molecularAlterationType'] == 'MUTATION_EXTENDED'), None)
    cna = next((p['molecularProfileId'] for p in data
                if p['molecularAlterationType'] == 'COPY_NUMBER_ALTERATION'
                and p.get('datatype') == 'DISCRETE'), None)
    if not mut or not cna:
        raise SystemExit(f'未找到 mutation/gistic 分子谱')
    print(f'[分子谱] mutation={mut}，gistic={cna}')
    return mut, cna


def get_sample_list(category='all_cases_in_study'):
    """取研究内全体样本列表 ID（mutations 端点必需 sampleListId）。"""
    lists = api_get(f'/studies/{C.STUDY_ID}/sample-lists')
    hit = next((l['sampleListId'] for l in lists if l.get('category') == category), None)
    if hit is None and lists:
        hit = lists[0]['sampleListId']
    print(f'[样本列表] {hit}')
    return hit


def fetch_mutations(profile, sample_list):
    data = api_get(f'/molecular-profiles/{profile}/mutations',
                   entrezGeneId=C.ENTREZ_ID, projection='DETAILED',
                   sampleListId=sample_list, pageSize=500, pageNumber=0)
    if not data:
        print('[突变] 该队列中未检索到 FSP1 突变记录')
        return pd.DataFrame()
    df = pd.DataFrame(data)
    cols = ['sampleId', 'patientId', 'proteinChange', 'mutationType',
            'variantType', 'proteinPosStart', 'proteinPosEnd']
    df = df[[c for c in cols if c in df.columns]].copy()
    df.to_csv(C.DATA_DIR / 'fsp1_mutations.csv', index=False, encoding='utf-8-sig')
    print(f'[突变] 记录数={len(df)}，样本数={df["sampleId"].nunique()}')
    return df


def fetch_cna(profile):
    data = api_post(f'/molecular-profiles/{profile}/molecular-data/fetch',
                    {'entrezGeneIds': [C.ENTREZ_ID],
                     'sampleListId': f'{C.STUDY_ID}_all'})
    if not data:
        print('[CNA] 未检索到拷贝数数据')
        return pd.Series(dtype=float)
    s = pd.Series({d['sampleId']: int(d['value']) for d in data}, name='cna')
    s.to_csv(C.DATA_DIR / 'fsp1_cna.csv', encoding='utf-8-sig')
    print(f'[CNA] 样本数={len(s)}，扩增(2)={int((s == 2).sum())}，'
          f'缺失(-2)={int((s == -2).sum())}，浅层缺失(-1)={int((s == -1).sum())}')
    return s


def fig_alteration(mut_df, cna):
    labels = ['FSP1 变异频率']
    amp = int((cna == 2).sum()); deep_del = int((cna == -2).sum())
    shallow = int((cna == -1).sum()); gain = int((cna == 1).sum())
    mut_n = 0 if mut_df.empty else int(mut_df['sampleId'].nunique())
    total = len(cna)
    data = {
        '突变': mut_n / total * 100,
        '扩增': amp / total * 100,
        '获得': gain / total * 100,
        '浅层缺失': shallow / total * 100,
        '深缺失': deep_del / total * 100,
    }
    fig, ax = plt.subplots(figsize=(5.6, 3.4), dpi=300)
    colors = ['#3F7F4D', '#C44E52', '#E3A19E', '#6A9EC4', '#2E5E8C']
    bottom = 0.0
    for (name, pct), color in zip(data.items(), colors):
        ax.bar(0, pct, bottom=bottom, color=color, width=0.5, label=f'{name} {pct:.1f}%')
        if pct > 1:
            ax.text(0, bottom + pct / 2, f'{pct:.1f}%', ha='center', va='center',
                    fontsize=8, color='white' if pct > 8 else 'black')
        bottom += pct
    ax.set_xlim(-0.45, 0.45)
    ax.set_xticks([0]); ax.set_xticklabels([labels[0]])
    ax.set_ylabel('样本占比（%）')
    ax.set_title(f'FSP1 在{C.CANCER_NAME}中的遗传变异谱\n（cBioPortal TCGA PanCancer Atlas，n={total}）', fontsize=10)
    ax.legend(fontsize=8, loc='upper left', bbox_to_anchor=(1.0, 1.0))
    fig.tight_layout()
    out = C.FIGURES_DIR / 'step3_alteration.png'
    fig.savefig(out, dpi=300, bbox_inches='tight'); plt.close(fig)
    print(f'[图] {out}')
    return {'total': total, 'mutated': mut_n, 'amplified': amp, 'gain': gain,
            'shallow_del': shallow, 'deep_del': deep_del}


def fig_lollipop(mut_df):
    if mut_df.empty or 'proteinPosStart' not in mut_df.columns:
        print('[lollipop] 无突变记录，跳过')
        return 0
    pos = pd.to_numeric(mut_df['proteinPosStart'], errors='coerce').dropna()
    if len(pos) == 0:
        print('[lollipop] 无有效蛋白位置，跳过')
        return 0
    counts = pos.value_counts().sort_index()
    fig, ax = plt.subplots(figsize=(6.8, 3.0), dpi=300)
    ax.stem(counts.index, counts.values, linefmt='#999999', markerfmt='o',
            basefmt='k-')
    ax.set_xlabel('蛋白位置（氨基酸）')
    ax.set_ylabel('突变次数')
    ax.set_title(f'FSP1 突变位点分布（TCGA BLCA，n={len(pos)} 次突变）', fontsize=10)
    fig.tight_layout()
    out = C.FIGURES_DIR / 'step3_lollipop.png'
    fig.savefig(out, dpi=300); plt.close(fig)
    print(f'[图] {out}')
    return int(len(pos))


def fig_cna(cna):
    counts = {v: int((cna == v).sum()) for v in (-2, -1, 0, 1, 2)}
    names = {-2: '深缺失', -1: '浅层缺失', 0: '二倍体', 1: '获得', 2: '扩增'}
    colors = ['#2E5E8C', '#6A9EC4', '#BBBBBB', '#E3A19E', '#C44E52']
    fig, ax = plt.subplots(figsize=(4.6, 3.4), dpi=300)
    ax.bar([names[v] for v in (-2, -1, 0, 1, 2)],
           [counts[v] for v in (-2, -1, 0, 1, 2)], color=colors, alpha=0.85)
    for i, v in enumerate((-2, -1, 0, 1, 2)):
        ax.text(i, counts[v] + 2, str(counts[v]), ha='center', fontsize=9)
    ax.set_ylabel('样本数')
    ax.set_title('FSP1 拷贝数状态分布（GISTIC）', fontsize=10)
    fig.tight_layout()
    out = C.FIGURES_DIR / 'step3_cna.png'
    fig.savefig(out, dpi=300); plt.close(fig)
    print(f'[图] {out}')


def main():
    mut_prof, cna_prof = get_profiles()
    all_list = get_sample_list()
    mut_df = fetch_mutations(mut_prof, all_list)
    cna = fetch_cna(cna_prof)
    summary = {'gene': 'FSP1', 'entrez': C.ENTREZ_ID, 'study': C.STUDY_ID,
               'alteration': fig_alteration(mut_df, cna),
               'mutation_sites': fig_lollipop(mut_df)}
    if len(cna):
        fig_cna(cna)
        summary['cna_counts'] = {str(k): int((cna == k).sum()) for k in (-2, -1, 0, 1, 2)}
    with open(C.DATA_DIR / 'step3_summary.json', 'w', encoding='utf-8') as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    print('[Step3 完成]', json.dumps(summary, ensure_ascii=False))


if __name__ == '__main__':
    main()
