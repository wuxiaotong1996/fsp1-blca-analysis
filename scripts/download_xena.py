# -*- coding: utf-8 -*-
"""Xena 数据并行分段下载工具：tcga.xenahubs.net 单连接限速（~18 KB/s），
用 8 线程 HTTP Range 并行下载后拼接，可将 22 MB 矩阵的下载时间从
约 20 分钟压缩到数分钟。下载结果存 raw/ 缓存，step1 检测到即跳过重复下载。"""
import os
import sys
from concurrent.futures import ThreadPoolExecutor

import requests

import config as C

N_PARTS = 8


def fetch_part(url, start, end, path, idx):
    headers = {'Range': f'bytes={start}-{end}', 'User-Agent': C.USER_AGENT}
    with requests.get(url, headers=headers, timeout=C.REQUEST_TIMEOUT, stream=True) as r:
        r.raise_for_status()
        part = str(path) + f'.part{idx}'
        with open(part, 'wb') as f:
            for chunk in r.iter_content(1 << 16):
                f.write(chunk)
        return part


def parallel_download(url, dest):
    with requests.get(url, stream=True, timeout=C.REQUEST_TIMEOUT,
                      headers={'User-Agent': C.USER_AGENT}) as r:
        r.raise_for_status()
        total = int(r.headers['Content-Length'])
    step = (total + N_PARTS - 1) // N_PARTS
    ranges = [(i * step, min((i + 1) * step - 1, total - 1)) for i in range(N_PARTS)]
    print(f'[并行下载] {url} 共 {total / 1e6:.1f} MB，{N_PARTS} 线程')
    with ThreadPoolExecutor(N_PARTS) as pool:
        parts = list(pool.map(lambda t: fetch_part(url, t[0], t[1], dest, t[2]),
                              [(s, e, i) for i, (s, e) in enumerate(ranges)]))
    with open(dest, 'wb') as out:
        for p in parts:
            with open(p, 'rb') as f:
                out.write(f.read())
            os.remove(p)
    print(f'[完成] {dest}（{os.path.getsize(dest) / 1e6:.1f} MB）')


def main():
    # 期望大小（来自服务器 Content-Length，精确校验避免残缺文件误判为缓存）
    EXPECTED = {
        C.EXPR_URL: 22827197,
        C.CLIN_URL: 485735,
    }
    jobs = [
        (C.EXPR_URL, C.RAW_DIR / 'TCGA.BLCA.HiSeqV2.gz'),
        (C.CLIN_URL, C.RAW_DIR / 'TCGA.BLCA.clinicalMatrix'),
    ]
    for url, dest in jobs:
        if dest.exists() and dest.stat().st_size == EXPECTED.get(url):
            print(f'[跳过] {dest} 已存在且完整')
            continue
        if dest.exists():
            print(f'[清理] 残缺文件 {dest}（{dest.stat().st_size} B）')
            dest.unlink()
        for attempt in range(1, C.RETRIES + 1):
            try:
                parallel_download(url, dest)
                break
            except Exception as e:
                print(f'[重试 {attempt}/{C.RETRIES}] {e}')
        else:
            sys.exit(1)


if __name__ == '__main__':
    main()
