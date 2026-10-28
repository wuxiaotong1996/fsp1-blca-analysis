# -*- coding: utf-8 -*-
"""一键运行全部步骤：Step1 表达+生存 → Step2 免疫浸润 → Step3 遗传变异。

单步失败不中断，结束后打印总览；任一步失败返回非零退出码。
用法：python run_all.py   （或 python run_all.py 1 只跑第 1 步）
"""
import subprocess
import sys
import os

BASE = os.path.dirname(os.path.abspath(__file__))  # 项目根目录（本脚本所在目录）

steps = [
    ("Step 1: GEPIA2 表达与生存分析", "step1_gepia2.py"),
    ("Step 2: TIMER2.0 免疫浸润分析", "step2_timer2.py"),
    ("Step 3: cBioPortal 遗传变异分析", "step3_cbioportal.py"),
]

# 支持只跑某一步：python run_all.py 1（或 2 / 3）
only = [int(x) for x in sys.argv[1:] if x.isdigit()]

print("=" * 60)
print("FSP1 膀胱癌生物信息学分析 — 方案一")
print("=" * 60)

results = {}
for i, (name, script) in enumerate(steps, 1):
    if only and i not in only:
        print(f"\n  – {name}（跳过）")
        continue
    print(f"\n{'─' * 60}")
    print(f"▶ {name}")
    print(f"{'─' * 60}")
    # 不捕获输出：网络下载期间需要实时进度（capture_output=True 会全程无输出，
    # 像卡死一样；且 Windows 管道按 cp936 编码，特殊符号易乱码）
    result = subprocess.run([sys.executable, script], cwd=BASE)
    ok = result.returncode == 0
    results[i] = ok
    print(f"◀ {'✅' if ok else '❌'} {name} 退出码 {result.returncode}")
    if not ok:
        print(f"请修复错误后重新运行此步骤（可只跑这一步：python run_all.py {i}）")

print(f"\n{'=' * 60}")
print("总览：")
for i, (name, _) in enumerate(steps, 1):
    if i in results:
        print(f"  {'✅' if results[i] else '❌'} {name}")
    else:
        print(f"  – {name}（跳过）")
print("全部完成。请检查 output/figures/ 和 output/data/ 目录。")
print("=" * 60)
sys.exit(0 if results and all(results.values()) else 1)  # 任一步失败则非零退出
