# -*- coding: utf-8 -*-
"""Step 1 浏览器自动化版（browser-use 0.13）：控制 GEPIA2 网页完成表达+生存分析。

与 Xena 版 step1_gepia2.py 的区别：本脚本走真实网页 UI，目标是拿到 GEPIA2
原版截图与数值用于交叉验证。实测（2026-10-15）：匿名接口 /api/json.boxplot
与 /api/json.survival 虽已 404，但网页 UI 可用——箱线图（num(T)=404,
num(N)=28）与 KM 曲线（HR(high)=0.86、logrank p=0.3）均成功出图，方向与
Xena 版一致。产物归档于 output/data/gepia2_browser/。**任何失败都会自动
回退运行 Xena 版**，保证 Step1 始终有标准结果。

依赖：pip install browser-use（playwright chromium 已装）
需要环境变量 ANTHROPIC_API_KEY（browser-use 用 Claude 驱动浏览器决策，
每次运行会产生少量 API 费用）。运行：python step1_gepia2_browser.py
"""
import asyncio
import subprocess
import sys
from pathlib import Path

import config as C

LLM_MODEL = 'claude-haiku-5-5'   # 浏览器决策用便宜快速的模型
MAX_STEPS = 60                   # 步数上限，防在坏页面上空转烧 token

FIG_DIR = str(C.FIGURES_DIR).replace('\\', '/')
DATA_DIR = str(C.DATA_DIR).replace('\\', '/')

TASK = f"""Go to {C.GEPIA2_URL}#analysis

Step A - Expression Box Plot:
1. Click 'Expression DIY', then 'Box Plot'
2. Enter gene '{C.GENE_SYMBOL}' (do NOT use 'FSP1')
3. Select '{C.CANCER_CODE}' Tumor and '{C.CANCER_CODE}' Normal datasets
4. Set 'Matched Normal data' to 'Match TCGA normal and GTEx data'
5. Click 'Plot', wait up to 60 seconds for the chart
6. If the chart appears: save the chart (screenshot or PDF) via your file
   system tool — Windows D: paths are NOT writable from the browser — and
   note the p-value shown on the page
7. If NO chart appears (blank, error, '404'): STOP Step A. Do not retry more
   than 3 times. Note exactly what the page shows.

Step B - Survival:
8. Open 'Survival' analysis: gene '{C.GENE_SYMBOL}', dataset '{C.CANCER_CODE}',
   Group Cutoff 'Median', click 'Plot', wait up to 60 seconds
9. If a Kaplan-Meier curve appears: save the chart (screenshot or PDF) via your
   file system tool and note HR and logrank p-value
10. If NO curve appears: STOP. Note exactly what the page shows.

Finally: write a JSON file named 'gepia2_browser_result.json' (via your file
system tool) with fields: expression_p_value, survival_hr, survival_logrank_p,
ok (true/false), notes. Describe in notes what you actually saw. NEVER invent
numbers - if a value was not shown, set it to null.
"""


async def run_browser():
    from browser_use import Agent, Browser
    from browser_use.llm import ChatAnthropic
    browser = Browser(headless=True, enable_default_extensions=False)  # 默认扩展 CDN 在本网络不通，必须关掉
    try:
        llm = ChatAnthropic(model=LLM_MODEL,
                            thinking={'type': 'disabled'})  # 5.x 模型服务端默认开思考模式，与强制 tool_choice 冲突（400），必须显式关
        agent = Agent(task=TASK, llm=llm,
                      browser=browser, use_vision=True,
                      use_thinking=False,  # 输出 schema 也不带 thinking 字段
                      use_judge=False,
                      file_system_path=str(C.BASE_DIR))
        return await agent.run(max_steps=MAX_STEPS)
    finally:
        try:
            await browser.stop()
        except Exception:
            pass


def main():
    print('[浏览器版] browser-use 控制 GEPIA2 网页（headless）…')
    try:
        result = asyncio.run(run_browser())
        print('[浏览器版] 结束：done=%s successful=%s' % (result.is_done(), result.is_successful()))
        print('[浏览器版] 最终文本：')
        print(result.final_result())
        if result.errors():
            print('[浏览器版] 错误：', result.errors())
    except Exception as e:
        print(f'[浏览器版失败] {type(e).__name__}: {e}')
    # 无论浏览器版结果如何，都回退跑 Xena 版，保证 Step1 出标准结果
    print('[回退] 运行 Xena 版 step1_gepia2.py …')
    r = subprocess.run([sys.executable, str(Path(C.BASE_DIR) / 'step1_gepia2.py')],
                       cwd=str(C.BASE_DIR))
    sys.exit(r.returncode)


if __name__ == '__main__':
    main()
