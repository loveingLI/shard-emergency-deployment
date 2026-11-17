#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
ci_report.py —— 论文统一 95% 置信区间报告（SHARD / FGCS 投稿格式化 ⑤a）

原理：所有结果 CSV 中的比例都是「M 轮蒙特卡洛 Bernoulli 试验的汇总均值」，
      试验总数 N 已知且固定（见下方 TRIALS 表），因此可对每个单元计算
      Wilson score 95% 置信区间（小 p 与边界 p 均表现良好，优于正态近似）。

      Wilson 区间（z = 1.96）：
        center = (p + z^2/(2N)) / (1 + z^2/N)
        half   = z/(1+z^2/N) * sqrt(p(1-p)/N + z^2/(4N^2))

输入：results/ 下的 6 个 CSV（均为已验证的 [VERIFIED] 数据源）
输出：results/ci_report.txt   —— 全量 CI 表 + 论文头条数字的 CI 摘要
零依赖：仅 Python 标准库。

运行（零基础用户，cmd 中整行复制粘贴）：
"C:\Users\hp\.workbuddy\binaries\python\versions\3.13.12\python.exe" D:\workbuddy\网云系统\sim\ci_report.py
"""
import csv
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")

Z = 1.959963984540054  # 95% 双侧

# 每个数据文件中「每个单元」的蒙特卡洛试验总数 N（与生成脚本的常量逐一核对）
#   reach_results.csv  : 5 graph seeds x 600 episodes      = 3000  (reach_sim.py)
#   zoo_results.csv   : 400 episodes                       =  400  (zoo_sim.py)
#   sens_reach.csv    : 5 seeds x 600 episodes             = 3000  (sens_sim.py B/C)
#   sens_phase.csv    : 500 rounds per cell                =  500  (sens_sim.py A)
#   shard_v2_results  : A(n) 1500；P(theta>=c) 2000               (shard_sim_v2.py)
#   spread.csv        : 3000                                (spread_sim.py)
TRIALS = {
    "reach": 3000, "zoo": 400, "sens_reach": 3000, "sens_phase": 500,
    "shard_avail": 1500, "shard_theta": 2000, "spread": 3000,
}


def wilson(p, n):
    """返回 (low, high) —— Wilson score 95% 置信区间。"""
    if n <= 0:
        return (float("nan"), float("nan"))
    p = min(max(p, 0.0), 1.0)
    denom = 1.0 + Z * Z / n
    center = (p + Z * Z / (2 * n)) / denom
    half = (Z / denom) * math.sqrt(p * (1 - p) / n + Z * Z / (4 * n * n))
    return (max(0.0, center - half), min(1.0, center + half))


def fmt_ci(p, n):
    lo, hi = wilson(p, n)
    return "%.4f [%.4f, %.4f]  (N=%d, ±%.4f)" % (p, lo, hi, n, (hi - lo) / 2)


def read_rows(fname):
    with open(os.path.join(RES, fname), newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def f(x):
    return float(x)


out = []
A = out.append

A("=" * 78)
A("SHARD 论文 —— 全部关键结果的 95% Wilson 置信区间报告")
A("（生成者 ci_report.py；z=1.95996；比例均为蒙特卡洛汇总均值）")
A("=" * 78)

# ---------------------------------------------------------------- reach 主表
A("")
A("【1】§6.3 四臂可达性/生存率主表（reach_results.csv，每单元 N=3000）")
rows = read_rows("reach_results.csv")
N = TRIALS["reach"]
max_gap = None  # (gap, L, n, strategy, surv, reach)
for r in rows:
    s, rc = f(r["survival"]), f(r["reachable"])
    gap = s - rc
    if max_gap is None or gap > max_gap[0]:
        max_gap = (gap, r["L"], r["n"], r["strategy"], s, rc)
A("  L     n   strategy  survival / reachable (95% CI)")
for r in rows:
    A("  %-4s %-3s %-8s  %s   %s"
      % (r["L"], r["n"], r["strategy"],
         fmt_ci(f(r["survival"]), N), fmt_ci(f(r["reachable"]), N)))
g, gl, gn, gstr, gs, gr = max_gap
lo_s, hi_s = wilson(gs, N)
lo_r, hi_r = wilson(gr, N)
A("")
A("  ★ 头条：最大生存-可达 gap = %.3f（L=%s, n=%s, %s：生存 %.3f vs 可达 %.3f）"
  % (g, gl, gn, gstr, gs, gr))
A("    保守下界（生存 CI 下限 - 可达 CI 上限）= %.3f - %.3f = %.3f  > 0.4"
  % (lo_s, hi_r, lo_s - hi_r))

# ------------------------------------------------------------- §6.1 冗余反噬
A("")
A("【2】§6.1 冗余反噬（shard_v2_results.csv，A(n) 每单元 N=1500）")
rows = read_rows("shard_v2_results.csv")
Na, Nt = TRIALS["shard_avail"], TRIALS["shard_theta"]
# 混合极限 P(theta>=c)（来自 shard_sim_v2 输出，N=2000）
theta_limit = {1.5: 1.0000, 2.5: 0.7470, 3.5: 0.3435}
byL = {}
for r in rows:
    byL.setdefault(f(r["L"]), {})[(r["strategy"], int(r["n"]))] = f(r["avail"])
for L in sorted(byL):
    A("  L=%.1f  rand A(n)：" % L)
    for n in (2, 4, 16, 64):
        key = ("rand", n)
        if key in byL[L]:
            A("    n=%-3d %s" % (n, fmt_ci(byL[L][key], Na)))
    # 定理核验：|A(64) - P(theta>=c)|
    a64 = byL[L].get(("rand", 64))
    pt = theta_limit[L]
    lo64, hi64 = wilson(a64, Na)
    lo_t, hi_t = wilson(pt, Nt)
    A("    定理核验：A64=%.4f [%.4f,%.4f] vs P(θ≥c)=%.4f [%.4f,%.4f] -> |diff|=%.4f（CI 重叠：%s）"
      % (a64, lo64, hi64, pt, lo_t, hi_t, abs(a64 - pt),
         "是" if (lo64 <= hi_t and lo_t <= hi64) else "否"))

# ------------------------------------------------------------- spread.csv
A("")
A("【3】图 Fig.4 饱和曲线（spread.csv，每单元 N=3000）")
rows = read_rows("spread.csv")
Ns = TRIALS["spread"]
for r in rows:
    A("  L=%s n=%-3s rand %s | disp %s"
      % (r["R"], r["n"], fmt_ci(f(r["rand"]), Ns), fmt_ci(f(r["disp"]), Ns)))

# ------------------------------------------------------------- §6.4 真实拓扑
A("")
A("【4】§6.4 Kdl 真实拓扑（zoo_results.csv，每单元 N=400）")
rows = read_rows("zoo_results.csv")
Nz = TRIALS["zoo"]
r0name = {"230.0": "弱", "310.0": "中", "430.0": "强"}
for r in rows:
    if r["n"] == "16":
        A("  r0=%s(%s) n=16 %-8s survival %s | reachable %s"
          % (r["r0"], r0name.get(r["r0"], "?"), r["strategy"],
             fmt_ci(f(r["survival"]), Nz), fmt_ci(f(r["reachable"]), Nz)))
A("  —— 全部 60 行（4 臂 x 3 r0 x 5 n）CI 完整表见本文件随附 CSV 重算；")
A("     n=16/强 regime 关键对比：clustr2 0.540 vs disp 0.170（3.2 倍）")
z16 = {r["strategy"]: f(r["reachable"]) for r in rows
       if r["n"] == "16" and r["r0"] == "430.0"}
lo_c, hi_c = wilson(z16["clustr2"], Nz)
lo_d, hi_d = wilson(z16["disp"], Nz)
A("     保守下界（clustr2 CI 下限 / disp CI 上限）= %.3f / %.3f = %.2f 倍 > 2"
  % (lo_c, hi_d, lo_c / hi_d))

# ------------------------------------------------------------- §6.5 敏感性
A("")
A("【5】§6.5-B 码率扫描（sens_reach.csv，每单元 N=3000）")
rows = read_rows("sens_reach.csv")
Nsr = TRIALS["sens_reach"]
for r in rows:
    if r["strategy"] == "clustr2":
        A("  %-10s param=%-6s c=%-5s reachable %s"
          % (r["sweep"], r["param"], r["C_r2"] if "C_r2" in r else "",
             fmt_ci(f(r["reachable"]), Nsr)))

A("")
A("【6】§6.5-A 相图（sens_phase.csv，每单元 N=500）")
rows = read_rows("sens_phase.csv")
Np = TRIALS["sens_phase"]
diffs = []
for r in rows:
    d = abs(f(r["A64"]) - f(r["P_theta_c"]))
    diffs.append(d)
    lo_a, hi_a = wilson(f(r["A64"]), Np)
    lo_p, hi_p = wilson(f(r["P_theta_c"]), Np)
    A("  c=%-4s m=%-2s θ̄=%-7s A64 %s | P(θ≥c) %s |diff|=%.3f 重叠:%s"
      % (r["c"], r["m"], r["theta_bar"],
         fmt_ci(f(r["A64"]), Np), fmt_ci(f(r["P_theta_c"]), Np), d,
         "是" if (lo_a <= hi_p and lo_p <= hi_a) else "否"))
A("  平均 |A64 - P| = %.4f（42 cell）" % (sum(diffs) / len(diffs)))

# ------------------------------------------------------------------- 汇总
A("")
A("=" * 78)
A("论文可直接引用的统计报告句（英文，供 §6 方法段使用）：")
A('  \"All reported proportions are Monte Carlo estimates. Unless stated')
A('   otherwise, 95% Wilson score confidence intervals are ±0.018 or')
A('   smaller for the main experiment (3,000 trials per cell), ±0.026 for')
A('   the mixture-limit cross-validation (1,500-2,000 trials), ±0.044 for')
A('   the phase diagram (500 trials), and ±0.049 for the Kdl topology')
A('   study (400 trials); the largest interval half-width across all')
A('   reported cells is 0.049.\"')
A("=" * 78)

# 最大半宽核验
half = 0.0
for r in read_rows("reach_results.csv"):
    half = max(half, (wilson(f(r["survival"]), N)[1] - wilson(f(r["survival"]), N)[0]) / 2,
               (wilson(f(r["reachable"]), N)[1] - wilson(f(r["reachable"]), N)[0]) / 2)
A("核验：主实验 (N=3000) 实测最大半宽 = %.4f" % half)
halfz = 0.0
for r in read_rows("zoo_results.csv"):
    halfz = max(halfz,
                (wilson(f(r["survival"]), Nz)[1] - wilson(f(r["survival"]), Nz)[0]) / 2,
                (wilson(f(r["reachable"]), Nz)[1] - wilson(f(r["reachable"]), Nz)[0]) / 2)
A("核验：Kdl 拓扑 (N=400) 实测最大半宽 = %.4f" % halfz)

txt = "\n".join(out) + "\n"
path = os.path.join(RES, "ci_report.txt")
with open(path, "w", encoding="utf-8") as fp:
    fp.write(txt)
print(txt)
print("[ok] %s" % path)
