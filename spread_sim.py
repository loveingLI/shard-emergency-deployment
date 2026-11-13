"""
SPREAD 原型 —— 抗空间相关失效的碎片化部署
验证：① 副本收益在相关失效下饱和 ② 几何构型（最大最小散布）优于随机 ③ 最优碎片度

纯标准库。运行：python spread_sim.py
"""

import math
import random
import csv
import os
import statistics

AREA = 10.0
EPISODES = 3000


def place_random(n, rng):
    pts = []
    while len(pts) < n:
        p = (rng.uniform(0, AREA), rng.uniform(0, AREA))
        if all(math.dist(p, q) > 0.15 for q in pts):
            pts.append(p)
    return pts


def place_dispersion(n, rng=None):
    """最大最小散布：贪心取距已有点集最远的候选点（网格离散化）"""
    grid = [(x * 0.5, y * 0.5) for x in range(21) for y in range(21)]
    start = max(grid, key=lambda p: (p[0] - 5) ** 2 + (p[1] - 5) ** 2)
    pts = [start]
    while len(pts) < n:
        best, bestd = None, -1
        for c in grid:
            d = min(math.dist(c, q) for q in pts)
            if d > bestd:
                best, bestd = c, d
        pts.append(best)
    return pts


def survival(pts, episodes, r0, seed=0, m=6):
    """簇过程毁伤：每轮生成 m 个半径 r0 的毁伤斑块，落在任一斑块内即失效。
    r0 表征灾害的空间相关长度。返回存活碎片数 >= k 的概率，k = ceil(n/2)。"""
    rng = random.Random(seed)
    n = len(pts)
    k = max(1, (n + 1) // 2)
    ok = 0
    for _ in range(episodes):
        centers = [(rng.uniform(-r0, AREA + r0), rng.uniform(-r0, AREA + r0))
                   for _ in range(m)]
        alive = sum(1 for p in pts if all(math.dist(c, p) > r0 for c in centers))
        if alive >= k:
            ok += 1
    return ok / episodes


def marginal_saturation(ns, vals, tol=0.003):
    """返回可用率取最大值时的碎片数（最优碎片度）"""
    best = max(range(len(vals)), key=lambda i: vals[i])
    return ns[best]


def main():
    os.makedirs("results", exist_ok=True)
    ns = [2, 4, 6, 8, 10, 12, 16, 20, 24]
    Rs = [1.5, 2.5, 3.5]
    rows = []
    for R in Rs:
        for n in ns:
            r_rand = statistics.mean(
                survival(place_random(n, random.Random(1000 + n * 13 + i * 97)),
                         EPISODES, R, seed=i, m=6) for i in range(5))
            r_disp = survival(place_dispersion(n), EPISODES, R, seed=7, m=6)
            rows.append({"R": R, "n": n, "rand": round(r_rand, 4),
                         "disp": round(r_disp, 4)})

    with open("results/spread.csv", "w", newline="", encoding="utf-8") as f:
        wr = csv.DictWriter(f, fieldnames=["R", "n", "rand", "disp"])
        wr.writeheader()
        wr.writerows(rows)

    lines = ["相关长度R   n   随机放置   最大最小散布   散布增益"]
    for R in Rs:
        sub = [x for x in rows if x["R"] == R]
        vals_d = [x["disp"] for x in sub]
        n_star = marginal_saturation(ns, vals_d)
        lines.append("  %.1f      -     -          -            (饱和点 n*=%d)" % (R, n_star))
        for x in sub:
            lines.append("           %-4d %.4f     %.4f        %+.2f" %
                         (x["n"], x["rand"], x["disp"], x["disp"] - x["rand"]))

    summ = []
    summ.append("%-7s %-8s %-10s %-11s %-10s %s" %
                ("相关长度", "饱和点n*", "n*可用率", "n=24可用率", "随机增益", "结论"))
    for R in Rs:
        sub = [x for x in rows if x["R"] == R]
        vals_d = [x["disp"] for x in sub]
        n_star = marginal_saturation(ns, vals_d)
        idx = ns.index(n_star) if n_star in ns else len(ns) - 1
        gain = vals_d[-1] - sub[-1]["rand"]
        summ.append("%-7.1f %-8d %-10.4f %-11.4f %-10.3f %s" %
                    (R, n_star, vals_d[idx], vals_d[-1], gain,
                     "碎片化收益明显" if n_star >= 6 else "收益提前饱和"))

    report = "\n".join(summ) + "\n\n" + "\n".join(lines)
    with open("results/spread_summary.txt", "w", encoding="utf-8") as f:
        f.write(report)
    print("\n".join(summ))

    svg = make_svg(ns, Rs, rows)
    with open("results/fig_saturation.svg", "w", encoding="utf-8") as f:
        f.write(svg)
    print("\n[ok] results/spread.csv, results/spread_summary.txt, results/fig_saturation.svg")


def make_svg(ns, Rs, rows):
    W, H, L, R_, T, B = 700, 400, 68, 175, 46, 58
    pw, ph = W - L - R_, H - T - B
    colors = ["#1D9E75", "#378ADD", "#BA7517", "#8A5CF6", "#C0392B"]

    def X(i):
        return L + pw * i / (len(ns) - 1)

    def Y(v):
        return T + ph * (1 - v)

    s = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" width="100%%">' % (W, H),
         '<rect width="%d" height="%d" fill="#ffffff"/>' % (W, H),
         '<text x="%d" y="24" font-size="14" font-family="sans-serif" fill="#2c2c2a">Service survival probability vs. fragment count</text>' % L]
    for k in range(5):
        yv = k / 4
        y = Y(yv)
        s.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#e8e8e3" stroke-width="0.6"/>' % (L, y, L + pw, y))
        s.append('<text x="%d" y="%.1f" font-size="11" text-anchor="end" font-family="sans-serif" fill="#5f5e5a">%.2f</text>' % (L - 8, y + 4, yv))
    for i, n in enumerate(ns):
        s.append('<text x="%.1f" y="%d" font-size="11" text-anchor="middle" font-family="sans-serif" fill="#5f5e5a">%d</text>' % (X(i), T + ph + 18, n))
    s.append('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#9a9a90" stroke-width="1"/>' % (L, T + ph, L + pw, T + ph))
    s.append('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#9a9a90" stroke-width="1"/>' % (L, T, L, T + ph))
    for ci, R in enumerate(Rs):
        sub = [x for x in rows if x["R"] == R]
        pts = " ".join("%.1f,%.1f" % (X(i), Y(x["disp"])) for i, x in enumerate(sub))
        s.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="2"/>' % (pts, colors[ci]))
        pts2 = " ".join("%.1f,%.1f" % (X(i), Y(x["rand"])) for i, x in enumerate(sub))
        s.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="1.2" stroke-dasharray="4 3"/>' % (pts2, colors[ci]))
    s.append('<text x="%d" y="%d" font-size="12" font-family="sans-serif" fill="#5f5e5a" text-anchor="middle">Fragment count n</text>' % (L + pw // 2, H - 20))
    s.append('<text x="14" y="%d" font-size="12" font-family="sans-serif" fill="#5f5e5a" transform="rotate(-90 14 %d)" text-anchor="middle">Survival probability</text>' % (T + ph // 2, T + ph // 2))
    ly = T + 6
    for ci, R in enumerate(Rs):
        s.append('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="%s" stroke-width="2"/>' % (L + pw + 14, ly, L + pw + 36, ly, colors[ci]))
        s.append('<text x="%d" y="%d" font-size="11" font-family="sans-serif" fill="#2c2c2a">R = %.1f km</text>' % (L + pw + 42, ly + 4, R))
        ly += 18
    s.append('<text x="%d" y="%d" font-size="10" font-family="sans-serif" fill="#7a7a70">实线 最大最小散布 / 虚线 随机</text>' % (L + pw + 14, ly + 8))
    s.append('</svg>')
    return "".join(s)


if __name__ == "__main__":
    main()
