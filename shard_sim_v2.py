"""
SHARD v2 — 论文主实验：碎片化放置与空间相关失效
====================================================
三臂对比：随机放置(rand) / 最大最小散布(disp) / 分簇+簇内冗余(clust, SHARD 设计规则)
regime：灾害相关长度 L ∈ {1.5, 2.5, 3.5} km（簇过程毁伤：m=6 个半径 L 的斑块）
指标：A(n) = P(存活碎片数 ≥ k)，k = ⌈n/2⌉（MDS 码率 c = 1/2）
理论交叉验证：混合极限定理  lim_{n→∞} A(n) = P(θ ≥ c)
              θ = 单次灾情的未覆盖面积比例（Monte Carlo 估计）

运行（在 sim/ 目录下）：python shard_sim_v2.py
依赖：仅 Python 3 标准库（无第三方依赖，可直接开源）

输出：results/shard_v2_results.csv    —— 全部原始数据（每行一个 (L, n, strategy) 单元）
      results/shard_v2_summary.txt    —— 汇总表（论文 Table 3 数据源）
      results/fig_regime_L15/25/35.svg —— 三个 regime 的 A(n) 曲线（论文 Fig.4）
      results/fig_theory.svg           —— 散布曲线 + 定理极限线（论文 Fig.5）
      results/fig_flip.svg             —— n=48 三策略对比（论文 Fig.6，几何反转）
"""

import math
import random
import csv
import os

# ----------------------------- 实验参数 -----------------------------
AREA = 10.0          # 部署区域边长 (km)
M_DISKS = 6          # 每轮毁伤斑块数（灾情强度）
K_RATE = 0.5         # MDS 码率 c = k/n
EPISODES = 1500      # 每单元蒙特卡洛轮数
TH_EPISODES = 2000   # θ 尾概率估计轮数
NS = [1, 2, 4, 6, 8, 10, 12, 16, 20, 24, 32, 48, 64]
R0S = [1.5, 2.5, 3.5]
RAND_DRAWS = 3       # 随机放置平均的独立抽样次数
CLUSTER_JITTER = 0.4 # 簇内碎片散布半径 (km)

STRATEGIES = ("rand", "disp", "clust")
GRID = [(x * 0.5, y * 0.5) for x in range(21) for y in range(21)]
TH_GRID = [(x * 0.5 + 0.25, y * 0.5 + 0.25) for x in range(20) for y in range(20)]


# ----------------------------- 灾害模型 -----------------------------
def disaster(rng, r0):
    """簇过程毁伤：m 个圆盘斑块，中心均匀分布于扩展区域（相关长度 = r0）"""
    lo, hi = -r0, AREA + r0
    return [(rng.uniform(lo, hi), rng.uniform(lo, hi)) for _ in range(M_DISKS)]


def survives(p, centers, r0):
    """点 p 是否被本轮毁伤覆盖"""
    r2 = r0 * r0
    for c in centers:
        dx = p[0] - c[0]
        dy = p[1] - c[1]
        if dx * dx + dy * dy <= r2:
            return False
    return True


# ----------------------------- 放置策略 -----------------------------
def place_rand(n, rng):
    return [(rng.uniform(0, AREA), rng.uniform(0, AREA)) for _ in range(n)]


def place_dispersion(n):
    """最大最小散布：贪心选距已选点最远的网格候选点"""
    pts = [max(GRID, key=lambda p: (p[0] - 5) ** 2 + (p[1] - 5) ** 2)]
    while len(pts) < n:
        best, best_d = None, -1.0
        for c in GRID:
            d = min((c[0] - q[0]) ** 2 + (c[1] - q[1]) ** 2 for q in pts)
            if d > best_d:
                best, best_d = c, d
        pts.append(best)
    return pts


def feasible_cluster_count(r0):
    """设计规则：两两间距 ≥ 2L 的独立失效单元数（区域几何上限）"""
    sep = 2.0 * r0
    centers = [max(GRID, key=lambda p: (p[0] - 5) ** 2 + (p[1] - 5) ** 2)]
    while len(centers) < 8:
        cands = [p for p in GRID
                 if all((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 >= sep * sep
                        for q in centers)]
        if not cands:
            break
        best = max(cands, key=lambda p: min((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2
                                            for q in centers))
        centers.append(best)
    return len(centers), centers


def place_clustered(n, r0, rng):
    """SHARD 布设：簇数 = 独立失效单元数（≤8），簇间距 > 2L，n 片均分入簇"""
    c_cap, _ = feasible_cluster_count(r0)
    C = max(1, min(c_cap, n))
    centers = feasible_cluster_count(r0)[1][:C]
    per = [n // C] * C
    for i in range(n - sum(per)):
        per[i] += 1
    pts = []
    for c, m in zip(centers, per):
        for _ in range(m):
            pts.append((c[0] + rng.uniform(-CLUSTER_JITTER, CLUSTER_JITTER),
                        c[1] + rng.uniform(-CLUSTER_JITTER, CLUSTER_JITTER)))
    return pts, C


# ----------------------------- 测度 -----------------------------
def availability(pts, r0, seed):
    """A(n) = P(存活碎片 ≥ ⌈n/2⌉)"""
    rng = random.Random(seed)
    n = len(pts)
    k = max(1, math.ceil(n * K_RATE))
    ok = 0
    for _ in range(EPISODES):
        cs = disaster(rng, r0)
        alive = 0
        for p in pts:
            if survives(p, cs, r0):
                alive += 1
        if alive >= k:
            ok += 1
    return ok / EPISODES


def theta_tail(r0, seed):
    """定理极限：P(θ ≥ c)，θ 为未覆盖面积比例（细网格估计）"""
    rng = random.Random(seed)
    ok = 0
    total = len(TH_GRID)
    for _ in range(TH_EPISODES):
        cs = disaster(rng, r0)
        alive = sum(1 for p in TH_GRID if survives(p, cs, r0))
        if alive / total >= K_RATE - 1e-9:
            ok += 1
    return ok / TH_EPISODES


# ----------------------------- SVG 绘图 -----------------------------
def polyline_svg(path, title, ylabel, xlabels, series, ymax=1.0):
    """series: [(label, color, values, dashed)]，values 与 xlabels 对齐"""
    W, H, L, R, T, B = 700, 400, 70, 200, 46, 60
    pw, ph = W - L - R, H - T - B

    def X(i):
        return L + pw * i / max(len(xlabels) - 1, 1)

    def Y(v):
        return T + ph * (1 - min(v, ymax) / ymax)

    s = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" width="100%%">' % (W, H),
         '<rect width="%d" height="%d" fill="#ffffff"/>' % (W, H),
         '<text x="%d" y="24" font-size="14" font-family="sans-serif" fill="#2c2c2a">%s</text>' % (L, title)]
    for g in range(5):
        yv = ymax * g / 4
        y = Y(yv)
        s.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#e8e8e3" stroke-width="0.6"/>' % (L, y, L + pw, y))
        s.append('<text x="%d" y="%.1f" font-size="11" text-anchor="end" font-family="sans-serif" fill="#5f5e5a">%.2f</text>' % (L - 8, y + 4, yv))
    for i, lab in enumerate(xlabels):
        s.append('<text x="%.1f" y="%d" font-size="11" text-anchor="middle" font-family="sans-serif" fill="#5f5e5a">%s</text>' % (X(i), T + ph + 18, lab))
    s.append('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#9a9a90" stroke-width="1"/>' % (L, T + ph, L + pw, T + ph))
    s.append('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#9a9a90" stroke-width="1"/>' % (L, T, L, T + ph))
    for label, color, vals, dashed in series:
        pts = " ".join("%.1f,%.1f" % (X(i), v if v is not None and not math.isnan(v) else T + ph)
                       for i, v in enumerate(vals))
        dash = ' stroke-dasharray="5 4"' if dashed else ""
        s.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="%d"%s/>' % (pts, color, 2 if not dashed else 1.5, dash))
        for i, v in enumerate(vals):
            if v is not None and not math.isnan(v):
                s.append('<circle cx="%.1f" cy="%.1f" r="3" fill="%s"/>' % (X(i), Y(v), color))
    s.append('<text x="%d" y="%d" font-size="12" font-family="sans-serif" fill="#5f5e5a" text-anchor="middle">Fragment count n</text>' % (L + pw // 2, H - 20))
    s.append('<text x="16" y="%d" font-size="12" font-family="sans-serif" fill="#5f5e5a" transform="rotate(-90 16 %d)" text-anchor="middle">%s</text>' % (T + ph // 2, T + ph // 2, ylabel))
    ly = T + 6
    for label, color, _, dashed in series:
        dash = ' stroke-dasharray="5 4"' if dashed else ""
        s.append('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="%s" stroke-width="2"%s/>' % (L + pw + 16, ly, L + pw + 40, ly, color, dash))
        s.append('<text x="%d" y="%d" font-size="11" font-family="sans-serif" fill="#2c2c2a">%s</text>' % (L + pw + 46, ly + 4, label))
        ly += 20
    s.append('</svg>')
    with open(path, "w", encoding="utf-8") as f:
        f.write("".join(s))


def bars_svg(path, title, ylabel, groups, series, ymax=1.0):
    """groups: [group_label]; series: [(label, color, values per group)]"""
    W, H, L, T = 700, 400, 70, 46
    pw, ph = W - L - 60, H - T - 70
    gw = pw / len(groups)
    bw = gw / (len(series) + 1)

    def Y(v):
        return T + ph * (1 - min(v, ymax) / ymax)

    s = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" width="100%%">' % (W, H),
         '<rect width="%d" height="%d" fill="#ffffff"/>' % (W, H),
         '<text x="%d" y="24" font-size="14" font-family="sans-serif" fill="#2c2c2a">%s</text>' % (L, title)]
    for g in range(5):
        yv = ymax * g / 4
        s.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#e8e8e3" stroke-width="0.6"/>' % (L, Y(yv), L + pw, Y(yv)))
        s.append('<text x="%d" y="%.1f" font-size="11" text-anchor="end" font-family="sans-serif" fill="#5f5e5a">%.2f</text>' % (L - 8, Y(yv) + 4, yv))
    for gi, glab in enumerate(groups):
        x0 = L + gi * gw
        s.append('<text x="%.1f" y="%d" font-size="12" text-anchor="middle" font-family="sans-serif" fill="#2c2c2a">%s</text>' % (x0 + gw / 2, T + ph + 20, glab))
        for si, (label, color, vals) in enumerate(series):
            bx = x0 + (si + 0.5) * bw
            v = vals[gi]
            hgt = ph * v / ymax
            s.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s" stroke="none"/>' % (bx, T + ph - hgt, bw * 0.8, hgt, color))
            s.append('<text x="%.1f" y="%.1f" font-size="10" text-anchor="middle" font-family="sans-serif" fill="#2c2c2a">%.2f</text>' % (bx + bw * 0.4, T + ph - hgt - 5, v))
    s.append('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#9a9a90" stroke-width="1"/>' % (L, T + ph, L + pw, T + ph))
    ly = T + 6
    for label, color, _ in series:
        s.append('<rect x="%d" y="%d" width="14" height="10" fill="%s"/>' % (L + pw + 16, ly - 8, color))
        s.append('<text x="%d" y="%d" font-size="11" font-family="sans-serif" fill="#2c2c2a">%s</text>' % (L + pw + 36, ly + 4, label))
        ly += 22
    s.append('<text x="16" y="%d" font-size="12" font-family="sans-serif" fill="#5f5e5a" transform="rotate(-90 16 %d)" text-anchor="middle">%s</text>' % (T + ph // 2, T + ph // 2, ylabel))
    s.append('</svg>')
    with open(path, "w", encoding="utf-8") as f:
        f.write("".join(s))


# ----------------------------- 主流程 -----------------------------
def main():
    os.makedirs("results", exist_ok=True)
    rows = []
    theta_limit = {}

    print(">> 估计混合极限 P(theta >= c) ...")
    for r0 in R0S:
        theta_limit[r0] = theta_tail(r0, seed=int(r0 * 1000))
        print("   L=%.1f  P(theta>=c) = %.4f" % (r0, theta_limit[r0]))

    print(">> 主实验（3 regime x %d 个 n x 3 策略）..." % len(NS))
    for r0 in R0S:
        c_cap, _ = feasible_cluster_count(r0)
        for n in NS:
            for si, strat in enumerate(STRATEGIES):
                base_seed = int(r0 * 1000) * 100 + n * 3 + si
                if strat == "rand":
                    vals = []
                    for d in range(RAND_DRAWS):
                        pts = place_rand(n, random.Random(base_seed * 7 + d))
                        vals.append(availability(pts, r0, base_seed * 13 + d))
                    a = sum(vals) / len(vals)
                elif strat == "disp":
                    pts = place_dispersion(n)
                    a = availability(pts, r0, base_seed * 13)
                else:
                    pts, C = place_clustered(n, r0, random.Random(base_seed * 5))
                    a = availability(pts, r0, base_seed * 13)
                    rows.append({"L": r0, "n": n, "strategy": strat,
                                 "clusters": C, "avail": round(a, 4)})
                    continue
                rows.append({"L": r0, "n": n, "strategy": strat,
                             "clusters": c_cap if strat == "disp" else 1,
                             "avail": round(a, 4)})
        print("   L=%.1f done (可行独立单元数 C_max=%d)" % (r0, c_cap))

    with open("results/shard_v2_results.csv", "w", newline="", encoding="utf-8") as f:
        wr = csv.DictWriter(f, fieldnames=["L", "n", "strategy", "clusters", "avail"])
        wr.writeheader()
        wr.writerows(rows)

    # 汇总表
    lines = []
    lines.append("可行独立失效单元数（区域上限，两两间距≥2L）：")
    for r0 in R0S:
        c_cap, _ = feasible_cluster_count(r0)
        lines.append("  L=%.1f km -> C_max=%d" % (r0, c_cap))
    lines.append("")
    lines.append("混合极限定理：lim A(n) = P(theta >= c), c = %.1f" % K_RATE)
    for r0 in R0S:
        lines.append("  L=%.1f -> %.4f" % (r0, theta_limit[r0]))
    lines.append("")
    header = "%-6s %-4s | %-8s %-8s %-8s" % ("L", "n", "rand", "disp", "clust")
    lines.append(header)
    lines.append("-" * len(header))
    for r0 in R0S:
        for n in NS:
            vals = []
            for strat in STRATEGIES:
                m = [x["avail"] for x in rows if x["L"] == r0 and x["n"] == n and x["strategy"] == strat]
                vals.append("%.4f" % m[0] if m else "-")
            lines.append("%-6.1f %-4d | %-8s %-8s %-8s" % (r0, n, vals[0], vals[1], vals[2]))
        lines.append("")
    report = "\n".join(lines)
    with open("results/shard_v2_summary.txt", "w", encoding="utf-8") as f:
        f.write(report + "\n")
    print(report)

    # 图
    colors = {"rand": "#B4B2A9", "disp": "#378ADD", "clust": "#1D9E75"}
    names = {"rand": "Random", "disp": "Dispersion", "clust": "Clustered (SHARD)"}
    for r0 in R0S:
        series = []
        for strat in STRATEGIES:
            vals = [next((x["avail"] for x in rows
                          if x["L"] == r0 and x["n"] == n and x["strategy"] == strat), None)
                    for n in NS]
            series.append((names[strat], colors[strat], vals, False))
        series.append(("Theorem limit P(theta>=c)", "#C0392B",
                       [theta_limit[r0]] * len(NS), True))
        polyline_svg("results/fig_regime_L%d.svg" % int(r0 * 10),
                     "Availability vs. fragment count  (L = %.1f km)" % r0,
                     "Availability A(n)", [str(x) for x in NS], series)

    series = []
    for r0 in R0S:
        vals = [next((x["avail"] for x in rows
                      if x["L"] == r0 and x["n"] == n and x["strategy"] == "disp"), None)
                for n in NS]
        series.append(("Dispersion, L=%.1f" % r0, colors_disp(r0), vals, False))
    for r0 in R0S:
        series.append(("limit L=%.1f" % r0, colors_disp(r0),
                       [theta_limit[r0]] * len(NS), True))
    polyline_svg("results/fig_theory.svg",
                 "Dispersion A(n) converges to the mixture limit P(theta >= c)",
                 "Availability A(n)", [str(x) for x in NS], series)

    n_bar = 48
    groups = ["L=%.1f km" % r0 for r0 in R0S]
    series = []
    for strat in STRATEGIES:
        vals = [next((x["avail"] for x in rows
                      if x["L"] == r0 and x["n"] == n_bar and x["strategy"] == strat), 0)
                for r0 in R0S]
        series.append((names[strat], colors[strat], vals))
    bars_svg("results/fig_flip.svg",
             "Geometry flip at n = %d: clustered dominates under strong correlation" % n_bar,
             "Availability", groups, series)

    print("\n[ok] results/shard_v2_results.csv")
    print("[ok] results/shard_v2_summary.txt")
    print("[ok] results/fig_regime_*.svg / fig_theory.svg / fig_flip.svg")


def colors_disp(r0):
    return {1.5: "#1D9E75", 2.5: "#378ADD", 3.5: "#C0392B"}[r0]


if __name__ == "__main__":
    main()
