"""
SHARD 连通层实验：生存 ≠ 可达（预审 P0-4 修复）+ R2 修正分簇臂（v2）
=====================================================================
在随机几何图上放置碎片，簇过程毁伤摧毁节点后：
  - 几何存活率 survival = P(存活碎片数 ≥ k)            （不考虑连通）
  - 可达服务率 reachable = P(同一连通分量内存活碎片 ≥ k)  （k-anycast 真正能凑齐）
两者的差 = 网络分区的服务惩罚，即论文 Fig.7 的数据源。

v2（2026-11-14）：新增第四臂 clustr2 —— R2 修正分簇。
  朴素 clust：簇数固定 ≤4，每簇 n/C 片，在 k=⌈n/2⌉ 下每簇 < k，违反 R2
  基数规则（每簇 ≥ k），分区后任何分量都凑不齐 k 片 → 可达崩塌。
  clustr2：簇数 C = max(1, ⌊n/k⌋)，保证每簇 ≥ k 片（R2 满足），
  簇心仍保持 ≥ 2L 分离。验证：R2 修正后分簇的可达率是否恢复。
  这把论文的 R2 从"架构规则"升级为"有实验证据的必要约束"。

运行（在 sim/ 目录下）：python reach_sim.py    （纯标准库）
输出：results/reach_results.csv
      results/reach_summary.txt
      results/fig_reach.svg
"""

import math
import random
import csv
import os
import statistics

AREA = 10.0
N_NODES = 80
COMM_RANGE = 1.8
M_DISKS = 6
EPISODES = 600
NS = [4, 8, 12, 16, 24]
R0S = [1.5, 2.5, 3.5]
STRATEGIES = ("rand", "disp", "clust", "clustr2")


def build_graph(seed):
    rng = random.Random(seed)
    pos = [(rng.uniform(0, AREA), rng.uniform(0, AREA)) for _ in range(N_NODES)]
    adj = [[] for _ in range(N_NODES)]
    for i in range(N_NODES):
        for j in range(i + 1, N_NODES):
            if math.dist(pos[i], pos[j]) <= COMM_RANGE:
                adj[i].append(j)
                adj[j].append(i)
    return pos, adj


def disaster(rng, r0):
    lo, hi = -r0, AREA + r0
    return [(rng.uniform(lo, hi), rng.uniform(lo, hi)) for _ in range(M_DISKS)]


def alive_nodes(pos, centers, r0):
    r2 = r0 * r0
    out = []
    for i, p in enumerate(pos):
        dead = False
        for c in centers:
            dx = p[0] - c[0]
            dy = p[1] - c[1]
            if dx * dx + dy * dy <= r2:
                dead = True
                break
        if not dead:
            out.append(i)
    return out


def components(adj, alive):
    seen = set()
    comps = []
    for s in alive:
        if s in seen:
            continue
        stack = [s]
        seen.add(s)
        comp = []
        while stack:
            u = stack.pop()
            comp.append(u)
            for v in adj[u]:
                if v in alive_set and v not in seen:
                    seen.add(v)
                    stack.append(v)
        comps.append(comp)
    return comps


def pick_shards(strategy, n, pos, adj, r0, rng):
    """返回承载碎片的节点下标列表"""
    if strategy == "rand":
        return rng.sample(range(N_NODES), n)
    if strategy == "clustr2":
        return pick_shards_r2(n, pos, r0)
    if strategy == "disp":
        # 节点上的最大最小散布
        first = max(range(N_NODES), key=lambda i: (pos[i][0] - 5) ** 2 + (pos[i][1] - 5) ** 2)
        chosen = [first]
        while len(chosen) < n:
            best, bd = None, -1
            for i in range(N_NODES):
                if i in chosen:
                    continue
                d = min((pos[i][0] - pos[j][0]) ** 2 + (pos[i][1] - pos[j][1]) ** 2
                        for j in chosen)
                if d > bd:
                    best, bd = i, d
            chosen.append(best)
        return chosen
    # clustered：先选彼此距离 ≥ 2L 的簇心节点（≤4），每簇用邻近节点填充
    sep = 2 * r0
    first = max(range(N_NODES), key=lambda i: (pos[i][0] - 5) ** 2 + (pos[i][1] - 5) ** 2)
    centers = [first]
    while len(centers) < min(4, n):
        cands = [i for i in range(N_NODES)
                 if i not in centers and
                 all(math.dist(pos[i], pos[c]) >= sep for c in centers)]
        if not cands:
            break
        best = max(cands, key=lambda i: min(math.dist(pos[i], pos[c]) for c in centers))
        centers.append(best)
    per = [n // len(centers)] * len(centers)
    for i in range(n - sum(per)):
        per[i] += 1
    shards = []
    used = set()
    for c, m in zip(centers, per):
        near = sorted((i for i in range(N_NODES) if i not in used),
                      key=lambda i: math.dist(pos[i], pos[c]))
        take = near[:m]
        if len(take) < m:  # 回退：补最近可用节点
            for i in near:
                if i not in used and i not in take:
                    take.append(i)
                if len(take) >= m:
                    break
        shards.extend(take)
        used.update(take)
    return shards[:n]


def k_of(n):
    return max(1, math.ceil(n * 0.5))


def pick_shards_r2(n, pos, r0):
    """R2 修正分簇：簇数 C = max(1, n // k)，保证每簇 ≥ k 片。
    簇心仍要求彼此分离 ≥ 2L（放不下时退化为最大距离簇心，R2 优先于分离性）。"""
    k = k_of(n)
    cmax = max(1, n // k)
    sep = 2 * r0
    first = max(range(N_NODES), key=lambda i: (pos[i][0] - 5) ** 2 + (pos[i][1] - 5) ** 2)
    centers = [first]
    while len(centers) < cmax:
        far = [i for i in range(N_NODES)
               if i not in centers and
               all(math.dist(pos[i], pos[c]) >= sep for c in centers)]
        if far:
            best = max(far, key=lambda i: min(math.dist(pos[i], pos[c]) for c in centers))
        else:  # 分离约束放不下：R2 优先于分离性
            best = max((i for i in range(N_NODES) if i not in centers),
                       key=lambda i: min(math.dist(pos[i], pos[c]) for c in centers))
        centers.append(best)
    per = [n // len(centers)] * len(centers)
    for i in range(n - sum(per)):
        per[i] += 1
    # R2 自检：除单簇情形外，每簇片数必须 ≥ k
    assert len(centers) == 1 or min(per) >= k, "R2 violated: per=%s k=%d" % (per, k)
    shards = []
    used = set()
    for c, m in zip(centers, per):
        near = sorted((i for i in range(N_NODES) if i not in used),
                      key=lambda i: math.dist(pos[i], pos[c]))
        shards.extend(near[:m])
        used.update(near[:m])
    return shards[:n]


def run_cell(strategy, n, r0, seed):
    pos, adj = build_graph(seed)
    rng = random.Random(seed * 31 + n)
    shards = pick_shards(strategy, n, pos, adj, r0, rng)
    shard_set = set(shards)
    k = k_of(n)
    rng2 = random.Random(seed * 977 + int(r0 * 100))
    surv_ok = 0
    reach_ok = 0
    global alive_set
    for _ in range(EPISODES):
        cs = disaster(rng2, r0)
        alive = alive_nodes(pos, cs, r0)
        alive_set = set(alive)
        alive_shards = [s for s in shards if s in alive_set]
        if len(alive_shards) >= k:
            surv_ok += 1
        # 可达：同一分量内凑够 k 片
        comps = components(adj, alive)
        best = 0
        for comp in comps:
            cnt = sum(1 for s in comp if s in shard_set)
            best = max(best, cnt)
        if best >= k:
            reach_ok += 1
    return surv_ok / EPISODES, reach_ok / EPISODES


def main():
    global alive_set
    alive_set = set()
    os.makedirs("results", exist_ok=True)
    rows = []
    print(">> 连通层实验（3 regime x %d 个 n x %d 策略，每单元 %d 轮，5 图种子取均值）..."
          % (len(NS), len(STRATEGIES), EPISODES))
    for r0 in R0S:
        for n in NS:
            for strat in STRATEGIES:
                sv, rc = [], []
                for gseed in range(5):
                    s, r = run_cell(strat, n, r0, gseed + 1)
                    sv.append(s)
                    rc.append(r)
                rows.append({"L": r0, "n": n, "strategy": strat,
                             "survival": round(statistics.mean(sv), 4),
                             "reachable": round(statistics.mean(rc), 4)})
        print("   L=%.1f done" % r0)

    with open("results/reach_results.csv", "w", newline="", encoding="utf-8") as f:
        wr = csv.DictWriter(f, fieldnames=["L", "n", "strategy", "survival", "reachable"])
        wr.writeheader()
        wr.writerows(rows)

    lines = ["%-6s %-4s %-6s | %-9s %-9s %-8s" % ("L", "n", "strat", "survival", "reachable", "gap")]
    for x in rows:
        lines.append("%-6.1f %-4d %-6s | %-9.4f %-9.4f %-8.4f"
                     % (x["L"], x["n"], x["strategy"],
                        x["survival"], x["reachable"],
                        x["survival"] - x["reachable"]))
    report = "\n".join(lines)
    with open("results/reach_summary.txt", "w", encoding="utf-8") as f:
        f.write(report + "\n")
    print(report)

    # 图：n=16 时各 regime/策略的 存活 vs 可达（4 策略）
    n_bar = 16
    groups = ["%s\n%s" % ("L=%.1f" % r0, st) for r0 in R0S for st in STRATEGIES]
    surv_vals, reach_vals = [], []
    for r0 in R0S:
        for strat in STRATEGIES:
            m = [x for x in rows if x["L"] == r0 and x["n"] == n_bar and x["strategy"] == strat][0]
            surv_vals.append(m["survival"])
            reach_vals.append(m["reachable"])
    W, H, L, T = 1160, 420, 70, 50
    pw, ph = 850, H - T - 80   # 画区宽度固定，右侧留白放图例（W 加宽不挤压图例）
    gw = pw / len(groups)
    bw = gw * 0.36
    s = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" width="100%%">' % (W, H),
         '<rect width="%d" height="%d" fill="#ffffff"/>' % (W, H),
         '<text x="%d" y="24" font-size="14" font-family="sans-serif" fill="#2c2c2a">Survival vs. reachable service at n = %d (k-anycast, k = n/2)</text>' % (L, n_bar)]
    for g in range(5):
        yv = g / 4
        y = T + ph * (1 - yv)
        s.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#e8e8e3" stroke-width="0.6"/>' % (L, y, L + pw, y))
        s.append('<text x="%d" y="%.1f" font-size="11" text-anchor="end" font-family="sans-serif" fill="#5f5e5a">%.2f</text>' % (L - 8, y + 4, yv))
    for gi, glab in enumerate(groups):
        x0 = L + gi * gw
        parts = glab.split("\n")
        s.append('<text x="%.1f" y="%d" font-size="10" text-anchor="middle" font-family="sans-serif" fill="#2c2c2a">%s</text>' % (x0 + gw / 2, T + ph + 16, parts[0]))
        s.append('<text x="%.1f" y="%d" font-size="10" text-anchor="middle" font-family="sans-serif" fill="#5f5e5a">%s</text>' % (x0 + gw / 2, T + ph + 30, parts[1]))
        for si, (vals, color) in enumerate(((surv_vals, "#B4B2A9"), (reach_vals, "#1D9E75"))):
            bx = x0 + si * (bw + 2) + (gw - 2 * bw - 4) / 2
            hgt = ph * vals[gi]
            s.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s"/>' % (bx, T + ph - hgt, bw, hgt, color))
            s.append('<text x="%.1f" y="%.1f" font-size="9" text-anchor="middle" font-family="sans-serif" fill="#2c2c2a">%.2f</text>' % (bx + bw / 2, T + ph - hgt - 4, vals[gi]))
    s.append('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#9a9a90" stroke-width="1"/>' % (L, T + ph, L + pw, T + ph))
    ly = T + 6
    for label, color in (("Survival (no connectivity)", "#B4B2A9"),
                         ("Reachable service (k-anycast)", "#1D9E75")):
        s.append('<rect x="%d" y="%d" width="14" height="10" fill="%s"/>' % (L + pw + 20, ly - 8, color))
        s.append('<text x="%d" y="%d" font-size="11" font-family="sans-serif" fill="#2c2c2a">%s</text>' % (L + pw + 40, ly + 4, label))
        ly += 22
    s.append('</svg>')
    with open("results/fig_reach.svg", "w", encoding="utf-8") as f:
        f.write("".join(s))

    print("\n[ok] results/reach_results.csv / reach_summary.txt / fig_reach.svg")


if __name__ == "__main__":
    main()
