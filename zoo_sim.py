"""
SHARD Topology Zoo 实验：真实 ISP 拓扑上的连通层验证（v1, 2026-11-17）
=====================================================================
目的：把合成随机几何图上的"二次反转 + R2 必要性"结论（reach_sim.py）
搬到真实 ISP 拓扑上检验外部效度。拓扑：Topology Zoo 的 Kdl
（Kentucky Datalink，美国光纤骨干，726 个带经纬度的节点，约 2000x2200 km）。

灾害模型：簇过程（M_DISKS 个半径 r0 的毁伤圆盘，圆心均匀落在拓扑包围盒内）。
r0 按 E[theta] 校准为三个 regime（弱/中/强相关），与合成实验语义对齐：
theta = 毁伤后未覆盖面积比例（定理 1 的混合变量）。

四臂放置（与 reach_sim v2 相同）：
  rand    —— 均匀随机节点
  disp    —— 贪心最大-最小散布（确定性）
  clust   —— 朴素 4 簇（每簇 n/4 < k，违反 R2，作反例）
  clustr2 —— R2 修正分簇（C = floor(n/k) 簇，每簇 >= k 片，簇间距 > 2*r0）

度量：
  survival  = P(存活碎片数 >= k)              （几何存活）
  reachable = P(同一连通分量内存活碎片 >= k)   （k-anycast 真正可聚）
  theta_node/theta_area = 每轮实际节点/面积存活比例（报均值）

运行（在 sim/ 目录下）：python zoo_sim.py   （纯标准库，零依赖）
输出：results/zoo_results.csv
      results/zoo_summary.txt
      results/fig_zoo.svg  （3 regime 面板：可达率 vs n，4 臂）
对应论文：§6.5 Real-topology validation；数据文件随仓库发布。

数据来源：The Internet Topology Zoo (Knight et al., IEEE JSAC 2011),
http://www.topology-zoo.org — Kdl.graphml（经 mroughan/InternetTopologyZoo
镜像获取）。
"""

import xml.etree.ElementTree as ET
import math
import random
import csv
import os
import statistics

TOPO_FILE = os.path.join("topozoo", "Kdl.graphml")
M_DISKS = 6
R0S = (230.0, 310.0, 430.0)     # km，按 E[theta]≈{0.80,0.67,0.42} 校准（与合成实验 regime 对齐）
EPISODES = 400
NS = (4, 8, 12, 16, 24)
STRATEGIES = ("rand", "disp", "clust", "clustr2")
NSAMPLE_THETA = 400             # theta_area 的面积采样点数


# ---------------------------------------------------------------- 拓扑加载
def load_topo():
    """解析 GraphML -> (pos, adj, bbox)，只保留有经纬度的最大连通分量。"""
    ns = {"g": "http://graphml.graphdrawing.org/xmlns"}
    root = ET.parse(TOPO_FILE).getroot()
    graph = root.find("g:graph", ns)
    raw = []
    idx = {}
    for n in graph.findall("g:node", ns):
        d = {kd.get("key"): kd.text for kd in n.findall("g:data", ns)}
        if d.get("d31") and d.get("d34"):
            idx[n.get("id")] = len(raw)
            raw.append((float(d["d31"]), float(d["d34"])))
    # 等距圆柱投影 -> 平面 km
    lat0 = statistics.mean(p[0] for p in raw)
    pos = [((lo - min(p[1] for p in raw)) * 111.32 * math.cos(math.radians(lat0)),
            (la - min(p[0] for p in raw)) * 110.57)
           for la, lo in raw]
    adj = [[] for _ in pos]
    for e in graph.findall("g:edge", ns):
        s, t = idx.get(e.get("source")), idx.get(e.get("target"))
        if s is not None and t is not None:
            adj[s].append(t)
            adj[t].append(s)
    # 最大连通分量
    seen = set()
    best = []
    for s0 in range(len(pos)):
        if s0 in seen:
            continue
        comp, stack = [], [s0]
        seen.add(s0)
        while stack:
            u = stack.pop()
            comp.append(u)
            for v in adj[u]:
                if v not in seen:
                    seen.add(v)
                    stack.append(v)
        if len(comp) > len(best):
            best = comp
    keep = set(best)
    remap = {u: i for i, u in enumerate(best)}
    pos2 = [pos[u] for u in best]
    adj2 = [[] for _ in best]
    for i, u in enumerate(best):
        adj2[i] = [remap[v] for v in adj[u] if v in keep]
    xs = [p[0] for p in pos2]
    ys = [p[1] for p in pos2]
    bbox = (min(xs), max(xs), min(ys), max(ys))
    return pos2, adj2, bbox


# ---------------------------------------------------------------- 放置策略
def place_rand(n, pos, rng):
    return rng.sample(range(len(pos)), n)


def place_disp(n, pos, _rng=None):
    """确定性贪心最大-最小散布（按 n 缓存）。"""
    key = ("disp", n)
    if key not in _DISP_CACHE:
        N = len(pos)
        # 起点：包围盒四角各自最近的节点中取离中心最远者
        cx = sum(p[0] for p in pos) / N
        cy = sum(p[1] for p in pos) / N
        first = max(range(N), key=lambda i: math.dist(pos[i], (cx, cy)))
        sel = [first]
        mind = [math.dist(pos[i], pos[first]) for i in range(N)]
        for _ in range(n - 1):
            j = max(range(N), key=lambda i: mind[i])
            sel.append(j)
            dj = math.dist(pos[j], pos[i]) if False else None  # 占位防误用
            for i in range(N):
                dji = math.dist(pos[j], pos[i])
                if dji < mind[i]:
                    mind[i] = dji
        _DISP_CACHE[key] = sel
    return list(_DISP_CACHE[key])


def _clusters(pos, centers, per, used_all=None):
    """每个中心取最近的 per[c] 个节点（互不重复）。"""
    shards, used = [], set()
    for c, m in zip(centers, per):
        near = sorted((i for i in range(len(pos)) if i not in used),
                      key=lambda i: math.dist(pos[i], pos[c]))
        take = near[:m]
        shards.extend(take)
        used.update(take)
    return shards


def place_clust(n, k, pos, rng, C, min_spacing):
    """C 簇；min_spacing>0 时簇间距须 > min_spacing（R1），否则朴素随机圆心。"""
    x0, x1, y0, y1 = _BBOX
    for _ in range(600):
        centers = [(rng.uniform(x0, x1), rng.uniform(y0, y1)) for _ in range(C)]
        if min_spacing <= 0 or all(
                math.dist(centers[i], centers[j]) > min_spacing
                for i in range(C) for j in range(i)):
            break
    # 圆心落在节点平面上：映射为最近节点索引即可（_clusters 按几何距离取成员）
    cidx = []
    for c in centers:
        cidx.append(min(range(len(pos)), key=lambda i: math.dist(pos[i], c)))
    q, r = divmod(n, C)
    per = [q + (1 if a < r else 0) for a in range(C)]
    return _clusters(pos, cidx, per)


_DISP_CACHE = {}


# ---------------------------------------------------------------- 毁伤与度量
def disaster(pos, r0, rng):
    x0, x1, y0, y1 = _BBOX
    return [(rng.uniform(x0, x1), rng.uniform(y0, y1)) for _ in range(M_DISKS)]


def alive_mask(pos, disks, r0):
    return [all(math.dist(p, c) > r0 for c in disks) for p in pos]


def components(alive, adj):
    """存活子图的连通分量 -> 每个节点的分量编号。"""
    comp = [-1] * len(adj)
    nc = 0
    for s0 in range(len(adj)):
        if not alive[s0] or comp[s0] >= 0:
            continue
        comp[s0] = nc
        stack = [s0]
        while stack:
            u = stack.pop()
            for v in adj[u]:
                if alive[v] and comp[v] < 0:
                    comp[v] = nc
                    stack.append(v)
        nc += 1
    return comp, nc


def theta_area_sample(disks, r0, rng):
    x0, x1, y0, y1 = _BBOX
    ok = 0
    for _ in range(NSAMPLE_THETA):
        p = (rng.uniform(x0, x1), rng.uniform(y0, y1))
        if all(math.dist(p, c) > r0 for c in disks):
            ok += 1
    return ok / NSAMPLE_THETA


# ---------------------------------------------------------------- 主实验
def run_cell(strategy, n, r0, seed):
    rng = random.Random(hash((strategy, n, r0, seed)) & 0xFFFFFFFF)
    k = max(1, math.ceil(n * 0.5))
    surv, reach, th_n, th_a = 0, 0, [], []
    for _ in range(EPISODES):
        if strategy == "rand":
            shards = place_rand(n, pos_all, rng)
        elif strategy == "disp":
            shards = place_disp(n, pos_all)
        elif strategy == "clust":
            shards = place_clust(n, k, pos_all, rng, C=4, min_spacing=-1)
        else:  # clustr2
            C = max(2, n // k)
            shards = place_clust(n, k, pos_all, rng, C=C, min_spacing=2 * r0)
        disks = disaster(pos_all, r0, rng)
        alive = alive_mask(pos_all, disks, r0)
        th_n.append(sum(alive) / len(alive))
        th_a.append(theta_area_sample(disks, r0, rng))
        a_shards = sum(1 for s in shards if alive[s])
        if a_shards >= k:
            surv += 1
            comp, _ = components(alive, adj_all)
            cnt = {}
            hit = False
            for s in shards:
                if alive[s]:
                    cnt[comp[s]] = cnt.get(comp[s], 0) + 1
                    if cnt[comp[s]] >= k:
                        hit = True
                        break
            if hit:
                reach += 1
    return (surv / EPISODES, reach / EPISODES,
            statistics.mean(th_n), statistics.mean(th_a))


def main():
    global pos_all, adj_all, _BBOX
    os.makedirs("results", exist_ok=True)
    pos_all, adj_all, _BBOX = load_topo()
    W_km = _BBOX[1] - _BBOX[0]
    H_km = _BBOX[3] - _BBOX[2]
    print(">> Topology Zoo 拓扑: Kdl (Kentucky Datalink)")
    print("   节点=%d  边=%d  包围盒 %.0f x %.0f km"
          % (len(pos_all), sum(len(a) for a in adj_all) // 2, W_km, H_km))

    rows = []
    for r0 in R0S:
        for n in NS:
            for strat in STRATEGIES:
                surv, reach, th_n, th_a = run_cell(strat, n, r0, seed=17)
                rows.append({"r0": r0, "n": n, "strategy": strat,
                             "k": max(1, math.ceil(n * 0.5)),
                             "survival": round(surv, 4),
                             "reachable": round(reach, 4),
                             "theta_node": round(th_n, 4),
                             "theta_area": round(th_a, 4)})
                print("   r0=%3.0fkm n=%2d %-7s surv=%.3f reach=%.3f thN=%.3f thA=%.3f"
                      % (r0, n, strat, surv, reach, th_n, th_a))

    with open("results/zoo_results.csv", "w", newline="", encoding="utf-8") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows)

    # ---------------- 汇总
    with open("results/zoo_summary.txt", "w", encoding="utf-8") as f:
        f.write("SHARD Topology Zoo 验证（Kdl, Kentucky Datalink, %d 节点, %s）\n"
                % (len(pos_all), "x".join("%.0fkm" % v for v in (W_km, H_km))))
        f.write("灾害: %d 圆盘/轮, %d 轮/单元, k=ceil(n/2)\n\n"
                % (M_DISKS, EPISODES))
        f.write("== regime 实测 theta（节点存活率/面积存活率, 各 n 平均）==\n")
        for r0 in R0S:
            sub = [r for r in rows if r["r0"] == r0]
            f.write("r0=%.0fkm: th_node=%.3f  th_area=%.3f\n"
                    % (r0, statistics.mean(r["theta_node"] for r in sub),
                       statistics.mean(r["theta_area"] for r in sub)))
        f.write("\n== n=16 关键对照（可达率）==\n")
        f.write("%-8s" % "r0(km)")
        for s in STRATEGIES:
            f.write(" %-9s" % s)
        f.write("\n")
        for r0 in R0S:
            f.write("%-8.0f" % r0)
            for s in STRATEGIES:
                m = [r for r in rows if r["r0"] == r0 and r["n"] == 16 and r["strategy"] == s][0]
                f.write(" %-9.3f" % m["reachable"])
            f.write("\n")
        f.write("\n== n=16 生存-可达 gap（生存率 - 可达率）==\n")
        for r0 in R0S:
            for s in STRATEGIES:
                m = [r for r in rows if r["r0"] == r0 and r["n"] == 16 and r["strategy"] == s][0]
                f.write("r0=%.0f %-7s gap=%.3f\n" % (r0, s, m["survival"] - m["reachable"]))

    # ---------------- 出图：3 面板（每 regime 一panel），可达率 vs n，4 臂
    PW, PH, ML, MT, GAP = 220, 190, 46, 40, 18
    W = ML * 2 + PW * 3 + GAP * 2
    H = MT + PH + 52
    COL = {"rand": "#4477aa", "disp": "#ee6677", "clust": "#aa3377", "clustr2": "#228833"}
    ymax = max(r["reachable"] for r in rows) or 1.0
    ymax = math.ceil(max(ymax, 0.6) * 10) / 10
    svg = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" '
           'font-family="Helvetica,Arial,sans-serif">' % (W, H),
           '<rect width="%d" height="%d" fill="#ffffff"/>' % (W, H),
           '<text x="%d" y="24" font-size="15" font-weight="bold" fill="#222">'
           'SHARD on real ISP topology (Kdl): reachable service vs fragment count'
           '</text>' % (W // 2)]
    for pi, r0 in enumerate(R0S):
        x0 = ML + pi * (PW + GAP)
        svg.append('<rect x="%d" y="%d" width="%d" height="%d" fill="none" '
                   'stroke="#999"/>' % (x0, MT, PW, PH))
        for gy in range(0, 6):
            y = MT + PH - gy * PH / 5
            svg.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#e6e6e6"/>'
                       % (x0, y, x0 + PW, y))
            svg.append('<text x="%d" y="%.1f" font-size="10" fill="#666" '
                       'text-anchor="end">%.1f</text>'
                       % (x0 - 4, y + 3, gy * ymax / 5))
        sub_th = statistics.mean(r["theta_area"] for r in rows if r["r0"] == r0)
        svg.append('<text x="%d" y="%d" font-size="12" font-weight="bold" fill="#222" '
                   'text-anchor="middle">r0=%.0f km (E[#]: %.2f)</text>'
                   % (x0 + PW // 2, MT - 12, r0, sub_th))
        for s in STRATEGIES:
            pts = []
            for n in NS:
                m = [r for r in rows if r["r0"] == r0 and r["n"] == n and r["strategy"] == s][0]
                px = x0 + (NS.index(n)) * (PW / (len(NS) - 1))
                py = MT + PH - m["reachable"] / ymax * PH
                pts.append("%.1f,%.1f" % (px, py))
            svg.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="2"/>'
                       % (" ".join(pts), COL[s]))
        svg.append('<text x="%d" y="%d" font-size="11" fill="#444" '
                   'text-anchor="middle">fragment count n</text>'
                   % (x0 + PW // 2, MT + PH + 18))
    lx = ML
    for s in STRATEGIES:
        label = {"rand": "random", "disp": "dispersion", "clust": "naive clustered",
                 "clustr2": "R2-clustered"}[s]
        svg.append('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="%s" stroke-width="3"/>'
                   % (lx, H - 22, lx + 18, H - 22, COL[s]))
        svg.append('<text x="%d" y="%d" font-size="11" fill="#333">%s</text>'
                   % (lx + 24, H - 18, label))
        lx += 130
    svg.append('<text x="14" y="%d" font-size="11" fill="#666">reachable service</text>'
               % (MT + PH // 2))
    svg.append('</svg>')
    with open("results/fig_zoo.svg", "w", encoding="utf-8") as f:
        f.write("\n".join(svg))

    print("\n>> 完成: results/zoo_results.csv / zoo_summary.txt / fig_zoo.svg")


if __name__ == "__main__":
    main()
