"""
GDDA 原型仿真器 —— 应急网云系统「部署梯度」架构
纯标准库实现，无第三方依赖，可直接运行与开源。

用途：验证「局部性」与「契约保持」两条性质，生成论文主图数据。

运行：python gdda_sim.py
输出：results/gdda_results.csv
      results/summary.txt
      results/stability.csv
      results/fig_locality.svg
      results/fig_sla.svg
"""

import math
import random
import csv
import os
import statistics

AREA = 10.0
N_DISTRICT = 6
PER_DISTRICT = 10
COMM_RANGE = 3.0
NODE_CAP = 40
MAX_HOPS = 3

# 梯度阈值：H 低于阈值即降级（L1/L2/L3 的门限）
TIER_TH = [0.92, 0.82, 0.72]


def build_world(seed):
    rng = random.Random(seed)
    centers = []
    for i in range(3):
        for j in range(2):
            centers.append((1.8 + i * 3.2 + rng.uniform(-0.4, 0.4),
                            2.6 + j * 4.8 + rng.uniform(-0.4, 0.4)))

    coords, district = [], []
    for d, (cx, cy) in enumerate(centers):
        for _ in range(PER_DISTRICT):
            x = min(AREA - 0.2, max(0.2, rng.gauss(cx, 0.8)))
            y = min(AREA - 0.2, max(0.2, rng.gauss(cy, 0.8)))
            coords.append((x, y))
            district.append(d)

    n = len(coords)
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if math.hypot(coords[i][0] - coords[j][0],
                          coords[i][1] - coords[j][1]) <= COMM_RANGE:
                adj[i].append(j)
                adj[j].append(i)

    deg0 = sum(len(a) for a in adj) / n
    cap0 = n * NODE_CAP
    load = [rng.uniform(50, 105) for _ in range(N_DISTRICT)]
    return {"coords": coords, "district": district, "adj": adj,
            "load": load, "centers": centers, "deg0": deg0, "cap0": cap0,
            "anchor0": [d * PER_DISTRICT for d in range(N_DISTRICT)]}


def damage(world, radius, seed, epicenter=None):
    rng = random.Random(seed * 7919 + int(radius * 100))
    if epicenter is None:
        ex = rng.uniform(1.5, AREA - 1.5)
        ey = rng.uniform(1.5, AREA - 1.5)
    else:
        ex, ey = epicenter
    alive, hot = set(), set()
    for i, (x, y) in enumerate(world["coords"]):
        r = math.hypot(x - ex, y - ey)
        p_kill = min(1.0, max(0.0, (1.0 - r / radius)) * 1.6) if radius > 0 else 0.0
        if rng.random() >= p_kill:
            alive.add(i)
        else:
            hot.add(i)
    return alive, hot


def components(adj, alive):
    seen, comps = set(), []
    for s in sorted(alive):
        if s in seen:
            continue
        stack, comp = [s], []
        seen.add(s)
        while stack:
            u = stack.pop()
            comp.append(u)
            for v in adj[u]:
                if v in alive and v not in seen:
                    seen.add(v)
                    stack.append(v)
        comps.append(comp)
    return comps


def bfs(adj, alive, src, max_hops=MAX_HOPS):
    if src not in alive:
        return {}
    dist = {src: 0}
    frontier = [src]
    for h in range(1, max_hops + 1):
        nxt = []
        for u in frontier:
            for v in adj[u]:
                if v in alive and v not in dist:
                    dist[v] = h
                    nxt.append(v)
        frontier = nxt
        if not frontier:
            break
    return dist


def health(adj, alive, cap0, deg0):
    """健康度 H = 0.45·相对度 + 0.40·相对容量 + 0.15·巨型连通分量占比。
    正式版第一项替换为归一化代数连通度 λ2（接口预留）。"""
    if not alive:
        return 0.0
    comps = components(adj, alive)
    giant = max(len(c) for c in comps)
    n = sum(len(a) for a in adj) and len(adj)
    deg = sum(1 for u in alive for v in adj[u] if v in alive) / max(len(alive), 1)
    h = 0.45 * min(1.0, deg / max(deg0, 1e-6))
    h += 0.40 * min(1.0, (len(alive) * NODE_CAP) / max(cap0, 1e-6))
    h += 0.15 * (giant / max(len(adj), 1))
    return h


def tier_of(h, prev=None, hyst=0.0):
    base = 0
    for i, th in enumerate(TIER_TH):
        if h < th:
            base = i + 1
    if prev is None or hyst <= 0:
        return base
    if base == prev:
        return prev
    if base > prev:
        th = TIER_TH[min(prev, len(TIER_TH) - 1)]
        return prev if h > th - hyst else base
    th = TIER_TH[max(0, min(base, len(TIER_TH) - 1))]
    return prev if h < th + hyst else base


def serve(world, alive, scheme, hot, tier):
    """返回 (需同步的节点数, 各区服务率, 各区时延违规率)"""
    adj, district, centers, load = (world["adj"], world["district"],
                                    world["centers"], world["load"])
    servo, viol = [], []

    if scheme == "static":
        sync = 0
        anchor = world["anchor0"]
    elif scheme == "gdda":
        sync = len(hot) + len({v for u in hot for v in adj[u] if v in alive})
        anchor = []
        for d in range(N_DISTRICT):
            cand = [i for i in district_alive(district, alive, d)]
            if not cand:
                anchor.append(None)
            else:
                anchor.append(max(cand, key=lambda i: (len(adj[i]), -i)))
    else:  # central：全局搜索，且需全网状态上报（含失效检测）
        sync = len(alive) + len(hot)
        anchor = []
        for d in range(N_DISTRICT):
            best, bestd = None, 99
            cx, cy = centers[d]
            for i in alive:
                dd = math.hypot(world["coords"][i][0] - cx,
                                world["coords"][i][1] - cy)
                if dd < bestd:
                    best, bestd = i, dd
            anchor.append(best)

    for d in range(N_DISTRICT):
        a = anchor[d]
        if a is None or a not in alive:
            servo.append(0.0)
            viol.append(1.0)
            continue
        dist = bfs(adj, alive, a)
        reach = sum(1 for v in dist if district[v] == d)
        # gdda/static 只能用本区节点；central 可用任意可达节点
        if scheme == "central":
            reach = len(dist)
        cap = reach * NODE_CAP
        servo.append(min(1.0, cap / max(load[d], 1e-6)))
        hops = dist.get(d * PER_DISTRICT, MAX_HOPS)
        viol.append(0.0 if hops <= MAX_HOPS else 1.0)
    return sync, servo, viol


def district_alive(district, alive, d):
    return [i for i in range(len(district)) if district[i] == d and i in alive]


def run_once(seed, radius):
    w = build_world(seed)
    alive, hot = damage(w, radius, seed)
    h = health(w["adj"], alive, w["cap0"], w["deg0"])
    tier = tier_of(h)

    row = {"seed": seed, "radius": radius, "alive": len(alive),
           "comps": len(components(w["adj"], alive)),
           "health": round(h, 4), "tier": tier}

    for sc in ("gdda", "central", "static"):
        sync, servo, viol = serve(w, alive, sc, hot, tier)
        row[sc + "_sync"] = sync
        row[sc + "_sla"] = round(statistics.mean(
            s * (1 - v) for s, v in zip(servo, viol)), 4)
        row[sc + "_serve"] = round(statistics.mean(servo), 4)
    return row


def experiment():
    radii = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
    seeds = list(range(1, 31))
    rows = [run_once(s, r) for r in radii for s in seeds]

    os.makedirs("results", exist_ok=True)
    with open("results/gdda_results.csv", "w", newline="", encoding="utf-8") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows)

    agg = {}
    keys = ("alive", "comps", "health", "tier", "gdda_sync", "gdda_sla",
            "central_sync", "central_sla", "static_sync", "static_sla")
    for r in radii:
        sub = [x for x in rows if x["radius"] == r]
        agg[r] = {k: statistics.mean(x[k] for x in sub) for k in keys}
        agg[r]["gdda_sd"] = statistics.stdev(x["gdda_sync"] for x in sub)
        agg[r]["sla_sd"] = statistics.stdev(x["gdda_sla"] for x in sub)

    lines = ["%-7s %-6s %-6s %-7s %-4s %-14s %-14s %-14s" %
             ("radius", "alive", "comps", "health", "tier",
              "GDDA msg/SLA", "CENT msg/SLA", "STATIC msg/SLA")]
    for r in radii:
        a = agg[r]
        lines.append("%-7.1f %-6.1f %-6.1f %-7.3f %-4d %-14s %-14s %-14s" % (
            r, a["alive"], a["comps"], a["health"], int(a["tier"]),
            "%.1f/%.3f" % (a["gdda_sync"], a["gdda_sla"]),
            "%.1f/%.3f" % (a["central_sync"], a["central_sla"]),
            "%.1f/%.3f" % (a["static_sync"], a["static_sla"])))
    report = "\n".join(lines)
    with open("results/summary.txt", "w", encoding="utf-8") as f:
        f.write(report + "\n")

    with open("results/stability.csv", "w", newline="", encoding="utf-8") as f:
        wr = csv.writer(f)
        wr.writerow(["hyst", "switches_mean", "switches_sd"])
        for hh in (0.0, 0.02, 0.04, 0.06, 0.10, 0.15):
            sw = []
            for s in seeds:
                w = build_world(s)
                rng = random.Random(s * 104729)
                epi = (rng.uniform(3.0, 7.0), rng.uniform(3.0, 7.0))
                prev, cnt = None, 0
                for r in [3.0, 3.8] * 6:
                    alive, hot = damage(w, r, s, epicenter=epi)
                    h = health(w["adj"], alive, w["cap0"], w["deg0"])
                    t = tier_of(h, prev, hh)
                    if prev is not None and t != prev:
                        cnt += 1
                    prev = t
                sw.append(cnt)
            wr.writerow([hh, round(statistics.mean(sw), 2),
                         round(statistics.stdev(sw), 2) if len(sw) > 1 else 0])

    write_svg("results/fig_locality.svg", radii, agg,
              [("gdda_sync", "#1D9E75", "GDDA  O(|N[R]|)"),
               ("central_sync", "#C0392B", "Centralized  O(|V|)")],
              "Control messages per damage event", "Messages")
    write_svg("results/fig_sla.svg", radii, agg,
              [("gdda_sla", "#1D9E75", "GDDA"),
               ("central_sla", "#C0392B", "Centralized"),
               ("static_sla", "#7F77DD", "Static 3-tier")],
              "SLA satisfaction ratio", "SLA ratio", ymax=1.0)

    print(report)
    print("\n[ok] results/")


def write_svg(path, radii, agg, series, title, ylabel, ymax=None):
    W, H, L, R, T, B = 700, 380, 70, 190, 44, 58
    pw, ph = W - L - R, H - T - B
    top = ymax if ymax else max(max(agg[r][k] for r in radii) for k, _, _ in series)
    top = max(top, 1e-6)

    def X(i):
        return L + pw * i / (len(radii) - 1)

    def Y(v):
        return T + ph * (1 - min(v, top) / top)

    s = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" width="100%%">' % (W, H),
         '<rect width="%d" height="%d" fill="#ffffff"/>' % (W, H),
         '<text x="%d" y="24" font-size="14" font-family="sans-serif" fill="#2c2c2a">%s</text>' % (L, title)]
    for k in range(5):
        yv = top * k / 4
        y = Y(yv)
        s.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#e8e8e3" stroke-width="0.6"/>' % (L, y, L + pw, y))
        s.append('<text x="%d" y="%.1f" font-size="11" text-anchor="end" font-family="sans-serif" fill="#5f5e5a">%s</text>'
                 % (L - 8, y + 4, ("%.0f" % yv) if top > 3 else ("%.1f" % yv)))
    for i, r in enumerate(radii):
        s.append('<text x="%.1f" y="%d" font-size="11" text-anchor="middle" font-family="sans-serif" fill="#5f5e5a">%.1f</text>'
                 % (X(i), T + ph + 18, r))
    s.append('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#9a9a90" stroke-width="1"/>' % (L, T + ph, L + pw, T + ph))
    s.append('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#9a9a90" stroke-width="1"/>' % (L, T, L, T + ph))
    for k, color, _ in series:
        pts = " ".join("%.1f,%.1f" % (X(i), Y(agg[r][k])) for i, r in enumerate(radii))
        s.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="2"/>' % (pts, color))
        for i, r in enumerate(radii):
            s.append('<circle cx="%.1f" cy="%.1f" r="3" fill="%s"/>' % (X(i), Y(agg[r][k]), color))
    s.append('<text x="%d" y="%d" font-size="12" font-family="sans-serif" fill="#5f5e5a" text-anchor="middle">Damage radius (km)</text>'
             % (L + pw // 2, H - 20))
    s.append('<text x="16" y="%d" font-size="12" font-family="sans-serif" fill="#5f5e5a" transform="rotate(-90 16 %d)" text-anchor="middle">%s</text>'
             % (T + ph // 2, T + ph // 2, ylabel))
    ly = T + 6
    for _, color, label in series:
        s.append('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="%s" stroke-width="2"/>' % (L + pw + 16, ly, L + pw + 40, ly, color))
        s.append('<text x="%d" y="%d" font-size="11" font-family="sans-serif" fill="#2c2c2a">%s</text>' % (L + pw + 46, ly + 4, label))
        ly += 20
    s.append('</svg>')
    with open(path, "w", encoding="utf-8") as f:
        f.write("".join(s))


if __name__ == "__main__":
    experiment()
