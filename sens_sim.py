"""
SHARD 敏感性扫描（论文 §6.5，TODO ④）
====================================================
三部分：
A. 相图（定理 1 的 regime 分类验证）：
   在 (码率 c, 灾害强度 m) 网格上，用随机 iid 放置系综（与混合极限定理
   假设严格一致）测 A(4) 与 A(64)，同时经验估计 θ̄ 与 P(θ≥c)。
   预测：(i) A(64) ≈ P(θ≥c)（定理极限，42 个 cell 全部对照）；
        (ii) 冗余收益 benefit=A(64)-A(4) 的符号在 θ̄=c 对角线附近翻转
             （E[θ] > c → 加碎片有益；E[θ] < c → 有害）。
   这是全文理论-实验闭环最硬的一张图：regime 边界由参数扫描实测。

B. 码率 × 可达性耦合（连通层）：
   c ∈ {1/3, 1/2, 2/3}，L=2.5，n=16，四臂（复用 reach_sim 机制）。
   验证 §6.3 的设计空间推论：R2 强制簇数 C=⌊n/k⌋ 随 c 增大而坍缩
   （c=2/3 时 C=1，分簇退化为单自足区域——码率与可部署几何耦合）。

C. 节点密度 ρ 扫描（连通层）：
   N ∈ {50, 80, 120}，L=2.5，n=16，k=8，四臂。
   验证：生存-可达 gap 随密度上升而收窄（分区惩罚由连通性主导）。

运行（在 sim/ 目录下）：python sens_sim.py   （纯标准库，复用 reach_sim）
输出：results/sens_phase.csv、results/sens_reach.csv、results/sens_summary.txt
      results/fig_sens_phase.svg、results/fig_sens_reach.svg
对应论文：§6.5（Fig. 9 相图 + Fig. 10 码率/密度扫描）
"""

import math
import random
import csv
import os
import statistics

import reach_sim as RS  # 复用连通层机制（build_graph / 四臂放置 / 毁伤）

# ---------------- Part A 参数（相图） ----------------
AREA = 10.0
R0_A = 2.5                       # 固定灾害半径，用 m 扫灾害强度
CS_A = (0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)   # 码率网格
MS_A = (1, 2, 3, 4, 6, 8)        # 每轮毁伤圆盘数网格
EPS_A = 500                      # 每 cell 蒙特卡洛轮数
NTH = 250                        # θ 面积采样点数
NS_A = (4, 64)                   # 小/大碎片数

# ---------------- Part B/C 参数（连通层） ----------------
CS_R = (1/3, 1/2, 2/3)           # 码率扫描（B）
NS_R = (50, 80, 120)             # 节点密度扫描（C）
L_BC = 2.5
N_BC = 16
GSEEDS = (1, 2, 3, 4, 5)          # 与 reach_sim 主实验一致（5 图种子）
EPS_R = 600                      # 与 reach_sim 主实验一致（600 轮）
STRATEGIES = RS.STRATEGIES


# ================================================================ Part A
def run_phase_cell(c, m, seed):
    """随机 iid 放置系综（定理 1 的假设）：返回 A(4), A(64), theta_bar, P(theta>=c)"""
    rng = random.Random(seed * 131 + int(c * 100) * 7 + m * 1009)
    r2 = R0_A * R0_A
    k4 = max(1, math.ceil(NS_A[0] * c))
    k64 = max(1, math.ceil(NS_A[1] * c))
    ok4 = ok64 = 0
    thetas = []
    for _ in range(EPS_A):
        disks = [(rng.uniform(0, AREA), rng.uniform(0, AREA)) for _ in range(m)]
        # 68 个 iid 碎片点：前 4 个作 n=4 臂，后 64 个作 n=64 臂
        alive = []
        for _ in range(68):
            px, py = rng.uniform(0, AREA), rng.uniform(0, AREA)
            dead = False
            for dx0, dy0 in disks:
                dx, dy = px - dx0, py - dy0
                if dx * dx + dy * dy <= r2:
                    dead = True
                    break
            alive.append(not dead)
        # θ：面积采样（定理的混合变量）
        okth = 0
        for _ in range(NTH):
            px, py = rng.uniform(0, AREA), rng.uniform(0, AREA)
            dead = False
            for dx0, dy0 in disks:
                dx, dy = px - dx0, py - dy0
                if dx * dx + dy * dy <= r2:
                    dead = True
                    break
            if not dead:
                okth += 1
        thetas.append(okth / NTH)
        if sum(alive[:4]) >= k4:
            ok4 += 1
        if sum(alive[4:]) >= k64:
            ok64 += 1
    a4, a64 = ok4 / EPS_A, ok64 / EPS_A
    pth = sum(1 for t in thetas if t >= c - 1e-9) / len(thetas)
    return a4, a64, statistics.mean(thetas), pth


def part_a():
    print(">> Part A: 相图（%d 个码率 x %d 个灾害强度 = %d cell，每 cell %d 轮）"
          % (len(CS_A), len(MS_A), len(CS_A) * len(MS_A), EPS_A))
    rows = []
    for m in MS_A:
        for c in CS_A:
            a4, a64, th, pth = run_phase_cell(c, m, seed=23)
            ben = a64 - a4
            rows.append({"c": c, "m": m, "theta_bar": round(th, 4),
                         "P_theta_c": round(pth, 4), "A4": round(a4, 4),
                         "A64": round(a64, 4), "benefit": round(ben, 4),
                         "regime_pred": "beneficial" if th > c else "harmful"})
            print("   c=%.2f m=%d  th=%.3f  P(th>=c)=%.3f  A4=%.3f  A64=%.3f  ben=%+.3f  [%s]"
                  % (c, m, th, pth, a4, a64, ben, rows[-1]["regime_pred"]))
    with open("results/sens_phase.csv", "w", newline="", encoding="utf-8") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows)
    return rows


# ================================================================ Part B/C
def reach_cell(strategy, n, r0, k_target, n_nodes, seed, episodes):
    """monkey-patch reach_sim 的码率/节点数/轮数后调用其 run_cell"""
    RS.N_NODES = n_nodes
    RS.EPISODES = episodes
    RS.k_of = (lambda kk: lambda x: kk)(k_target)   # 固定码率的 k
    return RS.run_cell(strategy, n, r0, seed)


def sweep_bc():
    print(">> Part B: 码率扫描（c ∈ {1/3,1/2,2/3}，L=%.1f，n=%d，%d 图种子 x %d 轮）"
          % (L_BC, N_BC, len(GSEEDS), EPS_R))
    rows_b = []
    for c in CS_R:
        k = max(1, math.ceil(N_BC * c))
        C_r2 = max(1, N_BC // k)
        for strat in STRATEGIES:
            sv, rc = [], []
            for gs in GSEEDS:
                s, r = reach_cell(strat, N_BC, L_BC, k, 80, gs, EPS_R)
                sv.append(s)
                rc.append(r)
            rows_b.append({"sweep": "code_rate", "param": round(c, 4), "k": k,
                           "C_r2": C_r2, "strategy": strat,
                           "survival": round(statistics.mean(sv), 4),
                           "reachable": round(statistics.mean(rc), 4)})
            print("   c=%.2f k=%2d C=%d %-7s surv=%.3f reach=%.3f"
                  % (c, k, C_r2, strat, rows_b[-1]["survival"], rows_b[-1]["reachable"]))
    print(">> Part C: 节点密度扫描（N ∈ %s，L=%.1f，n=%d，k=8）"
          % (list(NS_R), L_BC, N_BC))
    rows_c = []
    for nn in NS_R:
        for strat in STRATEGIES:
            sv, rc = [], []
            for gs in GSEEDS:
                s, r = reach_cell(strat, N_BC, L_BC, 8, nn, gs, EPS_R)
                sv.append(s)
                rc.append(r)
            rows_c.append({"sweep": "node_density", "param": nn, "k": 8,
                           "C_r2": 2, "strategy": strat,
                           "survival": round(statistics.mean(sv), 4),
                           "reachable": round(statistics.mean(rc), 4)})
            print("   N=%3d %-7s surv=%.3f reach=%.3f"
                  % (nn, strat, rows_c[-1]["survival"], rows_c[-1]["reachable"]))
    rows = rows_b + rows_c
    with open("results/sens_reach.csv", "w", newline="", encoding="utf-8") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows)
    return rows_b, rows_c


# ================================================================ 出图
def fig_phase(rows):
    """两面板：左=(c, θ̄) 相图按 benefit 着色 + 对角线；右=A(64) vs P(θ≥c) 定理验证"""
    W, H = 980, 400
    PW1, PW2 = 420, 400
    ML, MT, GAP = 70, 70, 60
    CH = H - MT - 50
    # benefit 颜色：正=绿系，负=红系，零=浅灰
    bmax = max(abs(r["benefit"]) for r in rows) or 1.0

    def bcol(b):
        t = b / bmax
        if t >= 0:
            v = 40 + int(165 * t)
            return "rgb(%d,%d,%d)" % (110 - int(70 * t), v, 130 - int(60 * t))
        v = 40 + int(165 * (-t))
        return "rgb(%d,%d,%d)" % (150 + int(70 * (-t)), v, 100)

    s = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" '
         'font-family="Helvetica,Arial,sans-serif">' % (W, H),
         '<rect width="%d" height="%d" fill="#ffffff"/>' % (W, H),
         '<text x="%d" y="26" font-size="15" font-weight="bold" fill="#222">'
         'Sensitivity: regime phase diagram and mixture-limit check</text>' % (W // 2)]
    # -------- 左：相图
    x0, y0 = ML, MT
    cmin, cmax = min(CS_A), max(CS_A)
    th_lo = min(r["theta_bar"] for r in rows) - 0.05
    th_hi = max(r["theta_bar"] for r in rows) + 0.05

    def px(c): return x0 + (c - cmin) / (cmax - cmin) * PW1

    def py(th): return y0 + CH * (1 - (th - th_lo) / (th_hi - th_lo))

    for gy in range(6):
        yv = th_lo + gy * (th_hi - th_lo) / 5
        s.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#ececec"/>'
                 % (x0, py(yv), x0 + PW1, py(yv)))
        s.append('<text x="%d" y="%.1f" font-size="10" fill="#777" text-anchor="end">%.2f</text>'
                 % (x0 - 6, py(yv) + 3, yv))
    for cx in range(0, 7):
        c = cmin + cx * (cmax - cmin) / 6
        s.append('<text x="%.1f" y="%d" font-size="10" fill="#777" text-anchor="middle">%.1f</text>'
                 % (px(c), y0 + CH + 16, c))
    # cell 方块
    cw = PW1 / (len(CS_A) - 1)
    for r in rows:
        cx, cy = px(r["c"]), py(r["theta_bar"])
        s.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s" '
                 'stroke="#fff" stroke-width="1" opacity="0.92"/>'
                 % (cx - cw / 2, cy - cw / 2, cw, cw, bcol(r["benefit"])))
    # 对角线 θ̄=c（裁剪到绘图区）
    pts = []
    for c in (cmin, cmax):
        if th_lo <= c <= th_hi:
            pts.append("%.1f,%.1f" % (px(c), py(c)))
    if len(pts) == 2:
        s.append('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="#333" '
                 'stroke-width="1.6" stroke-dasharray="6 4"/>' % tuple(
                     [float(v.split(",")[0]) for v in pts] + [float(v.split(",")[1]) for v in pts]))
    s.append('<text x="%d" y="%d" font-size="12" fill="#222" text-anchor="middle">'
             'code rate c</text>' % (x0 + PW1 // 2, y0 + CH + 34))
    s.append('<text x="18" y="%d" font-size="11" fill="#666" '
             'transform="rotate(-90 24 %d)">E[theta]</text>' % (y0 + CH // 2, y0 + CH // 2))
    s.append('<text x="%d" y="%d" font-size="12" font-weight="bold" fill="#222">Phase diagram: A(64)-A(4)</text>'
             % (x0 + 10, y0 - 14))
    s.append('<text x="%d" y="%d" font-size="10" fill="#555">dashed: E[theta] = c</text>'
             % (x0 + 10, y0 - 2))
    # -------- 右：定理验证散点
    x1 = ML + PW1 + GAP
    s.append('<rect x="%d" y="%d" width="%d" height="%d" fill="none" stroke="#999"/>'
             % (x1, y0, PW2, CH))
    for gy in range(6):
        yv = gy / 5
        yy = y0 + CH * (1 - yv)
        s.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#ececec"/>'
                 % (x1, yy, x1 + PW2, yy))
        s.append('<text x="%d" y="%.1f" font-size="10" fill="#777" text-anchor="end">%.2f</text>'
                 % (x1 - 6, yy + 3, yv))
    s.append('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#333" stroke-width="1.4"/>'
             % (x1, y0 + CH, x1 + PW2, y0))
    for r in rows:
        xx = x1 + r["P_theta_c"] / 1.0 * PW2
        yy = y0 + CH * (1 - r["A64"])
        s.append('<circle cx="%.1f" cy="%.1f" r="4.5" fill="%s" fill-opacity="0.85"/>'
                 % (xx, yy, bcol(r["benefit"])))
    s.append('<text x="%d" y="%d" font-size="10" fill="#777" text-anchor="middle">P(theta &gt;= c) (measured)</text>'
             % (x1 + PW2 // 2, y0 + CH + 16))
    s.append('<text x="%d" y="%d" font-size="12" font-weight="bold" fill="#222">'
             'Mixture-limit check: A(64) vs P(theta&gt;=c)</text>' % (x1 + 10, y0 - 14))
    s.append('<text x="%d" y="%d" font-size="10" fill="#555">solid: y = x (Theorem 1)</text>'
             % (x1 + 10, y0 - 2))
    s.append('<text x="14" y="%d" font-size="11" fill="#666" '
             'transform="rotate(-90 %d %d)">A(64)</text>' % (y0 + CH // 2, x1 - 46, y0 + CH // 2))
    # 图例
    lx = ML
    for lab, col in (("harmful (ben&lt;0)", bcol(-bmax)),
                     ("neutral", bcol(0)),
                     ("beneficial (ben&gt;0)", bcol(bmax))):
        s.append('<rect x="%d" y="%d" width="14" height="10" fill="%s"/>' % (lx, H - 26, col))
        s.append('<text x="%d" y="%d" font-size="11" fill="#333">%s</text>' % (lx + 18, H - 17, lab))
        lx += 170
    s.append('</svg>')
    with open("results/fig_sens_phase.svg", "w", encoding="utf-8") as f:
        f.write("\n".join(s))


def fig_reach(rows_b, rows_c):
    """两面板柱状图：左=码率扫描（标 k 与 C），右=节点密度扫描"""
    W, H = 980, 400
    PW1, PW2 = 400, 400
    ML, MT, GAP = 60, 70, 60
    CH = H - MT - 60
    COL = {"rand": "#4477aa", "disp": "#ee6677", "clust": "#aa3377", "clustr2": "#228833"}
    s = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" '
         'font-family="Helvetica,Arial,sans-serif">' % (W, H),
         '<rect width="%d" height="%d" fill="#ffffff"/>' % (W, H),
         '<text x="%d" y="26" font-size="15" font-weight="bold" fill="#222">'
         'Sensitivity in the connectivity layer (L=2.5, n=16)</text>' % (W // 2)]

    def panel(x0, groups, rows, xlabel, title, subfmt):
        for gy in range(6):
            yv = gy / 5
            yy = MT + CH * (1 - yv)
            s.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#ececec"/>'
                     % (x0, yy, x0 + PW1, yy))
            s.append('<text x="%d" y="%.1f" font-size="10" fill="#777" text-anchor="end">%.2f</text>'
                     % (x0 - 6, yy + 3, yv))
        gw = PW1 / len(groups)
        bw = gw * 0.18
        for gi, g in enumerate(groups):
            cell = [r for r in rows if r["param"] == g]
            for si, strat in enumerate(STRATEGIES):
                r = [x for x in cell if x["strategy"] == strat][0]
                bx = x0 + gi * gw + gw / 2 - 2 * bw + si * (bw + 2)
                hgt = CH * r["reachable"]
                s.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s"/>'
                         % (bx, MT + CH - hgt, bw, hgt, COL[strat]))
                s.append('<text x="%.1f" y="%.1f" font-size="8" text-anchor="middle" fill="#222">%.2f</text>'
                         % (bx + bw / 2, MT + CH - hgt - 3, r["reachable"]))
            s.append('<text x="%.1f" y="%d" font-size="11" fill="#222" text-anchor="middle">%s</text>'
                     % (x0 + gi * gw + gw / 2, MT + CH + 16, subfmt(g)))
        s.append('<text x="%d" y="%d" font-size="12" fill="#222" text-anchor="middle">%s</text>'
                 % (x0 + PW1 // 2, MT + CH + 36, xlabel))
        s.append('<text x="%d" y="%d" font-size="12" font-weight="bold" fill="#222">%s</text>'
                 % (x0 + 10, MT - 16, title))

    panel(ML, sorted(set(r["param"] for r in rows_b)), rows_b, "code rate c",
          "Code-rate sweep (k-anycast reachability)",
          lambda c: "c=%.2f (k=%d, C=%d)" % (c, max(1, math.ceil(N_BC * c)),
                                             max(1, N_BC // max(1, math.ceil(N_BC * c)))))
    panel(ML + PW1 + GAP, sorted(set(r["param"] for r in rows_c)), rows_c,
          "node count N", "Node-density sweep (k=8)", lambda n: "N=%d" % n)
    lx = ML
    for strat in STRATEGIES:
        lab = {"rand": "random", "disp": "dispersion", "clust": "naive clustered",
               "clustr2": "R2-clustered"}[strat]
        s.append('<rect x="%d" y="%d" width="14" height="10" fill="%s"/>' % (lx, H - 26, COL[strat]))
        s.append('<text x="%d" y="%d" font-size="11" fill="#333">%s</text>' % (lx + 18, H - 17, lab))
        lx += 150
    s.append('</svg>')
    with open("results/fig_sens_reach.svg", "w", encoding="utf-8") as f:
        f.write("\n".join(s))


# ================================================================ 汇总
def write_summary(rows_a, rows_b, rows_c):
    with open("results/sens_summary.txt", "w", encoding="utf-8") as f:
        f.write("SHARD 敏感性扫描（%d 轮/cell，种子固定）\n" % EPS_A)
        f.write("=" * 64 + "\n\n")
        f.write("[A] 相图：regime 符号判定 vs 定理预测（E[θ] vs c）\n")
        agree = sum(1 for r in rows_a
                    if (r["benefit"] > 0) == (r["theta_bar"] > r["c"]))
        f.write("   符号一致率： %d / %d cells\n" % (agree, len(rows_a)))
        for r in rows_a:
            f.write("   c=%.2f m=%d  E[th]=%.3f  A4=%.3f  A64=%.3f  ben=%+.3f  pred=%s\n"
                    % (r["c"], r["m"], r["theta_bar"], r["A4"], r["A64"],
                       r["benefit"], r["regime_pred"]))
        f.write("\n[A] 定理极限对照：A(64) vs 经验 P(θ≥c)\n")
        err = [abs(r["A64"] - r["P_theta_c"]) for r in rows_a]
        f.write("   平均 |A64 - P(θ≥c)| = %.4f  最大 = %.4f\n"
                % (statistics.mean(err), max(err)))
        f.write("\n[B] 码率扫描（L=2.5, n=16, 可达率）\n")
        for r in rows_b:
            f.write("   c=%.2f k=%2d C=%d %-7s surv=%.3f reach=%.3f\n"
                    % (r["param"], r["k"], r["C_r2"], r["strategy"],
                       r["survival"], r["reachable"]))
        f.write("\n[C] 节点密度扫描（L=2.5, n=16, k=8, 可达率）\n")
        for r in rows_c:
            f.write("   N=%3d %-7s surv=%.3f reach=%.3f gap=%.3f\n"
                    % (r["param"], r["strategy"], r["survival"],
                       r["reachable"], r["survival"] - r["reachable"]))


# ================================================================ main
def main():
    os.makedirs("results", exist_ok=True)
    rows_a = part_a()
    rows_b, rows_c = sweep_bc()
    write_summary(rows_a, rows_b, rows_c)
    fig_phase(rows_a)
    fig_reach(rows_b, rows_c)
    print("\n[ok] results/sens_phase.csv / sens_reach.csv / sens_summary.txt / "
          "fig_sens_phase.svg / fig_sens_reach.svg")


if __name__ == "__main__":
    main()
