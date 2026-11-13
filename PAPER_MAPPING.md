# SHARD 仿真套件 — 使用说明

论文《When Redundancy Backfires: Correlation-Aware Service Fragmentation for
Post-Disaster Network-Cloud Systems》的全部实验代码。

## 快速开始

```bash
cd sim
python shard_sim_v2.py   # 主实验：三臂放置 × 三 regime + 混合极限交叉验证（约 2 分钟）
python reach_sim.py      # 连通层实验：生存 vs 可达（约 2 分钟）
```

依赖：**仅 Python 3 标准库**，无第三方依赖，可直接开源（投稿时作为 artifact）。
所有随机种子固定，结果可复现。

## 文件 → 论文图表映射

| 脚本 | 输出 | 论文位置 |
|---|---|---|
| shard_sim_v2.py | results/fig_regime_L15/25/35.svg | Fig. 4（三个 regime 的 A(n) 曲线） |
| shard_sim_v2.py | results/fig_theory.svg | Fig. 5（散布曲线收敛于定理极限 P(θ≥c)） |
| shard_sim_v2.py | results/fig_flip.svg | Fig. 6（n=48 几何反转柱状图） |
| shard_sim_v2.py | results/shard_v2_results.csv / summary.txt | Table 3 数据源 |
| reach_sim.py | results/fig_reach.svg | Fig. 7（生存 vs 可达缺口） |
| reach_sim.py | results/reach_results.csv / summary.txt | Table 4 数据源 |

## 核心结果速览（30/15 个种子、2026-11-12 运行）

**定理极限 P(θ≥c)**（c=1/2）：L=1.5 → 1.000；L=2.5 → 0.747；L=3.5 → 0.343

| 现象 | 数据 |
|---|---|
| 冗余反噬 | L=3.5 时散布放置：n=2 处 A=0.70，n=64 处 A=0.29 |
| 峰值后衰减 | L=2.5：峰值 0.90 (n=8) → 0.81 (n=64)，趋向极限 0.747 |
| 几何反转 | n=48 时：L=1.5 散布胜 (1.000 vs 0.994)；L=2.5 分簇胜 (0.821 vs 0.773)；L=3.5 分簇胜 (0.337 vs 0.308) |
| 可达缺口 | L=2.5, n=4, 散布：几何存活 0.839 → 可达服务仅 0.323（缺口 0.52） |
| 二次反转 | 连通层下随机放置可达率最高，散布最差（几何存活与网络可聚性互相冲突） |

## 已知诚实边界（写进论文 Discussion 的内容）

1. 区域内可行独立失效单元数有几何上限（本设定 C_max = 8/8/5），
   "簇数 = 1/ρ" 的上界受 ⌊Area/πL²⌋ 约束。
2. reach 实验中的朴素分簇臂（每簇碎片数 < k）违反基数规则 R2，
   其可达率崩塌正是 R2 必要性的实验证据，不是 bug。
3. 定理对"随机放置系综"严格成立；固定放置在有限 n 有波动
   （L=3.5 时散布 n=64 = 0.294 略低于极限 0.343，属有限 n 效应）。

## 待办（论文 TODO 清单对应）

- [ ] 接入 Topology Zoo 真实拓扑替换合成图（`build_graph` 一处改动）
- [ ] 置信区间：每单元已多种子，成稿前统一报 95% CI
- [ ] 分簇的修正臂（每簇 ≥ k 片，基数规则 R2 的对照验证）
- [ ] 参数敏感性：m（灾情强度）、c（码率）、节点密度扫描

## Legacy（已被取代，保留供参考）

- `gdda_sim.py` — GDDA 自适应控制原型。已知问题：全局健康度与分量级
  不变量不一致、baseline 不公平（V5–V7），修复前**不得**用于论文数据。
- `spread_sim.py` — v1 预实验（二臂），已被 shard_sim_v2.py 完全取代。
