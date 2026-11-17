# SHARD: Correlation-Aware Service Fragmentation for Geo-Distributed Edge–Cloud Infrastructures under Spatially Correlated Failures

Research code accompanying a manuscript under preparation (target: *Future
Generation Computer Systems*, Elsevier). All experiments are pure-Python,
standard library only, with fixed random seeds — fully reproducible on any
machine with **Python 3.8+** and **zero third-party dependencies**.

The thesis of this work: under spatially correlated hazards (clustered
destruction disks parameterized by a correlation length *L*), adding coded
redundancy to service deployment can *reduce* availability — and the optimal
placement geometry flips twice depending on whether you measure survival or
network-level reachability. All findings replicate on a real ISP topology
(Topology Zoo / Kdl).

## Quick start

```bash
python shard_sim_v2.py   # main experiment: 3 placement arms x 3 hazard regimes (~2-5 min)
python reach_sim.py      # connectivity layer: survival vs. k-anycast reachability, 4 arms (~2-5 min)
python zoo_sim.py        # real-topology validation on Kdl (Topology Zoo) (~2-5 min)
python gdda_sim.py       # adaptive-control material (see caveats below)
```

All outputs are written to `results/` as CSV tables and native SVG figures.

## Repository layout

```
shard_sim_v2.py   main experiment (three arms: random / dispersion / clustered)
reach_sim.py      connectivity-layer experiment, four arms (A_R metric)
zoo_sim.py        real-ISP-topology validation (Topology Zoo, Kdl graph)
topozoo/          Kdl.graphml — Topology Zoo dataset file (Knight et al., JSAC 2011)
gdda_sim.py       legacy adaptive-control prototype (see caveats)
spread_sim.py     superseded v1 prototype (archived)
results/          all raw data (CSV), summaries (TXT) and figures (SVG)
PAPER_MAPPING.md  output-to-paper figure/table mapping (Chinese)
```

## Outputs and paper mapping

| Script | Output | Paper location |
|---|---|---|
| `shard_sim_v2.py` | `results/fig_regime_L15/25/35.svg` | Fig. 4 — A(n) curves in three hazard regimes |
| `shard_sim_v2.py` | `results/fig_theory.svg` | Fig. 5 — simulated dispersion curve vs. mixture-limit bound P(θ ≥ c) |
| `shard_sim_v2.py` | `results/fig_flip.svg` | Fig. 6 — geometry reversal at n=48 |
| `shard_sim_v2.py` | `results/shard_v2_results.csv`, `shard_v2_summary.txt` | Table 3 data source |
| `reach_sim.py` | `results/fig_reach.svg` | Fig. 7 — survival vs. reachability gap, four arms |
| `reach_sim.py` | `results/reach_results.csv`, `reach_summary.txt` | Table 4 data source |
| `zoo_sim.py` | `results/fig_zoo.svg` | Fig. 8 — reachability on the real Kdl topology |
| `zoo_sim.py` | `results/zoo_results.csv`, `zoo_summary.txt` | §6.4 real-topology table data source |
| `gdda_sim.py` | `results/fig_locality.svg`, `fig_sla.svg`, `stability.csv` | §4.4 adaptive control (draft material) |

## Headline results (Monte Carlo, fixed seeds, run 2026-11-12)

Mixture-limit bound P(θ ≥ c) with rate c = 1/2: **L=1.5 → 1.000; L=2.5 → 0.747; L=3.5 → 0.343**

| Phenomenon | Data |
|---|---|
| Redundancy backfire (L=3.5) | dispersion: A = 0.70 at n=2 → **0.29 at n=64** |
| Post-peak decay (L=2.5) | peak 0.90 (n=8) → 0.81 (n=64), converging to the 0.747 bound |
| Geometry reversal (n=48) | L=1.5: dispersion wins (1.000 vs 0.994); L=2.5: clustered wins (0.821 vs 0.773); L=3.5: clustered wins (0.337 vs 0.308) |
| Reachability gap (L=2.5, n=4) | geometric survival 0.839 → reachable service only 0.323 |
| Second reversal | under partitioning, random placement achieves the *highest* reachability and dispersion the lowest |
| R2-corrected clustering | best reachability at every tested n under moderate/strong correlation (0.517 vs 0.272 naive, n=16, L=2.5); gap shrinks 0.487 → 0.195 |

**Real-topology validation (zoo_sim.py, Kdl / Kentucky Datalink, 709 nodes,
run 2026-11-17)** — hazard radii calibrated so E[θ_node] = {0.78, 0.64, 0.43}
matches the synthetic regimes; all four findings replicate, some sharper:
reachable service at n=16 (rand/disp/naive-clust/R2-clust) — weak regime
0.860/0.860/0.865/**0.907**; moderate 0.588/0.487/0.645/**0.782**; strong
0.235/0.170/0.360/**0.540** (R2-clustering beats dispersion by 3.2×).
Dispersion's survival–reachability gap reaches 0.51 while the
R2-corrected arm's stays at 0.02–0.05 in every regime.

## Model summary

- Deployment region: 10 × 10 km; hazard: m = 6 destruction disks of radius L per episode (cluster process; L is the spatial correlation length).
- Placement arms: uniform random (`rand`), max–min dispersion (`disp`), clustered with intra-cluster redundancy (`clust`, cluster count bounded by geometric feasibility ⌊Area / πL²⌋).
- Service model: (n, k) MDS fragmentation, k = ⌈n/2⌉; service available if ≥ k fragments survive (A) and, in the connectivity layer, gatherable within one connected component (A_R).
- Theoretical cross-check: mixture-limit theorem lim A(n) = P(θ ≥ c), where θ is the uncovered area fraction of a hazard realization.

## Known honest limitations

1. The number of feasible independent failure units is geometrically bounded
   (C_max = 8/8/5 in our setting); the "clusters ≈ 1/ρ" heuristic is capped by
   ⌊Area / πL²⌋.
2. The naive clustered arm in the reach experiment places < k fragments per
   cluster, violating cardinality rule R2 by design — its reachability
   collapse is the *experimental evidence for R2's necessity*, not a bug.
3. The mixture-limit theorem holds for the random-placement ensemble; fixed
   placements exhibit finite-n fluctuation (e.g., dispersion at L=3.5, n=64:
   0.294 vs. bound 0.343).
4. `gdda_sim.py` uses mean degree as a proxy for algebraic connectivity and
   non-fair baselines; it is included for the adaptive-control discussion
   only and must not be cited as paper data before repair.
5. The real-topology study covers a single ISP footprint (Kdl); the hazard
   radii are continental-scale (hurricane footprints, 230–430 km), matched
   to the synthetic regimes by measured E[θ], not by absolute length.

## Dataset attribution

`topozoo/Kdl.graphml` is from The Internet Topology Zoo (S. Knight, H. X.
Nguyen, N. Falkner, R. Bowden, M. Roughan, "The Internet Topology Zoo,"
*IEEE JSAC*, 29(9):1765–1775, 2011; http://www.topology-zoo.org), obtained
from the long-term mirror github.com/mroughan/InternetTopologyZoo.

## Parameters

All tunable parameters are at the top of each script with inline comments
(`EPISODES`, `R0S`, `K_RATE`, `M_DISKS`, `NS`, seed counts). Changing any
parameter invalidates the corresponding numbers in the manuscript — re-run
and update all quoted figures accordingly.

## Citation

Citation information will be added upon publication.

## License

MIT — see [LICENSE](LICENSE).
