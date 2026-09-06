# cofold-bench report

Targets: 22 heterodimers; predictors: Boltz-2, AF2-Multimer v3 (ColabFold)

## Per-predictor summary (top-ranked model per target)

| predictor | n | median DockQ | mean DockQ | acceptable+ | medium+ | high | incorrect | Spearman(ipTM, DockQ) | Spearman(pDockQ, DockQ) |
|---|---|---|---|---|---|---|---|---|---|
| Boltz-2 | 22 | 0.830 | 0.750 | 21 | 19 | 13 | 1 | 0.70 | 0.68 |
| AF2-Multimer v3 (ColabFold) | 22 | 0.801 | 0.679 | 19 | 17 | 11 | 3 | 0.79 | 0.86 |

## Paired comparison (same targets)

- n = 22; Boltz-2 higher DockQ on 15, AF2-Multimer higher on 7, ties 0.
- median dDockQ (Boltz-2 - AF2-Multimer) = +0.023; mean = +0.070; IQR [-0.018, +0.067].
- Wilcoxon signed-rank p = 0.156.
- Both acceptable (>=0.23): 19; only Boltz-2: 2; only AF2-Multimer: 0; neither: 1.

## Per-target results (top-ranked model)

| PDB | subset | length | DockQ Boltz-2 | DockQ AF2-M | ipTM Boltz-2 | ipTM AF2-M | title |
|---|---|---|---|---|---|---|---|
| 9FWR | main | 421 | 0.97 | 0.90 | 0.97 | 0.93 | Crystal Structure of SARS-CoV-2 NSP10-NSP14 (ExoN) in comple |
| 9FL4 | main | 329 | 0.97 | 0.91 | 0.97 | 0.93 | compound 5b bound KMT9 crystal structure |
| 9HPX | main | 485 | 0.95 | 0.91 | 0.96 | 0.94 | [FeFe]-hydrogenase from D. desulfuricans with synthetic acti |
| 9O24 | main | 458 | 0.94 | 0.92 | 0.97 | 0.95 | Heparanase P6 in complex with fragment J40 |
| 9RS3 | main | 535 | 0.92 | 0.91 | 0.95 | 0.87 | Crystal structure of the human METTL3-METTL14 in complex wit |
| 9D5P | main | 400 | 0.90 | 0.88 | 0.93 | 0.91 | Crystal structure of the ILK/alpha-parvin core complex bound |
| 9ASS | main | 326 | 0.90 | 0.87 | 0.90 | 0.91 | Crystal Structure of Neutrophil Elastase Inhibited by Eap4 f |
| 9IRO | main | 358 | 0.86 | 0.85 | 0.91 | 0.88 | Crystal structure of SeUGI and SAUDG |
| 9J9W | main | 314 | 0.84 | 0.90 | 0.95 | 0.94 | A broad-spectrum anti-fungal effector dictates bacterial-fun |
| 9Q0B | main | 425 | 0.84 | 0.61 | 0.87 | 0.87 | CTX-M-15 WT in complex with BLIP E73W |
| 9E9G | main | 210 | 0.83 | 0.19 | 0.87 | 0.77 | Heligmosomoides polygyrus TGF-beta Mimic 6 Domain 3 (TGM6-D3 |
| 9G0D | main | 330 | 0.83 | 0.85 | 0.89 | 0.88 | Structure of human Mical1 bMERB_V978A domain:Rab10 complex. |
| 8RO8 | main | 537 | 0.81 | 0.81 | 0.94 | 0.92 | Human cohesin SMC1A-HD(longCC-EQ)/RAD21-C complex - Apo clos |
| 9YGS | main | 251 | 0.78 | 0.73 | 0.88 | 0.87 | Crystal structure of GMPPNP bound KRAS-Y71H in complex with  |
| 9I5Y | main | 442 | 0.75 | 0.77 | 0.83 | 0.85 | Crystal structure of ADP-bound BiP ATPase domain in complex  |
| 8CM0 | main | 324 | 0.71 | 0.79 | 0.83 | 0.84 | Rhs2-CT endonuclease toxin in complex with cognate immunity  |
| 9UO0 | antibody | 343 | 0.61 | 0.01 | 0.90 | 0.53 | Crystal structure of nanobody Tnb165 with MERS-CoV RBD |
| 9OHL | main | 334 | 0.57 | 0.60 | 0.96 | 0.94 | TRMT112-METTL5 bound to SAM and FWG-33B |
| 9EZV | antibody | 261 | 0.57 | 0.35 | 0.81 | 0.60 | Structure of the single-domain antibody VHH_h2 in complex wi |
| 9S61 | main | 347 | 0.47 | 0.69 | 0.94 | 0.86 | Crystal structure of Aurora-A bound to DBL8 |
| 9HMX | main | 493 | 0.46 | 0.38 | 0.69 | 0.78 | Structure of SteB-RipA complex from Mycobacterium tuberculos |
| 9SPR | main | 328 | 0.01 | 0.12 | 0.81 | 0.19 | p53 cancer mutant R282W in complex with DARPin C10 |

## Confident failures (ipTM >= 0.7 but DockQ < 0.23)

- 9E9G / AF2-Multimer v3 (ColabFold): ipTM 0.77, pDockQ 0.21, DockQ 0.19
- 9SPR / Boltz-2: ipTM 0.81, pDockQ 0.24, DockQ 0.01

## MSA-depth side comparison (Boltz-2, 7 targets that fit in 8 GB with default MSA)

| PDB | DockQ 512/128 | DockQ 8192/1024 | dDockQ | ipTM 512/128 | ipTM 8192/1024 |
|---|---|---|---|---|---|
| 9E9G | 0.831 | 0.833 | -0.002 | 0.871 | 0.882 |
| 9Q0B | 0.836 | 0.832 | +0.003 | 0.874 | 0.942 |
| 8CM0 | 0.712 | 0.712 | +0.000 | 0.834 | 0.834 |
| 9J9W | 0.843 | 0.843 | +0.000 | 0.953 | 0.953 |
| 9O24 | 0.942 | 0.931 | +0.011 | 0.973 | 0.976 |
| 9FWR | 0.969 | 0.961 | +0.008 | 0.966 | 0.964 |
| 9HMX | 0.458 | 0.436 | +0.023 | 0.689 | 0.740 |

Mean dDockQ (512/128 minus 8192/1024) = +0.006; max |dDockQ| = 0.023.

## Wall time per target from `boltz_uniform.log` (RTX 4060 Laptop 8 GB; includes MSA server wait)

- n = 21 completed targets timed; median 2.0 min, min 1.1, max 4.3 (last target of a batch is not timed).

## Wall time per target from `colabfold_all_attempt2.log` (RTX 4060 Laptop 8 GB; includes MSA server wait)

- n = 15 completed targets timed; median 3.3 min, min 2.5, max 6.6 (last target of a batch is not timed).

## Wall time per target from `colabfold_all.log` (RTX 4060 Laptop 8 GB; includes MSA server wait)

- n = 5 completed targets timed; median 6.7 min, min 5.5, max 9.6 (last target of a batch is not timed).
