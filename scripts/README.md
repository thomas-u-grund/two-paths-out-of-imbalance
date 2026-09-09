# Analysis pipeline

25 R scripts, run in numeric order. Each writes its outputs to `../data/`
(intermediate panels) or `../results/` (models, figures, tables) and, for
most scripts, a matching log to `../logs/`. Unlike the working project's
copies, these scripts use paths relative to the replication pack root
rather than `setwd()`, so run them with the pack root (`replication_pack/`)
as the working directory: `cd replication_pack && Rscript scripts/NN_name.R`.

Core pipeline (01–03) builds the data; 04–25 are estimation, figures, and
robustness/diagnostic checks.

| # | Script | What it does | Paper section |
|---|--------|---------------|----------------|
| 01 | `01_build_dyad_year.R` | Builds the signed dyad-year panel (positive = alliance, negative = MID) from raw Correlates of War data | §3.1 |
| 02 | `02_triads_and_load.R` | Enumerates closed triads per year, classifies balance, computes embeddedness/betweenness/actor constraint | §3.1–3.2 |
| 03 | `03_resolution_panel.R` | Builds the fixed 5-year-window resolution panel (superseded as primary design by 07, kept as the appendix's fixed-window benchmark) | Appendix A5 |
| 04 | `04_analysis.R` | Fixed-window logistic models on the panel from 03 | Appendix A5 |
| 05 | `05_abm_robustness.R` | Compact agent-based simulation: does load-constrained tie dynamics generate persistent imbalance? | §5 discussion |
| 06 | `06_figures.R` | Early descriptive figure (balanced vs. unbalanced triads over time) | Appendix A7 (Figure A1) |
| 07 | `07_competing_risks.R` | **Primary specification.** Discrete-time cause-specific competing-risks hazard models (dissolution vs. realignment) | §3.3, §4.1, Table 2 |
| 08 | `08_ir_controls.R` | Adds contiguity and Polity5 regime-type controls to the primary hazards | §4.5 |
| 09 | `09_h2_tie_choice.R` | H2, first pass: identifies which specific tie changed sign at each realignment event | §4.3 |
| 10 | `10_dissolution_disaggregation.R` | Disaggregates dissolution events by tie type (alliance vs. MID) and state-exit coincidence | §4.4 |
| 11 | `11_abm_threshold_robustness.R` | ABM robustness across different dissolution-eligibility thresholds | §5 discussion |
| 12 | `12_h2_unbalanced_load.R` | H2, revised: replaces total embeddedness with *unbalanced load* (only unbalanced overlapping triangles) as the tie-choice predictor | §4.3, Table 3 |
| 13 | `13_h2_heterogeneity.R` | Splits H2 by direction of change (alliance→hostility vs. hostility→alliance) | §4.3, Table 4–5 |
| 14 | `14_duration_dependence.R` | Adds spell-duration fixed effects to the primary hazards (the actual Table 2 specification) and re-runs trade/contiguity/regime robustness with it | Table 2, Appendix A1 |
| 15 | `15_clique_diagnostic.R` | Diagnoses the pre-1945 volatility in Figure A1: traced to a handful of multilateral alliances dominating the closed-triad population | Appendix A7 |
| 16 | `16_bloc_robustness.R` | Extends 15 to the full set of post-1945 blocs (NATO, Warsaw Pact, Rio Pact, etc.); two-way triad/bloc clustering check | §4.5, Appendix A2 |
| 17 | `17_twomode_proxy.R` | Cheap two-mode (state × organization) proxy check for the one-mode alliance-clique projection | Appendix A3 |
| 18 | `18_regenerate_fig2.R` | Generates the hazard-coefficient plot (Figure 3) | Figure 3 |
| 19 | `19_schematic_figure.R` | Generates the schematic figure explaining embeddedness vs. unbalanced load (Figure 1) | Figure 1 |
| 20 | `20_endogeneity_checks.R` | Lagged-covariate and omitted-variable-bias sensitivity checks | §3.4, Appendix A4 |
| 21 | `21_fig_cif_greyscale.R` | Generates the cumulative-incidence plot (Figure 2) | Figure 2 |
| 22 | `22_h2_sign_and_dependence.R` | H2 diagnostics: sign-composition test, dyad-year dependence check, same-sign comparison | §4.3, Appendix A12 |
| 23 | `23_h2_dependence_and_multinomial.R` | H2 dependence-correction follow-up, plus a multinomial discrete-time robustness model | Appendix A12 |
| 24 | `24_h1_unbalanced_load_control.R` | Tests whether H1's embeddedness–realignment effect survives controlling for the triad's current unbalanced load (exposure vs. pressure) | §4.1, Appendix A13 |
| 25 | `25_weak_balance_robustness.R` | Reclassifies imbalance under Davis's (1967) weak-balance criterion (all-negative triads treated as balanced), rebuilds the spell structure, and re-estimates the primary hazard models | §3.1, Appendix A14 |
| 26 | `26_cloglog_robustness.R` | Re-estimates the primary hazard specification (Table 2) with a complementary log-log link instead of logit, to check whether the link choice drives the results | §3.3, Appendix A15 |

## Legacy scripts

The project's original prototype pipeline (`balance_A_generate.R` through
`balance_D_analysis.R`, plus `balance_main.R`) predates this rebuild and
is **not** part of the current analysis; none of the paper's reported
results depend on it (the rebuild re-downloaded all data fresh rather
than reusing its outputs). It is not included in this replication pack.
