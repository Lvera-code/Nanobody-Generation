# Analysis of the RBD campaign (run1)

Products of the analysis of the 9,988 designs of the RBD campaign. They are kept
here because they are not reproducible without re-running the complete analysis;
the heavy inputs are reproducible and are not versioned.

## Contents

- `RESULTADOS.txt`: the complete numerical output of the analysis.
- `results_overview.pdf`: visual summary.
- `contactos_final_100_designs.csv` and `..._before_refolding.csv`: heavy-atom
  contacts closer than 5 Å recomputed over the top-100 `.cif` files,
  before and after refolding.
- `contactos_muestra_poblacion_500.csv`: the same over a random sample
  of 500 designs from the population (seed 0).
- `correlaciones_spearman.csv`: matrix of correlations between metrics.
- `02_log_clean.txt`, `01_step_configs.txt`, `00_config_dump.txt`: the
  resolved configuration of the 5 steps and the cleaned log of the run.
- `figuras/`: 8 figures of the analysis.
- `scripts/`: the scripts that generate all of the above.

## Non-versioned inputs

The metrics CSVs (`all_designs_metrics.csv`, `aggregate_metrics_analyze.csv`,
`final_designs_metrics_100.csv`, `per_target_metrics_analyze.csv`) and the
top-100 structures live on the compute machine, under
`rbd_campaign/run1/`. They add up to about 90 MB and are downloaded again when
needed.

## Warnings when reusing these data

- `bindsite_under_Xrmsd` is not an RMSD: it is the fraction of epitope residues
  whose nearest designed residue falls below X Å centre to centre,
  measured over the generated structure. The canonical threshold is 8 Å.
- `native_rmsd` and `native_rmsd_bb` are 0.0 in all 9,988 rows because the
  campaign ran with `native: false`. They are not usable.
- `per_target_metrics_analyze.csv` has a single row: the campaign has one
  target and four scaffolds, not four targets.
