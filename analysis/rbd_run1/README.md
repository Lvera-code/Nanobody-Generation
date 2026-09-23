# Análisis de la campaña RBD (run1)

Productos del análisis de los 9.988 diseños de la campaña de RBD. Se guardan
aquí porque no son reproducibles sin volver a ejecutar el análisis completo;
las entradas pesadas sí lo son y no se versionan.

## Contenido

- `RESULTADOS.txt`: salida numérica completa del análisis.
- `results_overview.pdf`: resumen visual.
- `contactos_final_100_designs.csv` y `..._before_refolding.csv`: contactos de
  átomos pesados a menos de 5 Å recalculados sobre los `.cif` del top-100,
  antes y después del replegado.
- `contactos_muestra_poblacion_500.csv`: lo mismo sobre una muestra aleatoria
  de 500 diseños de la población (semilla 0).
- `correlaciones_spearman.csv`: matriz de correlaciones entre métricas.
- `02_log_clean.txt`, `01_step_configs.txt`, `00_config_dump.txt`:
  configuración resuelta de los 5 pasos y log depurado de la corrida.
- `figuras/`: 8 figuras del análisis.
- `scripts/`: los scripts que generan todo lo anterior.

## Entradas no versionadas

Los CSV de métricas (`all_designs_metrics.csv`, `aggregate_metrics_analyze.csv`,
`final_designs_metrics_100.csv`, `per_target_metrics_analyze.csv`) y las
estructuras del top-100 viven en la máquina de cómputo, bajo
`rbd_campaign/run1/`. Suman unos 90 MB y se vuelven a descargar cuando hagan
falta.

## Avisos al reutilizar estos datos

- `bindsite_under_Xrmsd` no es un RMSD: es la fracción de residuos del epítopo
  cuyo residuo diseñado más cercano queda por debajo de X Å centro a centro,
  medida sobre la estructura generada. El umbral canónico es 8 Å.
- `native_rmsd` y `native_rmsd_bb` valen 0,0 en las 9.988 filas porque la
  campaña corrió con `native: false`. No son utilizables.
- `per_target_metrics_analyze.csv` tiene una sola fila: la campaña tiene una
  diana y cuatro andamios, no cuatro dianas.
