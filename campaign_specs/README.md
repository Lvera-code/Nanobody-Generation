# Especificaciones de campaña

Entradas exactas con las que se lanzaron las campañas de diseño en la máquina
de cómputo, guardadas aquí para que una campaña pueda repetirse sin depender
del sistema de ficheros de esa máquina.

## rbd/

Campaña de nanocuerpos contra el RBD de SARS-CoV-2.

- `7BEP_chainE_boltzgen_target.cif`: la diana tal y como la ingiere BoltzGen.
  Cadena `Exp`, 197 residuos, numeración 1..197.
- `design_spec.yaml`: especificación de diseño. La ventana de epítopo
  `binding: 148..163` se refiere a esa numeración 1..197 y corresponde a los
  residuos 479..494 de la proteína Spike (`PCNGVEGFNCYFPLQS`), el núcleo del
  RBM que incluye E484, F486, Y489 y Q493.

Comando de lanzamiento:

```
boltzgen run --protocol nanobody-anything --output run1 --devices 4 \
             --num_designs 10000 --diffusion_batch_size 1 --reuse \
             --budget 100 design_spec.yaml
```

### Relación con `validation_panel/reference_targets/`

`validation_panel/reference_targets/7BEP_chainE_boltzgen_target.cif` contiene
**la misma estructura** (197 residuos, 1556 átomos, desviación máxima de
coordenadas 0.0000 Å) pero con la cadena nombrada `E` y la numeración
original de Spike (332..528). Esa variante sirve para comparar contra
estructuras del PDB; la de `campaign_specs/rbd/` es la que consume
`design_spec.yaml`. No son intercambiables: con la variante del panel, el
`chain: id: Exp` no resuelve y la ventana `148..163` seleccionaría residuos
distintos.

### Procedencia de los andamios

Las cuatro especificaciones de andamio de nanocuerpo (`7eow.yaml`,
`7xl0.yaml`, `8coh.yaml`, `8z8v.yaml`) que referencia `design_spec.yaml` **no
son artefactos de este proyecto**: son ficheros de ejemplo del propio BoltzGen,
sin modificar. Viven en `boltzgen/example/nanobody_scaffolds/`, que este repo
ignora a propósito (ver `.gitignore`) por ser un clon del repositorio ajeno.

Se reproducen clonando el fork en el commit `7ae564b`. Verificado: las copias
locales son idénticas byte a byte a las que consumió la campaña en la máquina
de cómputo, y `git status` sobre ese directorio está limpio respecto al fork.
