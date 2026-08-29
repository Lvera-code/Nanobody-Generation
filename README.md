# Proyecto3-BoltzGen

Tercer proyecto del CNB (encargado por Carlos Oscar Sorzano): diseno de novo
de VHH (nanobodies) dirigido por epitopo, envolviendo
[BoltzGen](https://github.com/HannesStark/boltzgen) (`HannesStark/boltzgen`,
MIT) como plugin de Scipion. Independiente de
[BCell-Epitope-Prediction](https://github.com/Lvera-code/BCell-Epitope-Prediction)
(proyecto 1) y [PTM-Prediction](https://github.com/Lvera-code/PTM-Prediction)
(proyecto 2), pero reutiliza codigo de ambos por import directo, nunca
copiado.

## Fase 0: seleccion de candidatos B-cell para BoltzGen

`fase0_seleccion_bcell.py` es un fork corregido/enriquecido de
`_select_bcell_candidates` (Fase 7 de BCell-Epitope-Prediction), sin
modificar ese repo (sus resultados alimentan un paper en escritura). Reusa
por import todo lo que es seguro reusar tal cual, y corrige/anade
localmente lo que necesitaba un cambio real:

1. Recorte de candidatos largos con re-chequeo de alergenicidad/N-
   glicosilacion sobre la secuencia YA recortada (el original solo lo
   re-chequeaba en el padding, no en el recorte).
2. `documented_region`/columnas `bnab_*`/`iedb_*` recalculadas sobre la
   secuencia FINAL (post pad/trim/autotolerancia), no sobre la region padre.
3. Filtro DURO por TMbed (helice/tira TM, peptido senal, intracelular): una
   region asi no es una diana valida para un nanobody contra el antigeno
   nativo plegado, a diferencia de un constructo lineal (donde solo se
   anota).
4. Anotacion de accesibilidad real de superficie (RSA/pLDDT de
   DiscoTope-3.0), retenida desde el residuo hasta la ventana final.
5. Anotacion de evidencia de via secretora (UniProt REST, reusando el
   cliente de PTM-Prediction) -- el antigeno completo debe ser plausible
   como diana de superficie, no solo el residuo.
6. Anotacion de estructura cuaternaria via el comentario `SUBUNIT` de
   UniProt (texto libre, sin parsear a un enum).
7. Anotacion de copias reales en el ensamblaje biologico via la API de
   estructuras de RCSB (generico para cualquier PDB ID, sin curacion manual
   por target -- necesario porque este pipeline se reusara tal cual para
   candidatos de competicion como Adaptyv, no solo el panel GP120/MonkeyPox).

Validado de principio a fin contra un caso real (HIV1 gp120, PDB `3NGB`,
`top_n=6`) incluyendo las 7 diferencias.

## Proximos pasos

- Scaffolding del plugin `scipion-chem-boltzgen`.
- Wrapper del protocolo Scipion para la Fase 1 de BoltzGen (generacion),
  alimentado por la salida de Fase 0.
