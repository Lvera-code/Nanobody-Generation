# Proyecto3-BoltzGen

Third CNB project (commissioned by Carlos Oscar Sorzano): epitope-directed
de novo design of VHH (nanobodies), wrapping
[BoltzGen](https://github.com/HannesStark/boltzgen) (`HannesStark/boltzgen`,
MIT) as a Scipion plugin. Independent of
[BCell-Epitope-Prediction](https://github.com/Lvera-code/BCell-Epitope-Prediction)
(project 1) and [PTM-Prediction](https://github.com/Lvera-code/PTM-Prediction)
(project 2), but it reuses code from both by direct import, never
copied.

## Phase 0: selection of B-cell candidates for BoltzGen

`fase0_seleccion_bcell.py` is a corrected/enriched fork of
`_select_bcell_candidates` (Phase 7 of BCell-Epitope-Prediction), without
modifying that repo (its results feed a paper in writing). It reuses
by import everything that is safe to reuse as-is, and corrects/adds
locally what needed a real change:

1. Trimming of long candidates with a re-check of allergenicity/N-
   glycosylation over the ALREADY trimmed sequence (the original only
   re-checked it in the padding, not in the trim).
2. `documented_region`/the `bnab_*`/`iedb_*` columns recomputed over the
   FINAL sequence (post pad/trim/self-tolerance), not over the parent region.
3. A HARD filter by TMbed (TM helix/strand, signal peptide, intracellular): a
   region like that is not a valid target for a nanobody against the native
   folded antigen, unlike a linear construct (where it is only
   annotated).
4. Annotation of real surface accessibility (RSA/pLDDT from
   DiscoTope-3.0), retained from the residue through to the final window.
5. Annotation of secretory pathway evidence (UniProt REST, reusing
   PTM-Prediction's client). The complete antigen must be plausible
   as a surface target, not just the residue.
6. Annotation of quaternary structure via UniProt's `SUBUNIT`
   comment (free text, not parsed into an enum).
7. Annotation of the real number of copies in the biological assembly via the RCSB
   structures API (generic for any PDB ID, without manual curation
   per target, which is necessary because this pipeline will be reused as-is for
   competition candidates such as Adaptyv, not only the GP120/MonkeyPox panel).

Validated end to end against a real case (HIV1 gp120, PDB `3NGB`,
`top_n=6`) including all 7 differences.

## Next steps

- Scaffolding of the `scipion-chem-boltzgen` plugin.
- Scipion protocol wrapper for Phase 1 of BoltzGen (generation),
  fed by the output of Phase 0.
