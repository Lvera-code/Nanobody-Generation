"""Fase 0 de Proyecto 3 (BoltzGen): seleccion de candidatos B-cell, version corregida.

**Por que este archivo existe por separado y no modifica Project 1**: Project 1
(``B-Cell-Epitope-Prediction``) es feature-complete y sus resultados ya
alimentan un paper en escritura -- cualquier cambio de comportamiento en su
``construct_assembly.py`` (no solo columnas nuevas sin consumidores, sino
valores que YA aparecen en ``source_score_note``/metadata persistida) podria
alterar numeros ya citados si el pipeline se vuelve a correr. Este modulo
COPIA la logica de seleccion B-cell de ``_select_bcell_candidates`` (Fase 7
de Project 1) y la corrige/enriquece aqui, sin tocar el original. Todo lo que
NO se modifica se REUTILIZA por import directo de Project 1 (no se duplica)
para no divergir en silencio si Project 1 corrige algo en esas partes.

Siete diferencias reales respecto al original (ver commit/sesion 2026-08-29):

1. **Bug de staleness en `allergen`/`glycosylated` para candidatos RECORTADOS**:
   el `_trim_long_bcell_candidates` original de Project 1 recorta la
   secuencia pero nunca re-corre AlgPred2/StackGlyEmbed sobre la ventana
   final -- dejaba el veredicto de la region PADRE mas ancha. La version
   `_trim_long_bcell_candidates_fresco` de aqui SI re-chequea, exactamente
   igual que el padding ya hace en el original (mismo patron, solo se
   extiende al caso que faltaba).
2. **`documented_region`/nuevas columnas `bnab_*`/`iedb_*` calculadas sobre
   la secuencia FINAL** (post pad/trim/self-tolerancia), no sobre la region
   padre -- reusa los mismos matches de referencia ya encontrados por Fase
   6/6c de Project 1 (``bnab_df``/``iedb_df``), re-verificando el
   solapamiento real contra la ventana final con la misma funcion pura
   (``_longest_common_substring_len``) que el motor original usa. Sin
   lectura de disco nueva: el universo de epitopos de referencia posible es
   el que ya vino en esos DataFrames.
3. **Filtro duro por `tmbed_masked` (TM helix/strand, signal peptide,
   intracellular)**: Project 1 solo anota esto (ADR 2026-08-13, un
   constructo lineal no necesita descartar la region). Proyecto 3 SI la
   descarta -- una region asi no es una diana de union fisicamente valida
   para un nanobody contra el antigeno nativo plegado. Se relee
   `{input_stem}_tmbed_regions.csv` fresco (ver `_load_tmbed_regions`) y se
   recalcula el solapamiento contra la ventana FINAL (post pad/trim), en vez
   de confiar en la columna `tmbed_masked` heredada de `safe_df`, que queda
   obsoleta en cuanto Fase 7 recorta/extiende `start`/`end`.
4. **Anotacion (no filtro) de accesibilidad real de superficie (RSA)**:
   Project 1 retiene `rsa`/`pLDDT`/`alphafold_struc_flag` de DiscoTope-3.0
   desde el 2026-08-29 (ADR en `discotope_engine.py`) pero se pierden al
   colapsar a region en Fase 3 -- nunca llegan a `safe_df`. Aqui se releen
   frescos del raw CSV cacheado y se exponen como `mean_rsa`/`mean_plddt`/
   `alphafold_struc_flag` por candidata (ventana FINAL). NO es un filtro
   duro como `tmbed_masked`: RSA es continuo y no hay umbral
   "enterrado/expuesto" publicado por los autores de DiscoTope-3.0 (a
   diferencia de `calibrated_score`), asi que fijar un corte aqui seria
   fabricar una constante biologica sin fuente -- queda como criterio de
   desempate explicito para quien arme el YAML final de BoltzGen.

5. **Anotacion (no filtro) de evidencia de VIA SECRETORA** (reuso cruzado
   con PTM-Prediction, ver sesion 2026-08-29 "que falta enriquecer"):
   `mean_rsa`/`alphafold_struc_flag` (diferencia 4) dicen si un residuo esta
   expuesto DENTRO del plegamiento que se modelo, pero no dicen si el
   antigeno completo es siquiera una diana de superficie plausible -- un
   residuo con RSA alto en una proteina en realidad citoplasmatica/nuclear
   no es una diana valida para un nanobody. `PTM-Prediction/src/structural/
   uniprot_localization_client.py` ya resuelve exactamente esto (evidencia
   real de UniProt, keywords tipo `Secreted`/`Signal`/`Cell membrane` o el
   comentario `SUBCELLULAR LOCATION`), construido alli para N-glicosilacion
   pero la senal es generica. Se carga por ruta de archivo (``importlib``,
   ver `_load_uniprot_localization_client`) en vez de `import src...`
   porque PTM-Prediction y Project 1 comparten nombre de paquete top-level
   (`src`) -- un import normal del segundo repo pisaria el modulo `src` ya
   registrado en `sys.modules` desde el primero. El `accession` de Project 1
   es el stem completo del fichero de entrada (`HIV1_gp41_full_P04578_512-
   856`), no un ID UniProt limpio -- se extrae con un regex (patron oficial
   de UniProt, `_UNIPROT_ACCESSION_TOKEN_RE`) sobre los tokens separados por
   `_`/`-`/`.` del stem; si ningun token calza el patron (comun cuando el
   accession no viene de UniProt), la anotacion queda en ``None`` sin
   fallar, mismo comportamiento no-fatal que en PTM-Prediction.
6. **Anotacion (no filtro) de estructura CUATERNARIA/oligomerica** (Nivel A
   de la robustez para targets arbitrarios, sesion 2026-08-29 -- este
   pipeline se reusara tal cual para candidatos de la competicion Adaptyv,
   no solo el panel GP120/MonkeyPox, asi que no puede depender de curacion
   manual por target): el comentario ``SUBUNIT`` real de UniProt (texto
   libre curado, p.ej. "consists of a homotrimer of non-covalently
   associated gp120-gp41 heterodimers" para P04578) se expone tal cual en
   ``quaternary_structure_note`` -- deliberadamente SIN parsearlo a un enum
   monomero/dimero/trimero (fabricaria una regla de lenguaje natural sin
   fuente verificada). Complementa a RSA/via secretora: ninguna de las dos
   dice si la ventana elegida queda enterrada en una interfaz de
   oligomerizacion real. El "Nivel B" (verificacion estructural real via la
   API de ensamblajes de RCSB quando el Camino es PDB) se investigo la
   misma sesion y es tecnicamente viable de forma generica
   (``rcsb_struct_symmetry[].oligomeric_state``/``stoichiometry`` por
   ``assembly_id``, verificado con datos reales: 4NCO/Env-SOSIP-trimer da
   ``Hetero 12-mer``, stoichiometry ``A3/B3/C3/D3``; 1LYZ/lisozima-monomero
   da ``Monomer``, ``A1``) pero exige mapear asym_id/entidad al chain
   elegido como antigeno y solo cubre estructuras PDB depositadas (no
   AlphaFold ni targets sin estructura resuelta) -- no implementado
   todavia, decision pendiente de si se construye antes del scaffolding.
7. **Anotacion (no filtro) de COPIAS REALES en el ensamblaje biologico**
   (Nivel B, implementado misma sesion 2026-08-29 tras confirmar
   viabilidad): ``assembly_copy_number``/``assembly_oligomeric_state``/
   ``assembly_symmetry_kind`` (ver `_annotate_assembly_copy_number`) usan la
   API de estructuras de RCSB para contar, de forma generica (cualquier PDB
   ID, sin curacion manual), cuantas veces aparece la cadena elegida por
   Fase 1.5 en el ensamblaje biologico real (``rcsb_struct_symmetry``).
   Requiere resolver primero auth_asym_id -> label_asym_id: verificado con
   datos reales que NO son arrays paralelos indexables (``7BEP``, entidad
   con 2 cadenas: ``auth_asym_ids=['C','E']``/``asym_ids=['A','F']``, pero
   el auth real de ``E`` es label ``'A'``, no ``'F'`` -- emparejar por
   indice habria dado un resultado silenciosamente incorrecto, ver
   ``_resolve_label_asym_id`` para el metodo correcto verificado, indexado
   por instancia). Con la resolucion correcta, ``7BEP`` (RBD de SARS-CoV-2,
   panel de Project 1) da ``copy_number=1`` -- coherente con que la RBD no
   es un oligomero obligado por si misma. CAVEAT que sigue siendo real de
   forma conceptual (no ligado a este caso concreto, ya corregido): el
   "ensamblaje biologico" que RCSB expone como preferido puede en otras
   entradas seguir incluyendo copias cristalograficas de un complejo
   completo en vez de auto-oligomerizacion real -- por eso esta anotacion es
   PURAMENTE INFORMATIVA y debe leerse junto al texto de UniProt del punto
   6, nunca de forma aislada. Solo aplica a Camino PDB con ID real
   extraible del stem (``_extract_pdb_id``); ``None`` para Camino FASTA,
   AlphaFold, o targets de competicion sin estructura resuelta.

`conservation_pct` NO se corrige aqui (se deja igual que en Project 1,
pre-recorte): arreglarlo de verdad exigiria correr BLASTp de nuevo contra el
panel de conservacion, que Project 1 corre deliberadamente ANTES de su Fase
7 para mantenerla "logica pura" -- ver ADR en `fase_6b_conservacion` de
`pipeline.py`. Puede abordarse en una sesion futura si hace falta.

Uso previsto: llamar a ``select_bcell_candidates_for_boltzgen`` con los
mismos DataFrames que ya produce el pipeline de Project 1 hasta Fase 6c
(``safe_df``, ``algpred_df``, ``stackgly_df``, y opcionalmente
``bnab_df``/``iedb_df`` de sus Fases 6/6c) -- el resultado es el candidato
final listo para construir el ``binding_types`` del YAML de BoltzGen.
"""

import importlib.util
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

_B_CELL_REPO = Path(__file__).resolve().parent.parent / "B-Cell-Epitope-Prediction"
if str(_B_CELL_REPO) not in sys.path:
    sys.path.insert(0, str(_B_CELL_REPO))

_PTM_PREDICTION_REPO = Path(__file__).resolve().parent.parent / "PTM-Prediction"
_UNIPROT_CLIENT_PATH = _PTM_PREDICTION_REPO / "src" / "structural" / "uniprot_localization_client.py"

import pandas as pd

from src.config.settings import Settings
from src.engines import discotope_engine
from src.engines.algpred_engine import predict_allergenicity
from src.engines.blast_engine import filter_self_tolerant
from src.engines.construct_assembly import (
    _annotate_conformational_coverage,
    _best_subwindow,
    _load_native_residues,
    _load_per_residue_scores,
    _pad_short_bcell_candidates,
)
from src.engines.iedb_engine import query_iedb_crossref
from src.engines.lanl_catnap_engine import _longest_common_substring_len, query_bnab_crossref
from src.engines.stackglyembed_engine import predict_nglycosylation
from src.utils.logger_config import setup_logger


def _load_uniprot_localization_client():
    """Carga ``uniprot_localization_client.py`` de PTM-Prediction por ruta de archivo.

    No usa ``import src.structural.uniprot_localization_client`` porque
    PTM-Prediction y Project 1 comparten el nombre de paquete top-level
    (``src``) -- Project 1 ya registro ``src`` en ``sys.modules`` (import de
    arriba), y un import normal del segundo repo buscaria ``structural``
    dentro del ``src`` de Project 1 (que no lo tiene) en vez de cargar el de
    PTM-Prediction. El modulo es autocontenido (solo ``json``/``urllib``/
    ``typing`` de la libreria estandar, ver su propio docstring), asi que
    cargarlo por ruta con un nombre unico es seguro -- no arrastra ningun
    otro modulo de PTM-Prediction ni su paquete ``src``.
    """
    spec = importlib.util.spec_from_file_location(
        "ptm_prediction_uniprot_localization_client", _UNIPROT_CLIENT_PATH
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_uniprot_client = _load_uniprot_localization_client()
UniProtLookupError = _uniprot_client.UniProtLookupError
lookup_secretory_pathway_evidence = _uniprot_client.lookup_secretory_pathway_evidence

logger = setup_logger(__name__)


def _trim_long_bcell_candidates_fresco(
    candidates: pd.DataFrame, output_dir: Path, input_stem: str, max_length: int,
) -> pd.DataFrame:
    """Copia corregida de ``construct_assembly._trim_long_bcell_candidates``.

    Identica al original de Project 1 en la eleccion de la mejor sub-ventana
    (``_best_subwindow``, reusada sin cambios), PERO ademas re-corre
    AlgPred2/StackGlyEmbed sobre la secuencia ya recortada -- el original
    dejaba ``allergen``/``glycosylated`` calculados sobre la region padre sin
    refrescar. Mismo patron de 2 pasadas (recortar todo primero, correr los
    motores una sola vez en lote, despues escribir resultados) que
    ``_pad_short_bcell_candidates`` ya usa en el original, para no reinventar
    un patron nuevo.
    """
    if candidates.empty:
        return candidates
    trimmed = candidates.copy()
    original_lengths = trimmed["sequence"].str.len()

    proposals: Dict[int, Tuple[str, int, int]] = {}
    for idx, row in trimmed.iterrows():
        if len(row["sequence"]) <= max_length:
            continue
        per_residue = _load_per_residue_scores(output_dir, input_stem, row["accession"])
        new_seq, new_start, new_end = _best_subwindow(row["sequence"], int(row["start"]), per_residue, max_length)
        proposals[idx] = (new_seq, new_start, new_end)

    trimmed["trimmed_from_length"] = original_lengths.where(original_lengths > max_length)
    if not proposals:
        return trimmed

    proposal_seqs = sorted({seq for seq, _, _ in proposals.values()})
    algpred_check = predict_allergenicity(proposal_seqs, output_dir, filename_prefix=f"{input_stem}_bcell_trim_")
    glyco_check = predict_nglycosylation(proposal_seqs, output_dir, filename_prefix=f"{input_stem}_bcell_trim_")
    allergen_seqs = set(algpred_check[algpred_check["algpred_veredicto"] == "Allergen"]["sequence"]) \
        if not algpred_check.empty else set()
    glyco_risky_seqs = set(glyco_check[glyco_check["stackglyembed_veredicto"] == "Glicosilado"]["sequence"]) \
        if not glyco_check.empty else set()

    for idx, (new_seq, new_start, new_end) in proposals.items():
        trimmed.at[idx, "sequence"] = new_seq
        trimmed.at[idx, "start"] = new_start
        trimmed.at[idx, "end"] = new_end
        trimmed.at[idx, "allergen"] = new_seq in allergen_seqs
        trimmed.at[idx, "glycosylated"] = new_seq in glyco_risky_seqs

    return trimmed


# Orden ordinal REAL de IEDB (no inventado): verificado contra los 4 valores
# realmente presentes en reference_db/iedb/bcell_protective_epitopes.csv
# (2026-08-29): 'Positive-High'=155, 'Positive'=2965, 'Positive-Low'=248,
# 'Positive-Intermediate'=61 filas. El orden High>Intermediate>Low>Positive
# es la lectura literal de los propios calificadores de IEDB, no un umbral
# biologico fabricado.
_QUALITATIVE_MEASURE_RANK = {
    "Positive-High": 3,
    "Positive-Intermediate": 2,
    "Positive-Low": 1,
    "Positive": 0,
}


def _strongest_qualitative_measure(values) -> str:
    """El ``qualitative_measure`` de IEDB con mas fuerza de evidencia entre los presentes.

    Valores no reconocidos (no deberian aparecer, ver panel real arriba) caen
    al fondo del orden (rank -1) en vez de lanzar una excepcion.
    """
    return max(values, key=lambda v: _QUALITATIVE_MEASURE_RANK.get(v, -1))


def _matches_within_overlap(sequence: str, reference_df: pd.DataFrame, min_overlap: int) -> pd.DataFrame:
    """Sub-tabla de ``reference_df`` (columna ``epitope_sequence``) que SI solapa con ``sequence``.

    Misma regla exacta que ``query_bnab_crossref``/``query_iedb_crossref``
    (``required_overlap = min(min_overlap, len(epitope_sequence))``,
    solapamiento via ``_longest_common_substring_len``) -- reutilizada aqui,
    no reinventada, para no divergir del criterio real del motor original.
    """
    if reference_df is None or reference_df.empty:
        return reference_df.iloc[0:0] if reference_df is not None else pd.DataFrame()
    seq_upper = sequence.upper()
    mask = reference_df["epitope_sequence"].apply(
        lambda ep: _longest_common_substring_len(seq_upper, ep) >= min(min_overlap, len(ep))
    )
    return reference_df[mask]


def _load_tmbed_regions(output_dir: Path, input_stem: str) -> pd.DataFrame:
    """Regiones TM/signal/intracelular reales de TMbed (Fase 3b de Project 1), leidas frescas.

    Lee ``{output_dir}/{input_stem}_tmbed_regions.csv`` (columnas
    ``accession,start,end,type``, 1-indexado inclusivo, misma convencion que
    ``start``/``end`` de las candidatas B-cell) -- el mismo archivo que
    ``fase_3b_tm_signal_masking`` escribe en Project 1
    (``src/engines/tmbed_engine.py``), pero que nada relee despues de esa
    fase. No existe helper reutilizable para esto en Project 1: la columna
    ``tmbed_masked`` que SI se propaga hasta ``safe_df`` queda obsoleta en
    cuanto Fase 7 recorta/extiende ``start``/``end`` (``_trim_long_bcell_candidates``/
    ``_pad_short_bcell_candidates`` no la refrescan) -- por eso aqui se relee
    el CSV crudo y se recalcula el solapamiento contra la ventana FINAL, en
    vez de confiar en la columna heredada.
    """
    path = output_dir / f"{input_stem}_tmbed_regions.csv"
    if not path.exists():
        return pd.DataFrame(columns=["accession", "start", "end", "type"])
    return pd.read_csv(path)


# Los 4 tipos que TMbed puede marcar (ver `_MASKED_CLASS_TYPES` en
# `tmbed_engine.py`) -- a diferencia de Project 1 (donde son solo
# informativos para un constructo lineal, ADR 2026-08-13), ninguno es una
# diana de union valida para un nanobody contra el antigeno nativo plegado:
# TM_beta_strand/TM_alpha_helix estan embebidos en la membrana,
# signal_peptide se escinde co-traduccionalmente, intracellular no es
# accesible en la superficie.
_TMBED_HARD_FILTER_TYPES = {"TM_beta_strand", "TM_alpha_helix", "signal_peptide", "intracellular"}


def _filter_tmbed_masked(candidates: pd.DataFrame, output_dir: Path, input_stem: str) -> pd.DataFrame:
    """Filtro DURO especifico de Proyecto 3: descarta candidatas cuya ventana FINAL
    solapa una region TM/signal/intracelular real de TMbed.

    A diferencia de Project 1 (donde ``tmbed_masked`` es solo anotacion), una
    region asi no es una diana fisicamente valida para un nanobody dirigido
    contra el antigeno nativo plegado -- de ahi que aqui SI se aplique como
    filtro. Recalcula el solapamiento contra ``start``/``end`` ya finales
    (post pad/trim), con la misma logica de interval-overlap que
    ``annotate_overlapping_regions`` usa en Project 1 (``tmbed_engine.py``),
    reusada aqui en vez de reinventada.
    """
    if candidates.empty:
        return candidates
    regions = _load_tmbed_regions(output_dir, input_stem)
    regions = regions[regions["type"].isin(_TMBED_HARD_FILTER_TYPES)]
    if regions.empty:
        return candidates

    def _is_masked(row) -> bool:
        acc_regions = regions[regions["accession"] == row["accession"]]
        overlaps = acc_regions[(acc_regions["start"] <= row["end"]) & (acc_regions["end"] >= row["start"])]
        return not overlaps.empty

    masked = candidates.apply(_is_masked, axis=1)
    return candidates[~masked]


def _load_rsa_scores(output_dir: Path, input_stem: str, accession: str) -> pd.DataFrame:
    """RSA/pLDDT/alphafold_struc_flag por residuo de DiscoTope-3.0, releidos frescos.

    Mismo patron que ``_load_per_residue_scores`` de Project 1: relee
    ``{output_dir}/{input_stem}_discotope_raw.csv`` (cacheado por Fase 2,
    ``_cached_structural_raw_scores`` en ``pipeline.py``), filtra por
    ``accession`` y usa el ORDEN de fila como posicion 1-indexada (el CSV no
    trae columna de posicion explicita). Estas 3 columnas las retiene
    Project 1 desde el ADR 2026-08-29 en ``discotope_engine.py``,
    especificamente por y para Proyecto 3: ``rsa`` distingue "enterrado en
    el propio plegamiento del dominio soluble" (fenomeno estructural) de
    ``tmbed_masked`` ("insertado en membrana/intracelular", contexto de
    dominio) -- un motor de antigenicidad de secuencia puede puntuar alto
    una ventana que en la estructura 3D real no esta expuesta.

    Returns:
        DataFrame con columnas ``position``/``rsa``/``plddt``/``alphafold_struc_flag``,
        vacio si DiscoTope-3.0 no corrio para este ``input_stem`` (input sin
        estructura, Camino FASTA) -- mismo comportamiento no-fatal que
        ``_load_tmbed_regions``.
    """
    empty = pd.DataFrame(columns=["position", "rsa", "plddt", "alphafold_struc_flag"])
    raw_path = output_dir / f"{input_stem}_discotope_raw.csv"
    if not raw_path.is_file():
        return empty
    raw_df = pd.read_csv(raw_path)
    required = {discotope_engine.ACCESSION_COLUMN, discotope_engine.RSA_COLUMN}
    if not required.issubset(raw_df.columns):
        return empty
    acc_rows = raw_df[raw_df[discotope_engine.ACCESSION_COLUMN] == accession].reset_index(drop=True)
    if acc_rows.empty:
        return empty
    return pd.DataFrame({
        "position": acc_rows.index + 1,
        "rsa": acc_rows[discotope_engine.RSA_COLUMN],
        "plddt": acc_rows[discotope_engine.PLDDT_COLUMN],
        "alphafold_struc_flag": acc_rows[discotope_engine.ALPHAFOLD_STRUC_FLAG_COLUMN],
    })


def _annotate_rsa_accessibility(candidates: pd.DataFrame, output_dir: Path, input_stem: str) -> pd.DataFrame:
    """Anota (sin filtrar) la accesibilidad de superficie real (RSA) de la ventana FINAL de cada candidata.

    Puramente informativo, mismo patron no-decisorio que
    ``conservation_pct``/``documented_region``/bnAb (ver decision de diseno
    de Fase 0, 2026-08-29): a diferencia de ``tmbed_masked`` (categorias
    fijas TM/signal/intracelular = imposibilidad fisica binaria, SI se
    filtra duro), RSA es una variable continua y no hay un umbral
    "enterrado vs expuesto" publicado por los autores de DiscoTope-3.0 (a
    diferencia de ``calibrated_score``, para el que si publican umbrales de
    referencia -- ver ADR "calibrated_score" en ``discotope_engine.py``).
    Inventar un corte aqui violaria la misma regla ya aplicada en el resto
    del pipeline: no fabricar constantes biologicas sin fuente verificada.
    Se expone como ``mean_rsa``/``mean_plddt``/``alphafold_struc_flag`` para
    que quien arme el YAML final de BoltzGen lo use como criterio de
    desempate explicito, igual que ``conservation_pct``/``documented_region``.

    Requiere el raw CSV de DiscoTope-3.0 ya cacheado (Camino PDB); si no
    esta disponible (input sin estructura), no anota nada para esa
    candidata en vez de fallar -- mismo comportamiento que
    ``_annotate_conformational_coverage``.
    """
    if candidates.empty:
        return candidates
    annotated = candidates.copy()
    annotated["mean_rsa"] = pd.NA
    annotated["mean_plddt"] = pd.NA
    annotated["alphafold_struc_flag"] = pd.NA

    cache: Dict[str, pd.DataFrame] = {}
    for idx, row in annotated.iterrows():
        accession = row["accession"]
        if accession not in cache:
            cache[accession] = _load_rsa_scores(output_dir, input_stem, accession)
        per_residue = cache[accession]
        if per_residue.empty:
            continue
        window = per_residue[(per_residue["position"] >= row["start"]) & (per_residue["position"] <= row["end"])]
        if window.empty:
            continue
        annotated.at[idx, "mean_rsa"] = window["rsa"].mean(skipna=True)
        annotated.at[idx, "mean_plddt"] = window["plddt"].mean(skipna=True)
        flags = window["alphafold_struc_flag"].dropna().unique()
        if len(flags) == 1:
            annotated.at[idx, "alphafold_struc_flag"] = flags[0]
        elif len(flags) > 1:
            annotated.at[idx, "alphafold_struc_flag"] = "mixed"

    return annotated


# Patron oficial de UniProt para accession numbers (formato de 6 o 10
# caracteres, ver https://www.uniprot.org/help/accession_numbers) -- no
# inventado, es el regex publicado por UniProt para validar sus propios IDs.
_UNIPROT_ACCESSION_TOKEN_RE = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$"
)


def _extract_uniprot_accession(accession_stem: str) -> Optional[str]:
    """Busca un token con forma de accession UniProt real dentro del stem de Project 1.

    El ``accession`` de Project 1 es el stem completo del fichero de entrada
    (p.ej. ``HIV1_gp41_full_P04578_512-856``), no un ID UniProt limpio.
    Parte el stem por cualquier separador no alfanumerico (``_``/``-``/``.``)
    y devuelve el primer token que calza EXACTAMENTE (ancla completa, no
    subcadena) el patron oficial de UniProt -- evita falsos positivos como
    matches parciales dentro de un token mas largo. ``None`` si ningun token
    calza (caso comun cuando el accession de entrada no viene de UniProt).
    """
    for token in re.split(r"[^A-Za-z0-9]+", accession_stem):
        if _UNIPROT_ACCESSION_TOKEN_RE.match(token):
            return token
    return None


def _annotate_from_uniprot(candidates: pd.DataFrame, column: str, lookup) -> pd.DataFrame:
    """Patron comun a las anotaciones informativas via UniProt de este modulo.

    Extraido tras comprobar que ``_annotate_secretory_pathway`` y la nueva
    anotacion de estructura cuaternaria (``_annotate_quaternary_structure``)
    repetian exactamente el mismo esqueleto: extraer el accession UniProt del
    stem (``_extract_uniprot_accession``), cachear por stem (una unica
    consulta HTTP por accession UniProt distinto, nunca por fila/candidata --
    es una propiedad de la proteina completa), y degradar a ``None`` sin
    abortar el resto de la tabla ante cualquier fallo (stem sin token
    UniProt, accession no reconocido, o fallo real de red --
    ``UniProtLookupError``, logueado como warning). No requiere
    ``output_dir``/``input_stem``: a diferencia de RSA/TMbed, no depende del
    cache por-fase de Project 1.

    Args:
        column: nombre de la columna nueva a anadir.
        lookup: funcion ``uniprot_id -> valor`` (puede lanzar
            ``UniProtLookupError``); la logica biologica especifica de cada
            anotacion vive aqui, no en este helper.
    """
    if candidates.empty:
        return candidates
    annotated = candidates.copy()
    annotated[column] = pd.NA

    cache: Dict[str, object] = {}
    for idx, row in annotated.iterrows():
        stem = row["accession"]
        if stem not in cache:
            uniprot_id = _extract_uniprot_accession(stem)
            if uniprot_id is None:
                cache[stem] = None
            else:
                try:
                    cache[stem] = lookup(uniprot_id)
                except UniProtLookupError as exc:
                    logger.warning(
                        "No se pudo consultar UniProt para '%s' (extraido de '%s') -- "
                        "'%s' queda sin determinar (no fatal): %s",
                        uniprot_id, stem, column, exc,
                    )
                    cache[stem] = None
        annotated.at[idx, column] = cache[stem]

    return annotated


def _annotate_secretory_pathway(candidates: pd.DataFrame) -> pd.DataFrame:
    """Anota (sin filtrar) evidencia real de UniProt de que el antigeno transita la via secretora.

    Puramente informativo, mismo patron no-decisorio que
    ``_annotate_rsa_accessibility``/``conservation_pct``/``documented_region``:
    complementa la RSA (accesibilidad DENTRO del plegamiento modelado) con
    una senal a nivel de proteina completa -- si el antigeno ni siquiera es
    plausible como diana de superficie (secretado/membrana/extracelular),
    ningun RSA local lo vuelve una diana valida para un nanobody. Reusa
    ``lookup_secretory_pathway_evidence`` de PTM-Prediction sin cambios (ver
    `_load_uniprot_localization_client`). Ver ``_annotate_from_uniprot`` para
    el patron de cacheo/degradacion compartido.
    """
    return _annotate_from_uniprot(candidates, "secretory_pathway_evidence", lookup_secretory_pathway_evidence)


def _fetch_uniprot_entry_json(accession: str) -> Optional[dict]:
    """GET minimo a UniProt REST, mismo endpoint/timeout/manejo de errores que el cliente de PTM-Prediction.

    No existe en ``uniprot_localization_client.py`` una funcion de fetch
    reutilizable por separado (la logica HTTP esta embebida dentro de
    ``lookup_secretory_pathway_evidence``) -- reescribir aqui el mismo
    mecanismo GET/timeout/manejo de 400-404 es duplicar ~10 lineas de
    mecanica HTTP sin logica biologica que pueda divergir en silencio (a
    diferencia del resto de este archivo, que SI reusa por import todo lo
    que tiene logica real que proteger), preferible a modificar el cliente
    de un proyecto hermano para una necesidad especifica de Proyecto 3.
    Reusa la constante real ``UNIPROT_URL_TEMPLATE``/``UNIPROT_TIMEOUT_SECONDS``
    del modulo cargado en vez de duplicar el string del endpoint.

    Returns:
        El JSON completo de la entrada, o ``None`` si UniProt no reconoce el
        accession (HTTP 400/404 -- mismo criterio que el cliente original).

    Raises:
        UniProtLookupError: fallo real de red o respuesta no parseable.
    """
    import json
    import urllib.error
    import urllib.request

    url = _uniprot_client.UNIPROT_URL_TEMPLATE.format(accession=accession)
    request = urllib.request.Request(url, headers={"Accept": "application/json"})
    try:
        with urllib.request.urlopen(request, timeout=_uniprot_client.UNIPROT_TIMEOUT_SECONDS) as response:
            return json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        if exc.code in (400, 404):
            return None
        raise UniProtLookupError(f"UniProt devolvio HTTP {exc.code} para '{accession}'.") from exc
    except urllib.error.URLError as exc:
        raise UniProtLookupError(f"No se pudo contactar a UniProt para '{accession}': {exc.reason}") from exc
    except (json.JSONDecodeError, UnicodeDecodeError) as exc:
        raise UniProtLookupError(f"Respuesta de UniProt no parseable como JSON para '{accession}': {exc}") from exc


def _extract_quaternary_structure_note(entry_json: Optional[dict]) -> Optional[str]:
    """Texto libre curado por UniProt sobre el ensamblaje real de la proteina (comentario ``SUBUNIT``).

    Verificado 2026-08-29 contra datos reales: P04578 (HIV-1 Env, gp160)
    devuelve DOS comentarios ``SUBUNIT`` (uno por ``molecule`` procesada,
    gp120 y gp41), ambos describiendo el mismo homotrimero real de
    heterodimeros no covalentes -- se concatenan (deduplicados, orden
    preservado) en vez de quedarse solo con el primero, porque no hay
    garantia de que el candidato final caiga en la molecula del primer
    comentario. Es texto libre, no un campo estructurado con vocabulario
    controlado (a diferencia de los ``keywords`` de la via secretora) --
    exponerlo tal cual, sin intentar clasificarlo en un enum
    monomero/dimero/trimero, evita fabricar una regla de parseo de lenguaje
    natural sin fuente verificada. ``None`` si la entrada no tiene ningun
    comentario ``SUBUNIT`` (caso real verificado con P01588/EPO).
    """
    if not entry_json:
        return None
    texts = []
    for comment in entry_json.get("comments", []):
        if comment.get("commentType") != "SUBUNIT":
            continue
        for text_entry in comment.get("texts", []):
            value = text_entry.get("value")
            if value:
                texts.append(value)
    if not texts:
        return None
    seen = set()
    unique_texts = []
    for value in texts:
        if value not in seen:
            seen.add(value)
            unique_texts.append(value)
    return " || ".join(unique_texts)


def _annotate_quaternary_structure(candidates: pd.DataFrame) -> pd.DataFrame:
    """Anota (sin filtrar) el comentario SUBUNIT real de UniProt sobre el ensamblaje del antigeno (Nivel A).

    Ni la RSA (diferencia 4, accesibilidad DENTRO del plegamiento modelado)
    ni la via secretora (diferencia 5) dicen si la region elegida queda
    enterrada en una interfaz de oligomerizacion real -- un antigeno
    documentado como homotrimero (p.ej. Env de HIV) puede ocultar en el
    complejo nativo una ventana que en una estructura/prediccion
    MONOMERICA parece perfectamente expuesta. Este es el "Nivel A" (barato,
    generico, sin curacion manual por target) de las opciones discutidas
    para robustecer Fase 0 de cara a targets arbitrarios (no solo el panel
    GP120/MonkeyPox, tambien candidatos de competicion como Adaptyv):
    texto libre de UniProt, informativo, funciona para cualquier accession
    con entrada UniProt sin conocer de antemano el estado oligomerico real.

    El "Nivel B" (recalcular RSA sobre el ensamblaje biologico real via la
    API de estructura de RCSB cuando el Camino es PDB) queda deliberadamente
    fuera de esta funcion -- da una verificacion estructural real en vez de
    una nota de texto, pero exige mapear asym_id/entidad al chain elegido
    como antigeno, solo aplica a estructuras PDB depositadas (no a modelos
    de AlphaFold ni a targets sin estructura resuelta), y no esta
    implementado todavia (ver investigacion de vialidad, sesion 2026-08-29).

    Ver ``_annotate_from_uniprot`` para el patron de cacheo/degradacion
    compartido con ``_annotate_secretory_pathway``.
    """
    return _annotate_from_uniprot(
        candidates,
        "quaternary_structure_note",
        lambda uniprot_id: _extract_quaternary_structure_note(_fetch_uniprot_entry_json(uniprot_id)),
    )


class RcsbLookupError(Exception):
    """Fallo al consultar la API de estructuras de RCSB: red, timeout, o respuesta no parseable."""


_RCSB_TIMEOUT_SECONDS = 20

# Formato oficial wwPDB de PDB ID (4 caracteres): 1 digito (1-9) + 3
# alfanumericos. Nunca colisiona con `_UNIPROT_ACCESSION_TOKEN_RE` (que
# siempre exige letra inicial) -- verificado contra los 17 stems reales del
# panel de validacion de Project 1 (convencion `Descripcion_PDBID`, p.ej.
# `SARSCoV2_RBD_7BEP`).
_PDB_ID_TOKEN_RE = re.compile(r"^[1-9][A-Za-z0-9]{3}$")


def _extract_pdb_id(accession_stem: str) -> Optional[str]:
    """Busca un token con forma de PDB ID real dentro del stem, mismo criterio que ``_extract_uniprot_accession``.

    Riesgo real aceptado (igual que con el accession UniProt): un stem de un
    target futuro sin convencion `Descripcion_PDBID` (p.ej. un candidato de
    Adaptyv nombrado de otra forma) podria, en teoria, contener un token de
    4 caracteres que calce el patron sin ser un PDB ID real -- la consulta a
    RCSB simplemente devolveria 404 en ese caso (degradacion no fatal, ver
    `_fetch_rcsb_json`), nunca un falso positivo silencioso.
    """
    for token in re.split(r"[^A-Za-z0-9]+", accession_stem):
        if _PDB_ID_TOKEN_RE.match(token):
            return token.upper()
    return None


def _fetch_rcsb_json(url: str) -> Optional[dict]:
    """GET minimo a la API de datos de RCSB, mismo criterio de errores que ``_fetch_uniprot_entry_json``."""
    import json
    import urllib.error
    import urllib.request

    request = urllib.request.Request(url, headers={"Accept": "application/json"})
    try:
        with urllib.request.urlopen(request, timeout=_RCSB_TIMEOUT_SECONDS) as response:
            return json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            return None
        raise RcsbLookupError(f"RCSB devolvio HTTP {exc.code} para '{url}'.") from exc
    except urllib.error.URLError as exc:
        raise RcsbLookupError(f"No se pudo contactar a RCSB para '{url}': {exc.reason}") from exc
    except (json.JSONDecodeError, UnicodeDecodeError) as exc:
        raise RcsbLookupError(f"Respuesta de RCSB no parseable como JSON para '{url}': {exc}") from exc


def _resolve_label_asym_id(pdb_id: str, auth_chain_id: str) -> Optional[str]:
    """Resuelve auth_asym_id (``chain_id`` de Project 1) -> label_asym_id (asym_id interno mmCIF de RCSB).

    Paso NECESARIO, no cosmetico: verificado con datos reales que ambos
    pueden diferir (``8FDD`` real, del propio panel de Project 1: la
    entidad 3 tiene auth ``'P'`` con label real ``'C'``).

    CORRECCION REAL sobre una primera version de esta funcion (misma sesion
    2026-08-29): se asumio que ``auth_asym_ids``/``asym_ids`` (las dos listas
    que devuelve ``polymer_entity`` a nivel de entidad) eran arrays
    PARALELOS del mismo orden -- verificado FALSO contra datos reales de
    ``7BEP`` (entidad con 2 cadenas: ``auth_asym_ids=['C','E']``,
    ``asym_ids=['A','F']``; emparejar por indice daria ``E -> F``, pero el
    endpoint de instancia -- fuente autoritativa, ver abajo -- confirma que
    el ``auth_asym_id`` real de ``E`` es ``E``, no ``F``). El esquema de
    RCSB no garantiza ese alineamiento posicional.

    Metodo correcto, verificado: ``polymer_entity_instance/{pdb_id}/{id}``
    esta indexado por ``label_asym_id`` (NO por auth, verificado con 8FDD:
    consultar con el label real ``'C'`` funciona y devuelve
    ``auth_asym_id='P'``; consultar directamente con el auth ``'P'`` da
    HTTP 404). Por eso aqui se enumera cada ``label_asym_id`` real de cada
    entidad polimerica (``asym_ids`` de ``polymer_entity``) y se consulta su
    instancia una por una hasta encontrar la que reporte el
    ``auth_chain_id`` buscado -- mas llamadas HTTP que un simple indexado,
    pero es el unico camino verificado como correcto.

    ``None`` si el PDB ID no existe o ninguna instancia reporta ese
    ``auth_chain_id``.
    """
    entry = _fetch_rcsb_json(f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}")
    if entry is None:
        return None
    entity_ids = entry.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids") or []
    for entity_id in entity_ids:
        entity = _fetch_rcsb_json(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}")
        if entity is None:
            continue
        ids = entity.get("rcsb_polymer_entity_container_identifiers", {})
        for label_asym_id in ids.get("asym_ids") or []:
            instance = _fetch_rcsb_json(
                f"https://data.rcsb.org/rest/v1/core/polymer_entity_instance/{pdb_id}/{label_asym_id}"
            )
            if instance is None:
                continue
            instance_ids = instance.get("rcsb_polymer_entity_instance_container_identifiers", {})
            if instance_ids.get("auth_asym_id") == auth_chain_id:
                return label_asym_id
    return None


def _lookup_assembly_copy_number(pdb_id: str, auth_chain_id: str) -> Optional[Tuple[int, str, str]]:
    """Cuenta cuantas veces aparece la cadena elegida en el ensamblaje biologico real que SI la contiene (RCSB).

    CAVEAT REAL (conceptual, ver ``_resolve_label_asym_id`` para el bug real
    de mapeo ya corregido en esta misma sesion): el "ensamblaje biologico"
    que RCSB expone puede seguir reflejando copias CRISTALOGRAFICAS de un
    complejo completo (p.ej. 2 copias independientes de un complejo
    anticuerpo-antigeno en la unidad asimetrica) en vez de auto-
    oligomerizacion biologica real del antigeno -- verificado que ``7BEP``
    (RBD de SARS-CoV-2, panel real de Project 1) da ``copy_number=1`` con la
    resolucion correcta (coherente con que la RBD no es un oligomero
    obligado), pero no hay garantia de que otra entrada con un mapeo igual
    de correcto no muestre igualmente redundancia cristalografica real. Por
    esto esta anotacion es PURAMENTE INFORMATIVA (nunca decide/filtra, mismo
    patron que el resto de Fase 0) y debe leerse junto a
    ``quaternary_structure_note`` (Nivel A, texto curado de UniProt) antes
    de concluir que hay un oligomero real -- ninguna de las dos senales por
    si sola es suficiente.

    CORRECCION REAL sobre una primera version (misma sesion 2026-08-29,
    encontrada corriendo Fase 0 de principio a fin contra ``3NGB``, el
    propio gp120 de Proyecto 3): se asumio que ``assembly_id='1'`` bastaba
    (RCSB lo describe como el ensamblaje "preferido"), pero ``3NGB`` tiene 4
    ensamblajes reales (``entry.assembly_ids``), cada uno una copia
    independiente de un complejo Fab-gp120 distinto -- la cadena elegida
    (auth ``'D'``, label ``'G'``) solo aparece en el ensamblaje 3, NO en el
    1. Asumir solo "1" habria devuelto ``None`` (degradacion no-fatal, no un
    resultado incorrecto, pero incompleta) para un caso real y relevante del
    propio target de Proyecto 3. Ahora se enumeran TODOS los
    ``assembly_ids`` de la entrada y se usa el primero donde la cadena
    aparezca en algun cluster. Prefiere, dentro de ese ensamblaje, el
    ``rcsb_struct_symmetry`` de ``kind='Global Symmetry'`` si existe
    (verificado que puede haber varios para el mismo ensamblaje, p.ej.
    tambien ``'Pseudo Symmetry'`` -- caso real ``4NCO``); si la cadena no
    aparece en ninguno de ese tipo, cae al primero disponible donde SI
    aparezca.

    Returns:
        ``(copy_number, oligomeric_state, symmetry_kind)``, o ``None`` si no
        se pudo determinar (PDB ID no reconocido, cadena no resuelta a
        label_asym_id, o la cadena no aparece en ningun cluster de NINGUN
        ensamblaje de la entrada).
    """
    label_asym_id = _resolve_label_asym_id(pdb_id, auth_chain_id)
    if label_asym_id is None:
        return None

    entry = _fetch_rcsb_json(f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}")
    if entry is None:
        return None
    assembly_ids = entry.get("rcsb_entry_container_identifiers", {}).get("assembly_ids") or ["1"]

    for assembly_id in assembly_ids:
        assembly = _fetch_rcsb_json(f"https://data.rcsb.org/rest/v1/core/assembly/{pdb_id}/{assembly_id}")
        if assembly is None:
            continue
        symmetry_entries = assembly.get("rcsb_struct_symmetry") or []
        ordered = sorted(symmetry_entries, key=lambda s: s.get("kind") != "Global Symmetry")
        found = _find_cluster_copy_number(ordered, label_asym_id)
        if found is not None:
            return found
    return None


def _find_cluster_copy_number(
    symmetry_entries: List[dict], label_asym_id: str
) -> Optional[Tuple[int, str, str]]:
    """Busca ``label_asym_id`` en los clusters de ``symmetry_entries`` (ya ordenados por preferencia de kind)."""
    for sym in symmetry_entries:
        for cluster in sym.get("clusters", []):
            members = [member.get("asym_id") for member in cluster.get("members", [])]
            if label_asym_id in members:
                return len(members), sym.get("oligomeric_state"), sym.get("kind")
    return None


def _load_chain_id(output_dir: Path, accession: str) -> Optional[str]:
    """Lee el ``chain_id`` (auth) elegido por Fase 1.5 para ``accession``, desde su position_mapping cacheado.

    ``{output_dir}/{accession}_position_mapping.csv`` (escrito por
    ``structure_parser.parse_structure``) tiene un unico ``chain_id`` real
    para todas sus filas (una cadena elegida por accession, ver docstring de
    ``StructureRecord``) -- basta leer la primera fila. ``None`` si el
    archivo no existe (Camino FASTA, o el accession no tiene estructura
    asociada) o esta vacio.
    """
    path = output_dir / f"{accession}_position_mapping.csv"
    if not path.is_file():
        return None
    mapping = pd.read_csv(path, nrows=1)
    if mapping.empty or "chain_id" not in mapping.columns:
        return None
    return str(mapping.at[0, "chain_id"])


def _annotate_assembly_copy_number(candidates: pd.DataFrame, output_dir: Path) -> pd.DataFrame:
    """Anota (sin filtrar) cuantas copias reales tiene la cadena elegida en el ensamblaje biologico (Nivel B).

    Complementa ``quaternary_structure_note`` (Nivel A, texto libre de
    UniProt) con una verificacion ESTRUCTURAL real cuando el Camino es PDB:
    RCSB expone la composicion del ensamblaje biologico real
    (``rcsb_struct_symmetry``) de forma estructurada y generica para
    cualquier PDB ID, sin curacion manual por target -- viabilidad
    investigada y verificada con datos reales (sesion 2026-08-29): ``4NCO``
    (trimero Env real, del tipo de estructura de este proyecto) da
    ``assembly_copy_number=3`` para la cadena del antigeno; ``1LYZ``
    (monomero conocido) da ``1``.

    Solo cubre estructuras PDB depositadas -- ``None`` para Camino FASTA,
    modelos de AlphaFold, o cualquier target sin ID de PDB real en el stem
    (esperable en un candidato de competicion sin estructura resuelta, p.ej.
    Adaptyv); en ese caso ``quaternary_structure_note`` queda como unica
    senal de estructura cuaternaria.

    CAVEAT REAL (conceptual), ver docstring de
    ``_lookup_assembly_copy_number``: ``assembly_copy_number`` puede, en
    principio, contar copias CRISTALOGRAFICAS de un complejo anticuerpo-
    antigeno en vez de auto-oligomerizacion biologica real -- PURAMENTE
    INFORMATIVO, cruzar siempre con ``quaternary_structure_note`` antes de
    descartar una ventana por esto.

    Una unica resolucion (label_asym_id + ensamblaje) por accession
    distinto, cacheada -- no se reusa ``_annotate_from_uniprot`` porque la
    fuente es RCSB (no UniProt) y la clave de busqueda necesita ademas el
    ``chain_id`` elegido por Fase 1.5, no solo el accession.
    """
    if candidates.empty:
        return candidates
    annotated = candidates.copy()
    annotated["assembly_copy_number"] = pd.NA
    annotated["assembly_oligomeric_state"] = pd.NA
    annotated["assembly_symmetry_kind"] = pd.NA

    cache: Dict[str, Optional[Tuple[int, str, str]]] = {}
    for idx, row in annotated.iterrows():
        stem = row["accession"]
        if stem not in cache:
            pdb_id = _extract_pdb_id(stem)
            chain_id = _load_chain_id(output_dir, stem) if pdb_id is not None else None
            if pdb_id is None or chain_id is None:
                cache[stem] = None
            else:
                try:
                    cache[stem] = _lookup_assembly_copy_number(pdb_id, chain_id)
                except RcsbLookupError as exc:
                    logger.warning(
                        "No se pudo consultar RCSB para '%s' (PDB ID '%s', cadena '%s') -- "
                        "'assembly_copy_number' queda sin determinar (no fatal): %s",
                        stem, pdb_id, chain_id, exc,
                    )
                    cache[stem] = None
        result = cache[stem]
        if result is not None:
            copy_number, oligomeric_state, symmetry_kind = result
            annotated.at[idx, "assembly_copy_number"] = copy_number
            annotated.at[idx, "assembly_oligomeric_state"] = oligomeric_state
            annotated.at[idx, "assembly_symmetry_kind"] = symmetry_kind

    return annotated


def _refresh_bnab_iedb_annotations(
    selected: pd.DataFrame,
    bnab_df: Optional[pd.DataFrame],
    iedb_df: Optional[pd.DataFrame],
    bnab_min_overlap: int,
    iedb_min_overlap: int,
) -> pd.DataFrame:
    """Recalcula bnAb/IEDB sobre la secuencia FINAL de ``selected`` (post pad/trim/self-tolerancia).

    No relee ningun archivo de referencia: reusa los matches YA encontrados
    por Fase 6/6c de Project 1 contra la region padre mas ancha (``bnab_df``/
    ``iedb_df``, columnas de ``lanl_catnap_engine.query_bnab_crossref``/
    ``iedb_engine.query_iedb_crossref``), y re-verifica el solapamiento real
    contra la ventana final ya recortada/extendida. Correcto por construccion
    para candidatos RECORTADOS (subcadena del padre: cualquier match posible
    en la sub-ventana ya estaba en el conjunto padre). Limitacion real y
    documentada para candidatos con PADDING (flancos nativos anadidos,
    contenido nuevo que Fase 6/6c nunca vio): un match que solo exista en el
    flanco anadido no se detectaria aqui -- arreglarlo exigiria releer los
    CSV de referencia completos, que es exactamente el costo de E/S que este
    diseño evita. Caso limite, no bloqueante: el padding solo aplica a
    candidatos <15 aa (``CONSTRUCT_BCELL_FLANK_THRESHOLD``).
    """
    if selected.empty:
        return selected
    selected = selected.copy()

    if bnab_df is not None:
        ref = bnab_df.drop_duplicates(subset="epitope_sequence") if not bnab_df.empty else bnab_df
        n_matches, best_ic50, max_nvir, any_neut = [], [], [], []
        for seq in selected["sequence"]:
            m = _matches_within_overlap(seq, ref, bnab_min_overlap)
            n_matches.append(len(m))
            best_ic50.append(m["catnap_mean_ic50"].min() if not m.empty else pd.NA)
            max_nvir.append(m["catnap_n_viruses"].max() if not m.empty else pd.NA)
            any_neut.append(bool(m["neutralizing"].fillna(False).any()) if not m.empty else False)
        selected["bnab_n_matches"] = n_matches
        selected["bnab_best_ic50"] = best_ic50
        selected["bnab_max_n_viruses"] = max_nvir
        selected["bnab_any_neutralizing"] = any_neut

    if iedb_df is not None:
        ref = iedb_df.drop_duplicates(subset="epitope_sequence") if not iedb_df.empty else iedb_df
        documented, n_matches, best_qual, pmids = [], [], [], []
        for seq in selected["sequence"]:
            m = _matches_within_overlap(seq, ref, iedb_min_overlap)
            documented.append(not m.empty)
            n_matches.append(len(m))
            best_qual.append(_strongest_qualitative_measure(m["qualitative_measure"]) if not m.empty else pd.NA)
            pmids.append(",".join(sorted({str(p) for p in m["pmid"].dropna()})) if not m.empty else "")
        selected["documented_region"] = documented
        selected["iedb_n_matches"] = n_matches
        selected["iedb_best_qualitative_measure"] = best_qual
        selected["iedb_pmids"] = pmids

    return selected


def select_bcell_candidates_for_boltzgen(
    safe_df: pd.DataFrame,
    algpred_df: pd.DataFrame,
    stackgly_df: pd.DataFrame,
    top_n: int,
    conservation_map: Optional[Dict[str, float]] = None,
    bnab_df: Optional[pd.DataFrame] = None,
    iedb_df: Optional[pd.DataFrame] = None,
    output_dir: Optional[Path] = None,
    input_stem: Optional[str] = None,
    max_length: Optional[int] = None,
    flank_threshold: Optional[int] = None,
    flank_padding: Optional[int] = None,
    blast_db: str = Settings.BLAST_HUMAN_DB,
    identity_threshold: float = Settings.BLAST_IDENTITY_THRESHOLD,
    bnab_min_overlap: int = Settings.LANL_CATNAP_MIN_OVERLAP,
    iedb_min_overlap: int = Settings.IEDB_MIN_OVERLAP,
) -> pd.DataFrame:
    """Fork corregido de ``construct_assembly._select_bcell_candidates`` para Fase 0 de Proyecto 3.

    Mismo ranking/top-N/padding que el original de Project 1 (reusado sin
    cambios via import). Difiere en 7 puntos, ver docstring del modulo:
    recorte con re-chequeo real de alergenicidad/glicosilacion, bnAb/IEDB/
    ``documented_region`` recalculados sobre la secuencia final, filtro duro
    por TMbed, anotacion de RSA/pLDDT, anotacion de evidencia de via
    secretora (UniProt, cruzado con PTM-Prediction), anotacion de estructura
    cuaternaria via texto libre de UniProt (Nivel A), y anotacion de copias
    reales en el ensamblaje biologico via RCSB (Nivel B).
    """
    if safe_df.empty:
        return safe_df

    allergen_seqs = set(algpred_df[algpred_df["algpred_veredicto"] == "Allergen"]["sequence"]) \
        if not algpred_df.empty else set()
    glyco_risky_seqs = set(stackgly_df[stackgly_df["stackglyembed_veredicto"] == "Glicosilado"]["sequence"]) \
        if not stackgly_df.empty else set()

    candidates = safe_df.copy()
    if candidates.empty:
        return candidates
    candidates["allergen"] = candidates["sequence"].isin(allergen_seqs)
    candidates["glycosylated"] = candidates["sequence"].isin(glyco_risky_seqs)
    if conservation_map:
        candidates["conservation_pct"] = candidates["sequence"].map(conservation_map)

    score_cols = [c for c in candidates.columns if c.endswith("_score")]
    if score_cols:
        percentiles = candidates[score_cols].rank(pct=True, na_option="keep")
        candidates["_rank_score"] = percentiles.mean(axis=1, skipna=True)
    else:
        candidates["_rank_score"] = 0.0
    candidates = candidates.sort_values("_rank_score", ascending=False)
    selected = candidates.head(top_n).drop(columns="_rank_score")

    if output_dir is not None and input_stem is not None:
        selected = _pad_short_bcell_candidates(
            selected, output_dir, input_stem,
            flank_threshold or Settings.CONSTRUCT_BCELL_FLANK_THRESHOLD,
            flank_padding or Settings.CONSTRUCT_BCELL_FLANK_PADDING,
        )
        selected = _trim_long_bcell_candidates_fresco(
            selected, output_dir, input_stem, max_length or Settings.CONSTRUCT_BCELL_MAX_LENGTH
        )
        selected = _filter_tmbed_masked(selected, output_dir, input_stem)
        selected = _annotate_conformational_coverage(selected, output_dir, input_stem)
        selected = _annotate_rsa_accessibility(selected, output_dir, input_stem)
        selected = _annotate_secretory_pathway(selected)
        selected = _annotate_quaternary_structure(selected)
        selected = _annotate_assembly_copy_number(selected, output_dir)
        selected = filter_self_tolerant(selected, "sequence", db_path=blast_db, identity_threshold=identity_threshold)

    selected = _refresh_bnab_iedb_annotations(selected, bnab_df, iedb_df, bnab_min_overlap, iedb_min_overlap)
    return selected
