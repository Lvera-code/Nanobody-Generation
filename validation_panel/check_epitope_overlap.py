"""Comprueba si el epitopo DiscoTope usado para generar (RBD ATMSEQ 148..163 sobre
7BEP cadena E; PfCSP ATMSEQ 41..54 sobre 7RXP cadena A) coincide con el sitio de
union real (contacto antigeno-nanobody) en estructuras VHH-antigeno resueltas.

Identificacion de cadenas: por la descripcion real de entidad del propio mmCIF de
RCSB (_entity.pdbx_description via _entity_poly.pdbx_strand_id), no por longitud.
La cadena nanobody es la entidad cuya descripcion contiene un termino de
nanocuerpo (nanobody/VHH/sdAb/single-domain/camel*); el antigeno es la otra cadena
polimerica relevante. Corrige el bug real de la heuristica anterior ("cadena mas
corta = nanobody"), que fallaba en PfCSP porque el antigeno (dominio alfaTSR,
~65-67 residuos) es MAS CORTO que el nanobody (~125 residuos) ahi -- al reves que
en RBD.

Correspondencia de numeracion real via alineamiento por pares: ATMSEQ es una
posicion 1-based sobre los residuos resueltos del propio target usado para
disenar (ver scipion-chem-boltzgen/boltzgen/utils/label_indexing.py), NO el
auth_seq_id del PDB. Comparar esa posicion directamente contra el indice
posicional de cada estructura del panel (como hacia la version anterior de este
script) asume que ambas cadenas empiezan en el mismo residuo real y sin huecos --
falso en general entre cristales distintos (limites de construccion, residuos no
resueltos). Aqui se alinea por pares (Bio.Align.PairwiseAligner, BLOSUM62) la
secuencia real del target de referencia (7BEP cadena E / 7RXP cadena A, mismo
archivo usado para construir el design_spec de la campana real) contra la
secuencia real del antigeno de cada estructura del panel, y se traduce la ventana
ATMSEQ a la numeracion auth_seq_id real de cada estructura a traves de los bloques
alineados -- una posicion de la ventana que cae en un hueco de la alineacion
(residuo no resuelto/ausente en ese cristal) se reporta como tal, no se ignora en
silencio.

Contacto real: cualquier residuo del antigeno con un atomo a <=5A de cualquier
atomo de la cadena nanobody, indexado por su auth_seq_id real (no por posicion).
"""
import re
from pathlib import Path

import gemmi
from Bio.Align import PairwiseAligner, substitution_matrices

REFERENCE_TARGETS = {
    "rbd": ("reference_targets/7BEP_chainE_boltzgen_target.cif", range(148, 164)),
    "pfcsp": ("reference_targets/7RXP_chainA_boltzgen_target.cif", range(41, 55)),
}

TARGETS = {
    "rbd": Path("rbd"),
    "pfcsp": Path("pfcsp"),
}

NANOBODY_DESCRIPTION_RE = re.compile(
    r"nanobody|vhh|sdab|single[- ]domain|camel", re.IGNORECASE
)

ALIGNER = PairwiseAligner()
ALIGNER.substitution_matrix = substitution_matrices.load("BLOSUM62")
ALIGNER.mode = "global"
ALIGNER.open_gap_score = -10
ALIGNER.extend_gap_score = -0.5


def chain_entity_descriptions(cif_path):
    """Devuelve {chain_auth_id: entity_description} leyendo _entity/_entity_poly
    reales del mmCIF (no inferido/adivinado)."""
    doc = gemmi.cif.read(str(cif_path))
    block = doc.sole_block()
    descriptions = {
        row[0]: gemmi.cif.as_string(row[1])
        for row in block.find("_entity.", ["id", "pdbx_description"])
    }
    chain_to_desc = {}
    for entity_id, strand_ids in block.find("_entity_poly.", ["entity_id", "pdbx_strand_id"]):
        desc = descriptions.get(entity_id, "")
        for chain_id in strand_ids.split(","):
            chain_to_desc[chain_id.strip()] = desc
    return chain_to_desc


def polymer_chains(model):
    chains = []
    for chain in model:
        poly = chain.get_polymer()
        if poly.length() >= 20:
            chains.append((chain, poly.length()))
    return chains


def chain_sequence_and_residues(chain):
    """Secuencia real (1 letra/residuo) y lista de objetos Residue reales del
    mismo indice, en el orden en que aparecen en la cadena."""
    residues = list(chain.get_polymer())
    seq = gemmi.one_letter_code(chain.get_polymer().extract_sequence())
    assert len(seq) == len(residues)
    return seq.upper(), residues


def contact_residue_indices(antigen_residues, nanobody_chain, cutoff=5.0):
    """Indices (0-based, sobre antigen_residues) de residuos del antigeno con
    algun atomo a <=cutoff A de la cadena nanobody."""
    ns_atoms = [atom.pos for residue in nanobody_chain.get_polymer() for atom in residue]
    contacts = set()
    for idx, residue in enumerate(antigen_residues):
        for atom in residue:
            if any(atom.pos.dist(p) <= cutoff for p in ns_atoms):
                contacts.add(idx)
                break
    return contacts


def map_window_through_alignment(ref_seq, window_0based, panel_seq):
    """Alinea ref_seq contra panel_seq y devuelve:
    - {indice_ventana_ref_0based: indice_panel_0based} solo para posiciones que
      caen dentro de un bloque alineado (no en un hueco real de la alineacion).
    - porcentaje de identidad real sobre las columnas alineadas (sin huecos).
    """
    alignment = ALIGNER.align(ref_seq, panel_seq)[0]
    ref_blocks, panel_blocks = alignment.aligned

    position_map = {}
    identical, aligned_cols = 0, 0
    for (rs, re_), (ps, pe) in zip(ref_blocks, panel_blocks):
        for offset in range(re_ - rs):
            ref_i, panel_i = rs + offset, ps + offset
            position_map[ref_i] = panel_i
            aligned_cols += 1
            if ref_seq[ref_i] == panel_seq[panel_i]:
                identical += 1

    pct_identity = 100.0 * identical / aligned_cols if aligned_cols else 0.0
    mapped_window = {i: position_map[i] for i in window_0based if i in position_map}
    unmapped = [i for i in window_0based if i not in position_map]
    return mapped_window, pct_identity, unmapped


for label, target_path in TARGETS.items():
    ref_path, atmseq_window = REFERENCE_TARGETS[label]
    ref_st = gemmi.read_structure(ref_path)
    ref_st.setup_entities()
    ref_chain, _ = polymer_chains(ref_st[0])[0]
    ref_seq, _ = chain_sequence_and_residues(ref_chain)
    window_0based = [p - 1 for p in atmseq_window]  # ATMSEQ es 1-based

    print(f"\n=== {label.upper()} (ventana DiscoTope ATMSEQ {atmseq_window.start}..{atmseq_window.stop - 1}, "
          f"ref={ref_path}, {len(ref_seq)} res) ===")

    for cif_path in sorted(target_path.glob("*.cif")):
        try:
            st = gemmi.read_structure(str(cif_path))
            st.setup_entities()
            model = st[0]
            chains = polymer_chains(model)
            if len(chains) < 2:
                print(f"{cif_path.name}: menos de 2 cadenas polimericas utiles, se omite")
                continue

            chain_desc = chain_entity_descriptions(cif_path)
            nanobody_candidates = [
                (chain, length) for chain, length in chains
                if NANOBODY_DESCRIPTION_RE.search(chain_desc.get(chain.name, ""))
            ]
            if len(nanobody_candidates) != 1:
                print(f"{cif_path.name}: {len(nanobody_candidates)} cadenas coinciden con el "
                      f"patron de nanobody en _entity.pdbx_description (se esperaba 1) -- "
                      f"revisar a mano, se omite. Descripciones: {chain_desc}")
                continue
            nanobody_chain, nanobody_len = nanobody_candidates[0]
            antigen_candidates = [(c, l) for c, l in chains if c.name != nanobody_chain.name]
            if len(antigen_candidates) != 1:
                print(f"{cif_path.name}: {len(antigen_candidates)} cadenas de antigeno candidatas "
                      f"(se esperaba 1) -- revisar a mano, se omite. Descripciones: {chain_desc}")
                continue
            antigen_chain, antigen_len = antigen_candidates[0]

            panel_seq, panel_residues = chain_sequence_and_residues(antigen_chain)
            mapped_window, pct_identity, unmapped = map_window_through_alignment(
                ref_seq, window_0based, panel_seq
            )

            contacts_idx = contact_residue_indices(panel_residues, nanobody_chain)
            overlap_idx = sorted(i for i in mapped_window.values() if i in contacts_idx)
            overlap_authseq = [panel_residues[i].seqid.num for i in overlap_idx]

            unmapped_note = f" ventana_no_resuelta={len(unmapped)}/{len(window_0based)}" if unmapped else ""
            print(f"{cif_path.name}: antigeno={antigen_chain.name}({antigen_len}res, "
                  f"'{chain_desc.get(antigen_chain.name, '')}') "
                  f"nanobody={nanobody_chain.name}({nanobody_len}res, "
                  f"'{chain_desc.get(nanobody_chain.name, '')}') "
                  f"identidad_alineamiento={pct_identity:.0f}% "
                  f"contactos_reales={len(contacts_idx)} "
                  f"solapan_con_ventana={len(overlap_authseq)}{unmapped_note} "
                  f"auth_seq_id={overlap_authseq if overlap_authseq else ''}")
        except Exception as exc:
            print(f"{cif_path.name}: ERROR real al procesar -- {exc}")
