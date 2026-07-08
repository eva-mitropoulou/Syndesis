from __future__ import annotations

from typing import Any

from Bio.PDB import NeighborSearch
from Bio.PDB.Polypeptide import is_aa

from egfr_dockingforge.stage1.ligand_extraction import LigandInstance, atom_element, heavy_atoms


def bool_reason(flag: bool, reason: str | None = None) -> str | None:
    return reason if flag else None


def smiles_has_warhead(smiles: str | None) -> tuple[bool, str | None]:
    if not smiles:
        return False, None
    normalized = smiles.replace("[", "").replace("]", "").lower()
    if "c=cc(=o)n" in normalized or "c=cc(=o)" in normalized:
        return True, "Michael-acceptor/acrylamide-like motif in ligand SMILES."
    if "c#cc(=o)n" in normalized or "c#cc(=o)" in normalized:
        return True, "Propiolamide-like motif in ligand SMILES."
    return False, None


def struct_conn_covalent_evidence(mmcif: dict[str, Any], ligand_comp_id: str) -> list[str]:
    conn_types = mmcif.get("_struct_conn.conn_type_id", [])
    ptnr1 = mmcif.get("_struct_conn.ptnr1_label_comp_id", [])
    ptnr2 = mmcif.get("_struct_conn.ptnr2_label_comp_id", [])
    if isinstance(conn_types, str):
        conn_types = [conn_types]
        ptnr1 = [ptnr1]
        ptnr2 = [ptnr2]
    evidence: list[str] = []
    for conn_type, comp1, comp2 in zip(conn_types, ptnr1, ptnr2, strict=False):
        if str(conn_type).lower().startswith("covale") and ligand_comp_id in {
            str(comp1).upper(),
            str(comp2).upper(),
        }:
            evidence.append(f"struct_conn covalent linkage involving {ligand_comp_id}")
    return evidence


def suspicious_cys_distance(chain: Any, ligand_residue_id: tuple[str, int, str]) -> str | None:
    ligand_residue = None
    cys_atoms = []
    for residue in chain.get_residues():
        if residue.id == ligand_residue_id:
            ligand_residue = residue
        if residue.get_resname().strip().upper() == "CYS":
            for atom in residue.get_atoms():
                if atom.get_name().strip().upper() == "SG":
                    cys_atoms.append(atom)
    if ligand_residue is None or not cys_atoms:
        return None
    ligand_atoms = heavy_atoms(ligand_residue)
    if not ligand_atoms:
        return None
    search = NeighborSearch(cys_atoms)
    for atom in ligand_atoms:
        neighbors = search.search(atom.coord, 2.2, level="A")
        if neighbors:
            return "Ligand heavy atom within 2.2 A of cysteine SG; likely covalent or requires review."
    return None


def active_site_completeness(chain: Any, _config: dict[str, Any]) -> tuple[list[str], list[str], float]:
    # KNOWN LIMITATION: Stage 1 only has a raw Bio.PDB chain here (deposited auth
    # numbering, no UniProt->auth residue-number mapping). Position-level verification
    # of the specific catalytic residues (e.g. hinge Met793, gatekeeper Thr790,
    # catalytic Lys745, DFG Asp855) requires the residue-number mapper from
    # stage2/pocket_mapping.py (resolve_uniprot_residue / residue_by_auth_seq), which is
    # not available at this stage. We therefore do NOT claim position-level completeness.
    #
    # Honest fallback: for each key residue CLASS we require that a standard amino acid
    # residue of the expected resname is present AND has a CA atom resolved (backbone
    # modelled), rather than merely that the residue type appears somewhere in the chain.
    # Completeness is the fraction of key classes with at least one such CA-resolved
    # residue. This is a lower-bound proxy and does not confirm the correct sequence
    # position; the position-accurate check lives in Stage 2 pocket mapping.
    resname_has_ca: dict[str, bool] = {}
    for residue in chain.get_residues():
        if not is_aa(residue, standard=True):
            continue
        resname = residue.get_resname().strip().upper()
        has_ca = any(atom.get_name().strip().upper() == "CA" for atom in residue.get_atoms())
        resname_has_ca[resname] = resname_has_ca.get(resname, False) or has_ca
    key_classes = {
        "Lys745": "LYS",
        "Glu762": "GLU",
        "Thr790": "THR",
        "Met793": "MET",
        "Cys797": "CYS",
        "Asp855": "ASP",
        "Phe856": "PHE",
        "Gly857": "GLY",
    }
    missing = [label for label, resname in key_classes.items() if not resname_has_ca.get(resname, False)]
    missing_atoms: list[str] = []
    score = (len(key_classes) - len(missing)) / len(key_classes)
    return missing, missing_atoms, score


def classify_complex(
    pdb_id: str,
    ligand: LigandInstance,
    chain: Any,
    mmcif: dict[str, Any],
    ligand_smiles: str | None,
    config: dict[str, Any],
) -> dict[str, Any]:
    filters = config["filters"]
    warnings: list[str] = []
    known_covalent = {str(value).upper() for value in filters.get("known_covalent_ligand_comp_ids", [])}
    allosteric_ids = {str(value).upper() for value in filters.get("allosteric_ligand_comp_ids", [])}
    excluded = {str(value).upper() for value in filters.get("excluded_ligand_comp_ids", [])}

    warhead_flag, warhead_reason = smiles_has_warhead(ligand_smiles)
    covalent_evidence = struct_conn_covalent_evidence(mmcif, ligand.comp_id)
    if ligand.comp_id in known_covalent:
        covalent_evidence.append(f"{ligand.comp_id} listed as known/likely EGFR covalent ligand in config.")
    if warhead_flag:
        covalent_evidence.append(warhead_reason or "Reactive warhead motif.")

    covalent_flag = bool(covalent_evidence)
    allosteric_flag = ligand.comp_id in allosteric_ids
    atp_site_flag = (
        ligand.min_distance_to_chain is not None
        and ligand.min_distance_to_chain <= float(filters.get("max_ligand_to_pocket_distance_angstrom", 5.0))
    )
    ligand_class = "small_molecule"
    if ligand.comp_id in excluded:
        ligand_class = "cofactor_or_excluded_additive"
    elif ligand.heavy_atom_count < int(filters.get("fragment_heavy_atom_cutoff", 12)):
        ligand_class = "fragment"

    missing_residues, missing_atoms, completeness = active_site_completeness(chain, config)
    reasons: list[str] = []
    if ligand.comp_id in excluded:
        reasons.append(f"Ligand {ligand.comp_id} is excluded by v1 scope.")
    if ligand_class == "fragment":
        reasons.append("Fragment-only ligand below heavy-atom cutoff.")
    if covalent_flag:
        reasons.append("Covalent or likely covalent ligand.")
    if allosteric_flag:
        reasons.append("Known or configured allosteric ligand.")
    if not atp_site_flag:
        reasons.append("Ligand is not within configured ATP-pocket distance cutoff.")
    if completeness < 0.75:
        reasons.append("Severe missing active-site residue coverage by fallback check.")

    return {
        "covalent_flag": covalent_flag,
        "covalent_evidence": "; ".join(covalent_evidence) if covalent_evidence else None,
        "warhead_flag": warhead_flag,
        "warhead_reason": warhead_reason,
        "atp_site_flag": atp_site_flag,
        "atp_site_evidence": (
            f"Minimum ligand-protein distance {ligand.min_distance_to_chain:.2f} A"
            if ligand.min_distance_to_chain is not None
            else None
        ),
        "allosteric_flag": allosteric_flag,
        "ligand_class": ligand_class,
        "inhibitor_type_if_known": "noncovalent ATP-site" if atp_site_flag and not covalent_flag else None,
        "missing_active_site_residues": ",".join(missing_residues),
        "missing_active_site_atoms": ",".join(missing_atoms),
        "active_site_completeness_score": completeness,
        "hard_exclusion_reasons": reasons,
        "classification_warnings": warnings,
    }
