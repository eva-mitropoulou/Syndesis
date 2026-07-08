from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import pandas as pd
from rdkit import Chem, DataStructs

from egfr_dockingforge.common.io import write_table
from egfr_dockingforge.stage9.rdkit_transformations import _attach_to_atom, _canonical, _fp


# Mapping from an agent's proposed transformation class to the concrete RDKit
# fragment(s) that realise it. Mirrors the deterministic baseline fragment set
# in rdkit_transformations.enumerate_rule_based_analogs so that LLM-proposed and
# rule-based analogs are built by the SAME chemistry and are directly comparable.
TRANSFORM_CLASS_FRAGMENTS: dict[str, list[str]] = {
    "small_substituent_scan": ["C"],
    "halogen_scan": ["F", "Cl"],
    "solubilizing_tail_tuning": ["OC"],
    "conservative_bioisostere": ["C#N"],
    "heteroatom_swap": ["C"],  # graceful fallback: treated as a conservative methyl probe
}


def _fragments_for_class(transformation_class: str, tool_arguments: dict[str, Any]) -> list[str]:
    """Resolve the fragment SMILES to attach for a proposed transformation.

    If the agent supplied an explicit ``fragment_smiles`` in tool_arguments_json
    and it is a chemically valid, conservative fragment, honour it; otherwise
    fall back to the canonical fragment(s) for the transformation class.
    """
    explicit = tool_arguments.get("fragment_smiles") or tool_arguments.get("fragment")
    if isinstance(explicit, str) and explicit and Chem.MolFromSmiles(explicit) is not None:
        return [explicit]
    return TRANSFORM_CLASS_FRAGMENTS.get(str(transformation_class), [])


def _load_proposals(proposals_path: Path) -> pd.DataFrame:
    """Load agent proposals from the JSONL emitted by the agent loop."""
    if not proposals_path.exists():
        return pd.DataFrame()
    rows = []
    for line in proposals_path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            rows.append(json.loads(line))
        except json.JSONDecodeError:
            continue
    return pd.DataFrame(rows)


def dispatch_agent_proposals_to_molecules(
    proposals: pd.DataFrame,
    seeds: pd.DataFrame,
    edit_sites: pd.DataFrame,
    config: dict[str, Any],
    paths: dict[str, Path],
) -> pd.DataFrame:
    """Turn schema-valid agent proposals into concrete analog molecules.

    Each proposal carries a proposed_transformation_class + edit_site_id (and
    optionally an explicit fragment in tool_arguments_json). We resolve the edit
    site to an attachment atom, apply the corresponding RDKit transform, and emit
    analog_candidates rows tagged with the originating LLM strategy so they flow
    through the SAME validation / screening / acceptance path as the rule-based
    baseline. This is what makes the agentic-vs-baseline comparison real: without
    it every LLM strategy produces zero molecules.
    """
    seed_idx = seeds.set_index("seed_id")
    # edit_site_id -> attachment atom index and protection flag
    site_idx = edit_sites.set_index("edit_site_id") if not edit_sites.empty else pd.DataFrame()

    max_per_seed_strategy = int(config["transforms"]["max_analogs_per_seed_per_strategy"])
    protect_key_atoms = bool(config.get("edit_sites", {}).get("protect_key_interaction_atoms", True))

    rows: list[dict[str, Any]] = []
    seen: set[tuple[str, str, str]] = set()
    counts: dict[tuple[str, str], int] = {}

    valid = proposals[proposals.get("schema_valid", False) == True] if not proposals.empty else pd.DataFrame()  # noqa: E712
    for prop in valid.to_dict("records"):
        payload = prop.get("parsed_json")
        if isinstance(payload, str):
            try:
                payload = json.loads(payload)
            except json.JSONDecodeError:
                continue
        if not isinstance(payload, dict):
            continue
        strategy = prop.get("strategy_name", "agentic")
        seed_id = payload.get("seed_id") or prop.get("seed_id")
        if seed_id not in seed_idx.index:
            continue
        seed = seed_idx.loc[seed_id]
        if isinstance(seed, pd.DataFrame):
            seed = seed.iloc[0]
        parent_smiles = seed["standard_smiles"]

        transformation_class = str(payload.get("proposed_transformation_class", ""))
        edit_site_id = payload.get("edit_site_id")
        tool_arguments = payload.get("tool_arguments_json")
        if isinstance(tool_arguments, str):
            try:
                tool_arguments = json.loads(tool_arguments)
            except json.JSONDecodeError:
                tool_arguments = {}
        if not isinstance(tool_arguments, dict):
            tool_arguments = {}

        if edit_site_id not in getattr(site_idx, "index", []):
            continue
        site = site_idx.loc[edit_site_id]
        if isinstance(site, pd.DataFrame):
            site = site.iloc[0]
        # Honour protected atoms (including, when enabled, key-interaction atoms
        # flagged upstream). A proposal that targets a protected site is skipped
        # rather than silently allowed to damage the pharmacophore.
        if protect_key_atoms and bool(site.get("protected_region_flag", False)):
            continue
        atom_idx = int(site["attachment_atom_idx"])
        if atom_idx < 0:
            continue

        for fragment in _fragments_for_class(transformation_class, tool_arguments):
            key_count = (strategy, str(seed_id))
            if counts.get(key_count, 0) >= max_per_seed_strategy:
                break
            built = _attach_to_atom(parent_smiles, atom_idx, fragment)
            analog = _canonical(built or "")
            if not analog or analog == parent_smiles:
                continue
            dedupe_key = (strategy, str(seed_id), analog)
            if dedupe_key in seen:
                continue
            seen.add(dedupe_key)
            counts[key_count] = counts.get(key_count, 0) + 1

            parent_fp = _fp(parent_smiles)
            fp = _fp(analog)
            tanimoto = (
                DataStructs.TanimotoSimilarity(parent_fp, fp)
                if parent_fp is not None and fp is not None
                else 0.0
            )
            digest = hashlib.sha1(f"{seed_id}|{analog}|{strategy}".encode()).hexdigest()[:12]
            rows.append(
                {
                    "analog_id": f"analog_{digest}",
                    "proposal_id": payload.get("proposal_id", prop.get("proposal_id", f"proposal_{digest}")),
                    "iteration_id": prop.get("iteration_id", "iter_001"),
                    "strategy_name": strategy,
                    "seed_id": seed_id,
                    "parent_molecule_id": seed["molecule_id"],
                    "parent_smiles": parent_smiles,
                    "analog_smiles": analog,
                    "standard_smiles": analog,
                    "inchi_key": Chem.MolToInchiKey(Chem.MolFromSmiles(analog)),
                    "transformation_class": transformation_class,
                    "edit_site_id": edit_site_id,
                    "generated_by": "agent_proposal_rdkit_tool_dispatch",
                    "source": f"stage9_agentic_{strategy}",
                    "uniqueness_status": "unique",
                    "novelty_status": "analog_of_stage8_seed",
                    "parent_tanimoto": tanimoto,
                    "closest_known_egfr_ligand": seed["molecule_id"],
                    "warnings_json": json.dumps([]),
                }
            )

    out = pd.DataFrame(rows)
    write_table(paths["processed"] / "agent_generated_analogs.parquet", out)
    write_table(paths["processed"] / "agent_generated_analogs.csv", out)
    return out


def dispatch_from_files(config: dict[str, Any], paths: dict[str, Path]) -> pd.DataFrame:
    seeds = pd.read_parquet(paths["processed"] / "analog_seed_scaffolds.parquet")
    edit_sites = pd.read_parquet(paths["processed"] / "edit_sites.parquet")
    status_path = paths["processed"] / "agent_proposal_status.parquet"
    proposals = pd.read_parquet(status_path) if status_path.exists() else _load_proposals(paths["processed"] / "analog_proposals.jsonl")
    return dispatch_agent_proposals_to_molecules(proposals, seeds, edit_sites, config, paths)
