"""Model states and polymer conformers retained by the canonical parser."""

from __future__ import annotations

from typing import Any

from structlens.core.parsing import InputSelection, SourceSnapshot, load_structure_evidence


def selected_coordinates(snapshot: SourceSnapshot, selection: InputSelection) -> tuple[int, list[dict[str, Any]]]:
    if len(selection.author_chain_ids) != 1:
        raise ValueError("PyMOL requires a single selected author chain per structure")
    parsed = load_structure_evidence(snapshot, selection=selection)
    chain_id = selection.author_chain_ids[0]
    chains = [chain for chain in parsed.protein_structure.chains if chain.chain_id == chain_id]
    if not chains:
        raise ValueError(f"Selected model {selection.model_id} / chain {chain_id} has no analyzed residues.")
    disordered = {component.residue_id for component in parsed.components if any(atom.altloc for atom in component.atoms)}
    conformers = []
    for chain in chains:
        for residue in chain.residue_records:
            rid = residue.residue_id
            if rid in disordered:
                conformers.append({
                    "resi": f"{rid.auth_seq_id}{rid.insertion_code or ''}", "resn": rid.residue_name,
                    "atoms": [[atom.name, atom.altloc or ""] for atom in residue.atoms],
                })
    return parsed.metadata.available_models.index(selection.model_id) + 1, conformers
