"""Ligand contact projections from observed coordinates and DockLens paired results."""

from __future__ import annotations

import hashlib
import json
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from scipy.spatial import cKDTree  # type: ignore[import-untyped]

from structlens.core.models import AtomRecord, ResidueId
from structlens.core.parsing import ParsedStructure
from structlens.core.parsing.limits import ParseLimits
from structlens.core.reports.safe_io import read_bounded_bytes

PROXIMITY_ANGSTROM = 4.5
_PROFILES = {"plip", "dsv", "luna", "luna_dsv"}
_MAX_RESULT_BYTES = 128 * 1024 * 1024  # DockLens docklens-result-v1 limit.


@dataclass(frozen=True, slots=True)
class Contact:
    residue: str
    family: str
    receptor_atom: str
    ligand_atom: str
    distance: float | None
    evidence: str
    residue_id: ResidueId | None = None


@dataclass(frozen=True, slots=True)
class ProximityMap:
    contacts: tuple[Contact, ...]
    outside: int
    unavailable: int


@dataclass(frozen=True, slots=True)
class PoseContacts:
    pose_id: str
    ligand_id: str
    rows: tuple[Contact, ...]


@dataclass(frozen=True, slots=True)
class DockLensContacts:
    profile: str
    version: str
    run_id: str
    pose_hash: str
    counting_unit: str
    errors: int
    poses: tuple[PoseContacts, ...]


def residue_tag(residue: ResidueId, *, insertion: bool = True) -> str:
    code = (residue.insertion_code or "") if insertion else ""
    return f"{residue.residue_name}{residue.auth_seq_id}{code}{residue.chain_id}"


def _observed_heavy(atom: AtomRecord) -> bool:
    return atom.element.upper() not in {"H", "D", "T"} and atom.occupancy != 0


def ligand_proximity(parsed: ParsedStructure, component_id: str) -> ProximityMap:
    """One closest heavy-atom pair per residue of the selected model/chain."""
    ligand = next((c for c in parsed.components if c.component_id == component_id and c.is_ligand), None)
    if ligand is None:
        raise ValueError("Select a deposited ligand from the completed reference selection.")
    atoms = [a for a in ligand.atoms if _observed_heavy(a)]
    if not atoms:
        raise ValueError("Ligand coordinates unavailable: no observed heavy atoms.")
    tree = cKDTree([a.coordinate for a in atoms])
    rows: list[Contact] = []
    outside = unavailable = 0
    for chain in parsed.protein_structure.chains:
        for record in chain.residue_records:
            protein = [a for a in record.atoms if _observed_heavy(a)]
            if not protein:
                unavailable += 1
                continue
            distances, indices = tree.query([a.coordinate for a in protein])
            index = int(distances.argmin())
            d = float(distances[index])
            if d <= PROXIMITY_ANGSTROM:
                residue = record.residue_id
                rows.append(
                    Contact(
                        residue_tag(residue),
                        "Heavy-atom proximity",
                        protein[index].name,
                        atoms[int(indices[index])].name,
                        d,
                        "Nearest observed pair",
                        residue,
                    )
                )
            else:
                outside += 1
    if not outside and not rows:
        raise ValueError("Reference coordinates unavailable: no observed protein heavy atoms.")
    return ProximityMap(tuple(sorted(rows, key=lambda row: (row.distance or 0, row.residue))), outside, unavailable)


def _text(value: Any, field: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"Invalid DockLens {field}.")
    return value


def _array(doc: dict[str, Any], field: str) -> list[dict[str, Any]]:
    value = doc.get(field)
    if not isinstance(value, list) or any(not isinstance(item, dict) for item in value):
        raise ValueError(f"Invalid DockLens {field} list.")
    return value


def _object(doc: dict[str, Any], key: str) -> dict[str, Any]:
    value = doc.get(key)
    if not isinstance(value, dict):
        raise ValueError(f"Invalid DockLens {key} object.")
    return value


def _measurement(value: Any, label: str, *, angle: bool = False) -> float | None:
    if value is None:
        return None
    if type(value) not in (float, int):
        raise ValueError(f"Invalid DockLens {label}.")
    try:
        number = float(value)
    except OverflowError as exc:
        raise ValueError(f"Invalid DockLens {label}.") from exc
    if not math.isfinite(number) or number < 0 or (angle and number > 180):
        raise ValueError(f"Invalid DockLens {label}.")
    return number


def _contact(row: dict[str, Any]) -> Contact:
    details = [
        str(row[key])
        for key in ("subtype", "chemistry_basis", "chemistry_confidence", "hydrogen_atom")
        if isinstance(row.get(key), str) and row[key]
    ]
    for key, label, unit in (
        ("hydrogen_acceptor_distance_A", "H···A", "Å"),
        ("donor_hydrogen_acceptor_angle_deg", "D–H···A", "°"),
        ("hydrogen_acceptor_base_angle_deg", "H···A–base", "°"),
        ("theta_deg", "θ", "°"),
    ):
        value = _measurement(row.get(key), label, angle=unit == "°")
        if value is not None:
            details.append(f"{label} {value:g}{unit}")
    return Contact(
        _text(row.get("receptor_residue"), "receptor residue"),
        _text(row.get("type"), "interaction type"),
        _text(row.get("receptor_atom"), "receptor atom"),
        _text(row.get("ligand_atom"), "ligand atom"),
        _measurement(row.get("distance_A"), "distance"),
        " · ".join(details) or "Not reported",
    )


def import_docklens_contacts(path: Path, receptor_sha256: str, poses_path: Path) -> DockLensContacts:
    """Bind a DockLens result to exact receptor and pose bytes; read no paths from its JSON."""
    data = read_bounded_bytes(path, max_bytes=_MAX_RESULT_BYTES, label="DockLens result")
    try:
        doc = json.loads(data)
    except (ValueError, UnicodeError, RecursionError) as exc:
        raise ValueError("Invalid DockLens JSON result.") from exc
    if not isinstance(doc, dict) or doc.get("schema") != "docklens-result-v1":
        raise ValueError("Expected a DockLens docklens-result-v1 JSON result.")
    meta, inputs, params, totals = (_object(doc, key) for key in ("generated_by", "inputs", "parameters", "totals"))
    if meta.get("application") != "DockLens":
        raise ValueError("Invalid DockLens result metadata.")
    version = _text(meta.get("version"), "version")
    run_id = _text(doc.get("run_id"), "run ID")
    profile = _text(params.get("hbond_preset"), "profile")
    if profile not in _PROFILES:
        raise ValueError(f"Unknown DockLens interaction profile: {profile}")
    if inputs.get("receptor_sha256") != receptor_sha256:
        raise ValueError("DockLens receptor SHA-256 differs from the completed reference input.")
    poses_bytes = read_bounded_bytes(poses_path, max_bytes=ParseLimits().max_raw_bytes, label="DockLens poses")
    pose_hash = hashlib.sha256(poses_bytes).hexdigest()
    if inputs.get("poses_sha256") != pose_hash:
        raise ValueError("DockLens poses SHA-256 differs from the selected poses file.")

    poses = _array(doc, "poses")
    contacts = _array(doc, "interactions")
    if totals.get("poses") != len(poses) or totals.get("interactions") != len(contacts):
        raise ValueError("DockLens result totals disagree with the pose or interaction list.")
    records: dict[str, list[Contact]] = defaultdict(list)
    for item in contacts:
        records[_text(item.get("pose_id"), "pose ID")].append(_contact(item))
    parsed = []
    known = set()
    for item in poses:
        key = _text(item.get("pose_id"), "pose ID")
        if key in known or len(records[key]) != item.get("total_interactions"):
            raise ValueError("Duplicate pose ID or inconsistent per-pose total in DockLens result.")
        known.add(key)
        parsed.append(PoseContacts(key, _text(item.get("ligand_id"), "ligand ID"), tuple(records[key])))
    if set(records) - known:
        raise ValueError("DockLens contact names an unknown pose.")
    errors = totals.get("errors")
    if type(errors) is not int or errors < 0:
        raise ValueError("Invalid DockLens error count.")
    return DockLensContacts(
        profile, version, run_id, pose_hash, _text(params.get("counting_unit"), "counting unit"), errors, tuple(parsed)
    )
