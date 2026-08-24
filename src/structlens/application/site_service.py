"""Active-site definition and descriptive geometry application service."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np

from structlens.core.models import AtomRecord, ResidueId, ResidueRecord
from structlens.core.sites import SiteDefinition, SiteDefinitionMode, SiteMetrics


def _coords(residue: ResidueRecord, names: set[str] | None = None) -> np.ndarray:
    atoms = [atom.coordinate for atom in residue.atoms if names is None or atom.name.upper() in names]
    return np.asarray(atoms, dtype=np.float64)


def define_site(
    definition: SiteDefinition,
    reference_residues: Sequence[ResidueRecord],
    *,
    ligand_atoms: Mapping[str, Sequence[AtomRecord]] | None = None,
) -> tuple[ResidueRecord, ...]:
    """Resolve a site from reference structure records without renumbering."""

    by_id = {residue.residue_id: residue for residue in reference_residues}
    if definition.mode is SiteDefinitionMode.KEY_RESIDUES:
        return tuple(by_id[item] for item in definition.reference_residues if item in by_id)
    if definition.mode is SiteDefinitionMode.RESIDUE_RADIUS:
        if definition.center_residue is None or definition.radius_angstrom is None:
            return ()
        center = by_id.get(definition.center_residue)
        if center is None:
            return ()
        center_atoms = _coords(center)
    else:
        atoms = tuple((ligand_atoms or {}).get(definition.ligand_id or "", ()))
        if not atoms or definition.radius_angstrom is None:
            return ()
        center_atoms = np.asarray([atom.coordinate for atom in atoms], dtype=np.float64)
    radius = float(definition.radius_angstrom or 0.0)
    selected: list[ResidueRecord] = []
    for residue in reference_residues:
        coordinates = _coords(residue)
        if len(coordinates) and np.min(np.linalg.norm(coordinates[:, None, :] - center_atoms[None, :, :], axis=2)) <= radius:
            selected.append(residue)
    return tuple(selected)


def _centroid(values: np.ndarray) -> np.ndarray | None:
    return np.mean(values, axis=0) if len(values) else None


def _rmsd(reference: np.ndarray, target: np.ndarray) -> float | None:
    if not len(reference) or reference.shape != target.shape:
        return None
    return float(np.sqrt(np.mean(np.sum((target - reference) ** 2, axis=1))))


def _envelope_volume(values: np.ndarray) -> float | None:
    if len(values) < 4:
        return 0.0 if len(values) else None
    try:
        from scipy.spatial import ConvexHull  # type: ignore[import-untyped]

        return float(ConvexHull(values).volume)
    except (ImportError, ValueError):
        return None


def calculate_site_metrics(
    definition: SiteDefinition,
    reference_residues: Sequence[ResidueRecord],
    target_residues: Sequence[ResidueRecord],
    correspondence: Mapping[ResidueId, ResidueId],
    *,
    target_structure_id: str,
    target_transform: np.ndarray | None = None,
    sasa_angstrom2: float | None = None,
    ligand_atoms: Mapping[str, Sequence[AtomRecord]] | None = None,
) -> SiteMetrics:
    """Calculate site metrics while preserving deleted/unmapped residues."""

    target_by_id = {item.residue_id: item for item in target_residues}
    selected = define_site(definition, reference_residues, ligand_atoms=ligand_atoms)
    mapped_pairs = [(item, target_by_id[correspondence[item.residue_id]]) for item in selected if item.residue_id in correspondence and correspondence[item.residue_id] in target_by_id]
    ref_ca: list[tuple[float, float, float]] = []
    tar_ca: list[tuple[float, float, float]] = []
    tar_atoms: list[tuple[float, float, float]] = []
    polar = charged = 0
    for reference, target in mapped_pairs:
        ref_coordinates = _coords(reference, {"CA"})
        target_coordinates = _coords(target, {"CA"})
        if len(ref_coordinates) and len(target_coordinates):
            ref_ca.append(tuple(ref_coordinates[0]))
            tar_ca.append(tuple(target_coordinates[0]))
        tar_atoms.extend(_coords(target).tolist())
        polar += int(reference.residue_name.upper() in {"SER", "THR", "ASN", "GLN", "TYR", "HIS"})
        charged += int(reference.residue_name.upper() in {"ARG", "LYS", "ASP", "GLU", "HIS"})
    ref_array = np.asarray(ref_ca, dtype=np.float64)
    tar_array = np.asarray(tar_ca, dtype=np.float64)
    if target_transform is not None and len(tar_array):
        rotation = np.asarray(target_transform, dtype=np.float64)
        if rotation.shape == (4, 4):
            homogeneous = np.c_[tar_array, np.ones(len(tar_array))]
            tar_array = (homogeneous @ rotation.T)[:, :3]
    ref_centroid = _centroid(ref_array)
    tar_centroid = _centroid(tar_array)
    centroid_displacement = float(np.linalg.norm(tar_centroid - ref_centroid)) if ref_centroid is not None and tar_centroid is not None else None
    rg = float(np.sqrt(np.mean(np.sum((ref_array - ref_centroid) ** 2, axis=1)))) if ref_centroid is not None else None
    global_rmsd = _rmsd(ref_array, tar_array)
    site_fit_rmsd = global_rmsd
    if ref_centroid is not None and tar_centroid is not None and len(ref_array):
        site_fit_rmsd = _rmsd(ref_array - ref_centroid, tar_array - tar_centroid)
    return SiteMetrics(
        definition.site_id,
        target_structure_id,
        len(mapped_pairs),
        len(mapped_pairs) / len(selected) if selected else 0.0,
        global_rmsd,
        site_fit_rmsd,
        centroid_displacement,
        rg,
        # Metrics describe the selected structure (the target here), so the
        # envelope is built from target heavy atoms; unavailable target atoms
        # remain unavailable rather than inheriting reference geometry.
        _envelope_volume(np.asarray(tar_atoms, dtype=np.float64)),
        sasa_angstrom2,
        polar / len(mapped_pairs) if mapped_pairs else None,
        charged / len(mapped_pairs) if mapped_pairs else None,
    )


__all__ = ["calculate_site_metrics", "define_site"]
