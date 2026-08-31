"""Active-site definition and descriptive geometry application service."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np

from structlens.core.geometry.kabsch import apply_transform, kabsch
from structlens.core.models import AtomRecord, ResidueId, ResidueRecord
from structlens.core.sites import SiteDefinition, SiteDefinitionMode, SiteMetrics

# A backbone RMSD is measured over the peptide backbone, not the C-alpha trace.
BACKBONE_ATOMS = ("N", "CA", "C", "O")


def _coords(residue: ResidueRecord, names: set[str] | None = None) -> np.ndarray:
    atoms = [atom.coordinate for atom in residue.atoms if names is None or atom.name.upper() in names]
    return np.asarray(atoms, dtype=np.float64)


def _backbone_coords(residue: ResidueRecord) -> np.ndarray | None:
    """Backbone coordinates in a fixed atom order, or None if any are absent.

    Returning None for an incomplete backbone keeps a partially resolved residue
    out of the RMSD rather than silently pairing mismatched atoms.
    """

    by_name = {atom.name.upper(): atom.coordinate for atom in residue.atoms}
    if any(name not in by_name for name in BACKBONE_ATOMS):
        return None
    return np.asarray([by_name[name] for name in BACKBONE_ATOMS], dtype=np.float64)


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
        if (
            len(coordinates)
            and np.min(np.linalg.norm(coordinates[:, None, :] - center_atoms[None, :, :], axis=2)) <= radius
        ):
            selected.append(residue)
    return tuple(selected)


def _centroid(values: np.ndarray) -> np.ndarray | None:
    return np.mean(values, axis=0) if len(values) else None


def _rmsd(reference: np.ndarray, target: np.ndarray) -> float | None:
    if not len(reference) or reference.shape != target.shape:
        return None
    return float(np.sqrt(np.mean(np.sum((target - reference) ** 2, axis=1))))


def _site_fitted_rmsd(reference: np.ndarray, target: np.ndarray) -> float | None:
    """RMSD after a proper Kabsch fit of the site onto its reference.

    Subtracting centroids removes translation only. A site that is rotated
    relative to its reference would still report a large residual, which is not
    what a site-fitted measurement means.
    """

    if not len(reference) or reference.shape != target.shape:
        return None
    try:
        rotation, translation = kabsch(reference, target)
    except ValueError:
        return None
    return _rmsd(reference, apply_transform(target, rotation, translation))


def _envelope_volume(values: np.ndarray) -> float | None:
    """Convex-hull envelope volume, or None when no 3D envelope exists.

    Fewer than four atoms, and atoms that are collinear or coplanar, enclose no
    measurable volume. Those cases are unavailable rather than zero, and a
    degenerate hull raises QhullError rather than returning a value.
    """

    if len(values) < 4:
        return None
    try:
        from scipy.spatial import ConvexHull, QhullError  # type: ignore[import-untyped]
    except ImportError:
        return None
    try:
        return float(ConvexHull(values).volume)
    except (QhullError, ValueError):
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
    mapped_pairs = [
        (item, target_by_id[correspondence[item.residue_id]])
        for item in selected
        if item.residue_id in correspondence and correspondence[item.residue_id] in target_by_id
    ]
    ref_ca: list[tuple[float, float, float]] = []
    tar_ca: list[tuple[float, float, float]] = []
    tar_atoms: list[tuple[float, float, float]] = []
    polar = charged = 0
    for reference, target in mapped_pairs:
        reference_backbone = _backbone_coords(reference)
        target_backbone = _backbone_coords(target)
        if reference_backbone is not None and target_backbone is not None:
            ref_ca.extend(tuple(row) for row in reference_backbone)
            tar_ca.extend(tuple(row) for row in target_backbone)
        tar_atoms.extend(_coords(target).tolist())
        # Composition belongs to the structure represented by this metric.
        # The target may have a chemically different residue after a mapped
        # substitution, so never inherit the reference residue name here.
        polar += int(target.residue_name.upper() in {"SER", "THR", "ASN", "GLN", "TYR", "HIS"})
        charged += int(target.residue_name.upper() in {"ARG", "LYS", "ASP", "GLU", "HIS"})
    ref_array = np.asarray(ref_ca, dtype=np.float64)
    tar_array = np.asarray(tar_ca, dtype=np.float64)
    if target_transform is not None and len(tar_array):
        rotation = np.asarray(target_transform, dtype=np.float64)
        if rotation.shape == (4, 4):
            homogeneous = np.c_[tar_array, np.ones(len(tar_array))]
            tar_array = (homogeneous @ rotation.T)[:, :3]
    ref_centroid = _centroid(ref_array)
    tar_centroid = _centroid(tar_array)
    centroid_displacement = (
        float(np.linalg.norm(tar_centroid - ref_centroid))
        if ref_centroid is not None and tar_centroid is not None
        else None
    )
    rg = float(np.sqrt(np.mean(np.sum((ref_array - ref_centroid) ** 2, axis=1)))) if ref_centroid is not None else None
    # A global-frame RMSD only means something once the target has been placed
    # in the reference frame. Without an authoritative transform the two
    # structures sit in unrelated frames, so the measurement is unavailable
    # rather than an arbitrarily large number.
    global_rmsd = _rmsd(ref_array, tar_array) if target_transform is not None else None
    site_fit_rmsd = _site_fitted_rmsd(ref_array, tar_array)
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
