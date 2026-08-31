"""Application orchestration for structural coordinate-quality evidence."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from importlib.metadata import version

from structlens.core.evidence import Availability, Diagnostic, DiagnosticSeverity
from structlens.core.models import ComponentKind, ResidueId
from structlens.core.parsing import ParsedStructure
from structlens.core.provenance import MethodProvenance
from structlens.core.quality.clashes import (
    CLASH_ALGORITHM_VERSION,
    VDW_RADII_VERSION,
    ClashAtom,
    ClashDiagnostic,
    ClashObservation,
    screen_heavy_atom_overlaps,
)
from structlens.core.quality.coordinates import check_coordinate_quality, diagnostics_for_atoms
from structlens.core.quality.geometry import screen_chain_geometry
from structlens.core.quality.models import (
    CoordinateAtom,
    CoordinateQCSettings,
    CoordinateResidue,
    StructureQualityReport,
)

_METHOD_ID = "structlens.coordinate_quality"
_METHOD_VERSION = "1"


@dataclass(frozen=True, slots=True)
class StructureQualityService:
    """Combine pure QC screens without placing scientific rules in the GUI."""

    settings: CoordinateQCSettings = CoordinateQCSettings()

    def __post_init__(self) -> None:
        if not isinstance(self.settings, CoordinateQCSettings):
            raise TypeError("settings must be CoordinateQCSettings")

    def analyze(self, parsed: ParsedStructure) -> StructureQualityReport:
        """Assess one selected, normalized structure deterministically."""

        if not isinstance(parsed, ParsedStructure):
            raise TypeError("parsed must be a ParsedStructure")

        source_id = parsed.selection.display_name
        coordinate_residues = tuple(
            CoordinateResidue(
                component.residue_id,
                tuple(CoordinateAtom.from_atom_record(atom) for atom in component.atoms),
                is_polymer=component.kind is ComponentKind.POLYMER_RESIDUE,
            )
            for component in parsed.components
            if component.residue_id is not None
        )
        atomic_report = check_coordinate_quality(
            coordinate_residues,
            source_id=source_id,
            settings=self.settings,
        )
        orphan_diagnostics = tuple(
            diagnostic
            for component in parsed.components
            if component.residue_id is None
            for diagnostic in diagnostics_for_atoms(
                component.atoms,
                source_id=source_id,
                settings=self.settings,
            )
        )
        geometry_diagnostics = screen_chain_geometry(
            parsed.protein_structure.chains,
            min_cn_distance_angstrom=self.settings.min_cn_distance_angstrom,
            max_cn_distance_angstrom=self.settings.max_cn_distance_angstrom,
            min_ca_distance_angstrom=self.settings.min_ca_distance_angstrom,
            max_ca_distance_angstrom=self.settings.max_ca_distance_angstrom,
        )
        clash_result = screen_heavy_atom_overlaps(
            _primary_polymer_atoms(parsed),
            overlap_tolerance_angstrom=self.settings.heavy_atom_overlap_tolerance_angstrom,
        )
        diagnostics = _merge_diagnostics(
            parsed.diagnostics,
            atomic_report.diagnostics,
            orphan_diagnostics,
            geometry_diagnostics,
            tuple(_clash_diagnostic(item, source_id) for item in clash_result.diagnostics),
            tuple(_overlap_diagnostic(item, source_id) for item in clash_result.observations),
        )

        analysis_atoms = sum(len(record.atoms) for chain in parsed.protein_structure.chains for record in chain.residue_records)
        errors = sum(item.severity is DiagnosticSeverity.ERROR for item in diagnostics)
        warnings = sum(item.severity is DiagnosticSeverity.WARNING for item in diagnostics)
        availability = (
            Availability.INVALID_INPUT
            if errors
            else Availability.AVAILABLE
            if analysis_atoms
            else Availability.NOT_APPLICABLE
        )
        component_counts = Counter(component.kind.value for component in parsed.components)
        counts = {
            "chains": len(parsed.protein_structure.chains),
            "polymer_residues": sum(
                len(chain.residue_records) for chain in parsed.protein_structure.chains
            ),
            "analysis_atoms": analysis_atoms,
            "source_components": len(parsed.components),
            "source_atoms": sum(len(component.atoms) for component in parsed.components),
            "screened_heavy_atoms": clash_result.screened_atom_count,
            "heavy_atom_overlaps": len(clash_result.observations),
            "errors": errors,
            "warnings": warnings,
            **{
                f"components_{kind.value}": component_counts.get(kind.value, 0)
                for kind in ComponentKind
            },
        }
        return StructureQualityReport(
            availability=availability,
            diagnostics=diagnostics,
            counts=counts,
            settings=self.settings,
            provenance=_provenance(parsed, self.settings),
        )


def assess_structure_quality(
    parsed: ParsedStructure,
    *,
    settings: CoordinateQCSettings | None = None,
) -> StructureQualityReport:
    """Functional entry point for the structural quality service."""

    return StructureQualityService(settings or CoordinateQCSettings()).analyze(parsed)


def _primary_polymer_atoms(parsed: ParsedStructure) -> tuple[ClashAtom, ...]:
    atoms: list[ClashAtom] = []
    for chain in parsed.protein_structure.chains:
        chain_id = f"{chain.structure_id}:{chain.model_id}:{chain.chain_id}"
        for residue_index, residue in enumerate(chain.residue_records):
            residue_id = _residue_label(residue.residue_id)
            for atom_index, atom in enumerate(residue.atoms):
                altloc = f":{atom.altloc}" if atom.altloc else ""
                atom_id = atom.source_atom_id or (
                    f"{residue_id}:{atom.name}{altloc}:{atom_index}"
                )
                atoms.append(
                    ClashAtom(
                        atom_id=atom_id,
                        element=atom.element,
                        coordinate=atom.coordinate,
                        atom_name=atom.name,
                        residue_id=residue_id,
                        component_id=residue_id,
                        chain_id=chain_id,
                        residue_index=residue_index,
                        residue_name=residue.residue_name,
                    )
                )
    return tuple(atoms)


def _residue_label(residue_id: ResidueId) -> str:
    insertion = residue_id.insertion_code or ""
    return ":".join(
        (
            residue_id.structure_id,
            residue_id.model_id,
            residue_id.chain_id,
            f"{residue_id.auth_seq_id}{insertion}",
            residue_id.residue_name,
        )
    )


def _clash_diagnostic(item: ClashDiagnostic, source_id: str) -> Diagnostic:
    return Diagnostic(
        code=f"quality.clashes.{item.code.casefold()}",
        severity=DiagnosticSeverity.WARNING,
        message=item.message,
        source_id=source_id,
        atom_id=item.atom_id,
        remediation="Review the source element and the declared van der Waals radius table.",
    )


def _overlap_diagnostic(item: ClashObservation, source_id: str) -> Diagnostic:
    return Diagnostic(
        code="quality.clashes.heavy_atom_overlap",
        severity=DiagnosticSeverity.WARNING,
        message=(
            f"Candidate heavy-atom overlap between {item.atom_a_id} and {item.atom_b_id}: "
            f"distance {item.distance_angstrom:.3f} Å; full van der Waals overlap "
            f"{item.overlap_angstrom:.3f} Å. This heavy-atom overlap screening does not add "
            "hydrogens and is not a crystallographic validation score."
        ),
        source_id=source_id,
        atom_id=f"{item.atom_a_id}|{item.atom_b_id}",
        residue_id=item.residue_a_id,
        remediation="Inspect the local coordinates, alternate conformer, and covalent context.",
    )


def _merge_diagnostics(*groups: tuple[Diagnostic, ...]) -> tuple[Diagnostic, ...]:
    unique: dict[tuple[object, ...], Diagnostic] = {}
    for diagnostic in (item for group in groups for item in group):
        key = (
            diagnostic.code,
            diagnostic.severity,
            diagnostic.message,
            diagnostic.source_id,
            diagnostic.atom_id,
            diagnostic.residue_id,
            diagnostic.remediation,
        )
        unique.setdefault(key, diagnostic)
    return tuple(
        sorted(
            unique.values(),
            key=lambda item: (
                item.code,
                item.source_id or "",
                item.residue_id or "",
                item.atom_id or "",
                item.message,
            ),
        )
    )


def _provenance(parsed: ParsedStructure, settings: CoordinateQCSettings) -> MethodProvenance:
    parameters = {
        **settings.to_json(),
        "selection_id": parsed.selection.selection_id,
        "model_id": parsed.selection.model_id,
        "author_chain_ids": parsed.selection.author_chain_ids,
        "label_chain_ids": parsed.selection.label_chain_ids,
        "altloc_policy": parsed.selection.altloc_policy.value,
        "assembly_scope": parsed.selection.assembly_scope.value,
        "atom_scope": "selected_primary_polymer_heavy_atoms",
        "component_evidence_scope": "all_selected_model_components",
        "connectivity_exclusion_policy": "same_residue_adjacent_residue_probable_disulfide",
        "source_connection_records": "not_retained_in_v0.4",
    }
    units = {
        "min_cn_distance_angstrom": "angstrom",
        "max_cn_distance_angstrom": "angstrom",
        "min_ca_distance_angstrom": "angstrom",
        "max_ca_distance_angstrom": "angstrom",
        "heavy_atom_overlap_tolerance_angstrom": "angstrom",
    }
    return MethodProvenance(
        method_id=_METHOD_ID,
        method_version=_METHOD_VERSION,
        parameters=parameters,
        units=units,
        backend_versions={
            "scipy": version("scipy"),
            "vdw_radii": VDW_RADII_VERSION,
            "overlap_screen": CLASH_ALGORITHM_VERSION,
        },
        input_hashes={
            "raw_source": parsed.raw_source_hash,
            "logical_content": parsed.selection.content_id,
        },
        analyzed_representation=parsed.selection.assembly_scope.value,
    )


__all__ = ["StructureQualityService", "assess_structure_quality"]
