"""Observed heavy-atom contacts between two protein chains of one deposited model."""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from collections.abc import Sequence
from dataclasses import dataclass, field
from threading import Event

from Bio.Data.PDBData import protein_letters_3to1_extended
from scipy.spatial import cKDTree  # type: ignore[import-untyped]

from structlens.core.errors import AnalysisCancelledError
from structlens.core.interactions import InteractionRecord, InteractionThresholds, InteractionType
from structlens.core.interactions.chemistry import atom_is_acceptor, atom_is_donor, residue_is_hydrophobic
from structlens.core.models import AtomRecord, ProteinChain, ResidueId, ResidueRecord
from structlens.core.parsing import InputSelection, ParsedStructure, SourceSnapshot, load_structure_evidence

_THRESHOLDS = InteractionThresholds()
_ACID = {"ASP": {"OD1", "OD2"}, "GLU": {"OE1", "OE2"}}
_BASE = {"LYS": {"NZ"}, "ARG": {"NE", "NH1", "NH2"}}
_BACKBONE = {"N", "CA", "C", "O", "OXT"}
_AMINO_ACIDS = frozenset(protein_letters_3to1_extended) | {"SEC", "PYL"}


@dataclass(frozen=True, slots=True)
class InterfaceContact:
    residue_a: ResidueId
    residue_b: ResidueId
    atom_a: str
    atom_b: str
    distance_angstrom: float
    atom_pair_count: int
    interactions: tuple[InteractionRecord, ...]


@dataclass(frozen=True, slots=True)
class InterfaceResidue:
    residue_id: ResidueId
    partner_count: int
    atom_pair_count: int
    min_distance_angstrom: float


@dataclass(frozen=True, slots=True)
class SubunitInterface:
    source_sha256: str
    model_id: str
    chain_a: str
    chain_b: str
    cutoff_angstrom: float
    contacts: tuple[InterfaceContact, ...]
    residues: tuple[InterfaceResidue, ...]
    observed_a: int
    observed_b: int
    unavailable_a: int
    unavailable_b: int
    excluded_a: int
    excluded_b: int


def _one_chain(parsed: ParsedStructure, chain_id: str) -> ProteinChain:
    matches = [chain for chain in parsed.protein_structure.chains if chain.chain_id == chain_id]
    if len(matches) != 1:
        raise ValueError(f"Author chain {chain_id!r} is absent or ambiguous in this model.")
    return matches[0]


def is_protein_chain(chain: ProteinChain) -> bool:
    return any(residue.residue_name in _AMINO_ACIDS for residue in chain.residue_records)


def eligible_interface_chains(chains: Sequence[ProteinChain], model_id: str) -> tuple[str, ...]:
    """Author IDs of unambiguous protein chains within a single model."""
    selected = [chain for chain in chains if chain.model_id == model_id]
    counts = Counter(chain.chain_id for chain in selected)
    return tuple(chain.chain_id for chain in selected if counts[chain.chain_id] == 1 and is_protein_chain(chain))


def _check_cancel(cancel_event: Event | None) -> None:
    if cancel_event is not None and cancel_event.is_set():
        raise AnalysisCancelledError("Interface analysis cancelled.")


def _eligible(
    chain: ProteinChain, cancel_event: Event | None
) -> tuple[list[tuple[ResidueRecord, AtomRecord]], int, int, int]:
    atoms: list[tuple[ResidueRecord, AtomRecord]] = []
    observed = unavailable = excluded = 0
    for residue in chain.residue_records:
        _check_cancel(cancel_event)
        if residue.residue_name not in _AMINO_ACIDS:
            excluded += 1
            continue
        eligible = [
            (residue, atom)
            for atom in residue.atoms
            if atom.element not in {"H", "D", "T"} and (atom.occupancy is None or atom.occupancy > 0)
        ]
        if eligible:
            observed += 1
            atoms.extend(eligible)
        else:
            unavailable += 1
    if not atoms:
        raise ValueError(f"Chain {chain.chain_id!r} has no observed protein heavy atoms in this model.")
    return atoms, observed, unavailable, excluded


def _charge(residue: ResidueRecord, atom: AtomRecord) -> int:
    if atom.element == "O" and atom.name in _ACID.get(residue.residue_name, ()):
        return -1
    if atom.element == "N" and atom.name in _BASE.get(residue.residue_name, ()):
        return 1
    return 0


def _donor(residue: ResidueRecord, atom: AtomRecord) -> bool:
    return atom_is_donor(atom) and (atom.name != "N" or residue.residue_name != "PRO")


def _acceptor(residue: ResidueRecord, atom: AtomRecord) -> bool:
    if atom.element == "N":
        return residue.residue_name == "HIS" and atom.name in {"ND1", "NE2"}
    return atom.element == "O" and (atom.name == "OXT" or atom_is_acceptor(atom))


def _nonpolar(residue: ResidueRecord, atom: AtomRecord) -> bool:
    return atom.element in {"C", "S"} and atom.name not in _BACKBONE and residue_is_hydrophobic(residue)


def _types(a: ResidueRecord, x: AtomRecord, b: ResidueRecord, y: AtomRecord, distance: float) -> list[InteractionType]:
    kinds = []
    if distance <= _THRESHOLDS.salt_bridge_distance_angstrom and _charge(a, x) * _charge(b, y) == -1:
        kinds.append(InteractionType.SALT_BRIDGE)
    if distance <= _THRESHOLDS.hbond_distance_angstrom and (
        (_donor(a, x) and _acceptor(b, y)) or (_donor(b, y) and _acceptor(a, x))
    ):
        kinds.append(InteractionType.HBOND_GEOMETRIC)
    if distance <= _THRESHOLDS.hydrophobic_distance_angstrom and _nonpolar(a, x) and _nonpolar(b, y):
        kinds.append(InteractionType.HYDROPHOBIC)
    return kinds


@dataclass
class _Pair:
    nearest: tuple[float, str, str] = (math.inf, "", "")
    count: int = 0
    typed: dict[InteractionType, tuple[float, str, str]] = field(default_factory=dict)


def _contacts(
    first: list[tuple[ResidueRecord, AtomRecord]],
    second: list[tuple[ResidueRecord, AtomRecord]],
    cutoff: float,
    cancel_event: Event | None,
) -> tuple[InterfaceContact, ...]:
    tree = cKDTree([atom.coordinate for _, atom in second])
    total = 0
    pairs: dict[tuple[ResidueId, ResidueId], _Pair] = {}
    for a, x in first:
        _check_cancel(cancel_event)
        total += int(tree.query_ball_point(x.coordinate, cutoff, return_length=True))
        if total > 2_000_000:
            raise ValueError("Interface exceeds two million candidate atom pairs; choose a smaller cutoff.")
        for n, j in enumerate(tree.query_ball_point(x.coordinate, cutoff)):
            if n % 1024 == 0:
                _check_cancel(cancel_event)
            b, y = second[j]
            d = math.dist(x.coordinate, y.coordinate)
            if d > cutoff:
                continue
            candidate = (d, x.name, y.name)
            key = a.residue_id, b.residue_id
            pair = pairs.get(key)
            if pair is None:
                if len(pairs) >= 100_000:
                    raise ValueError("Interface exceeds 100,000 residue pairs; choose a smaller cutoff.")
                pairs[key] = pair = _Pair()
            pair.count += 1
            pair.nearest = min(pair.nearest, candidate)
            for kind in _types(a, x, b, y, d):
                pair.typed[kind] = min(pair.typed.get(kind, candidate), candidate)

    rows = []
    for (ra, rb), pair in pairs.items():
        distance, atom_a, atom_b = pair.nearest
        interactions = tuple(
            InteractionRecord(ra.structure_id, kind, ra, rb, xname, yname, d, evidence_mode="heavy_atom_geometry")
            for kind, (d, xname, yname) in sorted(pair.typed.items())
        )
        rows.append(InterfaceContact(ra, rb, atom_a, atom_b, distance, pair.count, interactions))
    return tuple(
        sorted(rows, key=lambda row: (row.distance_angstrom, _identifier(row.residue_a), _identifier(row.residue_b)))
    )


def _identifier(residue: ResidueId) -> tuple[str, str, str, str]:
    return residue.chain_id, residue.auth_seq_id, residue.insertion_code or "", residue.residue_name


def _summaries(contacts: tuple[InterfaceContact, ...]) -> tuple[InterfaceResidue, ...]:
    by_residue: dict[ResidueId, list[InterfaceContact]] = defaultdict(list)
    for contact in contacts:
        by_residue[contact.residue_a].append(contact)
        by_residue[contact.residue_b].append(contact)
    rows = (
        InterfaceResidue(
            residue, len(partners), sum(p.atom_pair_count for p in partners), min(p.distance_angstrom for p in partners)
        )
        for residue, partners in by_residue.items()
    )
    return tuple(sorted(rows, key=lambda r: (-r.partner_count, r.min_distance_angstrom, _identifier(r.residue_id))))


def analyze_subunit_interface(
    snapshot: SourceSnapshot,
    model_id: str,
    chain_a: str,
    chain_b: str,
    cutoff_angstrom: float = 5.0,
    *,
    cancel_event: Event | None = None,
) -> SubunitInterface:
    """Describe deposited contacts for the chosen author chains, without applying a fit or assembly operators."""
    if not isinstance(snapshot, SourceSnapshot):
        raise TypeError("Interface analysis requires a captured source snapshot.")
    _check_cancel(cancel_event)
    if not isinstance(model_id, str) or not model_id.strip():
        raise ValueError("Choose one deposited model.")
    if not isinstance(chain_a, str) or not isinstance(chain_b, str) or chain_a.strip() == chain_b.strip():
        raise ValueError("Choose two distinct author chains.")
    try:
        cutoff = float(cutoff_angstrom)
    except (OverflowError, TypeError, ValueError) as exc:
        raise ValueError("Interface cutoff must be a finite positive distance.") from exc
    if isinstance(cutoff_angstrom, bool) or not math.isfinite(cutoff) or cutoff <= 0:
        raise ValueError("Interface cutoff must be a finite positive distance.")
    selection = InputSelection(
        snapshot.content_id,
        snapshot.display_name,
        snapshot.logical_format,
        model_id,
        author_chain_ids=(chain_a.strip(), chain_b.strip()),
    )
    parsed = load_structure_evidence(snapshot, selection=selection)
    first, observed_a, unavailable_a, excluded_a = _eligible(_one_chain(parsed, chain_a.strip()), cancel_event)
    second, observed_b, unavailable_b, excluded_b = _eligible(_one_chain(parsed, chain_b.strip()), cancel_event)
    contacts = _contacts(first, second, cutoff, cancel_event)
    residues = _summaries(contacts)
    _check_cancel(cancel_event)
    return SubunitInterface(
        snapshot.raw_sha256,
        selection.model_id,
        chain_a.strip(),
        chain_b.strip(),
        cutoff,
        contacts,
        residues,
        observed_a,
        observed_b,
        unavailable_a,
        unavailable_b,
        excluded_a,
        excluded_b,
    )
