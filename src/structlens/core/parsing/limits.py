"""Resource limits applied at the structure parsing boundary."""

from __future__ import annotations

from dataclasses import dataclass


class SnapshotError(ValueError):
    """A source could not be captured as a valid supported snapshot."""


class SnapshotLimitError(SnapshotError):
    """A source or parser count exceeded an explicit input limit."""

    code = "input_limit_exceeded"

    def __init__(self, limit_name: str, limit: int, observed: int) -> None:
        self.limit_name = limit_name
        self.limit = limit
        self.observed = observed
        super().__init__(
            f"input limit exceeded for {limit_name}: observed {observed}, limit {limit}"
        )

    @property
    def diagnostic(self) -> str:
        """Stable human-readable diagnostic suitable for a parser report."""

        return f"{self.code}: {self}"


class StructureParseError(SnapshotError):
    """A bounded PDB/mmCIF source is malformed at the parser boundary."""

    code = "structure_parse_error"


@dataclass(frozen=True, slots=True)
class ParseLimits:
    """Hard caps for source capture and parser-facing record counts.

    Byte limits are enforced while reading the source.  The record limits are
    exposed here so parser adapters can validate counts before constructing
    large object graphs.
    """

    max_raw_bytes: int = 64 * 1024 * 1024
    max_decompressed_bytes: int = 256 * 1024 * 1024
    max_models: int = 64
    max_chains: int = 1_024
    max_residues: int = 1_000_000
    max_atoms: int = 10_000_000

    def __post_init__(self) -> None:
        for name in (
            "max_raw_bytes",
            "max_decompressed_bytes",
            "max_models",
            "max_chains",
            "max_residues",
            "max_atoms",
        ):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
                raise ValueError(f"{name} must be a positive integer")

    def validate_counts(
        self,
        *,
        models: int | None = None,
        chains: int | None = None,
        residues: int | None = None,
        atoms: int | None = None,
    ) -> None:
        """Reject parser counts that exceed the configured hard caps."""

        for count_name, observed, limit_name in (
            ("models", models, "max_models"),
            ("chains", chains, "max_chains"),
            ("residues", residues, "max_residues"),
            ("atoms", atoms, "max_atoms"),
        ):
            if observed is None:
                continue
            if isinstance(observed, bool) or not isinstance(observed, int) or observed < 0:
                raise ValueError(f"{count_name} must be a non-negative integer")
            limit = getattr(self, limit_name)
            if observed > limit:
                raise SnapshotLimitError(count_name, limit, observed)


__all__ = ["ParseLimits", "SnapshotError", "SnapshotLimitError", "StructureParseError"]
