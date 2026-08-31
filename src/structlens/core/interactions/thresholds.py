"""Centralized geometric interaction thresholds in Angstroms."""

import math
from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class InteractionThresholds:
    hbond_distance_angstrom: float = 3.5
    salt_bridge_distance_angstrom: float = 4.0
    hydrophobic_distance_angstrom: float = 4.5
    pi_centroid_distance_angstrom: float = 5.5
    cation_pi_distance_angstrom: float = 6.0
    metal_distance_angstrom: float = 3.0
    pi_parallel_angle_tolerance_degrees: float = 30.0
    pi_t_shape_min_angle_degrees: float = 60.0
    pi_t_shape_max_angle_degrees: float = 120.0
    cation_pi_normal_tolerance_degrees: float = 45.0

    def __post_init__(self) -> None:
        distances = (
            "hbond_distance_angstrom",
            "salt_bridge_distance_angstrom",
            "hydrophobic_distance_angstrom",
            "pi_centroid_distance_angstrom",
            "cation_pi_distance_angstrom",
            "metal_distance_angstrom",
        )
        angles = (
            "pi_parallel_angle_tolerance_degrees",
            "pi_t_shape_min_angle_degrees",
            "pi_t_shape_max_angle_degrees",
            "cation_pi_normal_tolerance_degrees",
        )
        for name in distances:
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be finite and positive")
            object.__setattr__(self, name, value)
        for name in angles:
            value = float(getattr(self, name))
            if not math.isfinite(value) or not 0.0 <= value <= 180.0:
                raise ValueError(f"{name} must be finite and between 0 and 180 degrees")
            object.__setattr__(self, name, value)
        if self.pi_t_shape_min_angle_degrees > self.pi_t_shape_max_angle_degrees:
            raise ValueError("pi_t_shape_min_angle_degrees must not exceed pi_t_shape_max_angle_degrees")


__all__ = ["InteractionThresholds"]
