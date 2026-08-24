"""Centralized geometric interaction thresholds in Angstroms."""

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


__all__ = ["InteractionThresholds"]
