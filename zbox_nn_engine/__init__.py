"""
zbox_nn_engine — Topology-aware thermodynamics for DNA probe design.

A standalone Python package implementing three scoring layers for
DNA probe selectivity prediction:

1. S_ECI   — uniform-weight box-counting score
2. N_box   — analytical sub-box count
3. Z_box,NN — NN-weighted partition function

Based on SantaLucia 2004 unified nearest-neighbor parameters.

Usage
-----
    from zbox_nn_engine import (
        identify_retained_islands,
        compute_S_ECI,
        compute_N_box,
        compute_Zbox_NN,
        compute_all_scores,
        compute_box_dG,
    )

    # Single probe
    scores = compute_all_scores(
        probe_seq="ctcttgcctacgccaca",
        mt_target_seq="ctcttgcctacgccaca",
        wt_target_seq="ctcttgcctacgccacc",
        temperature_C=68,
        salt_mM=1000,
        threshold=2,
    )
    print(scores["Zbox_ratio"])
"""

from .nn_parameters import NN_PARAMS, INIT_GC, INIT_AT, R_GAS, get_nn_params
from .island_detection import identify_retained_islands, count_mismatches
from .scores import (
    compute_box_dG,
    compute_S_ECI,
    compute_N_box,
    compute_Zbox_NN,
    compute_all_scores,
)

__version__ = "1.0.0"
__author__ = "Thijs Van der Snickt et al."

__all__ = [
    # Parameters
    "NN_PARAMS",
    "INIT_GC",
    "INIT_AT",
    "R_GAS",
    "get_nn_params",
    # Island detection
    "identify_retained_islands",
    "count_mismatches",
    # Scoring
    "compute_box_dG",
    "compute_S_ECI",
    "compute_N_box",
    "compute_Zbox_NN",
    "compute_all_scores",
]
