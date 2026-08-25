"""
Topology-aware thermodynamic scores for DNA probe design.

Three scoring layers are provided:

1. **S_ECI** — uniform-weight box-counting score:
       S_ECI = sum_i c_i^2
   where c_i is the retained-core length of island i.

2. **N_box** — analytical sub-box count:
       N_box(L, k) = (L - k + 1)(L - k + 2) / 2
   giving the number of contiguous sub-boxes of length >= k in an island
   of length L.

3. **Z_box,NN** — NN-weighted partition function:
       Z_box,NN = sum_b exp[-dG_NN(b, T) / RT]
   where the sum runs over all retained islands (boxes) and dG_NN is the
   nearest-neighbor free energy of the perfect duplex segment at the
   specified temperature and salt concentration.
"""

import numpy as np

from .nn_parameters import R_GAS, get_nn_params, INIT_GC, INIT_AT


# ---------------------------------------------------------------------------
# NN free energy of a single box (perfect duplex segment)
# ---------------------------------------------------------------------------

def compute_box_dG(box_seq, temperature_C=37, salt_mM=1000):
    """
    Compute the NN free energy of a contiguous matched box
    (perfect duplex segment) at the specified temperature and salt.

    Uses SantaLucia 2004 unified NN parameters with:
    - Terminal initiation corrections (AT or GC ends)
    - Temperature correction:  dG(T) = dH - T * dS
    - Salt correction (SantaLucia 1998 entropy form):
          dS_salt = 0.368 * (N-1) * ln([Na+])

    Parameters
    ----------
    box_seq : str
        Top-strand sequence of the box (5'->3').
    temperature_C : float
        Temperature in degrees Celsius.
    salt_mM : float
        Na+ concentration in mM.

    Returns
    -------
    float
        Free energy dG in kcal/mol.
    """
    box_seq = box_seq.upper()

    if len(box_seq) < 2:
        # Single base pair: only initiation
        base = box_seq[0]
        init = INIT_GC if base in 'GC' else INIT_AT
        dG37 = init['dG37']
        dH = init['dH']
        dS = init['dS']
    else:
        dG37 = 0.0
        dH = 0.0
        dS = 0.0

        # Terminal initiation at 5' end
        init5 = INIT_GC if box_seq[0] in 'GC' else INIT_AT
        dG37 += init5['dG37']
        dH += init5['dH']
        dS += init5['dS']

        # Terminal initiation at 3' end
        init3 = INIT_GC if box_seq[-1] in 'GC' else INIT_AT
        dG37 += init3['dG37']
        dH += init3['dH']
        dS += init3['dS']

        # NN steps
        for i in range(len(box_seq) - 1):
            p = get_nn_params(box_seq[i:i + 2])
            dG37 += p['dG37']
            dH += p['dH']
            dS += p['dS']

    # Temperature correction: dG(T) = dH - T * dS  (dS in cal, dH in kcal)
    T_K = temperature_C + 273.15
    dG_T = dH - T_K * dS / 1000.0

    # Salt correction (SantaLucia 1998)
    N_phosphates = max(len(box_seq) - 1, 1)
    Na_M = salt_mM / 1000.0
    if Na_M > 0:
        dS_salt = 0.368 * N_phosphates * np.log(Na_M)
        dG_T += -T_K * dS_salt / 1000.0

    return dG_T


# ---------------------------------------------------------------------------
# Scoring functions
# ---------------------------------------------------------------------------

def compute_S_ECI(retained_islands):
    """
    Uniform-weight box-counting score.

    S_ECI = sum_i c_i^2

    where c_i is the retained-core length of each retained island.

    Parameters
    ----------
    retained_islands : list of dict
        Output from ``island_detection.identify_retained_islands``.

    Returns
    -------
    float
    """
    return sum(island['retained_core'] ** 2 for island in retained_islands)


def compute_N_box(L, k):
    """
    Analytical sub-box count for an island of length L with threshold k.

    N_box(L, k) = (L - k + 1)(L - k + 2) / 2

    This counts all contiguous sub-boxes of length >= k within an island
    of length L.

    Parameters
    ----------
    L : int
        Island length (raw, before edge correction).
    k : int
        Threshold (minimum sub-box length).

    Returns
    -------
    int
    """
    if L < k:
        return 0
    return (L - k + 1) * (L - k + 2) // 2


def compute_Zbox_NN(retained_islands, temperature_C=68, salt_mM=1000,
                    threshold=2, RT=None):
    """
    NN-weighted partition function over retained islands.

    Z_box,NN = sum_b exp[-dG_NN(b, T) / RT]

    Each retained island contributes one term.  The NN free energy is
    computed for the perfect duplex segment (box_seq) at the specified
    temperature and salt.

    Parameters
    ----------
    retained_islands : list of dict
        Output from ``island_detection.identify_retained_islands``.
    temperature_C : float
        Temperature in degrees Celsius.
    salt_mM : float
        Na+ concentration in mM.
    threshold : int
        Minimum retained-core length (passed for documentation; the islands
        are already filtered).
    RT : float or None
        Pre-computed RT in kcal/mol.  If None, computed from temperature.

    Returns
    -------
    float
    """
    if RT is None:
        T_K = temperature_C + 273.15
        RT = R_GAS * T_K / 1000.0   # kcal/mol

    Z = 0.0
    for island in retained_islands:
        dG = compute_box_dG(island['box_seq'], temperature_C, salt_mM)
        Z += np.exp(-dG / RT)

    return Z


# ---------------------------------------------------------------------------
# Convenience: compute all scores for one probe against MT and WT targets
# ---------------------------------------------------------------------------

def compute_all_scores(probe_seq, mt_target_seq, wt_target_seq,
                       temperature_C=68, salt_mM=1000, threshold=2):
    """
    Compute S_ECI and Z_box,NN for a probe against both MT and WT targets.

    Parameters
    ----------
    probe_seq : str
        Probe sequence (5'->3').
    mt_target_seq : str
        Mutant-type target sequence (5'->3'), same length as probe.
    wt_target_seq : str
        Wild-type target sequence (5'->3'), same length as probe.
    temperature_C : float
        Hybridization temperature in °C.
    salt_mM : float
        Na+ concentration in mM.
    threshold : int
        Minimum retained-core length for island retention.

    Returns
    -------
    dict with keys:
        S_ECI_MT, S_ECI_WT, S_ECI_ratio,
        Zbox_MT, Zbox_WT, Zbox_ratio,
        n_islands_MT, n_islands_WT,
        islands_MT, islands_WT
    """
    from .island_detection import identify_retained_islands

    # Against MT target
    islands_mt = identify_retained_islands(probe_seq, mt_target_seq, threshold)
    S_ECI_MT = compute_S_ECI(islands_mt)
    Zbox_MT = compute_Zbox_NN(islands_mt, temperature_C, salt_mM, threshold)

    # Against WT target
    islands_wt = identify_retained_islands(probe_seq, wt_target_seq, threshold)
    S_ECI_WT = compute_S_ECI(islands_wt)
    Zbox_WT = compute_Zbox_NN(islands_wt, temperature_C, salt_mM, threshold)

    # Selectivity ratios
    S_ECI_ratio = S_ECI_MT / S_ECI_WT if S_ECI_WT > 0 else float('inf')
    Zbox_ratio = Zbox_MT / Zbox_WT if Zbox_WT > 0 else float('inf')

    return {
        'S_ECI_MT': S_ECI_MT,
        'S_ECI_WT': S_ECI_WT,
        'S_ECI_ratio': S_ECI_ratio,
        'Zbox_MT': Zbox_MT,
        'Zbox_WT': Zbox_WT,
        'Zbox_ratio': Zbox_ratio,
        'n_islands_MT': len(islands_mt),
        'n_islands_WT': len(islands_wt),
        'islands_MT': islands_mt,
        'islands_WT': islands_wt,
    }
