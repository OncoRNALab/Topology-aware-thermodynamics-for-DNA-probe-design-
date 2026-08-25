"""
Retained-island detection for topology-aware probe scoring.

A *retained island* is a maximal run of Watson-Crick matched positions
between a probe and its target.  Each island is trimmed by an edge
correction that removes positions exposed to a neighbouring mismatch:

    c_i = max(0, L_i - e_i)

where L_i is the raw island length and e_i is the number of mismatch-exposed
edges (0, 1, or 2).  Only islands with c_i >= k (the threshold) are retained.
"""

import numpy as np


def identify_retained_islands(probe_seq, target_seq, threshold=2):
    """
    Identify retained islands in a probe-target alignment.

    The probe and target are aligned position-by-position (no gaps).
    A position is *matched* when the probe base equals the target base
    (Watson-Crick pairing).  Consecutive matches form islands; each island
    is trimmed by edge correction and filtered by the threshold.

    Parameters
    ----------
    probe_seq : str
        Probe sequence (5'->3'), same length as target_seq.
    target_seq : str
        Target sequence (5'->3'), same length as probe_seq.
    threshold : int, optional
        Minimum retained-core length for an island to be kept (default 2).

    Returns
    -------
    list of dict
        Each dict has keys:
        - start    : 1-indexed start position
        - end      : 1-indexed end position
        - length   : raw island length (L_i)
        - edges    : number of mismatch-exposed edges (0, 1, or 2)
        - retained_core : c_i = max(0, length - edges)
        - box_seq  : the probe subsequence for this island
    """
    probe_seq = probe_seq.lower()
    target_seq = target_seq.lower()
    L = len(probe_seq)

    # Find matched positions
    matches = []
    for i in range(L):
        target_base = target_seq[i] if i < len(target_seq) else 'n'
        matches.append(probe_seq[i] == target_base)

    # Find islands of consecutive matches
    islands = []
    start = None
    for i in range(L):
        if matches[i]:
            if start is None:
                start = i
        else:
            if start is not None:
                islands.append((start, i - 1))
                start = None
    if start is not None:
        islands.append((start, L - 1))

    # Apply edge correction and threshold filter
    retained_islands = []
    for s, e in islands:
        length = e - s + 1
        edges = 0
        if s > 0:          # mismatch before the island
            edges += 1
        if e < L - 1:      # mismatch after the island
            edges += 1
        retained_core = max(0, length - edges)

        if retained_core >= threshold:
            retained_islands.append({
                'start': s + 1,          # 1-indexed
                'end': e + 1,
                'length': length,
                'edges': edges,
                'retained_core': retained_core,
                'box_seq': probe_seq[s:e + 1]
            })

    return retained_islands


def count_mismatches(probe_seq, target_seq):
    """
    Count and locate mismatches between probe and target.

    Returns
    -------
    dict with keys:
    - n_mismatches : int
    - positions    : list of 1-indexed mismatch positions
    """
    probe_seq = probe_seq.lower()
    target_seq = target_seq.lower()
    positions = []
    for i in range(len(probe_seq)):
        target_base = target_seq[i] if i < len(target_seq) else 'n'
        if probe_seq[i] != target_base:
            positions.append(i + 1)   # 1-indexed
    return {'n_mismatches': len(positions), 'positions': positions}
