# zbox_nn_engine

Topology-aware thermodynamics for DNA probe design.

A standalone Python package implementing three scoring layers for DNA probe
selectivity prediction, as described in Van der Snickt et al.,
*"Topology-aware thermodynamics for DNA probe design under fixed stringency."*

## Scoring layers

| Score | Formula | Description |
|-------|---------|-------------|
| **S_ECI** | `sum_i c_i^2` | Uniform-weight box-counting score. `c_i` is the retained-core length of island `i` after edge correction. |
| **N_box** | `(L-k+1)(L-k+2)/2` | Analytical count of contiguous sub-boxes of length >= k in an island of length L. |
| **Z_box,NN** | `sum_b exp[-dG_NN(b,T)/RT]` | NN-weighted partition function over retained islands. `dG_NN` is the SantaLucia 2004 nearest-neighbor free energy at temperature T and salt [Na+]. |

## Installation

```bash
# From the repository root:
pip install -e ./zbox_nn_engine

# Or simply add the directory to your Python path:
# import sys; sys.path.insert(0, "path/to/zbox_nn_engine")
```

Requires Python >= 3.8 and NumPy >= 1.20.

## Quick start

```python
from zbox_nn_engine import compute_all_scores

# Compute selectivity scores for a 17-mer probe
scores = compute_all_scores(
    probe_seq="ctcttgcctacgccaca",
    mt_target_seq="ctcttgcctacgccaca",   # mutant target (perfect match)
    wt_target_seq="ctcttgcctacgccacc",   # wild-type target (1 mismatch)
    temperature_C=68,
    salt_mM=1000,
    threshold=2,
)

print(f"S_ECI ratio:  {scores['S_ECI_ratio']:.4f}")
print(f"Zbox ratio:   {scores['Zbox_ratio']:.4f}")
print(f"Islands (MT): {scores['n_islands_MT']}")
print(f"Islands (WT): {scores['n_islands_WT']}")
```

## API reference

### `identify_retained_islands(probe_seq, target_seq, threshold=2)`

Find contiguous matched regions (retained islands) in a probe-target alignment.
Each island is trimmed by edge correction: `c_i = max(0, L_i - e_i)` where
`e_i` is the number of mismatch-exposed edges (0, 1, or 2). Only islands with
`c_i >= threshold` are returned.

**Returns:** list of dicts with `start`, `end`, `length`, `edges`, `retained_core`, `box_seq`.

### `compute_S_ECI(retained_islands)`

Uniform-weight box-counting score: `S_ECI = sum_i c_i^2`.

### `compute_N_box(L, k)`

Analytical sub-box count: `N_box(L,k) = (L-k+1)(L-k+2)/2`.

### `compute_Zbox_NN(retained_islands, temperature_C=68, salt_mM=1000, threshold=2, RT=None)`

NN-weighted partition function: `Z_box,NN = sum_b exp[-dG_NN(b,T)/RT]`.

### `compute_box_dG(box_seq, temperature_C=37, salt_mM=1000)`

NN free energy of a perfect duplex segment, using SantaLucia 2004 parameters
with terminal initiation, temperature correction, and salt correction.

### `compute_all_scores(probe_seq, mt_target_seq, wt_target_seq, temperature_C=68, salt_mM=1000, threshold=2)`

Convenience function that computes S_ECI and Z_box,NN for a probe against both
MT and WT targets, plus selectivity ratios.

**Returns:** dict with `S_ECI_MT`, `S_ECI_WT`, `S_ECI_ratio`, `Zbox_MT`, `Zbox_WT`, `Zbox_ratio`, `n_islands_MT`, `n_islands_WT`, `islands_MT`, `islands_WT`.

### `count_mismatches(probe_seq, target_seq)`

Count and locate mismatches between probe and target.

**Returns:** dict with `n_mismatches` and `positions` (1-indexed).

## Thermodynamic parameters

All NN parameters are from SantaLucia (2004), *Annu. Rev. Biophys. Biomol. Struct.*
33:415-440, Table 1 (unified NN parameters at 1 M NaCl, pH 7, 37 deg C).

- Gas constant: R = 1.987 cal/(mol K)
- 10 Watson-Crick NN dinucleotides + terminal initiation (AT/GC)
- Temperature correction: dG(T) = dH - T * dS
- Salt correction (SantaLucia 1998): dS_salt = 0.368 * (N-1) * ln([Na+])

## Future work

This engine is designed to be reusable beyond the manuscript's benchmark data.
Potential applications include:
- Web interface for interactive probe design
- Integration with high-throughput probe selection pipelines
- Extension to other hybridization models (e.g., including mismatch penalties)
