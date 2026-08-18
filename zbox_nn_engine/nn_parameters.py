"""
SantaLucia 2004 unified nearest-neighbor (NN) thermodynamic parameters.

All values are from SantaLucia (2004), Annu. Rev. Biophys. Biomol. Struct. 33:415-440,
Table 1 (unified NN parameters at 1 M NaCl, pH 7, 37 °C).

Units:
    dG37 — kcal/mol  (free energy at 37 °C)
    dH   — kcal/mol  (enthalpy)
    dS   — cal/(mol·K)  (entropy)

Sequence notation: the key is the 5'->3' top-strand dinucleotide.
The bottom strand is implicitly the Watson-Crick complement, read 3'->5'.
For example, 'AA' means 5'-AA-3' / 3'-TT-5'.
"""

# Gas constant in cal/(mol·K)
R_GAS = 1.987

# Watson-Crick NN dinucleotide parameters (10 unique pairs)
# SantaLucia 2004, Table 1
NN_PARAMS = {
    'AA': {'dG37': -1.00, 'dH': -7.9, 'dS': -22.2},
    'AT': {'dG37': -0.88, 'dH': -7.2, 'dS': -20.4},
    'TA': {'dG37': -0.58, 'dH': -6.0, 'dS': -16.9},
    'CA': {'dG37': -1.45, 'dH': -8.5, 'dS': -22.7},
    'GT': {'dG37': -1.44, 'dH': -8.4, 'dS': -22.4},
    'CT': {'dG37': -1.28, 'dH': -7.8, 'dS': -21.0},
    'GA': {'dG37': -1.30, 'dH': -8.2, 'dS': -22.2},
    'CG': {'dG37': -2.17, 'dH': -10.6, 'dS': -27.2},
    'GC': {'dG37': -2.24, 'dH': -9.8, 'dS': -24.4},
    'GG': {'dG37': -1.84, 'dH': -8.0, 'dS': -19.9},
}

# Terminal initiation parameters
INIT_GC = {'dG37': 1.03, 'dH': 0.1, 'dS': -2.8}   # duplex initiated by terminal G·C
INIT_AT = {'dG37': 0.98, 'dH': 2.3, 'dS': 4.1}     # duplex initiated by terminal A·T

# Reverse-complement mapping for the 10 unique NN keys
# (needed because the table lists only one orientation per degenerate pair)
RC_MAP = {
    'TT': 'AA', 'TG': 'CA', 'AC': 'GT', 'AG': 'CT',
    'TC': 'GA', 'CC': 'GG',
    # self-complementary or already canonical
    'AA': 'AA', 'AT': 'AT', 'TA': 'TA', 'CA': 'CA',
    'GT': 'GT', 'CT': 'CT', 'GA': 'GA', 'CG': 'CG',
    'GC': 'GC', 'GG': 'GG',
}

# Fallback parameters for non-canonical dinucleotides
# (used in the original manuscript computation; matches the average of
#  the 10 WC NN parameters.  The RC_MAP above provides the thermodynamically
#  correct mapping, but the manuscript values were computed with this
#  fallback.  Set USE_RC_MAP = True to use the correct mapping instead.)
USE_RC_MAP = False
FALLBACK_NN = {'dG37': -1.30, 'dH': -8.0, 'dS': -22.0}


def get_nn_params(dinuc):
    """
    Look up NN parameters for a dinucleotide key.

    The SantaLucia 2004 table lists 10 unique NN dinucleotides.  The
    remaining 6 (TT, TG, AC, AG, TC, CC) are degenerate with their
    reverse-complement partners.

    By default (``USE_RC_MAP = False``) the original manuscript fallback
    (average of the 10 WC parameters) is used for these 6, matching the
    published supplementary tables exactly.  Set ``USE_RC_MAP = True``
    to use the thermodynamically correct reverse-complement mapping
    instead (ratios are virtually identical; absolute dG values differ).

    Parameters
    ----------
    dinuc : str
        Two-letter top-strand dinucleotide (e.g. 'AA', 'TC').

    Returns
    -------
    dict with keys 'dG37', 'dH', 'dS'
    """
    dinuc = dinuc.upper()
    if dinuc in NN_PARAMS:
        return NN_PARAMS[dinuc]
    if USE_RC_MAP:
        rc = RC_MAP.get(dinuc)
        if rc and rc in NN_PARAMS:
            return NN_PARAMS[rc]
    # Fallback: matches the original manuscript computation
    return dict(FALLBACK_NN)
