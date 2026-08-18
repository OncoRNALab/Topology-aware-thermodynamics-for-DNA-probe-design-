import numpy as np
import itertools

def DG_model(probe, target, T, Salt=(1,0), init=False, ATend_correction=False, edges=False, l_e=4, MMinteraction=False, l_a=4):
# calculate the DG associated with hybridisation of probe to target, using the Nearest Neighbor model with two emperical corrections
    '''
    probe 5'->3', target 3'->5'
    Salt concentration in M (Na, Mg)

    Approach:
        1. first calculate PM DG
        2. then determine MM's and penalties (DG(MM)-DG(PM)) of the mismatched triplets
        3. then add Edge effects (penalties for MM are weighted based on distance to the edge of the probe)
        4. then add MMinteraction (penalties for MM are weighted based on distance to each other)
        5. return total DG
    '''

    dg = 0
    MM_locs = MMindex(probe, target)

    if any(MM_locs[i + 1] - MM_locs[i] == 1 for i in range(len(MM_locs) - 1)):
        return np.nan # MM next to another MM situations are not considered

    if (0 in MM_locs) or (len(probe)-1 in MM_locs):
        return np.nan # MM at the edge of the probe are not considered

    # dg if the target were fully complementary to the probe
    dg += DG_PM(probe, T, Salt)

    # select the part of the duplex with one nt left and right of the MM, for each MM.
    MM_triplets = [(probe[i-1:i+2], target[i-1:i+2]) for i in MM_locs]
    MM_penalties = []
    for top, bottom in MM_triplets:
        MM_pair1 = top[0]+top[1]+bottom[1]
        MM_pair2 = bottom[2]+bottom[1]+top[1]

        # PM dg (to be subtracted)
        PM_dg = DG_PM(top, T, Salt)

        # add MM dg
        MM_dg = DG_MM(MM_pair2, T, Salt) + DG_MM(MM_pair1, T, Salt)
        MM_penalties += [MM_dg - PM_dg]

    if edges:
        # edge effect (penalties for MM are weighted based on distance to the edge of the probe)
        # l_e = persistence length of edge effects
        MM_penalties = [penalty*(1-np.exp(-(i+1)/l_e))*(1-np.exp(-(len(probe)-i)/l_e)) for i, penalty in zip(MM_locs, MM_penalties)]

    if MMinteraction and len(MM_locs)>=2:
        # MM interaction (penalties for MM are weighted based on distance to each other)
        # l_a = MM interaction persistence length
        MM_penalties2 = [0]
        for i, j in itertools.combinations(range(len(MM_locs)), 2):
            alpha_ij = 0.5*np.exp(-np.abs(MM_locs[i]-MM_locs[j])/l_a)
            MM_penalties2 += [(MM_penalties[i] + MM_penalties[j]) * (1 - alpha_ij) / (len(MM_locs)-1)]
        MM_penalties = MM_penalties2
    
    dg += sum(MM_penalties)


    if init: dg += 0.2-T*(-5.7)/1000
    if ATend_correction:
        n_AT_ends = 0
        if probe[0]=='A' and target[0]=='T': n_AT_ends += 1
        if probe[-1]=='A' and target[-1]=='T': n_AT_ends += 1
        dg += n_AT_ends * (2.2 - T*6.9/1000)
    return dg

def DG_PM(probe, T, Salt=(1,0)):
    probe=probe.upper()
    dg = 0
    for i in range(len(probe)-1):
        pair = probe[i:i+2] # inspect pairs of consecutive probe letters
        dg += DG_pair(pair, T, Salt)
    return dg

# Parameters from SantaLucia 2004
def DG_pair(pair, T, Salt=(1,0)):
    # pair is neighboring bases on same strand bound to their complement
    # e.g. 'AG' means 5'-AG-3'/3'-TC-5'
    dg = 0
    dh = 0
    ds = 0
    Na, Mg = Salt
    pair = pair.upper()
    match pair:
        case 'AA': dh =  7.6; ds = 21.3
        case 'AC': dh =  8.4; ds = 22.4
        case 'AG': dh =  7.8; ds = 21.0
        case 'AT': dh =  7.2; ds = 20.4
        case 'CA': dh =  8.5; ds = 22.7
        case 'CC': dh =  8.0; ds = 19.9
        case 'CG': dh = 10.6; ds = 27.2
        case 'CT': dh =  7.8; ds = 21.0
        case 'GA': dh =  8.2; ds = 22.2
        case 'GC': dh =  9.8; ds = 24.4
        case 'GG': dh =  8.0; ds = 19.9
        case 'GT': dh =  8.4; ds = 22.4
        case 'TA': dh =  7.2; ds = 21.3
        case 'TC': dh =  8.2; ds = 22.2
        case 'TG': dh =  8.5; ds = 22.7
        case 'TT': dh =  7.6; ds = 21.3
    # salt correction
    ds = ds + 0.368 * 0.5 * np.log(Na + Mg*2)
    return -dh+T*ds/1000

def DG_MM(MM, T, Salt=(1,0)):
    # MM is a three letter tuple with first two on top strand, and last two forming the mismatch
    # e.g. 'GGT' means 5'-GG-3'/3'-CT-5' (where G-T is the mismatch and the bottom C is complementary to the first G)
    dg = 0
    dh = 0
    ds = 0
    Na, Mg = Salt
    MM = MM.upper()
    match MM:
        case 'AAA': dh =  1.2; ds =   1.7
        case 'ACA': dh =  5.3; ds =  14.6
        case 'AGA': dh = -0.7; ds =  -2.3
        case 'CAA': dh = -0.9; ds =  -4.2
        case 'CCA': dh =  0.6; ds =  -0.6
        case 'CGA': dh = -4.0; ds = -13.2
        case 'GAA': dh = -2.9; ds =  -9.8
        case 'GCA': dh = -0.7; ds =  -3.8
        case 'GGA': dh =  0.5; ds =   3.2
        case 'TAA': dh =  4.7; ds =  12.9
        case 'TCA': dh =  7.6; ds =  20.2
        case 'TGA': dh =  3.0; ds =   7.4

        case 'AAC': dh =  2.3; ds =   4.6
        case 'ACC': dh =  0.0; ds =  -4.4
        case 'ATC': dh = -1.2; ds =  -6.2
        case 'CAC': dh =  1.9; ds =   3.7
        case 'CCC': dh = -1.5; ds =  -7.2
        case 'CTC': dh = -1.5; ds =  -6.1
        case 'GAC': dh =  5.2; ds =  14.2
        case 'GCC': dh =  3.6; ds =   8.9
        case 'GTC': dh =  5.2; ds =  13.5
        case 'TAC': dh =  3.4; ds =   8.0
        case 'TCC': dh =  6.1; ds =  16.4
        case 'TTC': dh =  1.0; ds =   0.7

        case 'AAG': dh = -0.6; ds =  -2.3
        case 'AGG': dh = -3.1; ds =  -9.5
        case 'ATG': dh = -2.5; ds =  -8.3
        case 'CAG': dh = -0.7; ds =  -2.3
        case 'CGG': dh = -4.9; ds = -15.3
        case 'CTG': dh = -2.8; ds =  -8.0
        case 'GAG': dh = -0.6; ds =  -1.0
        case 'GGG': dh = -6.0; ds = -15.8
        case 'GTG': dh = -4.4; ds = -12.3
        case 'TAG': dh =  0.7; ds =   0.7
        case 'TGG': dh =  1.6; ds =   3.6
        case 'TTG': dh = -1.3; ds =  -5.3

        case 'ACT': dh =  0.7; ds =   0.2
        case 'AGT': dh =  1.0; ds =   0.9
        case 'ATT': dh = -2.7; ds = -10.8
        case 'CCT': dh = -0.8; ds =  -4.5
        case 'CGT': dh = -4.1; ds = -11.7
        case 'CTT': dh = -5.0; ds = -15.8
        case 'GCT': dh =  2.3; ds =   5.4
        case 'GGT': dh =  3.3; ds =  10.4
        case 'GTT': dh = -2.2; ds =  -8.4
        case 'TCT': dh =  1.2; ds =   0.7
        case 'TGT': dh = -0.1; ds =  -1.7
        case 'TTT': dh =  0.2; ds =  -1.5
    # salt correction
    ds = ds + 0.368 * 0.5 * np.log(Na + Mg*2)
    return dh-T*ds/1000

def MMindex(probe, target):
    # returns array with indices of MM locations
    assert(len(probe) == len(target))
    return [i for i in range(len(target)) if compl(probe[i]) != target[i]]

def compl(seq, reverse=False):
    # returns WC complement of given sequence string
    seq_no_spaces=''.join(seq.split()).upper()
    WC = list(seq_no_spaces)
    for i, l in enumerate(seq_no_spaces):
        match l:
            case 'A': WC[i] = 'T'
            case 'T': WC[i] = 'A'
            case 'G': WC[i] = 'C'
            case 'C': WC[i] = 'G'
            case _: raise ValueError
    return ''.join(WC)[::-1] if reverse else ''.join(WC)
