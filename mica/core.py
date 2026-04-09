# ════════════════════════════════════════════════════════════════════
#  MICA-Kernel v1.0 — Core Physics Engine
#  kyriosMICA Research Institute · Benin, West Africa
#  © 2026 Cyrille Egnon Davoh
#
#  Molecular Information through Coherent Architecture
#  Implementation of TQIM-Davoh (Théorie Quantique de l'Information
#  Moléculaire de Davoh) — Qudits-36 Formalism
#
#  The four Postulates:
#    I.  Indétermination  — Quantum superposition before measurement
#    II. Cohérence        — Heptadic entanglement MIH-21
#    III.Résonance        — Vibrational frequency ν_TIV = E_res/h
#    IV. Mutabilité       — Mutation probability P(-1)×P(+1)
# ════════════════════════════════════════════════════════════════════

import math
import numpy as np
from typing import Dict, List, Tuple, Optional, Union

# ── Physical Constants ────────────────────────────────────────────
PLANCK_H    = 6.626e-34    # J·s
KB          = 1.38e-23     # J/K
AVOGADRO    = 6.022e23

# ── CHARMM36 Force Field Parameters ──────────────────────────────
K_HB        = 20.0         # kcal/mol/Å²  — H-bond spring constant
DELTA       = 0.20         # Å            — discretization width
E_PAULI     = 0.40         # kcal/mol     — Pauli repulsion
KAP_AT      = 0.30         # intra-pair coupling AT
KAP_CG      = 0.50         # intra-pair coupling CG (stiffer)
LAM         = 0.12         # inter-base coupling (π-stacking)
KCAL_TO_J   = 4184.0       # 1 kcal/mol = 4184 J/mol

# ── Ternary State Definitions ─────────────────────────────────────
STATES_MAP = {
    'A': [-1,  0],          # Adenine  — 2 H-bonds
    'T': [ 0, +1],          # Thymine  — 2 H-bonds
    'C': [+1, -1,  0],      # Cytosine — 3 H-bonds
    'G': [-1, +1,  0],      # Guanine  — 3 H-bonds
}
BASE_PAIR  = {'A':'T','T':'A','C':'G','G':'C'}
N_BONDS    = {'A':2,'T':2,'C':3,'G':3}


# ════════════════════════════════════════════════════════════════════
#  MODULE 1 — THERMODYNAMICS
#  Postulat de l'Indétermination
# ════════════════════════════════════════════════════════════════════

def kbt(T_K: float) -> float:
    """Thermal energy k_B·T in kcal/mol"""
    return KB * T_K * AVOGADRO / KCAL_TO_J

def delta_E(state: int) -> float:
    """Energy of H-bond state s ∈ {-1, 0, +1} in kcal/mol"""
    if state == 0:  return 0.0
    if state == -1: return 0.5 * K_HB * DELTA**2 + E_PAULI
    return 0.5 * K_HB * DELTA**2   # state == +1

def boltzmann_populations(T_K: float = 300.0,
                          pH: float = 7.4,
                          ionic_strength: float = 0.15) -> Dict[int,float]:
    """
    Boltzmann populations P(s) for s ∈ {-1, 0, +1}
    at temperature T_K, pH and ionic strength I.

    pH correction: acid pH shifts population toward state -1
      (compression — proton tunneling facilitated)
    Ionic strength: modifies coupling constants via Debye-Hückel
      (not applied to populations directly but passed to MRK)
    """
    kT = kbt(T_K)

    # pH correction: ΔE(-1) decreases with acidity
    # At pH 7.4 (physiological): no correction
    # At pH 5.0 (stress/acidosis): -0.15 kcal/mol on state -1
    pH_correction = 0.0
    if pH < 7.4:
        pH_correction = -0.15 * (7.4 - pH) / 2.4  # linear interpolation

    energies = {
        -1: delta_E(-1) + pH_correction,
         0: delta_E(0),
        +1: delta_E(+1)
    }

    Z = sum(math.exp(-e/kT) for e in energies.values())
    return {s: math.exp(-e/kT)/Z for s, e in energies.items()}

def quantum_amplitudes(T_K: float = 300.0,
                       pH: float = 7.4,
                       ionic_strength: float = 0.15) -> Dict[int,float]:
    """
    Quantum amplitude vector for H-bond state (Postulat I)
    α(s) = √P(s)  such that Σ|α(s)|² = 1

    |ψ_i⟩ = α(-1)|−1⟩ + α(0)|0⟩ + α(+1)|+1⟩
    """
    P = boltzmann_populations(T_K, pH, ionic_strength)
    return {s: math.sqrt(p) for s, p in P.items()}


# ════════════════════════════════════════════════════════════════════
#  MODULE 2 — EVK (Espace Vectoriel Kyriosmica)
# ════════════════════════════════════════════════════════════════════

def evk_vector(codon: str) -> List[int]:
    """
    EVK vector of a codon in T₃ⁿ space.
    Uses default Watson-Crick states (most probable configuration).

    V_codon = V_B1 ⊕ V_B2 ⊕ V_B3  (direct sum)
    """
    v = []
    for base in codon.upper():
        v.extend(STATES_MAP.get(base, [0, 0]))
    return v

def evk_dimension(codon: str) -> int:
    """Dimension of EVK space for this codon"""
    return sum(N_BONDS.get(b, 2) for b in codon.upper())

def evk_states_count(codon: str) -> int:
    """Total number of accessible states: 3^n"""
    return 3 ** evk_dimension(codon)

def evk_composition(codon: str) -> Tuple[int,int]:
    """Returns (n_AT_bases, n_CG_bases) for a codon"""
    at = sum(1 for b in codon.upper() if b in 'AT')
    cg = sum(1 for b in codon.upper() if b in 'CG')
    return at, cg

def evk_label(codon: str) -> str:
    """Human-readable EVK space label, e.g. 'AT×2+CG'"""
    at, cg = evk_composition(codon)
    parts = []
    if at > 0: parts.append(f'AT×{at}' if at > 1 else 'AT')
    if cg > 0: parts.append(f'CG×{cg}' if cg > 1 else 'CG')
    return '+'.join(parts)


# ════════════════════════════════════════════════════════════════════
#  MODULE 3 — ERK (Espace de Résonance Kyriosmica)
# ════════════════════════════════════════════════════════════════════

def erk_energy(codon: str,
               T_K: float = 300.0,
               pH: float = 7.4,
               ionic_strength: float = 0.15,
               mode: str = 'classical') -> Dict:
    """
    ERK resonance energy for a codon.

    Modes:
      'classical' : uses default Watson-Crick states
      'quantum'   : uses expected value over Boltzmann distribution
    """
    v = evk_vector(codon)

    if mode == 'classical':
        E = sum(delta_E(t) for t in v)
        return {'E_ERK': round(E, 4), 'mode': 'classical'}

    elif mode == 'quantum':
        P = boltzmann_populations(T_K, pH, ionic_strength)
        # Expected ERK over quantum distribution
        E_mean = 0.0
        E_sq   = 0.0
        n = evk_dimension(codon)
        for s in [-1, 0, 1]:
            e = delta_E(s)
            # Each bond independently draws state s
            for bond_idx in range(n):
                E_mean += P[s] * e
                E_sq   += P[s] * e**2
        sigma_sq = E_sq/n - (E_mean/n)**2
        return {
            'E_ERK_mean':  round(E_mean, 4),
            'E_ERK_sigma': round(math.sqrt(max(0, sigma_sq)), 4),
            'mode': 'quantum'
        }


# ════════════════════════════════════════════════════════════════════
#  MODULE 4 — MRK (Matrice de Résonance Kyriosmica)
# ════════════════════════════════════════════════════════════════════

def mrk_matrix(codon: str,
               ionic_strength: float = 0.15) -> np.ndarray:
    """
    MRK matrix M ∈ ℝⁿˣⁿ for a codon.

    Diagonal blocks: intra-pair coupling
      M[a,a] = t_a²
      M[a,b] = κ·t_a·t_b  (a≠b, same base)

    Off-diagonal blocks: inter-base coupling (π-stacking)
      M[i,j] = λ·t_i·t_j  (different bases)

    Ionic strength modifies λ via Debye-Hückel screening:
      λ_eff = λ · exp(-κ_D · d)  where κ_D ∝ √I
    """
    v   = evk_vector(codon)
    n   = len(v)
    M   = np.zeros((n, n))

    # Debye-Hückel correction on inter-base coupling
    kappa_D = 3.28e9 * math.sqrt(ionic_strength)  # m⁻¹
    d_stack = 3.4e-10  # m — base stacking distance
    lam_eff = LAM * math.exp(-kappa_D * d_stack)

    bases = [(b, len(STATES_MAP.get(b, [0,0]))) for b in codon.upper()]
    offsets = [0]
    for _, nb in bases:
        offsets.append(offsets[-1] + nb)

    # Diagonal blocks (intra-pair)
    for bi, (base, nb) in enumerate(bases):
        oi  = offsets[bi]
        kap = KAP_CG if base in 'CG' else KAP_AT
        for a in range(nb):
            for b in range(nb):
                ta, tb = v[oi+a], v[oi+b]
                M[oi+a, oi+b] = ta**2 if a==b else kap*ta*tb

    # Off-diagonal blocks (inter-base, π-stacking)
    for bi, (_, nb_i) in enumerate(bases):
        for bj, (_, nb_j) in enumerate(bases):
            if bj <= bi:
                continue
            oi, oj = offsets[bi], offsets[bj]
            for a in range(nb_i):
                for b in range(nb_j):
                    val = lam_eff * v[oi+a] * v[oj+b]
                    M[oi+a, oj+b] = val
                    M[oj+b, oi+a] = val

    # Tikhonov regularization (semi-definite positive)
    min_eig = np.linalg.eigvalsh(M).min()
    if min_eig < 0:
        M += (-min_eig + 1e-8) * np.eye(n)

    return M

def mrk_invariants(codon: str,
                   ionic_strength: float = 0.15) -> Dict:
    """
    Compute all MRK invariants for a codon:
    - Rank (conformational signature)
    - Trace (total coupling energy)
    - Frobenius norm (spectral intensity)
    - Eigenvalues {λ_i} (vibrational modes)
    - E_res = λ_max (dominant resonance frequency)
    """
    M       = mrk_matrix(codon, ionic_strength)
    eigvals = np.linalg.eigvalsh(M)
    rank    = int(np.linalg.matrix_rank(M, tol=1e-6))
    E_res   = float(np.max(eigvals))

    return {
        'rank':       rank,
        'trace':      round(float(np.trace(M)), 4),
        'frobenius':  round(float(np.linalg.norm(M, 'fro')), 4),
        'E_res':      round(E_res, 4),
        'eigenvalues': [round(float(e), 4) for e in sorted(eigvals, reverse=True)],
        'nu_TIV_THz': round(E_res * KCAL_TO_J / AVOGADRO / PLANCK_H / 1e12, 4)
    }


# ════════════════════════════════════════════════════════════════════
#  MODULE 5 — POSTULAT DE LA MUTABILITÉ
# ════════════════════════════════════════════════════════════════════

def mutation_probability(codon: str,
                         T_K: float = 300.0,
                         pH: float = 7.4,
                         ionic_strength: float = 0.15) -> Dict:
    """
    Postulat de la Mutabilité (TQIM-Davoh):

    P(mutation_i) ∝ P_i(-1) · P_i(+1)

    A spontaneous mutation requires the proton to be accessible
    to both extreme states — compressed (-1, enabling tunneling)
    AND extended (+1, stabilizing the tautomer).

    Returns individual bond probabilities and codon total.
    """
    P  = boltzmann_populations(T_K, pH, ionic_strength)
    at, cg = evk_composition(codon)
    n_total_bonds = at * N_BONDS['A'] + cg * N_BONDS['C']

    # Per-bond mutation probability
    p_mut_per_bond = P[-1] * P[1]

    # Codon mutation probability (product over all bonds)
    # AT bonds and CG bonds have same P in Phase 1 (CHARMM36 additive)
    # Phase 2 (Drude): will differentiate AT vs CG
    p_mut_codon = p_mut_per_bond ** n_total_bonds

    # Mutability score: normalized rank-based risk
    inv = mrk_invariants(codon, ionic_strength)
    risk = 'HIGH' if inv['rank'] <= 3 else 'MEDIUM' if inv['rank'] == 4 else 'LOW'

    return {
        'P_mut_per_bond':  round(p_mut_per_bond, 6),
        'P_mut_codon':     round(p_mut_codon, 10),
        'P_mut_pct':       round(p_mut_per_bond * 100, 4),
        'n_bonds':         n_total_bonds,
        'mrk_rank':        inv['rank'],
        'risk_level':      risk,
        'note': 'Phase 1: AT/CG identical. Phase 2 (Drude) will differentiate.'
    }


# ════════════════════════════════════════════════════════════════════
#  MODULE 6 — MIH-21 (Postulat de la Cohérence)
# ════════════════════════════════════════════════════════════════════

def mih21_analysis(codons_7: List[str],
                   T_K: float = 300.0,
                   pH: float = 7.4,
                   ionic_strength: float = 0.15) -> Dict:
    """
    MIH-21 Heptadic Entanglement Analysis (Postulat II + III)

    Input: exactly 7 consecutive codons (21 base pairs = 2 B-DNA turns)
    Output: resonance cavity signature

    The MIH-21 cavity forms a topologically closed resonant cavity
    because 21 bp = 2 × 10.5 bp/turn (B-DNA geometry).
    Standing waves can establish inside this cavity.
    """
    if len(codons_7) != 7:
        raise ValueError(f"MIH-21 requires exactly 7 codons, got {len(codons_7)}")

    # Compute MRK invariants for each codon
    codon_signatures = []
    for codon in codons_7:
        inv = mrk_invariants(codon, ionic_strength)
        erk = erk_energy(codon, T_K, pH, ionic_strength, mode='quantum')
        mut = mutation_probability(codon, T_K, pH, ionic_strength)
        codon_signatures.append({
            'codon':    codon,
            'evk_dim':  evk_dimension(codon),
            'rank':     inv['rank'],
            'E_res':    inv['E_res'],
            'frobenius':inv['frobenius'],
            'P_mut_pct':mut['P_mut_pct'],
            'risk':     mut['risk_level']
        })

    # Cavity resonance frequency (dominant mode)
    # ν_TIV = Σ(E_res_i) / (h · n_codons) — mean field approximation
    E_res_values = [s['E_res'] for s in codon_signatures]
    E_res_total  = sum(E_res_values)
    E_res_mean   = E_res_total / 7

    # Standing wave analysis
    # L = 21 × 3.4 Å = 71.4 Å (B-DNA inter-base distance)
    L_angstrom   = 21 * 3.4   # Å
    lambda_fund  = 2 * L_angstrom  # fundamental mode wavelength

    # Convert E_res to THz
    nu_TIV = (E_res_mean * KCAL_TO_J / AVOGADRO) / PLANCK_H / 1e12

    # Conformational stability of the cavity
    ranks = [s['rank'] for s in codon_signatures]
    high_risk = [s for s in codon_signatures if s['risk'] == 'HIGH']

    return {
        'codons':          codons_7,
        'L_angstrom':      L_angstrom,
        'lambda_fund_A':   round(lambda_fund, 1),
        'E_res_mean_kcal': round(E_res_mean, 4),
        'nu_TIV_THz':      round(nu_TIV, 4),
        'mean_rank':       round(sum(ranks)/7, 3),
        'min_rank':        min(ranks),
        'max_rank':        max(ranks),
        'high_risk_codons': len(high_risk),
        'cavity_stability': 'UNSTABLE' if len(high_risk) >= 3
                            else 'MODERATE' if len(high_risk) >= 1
                            else 'STABLE',
        'codon_signatures': codon_signatures
    }


# ════════════════════════════════════════════════════════════════════
#  MODULE 7 — MONTE CARLO QUANTUM COLLAPSE
# ════════════════════════════════════════════════════════════════════

def quantum_collapse_mc(codon: str,
                        N: int = 2000,
                        T_K: float = 300.0,
                        pH: float = 7.4,
                        ionic_strength: float = 0.15) -> Dict:
    """
    Monte Carlo simulation of quantum collapses (Postulat I + Measurement).

    Simulates N independent biological measurement events.
    Each event collapses each H-bond independently according to P(s).

    Returns: distribution of MRK ranks, conformational variance,
    conformational entropy.
    """
    P   = boltzmann_populations(T_K, pH, ionic_strength)
    probs = [P[-1], P[0], P[1]]
    n   = evk_dimension(codon)

    np.random.seed(42)
    ranks  = []
    frobs  = []
    eigmax = []

    for _ in range(N):
        # Quantum collapse: each bond draws state from P(s)
        v_collapsed = list(np.random.choice([-1, 0, 1], size=n, p=probs))
        # Build MRK from collapsed state
        M = mrk_matrix_from_vector(v_collapsed, codon, ionic_strength)
        eigvals = np.linalg.eigvalsh(M)
        ranks.append(int(np.linalg.matrix_rank(M, tol=1e-6)))
        frobs.append(float(np.linalg.norm(M, 'fro')))
        eigmax.append(float(np.max(eigvals)))

    ranks_arr = np.array(ranks)
    frobs_arr = np.array(frobs)

    # Conformational entropy H(rank) in bits
    unique, counts = np.unique(ranks_arr, return_counts=True)
    probs_rank = counts / N
    H_rank = float(-np.sum(probs_rank * np.log2(probs_rank + 1e-12)))

    return {
        'N_simulations':   N,
        'rank_mean':       round(float(np.mean(ranks_arr)), 3),
        'rank_std':        round(float(np.std(ranks_arr)), 3),
        'rank_min':        int(np.min(ranks_arr)),
        'rank_max':        int(np.max(ranks_arr)),
        'rank_classical':  mrk_invariants(codon, ionic_strength)['rank'],
        'var_M':           round(float(np.var(frobs_arr)), 6),
        'E_res_mean':      round(float(np.mean(eigmax)), 4),
        'E_res_std':       round(float(np.std(eigmax)), 4),
        'H_rank_bits':     round(H_rank, 3),
        'n_eff_conformations': round(2**H_rank, 2),
        'rank_distribution': {
            int(r): round(float(c/N), 4)
            for r, c in zip(unique, counts)
        }
    }


def mrk_matrix_from_vector(v: List[int],
                            codon: str,
                            ionic_strength: float = 0.15) -> np.ndarray:
    """Build MRK matrix from a given collapsed vector (used in MC)"""
    n   = len(v)
    M   = np.zeros((n, n))

    kappa_D = 3.28e9 * math.sqrt(ionic_strength)
    d_stack = 3.4e-10
    lam_eff = LAM * math.exp(-kappa_D * d_stack)

    bases   = [(b, len(STATES_MAP.get(b,[0,0]))) for b in codon.upper()]
    offsets = [0]
    for _, nb in bases: offsets.append(offsets[-1] + nb)

    for bi, (base, nb) in enumerate(bases):
        oi  = offsets[bi]
        kap = KAP_CG if base in 'CG' else KAP_AT
        for a in range(nb):
            for b in range(nb):
                ta, tb = v[oi+a], v[oi+b]
                M[oi+a, oi+b] = ta**2 if a==b else kap*ta*tb

    for bi, (_, nb_i) in enumerate(bases):
        for bj, (_, nb_j) in enumerate(bases):
            if bj <= bi: continue
            oi, oj = offsets[bi], offsets[bj]
            for a in range(nb_i):
                for b in range(nb_j):
                    val = lam_eff * v[oi+a] * v[oj+b]
                    M[oi+a, oj+b] = val
                    M[oj+b, oi+a] = val

    min_eig = np.linalg.eigvalsh(M).min()
    if min_eig < 0:
        M += (-min_eig + 1e-8) * np.eye(n)
    return M


# ════════════════════════════════════════════════════════════════════
#  MODULE 8 — FULL SEQUENCE ANALYSIS
# ════════════════════════════════════════════════════════════════════

def analyze_sequence(sequence:        str,
                     T_K:             float = 300.0,
                     pH:              float = 7.4,
                     ionic_strength:  float = 0.15,
                     pressure_atm:    float = 1.0,
                     run_mc:          bool  = False,
                     mc_samples:      int   = 1000) -> Dict:
    """
    Full TQIM-Davoh analysis of a DNA/RNA sequence.

    INPUT:
      sequence       : DNA string (FASTA content or raw sequence)
      T_K            : Temperature in Kelvin (default 310K body temp)
      pH             : pH of biological medium (default 7.4)
      ionic_strength : Ionic strength in mol/L (default 0.15 = serum)
      pressure_atm   : Pressure in atm (default 1.0)
      run_mc         : Run Monte Carlo quantum collapse simulation
      mc_samples     : Number of MC samples per codon (if run_mc)

    OUTPUT:
      Complete TQIM-Davoh signature including:
      - EVK, ERK, MRK per codon
      - Quantum amplitudes
      - Mutation probability map (Postulat IV)
      - MIH-21 cavity analysis (Postulat II + III)
      - Summary statistics
      - High-risk site predictions
    """

    # ── Parse sequence ───────────────────────────────────────────
    seq = _parse_sequence(sequence)
    if len(seq) < 3:
        raise ValueError("Sequence must be at least 3 nucleotides (1 codon)")

    # Extract codons
    codons = [seq[i:i+3] for i in range(0, len(seq)-2, 3)
              if len(seq[i:i+3]) == 3
              and all(b in 'ATCGN' for b in seq[i:i+3])]

    # Remove stop codons
    stop_codons = {'TAA','TAG','TGA'}
    codons_cds  = []
    for c in codons:
        if c.upper() in stop_codons:
            break
        codons_cds.append(c.upper())

    if not codons_cds:
        raise ValueError("No valid codons found in sequence")

    # ── Compute quantum amplitudes (same for all bonds, Phase 1) ─
    P   = boltzmann_populations(T_K, pH, ionic_strength)
    amp = quantum_amplitudes(T_K, pH, ionic_strength)

    # ── Analyze each codon ───────────────────────────────────────
    codon_results = []
    for i, codon in enumerate(codons_cds):
        inv = mrk_invariants(codon, ionic_strength)
        erk = erk_energy(codon, T_K, pH, ionic_strength, 'quantum')
        mut = mutation_probability(codon, T_K, pH, ionic_strength)

        result = {
            'position':     i + 1,
            'codon':        codon,
            'evk_label':    evk_label(codon),
            'evk_dim':      evk_dimension(codon),
            'evk_states':   evk_states_count(codon),
            'erk_mean':     erk['E_ERK_mean'],
            'rank':         inv['rank'],
            'trace':        inv['trace'],
            'frobenius':    inv['frobenius'],
            'E_res':        inv['E_res'],
            'nu_TIV_THz':   inv['nu_TIV_THz'],
            'P_mut_pct':    mut['P_mut_pct'],
            'risk_level':   mut['risk_level'],
        }

        if run_mc:
            mc = quantum_collapse_mc(codon, mc_samples, T_K, pH, ionic_strength)
            result['rank_quantum']  = mc['rank_mean']
            result['var_M']         = mc['var_M']
            result['H_rank_bits']   = mc['H_rank_bits']
            result['n_eff_conf']    = mc['n_eff_conformations']

        codon_results.append(result)

    # ── MIH-21 analysis ─────────────────────────────────────────
    mih21_windows = []
    for i in range(0, len(codons_cds) - 6, 7):
        window = codons_cds[i:i+7]
        if len(window) == 7:
            mih = mih21_analysis(window, T_K, pH, ionic_strength)
            mih['window_start'] = i + 1
            mih['window_end']   = i + 7
            mih21_windows.append(mih)

    # ── Summary statistics ───────────────────────────────────────
    ranks    = [r['rank']    for r in codon_results]
    erks     = [r['erk_mean'] for r in codon_results]
    E_res_v  = [r['E_res']   for r in codon_results]
    high_risk = [r for r in codon_results if r['risk_level'] == 'HIGH']
    med_risk  = [r for r in codon_results if r['risk_level'] == 'MEDIUM']

    summary = {
        'n_codons':         len(codons_cds),
        'n_nucleotides':    len(seq),
        'conditions': {
            'T_K':            T_K,
            'T_C':            round(T_K - 273.15, 1),
            'pH':             pH,
            'ionic_strength': ionic_strength,
            'pressure_atm':   pressure_atm,
        },
        'amplitudes': {
            'alpha_minus1': round(amp[-1], 4),
            'alpha_0':      round(amp[0],  4),
            'alpha_plus1':  round(amp[1],  4),
            'P_minus1':     round(P[-1],   4),
            'P_0':          round(P[0],    4),
            'P_plus1':      round(P[1],    4),
        },
        'mrk_rank_mean':    round(np.mean(ranks),  3),
        'mrk_rank_std':     round(np.std(ranks),   3),
        'erk_mean':         round(np.mean(erks),   4),
        'E_res_mean':       round(np.mean(E_res_v),4),
        'high_risk_sites':  len(high_risk),
        'medium_risk_sites':len(med_risk),
        'low_risk_sites':   len(codons_cds) - len(high_risk) - len(med_risk),
        'n_mih21_windows':  len(mih21_windows),
    }

    return {
        'metadata': {
            'tool':        'MICA-Kernel v1.0',
            'institution': 'kyriosMICA · Benin, West Africa',
            'framework':   'TQIM-Davoh · Qudits-36',
            'author':      'Cyrille Egnon Davoh',
            'version':     '1.0.0',
            'doi':         'https://zenodo.org/kyriosmica',
        },
        'summary':         summary,
        'codons':          codon_results,
        'mih21_windows':   mih21_windows,
        'high_risk_sites': [
            {'pos': r['position'], 'codon': r['codon'],
             'rank': r['rank'], 'E_res': r['E_res']}
            for r in high_risk
        ],
        'predictions': {
            'method': 'Postulat de la Mutabilité (TQIM-Davoh)',
            'rule':   'MRK rank ≤ 3 → HIGH mutation risk',
            'sites':  [r['position'] for r in high_risk],
            'note':   'Pre-registered prediction April 2026 — kyriosMICA'
        }
    }


# ── Utility: parse FASTA or raw sequence ─────────────────────────
def _parse_sequence(raw: str) -> str:
    """Accept FASTA format or raw nucleotide string"""
    lines = raw.strip().split('\n')
    seq   = ''
    for line in lines:
        if line.startswith('>'):
            continue  # Skip FASTA header
        seq += line.strip().upper()
    # Keep only valid nucleotides
    return ''.join(b for b in seq if b in 'ATCGNU')
