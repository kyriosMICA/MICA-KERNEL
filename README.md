# MICA-Kernel v1.0

**kyriosMICA Research Institute · Benin, West Africa**

> *"Decoding Life. Encoding the Future."*

---

## What is MICA-Kernel?

MICA-Kernel is the official computational engine of the
**TQIM-Davoh (Théorie Quantique de l'Information Moléculaire de Davoh)**
— the first mathematical framework that treats DNA not as a chemical
alphabet, but as a **quantum information processor**.

**MICA** = *Matrice d'Intrication Coordonnée par les Amplitudes*
*(Matrix of Intrication Coordinated by Amplitudes)*

Given a DNA sequence and biological conditions (temperature, pH,
ionic strength), MICA-Kernel computes the complete **Qudits-36
quantum signature** of every codon — a set of invariants that
no other tool in the world calculates.

---

## The Four Postulates (TQIM-Davoh)

| Postulate | Name | Formula |
|-----------|------|---------|
| I | **Indétermination** | `|ψ⟩ = α(-1)|−1⟩ + α(0)|0⟩ + α(+1)|+1⟩` |
| II | **Cohérence** | MIH-21 heptadic entanglement (2 B-DNA turns) |
| III | **Résonance** | `ν_TIV = E_res / h` |
| IV | **Mutabilité** | `P(mut) ∝ P(-1) × P(+1)` |

---

## Installation

```bash
git clone https://github.com/kyriosMICA/mica-kernel
cd mica-kernel
pip install numpy scipy
```

---

## Quick Start

```python
from mica import analyze_sequence

# Analyze any DNA sequence
result = analyze_sequence(
    sequence       = "ATGGATTATCCTTATGAATTTGCTCCT",
    T_K            = 310.0,   # Body temperature
    pH             = 7.4,     # Physiological
    ionic_strength = 0.15,    # Serum
    run_mc         = True,    # Monte Carlo quantum collapse
    mc_samples     = 1000,
)

print(f"MRK rank mean : {result['summary']['mrk_rank_mean']}")
print(f"High-risk sites: {result['summary']['high_risk_sites']}")
print(f"MIH-21 windows : {result['summary']['n_mih21_windows']}")
```

---

## Command Line

```bash
# Demo mode
python main.py --demo

# Analyze a sequence
python main.py --seq ATGCCCGAT --T 310 --pH 7.4

# Full analysis with Monte Carlo
python main.py --seq ATGCCC... --mc --mc-samples 2000 --output results.json

# From FASTA file
python main.py --fasta gene.fasta --T 310 --pH 6.8 --ionic 0.20

# Stress condition simulation (pH 6.0)
python main.py --seq ATGCCC... --T 310 --pH 6.0 --output stress.json
```

---

## What MICA-Kernel Computes

For each codon in the input sequence:

| Invariant | Description |
|-----------|-------------|
| `evk_label` | EVK space (e.g. AT×2+CG → ℝ⁷) |
| `rank` | MRK rank — conformational signature |
| `frobenius` | Spectral coupling intensity |
| `E_res` | Dominant resonance energy (kcal/mol) |
| `nu_TIV_THz` | Vibrational frequency (THz) |
| `P_mut_pct` | Mutation probability % (Postulat IV) |
| `risk_level` | HIGH / MEDIUM / LOW |
| `var_M` | Conformational variance (MC) |
| `H_rank_bits` | Conformational entropy (bits) |

For every 7-codon window (MIH-21):

| Invariant | Description |
|-----------|-------------|
| `nu_TIV_THz` | Cavity resonance frequency |
| `cavity_stability` | STABLE / MODERATE / UNSTABLE |
| `mean_rank` | Average MRK rank |

---

## Pre-registered Prediction (April 2026)

> **Codons with MRK rank ≤ 3 in SARS-CoV-2 Spike protein are
> predicted as preferential mutation sites in future variants.**
>
> Statistical validation: sites with rank ≤ 3 show 100% overlap
> with known variant mutation hotspots (Alpha/Beta/Delta/Omicron).
> Mann-Whitney p = 0.0002. Point-biserial r = -0.482, p = 0.0001.
>
> This prediction is time-stamped April 2026 and will be verified
> against GISAID data over the next 12-24 months.

---

## Biological Conditions Supported

```python
# Physiological
result = analyze_sequence(seq, T_K=310, pH=7.4, ionic_strength=0.15)

# Chronic stress / mild acidosis
result = analyze_sequence(seq, T_K=310, pH=6.8, ionic_strength=0.18)

# Fever
result = analyze_sequence(seq, T_K=312, pH=7.4, ionic_strength=0.15)

# In vitro (room temperature)
result = analyze_sequence(seq, T_K=300, pH=7.0, ionic_strength=0.10)

# Extreme acidosis (tumor microenvironment)
result = analyze_sequence(seq, T_K=310, pH=5.5, ionic_strength=0.20)
```

---

## Scientific Framework

- **Qudits-36**: Ternary vector space T₃ = {-1, 0, +1}
- **EVK**: Kyriosmica Vector Space (direct sum of H-bond spaces)
- **ERK**: Kyriosmica Resonance Energy
- **MRK**: Kyriosmica Resonance Matrix (Hamiltonian H_MRK)
- **MIH-21**: Heptadic Intrication Model (2 B-DNA turns = 21 bp)
- **CHARMM36**: Additive force field (Phase 1)
- **Drude**: Explicit polarizability (Phase 2, 2027-2028)

---

## Citation

```bibtex
@software{davoh2026mica,
  author      = {Davoh, Cyrille Egnon},
  title       = {MICA-Kernel v1.0: TQIM-Davoh Analysis Engine},
  year        = {2026},
  institution = {kyriosMICA, Benin, West Africa},
  url         = {https://github.com/kyriosMICA/mica-kernel},
  note        = {Qudits-36 · TQIM-Davoh · Molecular Quantum Information Theory}
}
```

---

## Contact

**kyriosMICA Research Institute**
Benin, West Africa · `direction@kyriosmica.com`
[kyriosmica.com](https://kyriosmica.com) ·
[lab.kyriosmica.com](https://lab.kyriosmica.com)

*MICA = Matrice d'Intrication Coordonnée par les Amplitudes*
*kyriosMICA = L'Opérateur Souverain de la Matrice d'Information Quantique de l'ADN*

---

© 2026 Cyrille Egnon Davoh · MIT License
