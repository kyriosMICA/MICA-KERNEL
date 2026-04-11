# MICA-Kernel v3.0 — TQIM-Davoh v3

**kyriosMICA Research Institute** · Bénin, West Africa
*Decoding Life. Encoding the Future.*

## API

**Production**: https://api.kyriosmica.com
**Documentation**: https://api.kyriosmica.com/docs

## TQIM-Davoh v3 — Complete Implementation

### 4 Postulats Fondateurs
- **Indétermination** — Per-bond Z_j, |ψ⟩ = √P(-1)|−1⟩ + √P(0)|0⟩ + √P(+1)|+1⟩
- **Cohérence** — MIH-21 cavity (21bp = 2 B-DNA turns), MPS Level 1+2
- **Résonance** — ν_TIV = E_res/h from MRK eigenvalue decomposition
- **Mutabilité** — P(mutation) ∝ P(-1) × P(+1), rank ≤ 3 → HIGH risk

### Physics Engine
- **Bell/GHZ filter**: AT P(|−1,−1⟩)→0 · CG triple extreme→0
- **MPS contraction**: Level 1 (codon) + Level 2 (7-codon MIH-21 cavity)
- **T₃-Net QNN**: Ternary neural network with STE backpropagation
- **Na⁺/K⁺/Mg²⁺**: Individual ion Debye-Hückel (I = 0.5 × Σ cᵢzᵢ²)
- **Von Neumann entropy**: −Tr(ρ log₂ ρ) on cavity density matrix
- **B/A/Stress classification**: Global conformational analysis
- **3-mode inference**: ARGMAX / BOLTZMANN / SAMPLE from ρ_MIH21
- **Epigenetics**: 5mC, 5hmC, m6A, m4C with bond-specific corrections
- **Thermo-classification**: tendu / perturbé / équilibre

### Force Field
Phase 1: CHARMM36 (additive) — Best et al. JCTC 2012, Šponer et al. 2018
Phase 2 (2027-2028): Drude polarisable

### Constants (FLAG_CALIBRATE)
| Constant | Value | Flag |
|----------|-------|------|
| E_PAULI | 0.400 kcal/mol | DFT |
| λ_BASE | 0.12 | NMR |
| α_Mg | 0.08 kcal/mol/mM | DFT+Mg²⁺ |
| α_σ | 0.5 | NMR |
| k_prot | 0.15 kcal/mol | ITC |
| κ_AT | 0.30 | Šponer 2018 |
| κ_CG | 0.50 | Šponer 2018 |

## Endpoints

```
GET  /            Institution info
GET  /health      Health check
POST /analyze     Full sequence analysis
POST /codon       Single codon + Monte Carlo + MPS
POST /mih21       MIH-21 cavity (7 codons)
POST /conditions  Per-bond quantum state
POST /mutability  Mutation risk map
GET  /examples    Example requests
```

## Example

```bash
curl -X POST https://api.kyriosmica.com/analyze \
  -H "Content-Type: application/json" \
  -d '{
    "sequence": "ATGGATTATCCTTATGAATTTGCTCCT",
    "conditions": {
      "temperature_K": 310.15,
      "pH": 7.4,
      "sodium_mM": 140,
      "potassium_mM": 5,
      "magnesium_mM": 2.5
    },
    "topology": {
      "supercoiling_density": -0.06,
      "form": "B-DNA"
    },
    "epigenetics": [
      {"position": 1, "type": "5mC"}
    ],
    "analysis": {
      "force_field": "CHARMM36",
      "run_mc": true,
      "mc_samples": 500
    }
  }'
```

## Validation
- **NR3C1**: 777 codons, ΔRank MRK pH 7.4→6.8 = +2.0 uniform
- **SARS-CoV-2**: 100% rank-3 sites = known mutation hotspots (p=0.0002)
- **Collaboration**: Julien Boblique · SET Theory · DOI: 10.23880/izab-16000673

## Deploy

```bash
# Render (auto-deploy from GitHub)
# Procfile: web: uvicorn api:app --host 0.0.0.0 --port $PORT

# Local
pip install -r requirements.txt
python main.py --demo
uvicorn api:app --host 0.0.0.0 --port 8000
```

## References

1. Watson & Crick (1953). Nature, 171, 737-738.
2. Shannon (1949). Bell Syst. Tech. J., 28(4), 656-715.
3. Best et al. (2012). JCTC, 8, 3257. [CHARMM36]
4. Šponer et al. (2018). Chem.Rev., 118, 4177.
5. Davoh, C.E. (2026). Qudits-36. kyriosMICA. DOI: 10.5281/zenodo.19454825.

---

**kyriosMICA** — Matrice d'Intrication Coordonnée par les Amplitudes
© 2026 Cyrille Egnon Davoh · direction@kyriosmica.com
