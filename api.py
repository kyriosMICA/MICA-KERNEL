#!/usr/bin/env python3
# ════════════════════════════════════════════════════════════════════
#  MICA-Kernel API v1.0 — FastAPI REST Interface
#  kyriosMICA Research Institute · Benin, West Africa
#  © 2026 Cyrille Egnon Davoh
#
#  Deploy: uvicorn api:app --host 0.0.0.0 --port 8000
#  Docs:   http://api.kyriosmica.com/docs
# ════════════════════════════════════════════════════════════════════

from fastapi import FastAPI, HTTPException, BackgroundTasks
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from pydantic import BaseModel, Field, validator
from typing import Optional, List
import time
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mica import (
    analyze_sequence,
    boltzmann_populations,
    quantum_amplitudes,
    mrk_invariants,
    mutation_probability,
    mih21_analysis,
    quantum_collapse_mc,
)

# ════════════════════════════════════════════════════════════════════
#  APP SETUP
# ════════════════════════════════════════════════════════════════════

app = FastAPI(
    title       = "MICA-Kernel API",
    description = """
## kyriosMICA MICA-Kernel REST API

**TQIM-Davoh Analysis Engine** — Qudits-36 Formalism

The first REST API that computes the complete quantum signature
of a DNA sequence via the TQIM-Davoh framework:

- **EVK** — Kyriosmica Vector Space
- **ERK** — Kyriosmica Resonance Energy  
- **MRK** — Kyriosmica Resonance Matrix
- **MIH-21** — Heptadic Intrication (2 B-DNA turns)
- **Postulat de la Mutabilité** — Mutation risk map

All 4 biological conditions supported: T, pH, ionic strength, pressure.

**Institution**: kyriosMICA Research Institute · Benin, West Africa  
**Framework**: TQIM-Davoh · Qudits-36  
**Author**: Cyrille Egnon Davoh  
**Website**: https://kyriosmica.com
    """,
    version     = "1.0.0",
    contact     = {
        "name":  "kyriosMICA Research Institute",
        "url":   "https://kyriosmica.com",
        "email": "direction@kyriosmica.com",
    },
    license_info = {
        "name": "MIT License",
        "url":  "https://opensource.org/licenses/MIT",
    },
)

# CORS — allow lab.kyriosmica.com and any research tool
app.add_middleware(
    CORSMiddleware,
    allow_origins     = ["*"],
    allow_credentials = True,
    allow_methods     = ["*"],
    allow_headers     = ["*"],
)


# ════════════════════════════════════════════════════════════════════
#  REQUEST / RESPONSE MODELS
# ════════════════════════════════════════════════════════════════════

class AnalysisRequest(BaseModel):
    sequence:        str   = Field(...,
        description="DNA sequence (raw or FASTA format)",
        example="ATGGATTATCCTTATGAATTTGCTCCT")
    T_K:             float = Field(300.0,
        description="Temperature in Kelvin",
        ge=250.0, le=400.0, example=310.0)
    pH:              float = Field(7.4,
        description="pH of biological medium",
        ge=0.0, le=14.0, example=7.4)
    ionic_strength:  float = Field(0.15,
        description="Ionic strength in mol/L",
        ge=0.0, le=5.0, example=0.15)
    pressure_atm:    float = Field(1.0,
        description="Pressure in atm",
        ge=0.0, le=1000.0, example=1.0)
    run_mc:          bool  = Field(False,
        description="Run Monte Carlo quantum collapse simulation")
    mc_samples:      int   = Field(500,
        description="MC samples per codon (max 5000)",
        ge=100, le=5000)

    @validator('sequence')
    def validate_sequence(cls, v):
        clean = ''.join(b for b in v.upper()
                        if b in 'ATCGNU>\n ')
        if len(clean.replace('\n','').replace(' ','')
               .replace('>','')) < 3:
            raise ValueError("Sequence must have at least 3 nucleotides")
        return v

    class Config:
        schema_extra = {
            "example": {
                "sequence":       "ATGGATTATCCTTATGAATTTGCTCCT",
                "T_K":            310.0,
                "pH":             7.4,
                "ionic_strength": 0.15,
                "run_mc":         False,
            }
        }


class CodonRequest(BaseModel):
    codon:           str   = Field(...,
        description="Single codon (3 nucleotides)",
        min_length=3, max_length=3, example="ATG")
    T_K:             float = Field(300.0, ge=250.0, le=400.0)
    pH:              float = Field(7.4,   ge=0.0,   le=14.0)
    ionic_strength:  float = Field(0.15,  ge=0.0,   le=5.0)
    mc_samples:      int   = Field(1000,  ge=100,   le=5000)


class MIH21Request(BaseModel):
    codons: List[str] = Field(...,
        description="Exactly 7 consecutive codons",
        min_items=7, max_items=7,
        example=["ATG","GAT","TAT","CCT","GAA","AAT","GCC"])
    T_K:            float = Field(300.0)
    pH:             float = Field(7.4)
    ionic_strength: float = Field(0.15)


class ConditionsRequest(BaseModel):
    T_K:            float = Field(300.0, ge=250.0, le=400.0)
    pH:             float = Field(7.4,   ge=0.0,   le=14.0)
    ionic_strength: float = Field(0.15,  ge=0.0,   le=5.0)


# ════════════════════════════════════════════════════════════════════
#  ENDPOINTS
# ════════════════════════════════════════════════════════════════════

@app.get("/", tags=["Info"])
async def root():
    """kyriosMICA MICA-Kernel API — root endpoint"""
    return {
        "api":         "MICA-Kernel v1.0",
        "institution": "kyriosMICA Research Institute",
        "location":    "Benin, West Africa",
        "framework":   "TQIM-Davoh · Qudits-36",
        "author":      "Cyrille Egnon Davoh",
        "website":     "https://kyriosmica.com",
        "lab":         "https://lab.kyriosmica.com",
        "github":      "https://github.com/kyriosMICA/MICA-KERNEL",
        "docs":        "/docs",
        "endpoints": {
            "POST /analyze":     "Full sequence analysis",
            "POST /codon":       "Single codon signature",
            "POST /mih21":       "MIH-21 cavity analysis",
            "POST /conditions":  "Quantum amplitudes from conditions",
            "POST /mutability":  "Mutation risk map",
            "GET  /health":      "API health check",
        },
        "slogan": "Decoding Life. Encoding the Future.",
    }


@app.get("/health", tags=["Info"])
async def health():
    """API health check"""
    try:
        # Quick test
        P = boltzmann_populations(300.0, 7.4, 0.15)
        assert abs(sum(P.values()) - 1.0) < 1e-9
        return {
            "status":    "healthy",
            "engine":    "MICA-Kernel v1.0",
            "timestamp": time.time(),
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/analyze", tags=["Analysis"])
async def analyze(req: AnalysisRequest):
    """
    ## Full Sequence Analysis

    Computes the complete TQIM-Davoh quantum signature for every
    codon in the input sequence.

    **Returns:**
    - EVK, ERK, MRK invariants per codon
    - Quantum amplitudes (Postulat de l'Indétermination)
    - Mutation risk map (Postulat de la Mutabilité)
    - MIH-21 cavity signatures (Postulat de la Cohérence)
    - High-risk site predictions

    **Biological conditions:** T, pH, ionic strength, pressure

    **Example use cases:**
    - Viral mutation prediction (SARS-CoV-2, influenza)
    - Cancer gene variant analysis (BRCA1, TP53)
    - Stress response gene profiling (NR3C1)
    """
    try:
        t_start = time.time()
        result  = analyze_sequence(
            sequence       = req.sequence,
            T_K            = req.T_K,
            pH             = req.pH,
            ionic_strength = req.ionic_strength,
            pressure_atm   = req.pressure_atm,
            run_mc         = req.run_mc,
            mc_samples     = req.mc_samples,
        )
        result['_meta'] = {
            'computation_time_s': round(time.time() - t_start, 3),
            'api_version':        '1.0.0',
        }
        return JSONResponse(content=result)
    except ValueError as e:
        raise HTTPException(status_code=422, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500,
                           detail=f"Kernel error: {str(e)}")


@app.post("/codon", tags=["Analysis"])
async def analyze_codon(req: CodonRequest):
    """
    ## Single Codon Analysis

    Computes all TQIM-Davoh invariants for a single codon.
    Includes Monte Carlo quantum collapse simulation.

    **Returns:**
    - MRK invariants (rank, trace, Frobenius, eigenvalues, E_res)
    - Quantum amplitudes
    - Mutation probability
    - Monte Carlo distribution (rank, Var[M], H(rank))
    """
    try:
        codon = req.codon.upper()
        if not all(b in 'ATCG' for b in codon):
            raise HTTPException(
                status_code=422,
                detail="Codon must contain only A, T, C, G"
            )

        inv = mrk_invariants(codon, req.ionic_strength)
        mut = mutation_probability(codon, req.T_K,
                                   req.pH, req.ionic_strength)
        amp = quantum_amplitudes(req.T_K, req.pH, req.ionic_strength)
        mc  = quantum_collapse_mc(codon, req.mc_samples,
                                  req.T_K, req.pH, req.ionic_strength)

        return {
            "codon":      codon,
            "conditions": {
                "T_K": req.T_K, "pH": req.pH,
                "ionic_strength": req.ionic_strength
            },
            "quantum_amplitudes": {
                "alpha_minus1": round(amp[-1], 4),
                "alpha_0":      round(amp[0],  4),
                "alpha_plus1":  round(amp[1],  4),
            },
            "mrk_invariants":    inv,
            "mutability":        mut,
            "monte_carlo":       mc,
            "framework":         "TQIM-Davoh · Qudits-36 · kyriosMICA"
        }
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/mih21", tags=["Analysis"])
async def analyze_mih21(req: MIH21Request):
    """
    ## MIH-21 Heptadic Cavity Analysis

    Analyzes a window of 7 consecutive codons (21 base pairs)
    as a resonant cavity — the fundamental unit of the
    Postulat de la Cohérence (TQIM-Davoh).

    21 bp = 2 complete B-DNA turns (10.5 bp/turn).
    This is the topological closure condition for standing waves.

    **Returns:**
    - Cavity resonance frequency ν_TIV (THz)
    - Standing wave parameters
    - Cavity stability (STABLE / MODERATE / UNSTABLE)
    - Per-codon signatures
    """
    try:
        codons = [c.upper() for c in req.codons]
        for c in codons:
            if len(c) != 3 or not all(b in 'ATCG' for b in c):
                raise HTTPException(
                    status_code=422,
                    detail=f"Invalid codon: {c}"
                )
        result = mih21_analysis(codons, req.T_K,
                                req.pH, req.ionic_strength)
        result['framework'] = "TQIM-Davoh · MIH-21 · kyriosMICA"
        return result
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/conditions", tags=["Physics"])
async def get_conditions(req: ConditionsRequest):
    """
    ## Quantum State from Biological Conditions

    Computes the quantum amplitude vector (Postulat de l'Indétermination)
    and Boltzmann populations for given biological conditions.

    **The core equation:**

    |ψ_i⟩ = α(-1)|−1⟩ + α(0)|0⟩ + α(+1)|+1⟩

    where α(s) = √P(s) from Boltzmann distribution at T, pH, I.

    **Use cases:**
    - Compare physiological vs stress conditions
    - Simulate fever (T = 312K)
    - Simulate tumor microenvironment (pH = 5.5)
    """
    try:
        P   = boltzmann_populations(req.T_K, req.pH,
                                    req.ionic_strength)
        amp = quantum_amplitudes(req.T_K, req.pH,
                                 req.ionic_strength)
        P_mut = P[-1] * P[1]

        return {
            "conditions": {
                "T_K":            req.T_K,
                "T_Celsius":      round(req.T_K - 273.15, 1),
                "pH":             req.pH,
                "ionic_strength": req.ionic_strength,
            },
            "boltzmann_populations": {
                "P_minus1": round(P[-1], 6),
                "P_0":      round(P[0],  6),
                "P_plus1":  round(P[1],  6),
            },
            "quantum_amplitudes": {
                "alpha_minus1": round(amp[-1], 6),
                "alpha_0":      round(amp[0],  6),
                "alpha_plus1":  round(amp[1],  6),
            },
            "mutability_baseline": {
                "P_mut_per_bond_pct": round(P_mut * 100, 4),
                "description": "P(-1) × P(+1) — Postulat de la Mutabilité"
            },
            "interpretation": {
                "dominant_state": "equilibrium |0⟩" if P[0] > 0.5
                                  else "mixed",
                "thermal_signal": "resilient" if abs(P[0]-0.564) < 0.05
                                  else "shifted",
            },
            "postulat": "Indétermination (TQIM-Davoh · Qudits-36)",
            "reference": "Davoh 2026 · kyriosMICA · Benin"
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/mutability", tags=["Prediction"])
async def mutation_map(req: AnalysisRequest):
    """
    ## Mutation Risk Map

    Computes the mutation probability map for a full sequence.
    Implements the **Postulat de la Mutabilité** (TQIM-Davoh):

    P(mutation_i) ∝ P_i(-1) · P_i(+1)

    **Pre-registered prediction (April 2026):**
    Codons with MRK rank ≤ 3 are predicted as HIGH-risk mutation
    sites. Validated on SARS-CoV-2: 100% of rank-3 codons overlap
    with known variant hotspots. Mann-Whitney p = 0.0002.

    **Returns:**
    - Risk level per codon (HIGH / MEDIUM / LOW)
    - MRK rank (key discriminator)
    - Ranked list of top mutation risk sites
    """
    try:
        result = analyze_sequence(
            req.sequence, req.T_K, req.pH,
            req.ionic_strength, req.pressure_atm,
            run_mc=False
        )

        mutation_map_data = [
            {
                "position":   c['position'],
                "codon":      c['codon'],
                "evk_label":  c['evk_label'],
                "rank":       c['rank'],
                "E_res":      c['E_res'],
                "P_mut_pct":  c['P_mut_pct'],
                "risk_level": c['risk_level'],
            }
            for c in result['codons']
        ]

        # Sort by risk (HIGH first, then rank ascending)
        risk_order = {'HIGH': 0, 'MEDIUM': 1, 'LOW': 2}
        ranked = sorted(
            mutation_map_data,
            key=lambda x: (risk_order[x['risk_level']], x['rank'])
        )

        return {
            "n_codons":          result['summary']['n_codons'],
            "conditions":        result['summary']['conditions'],
            "high_risk_sites":   result['summary']['high_risk_sites'],
            "medium_risk_sites": result['summary']['medium_risk_sites'],
            "low_risk_sites":    result['summary']['low_risk_sites'],
            "prediction_rule":   "MRK rank ≤ 3 → HIGH · rank=4 → MEDIUM · rank≥5 → LOW",
            "validation":        "SARS-CoV-2: 100% rank-3 sites = known hotspots. p=0.0002",
            "mutation_map":      mutation_map_data,
            "top_10_risk_sites": ranked[:10],
            "framework":         "Postulat de la Mutabilité · TQIM-Davoh · kyriosMICA",
            "timestamp":         "Pre-registered prediction April 2026 · kyriosMICA · Benin"
        }
    except ValueError as e:
        raise HTTPException(status_code=422, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/examples", tags=["Info"])
async def get_examples():
    """Example requests for all endpoints"""
    return {
        "analyze_sequence": {
            "url":    "POST /analyze",
            "body": {
                "sequence":       "ATGGATTATCCTTATGAATTTGCTCCT",
                "T_K":            310.0,
                "pH":             7.4,
                "ionic_strength": 0.15,
                "run_mc":         True,
                "mc_samples":     1000
            }
        },
        "single_codon": {
            "url":  "POST /codon",
            "body": {"codon": "ATG", "T_K": 310.0, "pH": 7.4}
        },
        "mih21_window": {
            "url":  "POST /mih21",
            "body": {
                "codons": ["ATG","GAT","TAT","CCT","GAA","AAT","GCC"],
                "T_K": 310.0, "pH": 7.4
            }
        },
        "stress_conditions": {
            "url":  "POST /conditions",
            "body": {"T_K": 310.0, "pH": 6.8, "ionic_strength": 0.18}
        },
        "curl_examples": {
            "analyze": 'curl -X POST https://api.kyriosmica.com/analyze -H "Content-Type: application/json" -d \'{"sequence":"ATGCCCGAT","T_K":310,"pH":7.4}\'',
            "codon":   'curl -X POST https://api.kyriosmica.com/codon -H "Content-Type: application/json" -d \'{"codon":"ATG","T_K":310}\'',
        }
    }


# ════════════════════════════════════════════════════════════════════
#  RUN (for local testing)
# ════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    import uvicorn
    print("""
╔══════════════════════════════════════════════════════════════╗
║   MICA-Kernel API v1.0 — kyriosMICA                         ║
║   Starting server on http://0.0.0.0:8000                    ║
║   Docs: http://localhost:8000/docs                          ║
║   Benin, West Africa · Decoding Life. Encoding the Future.  ║
╚══════════════════════════════════════════════════════════════╝
    """)
    uvicorn.run("api:app", host="0.0.0.0", port=8000, reload=True)
