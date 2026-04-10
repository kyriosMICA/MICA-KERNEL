#!/usr/bin/env python3
# MICA-Kernel API v3.0 — kyriosMICA · TQIM-Davoh v3 · Bell/GHZ · MPS
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field, validator
from typing import Optional, List
import time, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mica import (analyze_sequence, boltzmann_populations, quantum_amplitudes, boltzmann_bond, mrk_invariants, mutation_probability, mih21_analysis, quantum_collapse_mc, erk_energy, evk_vector, evk_dimension, evk_states_count, evk_label, ionic_strength_from_ions, lambda_eff, build_pair_density, mps_contract_codon)
from mica.core import N_BONDS

app = FastAPI(title="MICA-Kernel API", description="kyriosMICA v3.0 · TQIM-Davoh v3 · Bell/GHZ · MPS · Na⁺/K⁺/Mg²⁺", version="3.0.0", contact={"name":"kyriosMICA","url":"https://kyriosmica.com","email":"direction@kyriosmica.com"})
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_credentials=True, allow_methods=["*"], allow_headers=["*"])

class EpigeneticMod(BaseModel):
    position: int; base: str = "C"; type: str = "5mC"; strand: str = "sense"
class Conditions(BaseModel):
    temperature_K: float = Field(310.15, ge=250, le=400)
    pH: float = Field(7.4, ge=0, le=14)
    sodium_mM: float = Field(140.0, ge=0, le=1000)
    potassium_mM: float = Field(5.0, ge=0, le=500)
    magnesium_mM: float = Field(0.0, ge=0, le=100)
    pressure_atm: float = Field(1.0, ge=0, le=1000)
class Topology(BaseModel):
    supercoiling_density: float = Field(0.0, ge=-0.15, le=0.10)
    form: str = "B-DNA"
class AnalysisConfig(BaseModel):
    force_field: str = "CHARMM36"; run_mc: bool = False; mc_samples: int = Field(500, ge=100, le=5000)
class AnalyzeRequest(BaseModel):
    sequence: str = Field(...); conditions: Optional[Conditions] = None; topology: Optional[Topology] = None; epigenetics: Optional[List[EpigeneticMod]] = None; analysis: Optional[AnalysisConfig] = None
class CodonRequest(BaseModel):
    codon: str = Field(..., min_length=3, max_length=3); conditions: Optional[Conditions] = None; topology: Optional[Topology] = None; epigenetics: Optional[List[EpigeneticMod]] = None; mc_samples: int = Field(1000, ge=100, le=5000)
class MIH21Request(BaseModel):
    codons: List[str] = Field(..., min_items=7, max_items=7); conditions: Optional[Conditions] = None; topology: Optional[Topology] = None
class ConditionsRequest(BaseModel):
    conditions: Conditions; base: str = "A"; bond_index: int = Field(0, ge=0, le=2); methylation_type: Optional[str] = None

def _c(c): return c or Conditions()
def _t(t): return t or Topology()
def _a(a): return a or AnalysisConfig()
def _m(e): return [{'position':m.position,'base':m.base,'type':m.type,'strand':m.strand} for m in e] if e else None
def _I(c): return ionic_strength_from_ions(c.sodium_mM, c.potassium_mM, c.magnesium_mM)

@app.get("/")
async def root():
    return {"api":"MICA-Kernel v3.0","spec":"TQIM-Davoh v3","institution":"kyriosMICA","features":["Bell/GHZ","MPS","Na⁺/K⁺/Mg²⁺ Debye-Hückel","per-bond Z_j","5mC/5hmC/m6A/m4C"],"endpoints":{"POST /analyze":"Full","POST /codon":"Single+MC","POST /mih21":"MIH-21","POST /conditions":"Per-bond","POST /mutability":"Risk map"},"slogan":"Decoding Life. Encoding the Future."}

@app.get("/health")
async def health():
    P = boltzmann_populations(300, 7.4, 0.15); assert abs(sum(P.values())-1)<1e-9
    return {"status":"healthy","engine":"MICA-Kernel v3.0","features":"Bell/GHZ·MPS·Na⁺/K⁺/Mg²⁺"}

@app.post("/analyze")
async def analyze(req: AnalyzeRequest):
    try:
        c,t,a,m = _c(req.conditions),_t(req.topology),_a(req.analysis),_m(req.epigenetics)
        t0=time.time()
        r = analyze_sequence(req.sequence, T_K=c.temperature_K, pH=c.pH, sodium_mM=c.sodium_mM, potassium_mM=c.potassium_mM, magnesium_mM=c.magnesium_mM, pressure_atm=c.pressure_atm, supercoiling=t.supercoiling_density, methylations=m, run_mc=a.run_mc, mc_samples=a.mc_samples, force_field=a.force_field)
        r['compute_time_s'] = round(time.time()-t0,3); return r
    except ValueError as e: raise HTTPException(422,str(e))
    except Exception as e: raise HTTPException(500,str(e))

@app.post("/codon")
async def single_codon(req: CodonRequest):
    try:
        codon=req.codon.upper(); c,t=_c(req.conditions),_t(req.topology); m=_m(req.epigenetics)
        I=_I(c); lam=lambda_eff(I,t.supercoiling_density)
        inv=mrk_invariants(codon,I,t.supercoiling_density)
        erk=erk_energy(codon,c.temperature_K,c.pH,I,'quantum',c.magnesium_mM,[{'position':mm.position,'type':mm.type} for mm in req.epigenetics] if req.epigenetics else None)
        mut=mutation_probability(codon,c.temperature_K,c.pH,I,t.supercoiling_density,c.magnesium_mM)
        mc=quantum_collapse_mc(codon,req.mc_samples,c.temperature_K,c.pH,I,t.supercoiling_density,c.magnesium_mM)
        mps=mps_contract_codon(list(codon),c.temperature_K,c.pH,c.magnesium_mM,m,lam)
        bonds=[]
        for bi,base in enumerate(codon):
            for j in range(N_BONDS.get(base,2)):
                bb=boltzmann_bond(base,j,c.temperature_K,c.pH,c.magnesium_mM)
                bonds.append({'base':base,'bond_index':j,'P_minus1':round(bb['probs'][-1],6),'P_0':round(bb['probs'][0],6),'P_plus1':round(bb['probs'][1],6),'dominant':bb['dominant'],'Z':round(bb['Z'],6)})
        from mica.core import _codon_to_aa
        return {"codon":codon,"amino_acid":_codon_to_aa(codon),"conditions":{"I_mol_L":round(I,4),"lambda_eff":round(lam,6),"sodium_mM":c.sodium_mM,"potassium_mM":c.potassium_mM,"magnesium_mM":c.magnesium_mM},"evk":{"label":evk_label(codon),"dim":evk_dimension(codon),"states":evk_states_count(codon),"vector":evk_vector(codon)},"erk":erk,"mrk_invariants":inv,"per_bond_boltzmann":bonds,"mutability":mut,"monte_carlo":mc,"mps":mps,"framework":"TQIM-Davoh v3 · Bell/GHZ · MPS · kyriosMICA"}
    except Exception as e: raise HTTPException(500,str(e))

@app.post("/mih21")
async def analyze_mih21(req: MIH21Request):
    try:
        codons=[co.upper() for co in req.codons]; c,t=_c(req.conditions),_t(req.topology); I=_I(c)
        r=mih21_analysis(codons,c.temperature_K,c.pH,I,t.supercoiling_density,c.magnesium_mM)
        r['conditions']={'I_mol_L':round(I,4),'sodium_mM':c.sodium_mM,'potassium_mM':c.potassium_mM,'magnesium_mM':c.magnesium_mM}
        return r
    except Exception as e: raise HTTPException(500,str(e))

@app.post("/conditions")
async def get_conditions(req: ConditionsRequest):
    try:
        c=req.conditions; I=_I(c)
        bb=boltzmann_bond(req.base,req.bond_index,c.temperature_K,c.pH,c.magnesium_mM,req.methylation_type)
        P,amp=bb['probs'],bb['amplitudes']
        return {"conditions":{**c.dict(),"I_mol_L":round(I,4)},"base":req.base,"bond_index":req.bond_index,"methylation_type":req.methylation_type,"boltzmann":{"P_minus1":round(P[-1],6),"P_0":round(P[0],6),"P_plus1":round(P[1],6)},"amplitudes":{"alpha_minus1":round(amp[-1],6),"alpha_0":round(amp[0],6),"alpha_plus1":round(amp[1],6)},"Z":round(bb['Z'],6),"dominant":bb['dominant'],"purity":round(bb['purity'],6)}
    except Exception as e: raise HTTPException(500,str(e))

@app.post("/mutability")
async def mutation_map(req: AnalyzeRequest):
    try:
        c,t,m=_c(req.conditions),_t(req.topology),_m(req.epigenetics)
        r=analyze_sequence(req.sequence,c.temperature_K,c.pH,sodium_mM=c.sodium_mM,potassium_mM=c.potassium_mM,magnesium_mM=c.magnesium_mM,supercoiling=t.supercoiling_density,methylations=m)
        mm=[{'position':co['position'],'codon':co['codon'],'amino_acid':co.get('amino_acid','?'),'rank':co['rank'],'risk_level':co['risk_level'],'P_mut_pct':co['P_mut_pct'],'mps_entropy':co.get('mps_entropy',0)} for co in r['codons']]
        return {"n_codons":r['summary']['n_codons'],"conditions":r['summary']['conditions'],"high_risk_sites":r['summary']['high_risk_sites'],"mutation_map":mm,"framework":"TQIM-Davoh v3 · Bell/GHZ · MPS"}
    except Exception as e: raise HTTPException(500,str(e))

@app.get("/examples")
async def examples():
    return {"v3":{"url":"POST /analyze","body":{"sequence":"ATGGATTATCCTTATGAATTTGCTCCT","conditions":{"temperature_K":310.15,"pH":7.4,"sodium_mM":140,"potassium_mM":5,"magnesium_mM":2.5},"topology":{"supercoiling_density":-0.06},"epigenetics":[{"position":1,"type":"5mC"}],"analysis":{"run_mc":True,"mc_samples":500}}},"note":"Na⁺/K⁺/Mg²⁺ → I via Debye-Hückel. Bell/GHZ + MPS automatic."}

if __name__ == "__main__":
    import uvicorn; uvicorn.run("api:app", host="0.0.0.0", port=8000, reload=True)
