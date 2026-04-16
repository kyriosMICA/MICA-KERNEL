# Déploiement MICA-Kernel v3.0

## Prérequis
- Python 3.11+
- NumPy, SciPy, FastAPI, Uvicorn, Pydantic

## Render (Production)

Le repo est connecté à Render via GitHub.
Push sur `main` → Render rebuild automatiquement.

```
Procfile: web: uvicorn api:app --host 0.0.0.0 --port $PORT
Runtime:  python-3.11.8
```

## Vérification post-deploy

```bash
curl https://api.kyriosmica.com/health
# → {"status":"healthy","engine":"MICA-Kernel v3.0","features":"Bell/GHZ·MPS·Na⁺/K⁺/Mg²⁺"}

curl https://api.kyriosmica.com/
# → liste des endpoints
```

## Test complet

```bash
curl -X POST https://api.kyriosmica.com/codon \
  -H "Content-Type: application/json" \
  -d '{"codon":"CGC","conditions":{"temperature_K":310,"pH":7.4,"sodium_mM":140,"potassium_mM":5,"magnesium_mM":2.5},"topology":{"supercoiling_density":-0.06}}'
```

## Fichiers modifiés v1 → v3

| Fichier | Action |
|---------|--------|
| api.py | REMPLACÉ — Pydantic v3, Na⁺/K⁺/Mg²⁺, Bell/GHZ |
| mica/core.py | REMPLACÉ — 5 nouveaux modules |
| mica/__init__.py | REMPLACÉ — exports mis à jour |
| requirements.txt | REMPLACÉ — identique mais vérifié |
| main.py | INCHANGÉ |
| setup.py | INCHANGÉ |
| Procfile | INCHANGÉ |
| runtime.txt | INCHANGÉ |
