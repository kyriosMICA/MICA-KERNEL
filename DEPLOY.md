# MICA-Kernel API — Guide de Déploiement
# api.kyriosmica.com
# kyriosMICA · © 2026 Cyrille Egnon Davoh

# ════════════════════════════════════════
# OPTION 1 — Railway.app (GRATUIT, 5 min)
# ════════════════════════════════════════

# 1. Créer un compte sur railway.app
# 2. New Project → Deploy from GitHub repo
# 3. Sélectionner kyriosMICA/MICA-KERNEL
# 4. Railway détecte automatiquement Python
# 5. Variables d'environnement : aucune requise
# 6. Domain → Custom domain → api.kyriosmica.com

# ════════════════════════════════════════
# OPTION 2 — Render.com (GRATUIT, 10 min)
# ════════════════════════════════════════

# 1. render.com → New Web Service
# 2. Connect GitHub → kyriosMICA/MICA-KERNEL
# 3. Build Command : pip install -r requirements.txt
# 4. Start Command : uvicorn api:app --host 0.0.0.0 --port $PORT
# 5. Free tier disponible
# 6. Custom domain → api.kyriosmica.com

# ════════════════════════════════════════
# OPTION 3 — Local (test immédiat)
# ════════════════════════════════════════

# pip install fastapi uvicorn numpy scipy
# uvicorn api:app --host 0.0.0.0 --port 8000
# → http://localhost:8000
# → http://localhost:8000/docs  (documentation interactive)

# ════════════════════════════════════════
# TESTS RAPIDES APRÈS DÉPLOIEMENT
# ════════════════════════════════════════

# Health check
# curl https://api.kyriosmica.com/health

# Analyser une séquence
# curl -X POST https://api.kyriosmica.com/analyze \
#   -H "Content-Type: application/json" \
#   -d '{"sequence":"ATGCCCGAT","T_K":310,"pH":7.4}'

# Codon unique
# curl -X POST https://api.kyriosmica.com/codon \
#   -H "Content-Type: application/json" \
#   -d '{"codon":"ATG","T_K":310}'

# Conditions stress
# curl -X POST https://api.kyriosmica.com/conditions \
#   -H "Content-Type: application/json" \
#   -d '{"T_K":310,"pH":6.0,"ionic_strength":0.18}'
