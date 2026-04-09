# ════════════════════════════════════════════════════════════════════
#  MICA-Kernel v1.0
#  kyriosMICA Research Institute · Benin, West Africa
#  © 2026 Cyrille Egnon Davoh
#  TQIM-Davoh · Qudits-36
# ════════════════════════════════════════════════════════════════════

__version__     = "1.0.0"
__author__      = "Cyrille Egnon Davoh"
__institution__ = "kyriosMICA · Benin, West Africa"
__framework__   = "TQIM-Davoh · Qudits-36"

from .core import (
    # Thermodynamics — Postulat de l'Indétermination
    boltzmann_populations,
    quantum_amplitudes,
    kbt,
    delta_E,

    # EVK — Espace Vectoriel Kyriosmica
    evk_vector,
    evk_dimension,
    evk_states_count,
    evk_composition,
    evk_label,

    # ERK — Espace de Résonance Kyriosmica
    erk_energy,

    # MRK — Matrice de Résonance Kyriosmica
    mrk_matrix,
    mrk_invariants,

    # Postulat de la Mutabilité
    mutation_probability,

    # MIH-21 — Postulat de la Cohérence
    mih21_analysis,

    # Monte Carlo quantum collapse
    quantum_collapse_mc,

    # Full sequence analysis
    analyze_sequence,
)

__all__ = [
    'boltzmann_populations', 'quantum_amplitudes', 'kbt', 'delta_E',
    'evk_vector', 'evk_dimension', 'evk_states_count',
    'evk_composition', 'evk_label',
    'erk_energy',
    'mrk_matrix', 'mrk_invariants',
    'mutation_probability',
    'mih21_analysis',
    'quantum_collapse_mc',
    'analyze_sequence',
]
