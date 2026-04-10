__version__ = "3.0.0"
__author__ = "Cyrille Egnon Davoh"
__institution__ = "kyriosMICA · Bénin"
__framework__ = "TQIM-Davoh v3 · Qudits-36 · Bell/GHZ · MPS · QNN"

from .core import (
    # Constants
    ionic_strength_from_ions,
    # Postulat I
    boltzmann_populations, quantum_amplitudes, boltzmann_bond,
    kbt_kcal, energy_total, delta_E,
    # EVK
    evk_vector, evk_dimension, evk_states_count, evk_composition, evk_label,
    # ERK
    erk_energy, lambda_eff,
    # MRK
    mrk_matrix, mrk_invariants, mrk_matrix_from_vector,
    # Bell/GHZ
    bell_ghz_filter_AT, bell_ghz_filter_CG, build_pair_density,
    # MPS Level 1
    mps_contract_codon,
    # MPS Level 2
    mps_contract_mih21,
    # Postulat IV
    mutation_probability,
    # MIH-21
    mih21_analysis,
    # Monte Carlo
    quantum_collapse_mc,
    # Von Neumann
    von_neumann_entropy,
    # Classification
    classify_conformation,
    # Inference
    infer_from_cavity, sample_from_mih21,
    # T₃-Net QNN
    T3Net, t3net_train_on_sequence, t3net_extract_features,
    # Full analysis
    analyze_sequence,
)
