#!/usr/bin/env python3
# ════════════════════════════════════════════════════════════════════
#  MICA-Kernel v1.0 — Command Line Interface
#  kyriosMICA · © 2026 Cyrille Egnon Davoh
#
#  Usage:
#    python main.py --seq ATGCCCGAT --T 310 --pH 7.4
#    python main.py --fasta sequence.fasta --output results.json
#    python main.py --demo
# ════════════════════════════════════════════════════════════════════

import json
import argparse
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from mica import analyze_sequence, boltzmann_populations, quantum_amplitudes
from mica import mrk_invariants, mutation_probability, mih21_analysis


def print_banner():
    print("""
╔══════════════════════════════════════════════════════════════╗
║   kyrios MICA  — Molecular Information through               ║
║                  Coherent Architecture                       ║
║                                                              ║
║   MICA-Kernel v1.0  ·  TQIM-Davoh  ·  Qudits-36             ║
║   © 2026 Cyrille Egnon Davoh · kyriosMICA · Benin            ║
║                                                              ║
║   "Decoding Life. Encoding the Future."                      ║
╚══════════════════════════════════════════════════════════════╝
""")


def print_summary(result: dict):
    """Pretty print analysis results"""
    meta = result['metadata']
    summ = result['summary']
    cond = summ['conditions']
    amp  = summ['amplitudes']

    print(f"\n{'═'*60}")
    print(f"  ANALYSE TQIM-Davoh — RÉSULTATS")
    print(f"{'═'*60}")
    print(f"\n  SÉQUENCE:")
    print(f"  Codons analysés    : {summ['n_codons']}")
    print(f"  Nucléotides        : {summ['n_nucleotides']}")

    print(f"\n  CONDITIONS:")
    print(f"  Température        : {cond['T_K']} K ({cond['T_C']} °C)")
    print(f"  pH                 : {cond['pH']}")
    print(f"  Force ionique      : {cond['ionic_strength']} mol/L")
    print(f"  Pression           : {cond['pressure_atm']} atm")

    print(f"\n  POSTULAT I — INDÉTERMINATION (amplitudes à {cond['T_K']}K):")
    print(f"  α(-1) = {amp['alpha_minus1']}  [P = {amp['P_minus1']*100:.1f}%]  "
          f"← compression, tunnel")
    print(f"  α( 0) = {amp['alpha_0']}   [P = {amp['P_0']*100:.1f}%]  "
          f"← équilibre dominant")
    print(f"  α(+1) = {amp['alpha_plus1']}  [P = {amp['P_plus1']*100:.1f}%]  "
          f"← extension, ouverture")

    print(f"\n  SIGNATURES MRK (moyenne sur {summ['n_codons']} codons):")
    print(f"  Rang MRK moyen     : {summ['mrk_rank_mean']} ± {summ['mrk_rank_std']}")
    print(f"  ERK moyen          : {summ['erk_mean']} kcal/mol")
    print(f"  E_res moyen        : {summ['E_res_mean']} kcal/mol")

    print(f"\n  POSTULAT IV — MUTABILITÉ (carte de risque):")
    print(f"  Sites HAUT risque  : {summ['high_risk_sites']} "
          f"(rang MRK ≤ 3)")
    print(f"  Sites MOYEN risque : {summ['medium_risk_sites']} "
          f"(rang MRK = 4)")
    print(f"  Sites BAS risque   : {summ['low_risk_sites']} "
          f"(rang MRK ≥ 5)")

    if result['high_risk_sites']:
        print(f"\n  ⚠ SITES À SURVEILLER (rang MRK ≤ 3):")
        for site in result['high_risk_sites']:
            print(f"    Pos {site['pos']:>4} · {site['codon']} · "
                  f"Rang={site['rank']} · E_res={site['E_res']} kcal/mol")

    if result['mih21_windows']:
        print(f"\n  POSTULAT II — COHÉRENCE (fenêtres MIH-21):")
        for w in result['mih21_windows'][:3]:
            print(f"  Fenêtre pos {w['window_start']}-{w['window_end']}:"
                  f" ν_TIV={w['nu_TIV_THz']} THz · "
                  f"stabilité={w['cavity_stability']}")

    print(f"\n{'═'*60}")
    print(f"  {meta['tool']} · {meta['institution']}")
    print(f"  Framework: {meta['framework']}")
    print(f"{'═'*60}\n")


def run_demo():
    """Run demonstration analysis on well-known sequences"""
    print_banner()
    print("  MODE DÉMONSTRATION\n")

    demos = [
        {
            'name': 'BRCA1 N-terminal (6 codons)',
            'seq':  'ATGGATTATCCTTATGAATTTGCTCCT',
            'T_K':  300.0, 'pH': 7.4,
        },
        {
            'name': 'SARS-CoV-2 Spike fragment (D614G site)',
            'seq':  'GATGCTGTTTTAAAAGGTTTGGAGTTTTCTTCTGGG',
            'T_K':  310.0, 'pH': 7.4,
        },
        {
            'name': 'NR3C1 DBD fragment (stress condition)',
            'seq':  'GCAGAGTGTGGTACCTGTAAAGACCTGGCC',
            'T_K':  310.0, 'pH': 6.8,  # mild acidosis under stress
        },
    ]

    for demo in demos:
        print(f"  ── {demo['name']}")
        result = analyze_sequence(
            demo['seq'],
            T_K=demo['T_K'],
            pH=demo['pH'],
            run_mc=False
        )
        print_summary(result)


def main():
    parser = argparse.ArgumentParser(
        description='MICA-Kernel v1.0 — TQIM-Davoh Analysis Engine',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python main.py --demo
  python main.py --seq ATGCCCGAT --T 310 --pH 7.4
  python main.py --seq ATGCCC --output results.json --mc
  python main.py --fasta my_gene.fasta --T 300 --pH 5.0 --ionic 0.15
        """
    )

    parser.add_argument('--seq',    type=str, help='DNA sequence (raw)')
    parser.add_argument('--fasta',  type=str, help='FASTA file path')
    parser.add_argument('--T',      type=float, default=300.0,
                        help='Temperature in Kelvin (default: 300)')
    parser.add_argument('--pH',     type=float, default=7.4,
                        help='pH (default: 7.4)')
    parser.add_argument('--ionic',  type=float, default=0.15,
                        help='Ionic strength mol/L (default: 0.15)')
    parser.add_argument('--pressure', type=float, default=1.0,
                        help='Pressure in atm (default: 1.0)')
    parser.add_argument('--mc',     action='store_true',
                        help='Run Monte Carlo quantum collapse')
    parser.add_argument('--mc-samples', type=int, default=1000,
                        help='MC samples per codon (default: 1000)')
    parser.add_argument('--output', type=str,
                        help='Output JSON file path')
    parser.add_argument('--demo',   action='store_true',
                        help='Run demonstration mode')
    parser.add_argument('--quiet',  action='store_true',
                        help='Suppress banner and verbose output')

    args = parser.parse_args()

    if not args.quiet:
        print_banner()

    if args.demo:
        run_demo()
        return

    # Get sequence
    seq = None
    if args.fasta:
        with open(args.fasta) as f:
            seq = f.read()
    elif args.seq:
        seq = args.seq
    else:
        parser.print_help()
        sys.exit(1)

    # Run analysis
    print(f"  Analyzing sequence ({len(seq)} chars) ...")
    print(f"  Conditions: T={args.T}K · pH={args.pH} · "
          f"I={args.ionic} mol/L")
    if args.mc:
        print(f"  Monte Carlo: {args.mc_samples} samples/codon")
    print()

    result = analyze_sequence(
        seq,
        T_K            = args.T,
        pH             = args.pH,
        ionic_strength = args.ionic,
        pressure_atm   = args.pressure,
        run_mc         = args.mc,
        mc_samples     = args.mc_samples,
    )

    if not args.quiet:
        print_summary(result)

    # Save JSON output
    if args.output:
        with open(args.output, 'w') as f:
            json.dump(result, f, indent=2)
        print(f"  Results saved to: {args.output}")
    else:
        print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
