#!/usr/bin/env python3
"""
01_preprocess.py
================
Pré-traitement de la matrice de counts RNA-seq :
  - Filtrage des gènes à faible expression
  - Statistiques descriptives
  - Visualisation de la distribution

Usage:
    python scripts/01_preprocess.py --input data/counts_matrix.csv --output data/counts_filtered.csv

Auteur: Roméo DJOMAN | AgroParisTech
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path


def load_counts(path: str) -> pd.DataFrame:
    """Charge la matrice de counts."""
    df = pd.read_csv(path, index_col=0)
    print(f"[INFO] Matrice chargée : {df.shape[0]} gènes x {df.shape[1]} échantillons")
    return df


def filter_low_counts(df: pd.DataFrame, min_count: int = 10, min_samples: int = 2) -> pd.DataFrame:
    """
    Filtre les gènes avec moins de `min_count` counts dans au moins `min_samples` échantillons.
    Règle standard HTSeq/DESeq2.
    """
    mask = (df >= min_count).sum(axis=1) >= min_samples
    filtered = df[mask]
    removed = df.shape[0] - filtered.shape[0]
    pct = removed / df.shape[0] * 100
    print(f"[INFO] Gènes filtrés : {removed} ({pct:.1f}%) | Conservés : {filtered.shape[0]}")
    return filtered


def compute_stats(df: pd.DataFrame) -> pd.DataFrame:
    """Calcule des statistiques descriptives par échantillon."""
    stats = pd.DataFrame({
        'total_counts': df.sum(),
        'median_counts': df.median(),
        'n_expressed_genes': (df > 0).sum(),
        'log10_library_size': np.log10(df.sum())
    })
    print("\n[INFO] Statistiques par échantillon :")
    print(stats.to_string())
    return stats


def plot_library_sizes(df: pd.DataFrame, output_dir: Path) -> None:
    """Visualise les tailles de librairies."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle('Contrôle Qualité — Distribution des Counts RNA-seq', fontsize=14, fontweight='bold')

    # Tailles de librairie
    lib_sizes = df.sum() / 1e6
    colors = ['#2196F3' if 'ctrl' in s.lower() else '#E91E63' for s in df.columns]
    axes[0].bar(range(len(lib_sizes)), lib_sizes, color=colors, edgecolor='white', linewidth=0.5)
    axes[0].set_xlabel('Échantillons', fontsize=11)
    axes[0].set_ylabel('Taille librairie (millions de reads)', fontsize=11)
    axes[0].set_title('Tailles des librairies', fontsize=12)
    axes[0].set_xticks(range(len(df.columns)))
    axes[0].set_xticklabels(df.columns, rotation=45, ha='right')
    patch1 = mpatches.Patch(color='#2196F3', label='Contrôle')
    patch2 = mpatches.Patch(color='#E91E63', label='Traité')
    axes[0].legend(handles=[patch1, patch2])

    # Distribution log2 counts
    for i, col in enumerate(df.columns):
        log_counts = np.log2(df[col] + 1)
        axes[1].hist(log_counts, bins=50, alpha=0.4,
                     color='#2196F3' if 'ctrl' in col.lower() else '#E91E63',
                     label=col)
    axes[1].set_xlabel('log2(counts + 1)', fontsize=11)
    axes[1].set_ylabel('Nombre de gènes', fontsize=11)
    axes[1].set_title('Distribution des counts (log2)', fontsize=12)
    axes[1].legend(fontsize=8)

    plt.tight_layout()
    out_path = output_dir / 'qc_library_sizes.png'
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f"[INFO] Figure sauvegardée : {out_path}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Pré-traitement matrice counts RNA-seq')
    parser.add_argument('--input', default='data/counts_matrix.csv', help='Matrice counts bruts')
    parser.add_argument('--output', default='data/counts_filtered.csv', help='Matrice filtrée')
    parser.add_argument('--min_count', type=int, default=10, help='Seuil de counts minimum')
    parser.add_argument('--min_samples', type=int, default=2, help='Nombre minimum d\'échantillons')
    args = parser.parse_args()

    output_dir = Path('results')
    output_dir.mkdir(exist_ok=True)

    df = load_counts(args.input)
    stats = compute_stats(df)
    df_filtered = filter_low_counts(df, args.min_count, args.min_samples)
    df_filtered.to_csv(args.output)
    plot_library_sizes(df_filtered, output_dir)
    stats.to_csv(output_dir / 'sample_stats.csv')
    print(f"\n[SUCCESS] Matrice filtrée exportée : {args.output}")


if __name__ == '__main__':
    main()
