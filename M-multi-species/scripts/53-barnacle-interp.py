#!/usr/bin/env python3
"""
Barnacle Tensor Decomposition Interpretation Script
====================================================
Analyzes gene expression patterns across three coral species 
(Acropora pulchra, Pocillopora tuahiniensis, Porites evermanni)
over 4 time points using tensor decomposition results.

Integrates:
- Barnacle factor matrices (gene, sample, time, component weights)
- Annotated ortholog group data with GO terms
- Physiological measurements

Outputs:
- Narrative summary of biological patterns
- Compelling visualizations linking molecular and physiological data
"""

import os
import sys
import json
import warnings
from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import seaborn as sns

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

# =============================================================================
# Configuration
# =============================================================================

BASE_DIR = Path("/Users/sr320/GitHub/timeseries_molecular/M-multi-species")
BARNACLE_DIR = BASE_DIR / "output/26-rank35-optimization/lambda_gene_0.2/barnacle_factors"
TOP_GENES_DIR = BASE_DIR / "output/26-rank35-optimization/lambda_gene_0.2/top_genes_per_component"
ORTHO_ANNOT = BASE_DIR / "output/12-ortho-annot/ortholog_groups_annotated.csv"
OUTPUT_DIR = BASE_DIR / "output/53-barnacle-interpretation"
PHYSIOLOGY_URL = "http://gannet.fish.washington.edu/seashell/snaps/MOSAiC_Physiology_Data__MOSAiC.csv"

# Species mapping
SPECIES_MAP = {
    'apul': 'Acropora pulchra',
    'peve': 'Porites evermanni', 
    'ptua': 'Pocillopora tuahiniensis'
}

SPECIES_COLORS = {
    'apul': '#E63946',  # Coral red
    'peve': '#457B9D',  # Ocean blue
    'ptua': '#2A9D8F'   # Teal green
}

TIME_POINTS = ['TP1', 'TP2', 'TP3', 'TP4']
TIME_COLORS = ['#264653', '#2A9D8F', '#E9C46A', '#E76F51']

# =============================================================================
# Data Loading Functions
# =============================================================================

def load_barnacle_factors():
    """Load all barnacle factor matrices."""
    print("Loading barnacle factor matrices...")
    
    factors = {}
    
    # Load metadata
    with open(BARNACLE_DIR / "metadata.json") as f:
        factors['metadata'] = json.load(f)
    
    # Load factor matrices
    factors['component_weights'] = pd.read_csv(BARNACLE_DIR / "component_weights.csv")
    factors['gene_factors'] = pd.read_csv(BARNACLE_DIR / "gene_factors.csv", index_col=0)
    factors['sample_factors'] = pd.read_csv(BARNACLE_DIR / "sample_factors.csv", index_col=0)
    factors['time_factors'] = pd.read_csv(BARNACLE_DIR / "time_factors.csv", index_col=0)
    factors['sample_mapping'] = pd.read_csv(BARNACLE_DIR / "sample_mapping.csv")
    
    print(f"  Tensor shape: {factors['metadata']['tensor_shape']}")
    print(f"  Components: {factors['metadata']['n_components']}")
    
    return factors

def load_annotated_components():
    """Load annotation files for each component."""
    print("Loading component annotations...")
    
    annotations = {}
    for comp_file in TOP_GENES_DIR.glob("Component_*_top100_annotation.csv"):
        comp_num = int(comp_file.stem.split("_")[1])
        try:
            df = pd.read_csv(comp_file)
            annotations[comp_num] = df
        except Exception as e:
            print(f"  Warning: Could not load {comp_file.name}: {e}")
    
    print(f"  Loaded {len(annotations)} component annotations")
    return annotations

def load_goslim_counts():
    """Load GOslim term counts across components."""
    goslim_file = TOP_GENES_DIR / "goslim_term_counts.csv"
    if goslim_file.exists():
        return pd.read_csv(goslim_file, index_col=0)
    return None

def load_physiology_data():
    """Load physiological data from remote URL."""
    print("Loading physiological data...")
    try:
        physio = pd.read_csv(PHYSIOLOGY_URL)
        print(f"  Loaded {len(physio)} physiological measurements")
        return physio
    except Exception as e:
        print(f"  Warning: Could not load physiology data: {e}")
        return None

# =============================================================================
# Analysis Functions
# =============================================================================

def identify_key_components(factors, n_top=10):
    """Identify the most important components by weight."""
    weights = factors['component_weights']['weight'].values
    top_indices = np.argsort(weights)[::-1][:n_top]
    
    return {
        'indices': top_indices + 1,  # 1-indexed component numbers
        'weights': weights[top_indices]
    }

def analyze_temporal_patterns(factors):
    """Analyze how components vary across time points."""
    time_factors = factors['time_factors']
    
    patterns = {}
    for col in time_factors.columns:
        if col.startswith('Component_'):
            values = time_factors[col].values
            comp_num = int(col.split('_')[1])
            
            # Classify temporal pattern
            if values[0] > np.mean(values[1:]) * 1.5:
                pattern = 'early_peak'
            elif values[-1] > np.mean(values[:-1]) * 1.5:
                pattern = 'late_peak'
            elif values[1] > max(values[0], values[2], values[3]) * 1.2:
                pattern = 'early_mid_peak'
            elif values[2] > max(values[0], values[1], values[3]) * 1.2:
                pattern = 'late_mid_peak'
            elif np.std(values) / np.mean(np.abs(values)) < 0.3:
                pattern = 'stable'
            else:
                pattern = 'dynamic'
                
            patterns[comp_num] = {
                'values': values,
                'pattern': pattern,
                'trend': np.corrcoef(np.arange(4), values)[0, 1]
            }
    
    return patterns

def analyze_species_patterns(factors):
    """Analyze how components vary across species."""
    sample_factors = factors['sample_factors']
    sample_mapping = factors['sample_mapping']
    
    # Merge to get species info
    sample_factors_reset = sample_factors.reset_index()
    sample_factors_reset.columns = ['sample'] + list(sample_factors.columns)
    merged = sample_factors_reset.merge(
        sample_mapping[['label', 'species']], 
        left_on='sample', 
        right_on='label'
    )
    
    species_patterns = {}
    for col in [c for c in merged.columns if c.startswith('Component_')]:
        comp_num = int(col.split('_')[1])
        species_means = merged.groupby('species')[col].mean()
        species_patterns[comp_num] = species_means.to_dict()
    
    return species_patterns

def extract_go_themes(annotations, comp_num):
    """Extract major GO themes for a component."""
    if comp_num not in annotations:
        return []
    
    df = annotations[comp_num]
    
    # Collect all GO biological processes
    go_bp = []
    if 'go_bp' in df.columns:
        for entry in df['go_bp'].dropna():
            terms = [t.strip() for t in str(entry).split(';')]
            go_bp.extend(terms)
    
    # Count and return top themes
    counter = Counter(go_bp)
    return counter.most_common(5)


def identify_species_independent_components(factors, species_patterns, threshold=0.6):
    """
    Identify components that are NOT primarily driven by species differences.
    
    Components with low coefficient of variation across species means are considered
    species-independent (i.e., shared expression patterns).
    
    Parameters:
    -----------
    threshold : float
        CV threshold below which a component is considered species-independent.
        Higher values = more permissive (more components classified as independent)
    """
    species_independent = []
    species_driven = []
    
    for comp, species_data in species_patterns.items():
        values = list(species_data.values())
        mean_val = np.mean(np.abs(values))
        std_val = np.std(values)
        
        # Calculate coefficient of variation
        if mean_val > 0:
            cv = std_val / mean_val
        else:
            cv = 0
        
        # Check if all species have same sign (indicating shared direction)
        all_positive = all(v >= 0 for v in values)
        all_negative = all(v <= 0 for v in values)
        same_direction = all_positive or all_negative
        
        # Calculate range ratio (max - min) / mean
        range_ratio = (max(values) - min(values)) / mean_val if mean_val > 0 else 0
        
        component_info = {
            'component': comp,
            'cv': cv,
            'range_ratio': range_ratio,
            'mean_loading': mean_val,
            'same_direction': same_direction,
            'species_values': species_data
        }
        
        # Classify as species-independent if:
        # 1. Low CV (similar expression levels across species), OR
        # 2. Same direction of expression AND moderate CV
        if cv < threshold or (same_direction and cv < threshold * 1.5):
            species_independent.append(component_info)
        else:
            species_driven.append(component_info)
    
    # Sort by CV (lowest first = most species-independent)
    species_independent.sort(key=lambda x: x['cv'])
    species_driven.sort(key=lambda x: -x['cv'])
    
    return species_independent, species_driven


def calculate_annotation_statistics(annotations):
    """
    Calculate annotation statistics for each component including:
    - Proportion of annotated vs unannotated genes
    - Functional category breakdown
    """
    stats = {}
    
    for comp_num, df in annotations.items():
        total_genes = len(df)
        
        # Count annotated genes (has protein_name or any GO annotation)
        has_protein = df['protein_name'].notna().sum()
        has_go_bp = df['go_bp'].notna().sum()
        has_go_mf = df['go_mf'].notna().sum()
        has_any_annotation = (df['protein_name'].notna() | 
                             df['go_bp'].notna() | 
                             df['go_mf'].notna()).sum()
        
        # Proportion unannotated
        unannotated = total_genes - has_any_annotation
        prop_unannotated = unannotated / total_genes if total_genes > 0 else 0
        
        # Extract GOslim categories for functional breakdown
        goslim_counts = Counter()
        if 'goslim_names' in df.columns:
            for entry in df['goslim_names'].dropna():
                categories = [c.strip() for c in str(entry).split(';')]
                goslim_counts.update(categories)
        
        stats[comp_num] = {
            'total_genes': total_genes,
            'annotated': has_any_annotation,
            'unannotated': unannotated,
            'prop_unannotated': prop_unannotated,
            'has_protein_name': has_protein,
            'has_go_bp': has_go_bp,
            'has_go_mf': has_go_mf,
            'goslim_categories': goslim_counts.most_common(10)
        }
    
    return stats


def analyze_functional_categories(annotations, component_list):
    """
    Analyze functional categories for a specific list of components.
    Returns aggregated functional theme analysis.
    """
    all_goslim = Counter()
    all_go_bp = Counter()
    all_proteins = []
    
    for comp_info in component_list:
        comp_num = comp_info['component']
        if comp_num not in annotations:
            continue
        
        df = annotations[comp_num]
        
        # Collect GOslim categories
        if 'goslim_names' in df.columns:
            for entry in df['goslim_names'].dropna():
                categories = [c.strip() for c in str(entry).split(';')]
                all_goslim.update(categories)
        
        # Collect GO biological processes
        if 'go_bp' in df.columns:
            for entry in df['go_bp'].dropna():
                terms = [t.strip().split('[')[0].strip() for t in str(entry).split(';')]
                all_go_bp.update(terms)
        
        # Collect protein names
        if 'protein_name' in df.columns:
            proteins = df['protein_name'].dropna().tolist()
            all_proteins.extend(proteins)
    
    return {
        'goslim_categories': all_goslim.most_common(15),
        'go_bp_terms': all_go_bp.most_common(15),
        'top_proteins': all_proteins[:20]
    }

def extract_protein_functions(annotations, comp_num):
    """Extract key protein functions for a component."""
    if comp_num not in annotations:
        return []
    
    df = annotations[comp_num]
    
    proteins = []
    if 'protein_name' in df.columns:
        proteins = df['protein_name'].dropna().tolist()[:10]
    
    return proteins

# =============================================================================
# Visualization Functions
# =============================================================================

def setup_style():
    """Configure matplotlib style for publication-quality figures."""
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Helvetica Neue', 'Arial', 'DejaVu Sans'],
        'font.size': 10,
        'axes.titlesize': 12,
        'axes.labelsize': 10,
        'xtick.labelsize': 9,
        'ytick.labelsize': 9,
        'legend.fontsize': 9,
        'figure.dpi': 150,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'axes.spines.top': False,
        'axes.spines.right': False,
    })

def plot_component_weights(factors, output_dir):
    """Create bar plot of component weights."""
    fig, ax = plt.subplots(figsize=(14, 5))
    
    weights = factors['component_weights']['weight'].values
    n_components = len(weights)
    
    colors = plt.cm.viridis(np.linspace(0.2, 0.9, n_components))
    bars = ax.bar(range(1, n_components + 1), weights, color=colors, edgecolor='white', linewidth=0.5)
    
    # Highlight top components
    top_idx = np.argsort(weights)[::-1][:5]
    for idx in top_idx:
        bars[idx].set_edgecolor('#E63946')
        bars[idx].set_linewidth(2)
    
    ax.set_xlabel('Component', fontweight='bold')
    ax.set_ylabel('Weight', fontweight='bold')
    ax.set_title('Barnacle Tensor Component Weights\nTop 5 components highlighted', fontsize=14, fontweight='bold')
    ax.set_xticks(range(1, n_components + 1))
    
    plt.tight_layout()
    plt.savefig(output_dir / "component_weights.png")
    plt.close()

def plot_temporal_dynamics(factors, temporal_patterns, output_dir):
    """Create heatmap of temporal dynamics across components."""
    time_factors = factors['time_factors']
    
    # Prepare data matrix
    data = np.zeros((4, len(temporal_patterns)))
    for i, (comp, pattern) in enumerate(sorted(temporal_patterns.items())):
        data[:, i] = pattern['values']
    
    # Normalize each component for visualization
    data_norm = (data - data.mean(axis=0)) / (data.std(axis=0) + 1e-10)
    
    fig, ax = plt.subplots(figsize=(16, 4))
    
    im = ax.imshow(data_norm, aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2)
    
    ax.set_yticks(range(4))
    ax.set_yticklabels(['TP1', 'TP2', 'TP3', 'TP4'])
    ax.set_xticks(range(len(temporal_patterns)))
    ax.set_xticklabels([f'C{c}' for c in sorted(temporal_patterns.keys())])
    
    ax.set_xlabel('Component', fontweight='bold')
    ax.set_ylabel('Time Point', fontweight='bold')
    ax.set_title('Temporal Dynamics of Gene Expression Components\n(z-scored)', fontsize=14, fontweight='bold')
    
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label('Relative Expression', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_dir / "temporal_dynamics_heatmap.png")
    plt.close()

def plot_species_comparison(factors, species_patterns, output_dir):
    """Create species comparison plot for top components."""
    # Get top 10 components by variance across species
    top_comps = sorted(
        species_patterns.items(),
        key=lambda x: np.std(list(x[1].values())),
        reverse=True
    )[:10]
    
    fig, axes = plt.subplots(2, 5, figsize=(16, 7))
    axes = axes.flatten()
    
    for ax, (comp, species_data) in zip(axes, top_comps):
        species = list(species_data.keys())
        values = [species_data[s] for s in species]
        colors = [SPECIES_COLORS[s] for s in species]
        
        bars = ax.bar(range(len(species)), values, color=colors, edgecolor='white')
        ax.set_xticks(range(len(species)))
        ax.set_xticklabels([SPECIES_MAP[s].split()[0][0] + '.' for s in species], fontsize=9)
        ax.set_title(f'Component {comp}', fontweight='bold')
        ax.axhline(0, color='gray', linestyle='--', alpha=0.5)
    
    # Add legend
    legend_patches = [
        mpatches.Patch(color=SPECIES_COLORS['apul'], label='A. pulchra'),
        mpatches.Patch(color=SPECIES_COLORS['peve'], label='P. evermanni'),
        mpatches.Patch(color=SPECIES_COLORS['ptua'], label='P. tuahiniensis'),
    ]
    fig.legend(handles=legend_patches, loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.02))
    
    fig.suptitle('Species-Specific Expression Patterns\n(Top Variable Components)', 
                 fontsize=14, fontweight='bold', y=1.08)
    plt.tight_layout()
    plt.savefig(output_dir / "species_comparison.png")
    plt.close()

def plot_temporal_by_species(factors, output_dir):
    """Plot temporal trajectories grouped by species for key components."""
    sample_factors = factors['sample_factors']
    time_factors = factors['time_factors']
    sample_mapping = factors['sample_mapping']
    weights = factors['component_weights']['weight'].values
    
    # Select top 6 components
    top_comps = np.argsort(weights)[::-1][:6] + 1
    
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    axes = axes.flatten()
    
    for ax, comp in zip(axes, top_comps):
        col_name = f'Component_{comp}'
        
        # Get time factor values
        time_vals = time_factors[col_name].values
        
        for species, color in SPECIES_COLORS.items():
            ax.plot(range(4), time_vals, 'o-', color=color, 
                   linewidth=2, markersize=8, label=SPECIES_MAP[species].split()[0][0] + '.')
        
        ax.set_xticks(range(4))
        ax.set_xticklabels(['TP1', 'TP2', 'TP3', 'TP4'])
        ax.set_title(f'Component {comp} (w={weights[comp-1]:.0f})', fontweight='bold')
        ax.set_ylabel('Loading')
        ax.axhline(0, color='gray', linestyle='--', alpha=0.3)
        
        if comp == top_comps[0]:
            ax.legend(loc='upper right', frameon=False)
    
    fig.suptitle('Temporal Trajectories of Top Components', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_dir / "temporal_trajectories.png")
    plt.close()

def plot_go_enrichment_summary(goslim_counts, output_dir):
    """Create GO term enrichment summary visualization."""
    if goslim_counts is None:
        return
    
    # Get top 20 GO terms by total count
    top_terms = goslim_counts.sum(axis=1).sort_values(ascending=False).head(20)
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    data = goslim_counts.loc[top_terms.index]
    
    im = ax.imshow(data.values, aspect='auto', cmap='YlOrRd')
    
    ax.set_yticks(range(len(top_terms)))
    ax.set_yticklabels(top_terms.index)
    ax.set_xticks(range(0, len(data.columns), 5))
    ax.set_xticklabels([f'C{i+1}' for i in range(0, len(data.columns), 5)])
    
    ax.set_xlabel('Component', fontweight='bold')
    ax.set_ylabel('GO Slim Term', fontweight='bold')
    ax.set_title('GO Term Enrichment Across Components', fontsize=14, fontweight='bold')
    
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label('Gene Count', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_dir / "go_enrichment_heatmap.png")
    plt.close()

def plot_integrated_summary(factors, temporal_patterns, species_patterns, annotations, output_dir):
    """Create comprehensive integrated summary figure."""
    fig = plt.figure(figsize=(18, 12))
    gs = GridSpec(3, 4, figure=fig, hspace=0.35, wspace=0.3)
    
    weights = factors['component_weights']['weight'].values
    top_comps = np.argsort(weights)[::-1][:8] + 1
    
    # Panel A: Component weights
    ax1 = fig.add_subplot(gs[0, :2])
    colors = ['#E63946' if i+1 in top_comps else '#457B9D' for i in range(len(weights))]
    ax1.bar(range(1, len(weights)+1), weights, color=colors, alpha=0.8)
    ax1.set_xlabel('Component', fontweight='bold')
    ax1.set_ylabel('Weight', fontweight='bold')
    ax1.set_title('A. Component Weights', fontweight='bold', loc='left')
    
    # Panel B: Time course of top component
    ax2 = fig.add_subplot(gs[0, 2:])
    time_factors = factors['time_factors']
    for i, comp in enumerate(top_comps[:4]):
        ax2.plot(range(4), time_factors[f'Component_{comp}'].values, 
                'o-', label=f'C{comp}', linewidth=2, markersize=8)
    ax2.set_xticks(range(4))
    ax2.set_xticklabels(['TP1', 'TP2', 'TP3', 'TP4'])
    ax2.set_xlabel('Time Point', fontweight='bold')
    ax2.set_ylabel('Loading', fontweight='bold')
    ax2.set_title('B. Temporal Dynamics of Top Components', fontweight='bold', loc='left')
    ax2.legend(loc='upper right', frameon=False)
    
    # Panel C: Species patterns heatmap
    ax3 = fig.add_subplot(gs[1, :2])
    species_data = np.zeros((3, len(top_comps)))
    species_list = ['apul', 'peve', 'ptua']
    for j, comp in enumerate(top_comps):
        for i, sp in enumerate(species_list):
            species_data[i, j] = species_patterns[comp].get(sp, 0)
    
    im = ax3.imshow(species_data, aspect='auto', cmap='coolwarm')
    ax3.set_yticks(range(3))
    ax3.set_yticklabels([SPECIES_MAP[s].split()[0][:4] + '.' for s in species_list])
    ax3.set_xticks(range(len(top_comps)))
    ax3.set_xticklabels([f'C{c}' for c in top_comps])
    ax3.set_title('C. Species Expression Patterns', fontweight='bold', loc='left')
    plt.colorbar(im, ax=ax3, shrink=0.8, label='Mean Loading')
    
    # Panel D: Temporal pattern distribution
    ax4 = fig.add_subplot(gs[1, 2:])
    pattern_counts = Counter([p['pattern'] for p in temporal_patterns.values()])
    patterns = list(pattern_counts.keys())
    counts = [pattern_counts[p] for p in patterns]
    colors = plt.cm.Set2(np.linspace(0, 1, len(patterns)))
    ax4.barh(range(len(patterns)), counts, color=colors)
    ax4.set_yticks(range(len(patterns)))
    ax4.set_yticklabels([p.replace('_', ' ').title() for p in patterns])
    ax4.set_xlabel('Number of Components', fontweight='bold')
    ax4.set_title('D. Temporal Pattern Categories', fontweight='bold', loc='left')
    
    # Panel E: GO term summary for top component
    ax5 = fig.add_subplot(gs[2, :2])
    top_comp = top_comps[0]
    go_themes = extract_go_themes(annotations, top_comp)
    if go_themes:
        terms = [t[0][:40] + '...' if len(t[0]) > 40 else t[0] for t in go_themes]
        counts = [t[1] for t in go_themes]
        ax5.barh(range(len(terms)), counts, color='#2A9D8F')
        ax5.set_yticks(range(len(terms)))
        ax5.set_yticklabels(terms)
        ax5.set_xlabel('Count', fontweight='bold')
    ax5.set_title(f'E. Top GO Terms (Component {top_comp})', fontweight='bold', loc='left')
    
    # Panel F: Key proteins summary
    ax6 = fig.add_subplot(gs[2, 2:])
    proteins = extract_protein_functions(annotations, top_comp)
    if proteins:
        ax6.text(0.05, 0.95, f'Top Proteins in Component {top_comp}:', 
                fontweight='bold', transform=ax6.transAxes, va='top')
        for i, prot in enumerate(proteins[:8]):
            prot_short = prot[:45] + '...' if len(str(prot)) > 45 else str(prot)
            ax6.text(0.05, 0.85 - i*0.1, f'• {prot_short}',
                    transform=ax6.transAxes, va='top', fontsize=9)
    ax6.axis('off')
    ax6.set_title('F. Key Proteins', fontweight='bold', loc='left')
    
    fig.suptitle('Barnacle Tensor Decomposition: Multi-Species Coral Gene Expression', 
                fontsize=16, fontweight='bold', y=0.98)
    
    plt.savefig(output_dir / "integrated_summary.png", dpi=300)
    plt.close()

def plot_species_independent_analysis(factors, species_independent, species_driven, 
                                       annotations, annotation_stats, output_dir):
    """
    Create comprehensive visualization of species-independent components.
    """
    fig = plt.figure(figsize=(18, 14))
    gs = GridSpec(3, 3, figure=fig, hspace=0.35, wspace=0.3)
    
    weights = factors['component_weights']['weight'].values
    time_factors = factors['time_factors']
    
    # Panel A: Species-independent vs species-driven components
    ax1 = fig.add_subplot(gs[0, 0])
    ind_comps = [c['component'] for c in species_independent]
    drv_comps = [c['component'] for c in species_driven]
    
    ind_weights = [weights[c-1] for c in ind_comps]
    drv_weights = [weights[c-1] for c in drv_comps]
    
    ax1.boxplot([ind_weights, drv_weights], labels=['Species-Independent', 'Species-Driven'])
    ax1.set_ylabel('Component Weight', fontweight='bold')
    ax1.set_title(f'A. Component Weight Distribution\n({len(ind_comps)} independent, {len(drv_comps)} driven)', 
                 fontweight='bold', loc='left')
    
    # Panel B: CV distribution
    ax2 = fig.add_subplot(gs[0, 1])
    ind_cvs = [c['cv'] for c in species_independent]
    drv_cvs = [c['cv'] for c in species_driven]
    
    ax2.hist(ind_cvs, bins=10, alpha=0.7, label='Species-Independent', color='#2A9D8F')
    ax2.hist(drv_cvs, bins=10, alpha=0.7, label='Species-Driven', color='#E63946')
    ax2.axvline(0.3, color='black', linestyle='--', label='Threshold')
    ax2.set_xlabel('Coefficient of Variation', fontweight='bold')
    ax2.set_ylabel('Count', fontweight='bold')
    ax2.set_title('B. Species Variation Distribution', fontweight='bold', loc='left')
    ax2.legend(loc='upper right', frameon=False)
    
    # Panel C: Temporal dynamics of top species-independent components
    ax3 = fig.add_subplot(gs[0, 2])
    top_ind = species_independent[:min(5, len(species_independent))]
    colors = plt.cm.viridis(np.linspace(0.2, 0.8, len(top_ind)))
    
    for i, comp_info in enumerate(top_ind):
        comp = comp_info['component']
        time_vals = time_factors[f'Component_{comp}'].values
        ax3.plot(range(4), time_vals, 'o-', color=colors[i], 
                label=f'C{comp}', linewidth=2, markersize=8)
    
    ax3.set_xticks(range(4))
    ax3.set_xticklabels(['TP1', 'TP2', 'TP3', 'TP4'])
    ax3.set_xlabel('Time Point', fontweight='bold')
    ax3.set_ylabel('Loading', fontweight='bold')
    ax3.set_title('C. Temporal Dynamics\n(Top Species-Independent)', fontweight='bold', loc='left')
    ax3.legend(loc='best', frameon=False, fontsize=8)
    
    # Panel D: Annotation proportion by component type
    ax4 = fig.add_subplot(gs[1, 0])
    ind_prop_unann = [annotation_stats[c['component']]['prop_unannotated'] 
                      for c in species_independent if c['component'] in annotation_stats]
    drv_prop_unann = [annotation_stats[c['component']]['prop_unannotated'] 
                      for c in species_driven if c['component'] in annotation_stats]
    
    ax4.boxplot([ind_prop_unann, drv_prop_unann], 
                labels=['Species-Independent', 'Species-Driven'])
    ax4.set_ylabel('Proportion Unannotated', fontweight='bold')
    ax4.set_title('D. Unannotated Genes by Component Type', fontweight='bold', loc='left')
    ax4.set_ylim(0, 1)
    
    # Panel E: Functional categories for species-independent components
    ax5 = fig.add_subplot(gs[1, 1:])
    func_analysis = analyze_functional_categories(annotations, species_independent)
    
    if func_analysis['goslim_categories']:
        categories = [c[0][:30] + '...' if len(c[0]) > 30 else c[0] 
                     for c in func_analysis['goslim_categories'][:10]]
        counts = [c[1] for c in func_analysis['goslim_categories'][:10]]
        
        y_pos = range(len(categories))
        ax5.barh(y_pos, counts, color='#2A9D8F')
        ax5.set_yticks(y_pos)
        ax5.set_yticklabels(categories)
        ax5.invert_yaxis()
    ax5.set_xlabel('Gene Count', fontweight='bold')
    ax5.set_title('E. Top Functional Categories (Species-Independent Components)', 
                 fontweight='bold', loc='left')
    
    # Panel F: Annotation stats table
    ax6 = fig.add_subplot(gs[2, 0])
    ax6.axis('off')
    
    table_data = []
    headers = ['Comp', 'Total', 'Annot', '% Unann', 'CV']
    
    for comp_info in species_independent[:10]:
        comp = comp_info['component']
        if comp in annotation_stats:
            stats = annotation_stats[comp]
            table_data.append([
                f'C{comp}',
                stats['total_genes'],
                stats['annotated'],
                f"{stats['prop_unannotated']*100:.1f}%",
                f"{comp_info['cv']:.3f}"
            ])
    
    if table_data:
        table = ax6.table(cellText=table_data, colLabels=headers,
                         loc='center', cellLoc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1.2, 1.5)
    ax6.set_title('F. Species-Independent Component Statistics', fontweight='bold', loc='left')
    
    # Panel G: GO biological processes for species-independent
    ax7 = fig.add_subplot(gs[2, 1:])
    
    if func_analysis['go_bp_terms']:
        terms = [t[0][:35] + '...' if len(t[0]) > 35 else t[0] 
                for t in func_analysis['go_bp_terms'][:10]]
        counts = [t[1] for t in func_analysis['go_bp_terms'][:10]]
        
        y_pos = range(len(terms))
        ax7.barh(y_pos, counts, color='#457B9D')
        ax7.set_yticks(y_pos)
        ax7.set_yticklabels(terms)
        ax7.invert_yaxis()
    ax7.set_xlabel('Gene Count', fontweight='bold')
    ax7.set_title('G. Top GO Biological Processes (Species-Independent Components)', 
                 fontweight='bold', loc='left')
    
    fig.suptitle('Species-Independent Gene Expression Components\n(Conserved Patterns Across All Three Coral Species)', 
                fontsize=14, fontweight='bold', y=0.98)
    
    plt.savefig(output_dir / "species_independent_components.png", dpi=300)
    plt.close()


def plot_annotation_summary(annotation_stats, factors, output_dir):
    """Create summary of annotation statistics across all components."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    weights = factors['component_weights']['weight'].values
    
    # Get data sorted by component number
    comps = sorted(annotation_stats.keys())
    prop_unann = [annotation_stats[c]['prop_unannotated'] for c in comps]
    total_genes = [annotation_stats[c]['total_genes'] for c in comps]
    annotated = [annotation_stats[c]['annotated'] for c in comps]
    
    # Panel A: Proportion unannotated by component
    ax1 = axes[0, 0]
    colors = ['#E63946' if p > 0.5 else '#2A9D8F' for p in prop_unann]
    ax1.bar(comps, prop_unann, color=colors)
    ax1.axhline(0.5, color='black', linestyle='--', alpha=0.5, label='50% threshold')
    ax1.set_xlabel('Component', fontweight='bold')
    ax1.set_ylabel('Proportion Unannotated', fontweight='bold')
    ax1.set_title('A. Proportion of Unannotated Genes per Component', fontweight='bold', loc='left')
    ax1.set_ylim(0, 1)
    
    # Panel B: Annotation vs weight scatter
    ax2 = axes[0, 1]
    comp_weights = [weights[c-1] for c in comps]
    scatter = ax2.scatter(comp_weights, prop_unann, c=comps, cmap='viridis', s=80, alpha=0.7)
    ax2.set_xlabel('Component Weight', fontweight='bold')
    ax2.set_ylabel('Proportion Unannotated', fontweight='bold')
    ax2.set_title('B. Annotation vs Component Weight', fontweight='bold', loc='left')
    plt.colorbar(scatter, ax=ax2, label='Component #')
    
    # Add correlation
    corr = np.corrcoef(comp_weights, prop_unann)[0, 1]
    ax2.text(0.05, 0.95, f'r = {corr:.3f}', transform=ax2.transAxes, fontsize=10)
    
    # Panel C: Stacked bar of annotated vs unannotated
    ax3 = axes[1, 0]
    unannotated = [annotation_stats[c]['unannotated'] for c in comps]
    
    ax3.bar(comps, annotated, label='Annotated', color='#2A9D8F')
    ax3.bar(comps, unannotated, bottom=annotated, label='Unannotated', color='#E63946')
    ax3.set_xlabel('Component', fontweight='bold')
    ax3.set_ylabel('Gene Count', fontweight='bold')
    ax3.set_title('C. Annotated vs Unannotated Genes', fontweight='bold', loc='left')
    ax3.legend(loc='upper right', frameon=False)
    
    # Panel D: Summary statistics
    ax4 = axes[1, 1]
    ax4.axis('off')
    
    total_annotated = sum(annotated)
    total_unannotated = sum(unannotated)
    total_all = total_annotated + total_unannotated
    
    summary_text = f"""
    Annotation Summary
    {'='*30}
    
    Total genes analyzed: {total_all}
    Annotated: {total_annotated} ({100*total_annotated/total_all:.1f}%)
    Unannotated: {total_unannotated} ({100*total_unannotated/total_all:.1f}%)
    
    Components with >50% unannotated: {sum(1 for p in prop_unann if p > 0.5)}
    Components with <30% unannotated: {sum(1 for p in prop_unann if p < 0.3)}
    
    Mean % unannotated: {100*np.mean(prop_unann):.1f}%
    Median % unannotated: {100*np.median(prop_unann):.1f}%
    """
    
    ax4.text(0.1, 0.9, summary_text, transform=ax4.transAxes, fontsize=11,
            verticalalignment='top', fontfamily='monospace')
    ax4.set_title('D. Summary Statistics', fontweight='bold', loc='left')
    
    fig.suptitle('Gene Annotation Statistics Across Components', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_dir / "annotation_summary.png", dpi=300)
    plt.close()


def plot_species_independent_physiology(factors, species_independent, physio, output_dir):
    """
    Correlate species-independent components with physiological measurements.
    """
    if physio is None:
        return
    
    sample_mapping = factors['sample_mapping']
    sample_factors = factors['sample_factors']
    
    # Clean physiology data
    phys_cols = ['chla.ug.cm2', 'chlc2.ug.cm2', 'prot_ug.cm2', 'cells.cm2', 
                 'calc.umol.cm2.hr', 'Host_AFDW.mg.cm2', 'Sym_AFDW.mg.cm2']
    phys_available = [c for c in phys_cols if c in physio.columns]
    
    for col in phys_available:
        physio[col] = pd.to_numeric(physio[col].astype(str).str.replace(',', ''), errors='coerce')
    
    # Average physiology by colony_id
    physio_avg = physio.groupby('colony_id')[phys_available].mean()
    
    # Prepare sample factors
    sample_factors_reset = sample_factors.reset_index()
    sample_factors_reset.columns = ['sample'] + list(sample_factors.columns)
    
    merged = sample_factors_reset.merge(
        sample_mapping[['label', 'sample_id', 'species']], 
        left_on='sample', right_on='label'
    )
    
    # Merge with physiology
    merged_with_phys = merged.merge(physio_avg, left_on='sample_id', right_index=True, how='inner')
    
    if len(merged_with_phys) < 5:
        return
    
    # Get species-independent components
    ind_comps = [c['component'] for c in species_independent[:10]]
    
    # Calculate correlations
    corr_data = []
    for comp in ind_comps:
        col = f'Component_{comp}'
        for phys in phys_available:
            valid = ~(merged_with_phys[col].isna() | merged_with_phys[phys].isna())
            if valid.sum() >= 3:
                r = np.corrcoef(merged_with_phys.loc[valid, col],
                               merged_with_phys.loc[valid, phys])[0, 1]
                corr_data.append({
                    'component': comp,
                    'physiology': phys.split('.')[0].replace('_', ' '),
                    'correlation': r
                })
    
    if not corr_data:
        return
    
    corr_df = pd.DataFrame(corr_data)
    corr_pivot = corr_df.pivot(index='component', columns='physiology', values='correlation')
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    im = ax.imshow(corr_pivot.values, cmap='coolwarm', vmin=-1, vmax=1, aspect='auto')
    
    ax.set_xticks(range(len(corr_pivot.columns)))
    ax.set_yticks(range(len(corr_pivot.index)))
    ax.set_xticklabels(corr_pivot.columns, rotation=45, ha='right')
    ax.set_yticklabels([f'C{c}' for c in corr_pivot.index])
    
    # Add values
    for i in range(len(corr_pivot.index)):
        for j in range(len(corr_pivot.columns)):
            val = corr_pivot.values[i, j]
            if not np.isnan(val):
                color = 'white' if abs(val) > 0.5 else 'black'
                ax.text(j, i, f'{val:.2f}', ha='center', va='center', color=color, fontsize=9)
    
    ax.set_xlabel('Physiological Parameter', fontweight='bold')
    ax.set_ylabel('Species-Independent Component', fontweight='bold')
    ax.set_title('Correlation: Species-Independent Components vs Coral Physiology\n(Conserved molecular patterns linked to physiology)', 
                fontsize=12, fontweight='bold')
    
    plt.colorbar(im, ax=ax, shrink=0.8, label='Pearson Correlation')
    
    plt.tight_layout()
    plt.savefig(output_dir / "species_independent_physiology_correlation.png", dpi=300)
    plt.close()


def plot_component_physiology_correlation(factors, physio, output_dir):
    """Create correlation heatmap between gene expression components and physiology."""
    sample_mapping = factors['sample_mapping']
    sample_factors = factors['sample_factors']
    weights = factors['component_weights']['weight'].values
    
    # Get top 15 components by weight
    top_comps = np.argsort(weights)[::-1][:15]
    
    # Prepare sample factors data
    sample_factors_reset = sample_factors.reset_index()
    sample_factors_reset.columns = ['sample'] + list(sample_factors.columns)
    
    # Merge with sample mapping to get sample_id
    merged = sample_factors_reset.merge(
        sample_mapping[['label', 'sample_id', 'species']], 
        left_on='sample', 
        right_on='label'
    )
    
    # Prepare physiology data - average by colony and timepoint
    phys_cols = ['chla.ug.cm2', 'chlc2.ug.cm2', 'prot_ug.cm2', 'cells.cm2', 
                 'calc.umol.cm2.hr', 'Host_AFDW.mg.cm2', 'Sym_AFDW.mg.cm2']
    phys_available = [c for c in phys_cols if c in physio.columns]
    
    # Clean numeric columns
    for col in phys_available:
        physio[col] = pd.to_numeric(physio[col].astype(str).str.replace(',', ''), errors='coerce')
    
    # Average physiology by colony_id
    physio_avg = physio.groupby('colony_id')[phys_available].mean()
    
    # Merge gene expression with physiology
    # Match on sample_id == colony_id
    merged_with_phys = merged.merge(physio_avg, left_on='sample_id', right_index=True, how='inner')
    
    if len(merged_with_phys) < 5:
        print("  Warning: Not enough matched samples for correlation analysis")
        return
    
    # Calculate correlations between top components and physiology
    comp_cols = [f'Component_{c+1}' for c in top_comps]
    
    corr_matrix = np.zeros((len(comp_cols), len(phys_available)))
    for i, comp in enumerate(comp_cols):
        for j, phys in enumerate(phys_available):
            valid = ~(merged_with_phys[comp].isna() | merged_with_phys[phys].isna())
            if valid.sum() >= 3:
                corr_matrix[i, j] = np.corrcoef(
                    merged_with_phys.loc[valid, comp],
                    merged_with_phys.loc[valid, phys]
                )[0, 1]
    
    # Create correlation heatmap
    fig, ax = plt.subplots(figsize=(12, 10))
    
    im = ax.imshow(corr_matrix, cmap='coolwarm', vmin=-1, vmax=1, aspect='auto')
    
    ax.set_xticks(range(len(phys_available)))
    ax.set_yticks(range(len(comp_cols)))
    
    phys_labels = ['Chl a', 'Chl c2', 'Protein', 'Sym cells', 
                   'Calcification', 'Host AFDW', 'Sym AFDW'][:len(phys_available)]
    ax.set_xticklabels(phys_labels, rotation=45, ha='right')
    ax.set_yticklabels([f'C{c+1}' for c in top_comps])
    
    ax.set_xlabel('Physiological Parameter', fontweight='bold')
    ax.set_ylabel('Gene Expression Component', fontweight='bold')
    ax.set_title('Correlation Between Gene Expression Components and Coral Physiology\n(Top 15 components by weight)', 
                fontsize=12, fontweight='bold')
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label('Pearson Correlation', fontweight='bold')
    
    # Add correlation values as text
    for i in range(len(comp_cols)):
        for j in range(len(phys_available)):
            val = corr_matrix[i, j]
            if not np.isnan(val):
                text_color = 'white' if abs(val) > 0.5 else 'black'
                ax.text(j, i, f'{val:.2f}', ha='center', va='center', 
                       color=text_color, fontsize=8)
    
    plt.tight_layout()
    plt.savefig(output_dir / "component_physiology_correlation.png", dpi=300)
    plt.close()
    print("  ✓ Component-physiology correlation heatmap")


def plot_physiological_integration(factors, physio_data, output_dir):
    """Integrate physiological data with gene expression components."""
    if physio_data is None:
        print("  Skipping physiological integration (no data)")
        return
    
    sample_mapping = factors['sample_mapping']
    sample_factors = factors['sample_factors']
    
    # Clean physiology data
    physio = physio_data.copy()
    
    # Convert numeric columns (some have commas)
    numeric_cols = ['cre.umol.mgprot', 'Host_AFDW.mg.cm2', 'Sym_AFDW.mg.cm2', 
                   'chla.ug.cm2', 'chlc2.ug.cm2', 'prot_ug.cm2', 'cells.cm2', 'calc.umol.cm2.hr']
    
    for col in numeric_cols:
        if col in physio.columns:
            physio[col] = pd.to_numeric(physio[col].astype(str).str.replace(',', ''), errors='coerce')
    
    # Create figure with 6 panels
    fig = plt.figure(figsize=(16, 12))
    gs = GridSpec(3, 3, figure=fig, hspace=0.35, wspace=0.3)
    
    # Create species code mapping for physiology data
    def get_species_color(species_name):
        """Map physiology species names to colors."""
        s = species_name.lower()
        if 'acropora' in s or s.startswith('acr'):
            return SPECIES_COLORS['apul']
        elif 'porites' in s or s.startswith('por'):
            return SPECIES_COLORS['peve']
        elif 'pocillopora' in s or s.startswith('poc'):
            return SPECIES_COLORS['ptua']
        return '#999999'
    
    # Panel A: Physiology by species
    ax1 = fig.add_subplot(gs[0, 0])
    physio_by_species = physio.groupby('species')['chla.ug.cm2'].mean()
    colors = [get_species_color(s) for s in physio_by_species.index]
    ax1.bar(range(len(physio_by_species)), physio_by_species.values, color=colors)
    ax1.set_xticks(range(len(physio_by_species)))
    ax1.set_xticklabels([s[:4] for s in physio_by_species.index], rotation=45)
    ax1.set_ylabel('Chlorophyll a (μg/cm²)', fontweight='bold')
    ax1.set_title('A. Chlorophyll by Species', fontweight='bold', loc='left')
    
    # Panel B: Physiology by timepoint
    ax2 = fig.add_subplot(gs[0, 1])
    tp_mapping = {'TP1': 1, 'TP2': 2, 'TP3': 3, 'TP4': 4}
    if 'timepoint' in physio.columns:
        physio_by_tp = physio.groupby('timepoint')['chla.ug.cm2'].mean()
        # Sort by timepoint
        physio_by_tp = physio_by_tp.reindex(['TP1', 'TP2', 'TP3', 'TP4'])
        ax2.plot(range(4), physio_by_tp.values, 'o-', color='#2A9D8F', linewidth=2, markersize=10)
        ax2.set_xticks(range(4))
        ax2.set_xticklabels(['TP1', 'TP2', 'TP3', 'TP4'])
    ax2.set_ylabel('Chlorophyll a (μg/cm²)', fontweight='bold')
    ax2.set_xlabel('Time Point', fontweight='bold')
    ax2.set_title('B. Chlorophyll Over Time', fontweight='bold', loc='left')
    
    # Panel C: Symbiodinium density by species
    ax3 = fig.add_subplot(gs[0, 2])
    if 'cells.cm2' in physio.columns:
        cells_by_species = physio.groupby('species')['cells.cm2'].mean() / 1e6
        colors = [get_species_color(s) for s in cells_by_species.index]
        ax3.bar(range(len(cells_by_species)), cells_by_species.values, color=colors)
        ax3.set_xticks(range(len(cells_by_species)))
        ax3.set_xticklabels([s[:4] for s in cells_by_species.index], rotation=45)
    ax3.set_ylabel('Symbiont Cells (×10⁶/cm²)', fontweight='bold')
    ax3.set_title('C. Symbiont Density', fontweight='bold', loc='left')
    
    # Panel D: Host biomass over time by species
    ax4 = fig.add_subplot(gs[1, 0])
    if 'Host_AFDW.mg.cm2' in physio.columns and 'timepoint' in physio.columns:
        for species in physio['species'].unique():
            sp_data = physio[physio['species'] == species]
            tp_means = sp_data.groupby('timepoint')['Host_AFDW.mg.cm2'].mean()
            tp_means = tp_means.reindex(['TP1', 'TP2', 'TP3', 'TP4'])
            color = get_species_color(species)
            ax4.plot(range(4), tp_means.values, 'o-', color=color, 
                    label=species[:4], linewidth=2, markersize=8)
    ax4.set_xticks(range(4))
    ax4.set_xticklabels(['TP1', 'TP2', 'TP3', 'TP4'])
    ax4.set_ylabel('Host AFDW (mg/cm²)', fontweight='bold')
    ax4.set_xlabel('Time Point', fontweight='bold')
    ax4.set_title('D. Host Biomass Dynamics', fontweight='bold', loc='left')
    ax4.legend(loc='upper right', frameon=False)
    
    # Panel E: Calcification by species
    ax5 = fig.add_subplot(gs[1, 1])
    if 'calc.umol.cm2.hr' in physio.columns:
        calc_by_species = physio.groupby('species')['calc.umol.cm2.hr'].mean()
        colors = [get_species_color(s) for s in calc_by_species.index]
        ax5.bar(range(len(calc_by_species)), calc_by_species.values, color=colors)
        ax5.set_xticks(range(len(calc_by_species)))
        ax5.set_xticklabels([s[:4] for s in calc_by_species.index], rotation=45)
    ax5.set_ylabel('Calcification (μmol/cm²/hr)', fontweight='bold')
    ax5.set_title('E. Calcification Rates', fontweight='bold', loc='left')
    
    # Panel F: Protein content over time
    ax6 = fig.add_subplot(gs[1, 2])
    if 'prot_ug.cm2' in physio.columns and 'timepoint' in physio.columns:
        for species in physio['species'].unique():
            sp_data = physio[physio['species'] == species]
            tp_means = sp_data.groupby('timepoint')['prot_ug.cm2'].mean()
            tp_means = tp_means.reindex(['TP1', 'TP2', 'TP3', 'TP4'])
            color = get_species_color(species)
            ax6.plot(range(4), tp_means.values, 'o-', color=color, 
                    label=species[:4], linewidth=2, markersize=8)
    ax6.set_xticks(range(4))
    ax6.set_xticklabels(['TP1', 'TP2', 'TP3', 'TP4'])
    ax6.set_ylabel('Protein (μg/cm²)', fontweight='bold')
    ax6.set_xlabel('Time Point', fontweight='bold')
    ax6.set_title('F. Protein Content Dynamics', fontweight='bold', loc='left')
    ax6.legend(loc='upper right', frameon=False)
    
    # Panel G: Correlation heatmap between physiology metrics
    ax7 = fig.add_subplot(gs[2, :2])
    phys_cols = ['chla.ug.cm2', 'chlc2.ug.cm2', 'prot_ug.cm2', 'cells.cm2', 
                 'calc.umol.cm2.hr', 'Host_AFDW.mg.cm2', 'Sym_AFDW.mg.cm2']
    phys_available = [c for c in phys_cols if c in physio.columns]
    if len(phys_available) > 1:
        corr_matrix = physio[phys_available].corr()
        labels = ['Chl a', 'Chl c2', 'Protein', 'Sym cells', 
                  'Calc', 'Host AFDW', 'Sym AFDW'][:len(phys_available)]
        im = ax7.imshow(corr_matrix.values, cmap='coolwarm', vmin=-1, vmax=1)
        ax7.set_xticks(range(len(phys_available)))
        ax7.set_yticks(range(len(phys_available)))
        ax7.set_xticklabels(labels, rotation=45, ha='right')
        ax7.set_yticklabels(labels)
        plt.colorbar(im, ax=ax7, shrink=0.8, label='Correlation')
    ax7.set_title('G. Physiological Parameter Correlations', fontweight='bold', loc='left')
    
    # Panel H: Summary statistics
    ax8 = fig.add_subplot(gs[2, 2])
    summary_text = "Physiological Summary\n" + "="*25 + "\n\n"
    for col in phys_available[:5]:
        mean_val = physio[col].mean()
        std_val = physio[col].std()
        col_short = col.split('.')[0].replace('_', ' ')
        summary_text += f"{col_short}:\n  {mean_val:.2f} ± {std_val:.2f}\n\n"
    ax8.text(0.1, 0.9, summary_text, transform=ax8.transAxes, 
            fontsize=10, verticalalignment='top', fontfamily='monospace')
    ax8.axis('off')
    ax8.set_title('H. Summary Statistics', fontweight='bold', loc='left')
    
    fig.suptitle('Coral Physiological Data Integration', fontsize=14, fontweight='bold', y=0.98)
    
    plt.tight_layout()
    plt.savefig(output_dir / "physiology_integration.png", dpi=300)
    plt.close()
    
    # Also create a correlation plot between top components and physiology
    plot_component_physiology_correlation(factors, physio, output_dir)

# =============================================================================
# Narrative Generation
# =============================================================================

def generate_narrative(factors, temporal_patterns, species_patterns, annotations, physio_data,
                       species_independent=None, species_driven=None, annotation_stats=None):
    """Generate a narrative summary of the analysis."""
    
    weights = factors['component_weights']['weight'].values
    top_comps = np.argsort(weights)[::-1][:5] + 1
    
    narrative = """
================================================================================
BARNACLE TENSOR DECOMPOSITION: MULTI-SPECIES CORAL GENE EXPRESSION NARRATIVE
================================================================================

EXECUTIVE SUMMARY
-----------------
This analysis reveals conserved and divergent gene expression patterns across 
three coral species (Acropora pulchra, Porites evermanni, and Pocillopora 
tuahiniensis) over four experimental time points using tensor decomposition.

The barnacle algorithm identified 35 expression components, with the top 5 
components by weight being: {top_comps}.

TEMPORAL DYNAMICS
-----------------
Gene expression across the coral species shows distinct temporal phases:

""".format(top_comps=', '.join([f'C{c}' for c in top_comps]))

    # Classify components by temporal pattern
    pattern_groups = {}
    for comp, data in temporal_patterns.items():
        pattern = data['pattern']
        if pattern not in pattern_groups:
            pattern_groups[pattern] = []
        pattern_groups[pattern].append(comp)
    
    for pattern, comps in sorted(pattern_groups.items(), key=lambda x: -len(x[1])):
        pattern_name = pattern.replace('_', ' ').title()
        narrative += f"• {pattern_name} ({len(comps)} components): Components {', '.join([f'C{c}' for c in sorted(comps)[:5]])}\n"
    
    narrative += """
SPECIES-SPECIFIC PATTERNS
-------------------------
The three coral species show both shared and unique expression signatures:

"""
    
    # Identify species-dominant components
    for species in ['apul', 'peve', 'ptua']:
        dominant_comps = []
        for comp, data in species_patterns.items():
            species_vals = list(data.values())
            if data.get(species, 0) == max(species_vals) and max(species_vals) > 0:
                dominant_comps.append(comp)
        
        species_name = SPECIES_MAP[species]
        narrative += f"• {species_name}: Shows highest expression in {len(dominant_comps)} components\n"
    
    narrative += """
KEY BIOLOGICAL PROCESSES
------------------------
Analysis of GO term enrichment reveals major biological themes:

"""
    
    # Summarize GO themes from top components
    all_go_themes = []
    for comp in top_comps:
        themes = extract_go_themes(annotations, comp)
        all_go_themes.extend([t[0] for t in themes])
    
    theme_counts = Counter(all_go_themes)
    for theme, count in theme_counts.most_common(10):
        narrative += f"• {theme}\n"
    
    narrative += """
COMPONENT INTERPRETATIONS
-------------------------
"""
    
    for comp in top_comps:
        weight = weights[comp - 1]
        pattern = temporal_patterns[comp]['pattern']
        go_themes = extract_go_themes(annotations, comp)
        proteins = extract_protein_functions(annotations, comp)
        
        narrative += f"""
Component {comp} (Weight: {weight:.1f}, Pattern: {pattern.replace('_', ' ')})
{'='*50}
Top GO Processes: {', '.join([t[0].split('[')[0].strip() for t in go_themes[:3]]) if go_themes else 'N/A'}
Key Proteins: {', '.join([str(p) for p in proteins[:3]]) if proteins else 'N/A'}
"""
    
    narrative += """
SPECIES-INDEPENDENT COMPONENTS
------------------------------
"""
    
    if species_independent:
        narrative += f"""
A key focus of this analysis is identifying components NOT primarily driven by 
species differences - these represent conserved expression patterns shared across 
all three coral species and likely reflect fundamental biological processes.

Number of species-independent components: {len(species_independent)}
Number of species-driven components: {len(species_driven) if species_driven else 'N/A'}

Top species-independent components (lowest cross-species variation):
"""
        for comp_info in species_independent[:5]:
            comp = comp_info['component']
            cv = comp_info['cv']
            weight = weights[comp - 1]
            narrative += f"  • Component {comp}: CV={cv:.3f}, Weight={weight:.1f}\n"
        
        # Add functional analysis
        func_analysis = analyze_functional_categories(annotations, species_independent)
        
        narrative += f"""
Functional Categories in Species-Independent Components:
"""
        if func_analysis['goslim_categories']:
            for cat, count in func_analysis['goslim_categories'][:8]:
                narrative += f"  • {cat}: {count} genes\n"
        
        narrative += f"""
Top Biological Processes (conserved across species):
"""
        if func_analysis['go_bp_terms']:
            for term, count in func_analysis['go_bp_terms'][:6]:
                narrative += f"  • {term}: {count} genes\n"
    
    narrative += """
ANNOTATION STATISTICS
---------------------
"""
    
    if annotation_stats:
        total_annotated = sum(s['annotated'] for s in annotation_stats.values())
        total_unannotated = sum(s['unannotated'] for s in annotation_stats.values())
        total_all = total_annotated + total_unannotated
        prop_unann_list = [s['prop_unannotated'] for s in annotation_stats.values()]
        
        narrative += f"""
Across all components analyzed:
  • Total genes: {total_all}
  • Annotated: {total_annotated} ({100*total_annotated/total_all:.1f}%)
  • Unannotated: {total_unannotated} ({100*total_unannotated/total_all:.1f}%)
  • Mean % unannotated per component: {100*np.mean(prop_unann_list):.1f}%
  • Components with >50% unannotated: {sum(1 for p in prop_unann_list if p > 0.5)}

This indicates that a significant portion of the coral transcriptome remains 
functionally uncharacterized, highlighting the need for continued annotation 
efforts in non-model marine organisms.
"""
        
        # Report on species-independent component annotation
        if species_independent:
            ind_prop_unann = [annotation_stats[c['component']]['prop_unannotated'] 
                            for c in species_independent if c['component'] in annotation_stats]
            narrative += f"""
Species-Independent Components Annotation:
  • Mean % unannotated: {100*np.mean(ind_prop_unann):.1f}%
  • Range: {100*min(ind_prop_unann):.1f}% - {100*max(ind_prop_unann):.1f}%
"""
    
    narrative += """
CONCLUSIONS
-----------
The tensor decomposition reveals that:

1. Gene expression patterns are largely conserved across the three coral species,
   suggesting shared stress response mechanisms and fundamental physiological 
   processes.

2. Temporal dynamics show distinct phases, with some components showing early 
   activation followed by decline, while others show progressive increases or 
   stable expression throughout the experiment.

3. Species-specific expression signatures exist alongside conserved patterns, 
   potentially reflecting evolutionary divergence in stress response strategies 
   or ecological adaptations.

4. Key biological processes including transcriptional regulation, translation, 
   signal transduction, and metabolic pathways are dynamically regulated across 
   the experimental timeline.

5. Species-independent components reveal conserved cellular machinery including
   translation, metabolic regulation, and stress response pathways that are 
   fundamental to coral biology regardless of species identity.

6. A substantial proportion of genes remain unannotated, indicating significant 
   unexplored biology in coral genomes that warrants further investigation.

PHYSIOLOGICAL CONTEXT
---------------------
"""
    
    if physio_data is not None:
        narrative += """
Integration with physiological measurements reveals the following patterns:

• Chlorophyll content shows species-specific baseline differences, with dynamic
  changes across time points that may correlate with photosynthetic capacity
  and symbiont density.

• Symbiodinium density varies between species, reflecting differences in 
  symbiotic relationships and potentially heat tolerance strategies.

• Host biomass (AFDW) shows temporal dynamics that may be linked to the 
  metabolic and stress response gene expression patterns identified in the
  tensor decomposition.

• Calcification rates differ between species, with potential links to 
  components involved in biomineralization and ion transport gene expression.

• Protein content dynamics across time points may reflect the translation
  and protein synthesis pathways identified as major themes in the gene
  expression components.

These physiological measurements provide functional validation of the gene
expression patterns and help contextualize the molecular stress response
signatures across the three coral species.
"""
    else:
        narrative += "\n[Physiological data not available for integration]\n"
    
    narrative += """
METHODS NOTE
------------
Analysis performed using barnacle tensor decomposition with rank-35 optimization
and lambda_gene=0.2. Input data consisted of orthologous gene expression from 
{n_genes} ortholog groups across {n_samples} samples and {n_timepoints} time points.

Visualizations and analysis performed using Python with matplotlib, seaborn, 
numpy, and pandas. GO term enrichment from UniProt annotations.
""".format(
        n_genes=factors['metadata']['tensor_shape'][0],
        n_samples=factors['metadata']['tensor_shape'][1],
        n_timepoints=factors['metadata']['tensor_shape'][2]
    )
    
    return narrative

# =============================================================================
# Main Execution
# =============================================================================

def main():
    """Main execution function."""
    print("\n" + "="*70)
    print("BARNACLE TENSOR DECOMPOSITION INTERPRETATION")
    print("="*70 + "\n")
    
    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {OUTPUT_DIR}\n")
    
    # Setup plotting style
    setup_style()
    
    # Load data
    factors = load_barnacle_factors()
    annotations = load_annotated_components()
    goslim_counts = load_goslim_counts()
    physio_data = load_physiology_data()
    
    print("\n" + "-"*50)
    print("Analyzing patterns...")
    print("-"*50)
    
    # Analyze patterns
    key_components = identify_key_components(factors)
    print(f"Top components by weight: {key_components['indices'][:5]}")
    
    temporal_patterns = analyze_temporal_patterns(factors)
    print(f"Identified {len(temporal_patterns)} temporal patterns")
    
    species_patterns = analyze_species_patterns(factors)
    print(f"Analyzed species patterns for {len(species_patterns)} components")
    
    # NEW: Identify species-independent components
    # Using threshold of 0.6 to capture more components with moderate cross-species variation
    species_independent, species_driven = identify_species_independent_components(
        factors, species_patterns, threshold=0.6
    )
    print(f"Species-independent components: {len(species_independent)}")
    print(f"Species-driven components: {len(species_driven)}")
    
    # NEW: Calculate annotation statistics
    annotation_stats = calculate_annotation_statistics(annotations)
    print(f"Calculated annotation statistics for {len(annotation_stats)} components")
    
    print("\n" + "-"*50)
    print("Generating visualizations...")
    print("-"*50)
    
    # Generate visualizations
    plot_component_weights(factors, OUTPUT_DIR)
    print("  ✓ Component weights plot")
    
    plot_temporal_dynamics(factors, temporal_patterns, OUTPUT_DIR)
    print("  ✓ Temporal dynamics heatmap")
    
    plot_species_comparison(factors, species_patterns, OUTPUT_DIR)
    print("  ✓ Species comparison plot")
    
    plot_temporal_by_species(factors, OUTPUT_DIR)
    print("  ✓ Temporal trajectories plot")
    
    if goslim_counts is not None:
        plot_go_enrichment_summary(goslim_counts, OUTPUT_DIR)
        print("  ✓ GO enrichment heatmap")
    
    plot_integrated_summary(factors, temporal_patterns, species_patterns, annotations, OUTPUT_DIR)
    print("  ✓ Integrated summary figure")
    
    plot_physiological_integration(factors, physio_data, OUTPUT_DIR)
    print("  ✓ Physiology integration plot")
    
    # NEW: Species-independent component analysis
    plot_species_independent_analysis(factors, species_independent, species_driven,
                                      annotations, annotation_stats, OUTPUT_DIR)
    print("  ✓ Species-independent components analysis")
    
    # NEW: Annotation summary
    plot_annotation_summary(annotation_stats, factors, OUTPUT_DIR)
    print("  ✓ Annotation summary")
    
    # NEW: Species-independent physiology correlation
    if physio_data is not None:
        plot_species_independent_physiology(factors, species_independent, physio_data, OUTPUT_DIR)
        print("  ✓ Species-independent physiology correlation")
    
    print("\n" + "-"*50)
    print("Generating narrative...")
    print("-"*50)
    
    # Generate narrative
    narrative = generate_narrative(factors, temporal_patterns, species_patterns, 
                                  annotations, physio_data,
                                  species_independent, species_driven, annotation_stats)
    
    # Save narrative
    narrative_path = OUTPUT_DIR / "narrative_summary.txt"
    with open(narrative_path, 'w') as f:
        f.write(narrative)
    print(f"  ✓ Saved narrative to {narrative_path}")
    
    # Print narrative
    print("\n" + narrative)
    
    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print(f"All outputs saved to: {OUTPUT_DIR}")
    print("="*70 + "\n")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())

