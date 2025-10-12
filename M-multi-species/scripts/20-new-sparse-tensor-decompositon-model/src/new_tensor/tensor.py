"""Tensor construction and manipulation utilities."""

import logging
from typing import Dict, List, Tuple, Optional
import json
from pathlib import Path

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def build_tensor(
    df: pd.DataFrame,
    species_list: List[str],
    timepoints_list: List[int]
) -> Tuple[np.ndarray, Dict, Dict, Dict]:
    """Build 3D tensor from aggregated gene expression data.

    Args:
        df: Aggregated DataFrame (genes × species_timepoint)
        species_list: List of species codes
        timepoints_list: List of timepoint numbers

    Returns:
        Tuple of (tensor, gene_mapping, species_mapping, timepoint_mapping)
    """
    logger.info("Building 3D tensor from aggregated data")

    # Extract species-timepoint column names
    st_columns = [col for col in df.columns if col != 'group_id']
    n_genes = len(df)
    n_species = len(species_list)
    n_timepoints = len(timepoints_list)

    logger.info(f"Tensor dimensions: {n_genes} genes × {n_species} species × {n_timepoints} timepoints")

    # Initialize tensor
    tensor = np.zeros((n_genes, n_species, n_timepoints))

    # Create mappings
    gene_mapping = {gene_id: idx for idx, gene_id in enumerate(df['group_id'])}
    species_mapping = {species: idx for idx, species in enumerate(species_list)}
    timepoint_mapping = {tp: idx for idx, tp in enumerate(timepoints_list)}

    # Fill tensor
    for col in st_columns:
        # Parse species and timepoint from column name
        # Format: SPECIES-TP{N} (e.g., ACR-TP1, POR-TP2)
        parts = col.split('-')
        if len(parts) != 2:
            logger.warning(f"Skipping column {col} - unexpected format")
            continue

        species = parts[0]
        tp_str = parts[1]

        # Parse timepoint number from TP{N} format
        if tp_str.startswith('TP'):
            tp = int(tp_str.replace('TP', ''))
        else:
            logger.warning(f"Skipping column {col} - unexpected timepoint format: {tp_str}")
            continue

        if species not in species_mapping or tp not in timepoint_mapping:
            logger.warning(f"Skipping column {col} - species or timepoint not in mappings")
            continue

        species_idx = species_mapping[species]
        timepoint_idx = timepoint_mapping[tp]

        # Get data for this species-timepoint combination
        tensor[:, species_idx, timepoint_idx] = df[col].values

    # Verify tensor construction
    non_zero_entries = np.count_nonzero(~np.isnan(tensor))
    total_entries = tensor.size
    logger.info(f"Tensor density: {non_zero_entries}/{total_entries} non-zero entries")

    return tensor, gene_mapping, species_mapping, timepoint_mapping


def save_tensor_mappings(
    output_dir: str,
    gene_mapping: Dict[str, int],
    species_mapping: Dict[str, int],
    timepoint_mapping: Dict[str, int],
    species_list: List[str],
    timepoints_list: List[int]
) -> None:
    """Save tensor mappings to disk.

    Args:
        output_dir: Output directory path
        gene_mapping: Gene ID to index mapping
        species_mapping: Species code to index mapping
        timepoint_mapping: Timepoint string to index mapping
        species_list: List of species in order
        timepoints_list: List of timepoints in order
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Save gene mapping
    gene_df = pd.DataFrame([
        {'gene_id': gene_id, 'index': idx}
        for gene_id, idx in gene_mapping.items()
    ])
    gene_df.to_csv(output_path / 'genes.csv', index=False)

    # Save species mapping
    species_df = pd.DataFrame([
        {'species': species, 'index': idx}
        for species, idx in species_mapping.items()
    ])
    species_df.to_csv(output_path / 'species.csv', index=False)

    # Save timepoint mapping
    tp_df = pd.DataFrame([
        {'timepoint': tp, 'index': idx}
        for tp, idx in timepoint_mapping.items()
    ])
    tp_df.to_csv(output_path / 'timepoints.csv', index=False)

    # Save tensor shape info
    shape_info = {
        'n_genes': len(gene_mapping),
        'n_species': len(species_list),
        'n_timepoints': len(timepoints_list),
        'species_order': species_list,
        'timepoints_order': [f"TP{tp}" for tp in timepoints_list]
    }

    with open(output_path / 'tensor_shapes.json', 'w') as f:
        json.dump(shape_info, f, indent=2)

    logger.info(f"Saved tensor mappings to {output_path}")


def load_tensor_mappings(output_dir: str) -> Tuple[Dict, Dict, Dict, Dict]:
    """Load tensor mappings from disk.

    Args:
        output_dir: Directory containing saved mappings

    Returns:
        Tuple of (gene_mapping, species_mapping, timepoint_mapping, shape_info)
    """
    output_path = Path(output_dir)

    # Load gene mapping
    gene_df = pd.read_csv(output_path / 'genes.csv')
    gene_mapping = dict(zip(gene_df['gene_id'], gene_df['index']))

    # Load species mapping
    species_df = pd.read_csv(output_path / 'species.csv')
    species_mapping = dict(zip(species_df['species'], species_df['index']))

    # Load timepoint mapping
    tp_df = pd.read_csv(output_path / 'timepoints.csv')
    timepoint_mapping = dict(zip(tp_df['timepoint'], tp_df['index']))

    # Load shape info
    with open(output_path / 'tensor_shapes.json', 'r') as f:
        shape_info = json.load(f)

    return gene_mapping, species_mapping, timepoint_mapping, shape_info


def save_tensor(tensor: np.ndarray, output_dir: str, filename: str = 'tensor.npz') -> None:
    """Save tensor to disk.

    Args:
        tensor: 3D tensor to save
        output_dir: Output directory
        filename: Output filename
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    filepath = output_path / filename
    np.savez_compressed(filepath, tensor=tensor)

    logger.info(f"Saved tensor to {filepath}")


def load_tensor(output_dir: str, filename: str = 'tensor.npz') -> np.ndarray:
    """Load tensor from disk.

    Args:
        output_dir: Directory containing saved tensor
        filename: Tensor filename

    Returns:
        Loaded 3D tensor
    """
    filepath = Path(output_dir) / filename
    data = np.load(filepath)
    tensor = data['tensor']

    logger.info(f"Loaded tensor from {filepath} with shape {tensor.shape}")
    return tensor
