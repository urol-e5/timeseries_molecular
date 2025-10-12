"""Preprocessing utilities for tensor construction."""

import logging
from typing import Dict, List, Optional, Union, Literal
import numpy as np
import pandas as pd
from scipy import stats

logger = logging.getLogger(__name__)


def aggregate_replicates(
    df: pd.DataFrame,
    sample_info: Dict[str, Dict[str, str]],
    method: Literal['mean', 'median', 'trimmed_mean'] = 'median',
    trim_fraction: float = 0.1
) -> pd.DataFrame:
    """Aggregate replicate samples across individuals for each species-timepoint combination.

    Args:
        df: Input DataFrame with gene expression data
        sample_info: Parsed sample metadata from parse_sample_labels
        method: Aggregation method ('mean', 'median', 'trimmed_mean')
        trim_fraction: Fraction to trim for trimmed mean (default 0.1)

    Returns:
        DataFrame with aggregated values (genes × species_timepoint)
    """
    logger.info(f"Aggregating replicates using {method} method")

    # Group samples by species and timepoint
    species_timepoint_groups = {}
    for col, info in sample_info.items():
        key = (info['species'], f"TP{info['timepoint']}")
        if key not in species_timepoint_groups:
            species_timepoint_groups[key] = []
        species_timepoint_groups[key].append(col)

    # Aggregate replicates for each group
    aggregated_data = []
    new_columns = []

    for (species, timepoint), columns in species_timepoint_groups.items():
        if len(columns) == 0:
            continue

        # Extract data for this species-timepoint combination
        group_data = df[columns].values

        if method == 'mean':
            aggregated_values = np.nanmean(group_data, axis=1)
        elif method == 'median':
            aggregated_values = np.nanmedian(group_data, axis=1)
        elif method == 'trimmed_mean':
            aggregated_values = stats.trim_mean(group_data, proportiontocut=trim_fraction, axis=1)
        else:
            raise ValueError(f"Unknown aggregation method: {method}")

        aggregated_data.append(aggregated_values)
        new_columns.append(f"{species}-{timepoint}")

    # Create result DataFrame
    result_df = pd.DataFrame(
        np.column_stack(aggregated_data),
        index=df.index,
        columns=new_columns
    )

    # Add group_id column back
    result_df.insert(0, 'group_id', df['group_id'])

    logger.info(f"Aggregated to {len(new_columns)} species-timepoint combinations")
    return result_df


def normalize_tensor(
    df: pd.DataFrame,
    method: Literal['none', 'zscore', 'minmax', 'log1p'] = 'none',
    axis: int = 0
) -> pd.DataFrame:
    """Normalize gene expression data.

    Args:
        df: DataFrame to normalize (excluding group_id column)
        method: Normalization method
        axis: Axis along which to normalize (0=genes, 1=samples)

    Returns:
        Normalized DataFrame
    """
    if method == 'none':
        logger.info("Skipping normalization")
        return df

    logger.info(f"Normalizing data using {method} method")

    # Work on numeric columns only
    numeric_cols = [col for col in df.columns if col != 'group_id']
    data = df[numeric_cols].values

    if method == 'zscore':
        # Z-score normalization per gene (axis=0) or per sample (axis=1)
        if axis == 0:
            # Normalize each gene across samples
            mean = np.nanmean(data, axis=1, keepdims=True)
            std = np.nanstd(data, axis=1, keepdims=True)
            normalized = (data - mean) / (std + 1e-8)
        else:
            # Normalize each sample across genes
            mean = np.nanmean(data, axis=0, keepdims=True)
            std = np.nanstd(data, axis=0, keepdims=True)
            normalized = (data - mean) / (std + 1e-8)

    elif method == 'minmax':
        # Min-max normalization
        if axis == 0:
            min_vals = np.nanmin(data, axis=1, keepdims=True)
            max_vals = np.nanmax(data, axis=1, keepdims=True)
            normalized = (data - min_vals) / (max_vals - min_vals + 1e-8)
        else:
            min_vals = np.nanmin(data, axis=0, keepdims=True)
            max_vals = np.nanmax(data, axis=0, keepdims=True)
            normalized = (data - min_vals) / (max_vals - min_vals + 1e-8)

    elif method == 'log1p':
        # Log transformation (log1p to handle zeros)
        normalized = np.log1p(data)

    else:
        raise ValueError(f"Unknown normalization method: {method}")

    # Create result DataFrame
    result_df = df.copy()
    result_df[numeric_cols] = normalized

    return result_df


def filter_genes(
    df: pd.DataFrame,
    min_expression: float = 1.0,
    min_variance_percentile: float = 10.0,
    max_genes: Optional[int] = None
) -> pd.DataFrame:
    """Filter genes based on expression levels and variance.

    Args:
        df: Input DataFrame
        min_expression: Minimum expression threshold
        min_variance_percentile: Percentile threshold for variance filtering
        max_genes: Maximum number of genes to keep (optional)

    Returns:
        Filtered DataFrame
    """
    logger.info("Filtering genes based on expression and variance")

    # Calculate gene-wise statistics
    numeric_cols = [col for col in df.columns if col != 'group_id']
    data = df[numeric_cols].values

    # Filter by minimum expression
    gene_means = np.mean(data, axis=1)
    gene_vars = np.var(data, axis=1)

    # Expression filter
    expr_mask = gene_means >= min_expression

    # Variance filter (keep genes above percentile)
    var_threshold = np.percentile(gene_vars, min_variance_percentile)
    var_mask = gene_vars >= var_threshold

    # Combine masks
    keep_mask = expr_mask & var_mask

    logger.info(f"Expression filter: {expr_mask.sum()}/{len(expr_mask)} genes")
    logger.info(f"Variance filter: {var_mask.sum()}/{len(var_mask)} genes")
    logger.info(f"Combined filter: {keep_mask.sum()}/{len(keep_mask)} genes")

    # Apply filters
    filtered_df = df[keep_mask].copy()

    # Optionally limit number of genes
    if max_genes is not None and len(filtered_df) > max_genes:
        logger.info(f"Limiting to top {max_genes} genes by variance")
        top_genes_idx = np.argsort(gene_vars[keep_mask])[-max_genes:]
        filtered_df = filtered_df.iloc[top_genes_idx]

    return filtered_df
