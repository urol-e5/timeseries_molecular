"""Input/Output utilities for loading and parsing molecular data."""

import re
import logging
from typing import Dict, List, Tuple, Optional
from pathlib import Path

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def load_csv(file_path: str) -> pd.DataFrame:
    """Load CSV file containing gene expression data.

    Args:
        file_path: Path to the CSV file

    Returns:
        DataFrame with gene expression data

    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If required columns are missing
    """
    file_path = Path(file_path)

    if not file_path.exists():
        raise FileNotFoundError(f"Input file not found: {file_path}")

    logger.info(f"Loading CSV from {file_path}")
    df = pd.read_csv(file_path)

    # Validate required columns
    if 'group_id' not in df.columns:
        raise ValueError("CSV must contain 'group_id' column")

    if df.shape[1] < 2:
        raise ValueError("CSV must contain at least one sample column")

    logger.info(f"Loaded {df.shape[0]} genes with {df.shape[1] - 1} samples")
    return df


def parse_sample_labels(sample_columns: List[str]) -> Dict[str, Dict[str, str]]:
    """Parse sample column names to extract species, individual, and timepoint.

    Expected format: SPECIES-INDIVIDUAL-TP{N}
    Example: ACR-139-TP1, POR-216-TP2

    Args:
        sample_columns: List of sample column names

    Returns:
        Dictionary mapping column names to parsed components

    Raises:
        ValueError: If sample names don't match expected pattern
    """
    pattern = r'^(?P<species>[A-Z]{3})-(?P<individual>[^-]+)-TP(?P<timepoint>\d+)$'

    parsed = {}
    species_codes = set()
    timepoints = set()

    for col in sample_columns:
        match = re.match(pattern, col)
        if not match:
            raise ValueError(
                f"Sample name '{col}' doesn't match expected pattern. "
                "Expected format: SPECIES-INDIVIDUAL-TP{N}"
            )

        components = match.groupdict()
        parsed[col] = components

        species_codes.add(components['species'])
        timepoints.add(int(components['timepoint']))

    logger.info(f"Found species codes: {sorted(species_codes)}")
    logger.info(f"Found timepoints: {sorted(timepoints)}")

    return parsed


def extract_sample_metadata(df: pd.DataFrame) -> Tuple[Dict, List, List]:
    """Extract sample metadata from DataFrame.

    Args:
        df: Input DataFrame with gene expression data

    Returns:
        Tuple of (sample_info, species_list, timepoints_list)
    """
    sample_columns = [col for col in df.columns if col != 'group_id']

    # Parse sample labels
    sample_info = parse_sample_labels(sample_columns)

    # Extract unique species and timepoints
    species_list = sorted(set(info['species'] for info in sample_info.values()))
    timepoints_list = sorted(set(int(info['timepoint']) for info in sample_info.values()))

    return sample_info, species_list, timepoints_list


def validate_input_data(df: pd.DataFrame) -> None:
    """Validate input DataFrame for tensor construction.

    Args:
        df: Input DataFrame to validate

    Raises:
        ValueError: If validation fails
    """
    # Check for missing values
    missing_count = df.isnull().sum().sum()
    if missing_count > 0:
        logger.warning(f"Found {missing_count} missing values in input data")

    # Check for infinite values
    inf_count = np.isinf(df.select_dtypes(include=[np.number])).sum().sum()
    if inf_count > 0:
        raise ValueError(f"Found {inf_count} infinite values in input data")

    # Check for negative values (expression data should be non-negative)
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    if not numeric_cols.empty:
        negative_count = (df[numeric_cols] < 0).sum().sum()
        if negative_count > 0:
            logger.warning(f"Found {negative_count} negative values in expression data")

    logger.info("Input data validation completed")
