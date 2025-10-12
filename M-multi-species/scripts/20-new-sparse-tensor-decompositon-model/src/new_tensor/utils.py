"""Utility functions for the tensor decomposition package."""

import logging
import random
from typing import Dict, Any
from pathlib import Path

import numpy as np

logger = logging.getLogger(__name__)


def setup_logging(level: str = 'INFO', log_file: str = None) -> None:
    """Set up logging configuration.

    Args:
        level: Logging level ('DEBUG', 'INFO', 'WARNING', 'ERROR')
        log_file: Optional log file path
    """
    numeric_level = getattr(logging, level.upper(), logging.INFO)

    # Create formatters
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    # Set up console handler
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    console_handler.setLevel(numeric_level)

    # Set up file handler if specified
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        file_handler.setLevel(numeric_level)

    # Configure root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(numeric_level)
    root_logger.addHandler(console_handler)

    if log_file:
        root_logger.addHandler(file_handler)

    logger.info(f"Logging configured with level {level}")


def set_random_seed(seed: int = 42) -> None:
    """Set random seed for reproducible results.

    Args:
        seed: Random seed value
    """
    np.random.seed(seed)
    random.seed(seed)

    # Set TensorFlow seed if available
    try:
        import tensorflow as tf
        tf.random.set_seed(seed)
    except ImportError:
        pass

    logger.info(f"Random seed set to {seed}")


def load_config(config_path: str) -> Dict[str, Any]:
    """Load configuration from YAML file.

    Args:
        config_path: Path to YAML config file

    Returns:
        Configuration dictionary
    """
    import yaml

    config_path = Path(config_path)
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    logger.info(f"Loaded configuration from {config_path}")
    return config


def save_config(config: Dict[str, Any], config_path: str) -> None:
    """Save configuration to YAML file.

    Args:
        config: Configuration dictionary
        config_path: Path to save config file
    """
    import yaml

    config_path = Path(config_path)
    config_path.parent.mkdir(parents=True, exist_ok=True)

    with open(config_path, 'w') as f:
        yaml.dump(config, f, default_flow_style=False, indent=2)

    logger.info(f"Saved configuration to {config_path}")


def get_tensor_info(tensor_path: str) -> Dict[str, Any]:
    """Get information about a saved tensor.

    Args:
        tensor_path: Path to tensor .npz file

    Returns:
        Dictionary with tensor information
    """
    tensor_path = Path(tensor_path)

    if not tensor_path.exists():
        raise FileNotFoundError(f"Tensor file not found: {tensor_path}")

    # Load tensor
    data = np.load(tensor_path)
    tensor = data['tensor']

    info = {
        'shape': tensor.shape,
        'dtype': str(tensor.dtype),
        'size': tensor.size,
        'nnz': np.count_nonzero(~np.isnan(tensor)),
        'density': np.count_nonzero(~np.isnan(tensor)) / tensor.size,
        'file_path': str(tensor_path),
        'file_size': tensor_path.stat().st_size
    }

    return info


def validate_output_directory(output_dir: str) -> Path:
    """Validate and create output directory if needed.

    Args:
        output_dir: Output directory path

    Returns:
        Path object for output directory
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    logger.info(f"Output directory: {output_path.absolute()}")
    return output_path
