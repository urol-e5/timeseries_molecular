"""New Tensor Decomposition Package for Multi-Species Time Series Analysis."""

__version__ = "0.1.0"

from .io import load_csv, parse_sample_labels
from .preprocess import aggregate_replicates, normalize_tensor
from .tensor import build_tensor, save_tensor_mappings
from .model import sparse_cp_decomposition
from .optimize import optimize_hyperparameters
from .cli import app

__all__ = [
    "load_csv",
    "parse_sample_labels",
    "aggregate_replicates",
    "normalize_tensor",
    "build_tensor",
    "save_tensor_mappings",
    "sparse_cp_decomposition",
    "optimize_hyperparameters",
    "app",
]
