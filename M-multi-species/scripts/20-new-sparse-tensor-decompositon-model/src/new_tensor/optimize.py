"""Hyperparameter optimization for sparse CP decomposition."""

import logging
import json
import pickle
from typing import Dict, List, Optional, Tuple, Callable
from pathlib import Path

import numpy as np
import optuna
from optuna.samplers import TPESampler
import pandas as pd

from .model import sparse_cp_decomposition

logger = logging.getLogger(__name__)


def masked_rmse_loss(
    tensor: np.ndarray,
    factors: List[np.ndarray],
    mask_fraction: float = 0.2
) -> float:
    """Compute masked RMSE for cross-validation.

    Args:
        tensor: Original tensor
        factors: CP decomposition factors
        mask_fraction: Fraction of entries to mask for CV

    Returns:
        Masked RMSE loss
    """
    # Create reconstruction
    reconstruction = np.einsum('ir,jr,kr->ijk', factors[0], factors[1], factors[2])

    # Create mask for held-out entries
    np.random.seed(42)
    mask = np.random.random(tensor.shape) > mask_fraction

    # Compute masked error
    masked_error = np.linalg.norm((tensor - reconstruction) * mask) ** 2
    n_masked = np.sum(mask)

    if n_masked == 0:
        return float('inf')

    rmse = np.sqrt(masked_error / n_masked)

    # Add sparsity penalty
    sparsity_penalty = 0.01 * (
        np.sum(np.abs(factors[0])) +
        np.sum(np.abs(factors[1])) +
        np.sum(np.abs(factors[2]))
    )

    return rmse + sparsity_penalty


def objective_function(
    trial: optuna.Trial,
    tensor: np.ndarray,
    rank_range: Tuple[int, int],
    lambda_bounds: Tuple[float, float],
    eval_metric: str = 'masked_rmse'
) -> float:
    """Optuna objective function for hyperparameter optimization.

    Args:
        trial: Optuna trial object
        tensor: Input tensor
        rank_range: (min_rank, max_rank) tuple
        lambda_bounds: (min_lambda, max_lambda) tuple for L1 penalties
        eval_metric: Evaluation metric ('masked_rmse' or 'reconstruction_error')

    Returns:
        Objective value to minimize
    """
    # Sample hyperparameters
    rank = trial.suggest_int('rank', rank_range[0], rank_range[1])
    lambda_A = trial.suggest_float('lambda_A', lambda_bounds[0], lambda_bounds[1], log=True)
    lambda_B = trial.suggest_float('lambda_B', lambda_bounds[0], lambda_bounds[1], log=True)
    lambda_C = trial.suggest_float('lambda_C', lambda_bounds[0], lambda_bounds[1], log=True)
    non_negative = trial.suggest_categorical('non_negative', [True, False])

    # Fit model
    factors, loss_history, metrics = sparse_cp_decomposition(
        tensor=tensor,
        rank=rank,
        lambda_A=lambda_A,
        lambda_B=lambda_B,
        lambda_C=lambda_C,
        non_negative=non_negative,
        max_iter=50,  # Fewer iterations for faster optimization
        tolerance=1e-4,
        verbose=False
    )

    # Compute evaluation metric
    if eval_metric == 'masked_rmse':
        loss = masked_rmse_loss(tensor, factors)
    elif eval_metric == 'reconstruction_error':
        loss = metrics['reconstruction_error']
    else:
        raise ValueError(f"Unknown evaluation metric: {eval_metric}")

    # Store additional metrics in trial
    trial.set_user_attr('explained_variance', metrics['explained_variance'])
    trial.set_user_attr('sparsity_A', metrics['sparsity_A'])
    trial.set_user_attr('sparsity_B', metrics['sparsity_B'])
    trial.set_user_attr('sparsity_C', metrics['sparsity_C'])
    trial.set_user_attr('final_loss', metrics['final_loss'])

    return loss


def optimize_hyperparameters(
    tensor: np.ndarray,
    n_trials: int = 100,
    rank_range: Tuple[int, int] = (2, 12),
    lambda_bounds: Tuple[float, float] = (1e-4, 1.0),
    eval_metric: str = 'masked_rmse',
    study_name: str = 'cp_optimization',
    storage_path: Optional[str] = None,
    verbose: bool = True
) -> Tuple[Dict, optuna.Study]:
    """Run hyperparameter optimization for sparse CP decomposition.

    Args:
        tensor: Input tensor for optimization
        n_trials: Number of optimization trials
        rank_range: Range for rank sampling
        lambda_bounds: Bounds for L1 penalty sampling
        eval_metric: Evaluation metric
        study_name: Name for the Optuna study
        storage_path: Path for storing study (optional)
        verbose: Verbosity level

    Returns:
        Tuple of (best_params, study)
    """
    logger.info(f"Starting hyperparameter optimization with {n_trials} trials")

    # Set up storage
    if storage_path:
        storage_path = Path(storage_path)
        storage_path.parent.mkdir(parents=True, exist_ok=True)
        storage_name = f"sqlite:///{storage_path}"
    else:
        storage_name = None

    # Create study
    sampler = TPESampler(seed=42)
    study = optuna.create_study(
        study_name=study_name,
        storage=storage_name,
        sampler=sampler,
        direction='minimize',
        load_if_exists=True
    )

    # Run optimization
    study.optimize(
        lambda trial: objective_function(
            trial, tensor, rank_range, lambda_bounds, eval_metric
        ),
        n_trials=n_trials,
        show_progress_bar=verbose
    )

    # Extract best parameters
    best_trial = study.best_trial
    best_params = {
        'rank': best_trial.params['rank'],
        'lambda_A': best_trial.params['lambda_A'],
        'lambda_B': best_trial.params['lambda_B'],
        'lambda_C': best_trial.params['lambda_C'],
        'non_negative': best_trial.params['non_negative'],
        'best_value': best_trial.value,
        'explained_variance': best_trial.user_attrs.get('explained_variance', 0),
        'sparsity_A': best_trial.user_attrs.get('sparsity_A', 0),
        'sparsity_B': best_trial.user_attrs.get('sparsity_B', 0),
        'sparsity_C': best_trial.user_attrs.get('sparsity_C', 0),
    }

    logger.info(f"Best parameters: rank={best_params['rank']}, "
                f"λ_A={best_params['lambda_A']:.6f}, "
                f"λ_B={best_params['lambda_B']:.6f}, "
                f"λ_C={best_params['lambda_C']:.6f}")

    return best_params, study


def save_optimization_results(
    study: optuna.Study,
    output_dir: str,
    best_params: Dict
) -> None:
    """Save optimization results to disk.

    Args:
        study: Optuna study object
        output_dir: Output directory
        best_params: Best parameters dictionary
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Save study object
    with open(output_path / 'optuna_study.pkl', 'wb') as f:
        pickle.dump(study, f)

    # Save best parameters
    with open(output_path / 'best_params.json', 'w') as f:
        json.dump(best_params, f, indent=2)

    # Save trials DataFrame
    trials_df = study.trials_dataframe()
    trials_df.to_csv(output_path / 'trials.csv', index=False)

    # Save top trials summary
    top_trials = []
    for trial in study.trials:
        if trial.state.name == 'COMPLETE':
            top_trials.append({
                'trial': trial.number,
                'rank': trial.params['rank'],
                'lambda_A': trial.params['lambda_A'],
                'lambda_B': trial.params['lambda_B'],
                'lambda_C': trial.params['lambda_C'],
                'non_negative': trial.params['non_negative'],
                'value': trial.value,
                'explained_variance': trial.user_attrs.get('explained_variance', 0),
                'sparsity_A': trial.user_attrs.get('sparsity_A', 0),
                'sparsity_B': trial.user_attrs.get('sparsity_B', 0),
                'sparsity_C': trial.user_attrs.get('sparsity_C', 0),
            })

    # Sort by objective value
    top_trials_df = pd.DataFrame(top_trials).sort_values('value')
    top_trials_df.to_csv(output_path / 'top_trials.csv', index=False)

    logger.info(f"Saved optimization results to {output_path}")


def load_optimization_results(output_dir: str) -> Tuple[Dict, optuna.Study]:
    """Load optimization results from disk.

    Args:
        output_dir: Directory containing saved results

    Returns:
        Tuple of (best_params, study)
    """
    output_path = Path(output_dir)

    # Load best parameters
    with open(output_path / 'best_params.json', 'r') as f:
        best_params = json.load(f)

    # Load study
    with open(output_path / 'optuna_study.pkl', 'rb') as f:
        study = pickle.load(f)

    return best_params, study
