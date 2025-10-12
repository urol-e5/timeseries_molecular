"""Sparse CP decomposition model with L1 regularization."""

import logging
from typing import Dict, List, Optional, Tuple, Union
import time

import numpy as np
from scipy import sparse
import tensorly as tl
from tensorly.decomposition import CP

logger = logging.getLogger(__name__)


def soft_thresholding(x: np.ndarray, threshold: float) -> np.ndarray:
    """Apply soft thresholding for L1 regularization.

    Args:
        x: Input array
        threshold: Threshold value

    Returns:
        Soft-thresholded array
    """
    return np.sign(x) * np.maximum(np.abs(x) - threshold, 0)


def project_nonnegative(x: np.ndarray) -> np.ndarray:
    """Project onto non-negative orthant.

    Args:
        x: Input array

    Returns:
        Non-negative projection
    """
    return np.maximum(x, 0)


def sparse_cp_als(
    tensor: np.ndarray,
    rank: int,
    lambda_A: float = 0.1,
    lambda_B: float = 0.1,
    lambda_C: float = 0.1,
    non_negative: bool = False,
    max_iter: int = 100,
    tolerance: float = 1e-6,
    verbose: bool = True
) -> Tuple[List[np.ndarray], List[float]]:
    """Fit sparse CP decomposition using ALS with L1 regularization.

    Args:
        tensor: Input tensor X of shape (n_genes, n_species, n_timepoints)
        rank: Number of components
        lambda_A: L1 penalty for gene factors
        lambda_B: L1 penalty for species factors
        lambda_C: L1 penalty for time factors
        non_negative: Whether to enforce non-negativity
        max_iter: Maximum number of iterations
        tolerance: Convergence tolerance
        verbose: Whether to print progress

    Returns:
        Tuple of (factors, loss_history)
    """
    logger.info(f"Fitting sparse CP decomposition with rank={rank}")
    logger.info(f"L1 penalties: λ_A={lambda_A}, λ_B={lambda_B}, λ_C={lambda_C}")
    logger.info(f"Non-negative: {non_negative}")

    # Initialize factors
    n_genes, n_species, n_timepoints = tensor.shape

    # Random initialization
    np.random.seed(42)
    A = np.random.randn(n_genes, rank)
    B = np.random.randn(n_species, rank)
    C = np.random.randn(n_timepoints, rank)

    # Normalize factors
    A = A / np.linalg.norm(A, axis=0, keepdims=True)
    B = B / np.linalg.norm(B, axis=0, keepdims=True)
    C = C / np.linalg.norm(C, axis=0, keepdims=True)

    loss_history = []
    start_time = time.time()

    for iteration in range(max_iter):
        # Update A (gene factors)
        for r in range(rank):
            # Khatri-Rao product of B and C
            BC = np.kron(B[:, r], C[:, r])

            # Update rule for A[:, r]
            numerator = np.sum(tensor.reshape(n_genes, -1) * BC.reshape(1, -1), axis=1)
            denominator = np.sum(BC.reshape(1, -1) ** 2, axis=1) + lambda_A

            A[:, r] = numerator / (denominator + 1e-8)

        # Update B (species factors)
        for r in range(rank):
            # Khatri-Rao product of A and C
            AC = np.kron(A[:, r], C[:, r])

            # Update rule for B[:, r]
            numerator = np.sum(tensor.transpose(1, 0, 2).reshape(n_species, -1) * AC.reshape(1, -1), axis=1)
            denominator = np.sum(AC.reshape(1, -1) ** 2, axis=1) + lambda_B

            B[:, r] = numerator / (denominator + 1e-8)

        # Update C (time factors)
        for r in range(rank):
            # Khatri-Rao product of A and B
            AB = np.kron(A[:, r], B[:, r])

            # Update rule for C[:, r]
            numerator = np.sum(tensor.transpose(2, 0, 1).reshape(n_timepoints, -1) * AB.reshape(1, -1), axis=1)
            denominator = np.sum(AB.reshape(1, -1) ** 2, axis=1) + lambda_C

            C[:, r] = numerator / (denominator + 1e-8)

        # Apply soft thresholding for sparsity
        A = soft_thresholding(A, lambda_A)
        B = soft_thresholding(B, lambda_B)
        C = soft_thresholding(C, lambda_C)

        # Apply non-negativity constraint if requested
        if non_negative:
            A = project_nonnegative(A)
            B = project_nonnegative(B)
            C = project_nonnegative(C)

        # Re-normalize factors
        for r in range(rank):
            norm_factor = np.linalg.norm(A[:, r]) * np.linalg.norm(B[:, r]) * np.linalg.norm(C[:, r])
            if norm_factor > 0:
                A[:, r] = A[:, r] / np.linalg.norm(A[:, r])
                B[:, r] = B[:, r] / np.linalg.norm(B[:, r])
                C[:, r] = C[:, r] / np.linalg.norm(C[:, r]) * norm_factor ** (1/3)

        # Compute reconstruction and loss
        reconstruction = tl.cp_to_tensor((None, [A, B, C]))
        loss = 0.5 * np.linalg.norm(tensor - reconstruction) ** 2

        # Add L1 penalties to loss
        loss += lambda_A * np.sum(np.abs(A))
        loss += lambda_B * np.sum(np.abs(B))
        loss += lambda_C * np.sum(np.abs(C))

        loss_history.append(loss)

        # Check convergence
        if iteration > 0 and loss_history[-2] != 0:
            rel_change = abs(loss_history[-1] - loss_history[-2]) / abs(loss_history[-2])
            if rel_change < tolerance:
                if verbose:
                    logger.info(f"Converged at iteration {iteration} (rel_change={rel_change:.6f})")
                break

        if verbose and iteration % 10 == 0:
            logger.info(f"Iteration {iteration}: loss={loss:.6f}")

    total_time = time.time() - start_time
    logger.info(f"CP decomposition completed in {total_time:.2f}s")

    return [A, B, C], loss_history


def sparse_cp_decomposition(
    tensor: np.ndarray,
    rank: int,
    lambda_A: float = 0.1,
    lambda_B: float = 0.1,
    lambda_C: float = 0.1,
    non_negative: bool = False,
    method: str = 'als',
    max_iter: int = 100,
    tolerance: float = 1e-6,
    verbose: bool = True
) -> Tuple[List[np.ndarray], List[float], Dict]:
    """Main interface for sparse CP decomposition.

    Args:
        tensor: Input tensor
        rank: Number of components
        lambda_A: L1 penalty for gene factors
        lambda_B: L1 penalty for species factors
        lambda_C: L1 penalty for time factors
        non_negative: Whether to enforce non-negativity
        method: Decomposition method ('als' only for now)
        max_iter: Maximum iterations
        tolerance: Convergence tolerance
        verbose: Verbosity level

    Returns:
        Tuple of (factors, loss_history, metrics)
    """
    if method != 'als':
        raise ValueError(f"Method {method} not supported yet")

    # Fit the model
    factors, loss_history = sparse_cp_als(
        tensor=tensor,
        rank=rank,
        lambda_A=lambda_A,
        lambda_B=lambda_B,
        lambda_C=lambda_C,
        non_negative=non_negative,
        max_iter=max_iter,
        tolerance=tolerance,
        verbose=verbose
    )

    # Compute metrics
    reconstruction = tl.cp_to_tensor((None, factors))
    reconstruction_error = np.linalg.norm(tensor - reconstruction) ** 2 / 2

    # Explained variance
    tensor_norm = np.linalg.norm(tensor) ** 2
    if tensor_norm > 0:
        explained_variance = 1 - (reconstruction_error / (tensor_norm / 2))
    else:
        explained_variance = 0.0

    # Sparsity metrics
    sparsity_A = np.mean(factors[0] == 0)
    sparsity_B = np.mean(factors[1] == 0)
    sparsity_C = np.mean(factors[2] == 0)

    metrics = {
        'reconstruction_error': reconstruction_error,
        'explained_variance': explained_variance,
        'final_loss': loss_history[-1] if loss_history else float('inf'),
        'n_iterations': len(loss_history),
        'sparsity_A': sparsity_A,
        'sparsity_B': sparsity_B,
        'sparsity_C': sparsity_C,
        'converged': len(loss_history) < max_iter if loss_history else False
    }

    logger.info(f"Final metrics: R²={explained_variance:.4f}, sparsity_A={sparsity_A:.4f}")

    return factors, loss_history, metrics
