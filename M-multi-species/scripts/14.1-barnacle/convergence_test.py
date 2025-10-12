#!/usr/bin/env python3
"""
Systematic convergence testing for Barnacle SparseCP decomposition.
Tests multiple parameter combinations to find settings that achieve convergence.
"""
import argparse
import json
import os
import subprocess
from datetime import datetime
from typing import Dict, List, Tuple
import pandas as pd
import numpy as np


def find_optimal_rank(input_file: str, base_output_dir: str) -> int:
    """
    Test different rank values (5-35 in increments of 5) to find the optimal rank.
    Uses baseline parameters and returns the rank with best convergence performance.
    """
    print("🔍 STEP 1: Finding optimal rank value...")

    # Define rank values to test
    rank_values = [30, 35, 40, 45, 50, 55, 60]

    # Baseline parameters for rank testing
    baseline_params = {
        'max_iter': 10000,
        'tol': 1e-4,
        'lambda_gene': 0.1,
        'lambda_sample': 0.1,
        'lambda_time': 0.05
    }

    rank_results = []

    for rank in rank_values:
        print(f"\n🧪 Testing rank {rank}...")

        # Create test parameters for this rank
        test_params = {**baseline_params, 'rank': rank}

        # Run single test
        result = run_single_test(test_params, input_file, base_output_dir, rank)

        # Store results with rank info
        rank_results.append({
            'rank': rank,
            'converged': result['converged'],
            'final_loss': result['final_loss'],
            'success': result['success'],
            'output_dir': result['output_dir']
        })

    # Analyze results to find best rank
    # Prioritize convergence, then lowest loss
    converged_results = [r for r in rank_results if r['converged'] and r['success']]

    if converged_results:
        # If multiple converged, choose the one with lowest loss
        best_result = min(converged_results, key=lambda x: x['final_loss'] if x['final_loss'] is not None else float('inf'))
        best_rank = best_result['rank']
        print(f"✅ Best rank found: {best_rank} (converged with loss: {best_result['final_loss']})")
    else:
        # If no convergence, choose the one with lowest loss among successful runs
        successful_results = [r for r in rank_results if r['success']]
        if successful_results:
            best_result = min(successful_results, key=lambda x: x['final_loss'] if x['final_loss'] is not None else float('inf'))
            best_rank = best_result['rank']
            print(f"⚠️  No convergence achieved. Best rank: {best_rank} (lowest loss: {best_result['final_loss']})")
        else:
            # Fallback to middle rank if nothing worked
            best_rank = 20
            print(f"❌ No successful runs. Using fallback rank: {best_rank}")

    print(f"\n📊 Rank testing summary:")
    for result in rank_results:
        status = "✅" if result['converged'] else "❌"
        loss = f" (loss: {result['final_loss']})" if result['final_loss'] else ""
        print(f"  Rank {result['rank']}: {status}{loss}")

    return best_rank


def define_parameter_grid(optimal_rank: int) -> List[Dict]:
    """
    Define a comprehensive grid of parameters to test for convergence for a specific rank.
    Returns a list of parameter combinations to test for the optimal rank.
    """
    print(f"🎯 STEP 2: Optimizing parameters for rank {optimal_rank}...")

    # Test different max_iter values (increased to 25000 as requested)
    max_iter_values = [1000, 2000, 5000, 10000, 15000, 20000, 25000]

    # Test different tolerance values
    tol_values = [1e-5, 1e-4, 1e-3]

    # Test different lambda combinations
    # Focus on gene regularization since that's typically most important for sparsity
    lambda_combinations = [
        # Original values
        {'lambda_gene': 0.1, 'lambda_sample': 0.1, 'lambda_time': 0.05},
        # Increase gene regularization
        {'lambda_gene': 0.2, 'lambda_sample': 0.1, 'lambda_time': 0.05},
        {'lambda_gene': 0.5, 'lambda_sample': 0.1, 'lambda_time': 0.05},
        {'lambda_gene': 1.0, 'lambda_sample': 0.1, 'lambda_time': 0.05},
        # Decrease gene regularization
        {'lambda_gene': 0.05, 'lambda_sample': 0.1, 'lambda_time': 0.05},
        {'lambda_gene': 0.01, 'lambda_sample': 0.1, 'lambda_time': 0.05},
        # Increase sample regularization
        {'lambda_gene': 0.1, 'lambda_sample': 0.2, 'lambda_time': 0.05},
        {'lambda_gene': 0.1, 'lambda_sample': 0.5, 'lambda_time': 0.05},
        # Increase time regularization
        {'lambda_gene': 0.1, 'lambda_sample': 0.1, 'lambda_time': 0.1},
        {'lambda_gene': 0.1, 'lambda_sample': 0.1, 'lambda_time': 0.2},
    ]

    # Generate all combinations for the optimal rank
    param_combinations = []
    for max_iter in max_iter_values:
        for tol in tol_values:
            for lambda_combo in lambda_combinations:
                combo = {
                    'rank': optimal_rank,
                    'max_iter': max_iter,
                    'tol': tol,
                    **lambda_combo
                }
                param_combinations.append(combo)

    # Also add some specific promising combinations for the optimal rank
    promising_combos = [
        # High iterations, relaxed tolerance, moderate regularization
        {'rank': optimal_rank, 'max_iter': 5000, 'tol': 1e-4, 'lambda_gene': 0.1, 'lambda_sample': 0.1, 'lambda_time': 0.05},
        {'rank': optimal_rank, 'max_iter': 10000, 'tol': 1e-4, 'lambda_gene': 0.1, 'lambda_sample': 0.1, 'lambda_time': 0.05},
        # High iterations, very relaxed tolerance, moderate regularization
        {'rank': optimal_rank, 'max_iter': 5000, 'tol': 1e-3, 'lambda_gene': 0.1, 'lambda_sample': 0.1, 'lambda_time': 0.05},
        # Very high iterations, relaxed tolerance, increased regularization
        {'rank': optimal_rank, 'max_iter': 10000, 'tol': 1e-4, 'lambda_gene': 0.5, 'lambda_sample': 0.1, 'lambda_time': 0.05},
        # Maximum iterations test
        {'rank': optimal_rank, 'max_iter': 25000, 'tol': 1e-4, 'lambda_gene': 0.1, 'lambda_sample': 0.1, 'lambda_time': 0.05},
    ]

    param_combinations.extend(promising_combos)

    print(f"Testing {len(param_combinations)} parameter combinations")
    return param_combinations


def run_single_test(params: Dict, input_file: str, base_output_dir: str, test_id: int) -> Dict:
    """
    Run a single parameter combination test.
    Returns test results including convergence status and final loss.
    """
    output_dir = os.path.join(base_output_dir, f"test_{test_id:03d}")

    # Build the command
    cmd = [
        'uv', 'run', 'python',
        '/Users/sr320/GitHub/timeseries_molecular/M-multi-species/scripts/14.1-barnacle/build_tensor_and_run.py',
        '--input-file', input_file,
        '--output-dir', output_dir,
        '--rank', str(params['rank']),
        '--lambda-gene', str(params['lambda_gene']),
        '--lambda-sample', str(params['lambda_sample']),
        '--lambda-time', str(params['lambda_time']),
        '--max-iter', str(params['max_iter']),
        '--tol', str(params['tol']),
        '--seed', '42'
    ]

    print(f"\n{'='*60}")
    print(f"Test {test_id:03d}: {params}")
    print(f"Command: {' '.join(cmd)}")
    print(f"{'='*60}")

    # Run the command and capture output
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=3600  # 1 hour timeout
        )

        # Check for convergence in metadata
        metadata_file = os.path.join(output_dir, 'barnacle_factors', 'metadata.json')
        if os.path.exists(metadata_file):
            with open(metadata_file, 'r') as f:
                metadata = json.load(f)

            converged = metadata.get('model_converged', False)
            final_loss = metadata.get('final_loss', None)
            n_components = metadata.get('n_components', 5)

            test_result = {
                'test_id': test_id,
                'converged': converged,
                'final_loss': final_loss,
                'n_components': n_components,
                'rank': params['rank'],
                'max_iter': params['max_iter'],
                'tol': params['tol'],
                'lambda_gene': params['lambda_gene'],
                'lambda_sample': params['lambda_sample'],
                'lambda_time': params['lambda_time'],
                'return_code': result.returncode,
                'success': result.returncode == 0,
                'output_dir': output_dir
            }

            if converged:
                print(f"✅ CONVERGED! Final loss: {final_loss}")
            else:
                print(f"❌ Did not converge. Final loss: {final_loss}")

            return test_result
        else:
            print("❌ No metadata file found - test may have failed")
            return {
                'test_id': test_id,
                'converged': False,
                'final_loss': None,
                'success': False,
                'return_code': result.returncode,
                'rank': params['rank'],
                'max_iter': params['max_iter'],
                'tol': params['tol'],
                'lambda_gene': params['lambda_gene'],
                'lambda_sample': params['lambda_sample'],
                'lambda_time': params['lambda_time'],
                'output_dir': output_dir
            }

    except subprocess.TimeoutExpired:
        print("❌ Test timed out after 1 hour")
        return {
            'test_id': test_id,
            'converged': False,
            'final_loss': None,
            'success': False,
            'return_code': -1,
            'rank': params['rank'],
            'max_iter': params['max_iter'],
            'tol': params['tol'],
            'lambda_gene': params['lambda_gene'],
            'lambda_sample': params['lambda_sample'],
            'lambda_time': params['lambda_time'],
            'output_dir': output_dir
        }


def save_results(results: List[Dict], output_file: str):
    """Save test results to CSV and JSON files."""
    # Convert to DataFrame for CSV
    df = pd.DataFrame(results)

    # Save as CSV
    df.to_csv(output_file, index=False)

    # Save as JSON for detailed analysis
    json_file = output_file.replace('.csv', '_detailed.json')
    with open(json_file, 'w') as f:
        json.dump(results, f, indent=2, default=str)

    print(f"\n📊 Results saved to: {output_file}")
    print(f"📋 Detailed results: {json_file}")


def print_summary(results: List[Dict]):
    """Print a summary of test results."""
    converged_tests = [r for r in results if r['converged']]
    successful_tests = [r for r in results if r['success']]

    print(f"\n{'='*80}")
    print("TEST SUMMARY")
    print(f"{'='*80}")
    print(f"Total tests run: {len(results)}")
    print(f"Successful runs: {len(successful_tests)}")
    print(f"Converged tests: {len(converged_tests)}")
    print(f"Success rate: {len(successful_tests)/len(results)*100:.1f}%")
    print(f"Convergence rate: {len(converged_tests)/len(successful_tests)*100:.1f}%" if successful_tests else "Convergence rate: 0%")

    if converged_tests:
        print("\n🏆 CONVERGED PARAMETER COMBINATIONS:")
        for test in converged_tests:
            print(f"  Test {test['test_id']:03d}: rank={test['rank']}, max_iter={test['max_iter']}, tol={test['tol']}, "
                  f"λ_gene={test['lambda_gene']}, λ_sample={test['lambda_sample']}, λ_time={test['lambda_time']}")
            print(f"    Final loss: {test['final_loss']}")
            print(f"    Output dir: {test['output_dir']}")
    else:
        print("\n❌ No parameter combinations achieved convergence.")
        print("💡 Consider trying even more relaxed tolerances or higher iterations.")


def main():
    parser = argparse.ArgumentParser(description='Two-step Barnacle convergence testing: rank optimization then parameter optimization')
    parser.add_argument('--input-file', required=True, help='Path to merged vst_counts_matrix.csv file')
    parser.add_argument('--output-dir', required=True, help='Base directory for test outputs')
    parser.add_argument('--results-file', default='convergence_test_results.csv',
                       help='Output file for results summary')
    args = parser.parse_args()

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    print(f"🚀 Starting two-step convergence testing for Barnacle")
    print(f"📁 Input file: {args.input_file}")
    print(f"📁 Output directory: {args.output_dir}")

    # Step 1: Find optimal rank
    optimal_rank = find_optimal_rank(args.input_file, args.output_dir)

    # Step 2: Optimize parameters for the optimal rank
    param_combinations = define_parameter_grid(optimal_rank)

    # Run parameter optimization tests
    results = []
    for i, params in enumerate(param_combinations):
        result = run_single_test(params, args.input_file, args.output_dir, i + 1)
        results.append(result)

        # Save intermediate results every 10 tests
        if (i + 1) % 10 == 0:
            intermediate_file = os.path.join(args.output_dir, f'intermediate_results_{i+1}.csv')
            save_results(results, intermediate_file)

    # Save final results
    results_file = os.path.join(args.output_dir, args.results_file)
    save_results(results, results_file)

    # Print summary
    print_summary(results)


if __name__ == '__main__':
    main()
