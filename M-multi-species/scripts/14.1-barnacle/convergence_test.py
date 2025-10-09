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


def define_parameter_grid() -> List[Dict]:
    """
    Define a comprehensive grid of parameters to test for convergence.
    Returns a list of parameter combinations to test.
    """
    # Current baseline parameters (from original script)
    baseline = {
        'max_iter': 1000,
        'tol': 1e-5,
        'lambda_gene': 0.1,
        'lambda_sample': 0.1,
        'lambda_time': 0.05
    }

    # Test different max_iter values
    max_iter_values = [1000, 2000, 5000, 10000]

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

    # Generate all combinations
    param_combinations = []
    for max_iter in max_iter_values:
        for tol in tol_values:
            for lambda_combo in lambda_combinations:
                combo = {
                    'max_iter': max_iter,
                    'tol': tol,
                    **lambda_combo
                }
                param_combinations.append(combo)

    # Also add some specific promising combinations
    promising_combos = [
        # High iterations, relaxed tolerance, moderate regularization
        {'max_iter': 5000, 'tol': 1e-4, 'lambda_gene': 0.1, 'lambda_sample': 0.1, 'lambda_time': 0.05},
        {'max_iter': 10000, 'tol': 1e-4, 'lambda_gene': 0.1, 'lambda_sample': 0.1, 'lambda_time': 0.05},
        # High iterations, very relaxed tolerance, moderate regularization
        {'max_iter': 5000, 'tol': 1e-3, 'lambda_gene': 0.1, 'lambda_sample': 0.1, 'lambda_time': 0.05},
        # Very high iterations, relaxed tolerance, increased regularization
        {'max_iter': 10000, 'tol': 1e-4, 'lambda_gene': 0.5, 'lambda_sample': 0.1, 'lambda_time': 0.05},
    ]

    param_combinations.extend(promising_combos)

    print(f"Testing {len(param_combinations)} parameter combinations")
    return param_combinations


def run_single_test(params: Dict, input_dir: str, base_output_dir: str, test_id: int) -> Dict:
    """
    Run a single parameter combination test.
    Returns test results including convergence status and final loss.
    """
    output_dir = os.path.join(base_output_dir, f"test_{test_id:03d}")

    # Build the command
    cmd = [
        'uv', 'run', 'python',
        '/Users/sr320/GitHub/timeseries_molecular/M-multi-species/scripts/14.1-barnacle/build_tensor_and_run.py',
        '--input-dir', input_dir,
        '--output-dir', output_dir,
        '--rank', '5',
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
            print(f"  Test {test['test_id']:03d}: max_iter={test['max_iter']}, tol={test['tol']}, "
                  f"λ_gene={test['lambda_gene']}, λ_sample={test['lambda_sample']}, λ_time={test['lambda_time']}")
            print(f"    Final loss: {test['final_loss']}")
            print(f"    Output dir: {test['output_dir']}")
    else:
        print("\n❌ No parameter combinations achieved convergence.")
        print("💡 Consider trying even more relaxed tolerances or higher iterations.")


def main():
    parser = argparse.ArgumentParser(description='Systematically test Barnacle convergence parameters')
    parser.add_argument('--input-dir', required=True, help='Directory with normalized expression CSV files')
    parser.add_argument('--output-dir', required=True, help='Base directory for test outputs')
    parser.add_argument('--results-file', default='convergence_test_results.csv',
                       help='Output file for results summary')
    args = parser.parse_args()

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Define parameter grid
    param_combinations = define_parameter_grid()

    # Run tests
    results = []
    for i, params in enumerate(param_combinations):
        result = run_single_test(params, args.input_dir, args.output_dir, i + 1)
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
