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
from typing import Dict, List
import pandas as pd
import numpy as np

def find_optimal_rank(input_file: str, base_output_dir: str) -> int:
    print("\U0001F50D STEP 1: Finding optimal rank value...")
    rank_values = [30, 35, 40, 45, 50, 55, 60]
    baseline_params = {'max_iter': 10000, 'tol': 1e-4, 'lambda_gene': 0.1, 'lambda_sample': 0.1, 'lambda_time': 0.05}
    rank_results = []

    for i, rank in enumerate(rank_values):
        print(f"\n\U0001F9EA Testing rank {rank}...")
        test_params = {**baseline_params, 'rank': rank}
        result = run_single_test(test_params, input_file, base_output_dir, i + 1)
        rank_results.append({**result, 'rank': rank})

    converged = [r for r in rank_results if r['converged'] and r['success']]
    if converged:
        best_result = min(converged, key=lambda x: x['final_loss'] if x['final_loss'] is not None else float('inf'))
    else:
        successful = [r for r in rank_results if r['success']]
        if successful:
            best_result = min(successful, key=lambda x: x['final_loss'] if x['final_loss'] is not None else float('inf'))
        else:
            best_result = {'rank': int(np.median(rank_values))}
            print(f"\u274C No successful runs. Using fallback rank: {best_result['rank']}")
            return best_result['rank']

    print(f"\n\U0001F4CA Rank testing summary:")
    for r in rank_results:
        status = "\u2705" if r['converged'] else "\u274C"
        loss = f" (loss: {r['final_loss']})" if r['final_loss'] is not None else ""
        print(f"  Rank {r['rank']}: {status}{loss}")

    print(f"\u2705 Best rank found: {best_result['rank']} (loss: {best_result.get('final_loss')})")
    return best_result['rank']

def define_parameter_grid(optimal_rank: int) -> List[Dict]:
    print(f"\U0001F3AF STEP 2: Optimizing parameters for rank {optimal_rank}...")
    max_iter_values = [1000, 2000, 5000, 10000, 15000, 20000, 25000]
    tol_values = [1e-5, 1e-4, 1e-3]
    lambda_combinations = [
        {'lambda_gene': 0.1, 'lambda_sample': 0.1, 'lambda_time': 0.05},
        {'lambda_gene': 0.2, 'lambda_sample': 0.1, 'lambda_time': 0.05},
        {'lambda_gene': 0.5, 'lambda_sample': 0.1, 'lambda_time': 0.05},
        {'lambda_gene': 1.0, 'lambda_sample': 0.1, 'lambda_time': 0.05},
        {'lambda_gene': 0.05, 'lambda_sample': 0.1, 'lambda_time': 0.05},
        {'lambda_gene': 0.01, 'lambda_sample': 0.1, 'lambda_time': 0.05},
        {'lambda_gene': 0.1, 'lambda_sample': 0.2, 'lambda_time': 0.05},
        {'lambda_gene': 0.1, 'lambda_sample': 0.5, 'lambda_time': 0.05},
        {'lambda_gene': 0.1, 'lambda_sample': 0.1, 'lambda_time': 0.1},
        {'lambda_gene': 0.1, 'lambda_sample': 0.1, 'lambda_time': 0.2},
    ]
    param_combinations = [
        {**l, 'rank': optimal_rank, 'max_iter': mi, 'tol': tol}
        for mi in max_iter_values for tol in tol_values for l in lambda_combinations
    ]
    return param_combinations

def run_single_test(params: Dict, input_file: str, base_output_dir: str, test_id: int) -> Dict:
    output_dir = os.path.join(base_output_dir, f"test_{test_id:03d}")
    script_path = os.environ.get("BARNACLE_SCRIPT_PATH") or "../scripts/14.1-barnacle/build_tensor_and_run.py"
    cmd = [
        '/Users/sr320/.local/bin/uv', 'run', 'python', script_path,
        '--input-file', input_file, '--output-dir', output_dir,
        '--rank', str(params['rank']), '--lambda-gene', str(params['lambda_gene']),
        '--lambda-sample', str(params['lambda_sample']), '--lambda-time', str(params['lambda_time']),
        '--max-iter', str(params['max_iter']), '--tol', str(params['tol']), '--seed', '42']

    print(f"\n{'='*60}\nTest {test_id:03d}: {params}\nCommand: {' '.join(cmd)}\n{'='*60}")

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        if result.returncode != 0:
            print(f"\u274C Subprocess failed (code {result.returncode})\n{result.stderr}")

        metadata_path = os.path.join(output_dir, 'barnacle_factors', 'metadata.json')
        if os.path.exists(metadata_path):
            with open(metadata_path) as f:
                metadata = json.load(f)
            return {
                'test_id': test_id,
                'converged': metadata.get('model_converged', False),
                'final_loss': metadata.get('final_loss'),
                'success': result.returncode == 0,
                'output_dir': output_dir,
                **params
            }
        else:
            print("\u274C No metadata file found")
            return {'test_id': test_id, 'converged': False, 'final_loss': None, 'success': False, 'output_dir': output_dir, **params}
    except subprocess.TimeoutExpired:
        print("\u274C Timeout after 1 hour")
        return {'test_id': test_id, 'converged': False, 'final_loss': None, 'success': False, 'output_dir': output_dir, **params}

def save_results(results: List[Dict], output_file: str):
    df = pd.DataFrame(results)
    df.to_csv(output_file, index=False)
    json_file = output_file.replace('.csv', f'_{datetime.now():%Y%m%d_%H%M}.json')
    with open(json_file, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\n\U0001F4CA Results saved: {output_file}\n\U0001F4CB Detailed JSON: {json_file}")

def print_summary(results: List[Dict]):
    print(f"\n{'='*80}\nTEST SUMMARY\n{'='*80}")
    if not results:
        print("No tests were run.")
        return
    converged = [r for r in results if r['converged']]
    successful = [r for r in results if r['success']]
    print(f"Total: {len(results)} | Successful: {len(successful)} | Converged: {len(converged)}")
    print(f"Success rate: {len(successful)/len(results)*100:.1f}%")
    print(f"Convergence rate: {len(converged)/len(successful)*100:.1f}%" if successful else "Convergence rate: 0%")
    if converged:
        print("\n\U0001F3C6 CONVERGED TESTS:")
        for r in converged:
            print(f"  Test {r['test_id']:03d}: rank={r['rank']}, iter={r['max_iter']}, tol={r['tol']}, λ_gene={r['lambda_gene']}, λ_sample={r['lambda_sample']}, λ_time={r['lambda_time']} | Loss: {r['final_loss']}")
    else:
        print("\n❌ No tests converged. Consider looser tolerance or more iterations.")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-file', required=True)
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--results-file', default='convergence_test_results.csv')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    print(f"\U0001F680 Starting Barnacle Convergence Test")
    print(f"Input: {args.input_file}\nOutput: {args.output_dir}")

    best_rank = find_optimal_rank(args.input_file, args.output_dir)
    param_grid = define_parameter_grid(best_rank)

    results = []
    for i, p in enumerate(param_grid):
        result = run_single_test(p, args.input_file, args.output_dir, i + 100)
        results.append(result)
        if (i + 1) % 10 == 0:
            save_results(results, os.path.join(args.output_dir, f"intermediate_{i+1}.csv"))

    save_results(results, os.path.join(args.output_dir, args.results_file))
    print_summary(results)

if __name__ == '__main__':
    main()
