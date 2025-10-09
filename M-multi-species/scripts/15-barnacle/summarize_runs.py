#!/usr/bin/env python3
import argparse
import json
import os
from typing import Dict, List

import pandas as pd


def load_metadata(run_dir: str) -> Dict:
    meta_path = os.path.join(run_dir, 'barnacle_factors', 'metadata.json')
    if not os.path.exists(meta_path):
        return {}
    with open(meta_path) as fh:
        return json.load(fh)


def sparsity(path: str) -> float:
    if not os.path.exists(path):
        return float('nan')
    df = pd.read_csv(path, index_col=0)
    total = df.size
    zeros = (df == 0).sum().sum()
    return float(zeros) / float(total) if total > 0 else float('nan')


def main() -> None:
    parser = argparse.ArgumentParser(description='Summarize multiple Barnacle runs into SUMMARY.md')
    parser.add_argument('--base-dir', required=True, help='Base output dir containing rank-* subdirs')
    args = parser.parse_args()

    base = args.base_dir
    runs = []
    for name in sorted(os.listdir(base)):
        if not name.startswith('rank-'):
            continue
        run_dir = os.path.join(base, name)
        if not os.path.isdir(run_dir):
            continue
        meta = load_metadata(run_dir)
        rank = int(name.split('-', 1)[1])
        gene_path = os.path.join(run_dir, 'barnacle_factors', 'gene_factors.csv')
        sample_path = os.path.join(run_dir, 'barnacle_factors', 'sample_factors.csv')
        time_path = os.path.join(run_dir, 'barnacle_factors', 'time_factors.csv')

        runs.append({
            'rank': rank,
            'tensor_shape': ' × '.join(map(str, meta.get('tensor_shape', []))) if meta else '',
            'components': meta.get('n_components', ''),
            'converged': meta.get('model_converged', ''),
            'final_loss': meta.get('final_loss', ''),
            'gene_sparsity': sparsity(gene_path),
            'sample_sparsity': sparsity(sample_path),
            'time_sparsity': sparsity(time_path),
            'component_weights_png': os.path.relpath(os.path.join(run_dir, 'figures', 'component_weights.png'), base),
            'time_loadings_png': os.path.relpath(os.path.join(run_dir, 'figures', 'time_loadings.png'), base),
        })

    runs = sorted(runs, key=lambda x: x['rank'])

    lines: List[str] = []
    lines.append('# Barnacle multi-rank comparison (15-barnacle)')
    lines.append('')
    if runs and runs[0]['tensor_shape']:
        lines.append(f"- Tensor shape (genes × samples × time): {runs[0]['tensor_shape']}")
    ranks_str = ', '.join(str(r['rank']) for r in runs)
    lines.append(f"- Ranks evaluated: {ranks_str}")
    lines.append('')

    # Table
    lines.append('| Rank | Components | Converged | Final loss | Gene sparsity | Sample sparsity | Time sparsity |')
    lines.append('|---:|---:|:---:|---:|---:|---:|---:|')
    for r in runs:
        lines.append('| {rank} | {components} | {converged} | {final_loss} | {gs:.3f} | {ss:.3f} | {ts:.3f} |'.format(
            rank=r['rank'], components=r['components'], converged=r['converged'], final_loss=r['final_loss'],
            gs=r['gene_sparsity'], ss=r['sample_sparsity'], ts=r['time_sparsity']
        ))
    lines.append('')

    # Figures per run
    for r in runs:
        lines.append(f"## Rank {r['rank']}")
        lines.append('')
        lines.append(f"![Component weights]({r['component_weights_png']})")
        lines.append('')
        lines.append(f"![Timepoint loadings]({r['time_loadings_png']})")
        lines.append('')

    out_path = os.path.join(base, 'SUMMARY.md')
    with open(out_path, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()


