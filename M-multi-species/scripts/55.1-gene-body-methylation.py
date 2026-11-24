#!/usr/bin/env python3
"""Aggregate per-CpG methylation to gene bodies.

Reads a merged CpG x sample matrix of percent methylation (rows named
CpG_<scaffold>_<position>, empty cells where a sample lacks coverage) and a GFF
with `gene` features. For every gene and sample it reports the mean percent
methylation over the CpGs inside the gene that have a value in that sample.
Cells supported by fewer than --min-cpgs CpGs are left empty (NOT set to zero).

Outputs, in --out-dir:
  <species>-gene-body-methylation.csv.gz   gene_id x sample, mean % methylation
  <species>-gene-body-CpG-counts.csv.gz    gene_id x sample, CpGs behind each mean

The matrix may be a local path or an http(s) URL, which is streamed.
Called from 55-variance-partitioning.Rmd.
"""
import argparse
import bisect
import csv
import gzip
import os
import subprocess
import sys

import numpy as np


def read_genes(gff):
    """Return {scaffold: (starts, ends, prefix_max_end, gene_index)} and gene ids."""
    per_scaffold, ids = {}, []
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            gid = None
            for attr in f[8].split(";"):
                if attr.startswith("ID="):
                    gid = attr[3:]
                    break
            if gid is None:
                continue
            per_scaffold.setdefault(f[0], []).append((int(f[3]), int(f[4]), len(ids)))
            ids.append(gid)
    index = {}
    for scaffold, genes in per_scaffold.items():
        genes.sort()
        starts = [g[0] for g in genes]
        ends = [g[1] for g in genes]
        prefix_max, running = [], 0
        for e in ends:
            running = max(running, e)
            prefix_max.append(running)
        index[scaffold] = (starts, ends, prefix_max, [g[2] for g in genes])
    return index, ids


def genes_at(index, scaffold, pos):
    """Indices of every gene whose span contains pos (genes can overlap)."""
    entry = index.get(scaffold)
    if entry is None:
        return ()
    starts, ends, prefix_max, gidx = entry
    hits = []
    j = bisect.bisect_right(starts, pos) - 1
    while j >= 0 and prefix_max[j] >= pos:
        if ends[j] >= pos:
            hits.append(gidx[j])
        j -= 1
    return hits


def open_matrix(src):
    if src.startswith("http://") or src.startswith("https://"):
        proc = subprocess.Popen(["curl", "-sfL", src], stdout=subprocess.PIPE, text=True)
        return proc.stdout, proc
    return open(src), None


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--species", required=True)
    ap.add_argument("--cpg", required=True, help="merged CpG matrix, path or URL")
    ap.add_argument("--gff", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--min-cpgs", type=int, default=10)
    args = ap.parse_args()

    index, gene_ids = read_genes(args.gff)
    fh, proc = open_matrix(args.cpg)
    samples = [s.strip('"') for s in fh.readline().rstrip("\n").split(",")[1:]]
    sums = np.zeros((len(gene_ids), len(samples)))
    counts = np.zeros((len(gene_ids), len(samples)), dtype=np.int32)

    n_rows = n_genic = 0
    for line in fh:
        n_rows += 1
        cpg_id, rest = line.split(",", 1)
        scaffold, pos = cpg_id.strip('"')[4:].rsplit("_", 1)   # drop "CpG_"
        hits = genes_at(index, scaffold, int(pos))
        if not hits:
            continue
        n_genic += 1
        vals = np.array([float(v) if v else np.nan for v in rest.rstrip("\n").split(",")])
        ok = ~np.isnan(vals)
        add = np.where(ok, vals, 0.0)
        for g in hits:
            sums[g] += add
            counts[g] += ok
    if proc is not None and proc.wait() != 0:
        sys.exit(f"download failed: {args.cpg}")
    if n_rows == 0:
        sys.exit("no CpG rows read")

    with np.errstate(invalid="ignore", divide="ignore"):
        means = np.where(counts >= args.min_cpgs, sums / counts, np.nan)
    keep = np.where(~np.isnan(means).all(axis=1))[0]

    os.makedirs(args.out_dir, exist_ok=True)
    base = os.path.join(args.out_dir, f"{args.species}-gene-body-")
    with gzip.open(base + "methylation.csv.gz", "wt", newline="") as m_out, \
         gzip.open(base + "CpG-counts.csv.gz", "wt", newline="") as c_out:
        mw, cw = csv.writer(m_out), csv.writer(c_out)
        mw.writerow(["gene_id"] + samples)
        cw.writerow(["gene_id"] + samples)
        for g in keep:
            mw.writerow([gene_ids[g]] + ["" if np.isnan(v) else f"{v:.4f}" for v in means[g]])
            cw.writerow([gene_ids[g]] + counts[g].tolist())

    complete = int((~np.isnan(means[keep])).all(axis=1).sum())
    print(f"{args.species}: {n_rows:,} CpGs read, {n_genic:,} in gene bodies; "
          f"{len(keep):,} of {len(gene_ids):,} genes have >= {args.min_cpgs} CpGs in at "
          f"least one sample, {complete:,} in all {len(samples)} samples")


if __name__ == "__main__":
    main()
