#!/bin/bash
# Full orthology analysis for coral species
# This script runs the complete analysis on full protein datasets
# WARNING: This will take several hours to complete

set -e  # Exit on any error

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_DIR="${SCRIPT_DIR}/../output/11-orthology-analysis"
APUL_PROTEINS="${SCRIPT_DIR}/../../D-Apul/data/Apulchra-genome.pep.faa"
PEVE_PROTEINS="${SCRIPT_DIR}/../../E-Peve/data/Porites_evermanni_v1.annot.pep.fa"
PTUA_PROTEINS="${SCRIPT_DIR}/../../F-Ptua/data/Pocillopora_meandrina_HIv1.genes.pep.faa"

# BLAST parameters
EVALUE="1e-5"
THREADS="8"
OUTFMT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"

echo "=== Coral Orthology Analysis - Full Dataset ==="
echo "Started at: $(date)"
echo ""

# Check input files
echo "Checking input files..."
for file in "$APUL_PROTEINS" "$PEVE_PROTEINS" "$PTUA_PROTEINS"; do
    if [[ ! -f "$file" ]]; then
        echo "ERROR: Input file not found: $file"
        exit 1
    fi
    echo "✓ $(basename "$file"): $(grep -c '^>' "$file") proteins"
done
echo ""

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Create BLAST databases
echo "Creating BLAST databases..."
makeblastdb -in "$APUL_PROTEINS" -dbtype prot -out "$OUTPUT_DIR/Apul_proteins_full" -title "Acropora_pulchra_proteins"
makeblastdb -in "$PEVE_PROTEINS" -dbtype prot -out "$OUTPUT_DIR/Peve_proteins_full" -title "Porites_evermanni_proteins"
makeblastdb -in "$PTUA_PROTEINS" -dbtype prot -out "$OUTPUT_DIR/Ptua_proteins_full" -title "Pocillopora_tuahiniensis_proteins"
echo "✓ BLAST databases created"
echo ""

# Run all-vs-all BLAST comparisons
echo "Running BLAST comparisons (this will take several hours)..."

echo "  Apul vs Peve..."
blastp -query "$APUL_PROTEINS" -db "$OUTPUT_DIR/Peve_proteins_full" \
  -out "$OUTPUT_DIR/Apul_vs_Peve_full.blastp" \
  -evalue "$EVALUE" -num_threads "$THREADS" -outfmt "$OUTFMT" -max_target_seqs 1

echo "  Peve vs Apul..."
blastp -query "$PEVE_PROTEINS" -db "$OUTPUT_DIR/Apul_proteins_full" \
  -out "$OUTPUT_DIR/Peve_vs_Apul_full.blastp" \
  -evalue "$EVALUE" -num_threads "$THREADS" -outfmt "$OUTFMT" -max_target_seqs 1

echo "  Apul vs Ptua..."
blastp -query "$APUL_PROTEINS" -db "$OUTPUT_DIR/Ptua_proteins_full" \
  -out "$OUTPUT_DIR/Apul_vs_Ptua_full.blastp" \
  -evalue "$EVALUE" -num_threads "$THREADS" -outfmt "$OUTFMT" -max_target_seqs 1

echo "  Ptua vs Apul..."
blastp -query "$PTUA_PROTEINS" -db "$OUTPUT_DIR/Apul_proteins_full" \
  -out "$OUTPUT_DIR/Ptua_vs_Apul_full.blastp" \
  -evalue "$EVALUE" -num_threads "$THREADS" -outfmt "$OUTFMT" -max_target_seqs 1

echo "  Peve vs Ptua..."
blastp -query "$PEVE_PROTEINS" -db "$OUTPUT_DIR/Ptua_proteins_full" \
  -out "$OUTPUT_DIR/Peve_vs_Ptua_full.blastp" \
  -evalue "$EVALUE" -num_threads "$THREADS" -outfmt "$OUTFMT" -max_target_seqs 1

echo "  Ptua vs Peve..."
blastp -query "$PTUA_PROTEINS" -db "$OUTPUT_DIR/Peve_proteins_full" \
  -out "$OUTPUT_DIR/Ptua_vs_Peve_full.blastp" \
  -evalue "$EVALUE" -num_threads "$THREADS" -outfmt "$OUTFMT" -max_target_seqs 1

echo "✓ All BLAST comparisons completed"
echo ""

# Check BLAST results
echo "BLAST results summary:"
for file in "$OUTPUT_DIR"/*_full.blastp; do
    if [[ -f "$file" ]]; then
        echo "  $(basename "$file"): $(wc -l < "$file") hits"
    fi
done
echo ""

# Create full orthology analysis script
cat > "$OUTPUT_DIR/analyze_full_orthologs.py" << 'EOF'
#!/usr/bin/env python3
"""
Full orthology analysis for coral species.
Processes the complete BLAST results to identify orthologous proteins.
"""

import pandas as pd
import sys
import os

def load_blast_results(blast_file):
    """Load BLAST results with proper column names."""
    columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
               'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen']
    
    if not os.path.exists(blast_file):
        print(f"Warning: {blast_file} not found")
        return pd.DataFrame(columns=columns)
    
    try:
        df = pd.read_csv(blast_file, sep='\t', header=None, names=columns)
        return df
    except Exception as e:
        print(f"Error reading {blast_file}: {e}")
        return pd.DataFrame(columns=columns)

def find_reciprocal_best_hits(blast1, blast2, min_identity=30, min_coverage=50):
    """Find reciprocal best hits between two BLAST results."""
    
    def get_best_hits(df):
        if df.empty:
            return df
        
        # Calculate coverage
        df['qcoverage'] = (df['qend'] - df['qstart'] + 1) / df['qlen'] * 100
        df['scoverage'] = (df['send'] - df['sstart'] + 1) / df['slen'] * 100
        
        # Filter by identity and coverage
        filtered = df[
            (df['pident'] >= min_identity) & 
            ((df['qcoverage'] >= min_coverage) | (df['scoverage'] >= min_coverage))
        ]
        
        if filtered.empty:
            return filtered
        
        # Get best hit per query (highest bitscore)
        best_hits = filtered.loc[filtered.groupby('qseqid')['bitscore'].idxmax()]
        return best_hits
    
    best1 = get_best_hits(blast1)
    best2 = get_best_hits(blast2)
    
    if best1.empty or best2.empty:
        return pd.DataFrame()
    
    # Find reciprocal best hits
    rbh = pd.merge(
        best1[['qseqid', 'sseqid', 'pident', 'evalue', 'bitscore']],
        best2[['qseqid', 'sseqid', 'pident', 'evalue', 'bitscore']],
        left_on=['qseqid', 'sseqid'],
        right_on=['sseqid', 'qseqid'],
        suffixes=('_1', '_2')
    )
    
    if rbh.empty:
        return rbh
    
    # Clean up column names
    rbh = rbh.rename(columns={
        'qseqid_1': 'gene1',
        'sseqid_1': 'gene2',
        'pident_1': 'identity1',
        'pident_2': 'identity2',
        'evalue_1': 'evalue1',
        'evalue_2': 'evalue2',
        'bitscore_1': 'bitscore1',
        'bitscore_2': 'bitscore2'
    })
    
    rbh['avg_identity'] = (rbh['identity1'] + rbh['identity2']) / 2
    rbh['min_evalue'] = rbh[['evalue1', 'evalue2']].min(axis=1)
    rbh['max_bitscore'] = rbh[['bitscore1', 'bitscore2']].max(axis=1)
    
    return rbh[['gene1', 'gene2', 'avg_identity', 'min_evalue', 'max_bitscore']]

def main():
    """Main analysis function."""
    
    output_dir = os.path.dirname(os.path.abspath(__file__))
    
    print("=== Full Coral Orthology Analysis ===")
    print("Loading BLAST results...")
    
    # Load full BLAST results
    apul_vs_peve = load_blast_results(f"{output_dir}/Apul_vs_Peve_full.blastp")
    peve_vs_apul = load_blast_results(f"{output_dir}/Peve_vs_Apul_full.blastp")
    apul_vs_ptua = load_blast_results(f"{output_dir}/Apul_vs_Ptua_full.blastp")
    ptua_vs_apul = load_blast_results(f"{output_dir}/Ptua_vs_Apul_full.blastp")
    peve_vs_ptua = load_blast_results(f"{output_dir}/Peve_vs_Ptua_full.blastp")
    ptua_vs_peve = load_blast_results(f"{output_dir}/Ptua_vs_Peve_full.blastp")
    
    print(f"BLAST hits loaded:")
    print(f"  Apul vs Peve: {len(apul_vs_peve):,}")
    print(f"  Peve vs Apul: {len(peve_vs_apul):,}")
    print(f"  Apul vs Ptua: {len(apul_vs_ptua):,}")
    print(f"  Ptua vs Apul: {len(ptua_vs_apul):,}")
    print(f"  Peve vs Ptua: {len(peve_vs_ptua):,}")
    print(f"  Ptua vs Peve: {len(ptua_vs_peve):,}")
    print()
    
    # Find reciprocal best hits
    print("Identifying reciprocal best hits...")
    apul_peve_rbh = find_reciprocal_best_hits(apul_vs_peve, peve_vs_apul)
    apul_ptua_rbh = find_reciprocal_best_hits(apul_vs_ptua, ptua_vs_apul)
    peve_ptua_rbh = find_reciprocal_best_hits(peve_vs_ptua, ptua_vs_peve)
    
    print(f"Reciprocal best hits:")
    print(f"  Apul-Peve: {len(apul_peve_rbh):,}")
    print(f"  Apul-Ptua: {len(apul_ptua_rbh):,}")
    print(f"  Peve-Ptua: {len(peve_ptua_rbh):,}")
    print()
    
    # Save results
    apul_peve_rbh.to_csv(f"{output_dir}/apul_peve_rbh_full.csv", index=False)
    apul_ptua_rbh.to_csv(f"{output_dir}/apul_ptua_rbh_full.csv", index=False)
    peve_ptua_rbh.to_csv(f"{output_dir}/peve_ptua_rbh_full.csv", index=False)
    
    print("Results saved!")
    print("Run ortholog group analysis to create comprehensive ortholog groups.")

if __name__ == "__main__":
    main()
EOF

# Run orthology analysis
echo "Running orthology analysis..."
cd "$OUTPUT_DIR"
python3 analyze_full_orthologs.py

echo ""
echo "=== Analysis Complete ==="
echo "Completed at: $(date)"
echo ""
echo "Output files:"
echo "  - BLAST results: *_full.blastp"
echo "  - Reciprocal best hits: *_rbh_full.csv"
echo "  - Analysis script: analyze_full_orthologs.py"
echo ""
echo "Next steps:"
echo "  1. Create ortholog groups from RBH results"
echo "  2. Functional enrichment analysis"
echo "  3. Expression correlation analysis"
echo "  4. Phylogenetic analysis"