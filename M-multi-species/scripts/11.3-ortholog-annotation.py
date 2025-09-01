#!/usr/bin/env python3
"""
Ortholog Group Functional Annotation - Python Utilities

This script provides additional functionality for annotating ortholog groups
with functional information from various databases and sources.

Author: Multi-species comparative analysis
Date: 2024
"""

import os
import sys
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import subprocess
import requests
import json
from collections import defaultdict, Counter
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

class OrthologAnnotator:
    """
    Class for annotating ortholog groups with functional information
    """
    
    def __init__(self, ortholog_file, output_dir):
        """
        Initialize the annotator
        
        Args:
            ortholog_file (str): Path to ortholog groups CSV file
            output_dir (str): Output directory for results
        """
        self.ortholog_file = ortholog_file
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Load ortholog groups
        self.ortholog_groups = pd.read_csv(ortholog_file)
        print(f"Loaded {len(self.ortholog_groups)} ortholog groups")
        
        # Initialize annotation results
        self.annotations = {}
        
    def extract_representative_sequences(self, protein_files):
        """
        Extract representative protein sequences for each ortholog group
        
        Args:
            protein_files (dict): Dictionary with species:file_path pairs
        """
        print("Extracting representative sequences...")
        
        # Load protein sequences
        protein_sequences = {}
        for species, file_path in protein_files.items():
            if os.path.exists(file_path):
                protein_sequences[species] = {seq.id: seq for seq in SeqIO.parse(file_path, "fasta")}
                print(f"Loaded {len(protein_sequences[species])} sequences for {species}")
            else:
                print(f"Warning: Protein file not found for {species}: {file_path}")
        
        # Extract representative sequences
        representative_sequences = []
        representative_ids = []
        
        for idx, row in self.ortholog_groups.iterrows():
            group_id = row['group_id']
            sequences = {}
            
            # Extract sequences for each species
            for species in ['apul', 'peve', 'ptua']:
                gene_id = row[species]
                if pd.notna(gene_id) and species in protein_sequences:
                    if gene_id in protein_sequences[species]:
                        sequences[species] = protein_sequences[species][gene_id]
            
            # Choose representative sequence (prefer three-way orthologs)
            if row['type'] == 'three_way' and len(sequences) == 3:
                # Use Apul as representative for three-way orthologs
                if 'apul' in sequences:
                    representative_sequences.append(sequences['apul'])
                    representative_ids.append(f"{group_id}_Apul")
            else:
                # Use available sequence for two-way orthologs
                for species in ['apul', 'peve', 'ptua']:
                    if species in sequences:
                        representative_sequences.append(sequences[species])
                        representative_ids.append(f"{group_id}_{species.capitalize()}")
                        break
        
        # Write representative sequences to FASTA
        output_file = self.output_dir / "representative_sequences.faa"
        SeqIO.write(representative_sequences, output_file, "fasta")
        print(f"Wrote {len(representative_sequences)} representative sequences to {output_file}")
        
        return output_file
    
    def run_interproscan(self, fasta_file):
        """
        Run InterProScan on protein sequences
        
        Args:
            fasta_file (str): Path to FASTA file with protein sequences
        """
        print("Running InterProScan...")
        
        output_file = self.output_dir / "interproscan_results.tsv"
        
        # Check if InterProScan is available
        try:
            cmd = [
                "interproscan.sh",
                "-i", str(fasta_file),
                "-o", str(output_file),
                "-f", "TSV",
                "-dp",
                "-pa",
                "-goterms",
                "-pathways",
                "-appl", "Pfam,ProDom,SMART,SUPERFAMILY,PRINTS,ProSiteProfiles,ProSitePatterns"
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode == 0:
                print(f"InterProScan completed successfully. Results saved to {output_file}")
                return output_file
            else:
                print(f"InterProScan failed: {result.stderr}")
                return None
                
        except FileNotFoundError:
            print("InterProScan not found. Please install InterProScan or use alternative methods.")
            return None
    
    def run_blast_swissprot(self, fasta_file):
        """
        Run BLAST against Swiss-Prot database
        
        Args:
            fasta_file (str): Path to FASTA file with protein sequences
        """
        print("Running BLAST against Swiss-Prot...")
        
        # Check if Swiss-Prot database exists
        swissprot_db = Path("../data/swissprot")
        if not swissprot_db.with_suffix(".pin").exists():
            print("Swiss-Prot database not found. Please run the R script first to download it.")
            return None
        
        output_file = self.output_dir / "swissprot_blast.tsv"
        
        try:
            cmd = [
                "blastp",
                "-query", str(fasta_file),
                "-db", str(swissprot_db),
                "-out", str(output_file),
                "-outfmt", "6 qseqid sseqid pident qcovs evalue bitscore stitle",
                "-evalue", "1e-5",
                "-max_target_seqs", "5"
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode == 0:
                print(f"BLAST completed successfully. Results saved to {output_file}")
                return output_file
            else:
                print(f"BLAST failed: {result.stderr}")
                return None
                
        except FileNotFoundError:
            print("BLAST not found. Please install BLAST+.")
            return None
    
    def parse_interproscan_results(self, interpro_file):
        """
        Parse InterProScan results
        
        Args:
            interpro_file (str): Path to InterProScan results file
        """
        if not os.path.exists(interpro_file):
            print(f"InterProScan results file not found: {interpro_file}")
            return None
        
        print("Parsing InterProScan results...")
        
        # Read InterProScan results
        columns = [
            "protein_id", "md5", "length", "analysis", "signature_accession",
            "signature_description", "start", "stop", "score", "status", "date",
            "interpro_accession", "interpro_description", "go_terms", "pathways"
        ]
        
        interpro_data = pd.read_csv(interpro_file, sep='\t', header=None, names=columns)
        
        # Group by protein ID
        protein_annotations = {}
        for protein_id, group in interpro_data.groupby('protein_id'):
            annotations = {
                'pfam_domains': [],
                'go_terms': [],
                'kegg_pathways': [],
                'interpro_terms': []
            }
            
            for _, row in group.iterrows():
                # Extract Pfam domains
                if row['analysis'] == 'Pfam':
                    annotations['pfam_domains'].append(row['signature_accession'])
                
                # Extract GO terms
                if pd.notna(row['go_terms']):
                    go_terms = row['go_terms'].split('|')
                    annotations['go_terms'].extend(go_terms)
                
                # Extract KEGG pathways
                if pd.notna(row['pathways']):
                    pathways = row['pathways'].split('|')
                    annotations['kegg_pathways'].extend(pathways)
                
                # Extract InterPro terms
                if pd.notna(row['interpro_accession']):
                    annotations['interpro_terms'].append(row['interpro_accession'])
            
            # Remove duplicates
            for key in annotations:
                annotations[key] = list(set(annotations[key]))
            
            protein_annotations[protein_id] = annotations
        
        return protein_annotations
    
    def parse_blast_results(self, blast_file):
        """
        Parse BLAST results
        
        Args:
            blast_file (str): Path to BLAST results file
        """
        if not os.path.exists(blast_file):
            print(f"BLAST results file not found: {blast_file}")
            return None
        
        print("Parsing BLAST results...")
        
        # Read BLAST results
        columns = ['query_id', 'subject_id', 'pident', 'qcovs', 'evalue', 'bitscore', 'subject_title']
        blast_data = pd.read_csv(blast_file, sep='\t', header=None, names=columns)
        
        # Get best hit for each query
        best_hits = blast_data.loc[blast_data.groupby('query_id')['evalue'].idxmin()]
        
        return best_hits.set_index('query_id')
    
    def integrate_annotations(self, interpro_annotations, blast_results):
        """
        Integrate annotations from multiple sources
        
        Args:
            interpro_annotations (dict): InterProScan annotations
            blast_results (pd.DataFrame): BLAST results
        """
        print("Integrating annotations...")
        
        # Create integrated annotation table
        integrated_annotations = []
        
        for idx, row in self.ortholog_groups.iterrows():
            group_id = row['group_id']
            protein_id = f"{group_id}_Apul"  # Assuming Apul as representative
            
            annotation = {
                'group_id': group_id,
                'type': row['type'],
                'avg_identity': row['avg_identity'],
                'apul': row['apul'],
                'peve': row['peve'],
                'ptua': row['ptua'],
                'swissprot_hit': None,
                'swissprot_description': None,
                'swissprot_evalue': None,
                'pfam_domains': None,
                'go_terms': None,
                'kegg_pathways': None,
                'interpro_terms': None
            }
            
            # Add Swiss-Prot annotations
            if blast_results is not None and protein_id in blast_results.index:
                blast_hit = blast_results.loc[protein_id]
                annotation['swissprot_hit'] = blast_hit['subject_id']
                annotation['swissprot_description'] = blast_hit['subject_title']
                annotation['swissprot_evalue'] = blast_hit['evalue']
            
            # Add InterPro annotations
            if interpro_annotations is not None and protein_id in interpro_annotations:
                interpro_data = interpro_annotations[protein_id]
                annotation['pfam_domains'] = ';'.join(interpro_data['pfam_domains'])
                annotation['go_terms'] = ';'.join(interpro_data['go_terms'])
                annotation['kegg_pathways'] = ';'.join(interpro_data['kegg_pathways'])
                annotation['interpro_terms'] = ';'.join(interpro_data['interpro_terms'])
            
            integrated_annotations.append(annotation)
        
        # Convert to DataFrame
        annotations_df = pd.DataFrame(integrated_annotations)
        
        # Save integrated annotations
        output_file = self.output_dir / "integrated_annotations.csv"
        annotations_df.to_csv(output_file, index=False)
        print(f"Saved integrated annotations to {output_file}")
        
        return annotations_df
    
    def analyze_functional_distribution(self, annotations_df):
        """
        Analyze functional distribution of annotations
        
        Args:
            annotations_df (pd.DataFrame): Integrated annotations
        """
        print("Analyzing functional distribution...")
        
        # Pfam domain analysis
        pfam_domains = []
        for domains in annotations_df['pfam_domains'].dropna():
            if domains:
                pfam_domains.extend(domains.split(';'))
        
        pfam_counts = Counter(pfam_domains)
        pfam_df = pd.DataFrame(list(pfam_counts.items()), columns=['domain', 'count'])
        pfam_df = pfam_df.sort_values('count', ascending=False).head(50)
        
        # GO term analysis
        go_terms = []
        for terms in annotations_df['go_terms'].dropna():
            if terms:
                go_terms.extend(terms.split(';'))
        
        go_counts = Counter(go_terms)
        go_df = pd.DataFrame(list(go_counts.items()), columns=['term', 'count'])
        go_df = go_df.sort_values('count', ascending=False).head(50)
        
        # Save analysis results
        pfam_df.to_csv(self.output_dir / "pfam_domain_analysis.csv", index=False)
        go_df.to_csv(self.output_dir / "go_term_analysis.csv", index=False)
        
        return pfam_df, go_df
    
    def create_visualizations(self, annotations_df, pfam_df, go_df):
        """
        Create visualization plots
        
        Args:
            annotations_df (pd.DataFrame): Integrated annotations
            pfam_df (pd.DataFrame): Pfam domain analysis
            go_df (pd.DataFrame): GO term analysis
        """
        print("Creating visualizations...")
        
        # Set up plotting style
        plt.style.use('default')
        sns.set_palette("husl")
        
        # 1. Pfam domain distribution
        plt.figure(figsize=(12, 8))
        top_pfam = pfam_df.head(20)
        plt.barh(range(len(top_pfam)), top_pfam['count'])
        plt.yticks(range(len(top_pfam)), top_pfam['domain'])
        plt.xlabel('Count')
        plt.title('Top 20 Pfam Domains in Ortholog Groups')
        plt.tight_layout()
        plt.savefig(self.output_dir / "pfam_domain_distribution.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        # 2. GO term distribution
        plt.figure(figsize=(12, 8))
        top_go = go_df.head(20)
        plt.barh(range(len(top_go)), top_go['count'])
        plt.yticks(range(len(top_go)), top_go['term'])
        plt.xlabel('Count')
        plt.title('Top 20 GO Terms in Ortholog Groups')
        plt.tight_layout()
        plt.savefig(self.output_dir / "go_term_distribution.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        # 3. Annotation coverage by ortholog type
        coverage_data = []
        for ortholog_type in annotations_df['type'].unique():
            subset = annotations_df[annotations_df['type'] == ortholog_type]
            total = len(subset)
            
            coverage_data.append({
                'type': ortholog_type,
                'total': total,
                'swissprot': sum(subset['swissprot_hit'].notna()),
                'pfam': sum(subset['pfam_domains'].notna() & (subset['pfam_domains'] != '')),
                'go': sum(subset['go_terms'].notna() & (subset['go_terms'] != '')),
                'kegg': sum(subset['kegg_pathways'].notna() & (subset['kegg_pathways'] != ''))
            })
        
        coverage_df = pd.DataFrame(coverage_data)
        
        # Plot annotation coverage
        fig, ax = plt.subplots(figsize=(10, 6))
        x = np.arange(len(coverage_df))
        width = 0.2
        
        ax.bar(x - width*1.5, coverage_df['swissprot'] / coverage_df['total'] * 100, 
               width, label='Swiss-Prot', alpha=0.8)
        ax.bar(x - width*0.5, coverage_df['pfam'] / coverage_df['total'] * 100, 
               width, label='Pfam', alpha=0.8)
        ax.bar(x + width*0.5, coverage_df['go'] / coverage_df['total'] * 100, 
               width, label='GO Terms', alpha=0.8)
        ax.bar(x + width*1.5, coverage_df['kegg'] / coverage_df['total'] * 100, 
               width, label='KEGG', alpha=0.8)
        
        ax.set_xlabel('Ortholog Type')
        ax.set_ylabel('Percentage (%)')
        ax.set_title('Annotation Coverage by Ortholog Type')
        ax.set_xticks(x)
        ax.set_xticklabels(coverage_df['type'])
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(self.output_dir / "annotation_coverage.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        print("Visualizations saved to output directory")
    
    def create_summary_report(self, annotations_df):
        """
        Create comprehensive summary report
        
        Args:
            annotations_df (pd.DataFrame): Integrated annotations
        """
        print("Creating summary report...")
        
        # Calculate statistics
        total_groups = len(annotations_df)
        three_way_groups = sum(annotations_df['type'] == 'three_way')
        two_way_groups = sum(annotations_df['type'] != 'three_way')
        
        annotated_groups = sum(
            annotations_df['swissprot_hit'].notna() |
            (annotations_df['pfam_domains'].notna() & (annotations_df['pfam_domains'] != '')) |
            (annotations_df['go_terms'].notna() & (annotations_df['go_terms'] != ''))
        )
        
        # Create summary data
        summary_data = {
            'Metric': [
                'Total Ortholog Groups',
                'Three-way Orthologs',
                'Two-way Orthologs',
                'Annotated Groups',
                'Swiss-Prot Annotations',
                'Pfam Domain Annotations',
                'GO Term Annotations',
                'KEGG Pathway Annotations'
            ],
            'Count': [
                total_groups,
                three_way_groups,
                two_way_groups,
                annotated_groups,
                sum(annotations_df['swissprot_hit'].notna()),
                sum(annotations_df['pfam_domains'].notna() & (annotations_df['pfam_domains'] != '')),
                sum(annotations_df['go_terms'].notna() & (annotations_df['go_terms'] != '')),
                sum(annotations_df['kegg_pathways'].notna() & (annotations_df['kegg_pathways'] != ''))
            ],
            'Percentage': [
                100,
                round(three_way_groups / total_groups * 100, 1),
                round(two_way_groups / total_groups * 100, 1),
                round(annotated_groups / total_groups * 100, 1),
                round(sum(annotations_df['swissprot_hit'].notna()) / total_groups * 100, 1),
                round(sum(annotations_df['pfam_domains'].notna() & (annotations_df['pfam_domains'] != '')) / total_groups * 100, 1),
                round(sum(annotations_df['go_terms'].notna() & (annotations_df['go_terms'] != '')) / total_groups * 100, 1),
                round(sum(annotations_df['kegg_pathways'].notna() & (annotations_df['kegg_pathways'] != '')) / total_groups * 100, 1)
            ]
        }
        
        summary_df = pd.DataFrame(summary_data)
        summary_df.to_csv(self.output_dir / "annotation_summary_report.csv", index=False)
        
        print("Summary report saved to output directory")
        print("\nAnnotation Summary:")
        print(summary_df.to_string(index=False))
        
        return summary_df
    
    def run_full_annotation_pipeline(self, protein_files):
        """
        Run the complete annotation pipeline
        
        Args:
            protein_files (dict): Dictionary with species:file_path pairs
        """
        print("Starting full annotation pipeline...")
        
        # Step 1: Extract representative sequences
        fasta_file = self.extract_representative_sequences(protein_files)
        
        # Step 2: Run InterProScan
        interpro_file = self.run_interproscan(fasta_file)
        
        # Step 3: Run BLAST against Swiss-Prot
        blast_file = self.run_blast_swissprot(fasta_file)
        
        # Step 4: Parse results
        interpro_annotations = self.parse_interproscan_results(interpro_file)
        blast_results = self.parse_blast_results(blast_file)
        
        # Step 5: Integrate annotations
        annotations_df = self.integrate_annotations(interpro_annotations, blast_results)
        
        # Step 6: Analyze functional distribution
        pfam_df, go_df = self.analyze_functional_distribution(annotations_df)
        
        # Step 7: Create visualizations
        self.create_visualizations(annotations_df, pfam_df, go_df)
        
        # Step 8: Create summary report
        summary_df = self.create_summary_report(annotations_df)
        
        print("Annotation pipeline completed successfully!")
        return annotations_df


def main():
    """
    Main function to run the annotation pipeline
    """
    # Define file paths
    ortholog_file = "../output/11-orthology-analysis/ortholog_groups.csv"
    output_dir = "../output/11.3-ortholog-annotation"
    
    protein_files = {
        'apul': "../../D-Apul/data/Apulchra-genome.pep.faa",
        'peve': "../../E-Peve/data/Porites_evermanni_v1.annot.pep.fa",
        'ptua': "../../F-Ptua/data/Pocillopora_meandrina_HIv1.genes.pep.faa"
    }
    
    # Check if ortholog file exists
    if not os.path.exists(ortholog_file):
        print(f"Error: Ortholog file not found: {ortholog_file}")
        print("Please run the orthology analysis first (script 11-orthology-analysis.Rmd)")
        sys.exit(1)
    
    # Initialize annotator
    annotator = OrthologAnnotator(ortholog_file, output_dir)
    
    # Run full pipeline
    annotations_df = annotator.run_full_annotation_pipeline(protein_files)
    
    print(f"\nAnnotation completed! Results saved to: {output_dir}")


if __name__ == "__main__":
    main()
