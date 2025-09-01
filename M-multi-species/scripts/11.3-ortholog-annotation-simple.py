#!/usr/bin/env python3
"""
Simplified Ortholog Group Functional Annotation
This script provides basic annotation functionality without requiring external tools.
"""

import os
import sys
import pandas as pd
import numpy as np
from Bio import SeqIO
from collections import defaultdict, Counter
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

class SimpleOrthologAnnotator:
    """
    Simplified class for annotating ortholog groups with basic functional information
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
    
    def create_basic_annotations(self):
        """
        Create basic annotations based on sequence properties and ortholog group characteristics
        """
        print("Creating basic annotations...")
        
        # Create basic annotation table
        basic_annotations = []
        
        for idx, row in self.ortholog_groups.iterrows():
            group_id = row['group_id']
            
            # Basic annotation based on ortholog type
            annotation = {
                'group_id': group_id,
                'type': row['type'],
                'avg_identity': row['avg_identity'],
                'apul': row['apul'],
                'peve': row['peve'],
                'ptua': row['ptua'],
                'annotation_source': 'basic_analysis',
                'confidence_level': 'medium' if row['type'] == 'three_way' else 'low',
                'species_count': 3 if row['type'] == 'three_way' else 2,
                'conservation_level': 'high' if row['avg_identity'] > 80 else 'medium' if row['avg_identity'] > 60 else 'low'
            }
            
            # Add basic functional predictions based on ortholog type
            if row['type'] == 'three_way':
                annotation['functional_category'] = 'core_conserved'
                annotation['evolutionary_significance'] = 'high'
            else:
                annotation['functional_category'] = 'lineage_specific'
                annotation['evolutionary_significance'] = 'medium'
            
            # Add sequence length information if available
            protein_id = f"{group_id}_Apul"
            annotation['sequence_length'] = None  # Will be filled if sequence data is available
            
            basic_annotations.append(annotation)
        
        # Convert to DataFrame
        annotations_df = pd.DataFrame(basic_annotations)
        
        # Save basic annotations
        output_file = self.output_dir / "basic_annotations.csv"
        annotations_df.to_csv(output_file, index=False)
        print(f"Saved basic annotations to {output_file}")
        
        return annotations_df
    
    def analyze_ortholog_distribution(self, annotations_df):
        """
        Analyze distribution of ortholog groups
        
        Args:
            annotations_df (pd.DataFrame): Basic annotations
        """
        print("Analyzing ortholog distribution...")
        
        # Ortholog type distribution
        type_counts = annotations_df['type'].value_counts()
        type_df = pd.DataFrame({
            'type': type_counts.index,
            'count': type_counts.values,
            'percentage': (type_counts.values / len(annotations_df) * 100).round(1)
        })
        
        # Conservation level distribution
        conservation_counts = annotations_df['conservation_level'].value_counts()
        conservation_df = pd.DataFrame({
            'conservation_level': conservation_counts.index,
            'count': conservation_counts.values,
            'percentage': (conservation_counts.values / len(annotations_df) * 100).round(1)
        })
        
        # Functional category distribution
        func_counts = annotations_df['functional_category'].value_counts()
        func_df = pd.DataFrame({
            'functional_category': func_counts.index,
            'count': func_counts.values,
            'percentage': (func_counts.values / len(annotations_df) * 100).round(1)
        })
        
        # Save analysis results
        type_df.to_csv(self.output_dir / "ortholog_type_distribution.csv", index=False)
        conservation_df.to_csv(self.output_dir / "conservation_level_distribution.csv", index=False)
        func_df.to_csv(self.output_dir / "functional_category_distribution.csv", index=False)
        
        return type_df, conservation_df, func_df
    
    def create_visualizations(self, annotations_df, type_df, conservation_df, func_df):
        """
        Create visualization plots
        
        Args:
            annotations_df (pd.DataFrame): Basic annotations
            type_df (pd.DataFrame): Ortholog type distribution
            conservation_df (pd.DataFrame): Conservation level distribution
            func_df (pd.DataFrame): Functional category distribution
        """
        print("Creating visualizations...")
        
        # Set up plotting style
        plt.style.use('default')
        sns.set_palette("husl")
        
        # 1. Ortholog type distribution
        plt.figure(figsize=(10, 6))
        plt.pie(type_df['count'], labels=type_df['type'], autopct='%1.1f%%', startangle=90)
        plt.title('Distribution of Ortholog Types')
        plt.axis('equal')
        plt.tight_layout()
        plt.savefig(self.output_dir / "ortholog_type_distribution.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        # 2. Conservation level distribution
        plt.figure(figsize=(10, 6))
        colors = ['#ff7f0e', '#2ca02c', '#d62728']
        plt.bar(conservation_df['conservation_level'], conservation_df['count'], color=colors)
        plt.xlabel('Conservation Level')
        plt.ylabel('Count')
        plt.title('Distribution of Conservation Levels')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(self.output_dir / "conservation_level_distribution.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        # 3. Functional category distribution
        plt.figure(figsize=(10, 6))
        plt.bar(func_df['functional_category'], func_df['count'])
        plt.xlabel('Functional Category')
        plt.ylabel('Count')
        plt.title('Distribution of Functional Categories')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(self.output_dir / "functional_category_distribution.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        # 4. Average identity distribution by ortholog type
        plt.figure(figsize=(10, 6))
        for ortholog_type in annotations_df['type'].unique():
            subset = annotations_df[annotations_df['type'] == ortholog_type]
            plt.hist(subset['avg_identity'], alpha=0.7, label=ortholog_type, bins=30)
        
        plt.xlabel('Average Identity (%)')
        plt.ylabel('Count')
        plt.title('Distribution of Average Identity by Ortholog Type')
        plt.legend()
        plt.tight_layout()
        plt.savefig(self.output_dir / "identity_distribution_by_type.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        print("Visualizations saved to output directory")
    
    def create_summary_report(self, annotations_df):
        """
        Create comprehensive summary report
        
        Args:
            annotations_df (pd.DataFrame): Basic annotations
        """
        print("Creating summary report...")
        
        # Calculate statistics
        total_groups = len(annotations_df)
        three_way_groups = sum(annotations_df['type'] == 'three_way')
        two_way_groups = sum(annotations_df['type'] != 'three_way')
        
        high_conservation = sum(annotations_df['conservation_level'] == 'high')
        medium_conservation = sum(annotations_df['conservation_level'] == 'medium')
        low_conservation = sum(annotations_df['conservation_level'] == 'low')
        
        core_conserved = sum(annotations_df['functional_category'] == 'core_conserved')
        lineage_specific = sum(annotations_df['functional_category'] == 'lineage_specific')
        
        # Create summary data
        summary_data = {
            'Metric': [
                'Total Ortholog Groups',
                'Three-way Orthologs',
                'Two-way Orthologs',
                'High Conservation (>80% identity)',
                'Medium Conservation (60-80% identity)',
                'Low Conservation (<60% identity)',
                'Core Conserved Functions',
                'Lineage-specific Functions',
                'Average Identity (Overall)',
                'Average Identity (Three-way)',
                'Average Identity (Two-way)'
            ],
            'Count': [
                total_groups,
                three_way_groups,
                two_way_groups,
                high_conservation,
                medium_conservation,
                low_conservation,
                core_conserved,
                lineage_specific,
                round(annotations_df['avg_identity'].mean(), 1),
                round(annotations_df[annotations_df['type'] == 'three_way']['avg_identity'].mean(), 1),
                round(annotations_df[annotations_df['type'] != 'three_way']['avg_identity'].mean(), 1)
            ],
            'Percentage': [
                100,
                round(three_way_groups / total_groups * 100, 1),
                round(two_way_groups / total_groups * 100, 1),
                round(high_conservation / total_groups * 100, 1),
                round(medium_conservation / total_groups * 100, 1),
                round(low_conservation / total_groups * 100, 1),
                round(core_conserved / total_groups * 100, 1),
                round(lineage_specific / total_groups * 100, 1),
                '-',
                '-',
                '-'
            ]
        }
        
        summary_df = pd.DataFrame(summary_data)
        summary_df.to_csv(self.output_dir / "annotation_summary_report.csv", index=False)
        
        print("Summary report saved to output directory")
        print("\nAnnotation Summary:")
        print(summary_df.to_string(index=False))
        
        return summary_df
    
    def create_annotation_database(self, annotations_df):
        """
        Create a comprehensive annotation database for downstream analysis
        
        Args:
            annotations_df (pd.DataFrame): Basic annotations
        """
        print("Creating annotation database...")
        
        # Create lookup tables
        annotation_db = list()
        
        # 1. Ortholog group to annotation mapping
        annotation_db.append(annotations_df[['group_id', 'type', 'avg_identity', 'annotation_source', 
                                           'confidence_level', 'species_count', 'conservation_level',
                                           'functional_category', 'evolutionary_significance']])
        
        # 2. Gene ID to ortholog group mapping
        gene_mapping = []
        for _, row in annotations_df.iterrows():
            for species in ['apul', 'peve', 'ptua']:
                gene_id = row[species]
                if pd.notna(gene_id):
                    gene_mapping.append({
                        'gene_id': gene_id,
                        'species': species.capitalize(),
                        'group_id': row['group_id'],
                        'type': row['type'],
                        'conservation_level': row['conservation_level']
                    })
        
        gene_mapping_df = pd.DataFrame(gene_mapping)
        annotation_db.append(gene_mapping_df)
        
        # Save annotation database components
        annotations_df.to_csv(self.output_dir / "ortholog_annotations_database.csv", index=False)
        gene_mapping_df.to_csv(self.output_dir / "gene_to_ortholog_mapping.csv", index=False)
        
        print("Created comprehensive annotation database")
        return annotation_db
    
    def run_simple_annotation_pipeline(self, protein_files):
        """
        Run the simplified annotation pipeline
        
        Args:
            protein_files (dict): Dictionary with species:file_path pairs
        """
        print("Starting simplified annotation pipeline...")
        
        # Step 1: Extract representative sequences
        fasta_file = self.extract_representative_sequences(protein_files)
        
        # Step 2: Create basic annotations
        annotations_df = self.create_basic_annotations()
        
        # Step 3: Analyze distribution
        type_df, conservation_df, func_df = self.analyze_ortholog_distribution(annotations_df)
        
        # Step 4: Create visualizations
        self.create_visualizations(annotations_df, type_df, conservation_df, func_df)
        
        # Step 5: Create summary report
        summary_df = self.create_summary_report(annotations_df)
        
        # Step 6: Create annotation database
        annotation_db = self.create_annotation_database(annotations_df)
        
        print("Simplified annotation pipeline completed successfully!")
        return annotations_df


def main():
    """
    Main function to run the simplified annotation pipeline
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
    annotator = SimpleOrthologAnnotator(ortholog_file, output_dir)
    
    # Run simplified pipeline
    annotations_df = annotator.run_simple_annotation_pipeline(protein_files)
    
    print(f"\nSimplified annotation completed! Results saved to: {output_dir}")


if __name__ == "__main__":
    main()
