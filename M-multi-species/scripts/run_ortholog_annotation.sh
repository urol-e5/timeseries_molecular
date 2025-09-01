#!/bin/bash

# Ortholog Annotation Pipeline Runner
# This script runs the complete ortholog annotation pipeline

# Set error handling
set -euo pipefail

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
OUTPUT_DIR="$PROJECT_ROOT/output/11.3-ortholog-annotation"
LOG_FILE="$OUTPUT_DIR/annotation_pipeline.log"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Logging function
log() {
    echo -e "${BLUE}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} $1" | tee -a "$LOG_FILE"
}

log_success() {
    echo -e "${GREEN}[$(date +'%Y-%m-%d %H:%M:%S')] SUCCESS:${NC} $1" | tee -a "$LOG_FILE"
}

log_warning() {
    echo -e "${YELLOW}[$(date +'%Y-%m-%d %H:%M:%S')] WARNING:${NC} $1" | tee -a "$LOG_FILE"
}

log_error() {
    echo -e "${RED}[$(date +'%Y-%m-%d %H:%M:%S')] ERROR:${NC} $1" | tee -a "$LOG_FILE"
}

# Function to check dependencies
check_dependencies() {
    log "Checking dependencies..."
    
    # Check for R
    if ! command -v R &> /dev/null; then
        log_error "R is not installed or not in PATH"
        exit 1
    fi
    
    # Check for Python
    if ! command -v python3 &> /dev/null; then
        log_error "Python3 is not installed or not in PATH"
        exit 1
    fi
    
    # Check for BLAST
    if ! command -v blastp &> /dev/null; then
        log_warning "BLAST+ is not installed. Some annotation steps may be skipped."
    fi
    
    # Check for InterProScan
    if ! command -v interproscan.sh &> /dev/null; then
        log_warning "InterProScan is not installed. Some annotation steps may be skipped."
    fi
    
    # Check for required R packages
    log "Checking R packages..."
    R --slave -e "
    required_packages <- c('tidyverse', 'Biostrings', 'ggplot2', 'VennDiagram', 'pheatmap')
    missing_packages <- required_packages[!required_packages %in% installed.packages()[,'Package']]
    if (length(missing_packages) > 0) {
        cat('Missing R packages:', paste(missing_packages, collapse=', '), '\n')
        cat('Please install them using: install.packages(c(', paste0('\"', missing_packages, '\"', collapse=', '), '))\n')
        quit(status = 1)
    } else {
        cat('All required R packages are installed\n')
    }
    "
    
    # Check for required Python packages
    log "Checking Python packages..."
    python3 -c "
import sys
required_packages = ['pandas', 'numpy', 'Bio', 'matplotlib', 'seaborn']
missing_packages = []
for package in required_packages:
    try:
        __import__(package)
    except ImportError:
        missing_packages.append(package)
if missing_packages:
    print(f'Missing Python packages: {missing_packages}')
    print('Please install them using: pip install ' + ' '.join(missing_packages))
    sys.exit(1)
else:
    print('All required Python packages are installed')
"
    
    log_success "Dependency check completed"
}

# Function to check input files
check_input_files() {
    log "Checking input files..."
    
    # Check ortholog groups file
    ORTHOLOG_FILE="$PROJECT_ROOT/output/11-orthology-analysis/ortholog_groups.csv"
    if [[ ! -f "$ORTHOLOG_FILE" ]]; then
        log_error "Ortholog groups file not found: $ORTHOLOG_FILE"
        log_error "Please run the orthology analysis first (script 11-orthology-analysis.Rmd)"
        exit 1
    fi
    
    # Check protein files
    PROTEIN_FILES=(
        "$PROJECT_ROOT/../D-Apul/data/Apulchra-genome.pep.faa"
        "$PROJECT_ROOT/../E-Peve/data/Porites_evermanni_v1.annot.pep.fa"
        "$PROJECT_ROOT/../F-Ptua/data/Pocillopora_meandrina_HIv1.genes.pep.faa"
    )
    
    for file in "${PROTEIN_FILES[@]}"; do
        if [[ ! -f "$file" ]]; then
            log_warning "Protein file not found: $file"
        else
            log "Found protein file: $file"
        fi
    done
    
    log_success "Input file check completed"
}

# Function to create output directory
create_output_dir() {
    log "Creating output directory..."
    mkdir -p "$OUTPUT_DIR"
    log_success "Output directory created: $OUTPUT_DIR"
}

# Function to run R annotation script
run_r_annotation() {
    log "Running R annotation script..."
    
    cd "$SCRIPT_DIR"
    
    # Run R Markdown script
    if Rscript -e "
    if (!require(rmarkdown)) {
        install.packages('rmarkdown', repos='https://cran.rstudio.com/')
    }
    rmarkdown::render('11.3-ortholog-annotation.Rmd', 
                     output_format = 'html_document',
                     output_file = '$OUTPUT_DIR/ortholog_annotation_report.html',
                     params = list(output_dir = '$OUTPUT_DIR'))
    "; then
        log_success "R annotation script completed"
    else
        log_error "R annotation script failed"
        return 1
    fi
}

# Function to run Python annotation script
run_python_annotation() {
    log "Running Python annotation script..."
    
    cd "$SCRIPT_DIR"
    
    if python3 "11.3-ortholog-annotation.py"; then
        log_success "Python annotation script completed"
    else
        log_error "Python annotation script failed"
        return 1
    fi
}

# Function to create summary report
create_summary_report() {
    log "Creating summary report..."
    
    cat > "$OUTPUT_DIR/README.md" << 'EOF'
# Ortholog Annotation Results

This directory contains the results of the ortholog group functional annotation analysis.

## Files Description

### Core Results
- `integrated_annotations.csv` - Complete annotation table with all functional information
- `annotation_summary_report.csv` - Summary statistics of annotation coverage
- `functional_enrichment_analysis.csv` - Enriched functional terms in three-way orthologs

### Database Files
- `ortholog_annotations_database.csv` - Annotation database for downstream analysis
- `gene_to_ortholog_mapping.csv` - Mapping between gene IDs and ortholog groups
- `pfam_terms_database.csv` - Database of Pfam domains found in orthologs
- `go_terms_database.csv` - Database of GO terms found in orthologs
- `kegg_pathways_database.csv` - Database of KEGG pathways found in orthologs

### Analysis Results
- `pfam_domain_analysis.csv` - Analysis of Pfam domain distribution
- `go_term_analysis.csv` - Analysis of GO term distribution

### Visualizations
- `pfam_domain_distribution.png` - Top Pfam domains in ortholog groups
- `go_term_distribution.png` - Top GO terms in ortholog groups
- `annotation_coverage.png` - Annotation coverage by ortholog type
- `functional_enrichment.png` - Enriched functional terms

### Raw Data
- `representative_sequences.faa` - Representative protein sequences used for annotation
- `interproscan_results.tsv` - Raw InterProScan results (if available)
- `swissprot_blast.tsv` - Raw BLAST results against Swiss-Prot (if available)

## Usage

The annotated ortholog groups can be used for:

1. **Comparative functional genomics** - Understanding functional conservation across coral species
2. **Pathway analysis** - Identifying conserved metabolic and signaling pathways
3. **Gene set enrichment analysis** - Finding overrepresented functional terms in expression studies
4. **Evolutionary analysis** - Studying functional evolution of gene families
5. **Cross-species expression analysis** - Comparing expression of functionally annotated orthologs

## Citation

If you use these annotations in your research, please cite:

- InterProScan: Mitchell et al. (2019) Nucleic Acids Research
- Swiss-Prot: Bairoch & Apweiler (2000) Nucleic Acids Research
- Pfam: El-Gebali et al. (2019) Nucleic Acids Research
- Gene Ontology: Ashburner et al. (2000) Nature Genetics
- KEGG: Kanehisa & Goto (2000) Nucleic Acids Research

EOF

    log_success "Summary report created"
}

# Function to validate results
validate_results() {
    log "Validating results..."
    
    # Check if key output files exist
    REQUIRED_FILES=(
        "integrated_annotations.csv"
        "annotation_summary_report.csv"
        "ortholog_annotations_database.csv"
        "gene_to_ortholog_mapping.csv"
    )
    
    for file in "${REQUIRED_FILES[@]}"; do
        if [[ -f "$OUTPUT_DIR/$file" ]]; then
            log "Found required file: $file"
        else
            log_warning "Missing required file: $file"
        fi
    done
    
    # Check annotation coverage
    if [[ -f "$OUTPUT_DIR/annotation_summary_report.csv" ]]; then
        log "Annotation coverage summary:"
        cat "$OUTPUT_DIR/annotation_summary_report.csv" | tail -n +2 | while IFS=',' read -r metric count percentage; do
            echo "  $metric: $count ($percentage%)"
        done
    fi
    
    log_success "Results validation completed"
}

# Main execution function
main() {
    log "Starting ortholog annotation pipeline..."
    
    # Create output directory and log file
    create_output_dir
    touch "$LOG_FILE"
    
    # Run pipeline steps
    check_dependencies
    check_input_files
    run_r_annotation
    run_python_annotation
    create_summary_report
    validate_results
    
    log_success "Ortholog annotation pipeline completed successfully!"
    log "Results are available in: $OUTPUT_DIR"
    log "Log file: $LOG_FILE"
}

# Run main function
main "$@"
