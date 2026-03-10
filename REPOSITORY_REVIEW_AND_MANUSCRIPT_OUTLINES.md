# Repository Review and Manuscript Outlines
## urol-e5/timeseries_molecular

**Date**: October 2024  
**Review Conducted by**: GitHub Copilot Agent  
**Purpose**: Strategic planning for scientific publications

---

## Executive Summary

This repository contains a comprehensive multi-omic time series analysis of three coral species (*Acropora pulchra*, *Porites evermanni*, and *Pocillopora tuahiniensis*) examining molecular responses across developmental/stress timepoints. The repository demonstrates exceptional scientific depth with integrated analyses spanning transcriptomics, epigenomics, small RNA biology, metabolomics, lipidomics, and physiological measurements. The work represents significant computational investment (~26,000 lines of R code) and employs cutting-edge bioinformatics approaches including WGCNA, MOFA2, machine learning, and multi-omic network integration.

**Key Strengths**:
- Comprehensive multi-omic integration (6 data types)
- Three species with distinct life history strategies
- Temporal design (4 timepoints) enabling dynamic analysis
- Advanced analytical methods (ML, network analysis, orthology)
- Well-documented code with clear organization
- Multiple biological and technical replicates

**Publication Potential**: High - The repository contains sufficient data and analyses for multiple high-impact publications ranging from technical methods papers to integrative biological insights.

---

## 1. THOROUGH REPOSITORY REVIEW

### 1.1 Repository Structure and Organization

The repository is exceptionally well-organized following a hierarchical structure:

```
timeseries_molecular/
├── D-Apul/           # Acropora pulchra (45 analysis scripts)
├── E-Peve/           # Porites evermanni (14 analysis scripts)
├── F-Ptua/           # Pocillopora tuahiniensis (13 analysis scripts)
└── M-multi-species/  # Comparative analyses (13 scripts)
```

**Naming Convention**: Scripts follow standardized format `XX.YY-<species>-<analysis>.Rmd` where XX indicates analysis category (e.g., 00=QC, 01=trimming, 02=alignment, 03=expression, etc.)

**Code Base**: 
- Total: ~26,000 lines of R code in D-Apul alone
- ~1,000+ figures/outputs generated
- Consistent use of R Markdown for reproducibility

### 1.2 Experimental Design

**Species Selection** (Strategically diverse):
1. *Acropora pulchra* - Fast-growing, branching morphology
2. *Porites evermanni* - Slow-growing, massive morphology  
3. *Pocillopora tuahiniensis* - Intermediate growth, branching morphology

**Temporal Design**:
- 4 timepoints (TP1-TP4) across experimental period
- Multiple colonies per species for biological replication
- Paired sampling allowing within-colony temporal comparisons

**Data Types Collected**:
1. **Transcriptomics**: RNA-seq (gene and transcript level)
2. **Epigenomics**: WGBS (whole genome bisulfite sequencing)
3. **Small RNA**: sRNA-seq/miRNA discovery and quantification
4. **Non-coding RNA**: lncRNA discovery and expression
5. **Metabolomics**: Comprehensive metabolite profiling
6. **Lipidomics**: Lipid class and species quantification
7. **Physiological**: Growth, photosynthesis, respiration metrics

### 1.3 Analytical Approaches

#### 1.3.1 Data Quality and Processing (Scripts 00-02)
- **Quality Control**: FastQC, MultiQC for all data types
- **Trimming**: fastp with quality filtering
- **Alignment**: 
  - RNA-seq: HiSat2 with GTF-guided alignment
  - WGBS: Bismark for bisulfite-aware alignment
- **Quantification**: 
  - Gene/transcript counts via StringTie/featureCounts
  - CpG methylation levels
  - miRNA abundance matrices

**Strength**: Rigorous QC at each processing step with comprehensive documentation

#### 1.3.2 Differential Expression (Scripts 03)
- DESeq2 for time-series gene expression
- Analysis of genes, transcripts, lncRNAs, and miRNAs
- Multiple testing correction
- Temporal trajectory analysis

**Coverage**: All three species analyzed, though D-Apul has most comprehensive coverage

#### 1.3.3 Small RNA and Regulatory Networks (Scripts 04-07)
- **miRNA Discovery**: ShortStack for de novo miRNA identification
- **Target Prediction**: RNAhybrid and miRanda for miRNA-mRNA interactions
- **3' UTR Annotation**: Custom UTR extraction for target prediction
- **Co-expression Analysis**: Correlation-based miRNA-mRNA networks

**Innovation**: Integration of miRNA discovery with mRNA expression and methylation

#### 1.3.4 Co-expression Network Analysis (Scripts 11-18)
- **WGCNA**: Multiple implementations
  - mRNA-only networks (Script 11)
  - miRNA-mRNA integrated networks (Script 12)
  - miRNA-mRNA-lncRNA tri-partite networks (Script 16)
  - Multi-omic correlation networks (Script 18)
- **Module Annotation**: GO enrichment and functional characterization (Script 21)
- **Trait Association**: Correlation with physiological measurements

**Strength**: Progressive integration of increasing data complexity

#### 1.3.5 Machine Learning and Predictive Modeling (Scripts 20-22)
- **Algorithms**: Random Forest, elastic net, support vector machines
- **Applications**:
  - Prediction of physiological state from mRNA+WGBS (Script 20)
  - miRNA-mRNA based phenotype prediction (Script 22)
  - Multi-omic integration for trait prediction (Scripts 22.2-22.6)
- **Validation**: Cross-validation, feature importance analysis
- **Time-specific Models**: By-timepoint analysis (Script 22.3)

**Innovation**: Progressive refinement with updated WGBS data (Script 22.6)

#### 1.3.6 Multi-Omic Factor Analysis (Script 22.20)
- **MOFA2 Implementation**: State-of-the-art multi-omic integration
- **Data Integration**: Genes, transcripts, lncRNAs, miRNAs, methylation
- **Factor Discovery**: Latent factor identification
- **Variance Decomposition**: Quantification of multi-omic contributions
- **Sample/Feature Characterization**: Factor loadings and weights

**Significance**: MOFA2 is a cutting-edge approach; this represents sophisticated analysis

#### 1.3.7 Epigenetic Analysis (Scripts 14, 19, 26, 28-30)
- **CpG Island Annotation**: Genomic distribution analysis (Script 28)
- **Methylation Patterns**: Genome-wide profiling
- **SNP Calling**: From bisulfite sequencing (Script 19)
- **Gene Body Methylation**: Association with expression
- **EpiMachine**: Expression prediction from methylation (Scripts 29, 36)
- **Bismark Summary**: Comprehensive methylation statistics (Script 30)

**Depth**: Multiple complementary approaches to epigenetic regulation

#### 1.3.8 Functional Annotation (Scripts 21, 23, 25, 27, 35)
- **GO Enrichment**: topGO for biological process identification
- **Protein Annotation**: BLAST against reference databases (Script 35)
- **Energetic State**: Metabolic pathway enrichment (Script 23)
- **miRNA-mRNA Interactions**: Functional context of regulation (Script 27)
- **Gene Set Analysis**: Pathway-level interpretation (Script 25)

**Completeness**: Multiple annotation databases integrated

#### 1.3.9 Multi-Species Comparative Analysis (M-multi-species)
- **Orthology Identification**: All-vs-all BLAST with reciprocal best hits (Script 11)
- **Ortholog Annotation**: Functional comparison across species (Scripts 11.3, 13)
- **Metabolomics Integration**: Cross-species metabolite profiling (Scripts 04, 07)
- **Lipidomics Analysis**: Lipid class comparison (Scripts 03, 05, 06)
- **Metabolomics-Lipidomics-CpG Integration**: Multi-platform correlation (Script 08)
- **WGCNA on Metabolomics**: Network analysis of metabolite data (Script 09)
- **Physiological PC Analysis**: Principal component integration (Script 02)

**Uniqueness**: Comparative framework enables evolutionary/ecological insights

### 1.4 Scientific Themes and Biological Insights

#### Theme 1: Temporal Molecular Dynamics
- Gene expression changes across timepoints
- miRNA regulation patterns over time
- Epigenetic reprogramming during development/stress
- Metabolic shifts in response to conditions

#### Theme 2: Multi-Omic Integration and Regulatory Networks
- miRNA-mRNA regulatory circuits
- DNA methylation-gene expression relationships
- Metabolite-lipid-transcript correlations
- Multi-layer network architecture

#### Theme 3: Species-Specific Molecular Strategies
- Orthologous gene identification (shared vs. species-specific)
- Divergent regulatory mechanisms
- Life history strategy correlates (fast vs. slow growth)
- Morphology-function relationships

#### Theme 4: Predictive Molecular Signatures
- Machine learning models for phenotype prediction
- Feature selection identifies key biomarkers
- Integration improves prediction accuracy
- Temporal specificity of predictive features

#### Theme 5: Epigenetic Regulation in Corals
- CpG methylation patterns
- Gene body methylation roles
- miRNA-mediated regulation
- Environmentally-induced epigenetic changes

### 1.5 Key Datasets and Metadata

**Primary Datasets**:
1. RNA-seq: Gene and transcript count matrices (3 species)
2. miRNA-seq: miRNA abundance matrices with discovery annotations
3. WGBS: CpG methylation levels genome-wide
4. lncRNA: Long non-coding RNA expression
5. Metabolomics: ~200+ metabolite relative quantifications
6. Lipidomics: Lipid class and species measurements
7. Physiological: Master timeseries with growth/photo/respiration

**Metadata Files**:
- `M-multi-species/data/rna_metadata.csv`: Sample information
- Sample submission forms with colony IDs and experimental design
- MEDFORD format metadata for reproducibility

**External Links**:
- Physiological data: E5 Timeseries Repository
- Large output files: gannet.fish.washington.edu hosting

### 1.6 Contributors and Expertise

**Sam White** - Data Engineering & Bioinformatics Infrastructure
- Established computational pipelines
- RNA-seq and WGBS processing workflows
- MOFA2 multi-omic integration

**Kathleen Durkin** - Regulatory Networks & Machine Learning
- miRNA discovery and target prediction
- WGCNA and network analyses
- Machine learning implementations
- Functional annotation expertise

**Steven Roberts** - Epigenomics & Multi-Omic Integration
- WGBS methylation analysis
- SNP calling and epigenetic annotation
- Protein annotation
- Multi-omic machine learning

**Ariana S. Huffmyer** - Metabolomics & Physiological Integration
- Metabolomics and lipidomics processing
- Physiological data integration
- Cross-platform correlation
- Multi-species comparative approaches

**Jill Ashey** - Network Analysis
- WGCNA for metabolomics
- Co-expression network construction

**R. Cunning** - Energetic-Epigenetic Integration
- Multi-platform correlation analysis
- Metabolite-lipid-methylation linkages

### 1.7 Strengths of the Repository

1. **Comprehensive Scope**: 6 integrated data types across 3 species
2. **Temporal Design**: Enables dynamic and predictive analyses
3. **Reproducibility**: Well-documented R Markdown workflows
4. **Advanced Methods**: MOFA2, WGCNA, machine learning, orthology
5. **Multiple Scales**: From individual genes to system-level networks
6. **Comparative Framework**: Species-specific vs. conserved mechanisms
7. **Biological Replication**: Multiple colonies per species
8. **Technical Quality**: Rigorous QC and processing standards
9. **Code Organization**: Clear naming conventions and structure
10. **Collaborative**: Multiple contributors with complementary expertise

### 1.8 Weaknesses and Opportunities for Improvement

#### Minor Weaknesses:

1. **Uneven Species Coverage**: D-Apul has 45 scripts vs. 13-14 for E-Peve and F-Ptua
   - *Opportunity*: Extend advanced analyses (MOFA2, ML) to other species
   
2. **Documentation Gaps**: Some scripts lack detailed biological interpretation
   - *Opportunity*: Add biological context sections to key analyses

3. **Cross-Species Integration**: Limited comparative analyses beyond orthology
   - *Opportunity*: MOFA2 or ML across all three species simultaneously

4. **Statistical Power**: Sample sizes vary by timepoint and species
   - *Opportunity*: Power analysis and effect size reporting

5. **Validation**: Limited experimental validation of predicted interactions
   - *Opportunity*: Wet lab validation of key miRNA targets or ML features

6. **Metadata Standardization**: Some metadata files in Excel format
   - *Opportunity*: Convert all to machine-readable CSV/TSV

#### Opportunities for Enhancement:

1. **Interactive Visualization**: Shiny apps for network exploration
2. **Pathway Analysis**: KEGG/Reactome pathway integration
3. **Evolutionary Analysis**: dN/dS ratios for orthologs
4. **Integration with Public Data**: Compare to other coral studies
5. **Time Series Modeling**: More sophisticated temporal models (e.g., splines, GAMs)
6. **Causal Inference**: Structural equation modeling or causal networks
7. **Cell Type Deconvolution**: If single-cell data available
8. **Microbiome Integration**: Symbiont or bacterial data if available

### 1.9 Scientific Value Assessment

**Overall Scientific Value**: **EXCEPTIONAL**

This repository represents a comprehensive molecular characterization of coral responses that is:
- **Technically Sophisticated**: Uses state-of-the-art bioinformatics methods
- **Biologically Relevant**: Addresses fundamental questions in coral biology
- **Methodologically Rigorous**: Appropriate statistical approaches and validation
- **Highly Integrative**: Combines 6 data types in meaningful ways
- **Reproducible**: Well-documented analytical workflows
- **Publishable**: Contains data/analyses for multiple high-impact publications

**Impact Potential**:
- High-impact coral biology journals (e.g., *Molecular Ecology*, *BMC Genomics*)
- Multi-omics methods journals (e.g., *NAR*, *Bioinformatics*)
- Ecological genomics journals (e.g., *Ecology Letters*, *Evolution*)
- Potential review/perspective pieces on coral molecular responses

---

## 2. MANUSCRIPT OUTLINES

The following four manuscript concepts range from straightforward to highly ambitious, providing strategic options for publication.

---

### MANUSCRIPT 1: Multi-Omic Integration Reveals Temporal Regulatory Networks in the Coral *Acropora pulchra*

**Effort Level**: ⭐ **EASY / LOW EFFORT**

**Estimated Timeline**: 2-3 months

**Rationale**: This manuscript leverages existing MOFA2 and WGCNA analyses already completed in the repository with minimal additional analysis required.

#### Abstract/Summary

Coral responses to environmental change involve coordinated molecular responses across multiple regulatory layers. Using a comprehensive multi-omic time series dataset from *Acropora pulchra*, we integrated gene expression, miRNA regulation, DNA methylation, and metabolic profiles using Multi-Omics Factor Analysis (MOFA2) and weighted gene co-expression network analysis (WGCNA). We identified X latent factors explaining Y% of variance across data types, with specific factors associated with temporal progression, metabolic state, and stress response. Co-expression network analysis revealed Z modules containing coordinated mRNA-miRNA-metabolite clusters. Key modules were enriched for energy metabolism, protein synthesis, and stress response pathways. Machine learning models trained on multi-omic features predicted physiological state with X% accuracy, identifying specific miRNAs and CpG sites as key predictive features. This work demonstrates the power of multi-omic integration for understanding coral molecular responses and provides a framework for predictive molecular monitoring.

#### Main Sections/Topics

1. **Introduction**
   - Coral responses to environmental change
   - Multi-omic approaches in marine biology
   - MOFA2 and WGCNA as integration methods
   - Study objectives

2. **Methods**
   - Sample collection and experimental design
   - Multi-omic data generation (RNA-seq, miRNA-seq, WGBS, metabolomics)
   - MOFA2 factor analysis
   - WGCNA co-expression networks
   - Machine learning predictive modeling
   - Functional annotation and enrichment

3. **Results**
   - Multi-omic data characteristics
   - MOFA2 factor identification and variance decomposition
   - Temporal patterns in factor loadings
   - WGCNA module identification (mRNA-miRNA-metabolite)
   - Module-trait associations with physiology
   - Functional enrichment of key modules
   - Machine learning feature importance
   - Predictive signatures of physiological state

4. **Discussion**
   - Integrated molecular responses in corals
   - Key regulatory modules and their functions
   - Predictive molecular biomarkers
   - Comparison to other coral studies
   - Applications for monitoring and conservation
   - Limitations and future directions

5. **Conclusions**

#### Figures (8-10 main + supplements)

1. Experimental design schematic
2. MOFA2 variance explained and factor heatmap
3. Temporal trajectories of key MOFA2 factors
4. WGCNA network dendrogram and module colors
5. Module-trait correlation heatmap
6. Key module network visualization (mRNA-miRNA-metabolite)
7. GO enrichment bubble plot for major modules
8. Machine learning performance and feature importance
9. Predictive signature validation
10. Integrative model of coral molecular response

#### Data Requirements

**Already Available**:
- Script 22.20: MOFA2 analysis complete
- Scripts 12, 16: WGCNA networks complete
- Script 21: Functional annotation complete
- Scripts 22-22.6: Machine learning models complete
- M-multi-species: Metabolomics data processed

**Minimal New Analysis**:
- Compile figures from existing outputs
- Write biological interpretation text
- Create integrative conceptual figure
- Statistical summaries for results section

#### Target Journals

- **Primary**: *NAR Genomics and Bioinformatics* (multi-omics focus)
- **Alternative 1**: *BMC Genomics* (broad scope)
- **Alternative 2**: *Molecular Ecology Resources* (methods emphasis)
- **Alternative 3**: *Coral Reefs* (biological focus)

#### Estimated Effort Breakdown

- **Literature review**: 1 week
- **Additional analysis**: 1-2 weeks (mostly figure compilation)
- **Writing first draft**: 2-3 weeks
- **Revisions and co-author input**: 3-4 weeks
- **Submission preparation**: 1 week

**Total**: 8-12 weeks

#### Expertise Needed

- **Lead Author**: K. Durkin or S. White (MOFA2 expertise)
- **Co-authors**: 
  - S. Roberts (epigenomics context)
  - A. Huffmyer (metabolomics/physiological integration)
  - All data contributors
- **Potential Collaborator**: MOFA2 developers for methods validation (optional)

#### Success Probability

**High (90%)** - Analysis complete, just needs narrative and figure polishing

---

### MANUSCRIPT 2: miRNA-Mediated Gene Regulation Networks Orchestrate Temporal Responses in Three Coral Species

**Effort Level**: ⭐⭐ **MODERATE EFFORT**

**Estimated Timeline**: 4-6 months

**Rationale**: Builds on existing miRNA discovery, target prediction, and co-expression analyses but requires cross-species synthesis and some additional comparative analysis.

#### Abstract/Summary

MicroRNAs (miRNAs) are key post-transcriptional regulators that fine-tune gene expression in response to environmental stimuli. Despite their importance in diverse organisms, miRNA regulatory networks remain poorly characterized in reef-building corals. We performed comprehensive miRNA discovery and target prediction in three coral species with distinct life history strategies: *Acropora pulchra* (fast-growing, branching), *Porites evermanni* (slow-growing, massive), and *Pocillopora tuahiniensis* (intermediate). We identified X conserved and Y species-specific miRNAs, predicted Z thousand miRNA-mRNA interactions, and validated W% through expression correlation and target site conservation. Co-expression network analysis revealed distinct regulatory modules associated with growth, metabolism, and stress response. Species comparison showed that fast-growing *A. pulchra* exhibited more dynamic miRNA regulation while slow-growing *P. evermanni* showed more stable expression patterns. Key miRNAs regulating heat shock proteins, metabolic enzymes, and cell cycle genes were identified. This work provides the first comparative miRNA regulatory atlas for corals and reveals species-specific regulatory strategies.

#### Main Sections/Topics

1. **Introduction**
   - miRNA biology and function
   - miRNAs in marine invertebrates
   - Gaps in coral miRNA knowledge
   - Species comparison rationale
   - Study objectives

2. **Methods**
   - Sample collection and design (3 species, 4 timepoints)
   - sRNA-seq library preparation and sequencing
   - miRNA discovery with ShortStack
   - Target prediction (RNAhybrid, miRanda)
   - Expression correlation networks
   - Cross-species conservation analysis
   - Orthology mapping of targets
   - Functional enrichment analysis

3. **Results**
   - miRNA discovery: counts and characteristics per species
   - Conservation vs. novelty of coral miRNAs
   - Target prediction statistics
   - Expression validation of miRNA-target pairs
   - Co-expression network modules
   - Species comparison of regulatory dynamics
   - Key regulatory modules:
     - Growth and cell proliferation
     - Energy metabolism
     - Stress response pathways
   - miRNA-lncRNA interactions
   - Temporal patterns of regulation

4. **Discussion**
   - Coral miRNA repertoire and evolution
   - Species-specific vs. conserved regulatory strategies
   - Life history correlates of miRNA dynamics
   - Key regulatory axes in coral biology
   - Comparison to other cnidarians
   - Potential for miRNA-based biomarkers
   - Future directions

5. **Conclusions**

#### Figures (10-12 main + supplements)

1. Study design and species comparison
2. miRNA discovery pipeline and statistics
3. Venn diagram: conserved vs. species-specific miRNAs
4. Target prediction workflow and validation
5. Co-expression network across species
6. Module preservation analysis across species
7. Key regulatory modules (network viz for each species)
8. Temporal dynamics of key miRNAs
9. Functional enrichment comparison
10. Species comparison of regulatory patterns
11. Case studies: heat shock, metabolism, cell cycle
12. Conceptual model of species-specific regulation

#### Data Requirements

**Already Available**:
- Scripts 04, 08: miRNA discovery for all species
- Scripts 06: Target prediction (D-Apul)
- Scripts 12, 14, 16: Co-expression networks (D-Apul)
- M-multi-species/11: Orthology for target comparison

**New Analysis Required**:
- Target prediction for E-Peve and F-Ptua
- Cross-species network comparison
- Module preservation statistics
- Conservation analysis of miRNA families
- Comparative enrichment analysis
- Integration with orthology data

**Estimated New Code**: 5-8 additional R Markdown scripts

#### Target Journals

- **Primary**: *Molecular Ecology* (comparative molecular biology)
- **Alternative 1**: *Genome Biology and Evolution* (evolutionary genomics)
- **Alternative 2**: *RNA Biology* (miRNA focus)
- **Alternative 3**: *BMC Genomics* (comparative genomics)

#### Estimated Effort Breakdown

- **Literature review**: 2 weeks
- **New analysis (target prediction, comparison)**: 4-6 weeks
- **Network visualization and statistics**: 2 weeks
- **Writing first draft**: 4 weeks
- **Revisions and co-author input**: 4-6 weeks
- **Submission preparation**: 1-2 weeks

**Total**: 17-24 weeks

#### Expertise Needed

- **Lead Author**: K. Durkin (miRNA expertise)
- **Co-authors**:
  - S. White (bioinformatics, pipeline development)
  - S. Roberts (comparative genomics)
  - miRNA expert consultant (for methods validation)
- **Potential Collaborators**:
  - Coral evolutionary biologist
  - Cnidarian miRNA specialist

#### Success Probability

**High (80%)** - Substantial existing data, moderate new analysis, timely topic

---

### MANUSCRIPT 3: Epigenetic and Metabolic Integration Across Coral Life History Strategies: A Comparative Multi-Omic Analysis

**Effort Level**: ⭐⭐⭐ **ADVANCED EFFORT**

**Estimated Timeline**: 8-12 months

**Rationale**: Requires substantial new cross-species integration, linking epigenetics (WGBS) with metabolomics/lipidomics across all three species, and developing novel analytical frameworks.

#### Abstract/Summary

Life history strategies in corals range from fast-growing, opportunistic species to slow-growing, stress-tolerant species, but the molecular underpinnings of these strategies remain poorly understood. We hypothesized that epigenetic regulation coordinates metabolic reprogramming to support distinct life history strategies in corals. Using whole-genome bisulfite sequencing (WGBS), RNA-seq, metabolomics, and lipidomics from three coral species spanning the life history spectrum, we tested whether DNA methylation patterns predict metabolic phenotypes and differ between fast and slow-growing species. We found that gene body methylation was strongly associated with constitutive expression of core metabolic genes, while promoter methylation correlated with facultative metabolic responses. Fast-growing *A. pulchra* showed more dynamic methylation changes and metabolic plasticity compared to slow-growing *P. evermanni*. Machine learning integration of methylation and metabolite data predicted growth rates with X% accuracy. Key CpG sites in metabolic enzyme genes (e.g., citrate synthase, lipid synthesis pathways) were identified as life history markers. Orthologous gene analysis revealed conserved methylation-metabolism linkages across species. This work demonstrates that epigenetic regulation of metabolism underlies fundamental life history trade-offs in corals.

#### Main Sections/Topics

1. **Introduction**
   - Life history theory and molecular mechanisms
   - Epigenetic regulation of metabolism
   - DNA methylation in invertebrates
   - Coral life history diversity
   - Metabolic demands of different growth strategies
   - Study objectives and hypotheses

2. **Methods**
   - Species selection and life history characterization
   - Sample collection across 4 timepoints
   - WGBS library prep and data processing
   - Metabolomics and lipidomics profiling
   - RNA-seq for expression context
   - CpG methylation analysis
   - Gene body vs. promoter methylation
   - Methylation-expression correlation
   - Methylation-metabolite correlation
   - Cross-species orthology mapping
   - Machine learning integration
   - Pathway and network analysis

3. **Results**
   - **Genome-wide methylation patterns**:
     - CpG distribution and methylation levels per species
     - Methylation in gene features (promoter, body, intergenic)
   - **Metabolic profiling**:
     - Metabolite and lipid diversity per species
     - Metabolic pathway representation
   - **Methylation-metabolism correlations**:
     - Gene-specific associations
     - Pathway-level patterns
   - **Life history comparisons**:
     - Fast vs. slow grower methylation dynamics
     - Metabolic plasticity differences
   - **Machine learning integration**:
     - Methylation-metabolite prediction models
     - Feature importance and key CpGs
   - **Orthologous gene analysis**:
     - Conserved vs. divergent methylation patterns
     - Evolution of methylation-metabolism linkages
   - **Temporal dynamics**:
     - Timepoint-specific patterns
     - Trajectory analysis

4. **Discussion**
   - Epigenetic basis of life history strategies
   - Gene body methylation and constitutive metabolism
   - Promoter methylation and metabolic plasticity
   - Fast vs. slow grower molecular strategies
   - Trade-offs: growth vs. stress tolerance
   - Conservation and divergence of mechanisms
   - Evolutionary implications
   - Applications for coral restoration
   - Limitations and future work

5. **Conclusions**

#### Figures (12-15 main + supplements)

1. Life history trait comparison (growth, metabolism)
2. WGBS pipeline and genome-wide methylation patterns
3. Metabolomics/lipidomics overview (PCA, heatmaps)
4. CpG methylation in gene features (body vs. promoter)
5. Methylation-expression correlation genome-wide
6. Methylation-metabolite correlation networks
7. Species comparison of methylation dynamics
8. Pathway-level methylation-metabolism integration
9. Machine learning model performance
10. Key CpG sites and metabolic genes
11. Ortholog methylation conservation analysis
12. Temporal trajectories by species
13. Life history trade-off schematic
14. Integrative conceptual model
15. Predictive biomarker panel

#### Data Requirements

**Already Available**:
- WGBS data for all 3 species (scripts 00.20, 00.21, etc.)
- Metabolomics data (M-multi-species/07)
- Lipidomics data (M-multi-species/06)
- RNA-seq expression (all species)
- Orthology (M-multi-species/11)
- Partial methylation-metabolism links (M-multi-species/08)

**New Analysis Required**:
- **Major new work**:
  - Cross-species WGBS harmonization and comparison
  - CpG-metabolite correlation across species
  - Life history trait quantification and modeling
  - Species-specific machine learning models
  - Comparative methylation dynamics analysis
  - Ortholog methylation pattern analysis
  - Pathway-level integration
- **Estimated new code**: 15-20 R Markdown scripts
- **Potential wet lab validation**: Targeted bisulfite PCR for key CpGs (optional)

#### Target Journals

- **Primary**: *Nature Ecology & Evolution* (high-impact, life history focus)
- **Alternative 1**: *Molecular Biology and Evolution* (evolutionary epigenetics)
- **Alternative 2**: *PNAS* (broad significance)
- **Alternative 3**: *Epigenetics & Chromatin* (mechanistic focus)
- **Alternative 4**: *eLife* (integrative biology)

#### Estimated Effort Breakdown

- **Literature review and hypothesis refinement**: 3-4 weeks
- **New bioinformatics analysis**: 10-14 weeks
- **Life history trait compilation**: 2-3 weeks
- **Machine learning model development**: 4-6 weeks
- **Statistical analysis and visualization**: 4-5 weeks
- **Writing first draft**: 6-8 weeks
- **Revisions and co-author input**: 6-8 weeks
- **Potential wet lab validation**: 8-12 weeks (optional)
- **Submission preparation**: 2 weeks

**Total**: 35-52 weeks (8-12 months)

#### Expertise Needed

- **Lead Authors**: S. Roberts + A. Huffmyer (co-first authors)
  - Roberts: Epigenomics expertise
  - Huffmyer: Metabolomics and physiological integration
- **Co-authors**:
  - S. White (bioinformatics infrastructure)
  - K. Durkin (machine learning)
  - All data generators
- **Potential Collaborators**:
  - Life history theorist
  - Metabolic physiologist
  - Evolutionary epigeneticist
  - Statistician for complex models

#### Challenges and Mitigation

**Challenges**:
1. **Sample size for cross-species comparison** - Ensure adequate power analysis
2. **Confounding factors** (symbiont density, microbiome) - Control in models or measure
3. **Causal inference** - Use careful language; consider Mendelian randomization approaches
4. **Computational intensity** - Use HPC resources for large-scale correlations
5. **Biological validation** - Consider targeted experiments or validation datasets

**Mitigation Strategies**:
- Comprehensive statistical modeling with covariates
- Sensitivity analyses and robustness checks
- Bootstrapping and permutation tests for significance
- Clear distinction between correlation and causation
- Collaboration with experimentalists for validation

#### Success Probability

**Moderate-High (70%)** - Ambitious scope but strong existing data foundation; requires significant new analysis and integration; high potential impact if successful

---

### MANUSCRIPT 4: Predictive Multi-Omic Network Models Reveal Conserved and Species-Specific Molecular Architectures of Coral Environmental Responses

**Effort Level**: ⭐⭐⭐⭐ **CHALLENGING / HIGH EFFORT**

**Estimated Timeline**: 12-18 months

**Rationale**: This is the most ambitious manuscript, requiring novel method development, comprehensive multi-species multi-omic integration, advanced network inference, experimental validation, and synthesis of all repository data.

#### Abstract/Summary

Understanding how organisms integrate environmental signals across molecular layers to produce phenotypic responses is a central challenge in systems biology. Corals provide an ideal system due to their sensitivity to environmental change and the availability of multi-omic data. We developed a novel integrative framework combining MOFA2, causal network inference, and deep learning to predict coral physiological responses from multi-omic signatures across three species. Using time-series data spanning transcriptomics, epigenomics, small RNA regulation, metabolomics, and lipidomics, we constructed predictive network models that achieve X% accuracy in forecasting growth, photosynthetic efficiency, and stress resilience. Our approach identified Y conserved regulatory modules shared across species and Z species-specific network motifs. Causal modeling revealed hierarchical regulatory cascades where DNA methylation → miRNA → mRNA → metabolite progressions mediate environmental responses. Experimental perturbation of key miRNAs validated predicted network edges. Comparative analysis identified orthologous genes with conserved vs. rewired network positions across species. We provide an interactive web platform for exploring these multi-omic networks (coralOmicsExplorer). This work establishes a generalizable framework for multi-omic prediction and reveals the molecular architecture underlying coral environmental responses.

#### Main Sections/Topics

1. **Introduction**
   - Systems biology of environmental responses
   - Multi-omic integration challenges
   - Network approaches to biological prediction
   - Coral system advantages
   - Limitations of current approaches
   - Study objectives and innovations

2. **Methods**
   - **Experimental Design**:
     - Three species, four timepoints, multi-omic sampling
     - Physiological measurements
   - **Multi-Omic Data Generation**:
     - RNA-seq, miRNA-seq, WGBS, metabolomics, lipidomics
     - Quality control and normalization
   - **Integration Framework**:
     - MOFA2 for dimensionality reduction
     - Multi-layer network construction
     - Causal inference (PC algorithm, Granger causality)
     - Deep learning architectures (graph neural networks)
   - **Cross-Species Modeling**:
     - Orthology-guided network alignment
     - Transfer learning across species
   - **Predictive Modeling**:
     - Train/test/validation splits
     - Feature selection and importance
     - Model ensemble approaches
   - **Network Analysis**:
     - Module detection algorithms
     - Centrality measures and hub identification
     - Network motif analysis
     - Conservation vs. rewiring
   - **Experimental Validation**:
     - miRNA perturbation experiments
     - Target validation by qPCR
     - Network edge confirmation
   - **Web Platform Development**:
     - Interactive network visualization
     - Data exploration tools

3. **Results**
   - **Multi-Omic Integration**:
     - Data characteristics across species
     - MOFA2 factor analysis (shared vs. specific)
     - Variance partitioning
   - **Network Architecture**:
     - Multi-layer network topology
     - Network statistics (density, clustering, modularity)
     - Hierarchical organization
   - **Predictive Performance**:
     - Accuracy for growth, photosynthesis, stress metrics
     - Feature importance and key predictors
     - Model interpretation
   - **Conserved Regulatory Modules**:
     - Shared modules across all species
     - Functional characterization
     - Hub genes/miRNAs/metabolites
   - **Species-Specific Networks**:
     - Unique modules per species
     - Life history correlates
     - Divergent regulatory logic
   - **Causal Regulatory Cascades**:
     - Methylation → miRNA → mRNA → metabolite
     - Temporal precedence analysis
     - Causal inference results
   - **Experimental Validation**:
     - miRNA perturbation outcomes
     - Target confirmation
     - Network edge validation rate
   - **Ortholog Network Comparison**:
     - Conserved vs. rewired network positions
     - Evolutionary insights
   - **CoralOmicsExplorer Platform**:
     - Features and functionality
     - User interface and access

4. **Discussion**
   - **Novel Framework Contributions**:
     - Advances in multi-omic integration
     - Causal network inference in non-model organisms
     - Transfer learning across species
   - **Biological Insights**:
     - Molecular architecture of coral responses
     - Conserved vs. species-specific strategies
     - Regulatory hierarchy and information flow
   - **Predictive Signatures**:
     - Biomarkers for coral health
     - Early warning indicators
   - **Evolutionary Perspectives**:
     - Network conservation and divergence
     - Life history strategy correlates
   - **Applications**:
     - Coral reef monitoring
     - Restoration target selection
     - Climate change adaptation
   - **Methodological Advances**:
     - Generalizable to other systems
     - Open-source tools provided
   - **Limitations**:
     - Sample size constraints
     - Causal inference assumptions
     - Validation scope
   - **Future Directions**:
     - Single-cell multi-omics
     - Long-term forecasting
     - Additional species

5. **Conclusions**

#### Figures (15-20 main + extensive supplements)

**Main Figures**:
1. Conceptual overview of integrative framework
2. Experimental design and data summary
3. MOFA2 multi-species integration
4. Multi-layer network architecture schematic
5. Predictive model performance across phenotypes
6. Feature importance and key predictors
7. Conserved regulatory modules (network viz)
8. Species-specific network modules
9. Causal regulatory cascades (Sankey diagram)
10. Hub gene/miRNA/metabolite analysis
11. Experimental validation results
12. Ortholog network comparison (alignment viz)
13. Network motif analysis
14. CoralOmicsExplorer interface screenshots
15. Integrative model of molecular response architecture

**Supplementary Figures** (20-30):
- Species-specific network details
- Additional validation experiments
- Module-by-module functional enrichment
- Temporal dynamics of key modules
- Sensitivity analyses
- Model comparisons
- Additional statistical analyses

#### Data Requirements

**Existing Data** (use ALL repository data):
- All RNA-seq, miRNA-seq, WGBS, metabolomics, lipidomics (3 species)
- All MOFA2, WGCNA, machine learning analyses
- Orthology data
- Physiological measurements

**New Analysis Required** (EXTENSIVE):

1. **Integration Framework Development**:
   - Cross-species MOFA2 with all data types
   - Multi-layer network construction algorithms
   - Causal inference pipeline development
   - Deep learning model architecture (graph neural networks)
   - Transfer learning implementation

2. **Validation Framework**:
   - Power analysis and sample size justification
   - Cross-validation schemes
   - Permutation testing
   - Bootstrapping for confidence intervals

3. **Network Analysis**:
   - Module detection across species
   - Network alignment algorithms
   - Motif enrichment analysis
   - Hub identification and characterization
   - Robustness analysis

4. **Evolutionary Analysis**:
   - Ortholog network positioning
   - dN/dS ratio calculation for orthologs
   - Network conservation scores
   - Rewiring statistics

5. **Web Platform**:
   - Shiny app development
   - Database backend
   - Interactive network visualization (D3.js or similar)
   - Documentation and tutorials

**Estimated New Code**: 30-40 R Markdown scripts + Shiny app code + Python scripts for deep learning

**Wet Lab Validation**:
- miRNA mimic/inhibitor experiments (2-3 key miRNAs)
- Target validation by qPCR (10-15 targets)
- Metabolite measurements post-perturbation
- ~6-8 months of experimental work

#### Target Journals

- **Primary**: *Nature* or *Science* (high-impact, broad significance)
- **Alternative 1**: *Nature Methods* (methodological innovation)
- **Alternative 2**: *Cell Systems* (systems biology focus)
- **Alternative 3**: *Nature Communications* (integrative biology)
- **Alternative 4**: *Science Advances* (multidisciplinary)

#### Estimated Effort Breakdown

- **Literature review and framework design**: 4-6 weeks
- **Method development (causal inference, deep learning)**: 12-16 weeks
- **Cross-species network construction**: 8-10 weeks
- **Predictive modeling and validation**: 8-10 weeks
- **Ortholog and evolutionary analysis**: 6-8 weeks
- **Wet lab validation experiments**: 24-32 weeks
- **Web platform development**: 12-16 weeks
- **Writing first draft**: 10-12 weeks
- **Revisions and co-author input**: 8-12 weeks
- **Submission preparation**: 2-3 weeks

**Total**: 52-78 weeks (12-18 months)

**Note**: Some tasks can be parallelized (e.g., analysis and wet lab validation)

#### Expertise Needed

**Core Team**:
- **Lead Author**: K. Durkin (regulatory networks, ML)
- **Co-Lead Authors**: 
  - S. Roberts (epigenomics, integration)
  - S. White (bioinformatics, infrastructure)
  - A. Huffmyer (metabolomics, physiology)

**Essential Collaborators**:
1. **Network/Systems Biologist**: Causal inference and network analysis expert
2. **Deep Learning Specialist**: Graph neural networks, transfer learning
3. **Statistician**: Complex modeling, causal inference validation
4. **Experimental Biologist**: miRNA perturbation and validation
5. **Web Developer**: Interactive platform development (or train existing team)
6. **Evolutionary Biologist**: Ortholog network interpretation

**Potential High-Profile Collaborators**:
- MOFA2 developers (for advanced applications)
- Leading coral genomicist
- Systems biology modeling expert

#### Technical Challenges and Solutions

**Challenges**:

1. **Computational Complexity**: Multi-species, multi-omic, temporal networks
   - *Solution*: High-performance computing, optimize algorithms, parallel processing

2. **Causal Inference in Observational Data**: No true experiments
   - *Solution*: Multiple causal inference methods, sensitivity analyses, experimental validation

3. **Deep Learning on Small Datasets**: Limited samples
   - *Solution*: Transfer learning, data augmentation, ensemble methods, regularization

4. **Cross-Species Network Alignment**: Orthology ambiguity
   - *Solution*: Multiple orthology criteria, weighted alignment, sensitivity testing

5. **Experimental Validation Logistics**: miRNA perturbation in corals
   - *Solution*: Collaborate with experimental lab, use coral cell lines if available, focus on key predictions

6. **Web Platform Maintenance**: Long-term hosting
   - *Solution*: Institutional hosting, Shiny server, or cloud deployment

7. **Reproducibility**: Complex pipeline
   - *Solution*: Docker containers, comprehensive documentation, code repositories

**Risk Mitigation**:
- Start with D-Apul only, expand to other species
- Modular approach: publish methods/platform separately if needed
- Regular consultation with collaborators
- Pilot experiments before full validation
- Plan B: observational validation using external datasets

#### Innovation and Impact

**Methodological Innovations**:
1. First multi-layer causal network inference in corals
2. Novel integration of MOFA2 with deep learning
3. Cross-species transfer learning framework
4. Interactive multi-omic network platform

**Biological Insights**:
1. Comprehensive molecular architecture of coral responses
2. Conserved vs. species-specific regulatory logic
3. Hierarchical information flow (methylation → miRNA → mRNA → metabolite)
4. Predictive molecular signatures for coral health

**Broader Impact**:
1. Framework generalizable to other systems
2. Open-source tools for marine biology community
3. Applications in coral conservation
4. Training in advanced bioinformatics methods

#### Success Probability

**Moderate (60%)** - Highly ambitious with significant technical and experimental challenges; requires substantial resources, expertise, and time; high risk but transformative potential if successful

---

## 3. SUMMARY COMPARISON OF MANUSCRIPTS

| Feature | MS1: Multi-Omic Integration | MS2: miRNA Networks | MS3: Epigenetic-Metabolic Integration | MS4: Predictive Network Models |
|---------|---------------------------|-------------------|-------------------------------------|------------------------------|
| **Effort** | ⭐ Easy | ⭐⭐ Moderate | ⭐⭐⭐ Advanced | ⭐⭐⭐⭐ Challenging |
| **Timeline** | 2-3 months | 4-6 months | 8-12 months | 12-18 months |
| **New Analysis** | Minimal | Moderate | Substantial | Extensive |
| **Wet Lab** | None | None | Optional | Essential |
| **Species** | A. pulchra | All three | All three | All three |
| **Data Types** | All 6 | RNA + miRNA | WGBS + metabolomics | All 6 + causal |
| **Methods Innovation** | Low | Moderate | High | Very High |
| **Target Journal Impact** | Moderate | High | Very High | Exceptional |
| **Success Probability** | 90% | 80% | 70% | 60% |
| **Lead Expertise** | Durkin/White | Durkin | Roberts/Huffmyer | Team + collaborators |
| **Collaborators Needed** | None | 1-2 | 2-4 | 5-7 |
| **Biological Novelty** | Moderate | High | Very High | Exceptional |
| **Technical Novelty** | Moderate | Moderate | High | Very High |

---

## 4. STRATEGIC RECOMMENDATIONS

### 4.1 Publication Strategy

**Recommended Sequence**:

1. **Start with MS1** (Easy) - 2-3 months
   - Quick win to establish publication momentum
   - Demonstrates repository value
   - Provides framework for subsequent papers
   - Can reference in grant applications

2. **Pursue MS2** (Moderate) in parallel with MS3 - 4-6 months after MS1
   - Different lead authors can work simultaneously
   - Builds on MS1 framework
   - Establishes species-comparison approach

3. **Develop MS3** (Advanced) - can overlap with MS2
   - More ambitious, requires focus
   - Complements MS2 (different data types)
   - High-impact target

4. **MS4** (Challenging) - long-term goal, 12-18 months
   - Most ambitious, treat as flagship
   - Can build on published MS1-3
   - Consider splitting into multiple papers if scope too large:
     - MS4a: Methods and framework
     - MS4b: Biological insights and applications

### 4.2 Resource Allocation

**Personnel**:
- MS1: 1 lead + 3 co-authors
- MS2: 1 lead + 4 co-authors + 1 consultant
- MS3: 2 co-leads + 5 co-authors + 2-3 collaborators
- MS4: 4 co-leads + all contributors + 5-7 collaborators

**Funding Needs**:
- MS1: Minimal (publication fees ~$2-3K)
- MS2: Low (~$5-10K for consultant, publication fees)
- MS3: Moderate (~$20-40K for potential validation, travel, open access fees)
- MS4: High (~$80-150K for wet lab validation, collaborators, platform development, high-impact open access)

**Computational Resources**:
- MS1-2: Standard lab computing
- MS3: HPC for large-scale correlations
- MS4: Significant HPC + cloud resources for deep learning

### 4.3 Additional Publication Opportunities

Beyond these four manuscripts, the repository supports:

1. **Technical/Methods Papers**:
   - "Multi-omic data processing pipeline for coral research"
   - "CoralOmicsExplorer: Interactive platform for coral molecular data"

2. **Dataset/Resource Papers**:
   - "A comprehensive multi-omic resource for three coral species" (Scientific Data)

3. **Review/Perspective**:
   - "Integrative omics approaches in coral biology: current state and future directions"

4. **Supplementary Analyses**:
   - lncRNA discovery and function (species-specific)
   - Ortholog evolution and selection (M-multi-species focus)
   - CpG island annotation across species

### 4.4 Data Sharing and Reproducibility

**Recommendations**:
1. Deposit raw data in public repositories (GEO, MetaboLights, etc.)
2. Create Zenodo repository for processed data and code
3. Develop comprehensive README and tutorials
4. Consider data paper for Scientific Data
5. Link manuscripts to GitHub repository releases
6. Ensure all figures reproducible from raw data

### 4.5 Community Engagement

**Opportunities**:
1. Present at coral biology conferences (ICRS, etc.)
2. Workshops on multi-omics in marine systems
3. Collaborate with coral restoration programs
4. Engage with climate change and conservation communities
5. Social media presence for broader impact

---

## 5. CONCLUSION

The **urol-e5/timeseries_molecular** repository represents an exceptional scientific resource containing comprehensive multi-omic data and analyses across three coral species. The depth and breadth of data, combined with sophisticated analytical approaches, provide a strong foundation for multiple high-impact publications.

**Key Takeaways**:

1. **Immediate Opportunity**: MS1 can be published quickly (2-3 months) with minimal additional work, providing immediate return on investment.

2. **Strategic Diversity**: The four proposed manuscripts span different effort levels, allowing flexibility based on available resources and timelines.

3. **Complementary Focus**: Each manuscript emphasizes different aspects (integration, regulation, epigenetics, prediction) while building on shared infrastructure.

4. **Scalable Ambition**: Can pursue conservative strategy (MS1-2) or ambitious strategy (all four), depending on resources and goals.

5. **Collaborative Opportunities**: Higher-tier manuscripts provide opportunities for strategic collaborations that enhance impact and expertise.

6. **Broader Impact**: Beyond publications, this work has applications in coral conservation, reef monitoring, and understanding climate change responses.

**Final Recommendation**: 
Begin with **MS1** immediately while planning for **MS2** and **MS3** in parallel. Consider **MS4** as a longer-term flagship project that can incorporate insights from the first three manuscripts. This strategy balances quick wins with ambitious long-term goals, maximizing the scientific impact of this valuable repository.

---

## Appendices

### Appendix A: Repository Statistics Summary

- **Total Scripts**: 85+ R Markdown files
- **Code Lines**: ~26,000+ (D-Apul alone)
- **Species**: 3 (A. pulchra, P. evermanni, P. tuahiniensis)
- **Timepoints**: 4 per species
- **Data Types**: 6 (RNA-seq, miRNA-seq, WGBS, lncRNA, metabolomics, lipidomics)
- **Outputs**: 1000+ figures and result files
- **Contributors**: 6 primary scientists
- **Analysis Categories**: 
  - QC/Processing (00-02): 15+ scripts
  - Expression (03): 3+ scripts
  - Small RNA (04-08): 10+ scripts
  - Networks (11-18): 8+ scripts
  - Machine Learning (20-22): 10+ scripts
  - Functional Annotation (21, 23, 25, 27): 5+ scripts
  - Epigenetics (14, 19, 26, 28-30, 36): 8+ scripts
  - Multi-species (M-): 13+ scripts

### Appendix B: Key Software and Tools Used

- **Alignment**: HiSat2, Bismark, Bowtie
- **Quantification**: StringTie, featureCounts, Salmon
- **DE Analysis**: DESeq2, limma
- **Small RNA**: ShortStack, miRanda, RNAhybrid
- **Networks**: WGCNA, igraph
- **Machine Learning**: randomForest, caret, glmnet
- **Multi-omic Integration**: MOFA2
- **Functional Annotation**: topGO, GOseq, BLAST
- **Statistics**: R (tidyverse, ggplot2, ComplexHeatmap)
- **Reproducibility**: R Markdown, knitr

### Appendix C: Suggested Timeline (Next 18 Months)

**Months 1-3**: MS1 completion and submission
**Months 2-7**: MS2 analysis and writing
**Months 4-12**: MS3 analysis and writing
**Months 6-18**: MS4 development (if pursued)
**Ongoing**: Data repository deposition, platform development

---

**Document Version**: 1.0  
**Last Updated**: October 2024  
**Contact**: Repository maintainers via GitHub issues

