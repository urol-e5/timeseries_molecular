27-Apul-mRNA-miRNA-interactions-topGO
================
Kathleen Durkin
2025-04-28

- <a href="#1-format-topgo-files" id="toc-1-format-topgo-files">1 Format
  topGO files</a>
  - <a href="#11-read-in-and-format-annotation-files"
    id="toc-11-read-in-and-format-annotation-files">1.1 Read in and format
    annotation files</a>
  - <a href="#12-set-up-gene2go-object"
    id="toc-12-set-up-gene2go-object">1.2 Set up gene2GO object</a>
  - <a href="#13-define-reference-set" id="toc-13-define-reference-set">1.3
    Define reference set</a>
  - <a href="#14-read-in-pccmiranda-data"
    id="toc-14-read-in-pccmiranda-data">1.4 Read in PCC/miranda data</a>
- <a href="#2-fa-of-all-mirna-targets"
  id="toc-2-fa-of-all-mirna-targets">2 FA of all miRNA targets</a>
- <a href="#3-fe-of-specific-mirnas-targets-all-targets"
  id="toc-3-fe-of-specific-mirnas-targets-all-targets">3 FE of specific
  miRNA’s targets (all targets)</a>
- <a href="#4-fe-of-specific-mirnas-targets-high-05-cor-targets"
  id="toc-4-fe-of-specific-mirnas-targets-high-05-cor-targets">4 FE of
  specific miRNA’s targets (high 0.5 cor targets)</a>
- <a href="#5-fe-of-specific-mirnas-targets-high-06-cor-targets"
  id="toc-5-fe-of-specific-mirnas-targets-high-06-cor-targets">5 FE of
  specific miRNA’s targets (high 0.6 cor targets)</a>
- <a href="#6-fe-of-specific-mirnas-targets-high-07-cor-targets"
  id="toc-6-fe-of-specific-mirnas-targets-high-07-cor-targets">6 FE of
  specific miRNA’s targets (high 0.7 cor targets)</a>

This script will use topGO to analyze functional enrichment of miRNA
targets for Apul

Code used below was created by `Jill Ashey`, modified for use with
A.pulchra datasets by `Kathleen Durkin`

# 1 Format topGO files

## 1.1 Read in and format annotation files

``` r
# Read in Apul annotations
annot_locations <- read.delim("https://raw.githubusercontent.com/urol-e5/deep-dive-expression/refs/heads/main/D-Apul/output/02-Apul-reference-annotation/Apulcra-genome-mRNA-IDmapping-2024_12_12.tab")
# Remove unneeded columns 
annot_locations <- annot_locations %>% dplyr::select(-X, -V13)
# Ensure there are no duplicate rows
annot_locations <- annot_locations %>% distinct()

head(annot_locations)
```

    ##                     V1     V3 Protein.names
    ## 1   ntLink_4:1155-1537 P35061   Histone H2A
    ## 2   ntLink_4:2660-3441 P84239    Histone H3
    ## 3   ntLink_4:4515-6830 P35061   Histone H2A
    ## 4   ntLink_4:7096-7859 P84239    Histone H3
    ## 5   ntLink_4:8474-9669 P35061   Histone H2A
    ## 6 ntLink_4:11162-11925 P84239    Histone H3
    ##                                     Organism Gene.Ontology..biological.process.
    ## 1          Acropora formosa (Staghorn coral)                                   
    ## 2 Urechis caupo (Innkeeper worm) (Spoonworm)                                   
    ## 3          Acropora formosa (Staghorn coral)                                   
    ## 4 Urechis caupo (Innkeeper worm) (Spoonworm)                                   
    ## 5          Acropora formosa (Staghorn coral)                                   
    ## 6 Urechis caupo (Innkeeper worm) (Spoonworm)                                   
    ##                                            Gene.Ontology.IDs
    ## 1 GO:0000786; GO:0003677; GO:0005634; GO:0030527; GO:0046982
    ## 2 GO:0000786; GO:0003677; GO:0005634; GO:0030527; GO:0046982
    ## 3 GO:0000786; GO:0003677; GO:0005634; GO:0030527; GO:0046982
    ## 4 GO:0000786; GO:0003677; GO:0005634; GO:0030527; GO:0046982
    ## 5 GO:0000786; GO:0003677; GO:0005634; GO:0030527; GO:0046982
    ## 6 GO:0000786; GO:0003677; GO:0005634; GO:0030527; GO:0046982

``` r
# Looks good!
```

This file shows each gene as it’s genomic location. We want to use gene
IDs to associate genes, so add gene IDs to this annotation table

Read in file that associates each mRNA genomic location with
corresponding gene ID

``` r
mRNA_FUNids <- read.table("../output/05-Apul-annotate-UTRs/Apul-mRNA-FUNids.txt", header=FALSE, col.names=c("location", "type", "mRNA_ID", "gene_ID", "product"), sep="\t")

# Remove unwanted text from parent column
mRNA_FUNids$gene_ID <- gsub("Parent=", "", mRNA_FUNids$gene_ID)
# Only need to keep mRNA location and gene ID
mRNA_FUNids <- mRNA_FUNids %>% dplyr::select(location, gene_ID)
```

join with annotation file

``` r
# join
annot <- left_join(annot_locations, mRNA_FUNids, by = c("V1" = "location"))

# ensure there are no duplicate rows
annot <- annot %>% distinct()
```

## 1.2 Set up gene2GO object

Want to isolate a list of GO terms per gene

``` r
gene2go <- annot %>% filter(!is.na(Gene.Ontology.IDs)) %>% dplyr::select(gene_ID, Gene.Ontology.IDs)
gene2go <- gene2go %>% dplyr::rename(GO.ID = Gene.Ontology.IDs)

gene2go_list <- setNames(
  strsplit(as.character(gene2go$GO.ID), ";"), 
  gene2go$gene_ID
)
```

Note: I think this means genes that had a Uniprot ID but no GO terms are
excluded from this analysis

## 1.3 Define reference set

Define reference set of genes. This should be all genes *found in our
samples*, NOT all genes in the A.pulchra genome. Some genes (e.g.,
reproduction pathways) may not be found/expected in our samples for
valid biological reasons.

``` r
# Read in counts matrix
Apul_counts <- read.csv("../output/02.20-D-Apul-RNAseq-alignment-HiSat2/apul-gene_count_matrix.csv")
# Exclude genes with all 0 counts
Apul_counts <- Apul_counts[rowSums(Apul_counts[, 2:6]) != 0, ]

# Select gene IDs of the genes present in our samples
all_genes <- Apul_counts$gene_id
length(all_genes)
```

    ## [1] 30094

So we have a reference set of 30094 genes present in our samples.

## 1.4 Read in PCC/miranda data

This is a table of all putative miRNA-mRNA binding predicted by miRanda,
plus Pearsons correlation coefficients for coexpression of each putative
binding pair. It only includes interactions with significant PCC (pval
\< 0.05)

``` r
data <- read.csv("../output/14.1-Apul-miRNA-mRNA-coexpression-additional_inputs/miRanda-PCC-significant-mRNA_3UTR_5UTR.csv") %>% filter(!is.na(mirna))
nrow(data)
```

    ## [1] 122274

``` r
head(data)
```

    ##           mirna                 Target Score Energy_Kcal_Mol Query_Aln
    ## 1 Cluster_10452   ntLink_0:61628-78213   147          -15.88      2 19
    ## 2 Cluster_10452   ntLink_0:61628-78213   147          -16.10      2 12
    ## 3 Cluster_10452 ntLink_4:151074-152269   149          -14.47      2 20
    ## 4 Cluster_10452 ntLink_4:151074-152269   159          -17.21      2 21
    ## 5 Cluster_10452 ntLink_4:187657-188852   149          -14.47      2 20
    ## 6 Cluster_10452 ntLink_4:187657-188852   159          -17.21      2 21
    ##   Subject_Aln Al_Len Subject_Identity Query_Identity         X4
    ## 1   1286 1308     19           63.16%         84.21% FUN_000006
    ## 2 13037 13058     10           90.00%         90.00% FUN_000006
    ## 3     900 920     18           66.67%         77.78% FUN_000113
    ## 4     211 234     21           66.67%         80.95% FUN_000113
    ## 5     900 920     18           66.67%         77.78% FUN_000128
    ## 6     211 234     21           66.67%         80.95% FUN_000128
    ##                  interaction      X         miRNA       mRNA   PCC.cor
    ## 1 Cluster_10452 _ FUN_000006 107896 Cluster_10452 FUN_000006 0.3314907
    ## 2 Cluster_10452 _ FUN_000006 107896 Cluster_10452 FUN_000006 0.3314907
    ## 3 Cluster_10452 _ FUN_000113 110497 Cluster_10452 FUN_000113 0.3960556
    ## 4 Cluster_10452 _ FUN_000113 110497 Cluster_10452 FUN_000113 0.3960556
    ## 5 Cluster_10452 _ FUN_000128 112792 Cluster_10452 FUN_000128 0.3471957
    ## 6 Cluster_10452 _ FUN_000128 112792 Cluster_10452 FUN_000128 0.3471957
    ##      p_value adjusted_p_value Source Significance
    ## 1 0.03665931        0.3965835   mRNA     p < 0.05
    ## 2 0.03665931        0.3965835   mRNA     p < 0.05
    ## 3 0.01141128        0.2339969   mRNA     p < 0.05
    ## 4 0.01141128        0.2339969   mRNA     p < 0.05
    ## 5 0.02816252        0.3539704   mRNA     p < 0.05
    ## 6 0.02816252        0.3539704   mRNA     p < 0.05

For the purposes of functional annotation and enrichment, I don’t care
if an miRNA may bind to several locations of a given gene, and I don’t
want those duplicates being overrepresented during enrichment analysis.

``` r
# Simplify target genes to gene IDs only (no binding coordinates)
data <- data %>% dplyr::select(-Target, -Score, -Energy_Kcal_Mol, -Subject_Aln, -Query_Aln, -Al_Len, -Subject_Identity, -Query_Identity) %>%
  distinct()
nrow(data)
```

    ## [1] 71634

# 2 FA of all miRNA targets

Functional annotation of all putative miRNA targets

``` r
cor_bind_FA <- left_join(data, annot, by = c("mRNA" = "gene_ID")) %>% distinct()

nrow(cor_bind_FA)
```

    ## [1] 71634

``` r
nrow(cor_bind_FA[!is.na(cor_bind_FA$Gene.Ontology.IDs),])
```

    ## [1] 20294

Of the 71634 putative miRNA targets predicted by miRanda with
significant PCC, 20294 have available annotations

``` r
high0.5_cor_bind_FA <- cor_bind_FA[abs(cor_bind_FA$PCC.cor) > 0.5,]

nrow(high0.5_cor_bind_FA)
```

    ## [1] 9480

``` r
nrow(high0.5_cor_bind_FA[!is.na(high0.5_cor_bind_FA$Gene.Ontology.IDs),])
```

    ## [1] 2711

Of the 9480 putative miRNA targets predicted by miRanda that that have
highly correlated expression (magnitude of correlation \> 0.5), 2711
have available annotations.

``` r
high0.6_cor_bind_FA <- cor_bind_FA[abs(cor_bind_FA$PCC.cor) > 0.6,]

nrow(high0.6_cor_bind_FA)
```

    ## [1] 2090

``` r
nrow(high0.6_cor_bind_FA[!is.na(high0.6_cor_bind_FA$Gene.Ontology.IDs),])
```

    ## [1] 601

2090 have correlation of at least 0.6, and of those 601 are annotated

``` r
high0.7_cor_bind_FA <- cor_bind_FA[abs(cor_bind_FA$PCC.cor) > 0.7,]

nrow(high0.7_cor_bind_FA)
```

    ## [1] 314

``` r
nrow(high0.7_cor_bind_FA[!is.na(high0.7_cor_bind_FA$Gene.Ontology.IDs),])
```

    ## [1] 86

314 have correlation of at least 0.7, and of those 86 are annotated

Save

``` r
write.csv(cor_bind_FA, "../output/27-Apul-mRNA-miRNA-interactions-topGO/miRNA_targets_FA.csv")
write.csv(high0.5_cor_bind_FA, "../output/27-Apul-mRNA-miRNA-interactions-topGO/miRNA_high0.5_cor_targets_FA.csv")
```

# 3 FE of specific miRNA’s targets (all targets)

Create topGO function for use with miRNA names

``` r
miRNA_topGO_FE <- function(miRNA.name, input_interactions) {

  #Isolate genes in our input module of interest
  interacting_genes <- input_interactions %>%
    filter(miRNA == miRNA.name) %>%
    pull(mRNA)
  
  if (length(interacting_genes) > 0) {
    # Create factor for all reference genes, where 1 represents module membership and 0 means the gene is not in module of interest
    gene_list <- factor(as.integer(all_genes %in% interacting_genes))
    names(gene_list) <- all_genes
    str(gene_list)

    ## Biological Process ##
    # Create topGO object
    GO_BP <- new("topGOdata", 
                ontology = "BP", # Biological Process
                allGenes = gene_list,
                annot = annFUN.gene2GO, 
                gene2GO = gene2go_list,
                geneSel=topDiffGenes)
    
    # Run GO enrichment test
    GO_BP_FE <- runTest(GO_BP, algorithm = "weight01", statistic = "fisher")
    # View the results
    GO_BP_results <- GenTable(GO_BP, Fisher = GO_BP_FE, orderBy = "Fisher",  topNodes = 100, numChar = 51)
    # Filter by significant results
    GO_BP_results$Fisher<-as.numeric(GO_BP_results$Fisher)
    GO_BP_results_sig<-GO_BP_results[GO_BP_results$Fisher<0.05,]
    
    
    ## Molecular Function ##
    # Create topGO object
    GO_MF <- new("topGOdata", 
                ontology = "MF", # Molecular Function
                allGenes = gene_list,
                annot = annFUN.gene2GO, 
                gene2GO = gene2go_list,
                geneSel=topDiffGenes)
    
    # Run GO enrichment test
    GO_MF_FE <- runTest(GO_MF, algorithm = "weight01", statistic = "fisher")
    # View the results
    GO_MF_results <- GenTable(GO_MF, Fisher = GO_MF_FE, orderBy = "Fisher",  topNodes = 100, numChar = 51)
    # Filter by significant results
    GO_MF_results$Fisher<-as.numeric(GO_MF_results$Fisher)
    GO_MF_results_sig<-GO_MF_results[GO_MF_results$Fisher<0.05,]
  
    # Return
    # Add type column only if results exist
    if (nrow(GO_BP_results_sig) > 0) {
      GO_BP_results_sig$type <- "Biological.Process"
    }
    if (nrow(GO_MF_results_sig) > 0) {
      GO_MF_results_sig$type <- "Molecular.Function"
    }
    GO_results <- rbind(GO_BP_results_sig, GO_MF_results_sig)
    print(GO_results)
  }
}

miRNA_topGO_FE("Cluster_9532", cor_bind_FA)
```

    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 161 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 9:   8 nodes to be scored    (11 eliminated genes)

    ## 
    ##   Level 8:   11 nodes to be scored   (28 eliminated genes)

    ## 
    ##   Level 7:   18 nodes to be scored   (33 eliminated genes)

    ## 
    ##   Level 6:   27 nodes to be scored   (175 eliminated genes)

    ## 
    ##   Level 5:   34 nodes to be scored   (587 eliminated genes)

    ## 
    ##   Level 4:   26 nodes to be scored   (661 eliminated genes)

    ## 
    ##   Level 3:   20 nodes to be scored   (953 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (1118 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1287 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 129 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   4 nodes to be scored    (100 eliminated genes)

    ## 
    ##   Level 7:   15 nodes to be scored   (407 eliminated genes)

    ## 
    ##   Level 6:   25 nodes to be scored   (491 eliminated genes)

    ## 
    ##   Level 5:   29 nodes to be scored   (828 eliminated genes)

    ## 
    ##   Level 4:   25 nodes to be scored   (1186 eliminated genes)

    ## 
    ##   Level 3:   19 nodes to be scored   (1816 eliminated genes)

    ## 
    ##   Level 2:   6 nodes to be scored    (2086 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2480 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0000724 double-strand break repair via homologous recombina...        13
    ## 2  GO:0000184 nuclear-transcribed mRNA catabolic process, nonsens...        10
    ## 3  GO:0000712      resolution of meiotic recombination intermediates         1
    ## 4  GO:0006813                                potassium ion transport         1
    ## 5  GO:0000423                                              mitophagy         2
    ## 6  GO:0000028                       ribosomal small subunit assembly         3
    ## 7  GO:0002933                                    lipid hydroxylation         3
    ## 8  GO:0005524                                            ATP binding       263
    ## 9  GO:0004108                         citrate (Si)-synthase activity         1
    ## 10 GO:0004013                        adenosylhomocysteinase activity         1
    ## 11 GO:0004185                  serine-type carboxypeptidase activity         2
    ## 12 GO:0004089                         carbonate dehydratase activity         2
    ## 13 GO:0004517                         nitric-oxide synthase activity         2
    ## 14 GO:0005201            extracellular matrix structural constituent        17
    ##    Significant Expected  Fisher               type
    ## 1            3     0.20 0.00085 Biological.Process
    ## 2            2     0.16 0.00969 Biological.Process
    ## 3            1     0.02 0.01559 Biological.Process
    ## 4            1     0.02 0.01559 Biological.Process
    ## 5            1     0.03 0.03095 Biological.Process
    ## 6            1     0.05 0.04608 Biological.Process
    ## 7            1     0.05 0.04608 Biological.Process
    ## 8           12     5.02 0.00310 Molecular.Function
    ## 9            1     0.02 0.01910 Molecular.Function
    ## 10           1     0.02 0.01910 Molecular.Function
    ## 11           1     0.04 0.03780 Molecular.Function
    ## 12           1     0.04 0.03780 Molecular.Function
    ## 13           1     0.04 0.03780 Molecular.Function
    ## 14           2     0.32 0.04050 Molecular.Function

Loop through all miRNA and run functional enrichment on the miRNA’s
targets (all predicted targets)

``` r
interacting_miRNAs <- unique(cor_bind_FA$miRNA) %>% na.omit
results_all_targets <- NULL  # initialize empty df

for(miRNA in interacting_miRNAs) {
  
  # Run topGO enrichment function
  miRNA_results <- miRNA_topGO_FE(miRNA, cor_bind_FA)
  
  # Only keep results if not empty
  if (nrow(miRNA_results) > 0) {
    
    # Add the miRNA source column
    miRNA_results$miRNA <- miRNA

    # Bind to the accumulating results data frame
    results_all_targets <- rbind(results_all_targets, miRNA_results)
  }
}
```

    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 300 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  3 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 12:  6 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 11:  11 nodes to be scored   (159 eliminated genes)

    ## 
    ##   Level 10:  15 nodes to be scored   (214 eliminated genes)

    ## 
    ##   Level 9:   24 nodes to be scored   (238 eliminated genes)

    ## 
    ##   Level 8:   22 nodes to be scored   (284 eliminated genes)

    ## 
    ##   Level 7:   30 nodes to be scored   (315 eliminated genes)

    ## 
    ##   Level 6:   44 nodes to be scored   (464 eliminated genes)

    ## 
    ##   Level 5:   56 nodes to be scored   (656 eliminated genes)

    ## 
    ##   Level 4:   44 nodes to be scored   (813 eliminated genes)

    ## 
    ##   Level 3:   32 nodes to be scored   (1077 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (1307 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1388 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 183 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   8 nodes to be scored    (9 eliminated genes)

    ## 
    ##   Level 8:   16 nodes to be scored   (144 eliminated genes)

    ## 
    ##   Level 7:   25 nodes to be scored   (444 eliminated genes)

    ## 
    ##   Level 6:   35 nodes to be scored   (527 eliminated genes)

    ## 
    ##   Level 5:   37 nodes to be scored   (829 eliminated genes)

    ## 
    ##   Level 4:   33 nodes to be scored   (1046 eliminated genes)

    ## 
    ##   Level 3:   17 nodes to be scored   (1763 eliminated genes)

    ## 
    ##   Level 2:   7 nodes to be scored    (2014 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2462 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0000082                  G1/S transition of mitotic cell cycle       104
    ## 2  GO:0001172                            RNA-templated transcription         1
    ## 3  GO:0002230 positive regulation of defense response to virus by...         1
    ## 4  GO:0001913                           T cell mediated cytotoxicity         1
    ## 5  GO:0004176                       ATP-dependent peptidase activity         2
    ## 6  GO:0003964                   RNA-directed DNA polymerase activity       122
    ## 7  GO:0004334                           fumarylacetoacetase activity         1
    ## 8  GO:0004862       cAMP-dependent protein kinase inhibitor activity         1
    ## 9  GO:0008137               NADH dehydrogenase (ubiquinone) activity         1
    ## 10 GO:0036381 pyridoxal 5'-phosphate synthase (glutamine hydrolys...         1
    ## 11 GO:0000254                      C-4 methylsterol oxidase activity         1
    ## 12 GO:0004129                          cytochrome-c oxidase activity         1
    ##    Significant Expected  Fisher               type
    ## 1            9     2.65 0.00080 Biological.Process
    ## 2            1     0.03 0.02550 Biological.Process
    ## 3            1     0.03 0.02550 Biological.Process
    ## 4            1     0.03 0.02550 Biological.Process
    ## 5            2     0.06 0.00099 Molecular.Function
    ## 6           11     3.87 0.00135 Molecular.Function
    ## 7            1     0.03 0.03170 Molecular.Function
    ## 8            1     0.03 0.03170 Molecular.Function
    ## 9            1     0.03 0.03170 Molecular.Function
    ## 10           1     0.03 0.03170 Molecular.Function
    ## 11           1     0.03 0.03170 Molecular.Function
    ## 12           1     0.03 0.03170 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 312 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  5 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 12:  6 nodes to be scored    (17 eliminated genes)

    ## 
    ##   Level 11:  10 nodes to be scored   (162 eliminated genes)

    ## 
    ##   Level 10:  17 nodes to be scored   (210 eliminated genes)

    ## 
    ##   Level 9:   24 nodes to be scored   (213 eliminated genes)

    ## 
    ##   Level 8:   26 nodes to be scored   (260 eliminated genes)

    ## 
    ##   Level 7:   39 nodes to be scored   (305 eliminated genes)

    ## 
    ##   Level 6:   51 nodes to be scored   (430 eliminated genes)

    ## 
    ##   Level 5:   59 nodes to be scored   (609 eliminated genes)

    ## 
    ##   Level 4:   36 nodes to be scored   (735 eliminated genes)

    ## 
    ##   Level 3:   25 nodes to be scored   (1003 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (1245 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1361 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 168 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 9:   10 nodes to be scored   (1 eliminated genes)

    ## 
    ##   Level 8:   11 nodes to be scored   (147 eliminated genes)

    ## 
    ##   Level 7:   21 nodes to be scored   (464 eliminated genes)

    ## 
    ##   Level 6:   30 nodes to be scored   (573 eliminated genes)

    ## 
    ##   Level 5:   33 nodes to be scored   (850 eliminated genes)

    ## 
    ##   Level 4:   26 nodes to be scored   (1175 eliminated genes)

    ## 
    ##   Level 3:   20 nodes to be scored   (1905 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (2100 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2542 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0006886                        intracellular protein transport         5
    ## 2 GO:0002225 positive regulation of antimicrobial peptide produc...         1
    ## 3 GO:0000289       nuclear-transcribed mRNA poly(A) tail shortening         2
    ## 4 GO:0001578                           microtubule bundle formation         2
    ## 5 GO:0003844             1,4-alpha-glucan branching enzyme activity         1
    ## 6 GO:0005220 inositol 1,4,5-trisphosphate-sensitive calcium-rele...         1
    ## 7 GO:0004335                                 galactokinase activity         1
    ## 8 GO:0004467              long-chain fatty acid-CoA ligase activity         1
    ## 9 GO:0004997        thyrotropin-releasing hormone receptor activity         1
    ##   Significant Expected Fisher               type
    ## 1           2     0.10 0.0037 Biological.Process
    ## 2           1     0.02 0.0198 Biological.Process
    ## 3           1     0.04 0.0393 Biological.Process
    ## 4           1     0.04 0.0393 Biological.Process
    ## 5           1     0.03 0.0300 Molecular.Function
    ## 6           1     0.03 0.0300 Molecular.Function
    ## 7           1     0.03 0.0300 Molecular.Function
    ## 8           1     0.03 0.0300 Molecular.Function
    ## 9           1     0.03 0.0300 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 206 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  2 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 11:  2 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 10:  3 nodes to be scored    (70 eliminated genes)

    ## 
    ##   Level 9:   9 nodes to be scored    (75 eliminated genes)

    ## 
    ##   Level 8:   14 nodes to be scored   (116 eliminated genes)

    ## 
    ##   Level 7:   24 nodes to be scored   (171 eliminated genes)

    ## 
    ##   Level 6:   33 nodes to be scored   (298 eliminated genes)

    ## 
    ##   Level 5:   40 nodes to be scored   (661 eliminated genes)

    ## 
    ##   Level 4:   37 nodes to be scored   (821 eliminated genes)

    ## 
    ##   Level 3:   29 nodes to be scored   (1014 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (1182 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1342 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 133 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   6 nodes to be scored    (101 eliminated genes)

    ## 
    ##   Level 7:   14 nodes to be scored   (389 eliminated genes)

    ## 
    ##   Level 6:   23 nodes to be scored   (475 eliminated genes)

    ## 
    ##   Level 5:   30 nodes to be scored   (781 eliminated genes)

    ## 
    ##   Level 4:   26 nodes to be scored   (1057 eliminated genes)

    ## 
    ##   Level 3:   18 nodes to be scored   (1724 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (2098 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2467 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0000349 generation of catalytic spliceosome for first trans...         1
    ## 2  GO:0000390                       spliceosomal complex disassembly         1
    ## 3  GO:0002218                   activation of innate immune response        15
    ## 4  GO:0001502                                 cartilage condensation         2
    ## 5  GO:0001895                                     retina homeostasis         2
    ## 6  GO:0001503                                           ossification        42
    ## 7  GO:0003964                   RNA-directed DNA polymerase activity       122
    ## 8  GO:0047869          dimethylpropiothetin dethiomethylase activity        19
    ## 9  GO:0004095              carnitine O-palmitoyltransferase activity         1
    ## 10 GO:0005436                    sodium:phosphate symporter activity         1
    ## 11 GO:0004013                        adenosylhomocysteinase activity         1
    ## 12 GO:0004467              long-chain fatty acid-CoA ligase activity         1
    ## 13 GO:0000703 oxidized pyrimidine nucleobase lesion DNA N-glycosy...         1
    ## 14 GO:0000048                           peptidyltransferase activity         1
    ##    Significant Expected  Fisher               type
    ## 1            1     0.02 0.02200 Biological.Process
    ## 2            1     0.02 0.02200 Biological.Process
    ## 3            2     0.33 0.04100 Biological.Process
    ## 4            1     0.04 0.04300 Biological.Process
    ## 5            1     0.04 0.04300 Biological.Process
    ## 6            3     0.92 0.04500 Biological.Process
    ## 7           11     3.65 0.00082 Molecular.Function
    ## 8            4     0.57 0.00204 Molecular.Function
    ## 9            1     0.03 0.02990 Molecular.Function
    ## 10           1     0.03 0.02990 Molecular.Function
    ## 11           1     0.03 0.02990 Molecular.Function
    ## 12           1     0.03 0.02990 Molecular.Function
    ## 13           1     0.03 0.02990 Molecular.Function
    ## 14           1     0.03 0.02990 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 249 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  5 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 11:  6 nodes to be scored    (144 eliminated genes)

    ## 
    ##   Level 10:  7 nodes to be scored    (211 eliminated genes)

    ## 
    ##   Level 9:   15 nodes to be scored   (218 eliminated genes)

    ## 
    ##   Level 8:   18 nodes to be scored   (259 eliminated genes)

    ## 
    ##   Level 7:   27 nodes to be scored   (332 eliminated genes)

    ## 
    ##   Level 6:   41 nodes to be scored   (459 eliminated genes)

    ## 
    ##   Level 5:   52 nodes to be scored   (648 eliminated genes)

    ## 
    ##   Level 4:   36 nodes to be scored   (783 eliminated genes)

    ## 
    ##   Level 3:   27 nodes to be scored   (1052 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (1177 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1342 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 135 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   5 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   6 nodes to be scored    (101 eliminated genes)

    ## 
    ##   Level 7:   13 nodes to be scored   (390 eliminated genes)

    ## 
    ##   Level 6:   25 nodes to be scored   (475 eliminated genes)

    ## 
    ##   Level 5:   32 nodes to be scored   (781 eliminated genes)

    ## 
    ##   Level 4:   27 nodes to be scored   (1075 eliminated genes)

    ## 
    ##   Level 3:   16 nodes to be scored   (1628 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (2159 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2443 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0007035                                 vacuolar acidification         1
    ## 2  GO:0000349 generation of catalytic spliceosome for first trans...         1
    ## 3  GO:0000390                       spliceosomal complex disassembly         1
    ## 4  GO:0002218                   activation of innate immune response        15
    ## 5  GO:0001895                                     retina homeostasis         2
    ## 6  GO:0001503                                           ossification        42
    ## 7  GO:0003964                   RNA-directed DNA polymerase activity       122
    ## 8  GO:0043531                                            ADP binding         1
    ## 9  GO:0005436                    sodium:phosphate symporter activity         1
    ## 10 GO:0004013                        adenosylhomocysteinase activity         1
    ## 11 GO:0004435          phosphatidylinositol phospholipase C activity         1
    ## 12 GO:0004467              long-chain fatty acid-CoA ligase activity         1
    ## 13 GO:0000703 oxidized pyrimidine nucleobase lesion DNA N-glycosy...         1
    ## 14 GO:0000048                           peptidyltransferase activity         1
    ## 15 GO:0003824                                     catalytic activity      1042
    ## 16 GO:0009034                                 tryptophanase activity         2
    ## 17 GO:0016829                                         lyase activity        58
    ##    Significant Expected Fisher               type
    ## 1            1     0.02 0.0230 Biological.Process
    ## 2            1     0.02 0.0230 Biological.Process
    ## 3            1     0.02 0.0230 Biological.Process
    ## 4            2     0.34 0.0440 Biological.Process
    ## 5            1     0.05 0.0450 Biological.Process
    ## 6            3     0.95 0.0480 Biological.Process
    ## 7            9     3.03 0.0028 Molecular.Function
    ## 8            1     0.02 0.0249 Molecular.Function
    ## 9            1     0.02 0.0249 Molecular.Function
    ## 10           1     0.02 0.0249 Molecular.Function
    ## 11           1     0.02 0.0249 Molecular.Function
    ## 12           1     0.02 0.0249 Molecular.Function
    ## 13           1     0.02 0.0249 Molecular.Function
    ## 14           1     0.02 0.0249 Molecular.Function
    ## 15          31    25.90 0.0301 Molecular.Function
    ## 16           1     0.05 0.0491 Molecular.Function
    ## 17           5     1.44 0.0493 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 209 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  5 nodes to be scored    (156 eliminated genes)

    ## 
    ##   Level 10:  6 nodes to be scored    (210 eliminated genes)

    ## 
    ##   Level 9:   11 nodes to be scored   (217 eliminated genes)

    ## 
    ##   Level 8:   15 nodes to be scored   (258 eliminated genes)

    ## 
    ##   Level 7:   23 nodes to be scored   (307 eliminated genes)

    ## 
    ##   Level 6:   35 nodes to be scored   (428 eliminated genes)

    ## 
    ##   Level 5:   42 nodes to be scored   (651 eliminated genes)

    ## 
    ##   Level 4:   30 nodes to be scored   (749 eliminated genes)

    ## 
    ##   Level 3:   25 nodes to be scored   (1013 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (1138 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1342 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 124 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   4 nodes to be scored    (100 eliminated genes)

    ## 
    ##   Level 7:   12 nodes to be scored   (384 eliminated genes)

    ## 
    ##   Level 6:   22 nodes to be scored   (464 eliminated genes)

    ## 
    ##   Level 5:   28 nodes to be scored   (721 eliminated genes)

    ## 
    ##   Level 4:   24 nodes to be scored   (998 eliminated genes)

    ## 
    ##   Level 3:   19 nodes to be scored   (1537 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (1927 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2560 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0000390                       spliceosomal complex disassembly         1
    ## 2  GO:0002218                   activation of innate immune response        15
    ## 3  GO:0001895                                     retina homeostasis         2
    ## 4  GO:0001503                                           ossification        42
    ## 5  GO:0003964                   RNA-directed DNA polymerase activity       122
    ## 6  GO:0000993                      RNA polymerase II complex binding         1
    ## 7  GO:0043531                                            ADP binding         1
    ## 8  GO:0004013                        adenosylhomocysteinase activity         1
    ## 9  GO:0004467              long-chain fatty acid-CoA ligase activity         1
    ## 10 GO:0000703 oxidized pyrimidine nucleobase lesion DNA N-glycosy...         1
    ## 11 GO:0000048                           peptidyltransferase activity         1
    ##    Significant Expected Fisher               type
    ## 1            1     0.02 0.0210 Biological.Process
    ## 2            2     0.32 0.0390 Biological.Process
    ## 3            1     0.04 0.0420 Biological.Process
    ## 4            3     0.89 0.0420 Biological.Process
    ## 5            9     3.12 0.0034 Molecular.Function
    ## 6            1     0.03 0.0256 Molecular.Function
    ## 7            1     0.03 0.0256 Molecular.Function
    ## 8            1     0.03 0.0256 Molecular.Function
    ## 9            1     0.03 0.0256 Molecular.Function
    ## 10           1     0.03 0.0256 Molecular.Function
    ## 11           1     0.03 0.0256 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 2 1 1 1 2 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 239 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  5 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  7 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 10:  9 nodes to be scored    (217 eliminated genes)

    ## 
    ##   Level 9:   13 nodes to be scored   (234 eliminated genes)

    ## 
    ##   Level 8:   19 nodes to be scored   (292 eliminated genes)

    ## 
    ##   Level 7:   21 nodes to be scored   (326 eliminated genes)

    ## 
    ##   Level 6:   37 nodes to be scored   (449 eliminated genes)

    ## 
    ##   Level 5:   51 nodes to be scored   (619 eliminated genes)

    ## 
    ##   Level 4:   39 nodes to be scored   (792 eliminated genes)

    ## 
    ##   Level 3:   26 nodes to be scored   (1061 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (1145 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1380 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 147 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   6 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 8:   9 nodes to be scored    (101 eliminated genes)

    ## 
    ##   Level 7:   20 nodes to be scored   (407 eliminated genes)

    ## 
    ##   Level 6:   27 nodes to be scored   (503 eliminated genes)

    ## 
    ##   Level 5:   31 nodes to be scored   (841 eliminated genes)

    ## 
    ##   Level 4:   29 nodes to be scored   (1172 eliminated genes)

    ## 
    ##   Level 3:   15 nodes to be scored   (1907 eliminated genes)

    ## 
    ##   Level 2:   6 nodes to be scored    (2185 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2426 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0002042      cell migration involved in sprouting angiogenesis         7
    ## 2 GO:0000103                                   sulfate assimilation         1
    ## 3 GO:0003964                   RNA-directed DNA polymerase activity       122
    ## 4 GO:0003676                                   nucleic acid binding       497
    ## 5 GO:0004829                         threonine-tRNA ligase activity         1
    ## 6 GO:0001002 RNA polymerase III type 1 promoter sequence-specifi...         1
    ## 7 GO:0003980 UDP-glucose:glycoprotein glucosyltransferase activi...         1
    ## 8 GO:0004134                    4-alpha-glucanotransferase activity         1
    ## 9 GO:0000703 oxidized pyrimidine nucleobase lesion DNA N-glycosy...         1
    ##   Significant Expected Fisher               type
    ## 1           2     0.16 0.0098 Biological.Process
    ## 2           1     0.02 0.0227 Biological.Process
    ## 3           9     3.16 0.0037 Molecular.Function
    ## 4          14    12.89 0.0223 Molecular.Function
    ## 5           1     0.03 0.0259 Molecular.Function
    ## 6           1     0.03 0.0259 Molecular.Function
    ## 7           1     0.03 0.0259 Molecular.Function
    ## 8           1     0.03 0.0259 Molecular.Function
    ## 9           1     0.03 0.0259 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 147 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  4 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 10:  6 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 9:   11 nodes to be scored   (151 eliminated genes)

    ## 
    ##   Level 8:   10 nodes to be scored   (163 eliminated genes)

    ## 
    ##   Level 7:   12 nodes to be scored   (179 eliminated genes)

    ## 
    ##   Level 6:   21 nodes to be scored   (265 eliminated genes)

    ## 
    ##   Level 5:   27 nodes to be scored   (487 eliminated genes)

    ## 
    ##   Level 4:   26 nodes to be scored   (567 eliminated genes)

    ## 
    ##   Level 3:   19 nodes to be scored   (878 eliminated genes)

    ## 
    ##   Level 2:   7 nodes to be scored    (1057 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1305 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 88 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   6 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 8:   5 nodes to be scored    (101 eliminated genes)

    ## 
    ##   Level 7:   12 nodes to be scored   (396 eliminated genes)

    ## 
    ##   Level 6:   17 nodes to be scored   (414 eliminated genes)

    ## 
    ##   Level 5:   17 nodes to be scored   (584 eliminated genes)

    ## 
    ##   Level 4:   12 nodes to be scored   (826 eliminated genes)

    ## 
    ##   Level 3:   11 nodes to be scored   (1432 eliminated genes)

    ## 
    ##   Level 2:   4 nodes to be scored    (1666 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2115 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0015969             guanosine tetraphosphate metabolic process         9
    ## 2 GO:0000122 negative regulation of transcription by RNA polymer...       140
    ## 3 GO:0003747                    translation release factor activity         1
    ## 4 GO:0001002 RNA polymerase III type 1 promoter sequence-specifi...         1
    ## 5 GO:0004315     3-oxoacyl-[acyl-carrier-protein] synthase activity         2
    ## 6 GO:0003676                                   nucleic acid binding       497
    ## 7 GO:0004622                             lysophospholipase activity         4
    ## 8 GO:0000340                      RNA 7-methylguanosine cap binding         4
    ##   Significant Expected Fisher               type
    ## 1           2     0.08 0.0027 Biological.Process
    ## 2           5     1.29 0.0060 Biological.Process
    ## 3           1     0.01 0.0110 Molecular.Function
    ## 4           1     0.01 0.0110 Molecular.Function
    ## 5           1     0.02 0.0220 Molecular.Function
    ## 6          14     5.55 0.0420 Molecular.Function
    ## 7           1     0.04 0.0440 Molecular.Function
    ## 8           1     0.04 0.0440 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 357 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  5 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  3 nodes to be scored    (11 eliminated genes)

    ## 
    ##   Level 11:  7 nodes to be scored    (160 eliminated genes)

    ## 
    ##   Level 10:  12 nodes to be scored   (203 eliminated genes)

    ## 
    ##   Level 9:   22 nodes to be scored   (213 eliminated genes)

    ## 
    ##   Level 8:   27 nodes to be scored   (308 eliminated genes)

    ## 
    ##   Level 7:   44 nodes to be scored   (359 eliminated genes)

    ## 
    ##   Level 6:   63 nodes to be scored   (524 eliminated genes)

    ## 
    ##   Level 5:   76 nodes to be scored   (806 eliminated genes)

    ## 
    ##   Level 4:   50 nodes to be scored   (968 eliminated genes)

    ## 
    ##   Level 3:   36 nodes to be scored   (1187 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (1334 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1402 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 288 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  6 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 9:   18 nodes to be scored   (13 eliminated genes)

    ## 
    ##   Level 8:   26 nodes to be scored   (151 eliminated genes)

    ## 
    ##   Level 7:   39 nodes to be scored   (500 eliminated genes)

    ## 
    ##   Level 6:   52 nodes to be scored   (619 eliminated genes)

    ## 
    ##   Level 5:   54 nodes to be scored   (942 eliminated genes)

    ## 
    ##   Level 4:   49 nodes to be scored   (1357 eliminated genes)

    ## 
    ##   Level 3:   28 nodes to be scored   (2029 eliminated genes)

    ## 
    ##   Level 2:   11 nodes to be scored   (2359 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2634 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0009734                      auxin-activated signaling pathway         4
    ## 2  GO:0000103                                   sulfate assimilation         1
    ## 3  GO:0001731         formation of translation preinitiation complex         1
    ## 4  GO:0001922                                 B-1 B cell homeostasis         1
    ## 5  GO:0000390                       spliceosomal complex disassembly         1
    ## 6  GO:0002230 positive regulation of defense response to virus by...         1
    ## 7  GO:0006884                                cell volume homeostasis         1
    ## 8  GO:0005543                                   phospholipid binding         8
    ## 9  GO:0000987    cis-regulatory region sequence-specific DNA binding       112
    ## 10 GO:0003676                                   nucleic acid binding       497
    ## 11 GO:0005044                            scavenger receptor activity        13
    ##    Significant Expected  Fisher               type
    ## 1            2     0.18 0.01100 Biological.Process
    ## 2            1     0.05 0.04500 Biological.Process
    ## 3            1     0.05 0.04500 Biological.Process
    ## 4            1     0.05 0.04500 Biological.Process
    ## 5            1     0.05 0.04500 Biological.Process
    ## 6            1     0.05 0.04500 Biological.Process
    ## 7            1     0.05 0.04500 Biological.Process
    ## 8            4     0.52 0.00099 Molecular.Function
    ## 9            8     7.22 0.03028 Molecular.Function
    ## 10          36    32.05 0.03676 Molecular.Function
    ## 11           3     0.84 0.04668 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 427 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  5 nodes to be scored    (2 eliminated genes)

    ## 
    ##   Level 11:  12 nodes to be scored   (146 eliminated genes)

    ## 
    ##   Level 10:  21 nodes to be scored   (207 eliminated genes)

    ## 
    ##   Level 9:   33 nodes to be scored   (224 eliminated genes)

    ## 
    ##   Level 8:   38 nodes to be scored   (306 eliminated genes)

    ## 
    ##   Level 7:   55 nodes to be scored   (380 eliminated genes)

    ## 
    ##   Level 6:   78 nodes to be scored   (613 eliminated genes)

    ## 
    ##   Level 5:   84 nodes to be scored   (917 eliminated genes)

    ## 
    ##   Level 4:   51 nodes to be scored   (1030 eliminated genes)

    ## 
    ##   Level 3:   32 nodes to be scored   (1221 eliminated genes)

    ## 
    ##   Level 2:   11 nodes to be scored   (1285 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1379 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 302 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   13 nodes to be scored   (10 eliminated genes)

    ## 
    ##   Level 8:   21 nodes to be scored   (145 eliminated genes)

    ## 
    ##   Level 7:   45 nodes to be scored   (464 eliminated genes)

    ## 
    ##   Level 6:   65 nodes to be scored   (596 eliminated genes)

    ## 
    ##   Level 5:   63 nodes to be scored   (925 eliminated genes)

    ## 
    ##   Level 4:   52 nodes to be scored   (1402 eliminated genes)

    ## 
    ##   Level 3:   25 nodes to be scored   (2053 eliminated genes)

    ## 
    ##   Level 2:   11 nodes to be scored   (2322 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2606 eliminated genes)

    ##        GO.ID                                          Term Annotated
    ## 1 GO:0000050                                    urea cycle         2
    ## 2 GO:0006487                protein N-linked glycosylation        19
    ## 3 GO:0000774    adenyl-nucleotide exchange factor activity         2
    ## 4 GO:0005302 L-tyrosine transmembrane transporter activity         9
    ## 5 GO:0004560                   alpha-L-fucosidase activity         6
    ##   Significant Expected Fisher               type
    ## 1           2     0.12 0.0036 Biological.Process
    ## 2           4     1.14 0.0237 Biological.Process
    ## 3           2     0.12 0.0037 Molecular.Function
    ## 4           3     0.55 0.0144 Molecular.Function
    ## 5           2     0.37 0.0475 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 762 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  10 nodes to be scored   (1 eliminated genes)

    ## 
    ##   Level 12:  14 nodes to be scored   (16 eliminated genes)

    ## 
    ##   Level 11:  25 nodes to be scored   (172 eliminated genes)

    ## 
    ##   Level 10:  40 nodes to be scored   (231 eliminated genes)

    ## 
    ##   Level 9:   71 nodes to be scored   (285 eliminated genes)

    ## 
    ##   Level 8:   83 nodes to be scored   (423 eliminated genes)

    ## 
    ##   Level 7:   104 nodes to be scored  (576 eliminated genes)

    ## 
    ##   Level 6:   132 nodes to be scored  (791 eliminated genes)

    ## 
    ##   Level 5:   133 nodes to be scored  (1077 eliminated genes)

    ## 
    ##   Level 4:   80 nodes to be scored   (1162 eliminated genes)

    ## 
    ##   Level 3:   52 nodes to be scored   (1319 eliminated genes)

    ## 
    ##   Level 2:   13 nodes to be scored   (1357 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1409 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 467 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  9 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 9:   28 nodes to be scored   (15 eliminated genes)

    ## 
    ##   Level 8:   37 nodes to be scored   (164 eliminated genes)

    ## 
    ##   Level 7:   76 nodes to be scored   (521 eliminated genes)

    ## 
    ##   Level 6:   105 nodes to be scored  (660 eliminated genes)

    ## 
    ##   Level 5:   80 nodes to be scored   (1086 eliminated genes)

    ## 
    ##   Level 4:   76 nodes to be scored   (1584 eliminated genes)

    ## 
    ##   Level 3:   38 nodes to be scored   (2183 eliminated genes)

    ## 
    ##   Level 2:   13 nodes to be scored   (2476 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2681 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0000398                         mRNA splicing, via spliceosome        63
    ## 2  GO:0000165                                           MAPK cascade        35
    ## 3  GO:0002218                   activation of innate immune response        15
    ## 4  GO:0002674     negative regulation of acute inflammatory response         4
    ## 5  GO:0003341                                        cilium movement        20
    ## 6  GO:0005524                                            ATP binding       263
    ## 7  GO:0004364                       glutathione transferase activity         4
    ## 8  GO:0000030                           mannosyltransferase activity        12
    ## 9  GO:0004675 transmembrane receptor protein serine/threonine kin...         3
    ## 10 GO:0003723                                            RNA binding       128
    ## 11 GO:0005198                           structural molecule activity        54
    ## 12 GO:0000146                           microfilament motor activity         9
    ## 13 GO:0005085             guanyl-nucleotide exchange factor activity        12
    ## 14 GO:0005229 intracellular calcium activated chloride channel ac...         4
    ## 15 GO:0000340                      RNA 7-methylguanosine cap binding         4
    ##    Significant Expected  Fisher               type
    ## 1           25    14.24 0.00170 Biological.Process
    ## 2           16     7.91 0.00190 Biological.Process
    ## 3            7     3.39 0.03350 Biological.Process
    ## 4            3     0.90 0.03820 Biological.Process
    ## 5            9     4.52 0.04420 Biological.Process
    ## 6           87    62.81 0.00024 Molecular.Function
    ## 7            4     0.96 0.00323 Molecular.Function
    ## 8            9     2.87 0.01020 Molecular.Function
    ## 9            3     0.72 0.01358 Molecular.Function
    ## 10          39    30.57 0.02325 Molecular.Function
    ## 11          21    12.90 0.03900 Molecular.Function
    ## 12           5     2.15 0.04048 Molecular.Function
    ## 13           6     2.87 0.04399 Molecular.Function
    ## 14           3     0.96 0.04461 Molecular.Function
    ## 15           3     0.96 0.04461 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 385 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  7 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  13 nodes to be scored   (156 eliminated genes)

    ## 
    ##   Level 10:  18 nodes to be scored   (209 eliminated genes)

    ## 
    ##   Level 9:   26 nodes to be scored   (239 eliminated genes)

    ## 
    ##   Level 8:   31 nodes to be scored   (333 eliminated genes)

    ## 
    ##   Level 7:   47 nodes to be scored   (408 eliminated genes)

    ## 
    ##   Level 6:   67 nodes to be scored   (597 eliminated genes)

    ## 
    ##   Level 5:   80 nodes to be scored   (903 eliminated genes)

    ## 
    ##   Level 4:   48 nodes to be scored   (1005 eliminated genes)

    ## 
    ##   Level 3:   33 nodes to be scored   (1173 eliminated genes)

    ## 
    ##   Level 2:   11 nodes to be scored   (1285 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1377 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 255 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  5 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 9:   11 nodes to be scored   (3 eliminated genes)

    ## 
    ##   Level 8:   16 nodes to be scored   (155 eliminated genes)

    ## 
    ##   Level 7:   32 nodes to be scored   (478 eliminated genes)

    ## 
    ##   Level 6:   50 nodes to be scored   (597 eliminated genes)

    ## 
    ##   Level 5:   50 nodes to be scored   (940 eliminated genes)

    ## 
    ##   Level 4:   47 nodes to be scored   (1411 eliminated genes)

    ## 
    ##   Level 3:   30 nodes to be scored   (2050 eliminated genes)

    ## 
    ##   Level 2:   11 nodes to be scored   (2354 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2642 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0000165                                           MAPK cascade        35
    ## 2  GO:0003341                                        cilium movement        20
    ## 3  GO:0019700                  organic phosphonate catabolic process        12
    ## 4  GO:0000462 maturation of SSU-rRNA from tricistronic rRNA trans...         5
    ## 5  GO:0003682                                      chromatin binding        13
    ## 6  GO:0001409 guanine nucleotide transmembrane transporter activi...        30
    ## 7  GO:0001965                        G-protein alpha-subunit binding        16
    ## 8  GO:0001786                             phosphatidylserine binding         4
    ## 9  GO:0003964                   RNA-directed DNA polymerase activity       122
    ## 10 GO:0004252                     serine-type endopeptidase activity        53
    ## 11 GO:0005524                                            ATP binding       263
    ## 12 GO:0003756                   protein disulfide isomerase activity         5
    ##    Significant Expected Fisher               type
    ## 1            7     2.01 0.0029 Biological.Process
    ## 2            5     1.15 0.0043 Biological.Process
    ## 3            3     0.69 0.0275 Biological.Process
    ## 4            2     0.29 0.0291 Biological.Process
    ## 5            4     0.95 0.0120 Molecular.Function
    ## 6            6     2.19 0.0190 Molecular.Function
    ## 7            4     1.17 0.0250 Molecular.Function
    ## 8            2     0.29 0.0290 Molecular.Function
    ## 9           15     8.92 0.0300 Molecular.Function
    ## 10           8     3.88 0.0360 Molecular.Function
    ## 11          27    19.23 0.0400 Molecular.Function
    ## 12           2     0.37 0.0460 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 470 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  4 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 12:  9 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 11:  16 nodes to be scored   (149 eliminated genes)

    ## 
    ##   Level 10:  22 nodes to be scored   (213 eliminated genes)

    ## 
    ##   Level 9:   36 nodes to be scored   (232 eliminated genes)

    ## 
    ##   Level 8:   43 nodes to be scored   (348 eliminated genes)

    ## 
    ##   Level 7:   60 nodes to be scored   (442 eliminated genes)

    ## 
    ##   Level 6:   77 nodes to be scored   (625 eliminated genes)

    ## 
    ##   Level 5:   88 nodes to be scored   (925 eliminated genes)

    ## 
    ##   Level 4:   61 nodes to be scored   (1010 eliminated genes)

    ## 
    ##   Level 3:   41 nodes to be scored   (1241 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (1358 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1402 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 231 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   10 nodes to be scored   (9 eliminated genes)

    ## 
    ##   Level 8:   14 nodes to be scored   (143 eliminated genes)

    ## 
    ##   Level 7:   28 nodes to be scored   (454 eliminated genes)

    ## 
    ##   Level 6:   48 nodes to be scored   (554 eliminated genes)

    ## 
    ##   Level 5:   48 nodes to be scored   (894 eliminated genes)

    ## 
    ##   Level 4:   39 nodes to be scored   (1296 eliminated genes)

    ## 
    ##   Level 3:   28 nodes to be scored   (2008 eliminated genes)

    ## 
    ##   Level 2:   11 nodes to be scored   (2260 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2659 eliminated genes)

    ##         GO.ID                                         Term Annotated
    ## 1  GO:0000165                                 MAPK cascade        35
    ## 2  GO:0008218                              bioluminescence        10
    ## 3  GO:0003341                              cilium movement        20
    ## 4  GO:0000038 very long-chain fatty acid metabolic process         9
    ## 5  GO:0006886              intracellular protein transport         5
    ## 6  GO:0005524                                  ATP binding       263
    ## 7  GO:0003964         RNA-directed DNA polymerase activity       122
    ## 8  GO:0004342   glucosamine-6-phosphate deaminase activity         2
    ## 9  GO:0005198                 structural molecule activity        54
    ## 10 GO:0000026       alpha-1,2-mannosyltransferase activity         3
    ## 11 GO:0001786                   phosphatidylserine binding         4
    ## 12 GO:0005201  extracellular matrix structural constituent        17
    ## 13 GO:0002020                             protease binding         5
    ##    Significant Expected   Fisher               type
    ## 1           10     2.68 0.000160 Biological.Process
    ## 2            5     0.77 0.000440 Biological.Process
    ## 3            5     1.53 0.014660 Biological.Process
    ## 4            3     0.69 0.026070 Biological.Process
    ## 5            2     0.38 0.049820 Biological.Process
    ## 6           36    18.95 0.000071 Molecular.Function
    ## 7           18     8.79 0.002300 Molecular.Function
    ## 8            2     0.14 0.005200 Molecular.Function
    ## 9           11     3.89 0.005400 Molecular.Function
    ## 10           2     0.22 0.014800 Molecular.Function
    ## 11           2     0.29 0.028100 Molecular.Function
    ## 12           4     1.22 0.029600 Molecular.Function
    ## 13           2     0.36 0.044700 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 318 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  4 nodes to be scored    (5 eliminated genes)

    ## 
    ##   Level 12:  7 nodes to be scored    (9 eliminated genes)

    ## 
    ##   Level 11:  11 nodes to be scored   (153 eliminated genes)

    ## 
    ##   Level 10:  14 nodes to be scored   (213 eliminated genes)

    ## 
    ##   Level 9:   26 nodes to be scored   (234 eliminated genes)

    ## 
    ##   Level 8:   28 nodes to be scored   (294 eliminated genes)

    ## 
    ##   Level 7:   34 nodes to be scored   (359 eliminated genes)

    ## 
    ##   Level 6:   48 nodes to be scored   (473 eliminated genes)

    ## 
    ##   Level 5:   61 nodes to be scored   (684 eliminated genes)

    ## 
    ##   Level 4:   40 nodes to be scored   (833 eliminated genes)

    ## 
    ##   Level 3:   29 nodes to be scored   (1068 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (1252 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1355 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 140 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   5 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   7 nodes to be scored    (142 eliminated genes)

    ## 
    ##   Level 7:   17 nodes to be scored   (438 eliminated genes)

    ## 
    ##   Level 6:   29 nodes to be scored   (482 eliminated genes)

    ## 
    ##   Level 5:   31 nodes to be scored   (832 eliminated genes)

    ## 
    ##   Level 4:   24 nodes to be scored   (1111 eliminated genes)

    ## 
    ##   Level 3:   17 nodes to be scored   (1744 eliminated genes)

    ## 
    ##   Level 2:   7 nodes to be scored    (1986 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2521 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0002064                            epithelial cell development         4
    ## 2  GO:0001945                               lymph vessel development         5
    ## 3  GO:0006836                             neurotransmitter transport         6
    ## 4  GO:0002314                 germinal center B cell differentiation         1
    ## 5  GO:0004252                     serine-type endopeptidase activity        53
    ## 6  GO:0001640 adenylate cyclase inhibiting G protein-coupled glut...         9
    ## 7  GO:0046982                    protein heterodimerization activity         1
    ## 8  GO:0010181                                            FMN binding         1
    ## 9  GO:0003975 UDP-N-acetylglucosamine-dolichyl-phosphate N-acetyl...         1
    ## 10 GO:0000978 RNA polymerase II cis-regulatory region sequence-sp...       100
    ## 11 GO:0004658                     propionyl-CoA carboxylase activity         2
    ## 12 GO:0004735             pyrroline-5-carboxylate reductase activity         2
    ##    Significant Expected Fisher               type
    ## 1            2     0.12 0.0050 Biological.Process
    ## 2            2     0.15 0.0082 Biological.Process
    ## 3            2     0.18 0.0120 Biological.Process
    ## 4            1     0.03 0.0298 Biological.Process
    ## 5            6     1.20 0.0011 Molecular.Function
    ## 6            2     0.20 0.0165 Molecular.Function
    ## 7            1     0.02 0.0227 Molecular.Function
    ## 8            1     0.02 0.0227 Molecular.Function
    ## 9            1     0.02 0.0227 Molecular.Function
    ## 10           6     2.27 0.0243 Molecular.Function
    ## 11           1     0.05 0.0449 Molecular.Function
    ## 12           1     0.05 0.0449 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 2 1 1 2 1 2 1 1 1 2 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 290 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  9 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 10:  13 nodes to be scored   (204 eliminated genes)

    ## 
    ##   Level 9:   24 nodes to be scored   (237 eliminated genes)

    ## 
    ##   Level 8:   21 nodes to be scored   (296 eliminated genes)

    ## 
    ##   Level 7:   29 nodes to be scored   (352 eliminated genes)

    ## 
    ##   Level 6:   52 nodes to be scored   (454 eliminated genes)

    ## 
    ##   Level 5:   54 nodes to be scored   (669 eliminated genes)

    ## 
    ##   Level 4:   42 nodes to be scored   (857 eliminated genes)

    ## 
    ##   Level 3:   29 nodes to be scored   (1078 eliminated genes)

    ## 
    ##   Level 2:   11 nodes to be scored   (1207 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1385 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 185 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  3 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 9:   10 nodes to be scored   (3 eliminated genes)

    ## 
    ##   Level 8:   10 nodes to be scored   (110 eliminated genes)

    ## 
    ##   Level 7:   19 nodes to be scored   (462 eliminated genes)

    ## 
    ##   Level 6:   32 nodes to be scored   (547 eliminated genes)

    ## 
    ##   Level 5:   39 nodes to be scored   (814 eliminated genes)

    ## 
    ##   Level 4:   36 nodes to be scored   (1132 eliminated genes)

    ## 
    ##   Level 3:   25 nodes to be scored   (1984 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (2264 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2611 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0002933                                    lipid hydroxylation         3
    ## 2  GO:0006368 transcription elongation by RNA polymerase II promo...         1
    ## 3  GO:0000335     negative regulation of transposition, DNA-mediated         1
    ## 4  GO:0015969             guanosine tetraphosphate metabolic process         9
    ## 5  GO:0004972                       NMDA glutamate receptor activity         3
    ## 6  GO:0005198                           structural molecule activity        54
    ## 7  GO:0005245                 voltage-gated calcium channel activity         7
    ## 8  GO:0005509                                    calcium ion binding        66
    ## 9  GO:0015101      organic cation transmembrane transporter activity         1
    ## 10 GO:0003909                                    DNA ligase activity         1
    ## 11 GO:0004441       inositol-1,4-bisphosphate 1-phosphatase activity         1
    ## 12 GO:0000400                          four-way junction DNA binding         1
    ## 13 GO:0005542                                     folic acid binding        10
    ##    Significant Expected Fisher               type
    ## 1            2     0.11 0.0038 Biological.Process
    ## 2            1     0.04 0.0361 Biological.Process
    ## 3            1     0.04 0.0361 Biological.Process
    ## 4            2     0.33 0.0392 Biological.Process
    ## 5            2     0.11 0.0037 Molecular.Function
    ## 6            6     1.93 0.0090 Molecular.Function
    ## 7            2     0.25 0.0235 Molecular.Function
    ## 8            6     2.35 0.0285 Molecular.Function
    ## 9            1     0.04 0.0357 Molecular.Function
    ## 10           1     0.04 0.0357 Molecular.Function
    ## 11           1     0.04 0.0357 Molecular.Function
    ## 12           1     0.04 0.0357 Molecular.Function
    ## 13           2     0.36 0.0470 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 503 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  6 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  9 nodes to be scored    (15 eliminated genes)

    ## 
    ##   Level 11:  14 nodes to be scored   (165 eliminated genes)

    ## 
    ##   Level 10:  20 nodes to be scored   (224 eliminated genes)

    ## 
    ##   Level 9:   36 nodes to be scored   (252 eliminated genes)

    ## 
    ##   Level 8:   48 nodes to be scored   (357 eliminated genes)

    ## 
    ##   Level 7:   65 nodes to be scored   (456 eliminated genes)

    ## 
    ##   Level 6:   91 nodes to be scored   (612 eliminated genes)

    ## 
    ##   Level 5:   93 nodes to be scored   (873 eliminated genes)

    ## 
    ##   Level 4:   62 nodes to be scored   (974 eliminated genes)

    ## 
    ##   Level 3:   44 nodes to be scored   (1280 eliminated genes)

    ## 
    ##   Level 2:   12 nodes to be scored   (1339 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1406 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 294 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   15 nodes to be scored   (9 eliminated genes)

    ## 
    ##   Level 8:   20 nodes to be scored   (147 eliminated genes)

    ## 
    ##   Level 7:   56 nodes to be scored   (494 eliminated genes)

    ## 
    ##   Level 6:   65 nodes to be scored   (605 eliminated genes)

    ## 
    ##   Level 5:   54 nodes to be scored   (1010 eliminated genes)

    ## 
    ##   Level 4:   43 nodes to be scored   (1430 eliminated genes)

    ## 
    ##   Level 3:   26 nodes to be scored   (2022 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (2291 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2587 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0002674     negative regulation of acute inflammatory response         4
    ## 2 GO:0001818             negative regulation of cytokine production        18
    ## 3 GO:0001732 formation of cytoplasmic translation initiation com...         6
    ## 4 GO:0005524                                            ATP binding       263
    ## 5 GO:0003964                   RNA-directed DNA polymerase activity       122
    ## 6 GO:0004177                                aminopeptidase activity         7
    ## 7 GO:0003682                                      chromatin binding        13
    ## 8 GO:0004623                              phospholipase A2 activity         4
    ## 9 GO:0000340                      RNA 7-methylguanosine cap binding         4
    ##   Significant Expected Fisher               type
    ## 1           3     0.42 0.0042 Biological.Process
    ## 2           6     1.89 0.0077 Biological.Process
    ## 3           3     0.63 0.0178 Biological.Process
    ## 4          39    24.82 0.0020 Molecular.Function
    ## 5          21    11.51 0.0041 Molecular.Function
    ## 6           3     0.66 0.0218 Molecular.Function
    ## 7           4     1.23 0.0279 Molecular.Function
    ## 8           2     0.38 0.0468 Molecular.Function
    ## 9           2     0.38 0.0468 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 602 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  6 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 12:  13 nodes to be scored   (3 eliminated genes)

    ## 
    ##   Level 11:  22 nodes to be scored   (147 eliminated genes)

    ## 
    ##   Level 10:  34 nodes to be scored   (229 eliminated genes)

    ## 
    ##   Level 9:   52 nodes to be scored   (258 eliminated genes)

    ## 
    ##   Level 8:   62 nodes to be scored   (368 eliminated genes)

    ## 
    ##   Level 7:   80 nodes to be scored   (492 eliminated genes)

    ## 
    ##   Level 6:   103 nodes to be scored  (692 eliminated genes)

    ## 
    ##   Level 5:   102 nodes to be scored  (986 eliminated genes)

    ## 
    ##   Level 4:   65 nodes to be scored   (1115 eliminated genes)

    ## 
    ##   Level 3:   47 nodes to be scored   (1276 eliminated genes)

    ## 
    ##   Level 2:   12 nodes to be scored   (1340 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1406 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 254 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 9:   12 nodes to be scored   (13 eliminated genes)

    ## 
    ##   Level 8:   19 nodes to be scored   (146 eliminated genes)

    ## 
    ##   Level 7:   35 nodes to be scored   (468 eliminated genes)

    ## 
    ##   Level 6:   52 nodes to be scored   (541 eliminated genes)

    ## 
    ##   Level 5:   49 nodes to be scored   (907 eliminated genes)

    ## 
    ##   Level 4:   41 nodes to be scored   (1323 eliminated genes)

    ## 
    ##   Level 3:   26 nodes to be scored   (1973 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (2269 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2622 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0000184 nuclear-transcribed mRNA catabolic process, nonsens...        10
    ## 2  GO:0002064                            epithelial cell development         4
    ## 3  GO:0001782                                     B cell homeostasis         7
    ## 4  GO:0003341                                        cilium movement        20
    ## 5  GO:0002244          hematopoietic progenitor cell differentiation         4
    ## 6  GO:0000723                                   telomere maintenance         4
    ## 7  GO:0003676                                   nucleic acid binding       497
    ## 8  GO:0001409 guanine nucleotide transmembrane transporter activi...        30
    ## 9  GO:0003964                   RNA-directed DNA polymerase activity       122
    ## 10 GO:0005509                                    calcium ion binding        66
    ## 11 GO:0003824                                     catalytic activity      1042
    ## 12 GO:0005085             guanyl-nucleotide exchange factor activity        12
    ## 13 GO:0004972                       NMDA glutamate receptor activity         3
    ## 14 GO:0005507                                     copper ion binding        18
    ## 15 GO:0003723                                            RNA binding       128
    ##    Significant Expected  Fisher               type
    ## 1            5     0.96 0.00130 Biological.Process
    ## 2            3     0.39 0.00330 Biological.Process
    ## 3            3     0.67 0.02290 Biological.Process
    ## 4            6     1.93 0.02890 Biological.Process
    ## 5            2     0.39 0.04860 Biological.Process
    ## 6            2     0.39 0.04860 Biological.Process
    ## 7           64    50.85 0.00057 Molecular.Function
    ## 8            8     3.07 0.00853 Molecular.Function
    ## 9           21    12.48 0.01051 Molecular.Function
    ## 10          13     6.75 0.01418 Molecular.Function
    ## 11          86   106.60 0.02051 Molecular.Function
    ## 12           4     1.23 0.02735 Molecular.Function
    ## 13           2     0.31 0.02918 Molecular.Function
    ## 14           5     1.84 0.03031 Molecular.Function
    ## 15          20    13.10 0.03108 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 712 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  5 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  8 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 12:  12 nodes to be scored   (22 eliminated genes)

    ## 
    ##   Level 11:  26 nodes to be scored   (173 eliminated genes)

    ## 
    ##   Level 10:  41 nodes to be scored   (234 eliminated genes)

    ## 
    ##   Level 9:   71 nodes to be scored   (294 eliminated genes)

    ## 
    ##   Level 8:   80 nodes to be scored   (400 eliminated genes)

    ## 
    ##   Level 7:   99 nodes to be scored   (513 eliminated genes)

    ## 
    ##   Level 6:   118 nodes to be scored  (721 eliminated genes)

    ## 
    ##   Level 5:   117 nodes to be scored  (996 eliminated genes)

    ## 
    ##   Level 4:   72 nodes to be scored   (1127 eliminated genes)

    ## 
    ##   Level 3:   48 nodes to be scored   (1317 eliminated genes)

    ## 
    ##   Level 2:   13 nodes to be scored   (1367 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1408 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 444 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  5 nodes to be scored    (2 eliminated genes)

    ## 
    ##   Level 10:  11 nodes to be scored   (7 eliminated genes)

    ## 
    ##   Level 9:   24 nodes to be scored   (18 eliminated genes)

    ## 
    ##   Level 8:   41 nodes to be scored   (166 eliminated genes)

    ## 
    ##   Level 7:   69 nodes to be scored   (509 eliminated genes)

    ## 
    ##   Level 6:   93 nodes to be scored   (675 eliminated genes)

    ## 
    ##   Level 5:   77 nodes to be scored   (1044 eliminated genes)

    ## 
    ##   Level 4:   69 nodes to be scored   (1501 eliminated genes)

    ## 
    ##   Level 3:   36 nodes to be scored   (2164 eliminated genes)

    ## 
    ##   Level 2:   13 nodes to be scored   (2486 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2679 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0003341                                        cilium movement        20
    ## 2  GO:0001933         negative regulation of protein phosphorylation         5
    ## 3  GO:0006886                        intracellular protein transport         5
    ## 4  GO:0001934         positive regulation of protein phosphorylation         4
    ## 5  GO:0002244          hematopoietic progenitor cell differentiation         4
    ## 6  GO:0002042      cell migration involved in sprouting angiogenesis         7
    ## 7  GO:0007224                           smoothened signaling pathway         2
    ## 8  GO:0001895                                     retina homeostasis         2
    ## 9  GO:0005524                                            ATP binding       263
    ## 10 GO:0004715 non-membrane spanning protein tyrosine kinase activ...         6
    ## 11 GO:0005388                    P-type calcium transporter activity         4
    ## 12 GO:0003964                   RNA-directed DNA polymerase activity       122
    ## 13 GO:0003676                                   nucleic acid binding       497
    ## 14 GO:0004518                                      nuclease activity        30
    ## 15 GO:0001965                        G-protein alpha-subunit binding        16
    ## 16 GO:0001786                             phosphatidylserine binding         4
    ## 17 GO:0005198                           structural molecule activity        54
    ## 18 GO:0003987                            acetate-CoA ligase activity         2
    ## 19 GO:0004100                               chitin synthase activity         2
    ## 20 GO:0004579 dolichyl-diphosphooligosaccharide-protein glycotran...         2
    ##    Significant Expected   Fisher               type
    ## 1           12     3.61 0.000031 Biological.Process
    ## 2            4     0.90 0.004500 Biological.Process
    ## 3            4     0.90 0.005800 Biological.Process
    ## 4            3     0.72 0.020200 Biological.Process
    ## 5            3     0.72 0.020200 Biological.Process
    ## 6            4     1.27 0.023100 Biological.Process
    ## 7            2     0.36 0.032600 Biological.Process
    ## 8            2     0.36 0.032600 Biological.Process
    ## 9           72    51.07 0.000610 Molecular.Function
    ## 10           5     1.16 0.001370 Molecular.Function
    ## 11           4     0.78 0.001410 Molecular.Function
    ## 12          37    23.69 0.002110 Molecular.Function
    ## 13         104    96.50 0.005160 Molecular.Function
    ## 14           9     5.82 0.015280 Molecular.Function
    ## 15           7     3.11 0.022520 Molecular.Function
    ## 16           3     0.78 0.024920 Molecular.Function
    ## 17          14    10.48 0.030040 Molecular.Function
    ## 18           2     0.39 0.037640 Molecular.Function
    ## 19           2     0.39 0.037640 Molecular.Function
    ## 20           2     0.39 0.037640 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 370 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  6 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 12:  8 nodes to be scored    (17 eliminated genes)

    ## 
    ##   Level 11:  16 nodes to be scored   (165 eliminated genes)

    ## 
    ##   Level 10:  21 nodes to be scored   (212 eliminated genes)

    ## 
    ##   Level 9:   37 nodes to be scored   (244 eliminated genes)

    ## 
    ##   Level 8:   32 nodes to be scored   (314 eliminated genes)

    ## 
    ##   Level 7:   42 nodes to be scored   (373 eliminated genes)

    ## 
    ##   Level 6:   58 nodes to be scored   (483 eliminated genes)

    ## 
    ##   Level 5:   67 nodes to be scored   (732 eliminated genes)

    ## 
    ##   Level 4:   42 nodes to be scored   (878 eliminated genes)

    ## 
    ##   Level 3:   26 nodes to be scored   (1042 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (1182 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1311 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 216 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  3 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 9:   10 nodes to be scored   (1 eliminated genes)

    ## 
    ##   Level 8:   14 nodes to be scored   (143 eliminated genes)

    ## 
    ##   Level 7:   31 nodes to be scored   (473 eliminated genes)

    ## 
    ##   Level 6:   43 nodes to be scored   (598 eliminated genes)

    ## 
    ##   Level 5:   47 nodes to be scored   (951 eliminated genes)

    ## 
    ##   Level 4:   36 nodes to be scored   (1298 eliminated genes)

    ## 
    ##   Level 3:   20 nodes to be scored   (1937 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (2179 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2564 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0002064                            epithelial cell development         4
    ## 2  GO:0000278                                     mitotic cell cycle       160
    ## 3  GO:0003179                              heart valve morphogenesis         4
    ## 4  GO:0001701                         in utero embryonic development        15
    ## 5  GO:0002218                   activation of innate immune response        15
    ## 6  GO:0001732 formation of cytoplasmic translation initiation com...         6
    ## 7  GO:0000712      resolution of meiotic recombination intermediates         1
    ## 8  GO:0001771                        immunological synapse formation         1
    ## 9  GO:0000737                 DNA catabolic process, endonucleolytic         1
    ## 10 GO:0002314                 germinal center B cell differentiation         1
    ## 11 GO:0002225 positive regulation of antimicrobial peptide produc...         1
    ## 12 GO:0001573                          ganglioside metabolic process         1
    ## 13 GO:0000472 endonucleolytic cleavage to generate mature 5'-end ...         1
    ## 14 GO:0006513                             protein monoubiquitination         1
    ## 15 GO:0000987    cis-regulatory region sequence-specific DNA binding       112
    ## 16 GO:0004325                                ferrochelatase activity         1
    ## 17 GO:0004458          D-lactate dehydrogenase (cytochrome) activity         1
    ## 18 GO:0003989                        acetyl-CoA carboxylase activity         1
    ## 19 GO:0004729 oxygen-dependent protoporphyrinogen oxidase activit...         1
    ## 20 GO:0004332                fructose-bisphosphate aldolase activity         1
    ## 21 GO:0008957              phenylacetaldehyde dehydrogenase activity         1
    ## 22 GO:0010181                                            FMN binding         1
    ## 23 GO:0004818                         glutamate-tRNA ligase activity         1
    ## 24 GO:0004368 glycerol-3-phosphate dehydrogenase (quinone) activi...         1
    ## 25 GO:0005220 inositol 1,4,5-trisphosphate-sensitive calcium-rele...         1
    ## 26 GO:0004408                     holocytochrome-c synthase activity         1
    ## 27 GO:0004802                                 transketolase activity         1
    ## 28 GO:0050218                         propionate-CoA ligase activity         1
    ## 29 GO:0003960                       NADPH:quinone reductase activity         1
    ##    Significant Expected  Fisher               type
    ## 1            3     0.17 0.00030 Biological.Process
    ## 2            9     6.92 0.00330 Biological.Process
    ## 3            2     0.17 0.01040 Biological.Process
    ## 4            3     0.65 0.02410 Biological.Process
    ## 5            3     0.65 0.02410 Biological.Process
    ## 6            2     0.26 0.02470 Biological.Process
    ## 7            1     0.04 0.04320 Biological.Process
    ## 8            1     0.04 0.04320 Biological.Process
    ## 9            1     0.04 0.04320 Biological.Process
    ## 10           1     0.04 0.04320 Biological.Process
    ## 11           1     0.04 0.04320 Biological.Process
    ## 12           1     0.04 0.04320 Biological.Process
    ## 13           1     0.04 0.04320 Biological.Process
    ## 14           1     0.04 0.04320 Biological.Process
    ## 15           7     4.20 0.00075 Molecular.Function
    ## 16           1     0.04 0.03746 Molecular.Function
    ## 17           1     0.04 0.03746 Molecular.Function
    ## 18           1     0.04 0.03746 Molecular.Function
    ## 19           1     0.04 0.03746 Molecular.Function
    ## 20           1     0.04 0.03746 Molecular.Function
    ## 21           1     0.04 0.03746 Molecular.Function
    ## 22           1     0.04 0.03746 Molecular.Function
    ## 23           1     0.04 0.03746 Molecular.Function
    ## 24           1     0.04 0.03746 Molecular.Function
    ## 25           1     0.04 0.03746 Molecular.Function
    ## 26           1     0.04 0.03746 Molecular.Function
    ## 27           1     0.04 0.03746 Molecular.Function
    ## 28           1     0.04 0.03746 Molecular.Function
    ## 29           1     0.04 0.03746 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 373 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  5 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  8 nodes to be scored    (15 eliminated genes)

    ## 
    ##   Level 11:  13 nodes to be scored   (164 eliminated genes)

    ## 
    ##   Level 10:  19 nodes to be scored   (216 eliminated genes)

    ## 
    ##   Level 9:   27 nodes to be scored   (235 eliminated genes)

    ## 
    ##   Level 8:   35 nodes to be scored   (320 eliminated genes)

    ## 
    ##   Level 7:   48 nodes to be scored   (386 eliminated genes)

    ## 
    ##   Level 6:   60 nodes to be scored   (502 eliminated genes)

    ## 
    ##   Level 5:   71 nodes to be scored   (780 eliminated genes)

    ## 
    ##   Level 4:   44 nodes to be scored   (893 eliminated genes)

    ## 
    ##   Level 3:   29 nodes to be scored   (1133 eliminated genes)

    ## 
    ##   Level 2:   11 nodes to be scored   (1266 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1360 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 194 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  1 nodes to be scored    (2 eliminated genes)

    ## 
    ##   Level 10:  5 nodes to be scored    (2 eliminated genes)

    ## 
    ##   Level 9:   9 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 8:   13 nodes to be scored   (157 eliminated genes)

    ## 
    ##   Level 7:   29 nodes to be scored   (471 eliminated genes)

    ## 
    ##   Level 6:   37 nodes to be scored   (570 eliminated genes)

    ## 
    ##   Level 5:   39 nodes to be scored   (900 eliminated genes)

    ## 
    ##   Level 4:   28 nodes to be scored   (1235 eliminated genes)

    ## 
    ##   Level 3:   20 nodes to be scored   (1946 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (2165 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2432 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0000045                                 autophagosome assembly         7
    ## 2  GO:0001558                              regulation of cell growth         3
    ## 3  GO:0001960 negative regulation of cytokine-mediated signaling ...         1
    ## 4  GO:0001573                          ganglioside metabolic process         1
    ## 5  GO:0007175 negative regulation of epidermal growth factor-acti...         1
    ## 6  GO:0003964                   RNA-directed DNA polymerase activity       122
    ## 7  GO:0005044                            scavenger receptor activity        13
    ## 8  GO:0005509                                    calcium ion binding        66
    ## 9  GO:0005245                 voltage-gated calcium channel activity         7
    ## 10 GO:0005324             long-chain fatty acid transporter activity         1
    ## 11 GO:0004305                           ethanolamine kinase activity         1
    ## 12 GO:0005536                                        glucose binding         1
    ## 13 GO:0004492               methylmalonyl-CoA decarboxylase activity         1
    ## 14 GO:0004321                       fatty-acyl-CoA synthase activity         1
    ## 15 GO:0002134                                            UTP binding         1
    ##    Significant Expected   Fisher               type
    ## 1            3     0.28 0.002000 Biological.Process
    ## 2            2     0.12 0.004700 Biological.Process
    ## 3            1     0.04 0.040400 Biological.Process
    ## 4            1     0.04 0.040400 Biological.Process
    ## 5            1     0.04 0.040400 Biological.Process
    ## 6           17     5.76 0.000038 Molecular.Function
    ## 7            3     0.61 0.021000 Molecular.Function
    ## 8            7     3.11 0.034000 Molecular.Function
    ## 9            2     0.33 0.040000 Molecular.Function
    ## 10           1     0.05 0.047000 Molecular.Function
    ## 11           1     0.05 0.047000 Molecular.Function
    ## 12           1     0.05 0.047000 Molecular.Function
    ## 13           1     0.05 0.047000 Molecular.Function
    ## 14           1     0.05 0.047000 Molecular.Function
    ## 15           1     0.05 0.047000 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 488 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  8 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 11:  14 nodes to be scored   (149 eliminated genes)

    ## 
    ##   Level 10:  24 nodes to be scored   (212 eliminated genes)

    ## 
    ##   Level 9:   39 nodes to be scored   (232 eliminated genes)

    ## 
    ##   Level 8:   46 nodes to be scored   (358 eliminated genes)

    ## 
    ##   Level 7:   71 nodes to be scored   (478 eliminated genes)

    ## 
    ##   Level 6:   82 nodes to be scored   (612 eliminated genes)

    ## 
    ##   Level 5:   87 nodes to be scored   (929 eliminated genes)

    ## 
    ##   Level 4:   60 nodes to be scored   (991 eliminated genes)

    ## 
    ##   Level 3:   41 nodes to be scored   (1148 eliminated genes)

    ## 
    ##   Level 2:   11 nodes to be scored   (1293 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1405 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 284 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  3 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 9:   10 nodes to be scored   (4 eliminated genes)

    ## 
    ##   Level 8:   26 nodes to be scored   (105 eliminated genes)

    ## 
    ##   Level 7:   39 nodes to be scored   (433 eliminated genes)

    ## 
    ##   Level 6:   62 nodes to be scored   (591 eliminated genes)

    ## 
    ##   Level 5:   57 nodes to be scored   (978 eliminated genes)

    ## 
    ##   Level 4:   45 nodes to be scored   (1426 eliminated genes)

    ## 
    ##   Level 3:   27 nodes to be scored   (2057 eliminated genes)

    ## 
    ##   Level 2:   12 nodes to be scored   (2319 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2625 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0000027                       ribosomal large subunit assembly         2
    ## 2 GO:0000423                                              mitophagy         2
    ## 3 GO:0002931                                   response to ischemia         3
    ## 4 GO:0000281                                    mitotic cytokinesis        15
    ## 5 GO:0000122 negative regulation of transcription by RNA polymer...       140
    ## 6 GO:0000038           very long-chain fatty acid metabolic process         9
    ## 7 GO:0005044                            scavenger receptor activity        13
    ## 8 GO:0004675 transmembrane receptor protein serine/threonine kin...         3
    ## 9 GO:0003964                   RNA-directed DNA polymerase activity       122
    ##   Significant Expected Fisher               type
    ## 1           2     0.15 0.0053 Biological.Process
    ## 2           2     0.15 0.0053 Biological.Process
    ## 3           2     0.22 0.0151 Biological.Process
    ## 4           4     1.09 0.0195 Biological.Process
    ## 5          17    10.22 0.0206 Biological.Process
    ## 6           3     0.66 0.0230 Biological.Process
    ## 7           4     1.01 0.0140 Molecular.Function
    ## 8           2     0.23 0.0170 Molecular.Function
    ## 9          15     9.45 0.0470 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 490 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  7 nodes to be scored    (2 eliminated genes)

    ## 
    ##   Level 11:  16 nodes to be scored   (147 eliminated genes)

    ## 
    ##   Level 10:  24 nodes to be scored   (225 eliminated genes)

    ## 
    ##   Level 9:   45 nodes to be scored   (261 eliminated genes)

    ## 
    ##   Level 8:   50 nodes to be scored   (367 eliminated genes)

    ## 
    ##   Level 7:   71 nodes to be scored   (459 eliminated genes)

    ## 
    ##   Level 6:   89 nodes to be scored   (665 eliminated genes)

    ## 
    ##   Level 5:   85 nodes to be scored   (952 eliminated genes)

    ## 
    ##   Level 4:   50 nodes to be scored   (1027 eliminated genes)

    ## 
    ##   Level 3:   37 nodes to be scored   (1167 eliminated genes)

    ## 
    ##   Level 2:   11 nodes to be scored   (1275 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1392 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 270 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  6 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 9:   15 nodes to be scored   (12 eliminated genes)

    ## 
    ##   Level 8:   19 nodes to be scored   (156 eliminated genes)

    ## 
    ##   Level 7:   38 nodes to be scored   (481 eliminated genes)

    ## 
    ##   Level 6:   55 nodes to be scored   (583 eliminated genes)

    ## 
    ##   Level 5:   50 nodes to be scored   (909 eliminated genes)

    ## 
    ##   Level 4:   43 nodes to be scored   (1375 eliminated genes)

    ## 
    ##   Level 3:   28 nodes to be scored   (2053 eliminated genes)

    ## 
    ##   Level 2:   12 nodes to be scored   (2289 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2565 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0001822                                     kidney development        65
    ## 2 GO:0000165                                           MAPK cascade        35
    ## 3 GO:0001764                                       neuron migration         8
    ## 4 GO:0000184 nuclear-transcribed mRNA catabolic process, nonsens...        10
    ## 5 GO:0005524                                            ATP binding       263
    ## 6 GO:0001409 guanine nucleotide transmembrane transporter activi...        30
    ## 7 GO:0003964                   RNA-directed DNA polymerase activity       122
    ## 8 GO:0001222                      transcription corepressor binding         3
    ## 9 GO:0003777                             microtubule motor activity        16
    ##   Significant Expected   Fisher               type
    ## 1          18     5.44 0.000002 Biological.Process
    ## 2           8     2.93 0.006500 Biological.Process
    ## 3           3     0.67 0.023400 Biological.Process
    ## 4           3     0.84 0.044300 Biological.Process
    ## 5          38    20.37 0.000066 Molecular.Function
    ## 6           8     2.32 0.001500 Molecular.Function
    ## 7          17     9.45 0.011400 Molecular.Function
    ## 8           2     0.23 0.017000 Molecular.Function
    ## 9           4     1.24 0.030400 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 294 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  5 nodes to be scored    (11 eliminated genes)

    ## 
    ##   Level 11:  9 nodes to be scored    (159 eliminated genes)

    ## 
    ##   Level 10:  15 nodes to be scored   (215 eliminated genes)

    ## 
    ##   Level 9:   23 nodes to be scored   (242 eliminated genes)

    ## 
    ##   Level 8:   24 nodes to be scored   (303 eliminated genes)

    ## 
    ##   Level 7:   30 nodes to be scored   (332 eliminated genes)

    ## 
    ##   Level 6:   44 nodes to be scored   (451 eliminated genes)

    ## 
    ##   Level 5:   59 nodes to be scored   (665 eliminated genes)

    ## 
    ##   Level 4:   45 nodes to be scored   (830 eliminated genes)

    ## 
    ##   Level 3:   26 nodes to be scored   (1039 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (1196 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1325 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 203 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 9:   9 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 8:   13 nodes to be scored   (152 eliminated genes)

    ## 
    ##   Level 7:   22 nodes to be scored   (462 eliminated genes)

    ## 
    ##   Level 6:   37 nodes to be scored   (576 eliminated genes)

    ## 
    ##   Level 5:   45 nodes to be scored   (836 eliminated genes)

    ## 
    ##   Level 4:   38 nodes to be scored   (1291 eliminated genes)

    ## 
    ##   Level 3:   22 nodes to be scored   (1996 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (2290 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2589 eliminated genes)

    ##         GO.ID                                               Term Annotated
    ## 1  GO:0002221     pattern recognition receptor signaling pathway        11
    ## 2  GO:0002674 negative regulation of acute inflammatory response         4
    ## 3  GO:0002040                             sprouting angiogenesis        12
    ## 4  GO:0001172                        RNA-templated transcription         1
    ## 5  GO:0000076               DNA replication checkpoint signaling         1
    ## 6  GO:0000492                            box C/D snoRNP assembly         1
    ## 7  GO:0005524                                        ATP binding       263
    ## 8  GO:0005201        extracellular matrix structural constituent        17
    ## 9  GO:0004017                          adenylate kinase activity         1
    ## 10 GO:0004505             phenylalanine 4-monooxygenase activity         1
    ##    Significant Expected   Fisher               type
    ## 1            5     0.47 0.000045 Biological.Process
    ## 2            3     0.17 0.000280 Biological.Process
    ## 3            3     0.51 0.016000 Biological.Process
    ## 4            1     0.04 0.042520 Biological.Process
    ## 5            1     0.04 0.042520 Biological.Process
    ## 6            1     0.04 0.042520 Biological.Process
    ## 7           20    11.27 0.007000 Molecular.Function
    ## 8            3     0.73 0.034000 Molecular.Function
    ## 9            1     0.04 0.043000 Molecular.Function
    ## 10           1     0.04 0.043000 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 407 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  9 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 10:  13 nodes to be scored   (204 eliminated genes)

    ## 
    ##   Level 9:   29 nodes to be scored   (218 eliminated genes)

    ## 
    ##   Level 8:   39 nodes to be scored   (289 eliminated genes)

    ## 
    ##   Level 7:   53 nodes to be scored   (406 eliminated genes)

    ## 
    ##   Level 6:   78 nodes to be scored   (625 eliminated genes)

    ## 
    ##   Level 5:   76 nodes to be scored   (906 eliminated genes)

    ## 
    ##   Level 4:   53 nodes to be scored   (1025 eliminated genes)

    ## 
    ##   Level 3:   38 nodes to be scored   (1215 eliminated genes)

    ## 
    ##   Level 2:   13 nodes to be scored   (1310 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1387 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 274 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   13 nodes to be scored   (1 eliminated genes)

    ## 
    ##   Level 8:   19 nodes to be scored   (150 eliminated genes)

    ## 
    ##   Level 7:   34 nodes to be scored   (479 eliminated genes)

    ## 
    ##   Level 6:   53 nodes to be scored   (594 eliminated genes)

    ## 
    ##   Level 5:   57 nodes to be scored   (938 eliminated genes)

    ## 
    ##   Level 4:   53 nodes to be scored   (1347 eliminated genes)

    ## 
    ##   Level 3:   27 nodes to be scored   (2008 eliminated genes)

    ## 
    ##   Level 2:   12 nodes to be scored   (2384 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2635 eliminated genes)

    ##         GO.ID                                           Term Annotated
    ## 1  GO:0003341                                cilium movement        20
    ## 2  GO:0002790                              peptide secretion         2
    ## 3  GO:0001822                             kidney development        65
    ## 4  GO:0002221 pattern recognition receptor signaling pathway        11
    ## 5  GO:0001558                      regulation of cell growth         3
    ## 6  GO:0015969     guanosine tetraphosphate metabolic process         9
    ## 7  GO:0001707                             mesoderm formation         4
    ## 8  GO:0005524                                    ATP binding       263
    ## 9  GO:0003964           RNA-directed DNA polymerase activity       122
    ## 10 GO:0000340              RNA 7-methylguanosine cap binding         4
    ## 11 GO:0005385    zinc ion transmembrane transporter activity         4
    ##    Significant Expected   Fisher               type
    ## 1            8     1.46 0.000037 Biological.Process
    ## 2            2     0.15 0.005300 Biological.Process
    ## 3           11     4.74 0.005800 Biological.Process
    ## 4            4     0.80 0.005900 Biological.Process
    ## 5            2     0.22 0.015100 Biological.Process
    ## 6            3     0.66 0.023000 Biological.Process
    ## 7            2     0.29 0.028700 Biological.Process
    ## 8           35    19.61 0.000330 Molecular.Function
    ## 9           17     9.10 0.007870 Molecular.Function
    ## 10           2     0.30 0.030030 Molecular.Function
    ## 11           2     0.30 0.030030 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 201 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  4 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 12:  5 nodes to be scored    (8 eliminated genes)

    ## 
    ##   Level 11:  7 nodes to be scored    (167 eliminated genes)

    ## 
    ##   Level 10:  10 nodes to be scored   (207 eliminated genes)

    ## 
    ##   Level 9:   15 nodes to be scored   (209 eliminated genes)

    ## 
    ##   Level 8:   14 nodes to be scored   (289 eliminated genes)

    ## 
    ##   Level 7:   16 nodes to be scored   (305 eliminated genes)

    ## 
    ##   Level 6:   25 nodes to be scored   (422 eliminated genes)

    ## 
    ##   Level 5:   39 nodes to be scored   (503 eliminated genes)

    ## 
    ##   Level 4:   32 nodes to be scored   (746 eliminated genes)

    ## 
    ##   Level 3:   21 nodes to be scored   (1019 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (1116 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1298 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 171 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   10 nodes to be scored   (10 eliminated genes)

    ## 
    ##   Level 8:   12 nodes to be scored   (143 eliminated genes)

    ## 
    ##   Level 7:   20 nodes to be scored   (457 eliminated genes)

    ## 
    ##   Level 6:   30 nodes to be scored   (556 eliminated genes)

    ## 
    ##   Level 5:   33 nodes to be scored   (805 eliminated genes)

    ## 
    ##   Level 4:   29 nodes to be scored   (1091 eliminated genes)

    ## 
    ##   Level 3:   21 nodes to be scored   (1852 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (2184 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2541 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0003341                                        cilium movement        20
    ## 2  GO:0001676                long-chain fatty acid metabolic process         3
    ## 3  GO:0001736                       establishment of planar polarity         5
    ## 4  GO:0003943             N-acetylgalactosamine-4-sulfatase activity         8
    ## 5  GO:0004587               ornithine-oxo-acid transaminase activity         1
    ## 6  GO:0004760                  serine-pyruvate transaminase activity         1
    ## 7  GO:0004149 dihydrolipoyllysine-residue succinyltransferase act...         1
    ## 8  GO:0042802                              identical protein binding         1
    ## 9  GO:0005337          nucleoside transmembrane transporter activity         1
    ## 10 GO:0004054                               arginine kinase activity         1
    ## 11 GO:0004591 oxoglutarate dehydrogenase (succinyl-transferring) ...         1
    ## 12 GO:0000048                           peptidyltransferase activity         1
    ## 13 GO:0003980 UDP-glucose:glycoprotein glucosyltransferase activi...         1
    ## 14 GO:0001002 RNA polymerase III type 1 promoter sequence-specifi...         1
    ## 15 GO:0004842                 ubiquitin-protein transferase activity        35
    ## 16 GO:0004190                   aspartic-type endopeptidase activity         2
    ## 17 GO:0005272                                sodium channel activity         2
    ## 18 GO:0004714 transmembrane receptor protein tyrosine kinase acti...         2
    ##    Significant Expected  Fisher               type
    ## 1            4     0.40 0.00048 Biological.Process
    ## 2            2     0.06 0.00113 Biological.Process
    ## 3            2     0.10 0.00366 Biological.Process
    ## 4            2     0.17 0.01200 Molecular.Function
    ## 5            1     0.02 0.02200 Molecular.Function
    ## 6            1     0.02 0.02200 Molecular.Function
    ## 7            1     0.02 0.02200 Molecular.Function
    ## 8            1     0.02 0.02200 Molecular.Function
    ## 9            1     0.02 0.02200 Molecular.Function
    ## 10           1     0.02 0.02200 Molecular.Function
    ## 11           1     0.02 0.02200 Molecular.Function
    ## 12           1     0.02 0.02200 Molecular.Function
    ## 13           1     0.02 0.02200 Molecular.Function
    ## 14           1     0.02 0.02200 Molecular.Function
    ## 15           3     0.76 0.03800 Molecular.Function
    ## 16           1     0.04 0.04300 Molecular.Function
    ## 17           1     0.04 0.04300 Molecular.Function
    ## 18           1     0.04 0.04300 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 234 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  2 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 12:  5 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 11:  9 nodes to be scored    (144 eliminated genes)

    ## 
    ##   Level 10:  12 nodes to be scored   (146 eliminated genes)

    ## 
    ##   Level 9:   19 nodes to be scored   (149 eliminated genes)

    ## 
    ##   Level 8:   18 nodes to be scored   (160 eliminated genes)

    ## 
    ##   Level 7:   26 nodes to be scored   (207 eliminated genes)

    ## 
    ##   Level 6:   34 nodes to be scored   (357 eliminated genes)

    ## 
    ##   Level 5:   47 nodes to be scored   (473 eliminated genes)

    ## 
    ##   Level 4:   32 nodes to be scored   (720 eliminated genes)

    ## 
    ##   Level 3:   18 nodes to be scored   (1000 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (1204 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1312 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 159 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   5 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 8:   10 nodes to be scored   (11 eliminated genes)

    ## 
    ##   Level 7:   17 nodes to be scored   (341 eliminated genes)

    ## 
    ##   Level 6:   32 nodes to be scored   (432 eliminated genes)

    ## 
    ##   Level 5:   34 nodes to be scored   (697 eliminated genes)

    ## 
    ##   Level 4:   30 nodes to be scored   (1026 eliminated genes)

    ## 
    ##   Level 3:   19 nodes to be scored   (1806 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (2121 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2447 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0008218                                        bioluminescence        10
    ## 2  GO:0002025 norepinephrine-epinephrine-mediated vasodilation in...         1
    ## 3  GO:0001731         formation of translation preinitiation complex         1
    ## 4  GO:0001832                                      blastocyst growth         1
    ## 5  GO:0001516                     prostaglandin biosynthetic process         2
    ## 6  GO:0002121                         inter-male aggressive behavior         2
    ## 7  GO:0004089                         carbonate dehydratase activity         2
    ## 8  GO:0005302          L-tyrosine transmembrane transporter activity         9
    ## 9  GO:0004142      diacylglycerol cholinephosphotransferase activity         1
    ## 10 GO:0004067                                  asparaginase activity         1
    ## 11 GO:0003913                                DNA photolyase activity         1
    ## 12 GO:0002161                        aminoacyl-tRNA editing activity         1
    ## 13 GO:0004017                              adenylate kinase activity         1
    ## 14 GO:0004252                     serine-type endopeptidase activity        53
    ## 15 GO:0003978                       UDP-glucose 4-epimerase activity         2
    ## 16 GO:0004081 bis(5'-nucleosyl)-tetraphosphatase (asymmetrical) a...         2
    ## 17 GO:0004720                      protein-lysine 6-oxidase activity         2
    ## 18 GO:0004176                       ATP-dependent peptidase activity         2
    ## 19 GO:0004566                            beta-glucuronidase activity         2
    ##    Significant Expected    Fisher               type
    ## 1            6     0.18 4.200e-09 Biological.Process
    ## 2            1     0.02 1.800e-02 Biological.Process
    ## 3            1     0.02 1.800e-02 Biological.Process
    ## 4            1     0.02 1.800e-02 Biological.Process
    ## 5            1     0.04 3.700e-02 Biological.Process
    ## 6            1     0.04 3.700e-02 Biological.Process
    ## 7            2     0.04 4.300e-04 Molecular.Function
    ## 8            2     0.19 1.406e-02 Molecular.Function
    ## 9            1     0.02 2.089e-02 Molecular.Function
    ## 10           1     0.02 2.089e-02 Molecular.Function
    ## 11           1     0.02 2.089e-02 Molecular.Function
    ## 12           1     0.02 2.089e-02 Molecular.Function
    ## 13           1     0.02 2.089e-02 Molecular.Function
    ## 14           4     1.11 2.348e-02 Molecular.Function
    ## 15           1     0.04 4.136e-02 Molecular.Function
    ## 16           1     0.04 4.136e-02 Molecular.Function
    ## 17           1     0.04 4.136e-02 Molecular.Function
    ## 18           1     0.04 4.136e-02 Molecular.Function
    ## 19           1     0.04 4.136e-02 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 564 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  5 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  7 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  17 nodes to be scored   (161 eliminated genes)

    ## 
    ##   Level 10:  25 nodes to be scored   (223 eliminated genes)

    ## 
    ##   Level 9:   51 nodes to be scored   (276 eliminated genes)

    ## 
    ##   Level 8:   63 nodes to be scored   (396 eliminated genes)

    ## 
    ##   Level 7:   83 nodes to be scored   (489 eliminated genes)

    ## 
    ##   Level 6:   96 nodes to be scored   (690 eliminated genes)

    ## 
    ##   Level 5:   100 nodes to be scored  (1010 eliminated genes)

    ## 
    ##   Level 4:   63 nodes to be scored   (1151 eliminated genes)

    ## 
    ##   Level 3:   42 nodes to be scored   (1285 eliminated genes)

    ## 
    ##   Level 2:   11 nodes to be scored   (1339 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1406 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 349 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  5 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   18 nodes to be scored   (12 eliminated genes)

    ## 
    ##   Level 8:   24 nodes to be scored   (153 eliminated genes)

    ## 
    ##   Level 7:   48 nodes to be scored   (488 eliminated genes)

    ## 
    ##   Level 6:   77 nodes to be scored   (632 eliminated genes)

    ## 
    ##   Level 5:   67 nodes to be scored   (1020 eliminated genes)

    ## 
    ##   Level 4:   62 nodes to be scored   (1506 eliminated genes)

    ## 
    ##   Level 3:   33 nodes to be scored   (2140 eliminated genes)

    ## 
    ##   Level 2:   12 nodes to be scored   (2448 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2683 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0001822                                     kidney development        65
    ## 2  GO:0003341                                        cilium movement        20
    ## 3  GO:0002221         pattern recognition receptor signaling pathway        11
    ## 4  GO:0000082                  G1/S transition of mitotic cell cycle       104
    ## 5  GO:0005524                                            ATP binding       263
    ## 6  GO:0000030                           mannosyltransferase activity        12
    ## 7  GO:0003964                   RNA-directed DNA polymerase activity       122
    ## 8  GO:0004622                             lysophospholipase activity         4
    ## 9  GO:0000340                      RNA 7-methylguanosine cap binding         4
    ## 10 GO:0004571 mannosyl-oligosaccharide 1,2-alpha-mannosidase acti...         2
    ## 11 GO:0004342             glucosamine-6-phosphate deaminase activity         2
    ## 12 GO:0001965                        G-protein alpha-subunit binding        16
    ## 13 GO:0005085             guanyl-nucleotide exchange factor activity        12
    ##    Significant Expected    Fisher               type
    ## 1           19     8.84 0.0005700 Biological.Process
    ## 2            9     2.72 0.0005800 Biological.Process
    ## 3            5     1.50 0.0102400 Biological.Process
    ## 4           22    14.15 0.0182400 Biological.Process
    ## 5           66    39.13 0.0000032 Molecular.Function
    ## 6            7     1.79 0.0005900 Molecular.Function
    ## 7           29    18.15 0.0052900 Molecular.Function
    ## 8            3     0.60 0.0116400 Molecular.Function
    ## 9            3     0.60 0.0116400 Molecular.Function
    ## 10           2     0.30 0.0220900 Molecular.Function
    ## 11           2     0.30 0.0220900 Molecular.Function
    ## 12           6     2.38 0.0223300 Molecular.Function
    ## 13           5     1.79 0.0228800 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 2 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 324 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  5 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 11:  10 nodes to be scored   (143 eliminated genes)

    ## 
    ##   Level 10:  12 nodes to be scored   (212 eliminated genes)

    ## 
    ##   Level 9:   22 nodes to be scored   (239 eliminated genes)

    ## 
    ##   Level 8:   30 nodes to be scored   (322 eliminated genes)

    ## 
    ##   Level 7:   39 nodes to be scored   (380 eliminated genes)

    ## 
    ##   Level 6:   52 nodes to be scored   (574 eliminated genes)

    ## 
    ##   Level 5:   65 nodes to be scored   (831 eliminated genes)

    ## 
    ##   Level 4:   45 nodes to be scored   (944 eliminated genes)

    ## 
    ##   Level 3:   29 nodes to be scored   (1139 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (1307 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1380 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 255 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   11 nodes to be scored   (1 eliminated genes)

    ## 
    ##   Level 8:   20 nodes to be scored   (101 eliminated genes)

    ## 
    ##   Level 7:   37 nodes to be scored   (414 eliminated genes)

    ## 
    ##   Level 6:   63 nodes to be scored   (600 eliminated genes)

    ## 
    ##   Level 5:   46 nodes to be scored   (952 eliminated genes)

    ## 
    ##   Level 4:   41 nodes to be scored   (1424 eliminated genes)

    ## 
    ##   Level 3:   24 nodes to be scored   (2014 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (2284 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2601 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0000226                  microtubule cytoskeleton organization        43
    ## 2 GO:0006813                                potassium ion transport         1
    ## 3 GO:0000052                           citrulline metabolic process         1
    ## 4 GO:0000349 generation of catalytic spliceosome for first trans...         1
    ## 5 GO:0000390                       spliceosomal complex disassembly         1
    ## 6 GO:0003985                acetyl-CoA C-acetyltransferase activity         2
    ## 7 GO:0005272                                sodium channel activity         2
    ## 8 GO:0005524                                            ATP binding       263
    ##   Significant Expected Fisher               type
    ## 1           6     1.80 0.0280 Biological.Process
    ## 2           1     0.04 0.0420 Biological.Process
    ## 3           1     0.04 0.0420 Biological.Process
    ## 4           1     0.04 0.0420 Biological.Process
    ## 5           1     0.04 0.0420 Biological.Process
    ## 6           2     0.11 0.0033 Molecular.Function
    ## 7           2     0.11 0.0033 Molecular.Function
    ## 8          26    15.06 0.0033 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 168 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  4 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 10:  6 nodes to be scored    (203 eliminated genes)

    ## 
    ##   Level 9:   9 nodes to be scored    (205 eliminated genes)

    ## 
    ##   Level 8:   10 nodes to be scored   (248 eliminated genes)

    ## 
    ##   Level 7:   10 nodes to be scored   (271 eliminated genes)

    ## 
    ##   Level 6:   22 nodes to be scored   (319 eliminated genes)

    ## 
    ##   Level 5:   34 nodes to be scored   (371 eliminated genes)

    ## 
    ##   Level 4:   33 nodes to be scored   (504 eliminated genes)

    ## 
    ##   Level 3:   25 nodes to be scored   (778 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (943 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1146 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 142 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   6 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   10 nodes to be scored   (102 eliminated genes)

    ## 
    ##   Level 7:   21 nodes to be scored   (391 eliminated genes)

    ## 
    ##   Level 6:   28 nodes to be scored   (550 eliminated genes)

    ## 
    ##   Level 5:   30 nodes to be scored   (840 eliminated genes)

    ## 
    ##   Level 4:   22 nodes to be scored   (1084 eliminated genes)

    ## 
    ##   Level 3:   15 nodes to be scored   (1782 eliminated genes)

    ## 
    ##   Level 2:   7 nodes to be scored    (2054 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2402 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0001558                              regulation of cell growth         3
    ## 2  GO:0003844             1,4-alpha-glucan branching enzyme activity         1
    ## 3  GO:0004325                                ferrochelatase activity         1
    ## 4  GO:0004066 asparagine synthase (glutamine-hydrolyzing) activit...         1
    ## 5  GO:0004618                       phosphoglycerate kinase activity         1
    ## 6  GO:0004777 succinate-semialdehyde dehydrogenase (NAD+) activit...         1
    ## 7  GO:0004095              carnitine O-palmitoyltransferase activity         1
    ## 8  GO:0005391 P-type sodium:potassium-exchanging transporter acti...         2
    ## 9  GO:0003730                                    mRNA 3'-UTR binding         2
    ## 10 GO:0003987                            acetate-CoA ligase activity         2
    ## 11 GO:0004029                 aldehyde dehydrogenase (NAD+) activity         3
    ## 12 GO:0004031                              aldehyde oxidase activity         3
    ##    Significant Expected Fisher               type
    ## 1            1     0.03  0.027 Biological.Process
    ## 2            1     0.02  0.015 Molecular.Function
    ## 3            1     0.02  0.015 Molecular.Function
    ## 4            1     0.02  0.015 Molecular.Function
    ## 5            1     0.02  0.015 Molecular.Function
    ## 6            1     0.02  0.015 Molecular.Function
    ## 7            1     0.02  0.015 Molecular.Function
    ## 8            1     0.03  0.030 Molecular.Function
    ## 9            1     0.03  0.030 Molecular.Function
    ## 10           1     0.03  0.030 Molecular.Function
    ## 11           1     0.05  0.045 Molecular.Function
    ## 12           1     0.05  0.045 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 656 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  7 nodes to be scored    (5 eliminated genes)

    ## 
    ##   Level 12:  11 nodes to be scored   (20 eliminated genes)

    ## 
    ##   Level 11:  23 nodes to be scored   (172 eliminated genes)

    ## 
    ##   Level 10:  34 nodes to be scored   (227 eliminated genes)

    ## 
    ##   Level 9:   63 nodes to be scored   (275 eliminated genes)

    ## 
    ##   Level 8:   75 nodes to be scored   (390 eliminated genes)

    ## 
    ##   Level 7:   92 nodes to be scored   (532 eliminated genes)

    ## 
    ##   Level 6:   112 nodes to be scored  (707 eliminated genes)

    ## 
    ##   Level 5:   108 nodes to be scored  (1018 eliminated genes)

    ## 
    ##   Level 4:   65 nodes to be scored   (1098 eliminated genes)

    ## 
    ##   Level 3:   46 nodes to be scored   (1240 eliminated genes)

    ## 
    ##   Level 2:   13 nodes to be scored   (1344 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1406 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 419 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  9 nodes to be scored    (2 eliminated genes)

    ## 
    ##   Level 9:   21 nodes to be scored   (9 eliminated genes)

    ## 
    ##   Level 8:   40 nodes to be scored   (160 eliminated genes)

    ## 
    ##   Level 7:   64 nodes to be scored   (502 eliminated genes)

    ## 
    ##   Level 6:   94 nodes to be scored   (661 eliminated genes)

    ## 
    ##   Level 5:   68 nodes to be scored   (1036 eliminated genes)

    ## 
    ##   Level 4:   67 nodes to be scored   (1518 eliminated genes)

    ## 
    ##   Level 3:   37 nodes to be scored   (2150 eliminated genes)

    ## 
    ##   Level 2:   12 nodes to be scored   (2458 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2682 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0003341                                        cilium movement        20
    ## 2  GO:0000184 nuclear-transcribed mRNA catabolic process, nonsens...        10
    ## 3  GO:0002064                            epithelial cell development         4
    ## 4  GO:0000045                                 autophagosome assembly         7
    ## 5  GO:0007224                           smoothened signaling pathway         2
    ## 6  GO:0002790                                      peptide secretion         2
    ## 7  GO:0001933         negative regulation of protein phosphorylation         5
    ## 8  GO:0000038           very long-chain fatty acid metabolic process         9
    ## 9  GO:0005388                    P-type calcium transporter activity         4
    ## 10 GO:0003964                   RNA-directed DNA polymerase activity       122
    ## 11 GO:0003676                                   nucleic acid binding       497
    ## 12 GO:0005044                            scavenger receptor activity        13
    ## 13 GO:0005524                                            ATP binding       263
    ## 14 GO:0004345             glucose-6-phosphate dehydrogenase activity         2
    ## 15 GO:0004100                               chitin synthase activity         2
    ## 16 GO:0003730                                    mRNA 3'-UTR binding         2
    ## 17 GO:0004573 Glc3Man9GlcNAc2 oligosaccharide glucosidase activit...         2
    ## 18 GO:0003756                   protein disulfide isomerase activity         5
    ## 19 GO:0003723                                            RNA binding       128
    ##    Significant Expected Fisher               type
    ## 1            9     3.23 0.0060 Biological.Process
    ## 2            5     1.62 0.0130 Biological.Process
    ## 3            3     0.65 0.0150 Biological.Process
    ## 4            4     1.13 0.0160 Biological.Process
    ## 5            2     0.32 0.0260 Biological.Process
    ## 6            2     0.32 0.0260 Biological.Process
    ## 7            3     0.81 0.0320 Biological.Process
    ## 8            4     1.45 0.0430 Biological.Process
    ## 9            4     0.73 0.0011 Molecular.Function
    ## 10          35    22.41 0.0029 Molecular.Function
    ## 11         108    91.31 0.0033 Molecular.Function
    ## 12           6     2.39 0.0199 Molecular.Function
    ## 13          61    48.32 0.0228 Molecular.Function
    ## 14           2     0.37 0.0337 Molecular.Function
    ## 15           2     0.37 0.0337 Molecular.Function
    ## 16           2     0.37 0.0337 Molecular.Function
    ## 17           2     0.37 0.0337 Molecular.Function
    ## 18           3     0.92 0.0460 Molecular.Function
    ## 19          29    23.52 0.0465 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 2 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 412 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  8 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  13 nodes to be scored   (144 eliminated genes)

    ## 
    ##   Level 10:  19 nodes to be scored   (222 eliminated genes)

    ## 
    ##   Level 9:   33 nodes to be scored   (259 eliminated genes)

    ## 
    ##   Level 8:   36 nodes to be scored   (354 eliminated genes)

    ## 
    ##   Level 7:   48 nodes to be scored   (427 eliminated genes)

    ## 
    ##   Level 6:   72 nodes to be scored   (591 eliminated genes)

    ## 
    ##   Level 5:   76 nodes to be scored   (920 eliminated genes)

    ## 
    ##   Level 4:   55 nodes to be scored   (1051 eliminated genes)

    ## 
    ##   Level 3:   38 nodes to be scored   (1197 eliminated genes)

    ## 
    ##   Level 2:   11 nodes to be scored   (1332 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1399 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 279 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 9:   13 nodes to be scored   (3 eliminated genes)

    ## 
    ##   Level 8:   22 nodes to be scored   (146 eliminated genes)

    ## 
    ##   Level 7:   36 nodes to be scored   (467 eliminated genes)

    ## 
    ##   Level 6:   61 nodes to be scored   (607 eliminated genes)

    ## 
    ##   Level 5:   54 nodes to be scored   (942 eliminated genes)

    ## 
    ##   Level 4:   48 nodes to be scored   (1368 eliminated genes)

    ## 
    ##   Level 3:   27 nodes to be scored   (2003 eliminated genes)

    ## 
    ##   Level 2:   11 nodes to be scored   (2320 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2606 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0001508                                       action potential         4
    ## 2 GO:0005524                                            ATP binding       263
    ## 3 GO:0001409 guanine nucleotide transmembrane transporter activi...        30
    ## 4 GO:0003964                   RNA-directed DNA polymerase activity       122
    ## 5 GO:0004888              transmembrane signaling receptor activity       224
    ## 6 GO:0004972                       NMDA glutamate receptor activity         3
    ## 7 GO:0003676                                   nucleic acid binding       497
    ## 8 GO:0001965                        G-protein alpha-subunit binding        16
    ##   Significant Expected   Fisher               type
    ## 1           2     0.28 0.026000 Biological.Process
    ## 2          41    21.22 0.000014 Molecular.Function
    ## 3           9     2.42 0.000390 Molecular.Function
    ## 4          19     9.84 0.003440 Molecular.Function
    ## 5          18    18.07 0.015910 Molecular.Function
    ## 6           2     0.24 0.018420 Molecular.Function
    ## 7          42    40.10 0.019940 Molecular.Function
    ## 8           4     1.29 0.034670 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 322 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  8 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  14 nodes to be scored   (142 eliminated genes)

    ## 
    ##   Level 10:  16 nodes to be scored   (216 eliminated genes)

    ## 
    ##   Level 9:   22 nodes to be scored   (235 eliminated genes)

    ## 
    ##   Level 8:   26 nodes to be scored   (279 eliminated genes)

    ## 
    ##   Level 7:   38 nodes to be scored   (359 eliminated genes)

    ## 
    ##   Level 6:   51 nodes to be scored   (462 eliminated genes)

    ## 
    ##   Level 5:   62 nodes to be scored   (694 eliminated genes)

    ## 
    ##   Level 4:   43 nodes to be scored   (830 eliminated genes)

    ## 
    ##   Level 3:   29 nodes to be scored   (1098 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (1293 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1351 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 161 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   6 nodes to be scored    (9 eliminated genes)

    ## 
    ##   Level 8:   8 nodes to be scored    (142 eliminated genes)

    ## 
    ##   Level 7:   20 nodes to be scored   (442 eliminated genes)

    ## 
    ##   Level 6:   33 nodes to be scored   (472 eliminated genes)

    ## 
    ##   Level 5:   34 nodes to be scored   (782 eliminated genes)

    ## 
    ##   Level 4:   30 nodes to be scored   (1211 eliminated genes)

    ## 
    ##   Level 3:   18 nodes to be scored   (1913 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (2169 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2455 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0003341                                        cilium movement        20
    ## 2  GO:0015969             guanosine tetraphosphate metabolic process         9
    ## 3  GO:0000972 transcription-dependent tethering of RNA polymerase...         1
    ## 4  GO:0001731         formation of translation preinitiation complex         1
    ## 5  GO:0007175 negative regulation of epidermal growth factor-acti...         1
    ## 6  GO:0000390                       spliceosomal complex disassembly         1
    ## 7  GO:0006884                                cell volume homeostasis         1
    ## 8  GO:0006513                             protein monoubiquitination         1
    ## 9  GO:0002291 T cell activation via T cell receptor contact with ...        10
    ## 10 GO:0000082                  G1/S transition of mitotic cell cycle       104
    ## 11 GO:0001409 guanine nucleotide transmembrane transporter activi...        30
    ## 12 GO:0004315     3-oxoacyl-[acyl-carrier-protein] synthase activity         2
    ## 13 GO:0003964                   RNA-directed DNA polymerase activity       122
    ## 14 GO:0005524                                            ATP binding       263
    ## 15 GO:0003676                                   nucleic acid binding       497
    ## 16 GO:0004142      diacylglycerol cholinephosphotransferase activity         1
    ## 17 GO:0005080                               protein kinase C binding         1
    ## 18 GO:0005338    nucleotide-sugar transmembrane transporter activity         1
    ## 19 GO:0004591 oxoglutarate dehydrogenase (succinyl-transferring) ...         1
    ## 20 GO:0003863 3-methyl-2-oxobutanoate dehydrogenase (2-methylprop...         1
    ## 21 GO:0004022                  alcohol dehydrogenase (NAD+) activity         1
    ## 22 GO:0004013                        adenosylhomocysteinase activity         1
    ## 23 GO:0001758                         retinal dehydrogenase activity         1
    ##    Significant Expected Fisher               type
    ## 1            3     0.62 0.0220 Biological.Process
    ## 2            2     0.28 0.0300 Biological.Process
    ## 3            1     0.03 0.0310 Biological.Process
    ## 4            1     0.03 0.0310 Biological.Process
    ## 5            1     0.03 0.0310 Biological.Process
    ## 6            1     0.03 0.0310 Biological.Process
    ## 7            1     0.03 0.0310 Biological.Process
    ## 8            1     0.03 0.0310 Biological.Process
    ## 9            2     0.31 0.0360 Biological.Process
    ## 10           7     3.24 0.0380 Biological.Process
    ## 11           6     1.25 0.0012 Molecular.Function
    ## 12           2     0.08 0.0017 Molecular.Function
    ## 13          12     5.10 0.0042 Molecular.Function
    ## 14          18    10.99 0.0227 Molecular.Function
    ## 15          21    20.77 0.0368 Molecular.Function
    ## 16           1     0.04 0.0418 Molecular.Function
    ## 17           1     0.04 0.0418 Molecular.Function
    ## 18           1     0.04 0.0418 Molecular.Function
    ## 19           1     0.04 0.0418 Molecular.Function
    ## 20           1     0.04 0.0418 Molecular.Function
    ## 21           1     0.04 0.0418 Molecular.Function
    ## 22           1     0.04 0.0418 Molecular.Function
    ## 23           1     0.04 0.0418 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 271 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  6 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  9 nodes to be scored    (152 eliminated genes)

    ## 
    ##   Level 10:  11 nodes to be scored   (222 eliminated genes)

    ## 
    ##   Level 9:   19 nodes to be scored   (239 eliminated genes)

    ## 
    ##   Level 8:   21 nodes to be scored   (273 eliminated genes)

    ## 
    ##   Level 7:   29 nodes to be scored   (295 eliminated genes)

    ## 
    ##   Level 6:   43 nodes to be scored   (455 eliminated genes)

    ## 
    ##   Level 5:   51 nodes to be scored   (729 eliminated genes)

    ## 
    ##   Level 4:   39 nodes to be scored   (828 eliminated genes)

    ## 
    ##   Level 3:   30 nodes to be scored   (1036 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (1217 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1371 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 102 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   10 nodes to be scored   (263 eliminated genes)

    ## 
    ##   Level 6:   16 nodes to be scored   (334 eliminated genes)

    ## 
    ##   Level 5:   22 nodes to be scored   (557 eliminated genes)

    ## 
    ##   Level 4:   21 nodes to be scored   (816 eliminated genes)

    ## 
    ##   Level 3:   20 nodes to be scored   (1641 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (1929 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2508 eliminated genes)

    ##         GO.ID                                               Term Annotated
    ## 1  GO:0001822                                 kidney development        65
    ## 2  GO:0001818         negative regulation of cytokine production        18
    ## 3  GO:0000335 negative regulation of transposition, DNA-mediated         1
    ## 4  GO:0001913                       T cell mediated cytotoxicity         1
    ## 5  GO:0001696                             gastric acid secretion        10
    ## 6  GO:0005524                                        ATP binding       263
    ## 7  GO:0004017                          adenylate kinase activity         1
    ## 8  GO:0003682                                  chromatin binding        13
    ## 9  GO:0004100                           chitin synthase activity         2
    ## 10 GO:0004660               protein farnesyltransferase activity         2
    ## 11 GO:0001965                    G-protein alpha-subunit binding        16
    ##    Significant Expected  Fisher               type
    ## 1            7     1.93 0.00240 Biological.Process
    ## 2            3     0.54 0.01470 Biological.Process
    ## 3            1     0.03 0.02980 Biological.Process
    ## 4            1     0.03 0.02980 Biological.Process
    ## 5            2     0.30 0.03350 Biological.Process
    ## 6           16     5.87 0.00013 Molecular.Function
    ## 7            1     0.02 0.02233 Molecular.Function
    ## 8            2     0.29 0.03268 Molecular.Function
    ## 9            1     0.04 0.04418 Molecular.Function
    ## 10           1     0.04 0.04418 Molecular.Function
    ## 11           2     0.36 0.04818 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 269 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  4 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 11:  8 nodes to be scored    (2 eliminated genes)

    ## 
    ##   Level 10:  9 nodes to be scored    (79 eliminated genes)

    ## 
    ##   Level 9:   17 nodes to be scored   (96 eliminated genes)

    ## 
    ##   Level 8:   22 nodes to be scored   (146 eliminated genes)

    ## 
    ##   Level 7:   30 nodes to be scored   (181 eliminated genes)

    ## 
    ##   Level 6:   43 nodes to be scored   (321 eliminated genes)

    ## 
    ##   Level 5:   59 nodes to be scored   (637 eliminated genes)

    ## 
    ##   Level 4:   42 nodes to be scored   (812 eliminated genes)

    ## 
    ##   Level 3:   23 nodes to be scored   (1054 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (1237 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1291 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 150 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   5 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   9 nodes to be scored    (7 eliminated genes)

    ## 
    ##   Level 7:   18 nodes to be scored   (445 eliminated genes)

    ## 
    ##   Level 6:   32 nodes to be scored   (560 eliminated genes)

    ## 
    ##   Level 5:   31 nodes to be scored   (785 eliminated genes)

    ## 
    ##   Level 4:   26 nodes to be scored   (1169 eliminated genes)

    ## 
    ##   Level 3:   17 nodes to be scored   (1813 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (2105 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2425 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0015969             guanosine tetraphosphate metabolic process         9
    ## 2  GO:0001731         formation of translation preinitiation complex         1
    ## 3  GO:0000349 generation of catalytic spliceosome for first trans...         1
    ## 4  GO:0016311                                      dephosphorylation         1
    ## 5  GO:0006091         generation of precursor metabolites and energy        12
    ## 6  GO:0000079 regulation of cyclin-dependent protein serine/threo...         2
    ## 7  GO:0001569       branching involved in blood vessel morphogenesis         2
    ## 8  GO:0002790                                      peptide secretion         2
    ## 9  GO:0005302          L-tyrosine transmembrane transporter activity         9
    ## 10 GO:0003824                                     catalytic activity      1042
    ## 11 GO:0004140                          dephospho-CoA kinase activity         1
    ## 12 GO:0000822                      inositol hexakisphosphate binding         1
    ## 13 GO:0004334                           fumarylacetoacetase activity         1
    ## 14 GO:0004095              carnitine O-palmitoyltransferase activity         1
    ## 15 GO:0001054                              RNA polymerase I activity         1
    ##    Significant Expected Fisher               type
    ## 1            2     0.19 0.0140 Biological.Process
    ## 2            1     0.02 0.0210 Biological.Process
    ## 3            1     0.02 0.0210 Biological.Process
    ## 4            1     0.02 0.0210 Biological.Process
    ## 5            2     0.26 0.0250 Biological.Process
    ## 6            1     0.04 0.0420 Biological.Process
    ## 7            1     0.04 0.0420 Biological.Process
    ## 8            1     0.04 0.0420 Biological.Process
    ## 9            3     0.24 0.0014 Molecular.Function
    ## 10          30    28.15 0.0122 Molecular.Function
    ## 11           1     0.03 0.0270 Molecular.Function
    ## 12           1     0.03 0.0270 Molecular.Function
    ## 13           1     0.03 0.0270 Molecular.Function
    ## 14           1     0.03 0.0270 Molecular.Function
    ## 15           1     0.03 0.0270 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 356 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  8 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 11:  11 nodes to be scored   (144 eliminated genes)

    ## 
    ##   Level 10:  17 nodes to be scored   (219 eliminated genes)

    ## 
    ##   Level 9:   31 nodes to be scored   (236 eliminated genes)

    ## 
    ##   Level 8:   34 nodes to be scored   (330 eliminated genes)

    ## 
    ##   Level 7:   46 nodes to be scored   (398 eliminated genes)

    ## 
    ##   Level 6:   58 nodes to be scored   (519 eliminated genes)

    ## 
    ##   Level 5:   61 nodes to be scored   (761 eliminated genes)

    ## 
    ##   Level 4:   45 nodes to be scored   (882 eliminated genes)

    ## 
    ##   Level 3:   30 nodes to be scored   (1085 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (1273 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1372 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 212 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   7 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   10 nodes to be scored   (100 eliminated genes)

    ## 
    ##   Level 7:   27 nodes to be scored   (429 eliminated genes)

    ## 
    ##   Level 6:   48 nodes to be scored   (532 eliminated genes)

    ## 
    ##   Level 5:   49 nodes to be scored   (884 eliminated genes)

    ## 
    ##   Level 4:   42 nodes to be scored   (1288 eliminated genes)

    ## 
    ##   Level 3:   19 nodes to be scored   (1941 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (2226 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2562 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0002218                   activation of innate immune response        15
    ## 2  GO:0006397                                        mRNA processing        64
    ## 3  GO:0006813                                potassium ion transport         1
    ## 4  GO:0006601                          creatine biosynthetic process         1
    ## 5  GO:0000349 generation of catalytic spliceosome for first trans...         1
    ## 6  GO:0001731         formation of translation preinitiation complex         1
    ## 7  GO:0002265       astrocyte activation involved in immune response         1
    ## 8  GO:0001113 transcription open complex formation at RNA polymer...         1
    ## 9  GO:0002230 positive regulation of defense response to virus by...         1
    ## 10 GO:0009071             serine family amino acid catabolic process         1
    ## 11 GO:0016311                                      dephosphorylation         1
    ## 12 GO:0003730                                    mRNA 3'-UTR binding         2
    ## 13 GO:0003676                                   nucleic acid binding       497
    ## 14 GO:0000016                                       lactase activity         4
    ## 15 GO:0004222                          metalloendopeptidase activity        15
    ## 16 GO:0004458          D-lactate dehydrogenase (cytochrome) activity         1
    ## 17 GO:0004334                           fumarylacetoacetase activity         1
    ## 18 GO:0002161                        aminoacyl-tRNA editing activity         1
    ## 19 GO:0015101      organic cation transmembrane transporter activity         1
    ## 20 GO:0004516          nicotinate phosphoribosyltransferase activity         1
    ## 21 GO:0008137               NADH dehydrogenase (ubiquinone) activity         1
    ## 22 GO:0004618                       phosphoglycerate kinase activity         1
    ## 23 GO:0004642    phosphoribosylformylglycinamidine synthase activity         1
    ## 24 GO:0004325                                ferrochelatase activity         1
    ## 25 GO:0008810                                     cellulase activity         1
    ## 26 GO:0004777 succinate-semialdehyde dehydrogenase (NAD+) activit...         1
    ## 27 GO:0003997                              acyl-CoA oxidase activity         1
    ## 28 GO:0004637            phosphoribosylamine-glycine ligase activity         1
    ## 29 GO:0005524                                            ATP binding       263
    ##    Significant Expected Fisher               type
    ## 1            4     0.52 0.0013 Biological.Process
    ## 2            4     2.22 0.0341 Biological.Process
    ## 3            1     0.03 0.0347 Biological.Process
    ## 4            1     0.03 0.0347 Biological.Process
    ## 5            1     0.03 0.0347 Biological.Process
    ## 6            1     0.03 0.0347 Biological.Process
    ## 7            1     0.03 0.0347 Biological.Process
    ## 8            1     0.03 0.0347 Biological.Process
    ## 9            1     0.03 0.0347 Biological.Process
    ## 10           1     0.03 0.0347 Biological.Process
    ## 11           1     0.03 0.0347 Biological.Process
    ## 12           2     0.08 0.0017 Molecular.Function
    ## 13          28    20.77 0.0097 Molecular.Function
    ## 14           2     0.17 0.0098 Molecular.Function
    ## 15           3     0.63 0.0224 Molecular.Function
    ## 16           1     0.04 0.0418 Molecular.Function
    ## 17           1     0.04 0.0418 Molecular.Function
    ## 18           1     0.04 0.0418 Molecular.Function
    ## 19           1     0.04 0.0418 Molecular.Function
    ## 20           1     0.04 0.0418 Molecular.Function
    ## 21           1     0.04 0.0418 Molecular.Function
    ## 22           1     0.04 0.0418 Molecular.Function
    ## 23           1     0.04 0.0418 Molecular.Function
    ## 24           1     0.04 0.0418 Molecular.Function
    ## 25           1     0.04 0.0418 Molecular.Function
    ## 26           1     0.04 0.0418 Molecular.Function
    ## 27           1     0.04 0.0418 Molecular.Function
    ## 28           1     0.04 0.0418 Molecular.Function
    ## 29          17    10.99 0.0431 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 519 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  8 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 12:  8 nodes to be scored    (15 eliminated genes)

    ## 
    ##   Level 11:  18 nodes to be scored   (173 eliminated genes)

    ## 
    ##   Level 10:  26 nodes to be scored   (217 eliminated genes)

    ## 
    ##   Level 9:   49 nodes to be scored   (253 eliminated genes)

    ## 
    ##   Level 8:   51 nodes to be scored   (353 eliminated genes)

    ## 
    ##   Level 7:   69 nodes to be scored   (494 eliminated genes)

    ## 
    ##   Level 6:   84 nodes to be scored   (630 eliminated genes)

    ## 
    ##   Level 5:   91 nodes to be scored   (908 eliminated genes)

    ## 
    ##   Level 4:   58 nodes to be scored   (1001 eliminated genes)

    ## 
    ##   Level 3:   42 nodes to be scored   (1202 eliminated genes)

    ## 
    ##   Level 2:   11 nodes to be scored   (1316 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1394 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 284 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  1 nodes to be scored    (2 eliminated genes)

    ## 
    ##   Level 10:  5 nodes to be scored    (2 eliminated genes)

    ## 
    ##   Level 9:   14 nodes to be scored   (4 eliminated genes)

    ## 
    ##   Level 8:   19 nodes to be scored   (151 eliminated genes)

    ## 
    ##   Level 7:   36 nodes to be scored   (480 eliminated genes)

    ## 
    ##   Level 6:   63 nodes to be scored   (579 eliminated genes)

    ## 
    ##   Level 5:   56 nodes to be scored   (929 eliminated genes)

    ## 
    ##   Level 4:   48 nodes to be scored   (1361 eliminated genes)

    ## 
    ##   Level 3:   29 nodes to be scored   (2104 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (2386 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2632 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0001732 formation of cytoplasmic translation initiation com...         6
    ## 2 GO:0000054                  ribosomal subunit export from nucleus         3
    ## 3 GO:0003676                                   nucleic acid binding       497
    ## 4 GO:0003727                            single-stranded RNA binding         2
    ## 5 GO:0003723                                            RNA binding       128
    ## 6 GO:0003964                   RNA-directed DNA polymerase activity       122
    ## 7 GO:0005198                           structural molecule activity        54
    ##   Significant Expected    Fisher               type
    ## 1           4     0.45 0.0004200 Biological.Process
    ## 2           2     0.23 0.0162500 Biological.Process
    ## 3          59    39.03 0.0000034 Molecular.Function
    ## 4           2     0.16 0.0061000 Molecular.Function
    ## 5          20    10.05 0.0074000 Molecular.Function
    ## 6          17     9.58 0.0130000 Molecular.Function
    ## 7           9     4.24 0.0343000 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 167 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  4 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 10:  6 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 9:   8 nodes to be scored    (148 eliminated genes)

    ## 
    ##   Level 8:   8 nodes to be scored    (201 eliminated genes)

    ## 
    ##   Level 7:   15 nodes to be scored   (223 eliminated genes)

    ## 
    ##   Level 6:   29 nodes to be scored   (272 eliminated genes)

    ## 
    ##   Level 5:   40 nodes to be scored   (485 eliminated genes)

    ## 
    ##   Level 4:   28 nodes to be scored   (577 eliminated genes)

    ## 
    ##   Level 3:   17 nodes to be scored   (878 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (1078 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1245 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 75 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   6 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 6:   10 nodes to be scored   (10 eliminated genes)

    ## 
    ##   Level 5:   16 nodes to be scored   (46 eliminated genes)

    ## 
    ##   Level 4:   18 nodes to be scored   (225 eliminated genes)

    ## 
    ##   Level 3:   15 nodes to be scored   (1356 eliminated genes)

    ## 
    ##   Level 2:   4 nodes to be scored    (1846 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2434 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0006884                                cell volume homeostasis         1
    ## 2 GO:0001508                                       action potential         4
    ## 3 GO:0001510                                        RNA methylation         4
    ## 4 GO:0004095              carnitine O-palmitoyltransferase activity         1
    ## 5 GO:0004494                      methylmalonyl-CoA mutase activity         1
    ## 6 GO:0001055                             RNA polymerase II activity         1
    ## 7 GO:0004435          phosphatidylinositol phospholipase C activity         1
    ## 8 GO:0004591 oxoglutarate dehydrogenase (succinyl-transferring) ...         1
    ## 9 GO:0016491                                oxidoreductase activity       115
    ##   Significant Expected Fisher               type
    ## 1           1     0.01 0.0085 Biological.Process
    ## 2           1     0.03 0.0336 Biological.Process
    ## 3           1     0.03 0.0336 Biological.Process
    ## 4           1     0.01 0.0090 Molecular.Function
    ## 5           1     0.01 0.0090 Molecular.Function
    ## 6           1     0.01 0.0090 Molecular.Function
    ## 7           1     0.01 0.0090 Molecular.Function
    ## 8           1     0.01 0.0090 Molecular.Function
    ## 9           5     1.04 0.0220 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 222 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  4 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 10:  5 nodes to be scored    (147 eliminated genes)

    ## 
    ##   Level 9:   10 nodes to be scored   (149 eliminated genes)

    ## 
    ##   Level 8:   18 nodes to be scored   (180 eliminated genes)

    ## 
    ##   Level 7:   24 nodes to be scored   (244 eliminated genes)

    ## 
    ##   Level 6:   41 nodes to be scored   (396 eliminated genes)

    ## 
    ##   Level 5:   50 nodes to be scored   (560 eliminated genes)

    ## 
    ##   Level 4:   33 nodes to be scored   (768 eliminated genes)

    ## 
    ##   Level 3:   23 nodes to be scored   (1004 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (1238 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1370 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 193 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  1 nodes to be scored    (2 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (2 eliminated genes)

    ## 
    ##   Level 9:   8 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 8:   12 nodes to be scored   (153 eliminated genes)

    ## 
    ##   Level 7:   23 nodes to be scored   (465 eliminated genes)

    ## 
    ##   Level 6:   37 nodes to be scored   (569 eliminated genes)

    ## 
    ##   Level 5:   40 nodes to be scored   (803 eliminated genes)

    ## 
    ##   Level 4:   36 nodes to be scored   (1160 eliminated genes)

    ## 
    ##   Level 3:   22 nodes to be scored   (2007 eliminated genes)

    ## 
    ##   Level 2:   7 nodes to be scored    (2281 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2579 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0002221         pattern recognition receptor signaling pathway        11
    ## 2  GO:0001666                                    response to hypoxia        14
    ## 3  GO:0006556              S-adenosylmethionine biosynthetic process         1
    ## 4  GO:0001937 negative regulation of endothelial cell proliferati...         1
    ## 5  GO:0005198                           structural molecule activity        54
    ## 6  GO:0003756                   protein disulfide isomerase activity         5
    ## 7  GO:0004332                fructose-bisphosphate aldolase activity         1
    ## 8  GO:0001594                          trace-amine receptor activity         1
    ## 9  GO:0004305                           ethanolamine kinase activity         1
    ## 10 GO:0004352                glutamate dehydrogenase (NAD+) activity         1
    ## 11 GO:0004561                 alpha-N-acetylglucosaminidase activity         1
    ## 12 GO:0004568                                     chitinase activity         1
    ## 13 GO:0004321                       fatty-acyl-CoA synthase activity         1
    ##    Significant Expected   Fisher               type
    ## 1            5     0.32 6.60e-06 Biological.Process
    ## 2            3     0.41 6.60e-03 Biological.Process
    ## 3            1     0.03 2.91e-02 Biological.Process
    ## 4            1     0.03 2.91e-02 Biological.Process
    ## 5            6     1.85 7.60e-03 Molecular.Function
    ## 6            2     0.17 1.08e-02 Molecular.Function
    ## 7            1     0.03 3.42e-02 Molecular.Function
    ## 8            1     0.03 3.42e-02 Molecular.Function
    ## 9            1     0.03 3.42e-02 Molecular.Function
    ## 10           1     0.03 3.42e-02 Molecular.Function
    ## 11           1     0.03 3.42e-02 Molecular.Function
    ## 12           1     0.03 3.42e-02 Molecular.Function
    ## 13           1     0.03 3.42e-02 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 317 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  3 nodes to be scored    (11 eliminated genes)

    ## 
    ##   Level 11:  7 nodes to be scored    (155 eliminated genes)

    ## 
    ##   Level 10:  13 nodes to be scored   (203 eliminated genes)

    ## 
    ##   Level 9:   23 nodes to be scored   (227 eliminated genes)

    ## 
    ##   Level 8:   25 nodes to be scored   (301 eliminated genes)

    ## 
    ##   Level 7:   40 nodes to be scored   (386 eliminated genes)

    ## 
    ##   Level 6:   55 nodes to be scored   (570 eliminated genes)

    ## 
    ##   Level 5:   63 nodes to be scored   (813 eliminated genes)

    ## 
    ##   Level 4:   45 nodes to be scored   (909 eliminated genes)

    ## 
    ##   Level 3:   29 nodes to be scored   (1097 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (1296 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1382 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 230 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   12 nodes to be scored   (9 eliminated genes)

    ## 
    ##   Level 8:   17 nodes to be scored   (151 eliminated genes)

    ## 
    ##   Level 7:   30 nodes to be scored   (478 eliminated genes)

    ## 
    ##   Level 6:   50 nodes to be scored   (548 eliminated genes)

    ## 
    ##   Level 5:   42 nodes to be scored   (844 eliminated genes)

    ## 
    ##   Level 4:   40 nodes to be scored   (1314 eliminated genes)

    ## 
    ##   Level 3:   22 nodes to be scored   (1866 eliminated genes)

    ## 
    ##   Level 2:   11 nodes to be scored   (2159 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2607 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0019700                  organic phosphonate catabolic process        12
    ## 2  GO:0003179                              heart valve morphogenesis         4
    ## 3  GO:0001889                                      liver development         8
    ## 4  GO:0006021                          inositol biosynthetic process         1
    ## 5  GO:0002230 positive regulation of defense response to virus by...         1
    ## 6  GO:0006884                                cell volume homeostasis         1
    ## 7  GO:0000014     single-stranded DNA endodeoxyribonuclease activity         3
    ## 8  GO:0005262                               calcium channel activity        24
    ## 9  GO:0001409 guanine nucleotide transmembrane transporter activi...        30
    ## 10 GO:0004252                     serine-type endopeptidase activity        53
    ## 11 GO:0003676                                   nucleic acid binding       497
    ## 12 GO:0005245                 voltage-gated calcium channel activity         7
    ##    Significant Expected  Fisher               type
    ## 1            5     0.43 0.00003 Biological.Process
    ## 2            2     0.14 0.00710 Biological.Process
    ## 3            2     0.28 0.03010 Biological.Process
    ## 4            1     0.04 0.03540 Biological.Process
    ## 5            1     0.04 0.03540 Biological.Process
    ## 6            1     0.04 0.03540 Biological.Process
    ## 7            2     0.16 0.00820 Molecular.Function
    ## 8            6     1.28 0.01030 Molecular.Function
    ## 9            5     1.60 0.01950 Molecular.Function
    ## 10           7     2.83 0.02080 Molecular.Function
    ## 11          23    26.50 0.03810 Molecular.Function
    ## 12           2     0.37 0.04970 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 819 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  6 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  12 nodes to be scored   (6 eliminated genes)

    ## 
    ##   Level 12:  19 nodes to be scored   (23 eliminated genes)

    ## 
    ##   Level 11:  32 nodes to be scored   (179 eliminated genes)

    ## 
    ##   Level 10:  47 nodes to be scored   (244 eliminated genes)

    ## 
    ##   Level 9:   83 nodes to be scored   (291 eliminated genes)

    ## 
    ##   Level 8:   92 nodes to be scored   (416 eliminated genes)

    ## 
    ##   Level 7:   116 nodes to be scored  (578 eliminated genes)

    ## 
    ##   Level 6:   142 nodes to be scored  (794 eliminated genes)

    ## 
    ##   Level 5:   128 nodes to be scored  (1090 eliminated genes)

    ## 
    ##   Level 4:   76 nodes to be scored   (1176 eliminated genes)

    ## 
    ##   Level 3:   50 nodes to be scored   (1311 eliminated genes)

    ## 
    ##   Level 2:   12 nodes to be scored   (1376 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1408 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 495 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  5 nodes to be scored    (2 eliminated genes)

    ## 
    ##   Level 10:  9 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 9:   25 nodes to be scored   (18 eliminated genes)

    ## 
    ##   Level 8:   45 nodes to be scored   (163 eliminated genes)

    ## 
    ##   Level 7:   85 nodes to be scored   (510 eliminated genes)

    ## 
    ##   Level 6:   117 nodes to be scored  (669 eliminated genes)

    ## 
    ##   Level 5:   84 nodes to be scored   (1077 eliminated genes)

    ## 
    ##   Level 4:   71 nodes to be scored   (1590 eliminated genes)

    ## 
    ##   Level 3:   36 nodes to be scored   (2207 eliminated genes)

    ## 
    ##   Level 2:   13 nodes to be scored   (2505 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2680 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0002064                            epithelial cell development         4
    ## 2  GO:0000184 nuclear-transcribed mRNA catabolic process, nonsens...        10
    ## 3  GO:0001732 formation of cytoplasmic translation initiation com...         6
    ## 4  GO:0000302                    response to reactive oxygen species         6
    ## 5  GO:0000398                         mRNA splicing, via spliceosome        63
    ## 6  GO:0000086                  G2/M transition of mitotic cell cycle         7
    ## 7  GO:0007224                           smoothened signaling pathway         2
    ## 8  GO:0003964                   RNA-directed DNA polymerase activity       122
    ## 9  GO:0003676                                   nucleic acid binding       497
    ## 10 GO:0005388                    P-type calcium transporter activity         4
    ## 11 GO:0005164                 tumor necrosis factor receptor binding         3
    ## 12 GO:0003735                     structural constituent of ribosome        13
    ## 13 GO:0005385            zinc ion transmembrane transporter activity         4
    ## 14 GO:0004623                              phospholipase A2 activity         4
    ## 15 GO:0004517                         nitric-oxide synthase activity         2
    ## 16 GO:0004365 glyceraldehyde-3-phosphate dehydrogenase (NAD+) (ph...         2
    ## 17 GO:0005112                                          Notch binding         2
    ## 18 GO:0003727                            single-stranded RNA binding         2
    ## 19 GO:0003723                                            RNA binding       128
    ##    Significant Expected  Fisher               type
    ## 1            4     0.86 0.00210 Biological.Process
    ## 2            6     2.15 0.00910 Biological.Process
    ## 3            4     1.29 0.02190 Biological.Process
    ## 4            4     1.29 0.02190 Biological.Process
    ## 5           19    13.57 0.03050 Biological.Process
    ## 6            4     1.51 0.04260 Biological.Process
    ## 7            2     0.43 0.04630 Biological.Process
    ## 8           41    26.11 0.00096 Molecular.Function
    ## 9          122   106.35 0.00113 Molecular.Function
    ## 10           4     0.86 0.00208 Molecular.Function
    ## 11           3     0.64 0.00976 Molecular.Function
    ## 12           7     2.78 0.01013 Molecular.Function
    ## 13           3     0.86 0.03279 Molecular.Function
    ## 14           3     0.86 0.03279 Molecular.Function
    ## 15           2     0.43 0.04573 Molecular.Function
    ## 16           2     0.43 0.04573 Molecular.Function
    ## 17           2     0.43 0.04573 Molecular.Function
    ## 18           2     0.43 0.04573 Molecular.Function
    ## 19          34    27.39 0.04878 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 302 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  5 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  11 nodes to be scored   (140 eliminated genes)

    ## 
    ##   Level 10:  12 nodes to be scored   (150 eliminated genes)

    ## 
    ##   Level 9:   19 nodes to be scored   (176 eliminated genes)

    ## 
    ##   Level 8:   21 nodes to be scored   (220 eliminated genes)

    ## 
    ##   Level 7:   29 nodes to be scored   (278 eliminated genes)

    ## 
    ##   Level 6:   45 nodes to be scored   (444 eliminated genes)

    ## 
    ##   Level 5:   64 nodes to be scored   (701 eliminated genes)

    ## 
    ##   Level 4:   48 nodes to be scored   (818 eliminated genes)

    ## 
    ##   Level 3:   34 nodes to be scored   (1115 eliminated genes)

    ## 
    ##   Level 2:   12 nodes to be scored   (1279 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1378 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 183 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   7 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   8 nodes to be scored    (148 eliminated genes)

    ## 
    ##   Level 7:   22 nodes to be scored   (444 eliminated genes)

    ## 
    ##   Level 6:   38 nodes to be scored   (469 eliminated genes)

    ## 
    ##   Level 5:   39 nodes to be scored   (720 eliminated genes)

    ## 
    ##   Level 4:   31 nodes to be scored   (1057 eliminated genes)

    ## 
    ##   Level 3:   23 nodes to be scored   (1811 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (2063 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2608 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0006886                        intracellular protein transport         5
    ## 2  GO:0001731         formation of translation preinitiation complex         1
    ## 3  GO:0000904         cell morphogenesis involved in differentiation         1
    ## 4  GO:0001937 negative regulation of endothelial cell proliferati...         1
    ## 5  GO:0000712      resolution of meiotic recombination intermediates         1
    ## 6  GO:0008218                                        bioluminescence        10
    ## 7  GO:0001965                        G-protein alpha-subunit binding        16
    ## 8  GO:0004149 dihydrolipoyllysine-residue succinyltransferase act...         1
    ## 9  GO:0005315 inorganic phosphate transmembrane transporter activ...         1
    ## 10 GO:0004748 ribonucleoside-diphosphate reductase activity, thio...         1
    ## 11 GO:0004017                              adenylate kinase activity         1
    ## 12 GO:0004108                         citrate (Si)-synthase activity         1
    ## 13 GO:0004152                  dihydroorotate dehydrogenase activity         1
    ## 14 GO:0004505                 phenylalanine 4-monooxygenase activity         1
    ## 15 GO:0004615                            phosphomannomutase activity         1
    ## 16 GO:0003874           6-pyruvoyltetrahydropterin synthase activity         1
    ## 17 GO:0004067                                  asparaginase activity         1
    ##    Significant Expected Fisher               type
    ## 1            2     0.15 0.0078 Biological.Process
    ## 2            1     0.03 0.0291 Biological.Process
    ## 3            1     0.03 0.0291 Biological.Process
    ## 4            1     0.03 0.0291 Biological.Process
    ## 5            1     0.03 0.0291 Biological.Process
    ## 6            2     0.29 0.0320 Biological.Process
    ## 7            4     0.50 0.0012 Molecular.Function
    ## 8            1     0.03 0.0313 Molecular.Function
    ## 9            1     0.03 0.0313 Molecular.Function
    ## 10           1     0.03 0.0313 Molecular.Function
    ## 11           1     0.03 0.0313 Molecular.Function
    ## 12           1     0.03 0.0313 Molecular.Function
    ## 13           1     0.03 0.0313 Molecular.Function
    ## 14           1     0.03 0.0313 Molecular.Function
    ## 15           1     0.03 0.0313 Molecular.Function
    ## 16           1     0.03 0.0313 Molecular.Function
    ## 17           1     0.03 0.0313 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 488 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  5 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 12:  10 nodes to be scored   (2 eliminated genes)

    ## 
    ##   Level 11:  19 nodes to be scored   (147 eliminated genes)

    ## 
    ##   Level 10:  22 nodes to be scored   (221 eliminated genes)

    ## 
    ##   Level 9:   46 nodes to be scored   (255 eliminated genes)

    ## 
    ##   Level 8:   48 nodes to be scored   (371 eliminated genes)

    ## 
    ##   Level 7:   61 nodes to be scored   (483 eliminated genes)

    ## 
    ##   Level 6:   81 nodes to be scored   (610 eliminated genes)

    ## 
    ##   Level 5:   83 nodes to be scored   (867 eliminated genes)

    ## 
    ##   Level 4:   57 nodes to be scored   (1040 eliminated genes)

    ## 
    ##   Level 3:   39 nodes to be scored   (1188 eliminated genes)

    ## 
    ##   Level 2:   13 nodes to be scored   (1297 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1405 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 335 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  6 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   16 nodes to be scored   (1 eliminated genes)

    ## 
    ##   Level 8:   26 nodes to be scored   (155 eliminated genes)

    ## 
    ##   Level 7:   42 nodes to be scored   (474 eliminated genes)

    ## 
    ##   Level 6:   76 nodes to be scored   (616 eliminated genes)

    ## 
    ##   Level 5:   66 nodes to be scored   (973 eliminated genes)

    ## 
    ##   Level 4:   60 nodes to be scored   (1437 eliminated genes)

    ## 
    ##   Level 3:   29 nodes to be scored   (2080 eliminated genes)

    ## 
    ##   Level 2:   12 nodes to be scored   (2399 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2635 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0019700                  organic phosphonate catabolic process        12
    ## 2 GO:0002064                            epithelial cell development         4
    ## 3 GO:0002218                   activation of innate immune response        15
    ## 4 GO:0001707                                     mesoderm formation         4
    ## 5 GO:0003179                              heart valve morphogenesis         4
    ## 6 GO:0008218                                        bioluminescence        10
    ## 7 GO:0003682                                      chromatin binding        13
    ## 8 GO:0001409 guanine nucleotide transmembrane transporter activi...        30
    ## 9 GO:0005245                 voltage-gated calcium channel activity         7
    ##   Significant Expected  Fisher               type
    ## 1           6     0.97 0.00015 Biological.Process
    ## 2           3     0.32 0.00194 Biological.Process
    ## 3           5     1.21 0.00491 Biological.Process
    ## 4           2     0.32 0.03485 Biological.Process
    ## 5           2     0.32 0.03485 Biological.Process
    ## 6           3     0.81 0.04051 Biological.Process
    ## 7           5     1.28 0.00600 Molecular.Function
    ## 8           7     2.96 0.02400 Molecular.Function
    ## 9           3     0.69 0.02500 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 253 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  9 nodes to be scored    (155 eliminated genes)

    ## 
    ##   Level 10:  10 nodes to be scored   (208 eliminated genes)

    ## 
    ##   Level 9:   16 nodes to be scored   (240 eliminated genes)

    ## 
    ##   Level 8:   19 nodes to be scored   (284 eliminated genes)

    ## 
    ##   Level 7:   29 nodes to be scored   (305 eliminated genes)

    ## 
    ##   Level 6:   40 nodes to be scored   (440 eliminated genes)

    ## 
    ##   Level 5:   50 nodes to be scored   (629 eliminated genes)

    ## 
    ##   Level 4:   39 nodes to be scored   (789 eliminated genes)

    ## 
    ##   Level 3:   25 nodes to be scored   (1062 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (1250 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1371 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 227 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  6 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   13 nodes to be scored   (10 eliminated genes)

    ## 
    ##   Level 8:   18 nodes to be scored   (155 eliminated genes)

    ## 
    ##   Level 7:   31 nodes to be scored   (470 eliminated genes)

    ## 
    ##   Level 6:   41 nodes to be scored   (603 eliminated genes)

    ## 
    ##   Level 5:   43 nodes to be scored   (861 eliminated genes)

    ## 
    ##   Level 4:   36 nodes to be scored   (1218 eliminated genes)

    ## 
    ##   Level 3:   24 nodes to be scored   (1918 eliminated genes)

    ## 
    ##   Level 2:   12 nodes to be scored   (2229 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2569 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0003341                                        cilium movement        20
    ## 2  GO:0000082                  G1/S transition of mitotic cell cycle       104
    ## 3  GO:0000972 transcription-dependent tethering of RNA polymerase...         1
    ## 4  GO:0002224                   toll-like receptor signaling pathway         1
    ## 5  GO:0006884                                cell volume homeostasis         1
    ## 6  GO:0000340                      RNA 7-methylguanosine cap binding         4
    ## 7  GO:0004252                     serine-type endopeptidase activity        53
    ## 8  GO:0005245                 voltage-gated calcium channel activity         7
    ## 9  GO:0004714 transmembrane receptor protein tyrosine kinase acti...         2
    ## 10 GO:0004458          D-lactate dehydrogenase (cytochrome) activity         1
    ## 11 GO:0002161                        aminoacyl-tRNA editing activity         1
    ## 12 GO:0003980 UDP-glucose:glycoprotein glucosyltransferase activi...         1
    ## 13 GO:0004022                  alcohol dehydrogenase (NAD+) activity         1
    ## 14 GO:0004485               methylcrotonoyl-CoA carboxylase activity         1
    ## 15 GO:0001758                         retinal dehydrogenase activity         1
    ## 16 GO:0004334                           fumarylacetoacetase activity         1
    ## 17 GO:0004427             inorganic diphosphate phosphatase activity         1
    ## 18 GO:0042802                              identical protein binding         1
    ## 19 GO:0005005                 transmembrane-ephrin receptor activity         1
    ## 20 GO:0004368 glycerol-3-phosphate dehydrogenase (quinone) activi...         1
    ## 21 GO:0005381            iron ion transmembrane transporter activity         1
    ## 22 GO:0004149 dihydrolipoyllysine-residue succinyltransferase act...         1
    ## 23 GO:0001002 RNA polymerase III type 1 promoter sequence-specifi...         1
    ## 24 GO:0003810 protein-glutamine gamma-glutamyltransferase activit...         1
    ## 25 GO:0050218                         propionate-CoA ligase activity         1
    ##    Significant Expected Fisher               type
    ## 1            4     0.60 0.0023 Biological.Process
    ## 2            8     3.10 0.0095 Biological.Process
    ## 3            1     0.03 0.0298 Biological.Process
    ## 4            1     0.03 0.0298 Biological.Process
    ## 5            1     0.03 0.0298 Biological.Process
    ## 6            2     0.17 0.0098 Molecular.Function
    ## 7            6     2.21 0.0216 Molecular.Function
    ## 8            2     0.29 0.0317 Molecular.Function
    ## 9            2     0.08 0.0414 Molecular.Function
    ## 10           1     0.04 0.0418 Molecular.Function
    ## 11           1     0.04 0.0418 Molecular.Function
    ## 12           1     0.04 0.0418 Molecular.Function
    ## 13           1     0.04 0.0418 Molecular.Function
    ## 14           1     0.04 0.0418 Molecular.Function
    ## 15           1     0.04 0.0418 Molecular.Function
    ## 16           1     0.04 0.0418 Molecular.Function
    ## 17           1     0.04 0.0418 Molecular.Function
    ## 18           1     0.04 0.0418 Molecular.Function
    ## 19           1     0.04 0.0418 Molecular.Function
    ## 20           1     0.04 0.0418 Molecular.Function
    ## 21           1     0.04 0.0418 Molecular.Function
    ## 22           1     0.04 0.0418 Molecular.Function
    ## 23           1     0.04 0.0418 Molecular.Function
    ## 24           1     0.04 0.0418 Molecular.Function
    ## 25           1     0.04 0.0418 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 2 1 2 1 2 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 314 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  3 nodes to be scored    (11 eliminated genes)

    ## 
    ##   Level 11:  7 nodes to be scored    (155 eliminated genes)

    ## 
    ##   Level 10:  13 nodes to be scored   (203 eliminated genes)

    ## 
    ##   Level 9:   20 nodes to be scored   (213 eliminated genes)

    ## 
    ##   Level 8:   23 nodes to be scored   (269 eliminated genes)

    ## 
    ##   Level 7:   38 nodes to be scored   (343 eliminated genes)

    ## 
    ##   Level 6:   62 nodes to be scored   (531 eliminated genes)

    ## 
    ##   Level 5:   65 nodes to be scored   (746 eliminated genes)

    ## 
    ##   Level 4:   38 nodes to be scored   (934 eliminated genes)

    ## 
    ##   Level 3:   29 nodes to be scored   (1180 eliminated genes)

    ## 
    ##   Level 2:   11 nodes to be scored   (1298 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1384 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 234 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   12 nodes to be scored   (10 eliminated genes)

    ## 
    ##   Level 8:   20 nodes to be scored   (50 eliminated genes)

    ## 
    ##   Level 7:   32 nodes to be scored   (469 eliminated genes)

    ## 
    ##   Level 6:   42 nodes to be scored   (586 eliminated genes)

    ## 
    ##   Level 5:   47 nodes to be scored   (874 eliminated genes)

    ## 
    ##   Level 4:   42 nodes to be scored   (1208 eliminated genes)

    ## 
    ##   Level 3:   24 nodes to be scored   (1993 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (2257 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2562 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0002218                   activation of innate immune response        15
    ## 2  GO:0000902                                     cell morphogenesis         7
    ## 3  GO:0001889                                      liver development         8
    ## 4  GO:0001832                                      blastocyst growth         1
    ## 5  GO:0006488 dolichol-linked oligosaccharide biosynthetic proces...         1
    ## 6  GO:0000076                   DNA replication checkpoint signaling         1
    ## 7  GO:0002230 positive regulation of defense response to virus by...         1
    ## 8  GO:0005245                 voltage-gated calcium channel activity         7
    ## 9  GO:0003964                   RNA-directed DNA polymerase activity       122
    ## 10 GO:0004104                                cholinesterase activity         1
    ## 11 GO:0001758                         retinal dehydrogenase activity         1
    ## 12 GO:0004736                          pyruvate carboxylase activity         1
    ## 13 GO:0003863 3-methyl-2-oxobutanoate dehydrogenase (2-methylprop...         1
    ## 14 GO:0004435          phosphatidylinositol phospholipase C activity         1
    ## 15 GO:0004844                      uracil DNA N-glycosylase activity         1
    ## 16 GO:0015101      organic cation transmembrane transporter activity         1
    ## 17 GO:0004494                      methylmalonyl-CoA mutase activity         1
    ## 18 GO:0005381            iron ion transmembrane transporter activity         1
    ## 19 GO:0004760                  serine-pyruvate transaminase activity         1
    ## 20 GO:0004683           calmodulin-dependent protein kinase activity         1
    ## 21 GO:0001002 RNA polymerase III type 1 promoter sequence-specifi...         1
    ## 22 GO:0005080                               protein kinase C binding         1
    ## 23 GO:0005324             long-chain fatty acid transporter activity         1
    ## 24 GO:0005246                     calcium channel regulator activity         1
    ## 25 GO:0000703 oxidized pyrimidine nucleobase lesion DNA N-glycosy...         1
    ## 26 GO:0005524                                            ATP binding       263
    ##    Significant Expected Fisher               type
    ## 1            3     0.55 0.0160 Biological.Process
    ## 2            2     0.26 0.0250 Biological.Process
    ## 3            2     0.29 0.0320 Biological.Process
    ## 4            1     0.04 0.0370 Biological.Process
    ## 5            1     0.04 0.0370 Biological.Process
    ## 6            1     0.04 0.0370 Biological.Process
    ## 7            1     0.04 0.0370 Biological.Process
    ## 8            3     0.30 0.0023 Molecular.Function
    ## 9           12     5.14 0.0045 Molecular.Function
    ## 10           1     0.04 0.0421 Molecular.Function
    ## 11           1     0.04 0.0421 Molecular.Function
    ## 12           1     0.04 0.0421 Molecular.Function
    ## 13           1     0.04 0.0421 Molecular.Function
    ## 14           1     0.04 0.0421 Molecular.Function
    ## 15           1     0.04 0.0421 Molecular.Function
    ## 16           1     0.04 0.0421 Molecular.Function
    ## 17           1     0.04 0.0421 Molecular.Function
    ## 18           1     0.04 0.0421 Molecular.Function
    ## 19           1     0.04 0.0421 Molecular.Function
    ## 20           1     0.04 0.0421 Molecular.Function
    ## 21           1     0.04 0.0421 Molecular.Function
    ## 22           1     0.04 0.0421 Molecular.Function
    ## 23           1     0.04 0.0421 Molecular.Function
    ## 24           1     0.04 0.0421 Molecular.Function
    ## 25           1     0.04 0.0421 Molecular.Function
    ## 26          17    11.08 0.0463 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 206 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  8 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 10:  11 nodes to be scored   (142 eliminated genes)

    ## 
    ##   Level 9:   16 nodes to be scored   (158 eliminated genes)

    ## 
    ##   Level 8:   17 nodes to be scored   (196 eliminated genes)

    ## 
    ##   Level 7:   24 nodes to be scored   (225 eliminated genes)

    ## 
    ##   Level 6:   30 nodes to be scored   (277 eliminated genes)

    ## 
    ##   Level 5:   40 nodes to be scored   (404 eliminated genes)

    ## 
    ##   Level 4:   29 nodes to be scored   (632 eliminated genes)

    ## 
    ##   Level 3:   19 nodes to be scored   (994 eliminated genes)

    ## 
    ##   Level 2:   7 nodes to be scored    (1130 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1258 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 134 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   6 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   8 nodes to be scored    (107 eliminated genes)

    ## 
    ##   Level 7:   12 nodes to be scored   (453 eliminated genes)

    ## 
    ##   Level 6:   24 nodes to be scored   (546 eliminated genes)

    ## 
    ##   Level 5:   29 nodes to be scored   (732 eliminated genes)

    ## 
    ##   Level 4:   28 nodes to be scored   (1043 eliminated genes)

    ## 
    ##   Level 3:   19 nodes to be scored   (1754 eliminated genes)

    ## 
    ##   Level 2:   5 nodes to be scored    (1982 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2525 eliminated genes)

    ##         GO.ID                                          Term Annotated
    ## 1  GO:0006527                    arginine catabolic process         1
    ## 2  GO:0000492                       box C/D snoRNP assembly         1
    ## 3  GO:0001825                          blastocyst formation         1
    ## 4  GO:0000076          DNA replication checkpoint signaling         1
    ## 5  GO:0001516            prostaglandin biosynthetic process         2
    ## 6  GO:0000050                                    urea cycle         2
    ## 7  GO:0000028              ribosomal small subunit assembly         3
    ## 8  GO:0003179                     heart valve morphogenesis         4
    ## 9  GO:0001708                       cell fate specification         4
    ## 10 GO:0004408            holocytochrome-c synthase activity         1
    ## 11 GO:0004516 nicotinate phosphoribosyltransferase activity         1
    ## 12 GO:0004129                 cytochrome-c oxidase activity         1
    ## 13 GO:0005381   iron ion transmembrane transporter activity         1
    ## 14 GO:0097602                 cullin family protein binding         1
    ## 15 GO:0003995               acyl-CoA dehydrogenase activity         1
    ## 16 GO:0004222                 metalloendopeptidase activity        15
    ## 17 GO:0004658            propionyl-CoA carboxylase activity         2
    ## 18 GO:0004474                      malate synthase activity         2
    ## 19 GO:0003985       acetyl-CoA C-acetyltransferase activity         2
    ## 20 GO:0004720             protein-lysine 6-oxidase activity         2
    ## 21 GO:0005319                    lipid transporter activity         3
    ##    Significant Expected Fisher               type
    ## 1            1     0.01  0.011 Biological.Process
    ## 2            1     0.01  0.011 Biological.Process
    ## 3            1     0.01  0.011 Biological.Process
    ## 4            1     0.01  0.011 Biological.Process
    ## 5            1     0.02  0.023 Biological.Process
    ## 6            1     0.02  0.023 Biological.Process
    ## 7            1     0.03  0.034 Biological.Process
    ## 8            1     0.05  0.045 Biological.Process
    ## 9            1     0.05  0.045 Biological.Process
    ## 10           1     0.02  0.016 Molecular.Function
    ## 11           1     0.02  0.016 Molecular.Function
    ## 12           1     0.02  0.016 Molecular.Function
    ## 13           1     0.02  0.016 Molecular.Function
    ## 14           1     0.02  0.016 Molecular.Function
    ## 15           1     0.02  0.016 Molecular.Function
    ## 16           2     0.24  0.024 Molecular.Function
    ## 17           1     0.03  0.032 Molecular.Function
    ## 18           1     0.03  0.032 Molecular.Function
    ## 19           1     0.03  0.032 Molecular.Function
    ## 20           1     0.03  0.032 Molecular.Function
    ## 21           1     0.05  0.048 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 247 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  5 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  7 nodes to be scored    (159 eliminated genes)

    ## 
    ##   Level 10:  10 nodes to be scored   (211 eliminated genes)

    ## 
    ##   Level 9:   17 nodes to be scored   (220 eliminated genes)

    ## 
    ##   Level 8:   20 nodes to be scored   (263 eliminated genes)

    ## 
    ##   Level 7:   28 nodes to be scored   (303 eliminated genes)

    ## 
    ##   Level 6:   42 nodes to be scored   (462 eliminated genes)

    ## 
    ##   Level 5:   51 nodes to be scored   (680 eliminated genes)

    ## 
    ##   Level 4:   32 nodes to be scored   (865 eliminated genes)

    ## 
    ##   Level 3:   21 nodes to be scored   (1047 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (1117 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1220 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 125 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   6 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   8 nodes to be scored    (142 eliminated genes)

    ## 
    ##   Level 7:   11 nodes to be scored   (455 eliminated genes)

    ## 
    ##   Level 6:   23 nodes to be scored   (561 eliminated genes)

    ## 
    ##   Level 5:   23 nodes to be scored   (763 eliminated genes)

    ## 
    ##   Level 4:   26 nodes to be scored   (1045 eliminated genes)

    ## 
    ##   Level 3:   17 nodes to be scored   (1681 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (2048 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2467 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0000077                        DNA damage checkpoint signaling         7
    ## 2  GO:0000972 transcription-dependent tethering of RNA polymerase...         1
    ## 3  GO:0003352                          regulation of cilium movement         1
    ## 4  GO:0000712      resolution of meiotic recombination intermediates         1
    ## 5  GO:0003964                   RNA-directed DNA polymerase activity       122
    ## 6  GO:0004609              phosphatidylserine decarboxylase activity        13
    ## 7  GO:0001716                          L-amino-acid oxidase activity         2
    ## 8  GO:0004658                     propionyl-CoA carboxylase activity         2
    ## 9  GO:0003987                            acetate-CoA ligase activity         2
    ## 10 GO:0003676                                   nucleic acid binding       497
    ##    Significant Expected  Fisher               type
    ## 1            2     0.15 0.00920 Biological.Process
    ## 2            1     0.02 0.02200 Biological.Process
    ## 3            1     0.02 0.02200 Biological.Process
    ## 4            1     0.02 0.02200 Biological.Process
    ## 5           10     2.55 0.00016 Molecular.Function
    ## 6            3     0.27 0.00213 Molecular.Function
    ## 7            1     0.04 0.04136 Molecular.Function
    ## 8            1     0.04 0.04136 Molecular.Function
    ## 9            1     0.04 0.04136 Molecular.Function
    ## 10          17    10.38 0.04918 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 345 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  6 nodes to be scored    (11 eliminated genes)

    ## 
    ##   Level 11:  11 nodes to be scored   (156 eliminated genes)

    ## 
    ##   Level 10:  17 nodes to be scored   (213 eliminated genes)

    ## 
    ##   Level 9:   25 nodes to be scored   (232 eliminated genes)

    ## 
    ##   Level 8:   29 nodes to be scored   (326 eliminated genes)

    ## 
    ##   Level 7:   38 nodes to be scored   (419 eliminated genes)

    ## 
    ##   Level 6:   58 nodes to be scored   (504 eliminated genes)

    ## 
    ##   Level 5:   70 nodes to be scored   (742 eliminated genes)

    ## 
    ##   Level 4:   41 nodes to be scored   (920 eliminated genes)

    ## 
    ##   Level 3:   32 nodes to be scored   (1190 eliminated genes)

    ## 
    ##   Level 2:   12 nodes to be scored   (1263 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1362 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 216 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   10 nodes to be scored   (1 eliminated genes)

    ## 
    ##   Level 8:   14 nodes to be scored   (150 eliminated genes)

    ## 
    ##   Level 7:   28 nodes to be scored   (462 eliminated genes)

    ## 
    ##   Level 6:   43 nodes to be scored   (575 eliminated genes)

    ## 
    ##   Level 5:   42 nodes to be scored   (884 eliminated genes)

    ## 
    ##   Level 4:   38 nodes to be scored   (1225 eliminated genes)

    ## 
    ##   Level 3:   23 nodes to be scored   (1976 eliminated genes)

    ## 
    ##   Level 2:   12 nodes to be scored   (2274 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2597 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0006813                                potassium ion transport         1
    ## 2  GO:0001825                                   blastocyst formation         1
    ## 3  GO:0001951                        intestinal D-glucose absorption         1
    ## 4  GO:0002265       astrocyte activation involved in immune response         1
    ## 5  GO:0000076                   DNA replication checkpoint signaling         1
    ## 6  GO:0000460                                maturation of 5.8S rRNA         1
    ## 7  GO:0016311                                      dephosphorylation         1
    ## 8  GO:0004560                            alpha-L-fucosidase activity         6
    ## 9  GO:0004325                                ferrochelatase activity         1
    ## 10 GO:0002161                        aminoacyl-tRNA editing activity         1
    ## 11 GO:0004777 succinate-semialdehyde dehydrogenase (NAD+) activit...         1
    ## 12 GO:0003980 UDP-glucose:glycoprotein glucosyltransferase activi...         1
    ## 13 GO:0004850                         uridine phosphorylase activity         1
    ## 14 GO:0001758                         retinal dehydrogenase activity         1
    ## 15 GO:0003923                       GPI-anchor transamidase activity         1
    ## 16 GO:0005381            iron ion transmembrane transporter activity         1
    ## 17 GO:0004134                    4-alpha-glucanotransferase activity         1
    ## 18 GO:0004067                                  asparaginase activity         1
    ## 19 GO:0001002 RNA polymerase III type 1 promoter sequence-specifi...         1
    ## 20 GO:0003960                       NADPH:quinone reductase activity         1
    ##    Significant Expected Fisher               type
    ## 1            1     0.03  0.034 Biological.Process
    ## 2            1     0.03  0.034 Biological.Process
    ## 3            1     0.03  0.034 Biological.Process
    ## 4            1     0.03  0.034 Biological.Process
    ## 5            1     0.03  0.034 Biological.Process
    ## 6            1     0.03  0.034 Biological.Process
    ## 7            1     0.03  0.034 Biological.Process
    ## 8            2     0.23  0.020 Molecular.Function
    ## 9            1     0.04  0.039 Molecular.Function
    ## 10           1     0.04  0.039 Molecular.Function
    ## 11           1     0.04  0.039 Molecular.Function
    ## 12           1     0.04  0.039 Molecular.Function
    ## 13           1     0.04  0.039 Molecular.Function
    ## 14           1     0.04  0.039 Molecular.Function
    ## 15           1     0.04  0.039 Molecular.Function
    ## 16           1     0.04  0.039 Molecular.Function
    ## 17           1     0.04  0.039 Molecular.Function
    ## 18           1     0.04  0.039 Molecular.Function
    ## 19           1     0.04  0.039 Molecular.Function
    ## 20           1     0.04  0.039 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 206 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  5 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 10:  6 nodes to be scored    (203 eliminated genes)

    ## 
    ##   Level 9:   12 nodes to be scored   (206 eliminated genes)

    ## 
    ##   Level 8:   14 nodes to be scored   (243 eliminated genes)

    ## 
    ##   Level 7:   23 nodes to be scored   (287 eliminated genes)

    ## 
    ##   Level 6:   37 nodes to be scored   (355 eliminated genes)

    ## 
    ##   Level 5:   40 nodes to be scored   (516 eliminated genes)

    ## 
    ##   Level 4:   33 nodes to be scored   (717 eliminated genes)

    ## 
    ##   Level 3:   22 nodes to be scored   (1032 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (1191 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1307 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 156 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   5 nodes to be scored    (100 eliminated genes)

    ## 
    ##   Level 7:   11 nodes to be scored   (389 eliminated genes)

    ## 
    ##   Level 6:   35 nodes to be scored   (407 eliminated genes)

    ## 
    ##   Level 5:   38 nodes to be scored   (628 eliminated genes)

    ## 
    ##   Level 4:   37 nodes to be scored   (976 eliminated genes)

    ## 
    ##   Level 3:   18 nodes to be scored   (1496 eliminated genes)

    ## 
    ##   Level 2:   6 nodes to be scored    (1831 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2301 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0019700                  organic phosphonate catabolic process        12
    ## 2  GO:0002933                                    lipid hydroxylation         3
    ## 3  GO:0001731         formation of translation preinitiation complex         1
    ## 4  GO:0000972 transcription-dependent tethering of RNA polymerase...         1
    ## 5  GO:0000398                         mRNA splicing, via spliceosome        63
    ## 6  GO:0000050                                             urea cycle         2
    ## 7  GO:0002790                                      peptide secretion         2
    ## 8  GO:0003746                 translation elongation factor activity         7
    ## 9  GO:0004342             glucosamine-6-phosphate deaminase activity         2
    ## 10 GO:0004142      diacylglycerol cholinephosphotransferase activity         1
    ## 11 GO:0004325                                ferrochelatase activity         1
    ## 12 GO:0004066 asparagine synthase (glutamine-hydrolyzing) activit...         1
    ## 13 GO:0004067                                  asparaginase activity         1
    ## 14 GO:0008137               NADH dehydrogenase (ubiquinone) activity         1
    ## 15 GO:0047830                      D-octopine dehydrogenase activity         1
    ## 16 GO:0002161                        aminoacyl-tRNA editing activity         1
    ## 17 GO:0004152                  dihydroorotate dehydrogenase activity         1
    ## 18 GO:0004334                           fumarylacetoacetase activity         1
    ## 19 GO:0046982                    protein heterodimerization activity         1
    ## 20 GO:0003960                       NADPH:quinone reductase activity         1
    ## 21 GO:0004458          D-lactate dehydrogenase (cytochrome) activity         1
    ## 22 GO:0004174 electron-transferring-flavoprotein dehydrogenase ac...         1
    ## 23 GO:0097602                          cullin family protein binding         1
    ## 24 GO:0000049                                           tRNA binding        30
    ## 25 GO:0004089                         carbonate dehydratase activity         2
    ##    Significant Expected    Fisher               type
    ## 1            5     0.26 0.0000022 Biological.Process
    ## 2            2     0.06 0.0013000 Biological.Process
    ## 3            1     0.02 0.0213000 Biological.Process
    ## 4            1     0.02 0.0213000 Biological.Process
    ## 5            4     1.34 0.0413000 Biological.Process
    ## 6            1     0.04 0.0421000 Biological.Process
    ## 7            1     0.04 0.0421000 Biological.Process
    ## 8            3     0.17 0.0004400 Molecular.Function
    ## 9            2     0.05 0.0005700 Molecular.Function
    ## 10           1     0.02 0.0241400 Molecular.Function
    ## 11           1     0.02 0.0241400 Molecular.Function
    ## 12           1     0.02 0.0241400 Molecular.Function
    ## 13           1     0.02 0.0241400 Molecular.Function
    ## 14           1     0.02 0.0241400 Molecular.Function
    ## 15           1     0.02 0.0241400 Molecular.Function
    ## 16           1     0.02 0.0241400 Molecular.Function
    ## 17           1     0.02 0.0241400 Molecular.Function
    ## 18           1     0.02 0.0241400 Molecular.Function
    ## 19           1     0.02 0.0241400 Molecular.Function
    ## 20           1     0.02 0.0241400 Molecular.Function
    ## 21           1     0.02 0.0241400 Molecular.Function
    ## 22           1     0.02 0.0241400 Molecular.Function
    ## 23           1     0.02 0.0241400 Molecular.Function
    ## 24           3     0.72 0.0342800 Molecular.Function
    ## 25           1     0.05 0.0477000 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 193 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  2 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 12:  6 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 11:  10 nodes to be scored   (144 eliminated genes)

    ## 
    ##   Level 10:  11 nodes to be scored   (209 eliminated genes)

    ## 
    ##   Level 9:   15 nodes to be scored   (223 eliminated genes)

    ## 
    ##   Level 8:   13 nodes to be scored   (252 eliminated genes)

    ## 
    ##   Level 7:   19 nodes to be scored   (258 eliminated genes)

    ## 
    ##   Level 6:   24 nodes to be scored   (379 eliminated genes)

    ## 
    ##   Level 5:   35 nodes to be scored   (555 eliminated genes)

    ## 
    ##   Level 4:   28 nodes to be scored   (627 eliminated genes)

    ## 
    ##   Level 3:   19 nodes to be scored   (838 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (1146 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1293 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 96 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   3 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 7:   10 nodes to be scored   (269 eliminated genes)

    ## 
    ##   Level 6:   15 nodes to be scored   (343 eliminated genes)

    ## 
    ##   Level 5:   19 nodes to be scored   (650 eliminated genes)

    ## 
    ##   Level 4:   21 nodes to be scored   (935 eliminated genes)

    ## 
    ##   Level 3:   16 nodes to be scored   (1340 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (1976 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2486 eliminated genes)

    ##        GO.ID                                    Term Annotated Significant
    ## 1 GO:0002933                     lipid hydroxylation         3           2
    ## 2 GO:0006836              neurotransmitter transport         6           2
    ## 3 GO:0006397                         mRNA processing        64           2
    ## 4 GO:0001503                            ossification        42           3
    ## 5 GO:0001516      prostaglandin biosynthetic process         2           1
    ## 6 GO:0003964    RNA-directed DNA polymerase activity       122           7
    ## 7 GO:0005436     sodium:phosphate symporter activity         1           1
    ## 8 GO:0004508 steroid 17-alpha-monooxygenase activity        17           2
    ## 9 GO:0004565             beta-galactosidase activity         3           1
    ##   Expected  Fisher               type
    ## 1     0.05 0.00082 Biological.Process
    ## 2     0.10 0.00399 Biological.Process
    ## 3     1.09 0.01706 Biological.Process
    ## 4     0.71 0.03213 Biological.Process
    ## 5     0.03 0.03374 Biological.Process
    ## 6     1.71 0.00130 Molecular.Function
    ## 7     0.01 0.01400 Molecular.Function
    ## 8     0.24 0.02290 Molecular.Function
    ## 9     0.04 0.04160 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 161 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 9:   8 nodes to be scored    (11 eliminated genes)

    ## 
    ##   Level 8:   11 nodes to be scored   (28 eliminated genes)

    ## 
    ##   Level 7:   18 nodes to be scored   (33 eliminated genes)

    ## 
    ##   Level 6:   27 nodes to be scored   (175 eliminated genes)

    ## 
    ##   Level 5:   34 nodes to be scored   (587 eliminated genes)

    ## 
    ##   Level 4:   26 nodes to be scored   (661 eliminated genes)

    ## 
    ##   Level 3:   20 nodes to be scored   (953 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (1118 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1287 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 129 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   4 nodes to be scored    (100 eliminated genes)

    ## 
    ##   Level 7:   15 nodes to be scored   (407 eliminated genes)

    ## 
    ##   Level 6:   25 nodes to be scored   (491 eliminated genes)

    ## 
    ##   Level 5:   29 nodes to be scored   (828 eliminated genes)

    ## 
    ##   Level 4:   25 nodes to be scored   (1186 eliminated genes)

    ## 
    ##   Level 3:   19 nodes to be scored   (1816 eliminated genes)

    ## 
    ##   Level 2:   6 nodes to be scored    (2086 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2480 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0000724 double-strand break repair via homologous recombina...        13
    ## 2  GO:0000184 nuclear-transcribed mRNA catabolic process, nonsens...        10
    ## 3  GO:0000712      resolution of meiotic recombination intermediates         1
    ## 4  GO:0006813                                potassium ion transport         1
    ## 5  GO:0000423                                              mitophagy         2
    ## 6  GO:0000028                       ribosomal small subunit assembly         3
    ## 7  GO:0002933                                    lipid hydroxylation         3
    ## 8  GO:0005524                                            ATP binding       263
    ## 9  GO:0004108                         citrate (Si)-synthase activity         1
    ## 10 GO:0004013                        adenosylhomocysteinase activity         1
    ## 11 GO:0004185                  serine-type carboxypeptidase activity         2
    ## 12 GO:0004089                         carbonate dehydratase activity         2
    ## 13 GO:0004517                         nitric-oxide synthase activity         2
    ## 14 GO:0005201            extracellular matrix structural constituent        17
    ##    Significant Expected  Fisher               type
    ## 1            3     0.20 0.00085 Biological.Process
    ## 2            2     0.16 0.00969 Biological.Process
    ## 3            1     0.02 0.01559 Biological.Process
    ## 4            1     0.02 0.01559 Biological.Process
    ## 5            1     0.03 0.03095 Biological.Process
    ## 6            1     0.05 0.04608 Biological.Process
    ## 7            1     0.05 0.04608 Biological.Process
    ## 8           12     5.02 0.00310 Molecular.Function
    ## 9            1     0.02 0.01910 Molecular.Function
    ## 10           1     0.02 0.01910 Molecular.Function
    ## 11           1     0.04 0.03780 Molecular.Function
    ## 12           1     0.04 0.03780 Molecular.Function
    ## 13           1     0.04 0.03780 Molecular.Function
    ## 14           2     0.32 0.04050 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 406 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  5 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  5 nodes to be scored    (5 eliminated genes)

    ## 
    ##   Level 11:  10 nodes to be scored   (167 eliminated genes)

    ## 
    ##   Level 10:  18 nodes to be scored   (213 eliminated genes)

    ## 
    ##   Level 9:   33 nodes to be scored   (237 eliminated genes)

    ## 
    ##   Level 8:   40 nodes to be scored   (336 eliminated genes)

    ## 
    ##   Level 7:   55 nodes to be scored   (413 eliminated genes)

    ## 
    ##   Level 6:   68 nodes to be scored   (597 eliminated genes)

    ## 
    ##   Level 5:   78 nodes to be scored   (882 eliminated genes)

    ## 
    ##   Level 4:   50 nodes to be scored   (988 eliminated genes)

    ## 
    ##   Level 3:   32 nodes to be scored   (1215 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (1330 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1387 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 217 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   6 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   12 nodes to be scored   (142 eliminated genes)

    ## 
    ##   Level 7:   26 nodes to be scored   (455 eliminated genes)

    ## 
    ##   Level 6:   46 nodes to be scored   (573 eliminated genes)

    ## 
    ##   Level 5:   49 nodes to be scored   (898 eliminated genes)

    ## 
    ##   Level 4:   43 nodes to be scored   (1339 eliminated genes)

    ## 
    ##   Level 3:   23 nodes to be scored   (2022 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (2305 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2597 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0002291 T cell activation via T cell receptor contact with ...        10
    ## 2 GO:0002674     negative regulation of acute inflammatory response         4
    ## 3 GO:0002221         pattern recognition receptor signaling pathway        11
    ## 4 GO:0005524                                            ATP binding       263
    ## 5 GO:0003676                                   nucleic acid binding       497
    ## 6 GO:0003824                                     catalytic activity      1042
    ## 7 GO:0004622                             lysophospholipase activity         4
    ## 8 GO:0003964                   RNA-directed DNA polymerase activity       122
    ##   Significant Expected  Fisher               type
    ## 1           3     0.52 0.01300 Biological.Process
    ## 2           2     0.21 0.01500 Biological.Process
    ## 3           3     0.58 0.01700 Biological.Process
    ## 4          29    15.25 0.00037 Molecular.Function
    ## 5          34    28.82 0.00101 Molecular.Function
    ## 6          64    60.43 0.01103 Molecular.Function
    ## 7           2     0.23 0.01856 Molecular.Function
    ## 8          12     7.08 0.04743 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 194 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  3 nodes to be scored    (11 eliminated genes)

    ## 
    ##   Level 11:  5 nodes to be scored    (159 eliminated genes)

    ## 
    ##   Level 10:  8 nodes to be scored    (203 eliminated genes)

    ## 
    ##   Level 9:   12 nodes to be scored   (209 eliminated genes)

    ## 
    ##   Level 8:   16 nodes to be scored   (249 eliminated genes)

    ## 
    ##   Level 7:   18 nodes to be scored   (294 eliminated genes)

    ## 
    ##   Level 6:   28 nodes to be scored   (370 eliminated genes)

    ## 
    ##   Level 5:   38 nodes to be scored   (483 eliminated genes)

    ## 
    ##   Level 4:   30 nodes to be scored   (650 eliminated genes)

    ## 
    ##   Level 3:   21 nodes to be scored   (917 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (1084 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1228 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 145 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   5 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   7 nodes to be scored    (10 eliminated genes)

    ## 
    ##   Level 7:   18 nodes to be scored   (348 eliminated genes)

    ## 
    ##   Level 6:   27 nodes to be scored   (360 eliminated genes)

    ## 
    ##   Level 5:   29 nodes to be scored   (579 eliminated genes)

    ## 
    ##   Level 4:   27 nodes to be scored   (931 eliminated genes)

    ## 
    ##   Level 3:   22 nodes to be scored   (1698 eliminated genes)

    ## 
    ##   Level 2:   7 nodes to be scored    (1934 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2276 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0006598                            polyamine catabolic process         3
    ## 2  GO:0042060                                          wound healing         1
    ## 3  GO:0006813                                potassium ion transport         1
    ## 4  GO:0001867                  complement activation, lectin pathway        13
    ## 5  GO:0001578                           microtubule bundle formation         2
    ## 6  GO:0003985                acetyl-CoA C-acetyltransferase activity         2
    ## 7  GO:0003676                                   nucleic acid binding       497
    ## 8  GO:0003923                       GPI-anchor transamidase activity         1
    ## 9  GO:0004060                 arylamine N-acetyltransferase activity         1
    ## 10 GO:0003913                                DNA photolyase activity         1
    ## 11 GO:0003997                              acyl-CoA oxidase activity         1
    ## 12 GO:0003944 N-acetylglucosamine-1-phosphodiester alpha-N-acetyl...         1
    ## 13 GO:0004777 succinate-semialdehyde dehydrogenase (NAD+) activit...         1
    ## 14 GO:0004748 ribonucleoside-diphosphate reductase activity, thio...         1
    ## 15 GO:0004176                       ATP-dependent peptidase activity         2
    ##    Significant Expected  Fisher               type
    ## 1            2     0.05 0.00089 Biological.Process
    ## 2            1     0.02 0.01772 Biological.Process
    ## 3            1     0.02 0.01772 Biological.Process
    ## 4            2     0.23 0.02087 Biological.Process
    ## 5            1     0.04 0.03513 Biological.Process
    ## 6            2     0.05 0.00057 Molecular.Function
    ## 7           19    12.00 0.00264 Molecular.Function
    ## 8            1     0.02 0.02414 Molecular.Function
    ## 9            1     0.02 0.02414 Molecular.Function
    ## 10           1     0.02 0.02414 Molecular.Function
    ## 11           1     0.02 0.02414 Molecular.Function
    ## 12           1     0.02 0.02414 Molecular.Function
    ## 13           1     0.02 0.02414 Molecular.Function
    ## 14           1     0.02 0.02414 Molecular.Function
    ## 15           1     0.05 0.04770 Molecular.Function

``` r
head(results_all_targets)
```

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0000082                  G1/S transition of mitotic cell cycle       104
    ## 2 GO:0001172                            RNA-templated transcription         1
    ## 3 GO:0002230 positive regulation of defense response to virus by...         1
    ## 4 GO:0001913                           T cell mediated cytotoxicity         1
    ## 5 GO:0004176                       ATP-dependent peptidase activity         2
    ## 6 GO:0003964                   RNA-directed DNA polymerase activity       122
    ##   Significant Expected  Fisher               type         miRNA
    ## 1           9     2.65 0.00080 Biological.Process Cluster_10452
    ## 2           1     0.03 0.02550 Biological.Process Cluster_10452
    ## 3           1     0.03 0.02550 Biological.Process Cluster_10452
    ## 4           1     0.03 0.02550 Biological.Process Cluster_10452
    ## 5           2     0.06 0.00099 Molecular.Function Cluster_10452
    ## 6          11     3.87 0.00135 Molecular.Function Cluster_10452

Save results

``` r
write.csv(results_all_targets, "../output/27-Apul-mRNA-miRNA-interactions-topGO/miRNA_all_targets_topGO_FE.csv")
```

# 4 FE of specific miRNA’s targets (high 0.5 cor targets)

Loop through all miRNA and run functional enrichment on the miRNA’s
highly correlated targets (PCC magnitude \> 0.5)

``` r
interacting_miRNAs_high0.5 <- unique(high0.5_cor_bind_FA$miRNA) %>% na.omit
results_high0.5_cor_targets <- NULL  # initialize empty df

for(miRNA in interacting_miRNAs_high0.5) {
  
  # Run topGO enrichment function
  miRNA_results <- miRNA_topGO_FE(miRNA, high0.5_cor_bind_FA)
  
  # Only keep results if not empty
  if (nrow(miRNA_results) > 0) {
    
    # Add the miRNA source column
    miRNA_results$miRNA <- miRNA

    # Bind to the accumulating results data frame
    results_high0.5_cor_targets <- rbind(results_high0.5_cor_targets, miRNA_results)
  }
}
```

    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 34 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   2 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 7:   3 nodes to be scored    (6 eliminated genes)

    ## 
    ##   Level 6:   4 nodes to be scored    (14 eliminated genes)

    ## 
    ##   Level 5:   5 nodes to be scored    (16 eliminated genes)

    ## 
    ##   Level 4:   6 nodes to be scored    (16 eliminated genes)

    ## 
    ##   Level 3:   5 nodes to be scored    (211 eliminated genes)

    ## 
    ##   Level 2:   4 nodes to be scored    (672 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (778 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 38 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 8:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   5 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   5 nodes to be scored    (12 eliminated genes)

    ## 
    ##   Level 5:   7 nodes to be scored    (160 eliminated genes)

    ## 
    ##   Level 4:   8 nodes to be scored    (166 eliminated genes)

    ## 
    ##   Level 3:   7 nodes to be scored    (472 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (1024 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1779 eliminated genes)

    ##        GO.ID                                               Term Annotated
    ## 1 GO:0006598                        polyamine catabolic process         3
    ## 2 GO:0001745                         compound eye morphogenesis         3
    ## 3 GO:0000014 single-stranded DNA endodeoxyribonuclease activity         3
    ## 4 GO:0005290     L-histidine transmembrane transporter activity         9
    ##   Significant Expected Fisher               type
    ## 1           1     0.00 0.0042 Biological.Process
    ## 2           1     0.00 0.0042 Biological.Process
    ## 3           1     0.01 0.0065 Molecular.Function
    ## 4           1     0.02 0.0193 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 42 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  1 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 12:  2 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 11:  3 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 10:  3 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 9:   5 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 8:   4 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 7:   4 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 6:   4 nodes to be scored    (72 eliminated genes)

    ## 
    ##   Level 5:   6 nodes to be scored    (125 eliminated genes)

    ## 
    ##   Level 4:   4 nodes to be scored    (138 eliminated genes)

    ## 
    ##   Level 3:   1 nodes to be scored    (209 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (239 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (246 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 57 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   5 nodes to be scored    (371 eliminated genes)

    ## 
    ##   Level 6:   9 nodes to be scored    (409 eliminated genes)

    ## 
    ##   Level 5:   12 nodes to be scored   (561 eliminated genes)

    ## 
    ##   Level 4:   10 nodes to be scored   (640 eliminated genes)

    ## 
    ##   Level 3:   11 nodes to be scored   (1150 eliminated genes)

    ## 
    ##   Level 2:   4 nodes to be scored    (1518 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2039 eliminated genes)

    ##        GO.ID                                             Term Annotated
    ## 1 GO:0001658 branching involved in ureteric bud morphogenesis         4
    ## 2 GO:0001945                         lymph vessel development         5
    ## 3 GO:0004252               serine-type endopeptidase activity        53
    ## 4 GO:0005201      extracellular matrix structural constituent        17
    ## 5 GO:0004771                         sterol esterase activity         5
    ##   Significant Expected Fisher               type
    ## 1           1     0.01 0.0057 Biological.Process
    ## 2           1     0.01 0.0071 Biological.Process
    ## 3           3     0.25 0.0016 Molecular.Function
    ## 4           2     0.08 0.0026 Molecular.Function
    ## 5           1     0.02 0.0232 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 70 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   2 nodes to be scored    (12 eliminated genes)

    ## 
    ##   Level 8:   3 nodes to be scored    (52 eliminated genes)

    ## 
    ##   Level 7:   7 nodes to be scored    (55 eliminated genes)

    ## 
    ##   Level 6:   10 nodes to be scored   (122 eliminated genes)

    ## 
    ##   Level 5:   14 nodes to be scored   (178 eliminated genes)

    ## 
    ##   Level 4:   11 nodes to be scored   (211 eliminated genes)

    ## 
    ##   Level 3:   11 nodes to be scored   (339 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (586 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (765 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 37 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   3 nodes to be scored    (263 eliminated genes)

    ## 
    ##   Level 6:   6 nodes to be scored    (267 eliminated genes)

    ## 
    ##   Level 5:   6 nodes to be scored    (397 eliminated genes)

    ## 
    ##   Level 4:   8 nodes to be scored    (522 eliminated genes)

    ## 
    ##   Level 3:   9 nodes to be scored    (611 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (1263 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1841 eliminated genes)

    ##        GO.ID                                 Term Annotated Significant
    ## 1 GO:0001895                   retina homeostasis         2           1
    ## 2 GO:0004013      adenosylhomocysteinase activity         1           1
    ## 3 GO:0005524                          ATP binding       263           5
    ## 4 GO:0009034               tryptophanase activity         2           1
    ## 5 GO:0003964 RNA-directed DNA polymerase activity       122           3
    ##   Expected Fisher               type
    ## 1     0.01 0.0085 Biological.Process
    ## 2     0.00 0.0047 Molecular.Function
    ## 3     1.23 0.0050 Molecular.Function
    ## 4     0.01 0.0093 Molecular.Function
    ## 5     0.57 0.0171 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 73 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   2 nodes to be scored    (12 eliminated genes)

    ## 
    ##   Level 8:   3 nodes to be scored    (52 eliminated genes)

    ## 
    ##   Level 7:   7 nodes to be scored    (55 eliminated genes)

    ## 
    ##   Level 6:   10 nodes to be scored   (122 eliminated genes)

    ## 
    ##   Level 5:   15 nodes to be scored   (178 eliminated genes)

    ## 
    ##   Level 4:   12 nodes to be scored   (211 eliminated genes)

    ## 
    ##   Level 3:   12 nodes to be scored   (359 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (605 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (827 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 36 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   3 nodes to be scored    (263 eliminated genes)

    ## 
    ##   Level 6:   6 nodes to be scored    (267 eliminated genes)

    ## 
    ##   Level 5:   6 nodes to be scored    (397 eliminated genes)

    ## 
    ##   Level 4:   7 nodes to be scored    (524 eliminated genes)

    ## 
    ##   Level 3:   9 nodes to be scored    (613 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (768 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1594 eliminated genes)

    ##        GO.ID                                 Term Annotated Significant
    ## 1 GO:0001895                   retina homeostasis         2           1
    ## 2 GO:0005524                          ATP binding       263           6
    ## 3 GO:0003964 RNA-directed DNA polymerase activity       122           4
    ## 4 GO:0009034               tryptophanase activity         2           1
    ## 5 GO:0004322                 ferroxidase activity         3           1
    ## 6 GO:0000287                magnesium ion binding        56           2
    ##   Expected Fisher               type
    ## 1     0.01 0.0099 Biological.Process
    ## 2     1.42 0.0016 Molecular.Function
    ## 3     0.66 0.0033 Molecular.Function
    ## 4     0.01 0.0108 Molecular.Function
    ## 5     0.02 0.0161 Molecular.Function
    ## 6     0.30 0.0355 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 70 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   2 nodes to be scored    (12 eliminated genes)

    ## 
    ##   Level 8:   3 nodes to be scored    (52 eliminated genes)

    ## 
    ##   Level 7:   7 nodes to be scored    (55 eliminated genes)

    ## 
    ##   Level 6:   10 nodes to be scored   (122 eliminated genes)

    ## 
    ##   Level 5:   14 nodes to be scored   (178 eliminated genes)

    ## 
    ##   Level 4:   11 nodes to be scored   (211 eliminated genes)

    ## 
    ##   Level 3:   11 nodes to be scored   (339 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (586 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (765 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 36 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   3 nodes to be scored    (263 eliminated genes)

    ## 
    ##   Level 6:   5 nodes to be scored    (267 eliminated genes)

    ## 
    ##   Level 5:   6 nodes to be scored    (397 eliminated genes)

    ## 
    ##   Level 4:   8 nodes to be scored    (458 eliminated genes)

    ## 
    ##   Level 3:   9 nodes to be scored    (613 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (1265 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1594 eliminated genes)

    ##        GO.ID                                 Term Annotated Significant
    ## 1 GO:0001895                   retina homeostasis         2           1
    ## 2 GO:0009034               tryptophanase activity         2           1
    ## 3 GO:0004322                 ferroxidase activity         3           1
    ## 4 GO:0003964 RNA-directed DNA polymerase activity       122           3
    ## 5 GO:0005524                          ATP binding       263           4
    ## 6 GO:0000287                magnesium ion binding        56           2
    ##   Expected Fisher               type
    ## 1     0.01 0.0085 Biological.Process
    ## 2     0.01 0.0086 Molecular.Function
    ## 3     0.01 0.0129 Molecular.Function
    ## 4     0.53 0.0136 Molecular.Function
    ## 5     1.14 0.0211 Molecular.Function
    ## 6     0.24 0.0232 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 24 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   1 nodes to be scored    (12 eliminated genes)

    ## 
    ##   Level 8:   3 nodes to be scored    (37 eliminated genes)

    ## 
    ##   Level 7:   3 nodes to be scored    (40 eliminated genes)

    ## 
    ##   Level 6:   4 nodes to be scored    (122 eliminated genes)

    ## 
    ##   Level 5:   4 nodes to be scored    (137 eliminated genes)

    ## 
    ##   Level 4:   3 nodes to be scored    (148 eliminated genes)

    ## 
    ##   Level 3:   1 nodes to be scored    (217 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (235 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (246 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 40 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   4 nodes to be scored    (263 eliminated genes)

    ## 
    ##   Level 6:   5 nodes to be scored    (270 eliminated genes)

    ## 
    ##   Level 5:   8 nodes to be scored    (400 eliminated genes)

    ## 
    ##   Level 4:   8 nodes to be scored    (534 eliminated genes)

    ## 
    ##   Level 3:   9 nodes to be scored    (1182 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (1465 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1794 eliminated genes)

    ##        GO.ID                                 Term Annotated Significant
    ## 1 GO:0002040               sprouting angiogenesis        12           1
    ## 2 GO:0001701       in utero embryonic development        15           1
    ## 3 GO:0004382             GDP phosphatase activity         3           1
    ## 4 GO:0003964 RNA-directed DNA polymerase activity       122           3
    ##   Expected Fisher               type
    ## 1     0.03  0.025 Biological.Process
    ## 2     0.03  0.032 Biological.Process
    ## 3     0.01  0.013 Molecular.Function
    ## 4     0.53  0.014 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 23 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   2 nodes to be scored    (7 eliminated genes)

    ## 
    ##   Level 6:   3 nodes to be scored    (7 eliminated genes)

    ## 
    ##   Level 5:   4 nodes to be scored    (7 eliminated genes)

    ## 
    ##   Level 4:   4 nodes to be scored    (29 eliminated genes)

    ## 
    ##   Level 3:   5 nodes to be scored    (181 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (471 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (532 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ##        GO.ID                                          Term Annotated
    ## 1 GO:0006511 ubiquitin-dependent protein catabolic process         7
    ##   Significant Expected Fisher               type
    ## 1           1        0  0.005 Biological.Process
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 115 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  3 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 9:   7 nodes to be scored    (142 eliminated genes)

    ## 
    ##   Level 8:   7 nodes to be scored    (154 eliminated genes)

    ## 
    ##   Level 7:   9 nodes to be scored    (161 eliminated genes)

    ## 
    ##   Level 6:   15 nodes to be scored   (254 eliminated genes)

    ## 
    ##   Level 5:   22 nodes to be scored   (378 eliminated genes)

    ## 
    ##   Level 4:   22 nodes to be scored   (437 eliminated genes)

    ## 
    ##   Level 3:   16 nodes to be scored   (830 eliminated genes)

    ## 
    ##   Level 2:   6 nodes to be scored    (1032 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1257 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 93 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   7 nodes to be scored    (263 eliminated genes)

    ## 
    ##   Level 6:   17 nodes to be scored   (267 eliminated genes)

    ## 
    ##   Level 5:   19 nodes to be scored   (540 eliminated genes)

    ## 
    ##   Level 4:   22 nodes to be scored   (797 eliminated genes)

    ## 
    ##   Level 3:   17 nodes to be scored   (1386 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (1734 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2475 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0006598                            polyamine catabolic process         3
    ## 2  GO:0001676                long-chain fatty acid metabolic process         3
    ## 3  GO:0001822                                     kidney development        65
    ## 4  GO:0003777                             microtubule motor activity        16
    ## 5  GO:0004591 oxoglutarate dehydrogenase (succinyl-transferring) ...         1
    ## 6  GO:0004494                      methylmalonyl-CoA mutase activity         1
    ## 7  GO:0005080                               protein kinase C binding         1
    ## 8  GO:0004190                   aspartic-type endopeptidase activity         2
    ## 9  GO:0005319                             lipid transporter activity         3
    ## 10 GO:0000016                                       lactase activity         4
    ##    Significant Expected Fisher               type
    ## 1            1     0.01 0.0130 Biological.Process
    ## 2            1     0.01 0.0130 Biological.Process
    ## 3            2     0.28 0.0280 Biological.Process
    ## 4            2     0.14 0.0087 Molecular.Function
    ## 5            1     0.01 0.0090 Molecular.Function
    ## 6            1     0.01 0.0090 Molecular.Function
    ## 7            1     0.01 0.0090 Molecular.Function
    ## 8            1     0.02 0.0179 Molecular.Function
    ## 9            1     0.03 0.0268 Molecular.Function
    ## 10           1     0.04 0.0356 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 115 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  3 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 11:  4 nodes to be scored    (142 eliminated genes)

    ## 
    ##   Level 10:  5 nodes to be scored    (203 eliminated genes)

    ## 
    ##   Level 9:   6 nodes to be scored    (205 eliminated genes)

    ## 
    ##   Level 8:   6 nodes to be scored    (207 eliminated genes)

    ## 
    ##   Level 7:   8 nodes to be scored    (211 eliminated genes)

    ## 
    ##   Level 6:   17 nodes to be scored   (257 eliminated genes)

    ## 
    ##   Level 5:   25 nodes to be scored   (285 eliminated genes)

    ## 
    ##   Level 4:   18 nodes to be scored   (369 eliminated genes)

    ## 
    ##   Level 3:   12 nodes to be scored   (527 eliminated genes)

    ## 
    ##   Level 2:   7 nodes to be scored    (962 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1103 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 113 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   6 nodes to be scored    (100 eliminated genes)

    ## 
    ##   Level 7:   17 nodes to be scored   (396 eliminated genes)

    ## 
    ##   Level 6:   22 nodes to be scored   (421 eliminated genes)

    ## 
    ##   Level 5:   25 nodes to be scored   (544 eliminated genes)

    ## 
    ##   Level 4:   20 nodes to be scored   (668 eliminated genes)

    ## 
    ##   Level 3:   12 nodes to be scored   (1272 eliminated genes)

    ## 
    ##   Level 2:   5 nodes to be scored    (1724 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2146 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0000349 generation of catalytic spliceosome for first trans...         1
    ## 2  GO:0001825                                   blastocyst formation         1
    ## 3  GO:0000050                                             urea cycle         2
    ## 4  GO:0002790                                      peptide secretion         2
    ## 5  GO:0001734         mRNA (N6-adenosine)-methyltransferase activity         1
    ## 6  GO:0004777 succinate-semialdehyde dehydrogenase (NAD+) activit...         1
    ## 7  GO:0004067                                  asparaginase activity         1
    ## 8  GO:0004035                          alkaline phosphatase activity         1
    ## 9  GO:0004408                     holocytochrome-c synthase activity         1
    ## 10 GO:0008137               NADH dehydrogenase (ubiquinone) activity         1
    ## 11 GO:0004748 ribonucleoside-diphosphate reductase activity, thio...         1
    ## 12 GO:0004325                                ferrochelatase activity         1
    ## 13 GO:0004176                       ATP-dependent peptidase activity         2
    ## 14 GO:0003730                                    mRNA 3'-UTR binding         2
    ## 15 GO:0003824                                     catalytic activity      1042
    ##    Significant Expected Fisher               type
    ## 1            1     0.00 0.0043 Biological.Process
    ## 2            1     0.00 0.0043 Biological.Process
    ## 3            1     0.01 0.0085 Biological.Process
    ## 4            1     0.01 0.0085 Biological.Process
    ## 5            1     0.01 0.0100 Molecular.Function
    ## 6            1     0.01 0.0100 Molecular.Function
    ## 7            1     0.01 0.0100 Molecular.Function
    ## 8            1     0.01 0.0100 Molecular.Function
    ## 9            1     0.01 0.0100 Molecular.Function
    ## 10           1     0.01 0.0100 Molecular.Function
    ## 11           1     0.01 0.0100 Molecular.Function
    ## 12           1     0.01 0.0100 Molecular.Function
    ## 13           1     0.02 0.0200 Molecular.Function
    ## 14           1     0.02 0.0200 Molecular.Function
    ## 15          17    10.51 0.0280 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 332 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  7 nodes to be scored    (155 eliminated genes)

    ## 
    ##   Level 10:  9 nodes to be scored    (206 eliminated genes)

    ## 
    ##   Level 9:   21 nodes to be scored   (219 eliminated genes)

    ## 
    ##   Level 8:   27 nodes to be scored   (260 eliminated genes)

    ## 
    ##   Level 7:   40 nodes to be scored   (345 eliminated genes)

    ## 
    ##   Level 6:   60 nodes to be scored   (478 eliminated genes)

    ## 
    ##   Level 5:   66 nodes to be scored   (753 eliminated genes)

    ## 
    ##   Level 4:   50 nodes to be scored   (905 eliminated genes)

    ## 
    ##   Level 3:   34 nodes to be scored   (1097 eliminated genes)

    ## 
    ##   Level 2:   11 nodes to be scored   (1271 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1374 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 173 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   9 nodes to be scored    (9 eliminated genes)

    ## 
    ##   Level 8:   12 nodes to be scored   (42 eliminated genes)

    ## 
    ##   Level 7:   19 nodes to be scored   (454 eliminated genes)

    ## 
    ##   Level 6:   34 nodes to be scored   (555 eliminated genes)

    ## 
    ##   Level 5:   36 nodes to be scored   (799 eliminated genes)

    ## 
    ##   Level 4:   29 nodes to be scored   (1180 eliminated genes)

    ## 
    ##   Level 3:   23 nodes to be scored   (1922 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (2171 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2543 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0002674     negative regulation of acute inflammatory response         4
    ## 2  GO:0001960 negative regulation of cytokine-mediated signaling ...         1
    ## 3  GO:0001967                                      suckling behavior         1
    ## 4  GO:0000183                         rDNA heterochromatin formation         1
    ## 5  GO:0000184 nuclear-transcribed mRNA catabolic process, nonsens...        10
    ## 6  GO:0001409 guanine nucleotide transmembrane transporter activi...        30
    ## 7  GO:0003682                                      chromatin binding        13
    ## 8  GO:0003723                                            RNA binding       128
    ## 9  GO:0004143                         diacylglycerol kinase activity         1
    ## 10 GO:0004441       inositol-1,4-bisphosphate 1-phosphatase activity         1
    ## 11 GO:0005289 high-affinity L-arginine transmembrane transporter ...         1
    ## 12 GO:0004568                                     chitinase activity         1
    ## 13 GO:0004033                    aldo-keto reductase (NADP) activity         1
    ## 14 GO:0003676                                   nucleic acid binding       497
    ##    Significant Expected Fisher               type
    ## 1            2     0.15 0.0076 Biological.Process
    ## 2            1     0.04 0.0369 Biological.Process
    ## 3            1     0.04 0.0369 Biological.Process
    ## 4            1     0.04 0.0369 Biological.Process
    ## 5            2     0.37 0.0496 Biological.Process
    ## 6            5     1.08 0.0038 Molecular.Function
    ## 7            3     0.47 0.0100 Molecular.Function
    ## 8           10     4.61 0.0154 Molecular.Function
    ## 9            1     0.04 0.0360 Molecular.Function
    ## 10           1     0.04 0.0360 Molecular.Function
    ## 11           1     0.04 0.0360 Molecular.Function
    ## 12           1     0.04 0.0360 Molecular.Function
    ## 13           1     0.04 0.0360 Molecular.Function
    ## 14          25    17.90 0.0412 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 117 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  3 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 10:  3 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 9:   4 nodes to be scored    (142 eliminated genes)

    ## 
    ##   Level 8:   4 nodes to be scored    (142 eliminated genes)

    ## 
    ##   Level 7:   9 nodes to be scored    (145 eliminated genes)

    ## 
    ##   Level 6:   20 nodes to be scored   (169 eliminated genes)

    ## 
    ##   Level 5:   25 nodes to be scored   (406 eliminated genes)

    ## 
    ##   Level 4:   21 nodes to be scored   (480 eliminated genes)

    ## 
    ##   Level 3:   17 nodes to be scored   (647 eliminated genes)

    ## 
    ##   Level 2:   7 nodes to be scored    (944 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1211 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 55 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   8 nodes to be scored    (263 eliminated genes)

    ## 
    ##   Level 6:   7 nodes to be scored    (277 eliminated genes)

    ## 
    ##   Level 5:   10 nodes to be scored   (463 eliminated genes)

    ## 
    ##   Level 4:   11 nodes to be scored   (471 eliminated genes)

    ## 
    ##   Level 3:   10 nodes to be scored   (842 eliminated genes)

    ## 
    ##   Level 2:   4 nodes to be scored    (1256 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1714 eliminated genes)

    ##        GO.ID                        Term Annotated Significant Expected Fisher
    ## 1 GO:0001755 neural crest cell migration         2           1     0.01 0.0110
    ## 2 GO:0000165                MAPK cascade        35           2     0.20 0.0150
    ## 3 GO:0004017   adenylate kinase activity         1           1     0.01 0.0061
    ## 4 GO:0005524                 ATP binding       263           5     1.61 0.0175
    ##                 type
    ## 1 Biological.Process
    ## 2 Biological.Process
    ## 3 Molecular.Function
    ## 4 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 24 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 7:   3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 5:   4 nodes to be scored    (154 eliminated genes)

    ## 
    ##   Level 4:   3 nodes to be scored    (170 eliminated genes)

    ## 
    ##   Level 3:   5 nodes to be scored    (215 eliminated genes)

    ## 
    ##   Level 2:   4 nodes to be scored    (408 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (514 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 46 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   6 nodes to be scored    (265 eliminated genes)

    ## 
    ##   Level 6:   6 nodes to be scored    (269 eliminated genes)

    ## 
    ##   Level 5:   8 nodes to be scored    (461 eliminated genes)

    ## 
    ##   Level 4:   6 nodes to be scored    (498 eliminated genes)

    ## 
    ##   Level 3:   10 nodes to be scored   (744 eliminated genes)

    ## 
    ##   Level 2:   4 nodes to be scored    (1278 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1816 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0000165                                           MAPK cascade        35
    ## 2 GO:0005524                                            ATP binding       263
    ## 3 GO:0004698            calcium-dependent protein kinase C activity         2
    ## 4 GO:0000900 mRNA regulatory element binding translation repress...         3
    ## 5 GO:0001640 adenylate cyclase inhibiting G protein-coupled glut...         9
    ##   Significant Expected   Fisher               type
    ## 1           4     0.15 4.60e-06 Biological.Process
    ## 2           5     0.95 1.20e-03 Molecular.Function
    ## 3           1     0.01 7.20e-03 Molecular.Function
    ## 4           1     0.01 1.08e-02 Molecular.Function
    ## 5           1     0.03 3.20e-02 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 96 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  1 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 12:  2 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 11:  4 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 9:   8 nodes to be scored    (16 eliminated genes)

    ## 
    ##   Level 8:   8 nodes to be scored    (16 eliminated genes)

    ## 
    ##   Level 7:   13 nodes to be scored   (19 eliminated genes)

    ## 
    ##   Level 6:   12 nodes to be scored   (89 eliminated genes)

    ## 
    ##   Level 5:   17 nodes to be scored   (293 eliminated genes)

    ## 
    ##   Level 4:   11 nodes to be scored   (314 eliminated genes)

    ## 
    ##   Level 3:   8 nodes to be scored    (433 eliminated genes)

    ## 
    ##   Level 2:   5 nodes to be scored    (616 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (907 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 50 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   7 nodes to be scored    (9 eliminated genes)

    ## 
    ##   Level 5:   14 nodes to be scored   (166 eliminated genes)

    ## 
    ##   Level 4:   10 nodes to be scored   (244 eliminated genes)

    ## 
    ##   Level 3:   9 nodes to be scored    (1003 eliminated genes)

    ## 
    ##   Level 2:   4 nodes to be scored    (1352 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2097 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0001945                               lymph vessel development         5
    ## 2 GO:0003333                     amino acid transmembrane transport         2
    ## 3 GO:0001658       branching involved in ureteric bud morphogenesis         4
    ## 4 GO:0001889                                      liver development         8
    ## 5 GO:0004252                     serine-type endopeptidase activity        53
    ## 6 GO:0046982                    protein heterodimerization activity         1
    ## 7 GO:0004771                               sterol esterase activity         5
    ## 8 GO:0001640 adenylate cyclase inhibiting G protein-coupled glut...         9
    ## 9 GO:0005302          L-tyrosine transmembrane transporter activity         9
    ##   Significant Expected  Fisher               type
    ## 1           2     0.03 0.00028 Biological.Process
    ## 2           1     0.01 0.01131 Biological.Process
    ## 3           1     0.02 0.02251 Biological.Process
    ## 4           1     0.05 0.04458 Biological.Process
    ## 5           3     0.29 0.00250 Molecular.Function
    ## 6           1     0.01 0.00540 Molecular.Function
    ## 7           1     0.03 0.02670 Molecular.Function
    ## 8           1     0.05 0.04770 Molecular.Function
    ## 9           1     0.05 0.04770 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 42 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   4 nodes to be scored    (9 eliminated genes)

    ## 
    ##   Level 8:   2 nodes to be scored    (9 eliminated genes)

    ## 
    ##   Level 7:   2 nodes to be scored    (9 eliminated genes)

    ## 
    ##   Level 6:   5 nodes to be scored    (9 eliminated genes)

    ## 
    ##   Level 5:   5 nodes to be scored    (21 eliminated genes)

    ## 
    ##   Level 4:   9 nodes to be scored    (36 eliminated genes)

    ## 
    ##   Level 3:   8 nodes to be scored    (321 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (457 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (603 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 38 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   3 nodes to be scored    (263 eliminated genes)

    ## 
    ##   Level 6:   3 nodes to be scored    (267 eliminated genes)

    ## 
    ##   Level 5:   6 nodes to be scored    (397 eliminated genes)

    ## 
    ##   Level 4:   9 nodes to be scored    (399 eliminated genes)

    ## 
    ##   Level 3:   11 nodes to be scored   (485 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (1152 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1737 eliminated genes)

    ##        GO.ID                                       Term Annotated Significant
    ## 1 GO:0015969 guanosine tetraphosphate metabolic process         9           1
    ## 2 GO:0019700      organic phosphonate catabolic process        12           1
    ## 3 GO:0001666                        response to hypoxia        14           1
    ## 4 GO:0005178                           integrin binding         4           1
    ## 5 GO:0003964       RNA-directed DNA polymerase activity       122           2
    ##   Expected Fisher               type
    ## 1     0.02  0.019 Biological.Process
    ## 2     0.03  0.025 Biological.Process
    ## 3     0.03  0.029 Biological.Process
    ## 4     0.01  0.011 Molecular.Function
    ## 5     0.35  0.045 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 141 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  3 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 10:  3 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 9:   6 nodes to be scored    (142 eliminated genes)

    ## 
    ##   Level 8:   9 nodes to be scored    (142 eliminated genes)

    ## 
    ##   Level 7:   15 nodes to be scored   (166 eliminated genes)

    ## 
    ##   Level 6:   22 nodes to be scored   (271 eliminated genes)

    ## 
    ##   Level 5:   28 nodes to be scored   (518 eliminated genes)

    ## 
    ##   Level 4:   23 nodes to be scored   (588 eliminated genes)

    ## 
    ##   Level 3:   20 nodes to be scored   (886 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (1058 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1251 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 108 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   7 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   7 nodes to be scored    (142 eliminated genes)

    ## 
    ##   Level 7:   14 nodes to be scored   (457 eliminated genes)

    ## 
    ##   Level 6:   18 nodes to be scored   (545 eliminated genes)

    ## 
    ##   Level 5:   21 nodes to be scored   (750 eliminated genes)

    ## 
    ##   Level 4:   19 nodes to be scored   (957 eliminated genes)

    ## 
    ##   Level 3:   14 nodes to be scored   (1639 eliminated genes)

    ## 
    ##   Level 2:   5 nodes to be scored    (2052 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2306 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0002674     negative regulation of acute inflammatory response         4
    ## 2 GO:0002221         pattern recognition receptor signaling pathway        11
    ## 3 GO:0000165                                           MAPK cascade        35
    ## 4 GO:0001818             negative regulation of cytokine production        18
    ## 5 GO:0001640 adenylate cyclase inhibiting G protein-coupled glut...         9
    ## 6 GO:0004143                         diacylglycerol kinase activity         1
    ## 7 GO:0005524                                            ATP binding       263
    ## 8 GO:0032450                     maltose alpha-glucosidase activity         2
    ## 9 GO:0000026                 alpha-1,2-mannosyltransferase activity         3
    ##   Significant Expected   Fisher               type
    ## 1           3     0.05 8.20e-06 Biological.Process
    ## 2           2     0.15 8.80e-03 Biological.Process
    ## 3           3     0.47 1.03e-02 Biological.Process
    ## 4           2     0.24 2.31e-02 Biological.Process
    ## 5           2     0.13 6.50e-03 Molecular.Function
    ## 6           1     0.01 1.40e-02 Molecular.Function
    ## 7           8     3.69 2.67e-02 Molecular.Function
    ## 8           1     0.03 2.79e-02 Molecular.Function
    ## 9           1     0.04 4.16e-02 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 245 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  7 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 10:  10 nodes to be scored   (203 eliminated genes)

    ## 
    ##   Level 9:   17 nodes to be scored   (222 eliminated genes)

    ## 
    ##   Level 8:   22 nodes to be scored   (271 eliminated genes)

    ## 
    ##   Level 7:   26 nodes to be scored   (339 eliminated genes)

    ## 
    ##   Level 6:   38 nodes to be scored   (412 eliminated genes)

    ## 
    ##   Level 5:   49 nodes to be scored   (503 eliminated genes)

    ## 
    ##   Level 4:   36 nodes to be scored   (615 eliminated genes)

    ## 
    ##   Level 3:   27 nodes to be scored   (1009 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (1217 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1346 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 126 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  3 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 9:   8 nodes to be scored    (12 eliminated genes)

    ## 
    ##   Level 8:   9 nodes to be scored    (145 eliminated genes)

    ## 
    ##   Level 7:   12 nodes to be scored   (460 eliminated genes)

    ## 
    ##   Level 6:   21 nodes to be scored   (486 eliminated genes)

    ## 
    ##   Level 5:   22 nodes to be scored   (704 eliminated genes)

    ## 
    ##   Level 4:   20 nodes to be scored   (949 eliminated genes)

    ## 
    ##   Level 3:   19 nodes to be scored   (1684 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (1935 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2389 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0000184 nuclear-transcribed mRNA catabolic process, nonsens...        10
    ## 2  GO:0006513                             protein monoubiquitination         1
    ## 3  GO:0001771                        immunological synapse formation         1
    ## 4  GO:0000183                         rDNA heterochromatin formation         1
    ## 5  GO:0002181                                cytoplasmic translation        21
    ## 6  GO:0000027                       ribosomal large subunit assembly         2
    ## 7  GO:0001409 guanine nucleotide transmembrane transporter activi...        30
    ## 8  GO:0003951                                   NAD+ kinase activity         1
    ## 9  GO:0004017                              adenylate kinase activity         1
    ## 10 GO:0005085             guanyl-nucleotide exchange factor activity        12
    ## 11 GO:0003964                   RNA-directed DNA polymerase activity       122
    ## 12 GO:0003727                            single-stranded RNA binding         2
    ## 13 GO:0003676                                   nucleic acid binding       497
    ##    Significant Expected  Fisher               type
    ## 1            3     0.20 0.00077 Biological.Process
    ## 2            1     0.02 0.01984 Biological.Process
    ## 3            1     0.02 0.01984 Biological.Process
    ## 4            1     0.02 0.01984 Biological.Process
    ## 5            3     0.42 0.02812 Biological.Process
    ## 6            1     0.04 0.03931 Biological.Process
    ## 7            3     0.56 0.01800 Molecular.Function
    ## 8            1     0.02 0.01900 Molecular.Function
    ## 9            1     0.02 0.01900 Molecular.Function
    ## 10           2     0.22 0.02000 Molecular.Function
    ## 11           6     2.29 0.02500 Molecular.Function
    ## 12           1     0.04 0.03700 Molecular.Function
    ## 13          14     9.31 0.03800 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 380 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  5 nodes to be scored    (2 eliminated genes)

    ## 
    ##   Level 11:  12 nodes to be scored   (143 eliminated genes)

    ## 
    ##   Level 10:  20 nodes to be scored   (211 eliminated genes)

    ## 
    ##   Level 9:   32 nodes to be scored   (249 eliminated genes)

    ## 
    ##   Level 8:   34 nodes to be scored   (342 eliminated genes)

    ## 
    ##   Level 7:   48 nodes to be scored   (418 eliminated genes)

    ## 
    ##   Level 6:   67 nodes to be scored   (550 eliminated genes)

    ## 
    ##   Level 5:   66 nodes to be scored   (765 eliminated genes)

    ## 
    ##   Level 4:   48 nodes to be scored   (915 eliminated genes)

    ## 
    ##   Level 3:   34 nodes to be scored   (1123 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (1278 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1396 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 163 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  3 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 9:   10 nodes to be scored   (10 eliminated genes)

    ## 
    ##   Level 8:   12 nodes to be scored   (44 eliminated genes)

    ## 
    ##   Level 7:   19 nodes to be scored   (356 eliminated genes)

    ## 
    ##   Level 6:   26 nodes to be scored   (375 eliminated genes)

    ## 
    ##   Level 5:   32 nodes to be scored   (579 eliminated genes)

    ## 
    ##   Level 4:   28 nodes to be scored   (847 eliminated genes)

    ## 
    ##   Level 3:   20 nodes to be scored   (1747 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (1998 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2568 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0001933         negative regulation of protein phosphorylation         5
    ## 2  GO:0001782                                     B cell homeostasis         7
    ## 3  GO:0001764                                       neuron migration         8
    ## 4  GO:0001731         formation of translation preinitiation complex         1
    ## 5  GO:0001960 negative regulation of cytokine-mediated signaling ...         1
    ## 6  GO:0001556                                      oocyte maturation         1
    ## 7  GO:0002092        positive regulation of receptor internalization         1
    ## 8  GO:0002026           regulation of the force of heart contraction         1
    ## 9  GO:0000038           very long-chain fatty acid metabolic process         9
    ## 10 GO:0005524                                            ATP binding       263
    ## 11 GO:0004197                   cysteine-type endopeptidase activity        18
    ## 12 GO:0004715 non-membrane spanning protein tyrosine kinase activ...         6
    ## 13 GO:0005436                    sodium:phosphate symporter activity         1
    ## 14 GO:0015175 neutral amino acid transmembrane transporter activi...         1
    ## 15 GO:0005220 inositol 1,4,5-trisphosphate-sensitive calcium-rele...         1
    ## 16 GO:0004745                     NAD-retinol dehydrogenase activity         1
    ## 17 GO:0003874           6-pyruvoyltetrahydropterin synthase activity         1
    ## 18 GO:0004862       cAMP-dependent protein kinase inhibitor activity         1
    ## 19 GO:0003989                        acetyl-CoA carboxylase activity         1
    ## 20 GO:0004423                         iduronate-2-sulfatase activity         1
    ## 21 GO:0005096                              GTPase activator activity        12
    ##    Significant Expected  Fisher               type
    ## 1            2     0.17 0.01100 Biological.Process
    ## 2            2     0.24 0.02200 Biological.Process
    ## 3            2     0.28 0.02900 Biological.Process
    ## 4            1     0.03 0.03500 Biological.Process
    ## 5            1     0.03 0.03500 Biological.Process
    ## 6            1     0.03 0.03500 Biological.Process
    ## 7            1     0.03 0.03500 Biological.Process
    ## 8            1     0.03 0.03500 Biological.Process
    ## 9            2     0.31 0.03600 Biological.Process
    ## 10          17     7.58 0.00099 Molecular.Function
    ## 11           4     0.52 0.00791 Molecular.Function
    ## 12           2     0.17 0.01141 Molecular.Function
    ## 13           1     0.03 0.02882 Molecular.Function
    ## 14           1     0.03 0.02882 Molecular.Function
    ## 15           1     0.03 0.02882 Molecular.Function
    ## 16           1     0.03 0.02882 Molecular.Function
    ## 17           1     0.03 0.02882 Molecular.Function
    ## 18           1     0.03 0.02882 Molecular.Function
    ## 19           1     0.03 0.02882 Molecular.Function
    ## 20           1     0.03 0.02882 Molecular.Function
    ## 21           2     0.35 0.04490 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 56 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  2 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 11:  3 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 10:  3 nodes to be scored    (5 eliminated genes)

    ## 
    ##   Level 9:   6 nodes to be scored    (5 eliminated genes)

    ## 
    ##   Level 8:   4 nodes to be scored    (42 eliminated genes)

    ## 
    ##   Level 7:   5 nodes to be scored    (55 eliminated genes)

    ## 
    ##   Level 6:   4 nodes to be scored    (120 eliminated genes)

    ## 
    ##   Level 5:   7 nodes to be scored    (293 eliminated genes)

    ## 
    ##   Level 4:   9 nodes to be scored    (369 eliminated genes)

    ## 
    ##   Level 3:   6 nodes to be scored    (621 eliminated genes)

    ## 
    ##   Level 2:   4 nodes to be scored    (693 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (859 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 12 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 7:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 5:   1 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 4:   2 nodes to be scored    (275 eliminated genes)

    ## 
    ##   Level 3:   5 nodes to be scored    (293 eliminated genes)

    ## 
    ##   Level 2:   1 nodes to be scored    (305 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1010 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0000472 endonucleolytic cleavage to generate mature 5'-end ...         1
    ## 2 GO:0001889                                      liver development         8
    ## 3 GO:0010181                                            FMN binding         1
    ##   Significant Expected  Fisher               type
    ## 1           1     0.00 0.00210 Biological.Process
    ## 2           1     0.02 0.01690 Biological.Process
    ## 3           1     0.00 0.00036 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 87 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   2 nodes to be scored    (10 eliminated genes)

    ## 
    ##   Level 8:   4 nodes to be scored    (14 eliminated genes)

    ## 
    ##   Level 7:   7 nodes to be scored    (15 eliminated genes)

    ## 
    ##   Level 6:   12 nodes to be scored   (144 eliminated genes)

    ## 
    ##   Level 5:   20 nodes to be scored   (373 eliminated genes)

    ## 
    ##   Level 4:   19 nodes to be scored   (420 eliminated genes)

    ## 
    ##   Level 3:   13 nodes to be scored   (737 eliminated genes)

    ## 
    ##   Level 2:   7 nodes to be scored    (884 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1041 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 64 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   7 nodes to be scored    (399 eliminated genes)

    ## 
    ##   Level 6:   9 nodes to be scored    (425 eliminated genes)

    ## 
    ##   Level 5:   13 nodes to be scored   (686 eliminated genes)

    ## 
    ##   Level 4:   11 nodes to be scored   (733 eliminated genes)

    ## 
    ##   Level 3:   13 nodes to be scored   (1282 eliminated genes)

    ## 
    ##   Level 2:   4 nodes to be scored    (1445 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2068 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0001960 negative regulation of cytokine-mediated signaling ...         1
    ## 2 GO:0000302                    response to reactive oxygen species         6
    ## 3 GO:0000184 nuclear-transcribed mRNA catabolic process, nonsens...        10
    ## 4 GO:0004492               methylmalonyl-CoA decarboxylase activity         1
    ## 5 GO:0005536                                        glucose binding         1
    ## 6 GO:0003964                   RNA-directed DNA polymerase activity       122
    ##   Significant Expected Fisher               type
    ## 1           1     0.00 0.0035 Biological.Process
    ## 2           1     0.02 0.0211 Biological.Process
    ## 3           1     0.04 0.0350 Biological.Process
    ## 4           1     0.01 0.0061 Molecular.Function
    ## 5           1     0.01 0.0061 Molecular.Function
    ## 6           3     0.75 0.0359 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 149 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  7 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 10:  9 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 9:   8 nodes to be scored    (12 eliminated genes)

    ## 
    ##   Level 8:   10 nodes to be scored   (37 eliminated genes)

    ## 
    ##   Level 7:   16 nodes to be scored   (42 eliminated genes)

    ## 
    ##   Level 6:   22 nodes to be scored   (106 eliminated genes)

    ## 
    ##   Level 5:   31 nodes to be scored   (319 eliminated genes)

    ## 
    ##   Level 4:   21 nodes to be scored   (431 eliminated genes)

    ## 
    ##   Level 3:   13 nodes to be scored   (577 eliminated genes)

    ## 
    ##   Level 2:   7 nodes to be scored    (640 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (762 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 44 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   3 nodes to be scored    (263 eliminated genes)

    ## 
    ##   Level 6:   7 nodes to be scored    (267 eliminated genes)

    ## 
    ##   Level 5:   11 nodes to be scored   (397 eliminated genes)

    ## 
    ##   Level 4:   8 nodes to be scored    (505 eliminated genes)

    ## 
    ##   Level 3:   10 nodes to be scored   (1009 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (1147 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1906 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0002230 positive regulation of defense response to virus by...         1
    ## 2 GO:0007175 negative regulation of epidermal growth factor-acti...         1
    ## 3 GO:0002933                                    lipid hydroxylation         3
    ## 4 GO:0000184 nuclear-transcribed mRNA catabolic process, nonsens...        10
    ## 5 GO:0004753                    saccharopine dehydrogenase activity         1
    ## 6 GO:0000016                                       lactase activity         4
    ## 7 GO:0008376               acetylgalactosaminyltransferase activity        14
    ##   Significant Expected Fisher               type
    ## 1           1     0.00 0.0035 Biological.Process
    ## 2           1     0.00 0.0035 Biological.Process
    ## 3           1     0.01 0.0106 Biological.Process
    ## 4           1     0.04 0.0350 Biological.Process
    ## 5           1     0.00 0.0032 Molecular.Function
    ## 6           1     0.01 0.0129 Molecular.Function
    ## 7           1     0.05 0.0445 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 134 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  6 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 10:  8 nodes to be scored    (203 eliminated genes)

    ## 
    ##   Level 9:   12 nodes to be scored   (215 eliminated genes)

    ## 
    ##   Level 8:   10 nodes to be scored   (216 eliminated genes)

    ## 
    ##   Level 7:   10 nodes to be scored   (219 eliminated genes)

    ## 
    ##   Level 6:   19 nodes to be scored   (374 eliminated genes)

    ## 
    ##   Level 5:   23 nodes to be scored   (433 eliminated genes)

    ## 
    ##   Level 4:   18 nodes to be scored   (473 eliminated genes)

    ## 
    ##   Level 3:   14 nodes to be scored   (699 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (871 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1031 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 63 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   6 nodes to be scored    (263 eliminated genes)

    ## 
    ##   Level 6:   8 nodes to be scored    (270 eliminated genes)

    ## 
    ##   Level 5:   12 nodes to be scored   (434 eliminated genes)

    ## 
    ##   Level 4:   14 nodes to be scored   (524 eliminated genes)

    ## 
    ##   Level 3:   14 nodes to be scored   (1210 eliminated genes)

    ## 
    ##   Level 2:   5 nodes to be scored    (1566 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2315 eliminated genes)

    ##        GO.ID                                               Term Annotated
    ## 1 GO:0001822                                 kidney development        65
    ## 2 GO:0000165                                       MAPK cascade        35
    ## 3 GO:0000335 negative regulation of transposition, DNA-mediated         1
    ## 4 GO:0004104                            cholinesterase activity         1
    ## 5 GO:0004382                           GDP phosphatase activity         3
    ## 6 GO:0000340                  RNA 7-methylguanosine cap binding         4
    ## 7 GO:0005178                                   integrin binding         4
    ## 8 GO:0005524                                        ATP binding       263
    ##   Significant Expected    Fisher               type
    ## 1           7     0.83 0.0000068 Biological.Process
    ## 2           4     0.45 0.0007600 Biological.Process
    ## 3           1     0.01 0.0127600 Biological.Process
    ## 4           1     0.01 0.0100000 Molecular.Function
    ## 5           1     0.03 0.0300000 Molecular.Function
    ## 6           1     0.04 0.0400000 Molecular.Function
    ## 7           1     0.04 0.0400000 Molecular.Function
    ## 8           6     2.65 0.0430000 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 94 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  4 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 10:  5 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 9:   7 nodes to be scored    (151 eliminated genes)

    ## 
    ##   Level 8:   5 nodes to be scored    (151 eliminated genes)

    ## 
    ##   Level 7:   8 nodes to be scored    (152 eliminated genes)

    ## 
    ##   Level 6:   15 nodes to be scored   (184 eliminated genes)

    ## 
    ##   Level 5:   18 nodes to be scored   (288 eliminated genes)

    ## 
    ##   Level 4:   16 nodes to be scored   (320 eliminated genes)

    ## 
    ##   Level 3:   8 nodes to be scored    (395 eliminated genes)

    ## 
    ##   Level 2:   4 nodes to be scored    (582 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (660 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 39 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   3 nodes to be scored    (263 eliminated genes)

    ## 
    ##   Level 6:   3 nodes to be scored    (267 eliminated genes)

    ## 
    ##   Level 5:   7 nodes to be scored    (397 eliminated genes)

    ## 
    ##   Level 4:   9 nodes to be scored    (399 eliminated genes)

    ## 
    ##   Level 3:   11 nodes to be scored   (921 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (1351 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1909 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0015969             guanosine tetraphosphate metabolic process         9
    ## 2 GO:0000122 negative regulation of transcription by RNA polymer...       140
    ## 3 GO:0003964                   RNA-directed DNA polymerase activity       122
    ## 4 GO:0005178                                       integrin binding         4
    ##   Significant Expected  Fisher               type
    ## 1           2     0.05 0.00099 Biological.Process
    ## 2           3     0.79 0.03687 Biological.Process
    ## 3           3     0.53 0.01400 Molecular.Function
    ## 4           1     0.02 0.01700 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 114 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  3 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 10:  3 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 9:   5 nodes to be scored    (142 eliminated genes)

    ## 
    ##   Level 8:   6 nodes to be scored    (142 eliminated genes)

    ## 
    ##   Level 7:   9 nodes to be scored    (169 eliminated genes)

    ## 
    ##   Level 6:   20 nodes to be scored   (293 eliminated genes)

    ## 
    ##   Level 5:   25 nodes to be scored   (380 eliminated genes)

    ## 
    ##   Level 4:   19 nodes to be scored   (517 eliminated genes)

    ## 
    ##   Level 3:   12 nodes to be scored   (720 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (897 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1078 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 80 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   4 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 8:   5 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 7:   9 nodes to be scored    (375 eliminated genes)

    ## 
    ##   Level 6:   14 nodes to be scored   (402 eliminated genes)

    ## 
    ##   Level 5:   15 nodes to be scored   (564 eliminated genes)

    ## 
    ##   Level 4:   12 nodes to be scored   (758 eliminated genes)

    ## 
    ##   Level 3:   13 nodes to be scored   (1127 eliminated genes)

    ## 
    ##   Level 2:   5 nodes to be scored    (1506 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2015 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0006511          ubiquitin-dependent protein catabolic process         7
    ## 2 GO:0001002 RNA polymerase III type 1 promoter sequence-specifi...         1
    ## 3 GO:0000703 oxidized pyrimidine nucleobase lesion DNA N-glycosy...         1
    ## 4 GO:0003980 UDP-glucose:glycoprotein glucosyltransferase activi...         1
    ## 5 GO:0000340                      RNA 7-methylguanosine cap binding         4
    ## 6 GO:0016836                                   hydro-lyase activity         4
    ## 7 GO:0005524                                            ATP binding       263
    ##   Significant Expected Fisher               type
    ## 1           1     0.03 0.0290 Biological.Process
    ## 2           1     0.01 0.0076 Molecular.Function
    ## 3           1     0.01 0.0076 Molecular.Function
    ## 4           1     0.01 0.0076 Molecular.Function
    ## 5           1     0.03 0.0299 Molecular.Function
    ## 6           1     0.03 0.0299 Molecular.Function
    ## 7           5     1.99 0.0422 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 5 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 5:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 4:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 3:   1 nodes to be scored    (20 eliminated genes)

    ## 
    ##   Level 2:   1 nodes to be scored    (20 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (63 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 38 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   3 nodes to be scored    (375 eliminated genes)

    ## 
    ##   Level 6:   5 nodes to be scored    (397 eliminated genes)

    ## 
    ##   Level 5:   7 nodes to be scored    (408 eliminated genes)

    ## 
    ##   Level 4:   7 nodes to be scored    (414 eliminated genes)

    ## 
    ##   Level 3:   8 nodes to be scored    (595 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (903 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1745 eliminated genes)

    ##        GO.ID                                          Term Annotated
    ## 1 GO:0003341                               cilium movement        20
    ## 2 GO:0000048                  peptidyltransferase activity         1
    ## 3 GO:0005337 nucleoside transmembrane transporter activity         1
    ##   Significant Expected Fisher               type
    ## 1           1     0.01 0.0140 Biological.Process
    ## 2           1     0.00 0.0018 Molecular.Function
    ## 3           1     0.00 0.0018 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 5 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 4:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 3:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (10 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (443 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 54 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   3 nodes to be scored    (7 eliminated genes)

    ## 
    ##   Level 7:   6 nodes to be scored    (332 eliminated genes)

    ## 
    ##   Level 6:   8 nodes to be scored    (341 eliminated genes)

    ## 
    ##   Level 5:   8 nodes to be scored    (499 eliminated genes)

    ## 
    ##   Level 4:   10 nodes to be scored   (578 eliminated genes)

    ## 
    ##   Level 3:   10 nodes to be scored   (945 eliminated genes)

    ## 
    ##   Level 2:   4 nodes to be scored    (1111 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1890 eliminated genes)

    ##        GO.ID                                   Term Annotated Significant
    ## 1 GO:0008218                        bioluminescence        10           6
    ## 2 GO:0004089         carbonate dehydratase activity         2           1
    ## 3 GO:0005245 voltage-gated calcium channel activity         7           1
    ##   Expected  Fisher               type
    ## 1     0.04 1.9e-14 Biological.Process
    ## 2     0.01 5.0e-03 Molecular.Function
    ## 3     0.02 1.8e-02 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 265 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  5 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  8 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 10:  11 nodes to be scored   (217 eliminated genes)

    ## 
    ##   Level 9:   19 nodes to be scored   (234 eliminated genes)

    ## 
    ##   Level 8:   23 nodes to be scored   (273 eliminated genes)

    ## 
    ##   Level 7:   28 nodes to be scored   (288 eliminated genes)

    ## 
    ##   Level 6:   41 nodes to be scored   (480 eliminated genes)

    ## 
    ##   Level 5:   50 nodes to be scored   (751 eliminated genes)

    ## 
    ##   Level 4:   38 nodes to be scored   (863 eliminated genes)

    ## 
    ##   Level 3:   30 nodes to be scored   (1067 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (1206 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1381 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 137 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   8 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   8 nodes to be scored    (142 eliminated genes)

    ## 
    ##   Level 7:   15 nodes to be scored   (465 eliminated genes)

    ## 
    ##   Level 6:   23 nodes to be scored   (501 eliminated genes)

    ## 
    ##   Level 5:   27 nodes to be scored   (769 eliminated genes)

    ## 
    ##   Level 4:   25 nodes to be scored   (968 eliminated genes)

    ## 
    ##   Level 3:   20 nodes to be scored   (1603 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (1850 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2395 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0001822                                     kidney development        65
    ## 2 GO:0000165                                           MAPK cascade        35
    ## 3 GO:0000335     negative regulation of transposition, DNA-mediated         1
    ## 4 GO:0001409 guanine nucleotide transmembrane transporter activi...        30
    ## 5 GO:0005524                                            ATP binding       263
    ## 6 GO:0005381            iron ion transmembrane transporter activity         1
    ## 7 GO:0004140                          dephospho-CoA kinase activity         1
    ## 8 GO:0005324             long-chain fatty acid transporter activity         1
    ## 9 GO:0003777                             microtubule motor activity        16
    ##   Significant Expected   Fisher               type
    ## 1          12     2.30 0.000001 Biological.Process
    ## 2           4     1.24 0.032000 Biological.Process
    ## 3           1     0.04 0.035000 Biological.Process
    ## 4           4     0.65 0.003600 Molecular.Function
    ## 5          12     5.68 0.009000 Molecular.Function
    ## 6           1     0.02 0.021600 Molecular.Function
    ## 7           1     0.02 0.021600 Molecular.Function
    ## 8           1     0.02 0.021600 Molecular.Function
    ## 9           2     0.35 0.045400 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 23 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   2 nodes to be scored    (30 eliminated genes)

    ## 
    ##   Level 7:   2 nodes to be scored    (38 eliminated genes)

    ## 
    ##   Level 6:   2 nodes to be scored    (40 eliminated genes)

    ## 
    ##   Level 5:   3 nodes to be scored    (44 eliminated genes)

    ## 
    ##   Level 4:   3 nodes to be scored    (111 eliminated genes)

    ## 
    ##   Level 3:   5 nodes to be scored    (154 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (468 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (532 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 35 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   5 nodes to be scored    (263 eliminated genes)

    ## 
    ##   Level 6:   5 nodes to be scored    (268 eliminated genes)

    ## 
    ##   Level 5:   6 nodes to be scored    (301 eliminated genes)

    ## 
    ##   Level 4:   6 nodes to be scored    (337 eliminated genes)

    ## 
    ##   Level 3:   7 nodes to be scored    (548 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (1068 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1777 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0000052                           citrulline metabolic process         1
    ## 2 GO:0000209                             protein polyubiquitination        30
    ## 3 GO:0003980 UDP-glucose:glycoprotein glucosyltransferase activi...         1
    ## 4 GO:0004559                             alpha-mannosidase activity         7
    ## 5 GO:0005507                                     copper ion binding        18
    ##   Significant Expected Fisher               type
    ## 1           1     0.00 0.0014 Biological.Process
    ## 2           1     0.04 0.0421 Biological.Process
    ## 3           1     0.00 0.0025 Molecular.Function
    ## 4           1     0.02 0.0175 Molecular.Function
    ## 5           1     0.05 0.0446 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 30 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   3 nodes to be scored    (2 eliminated genes)

    ## 
    ##   Level 7:   3 nodes to be scored    (2 eliminated genes)

    ## 
    ##   Level 6:   2 nodes to be scored    (60 eliminated genes)

    ## 
    ##   Level 5:   5 nodes to be scored    (90 eliminated genes)

    ## 
    ##   Level 4:   5 nodes to be scored    (102 eliminated genes)

    ## 
    ##   Level 3:   4 nodes to be scored    (557 eliminated genes)

    ## 
    ##   Level 2:   4 nodes to be scored    (938 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1000 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0005391 P-type sodium:potassium-exchanging transporter acti...         2
    ##   Significant Expected Fisher               type
    ## 1           1        0 0.0036 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 372 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  7 nodes to be scored    (11 eliminated genes)

    ## 
    ##   Level 11:  12 nodes to be scored   (156 eliminated genes)

    ## 
    ##   Level 10:  18 nodes to be scored   (215 eliminated genes)

    ## 
    ##   Level 9:   31 nodes to be scored   (235 eliminated genes)

    ## 
    ##   Level 8:   33 nodes to be scored   (334 eliminated genes)

    ## 
    ##   Level 7:   48 nodes to be scored   (422 eliminated genes)

    ## 
    ##   Level 6:   61 nodes to be scored   (527 eliminated genes)

    ## 
    ##   Level 5:   68 nodes to be scored   (836 eliminated genes)

    ## 
    ##   Level 4:   44 nodes to be scored   (929 eliminated genes)

    ## 
    ##   Level 3:   34 nodes to be scored   (1153 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (1279 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1372 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 261 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  5 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   17 nodes to be scored   (3 eliminated genes)

    ## 
    ##   Level 8:   17 nodes to be scored   (153 eliminated genes)

    ## 
    ##   Level 7:   34 nodes to be scored   (489 eliminated genes)

    ## 
    ##   Level 6:   54 nodes to be scored   (583 eliminated genes)

    ## 
    ##   Level 5:   51 nodes to be scored   (909 eliminated genes)

    ## 
    ##   Level 4:   43 nodes to be scored   (1332 eliminated genes)

    ## 
    ##   Level 3:   28 nodes to be scored   (2045 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (2278 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2570 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0002064                            epithelial cell development         4
    ## 2  GO:0000122 negative regulation of transcription by RNA polymer...       140
    ## 3  GO:0001732 formation of cytoplasmic translation initiation com...         6
    ## 4  GO:0005388                    P-type calcium transporter activity         4
    ## 5  GO:0005198                           structural molecule activity        54
    ## 6  GO:0003756                   protein disulfide isomerase activity         5
    ## 7  GO:0001786                             phosphatidylserine binding         4
    ## 8  GO:0005539                              glycosaminoglycan binding         5
    ## 9  GO:0002020                                       protease binding         5
    ## 10 GO:0003682                                      chromatin binding        13
    ## 11 GO:0004843                  cysteine-type deubiquitinase activity         6
    ##    Significant Expected  Fisher               type
    ## 1            3     0.24 0.00078 Biological.Process
    ## 2           14     8.33 0.03213 Biological.Process
    ## 3            2     0.36 0.04492 Biological.Process
    ## 4            4     0.23 0.00001 Molecular.Function
    ## 5           10     3.09 0.00025 Molecular.Function
    ## 6            3     0.29 0.00169 Molecular.Function
    ## 7            2     0.23 0.01812 Molecular.Function
    ## 8            2     0.29 0.02907 Molecular.Function
    ## 9            2     0.29 0.02907 Molecular.Function
    ## 10           3     0.74 0.03449 Molecular.Function
    ## 11           2     0.34 0.04199 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 52 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   4 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 7:   6 nodes to be scored    (5 eliminated genes)

    ## 
    ##   Level 6:   10 nodes to be scored   (73 eliminated genes)

    ## 
    ##   Level 5:   10 nodes to be scored   (99 eliminated genes)

    ## 
    ##   Level 4:   7 nodes to be scored    (173 eliminated genes)

    ## 
    ##   Level 3:   6 nodes to be scored    (273 eliminated genes)

    ## 
    ##   Level 2:   5 nodes to be scored    (337 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (848 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 51 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   6 nodes to be scored    (269 eliminated genes)

    ## 
    ##   Level 6:   7 nodes to be scored    (273 eliminated genes)

    ## 
    ##   Level 5:   8 nodes to be scored    (430 eliminated genes)

    ## 
    ##   Level 4:   11 nodes to be scored   (517 eliminated genes)

    ## 
    ##   Level 3:   11 nodes to be scored   (729 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (1292 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2073 eliminated genes)

    ##        GO.ID                                  Term Annotated Significant
    ## 1 GO:0001947                         heart looping         4           1
    ## 2 GO:0000902                    cell morphogenesis         7           1
    ## 3 GO:0019700 organic phosphonate catabolic process        12           1
    ## 4 GO:0001666                   response to hypoxia        14           1
    ## 5 GO:0000049                          tRNA binding        30           2
    ## 6 GO:0004565           beta-galactosidase activity         3           1
    ## 7 GO:0005525                           GTP binding         6           1
    ##   Expected Fisher               type
    ## 1     0.01  0.014 Biological.Process
    ## 2     0.02  0.025 Biological.Process
    ## 3     0.04  0.042 Biological.Process
    ## 4     0.05  0.049 Biological.Process
    ## 5     0.16  0.011 Molecular.Function
    ## 6     0.02  0.016 Molecular.Function
    ## 7     0.03  0.032 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 45 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   4 nodes to be scored    (9 eliminated genes)

    ## 
    ##   Level 8:   2 nodes to be scored    (9 eliminated genes)

    ## 
    ##   Level 7:   2 nodes to be scored    (9 eliminated genes)

    ## 
    ##   Level 6:   4 nodes to be scored    (9 eliminated genes)

    ## 
    ##   Level 5:   7 nodes to be scored    (113 eliminated genes)

    ## 
    ##   Level 4:   10 nodes to be scored   (121 eliminated genes)

    ## 
    ##   Level 3:   8 nodes to be scored    (486 eliminated genes)

    ## 
    ##   Level 2:   4 nodes to be scored    (661 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (779 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 30 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   3 nodes to be scored    (263 eliminated genes)

    ## 
    ##   Level 6:   4 nodes to be scored    (267 eliminated genes)

    ## 
    ##   Level 5:   4 nodes to be scored    (305 eliminated genes)

    ## 
    ##   Level 4:   5 nodes to be scored    (337 eliminated genes)

    ## 
    ##   Level 3:   8 nodes to be scored    (360 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (485 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1723 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0015969             guanosine tetraphosphate metabolic process         9
    ## 2 GO:0004180                              carboxypeptidase activity        28
    ## 3 GO:0001409 guanine nucleotide transmembrane transporter activi...        30
    ##   Significant Expected  Fisher               type
    ## 1           2     0.03 0.00022 Biological.Process
    ## 2           1     0.03 0.03000 Molecular.Function
    ## 3           1     0.03 0.03200 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 55 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   4 nodes to be scored    (9 eliminated genes)

    ## 
    ##   Level 8:   3 nodes to be scored    (9 eliminated genes)

    ## 
    ##   Level 7:   4 nodes to be scored    (9 eliminated genes)

    ## 
    ##   Level 6:   6 nodes to be scored    (74 eliminated genes)

    ## 
    ##   Level 5:   6 nodes to be scored    (121 eliminated genes)

    ## 
    ##   Level 4:   11 nodes to be scored   (131 eliminated genes)

    ## 
    ##   Level 3:   10 nodes to be scored   (569 eliminated genes)

    ## 
    ##   Level 2:   7 nodes to be scored    (743 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (907 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 30 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   3 nodes to be scored    (263 eliminated genes)

    ## 
    ##   Level 6:   3 nodes to be scored    (267 eliminated genes)

    ## 
    ##   Level 5:   4 nodes to be scored    (397 eliminated genes)

    ## 
    ##   Level 4:   6 nodes to be scored    (399 eliminated genes)

    ## 
    ##   Level 3:   8 nodes to be scored    (462 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (1094 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1588 eliminated genes)

    ##        GO.ID                                       Term Annotated Significant
    ## 1 GO:0001822                         kidney development        65           3
    ## 2 GO:0000165                               MAPK cascade        35           2
    ## 3 GO:0015969 guanosine tetraphosphate metabolic process         9           1
    ##   Expected Fisher               type
    ## 1     0.32 0.0029 Biological.Process
    ## 2     0.17 0.0116 Biological.Process
    ## 3     0.04 0.0439 Biological.Process
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 29 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  1 nodes to be scored    (7 eliminated genes)

    ## 
    ##   Level 9:   1 nodes to be scored    (7 eliminated genes)

    ## 
    ##   Level 8:   3 nodes to be scored    (8 eliminated genes)

    ## 
    ##   Level 7:   4 nodes to be scored    (8 eliminated genes)

    ## 
    ##   Level 6:   5 nodes to be scored    (8 eliminated genes)

    ## 
    ##   Level 5:   6 nodes to be scored    (36 eliminated genes)

    ## 
    ##   Level 4:   3 nodes to be scored    (99 eliminated genes)

    ## 
    ##   Level 3:   1 nodes to be scored    (228 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (239 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (246 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 34 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   6 nodes to be scored    (9 eliminated genes)

    ## 
    ##   Level 5:   9 nodes to be scored    (35 eliminated genes)

    ## 
    ##   Level 4:   6 nodes to be scored    (136 eliminated genes)

    ## 
    ##   Level 3:   6 nodes to be scored    (838 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (1062 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1705 eliminated genes)

    ##        GO.ID                                          Term Annotated
    ## 1 GO:0001843                           neural tube closure         7
    ## 2 GO:0003684                           damaged DNA binding         4
    ## 3 GO:0003943    N-acetylgalactosamine-4-sulfatase activity         8
    ## 4 GO:0005302 L-tyrosine transmembrane transporter activity         9
    ## 5 GO:0004197          cysteine-type endopeptidase activity        18
    ##   Significant Expected Fisher               type
    ## 1           1     0.00  0.005 Biological.Process
    ## 2           1     0.01  0.010 Molecular.Function
    ## 3           1     0.02  0.020 Molecular.Function
    ## 4           1     0.02  0.022 Molecular.Function
    ## 5           1     0.05  0.045 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 59 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   4 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 7:   8 nodes to be scored    (9 eliminated genes)

    ## 
    ##   Level 6:   8 nodes to be scored    (25 eliminated genes)

    ## 
    ##   Level 5:   9 nodes to be scored    (61 eliminated genes)

    ## 
    ##   Level 4:   11 nodes to be scored   (82 eliminated genes)

    ## 
    ##   Level 3:   10 nodes to be scored   (171 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (318 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (729 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 42 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 7:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   7 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 5:   9 nodes to be scored    (30 eliminated genes)

    ## 
    ##   Level 4:   11 nodes to be scored   (185 eliminated genes)

    ## 
    ##   Level 3:   10 nodes to be scored   (626 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (1241 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1867 eliminated genes)

    ##        GO.ID                                  Term Annotated Significant
    ## 1 GO:0000054 ribosomal subunit export from nucleus         3           1
    ## 2 GO:0000028      ribosomal small subunit assembly         3           1
    ## 3 GO:0001561            fatty acid alpha-oxidation         3           1
    ## 4 GO:0002933                   lipid hydroxylation         3           1
    ## 5 GO:0002161       aminoacyl-tRNA editing activity         1           1
    ## 6 GO:0003697           single-stranded DNA binding         8           1
    ## 7 GO:0005542                    folic acid binding        10           1
    ## 8 GO:0004222         metalloendopeptidase activity        15           1
    ##   Expected Fisher               type
    ## 1     0.01 0.0085 Biological.Process
    ## 2     0.01 0.0085 Biological.Process
    ## 3     0.01 0.0085 Biological.Process
    ## 4     0.01 0.0085 Biological.Process
    ## 5     0.00 0.0029 Molecular.Function
    ## 6     0.02 0.0229 Molecular.Function
    ## 7     0.03 0.0285 Molecular.Function
    ## 8     0.04 0.0425 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 144 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  3 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 10:  3 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 9:   6 nodes to be scored    (142 eliminated genes)

    ## 
    ##   Level 8:   7 nodes to be scored    (142 eliminated genes)

    ## 
    ##   Level 7:   17 nodes to be scored   (181 eliminated genes)

    ## 
    ##   Level 6:   26 nodes to be scored   (213 eliminated genes)

    ## 
    ##   Level 5:   28 nodes to be scored   (524 eliminated genes)

    ## 
    ##   Level 4:   23 nodes to be scored   (599 eliminated genes)

    ## 
    ##   Level 3:   19 nodes to be scored   (797 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (872 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1039 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 94 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   5 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   4 nodes to be scored    (100 eliminated genes)

    ## 
    ##   Level 7:   9 nodes to be scored    (409 eliminated genes)

    ## 
    ##   Level 6:   18 nodes to be scored   (427 eliminated genes)

    ## 
    ##   Level 5:   20 nodes to be scored   (663 eliminated genes)

    ## 
    ##   Level 4:   16 nodes to be scored   (907 eliminated genes)

    ## 
    ##   Level 3:   13 nodes to be scored   (1667 eliminated genes)

    ## 
    ##   Level 2:   7 nodes to be scored    (1873 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2338 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0000050                                             urea cycle         2
    ## 2 GO:0003676                                   nucleic acid binding       497
    ## 3 GO:0005085             guanyl-nucleotide exchange factor activity        12
    ## 4 GO:0004143                         diacylglycerol kinase activity         1
    ## 5 GO:0004568                                     chitinase activity         1
    ## 6 GO:0004571 mannosyl-oligosaccharide 1,2-alpha-mannosidase acti...         2
    ##   Significant Expected Fisher               type
    ## 1           1     0.02 0.0170 Biological.Process
    ## 2          12     6.98 0.0076 Molecular.Function
    ## 3           2     0.17 0.0116 Molecular.Function
    ## 4           1     0.01 0.0140 Molecular.Function
    ## 5           1     0.01 0.0140 Molecular.Function
    ## 6           1     0.03 0.0279 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 87 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  3 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 9:   5 nodes to be scored    (142 eliminated genes)

    ## 
    ##   Level 8:   4 nodes to be scored    (172 eliminated genes)

    ## 
    ##   Level 7:   6 nodes to be scored    (199 eliminated genes)

    ## 
    ##   Level 6:   15 nodes to be scored   (214 eliminated genes)

    ## 
    ##   Level 5:   20 nodes to be scored   (304 eliminated genes)

    ## 
    ##   Level 4:   13 nodes to be scored   (414 eliminated genes)

    ## 
    ##   Level 3:   9 nodes to be scored    (457 eliminated genes)

    ## 
    ##   Level 2:   4 nodes to be scored    (554 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (691 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ##        GO.ID                Term Annotated Significant Expected Fisher
    ## 1 GO:0001666 response to hypoxia        14           1     0.04  0.039
    ##                 type
    ## 1 Biological.Process
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 13 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 7:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 5:   1 nodes to be scored    (12 eliminated genes)

    ## 
    ##   Level 4:   3 nodes to be scored    (13 eliminated genes)

    ## 
    ##   Level 3:   3 nodes to be scored    (23 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (94 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (532 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 69 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   4 nodes to be scored    (7 eliminated genes)

    ## 
    ##   Level 7:   8 nodes to be scored    (332 eliminated genes)

    ## 
    ##   Level 6:   12 nodes to be scored   (346 eliminated genes)

    ## 
    ##   Level 5:   12 nodes to be scored   (534 eliminated genes)

    ## 
    ##   Level 4:   13 nodes to be scored   (707 eliminated genes)

    ## 
    ##   Level 3:   11 nodes to be scored   (1115 eliminated genes)

    ## 
    ##   Level 2:   4 nodes to be scored    (1175 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2055 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0008218                                        bioluminescence        10
    ## 2 GO:0019700                  organic phosphonate catabolic process        12
    ## 3 GO:0001409 guanine nucleotide transmembrane transporter activi...        30
    ## 4 GO:0004089                         carbonate dehydratase activity         2
    ## 5 GO:0004758                 serine C-palmitoyltransferase activity         5
    ## 6 GO:0005245                 voltage-gated calcium channel activity         7
    ##   Significant Expected Fisher               type
    ## 1           1     0.01 0.0140 Biological.Process
    ## 2           1     0.02 0.0170 Biological.Process
    ## 3           3     0.14 0.0003 Molecular.Function
    ## 4           1     0.01 0.0093 Molecular.Function
    ## 5           1     0.02 0.0232 Molecular.Function
    ## 6           1     0.03 0.0324 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 380 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  6 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  8 nodes to be scored    (11 eliminated genes)

    ## 
    ##   Level 11:  12 nodes to be scored   (158 eliminated genes)

    ## 
    ##   Level 10:  17 nodes to be scored   (211 eliminated genes)

    ## 
    ##   Level 9:   28 nodes to be scored   (233 eliminated genes)

    ## 
    ##   Level 8:   31 nodes to be scored   (321 eliminated genes)

    ## 
    ##   Level 7:   44 nodes to be scored   (385 eliminated genes)

    ## 
    ##   Level 6:   62 nodes to be scored   (444 eliminated genes)

    ## 
    ##   Level 5:   75 nodes to be scored   (624 eliminated genes)

    ## 
    ##   Level 4:   48 nodes to be scored   (857 eliminated genes)

    ## 
    ##   Level 3:   36 nodes to be scored   (1190 eliminated genes)

    ## 
    ##   Level 2:   11 nodes to be scored   (1312 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1404 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 223 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   8 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   15 nodes to be scored   (42 eliminated genes)

    ## 
    ##   Level 7:   33 nodes to be scored   (356 eliminated genes)

    ## 
    ##   Level 6:   48 nodes to be scored   (575 eliminated genes)

    ## 
    ##   Level 5:   43 nodes to be scored   (868 eliminated genes)

    ## 
    ##   Level 4:   37 nodes to be scored   (1249 eliminated genes)

    ## 
    ##   Level 3:   26 nodes to be scored   (2019 eliminated genes)

    ## 
    ##   Level 2:   11 nodes to be scored   (2271 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2620 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0001508                                       action potential         4
    ## 2  GO:0001732 formation of cytoplasmic translation initiation com...         6
    ## 3  GO:0009058                                   biosynthetic process       207
    ## 4  GO:0001771                        immunological synapse formation         1
    ## 5  GO:0006777         Mo-molybdopterin cofactor biosynthetic process         1
    ## 6  GO:0000183                         rDNA heterochromatin formation         1
    ## 7  GO:0000971 tRNA exon ligation utilizing 2',3' cyclic phosphate...         1
    ## 8  GO:0001706                                     endoderm formation         1
    ## 9  GO:0001967                                      suckling behavior         1
    ## 10 GO:0000454             snoRNA guided rRNA pseudouridine synthesis         1
    ## 11 GO:0002943                          tRNA dihydrouridine synthesis         1
    ## 12 GO:0007601                                      visual perception         1
    ## 13 GO:0000904         cell morphogenesis involved in differentiation         1
    ## 14 GO:0000105                         histidine biosynthetic process         1
    ## 15 GO:0001937 negative regulation of endothelial cell proliferati...         1
    ## 16 GO:0002100               tRNA wobble adenosine to inosine editing         1
    ## 17 GO:0003735                     structural constituent of ribosome        13
    ## 18 GO:0001786                             phosphatidylserine binding         4
    ## 19 GO:0003964                   RNA-directed DNA polymerase activity       122
    ## 20 GO:0000009                 alpha-1,6-mannosyltransferase activity         1
    ## 21 GO:0005315 inorganic phosphate transmembrane transporter activ...         1
    ## 22 GO:0004108                         citrate (Si)-synthase activity         1
    ## 23 GO:0004856                                  xylulokinase activity         1
    ## 24 GO:0004638      phosphoribosylaminoimidazole carboxylase activity         1
    ## 25 GO:0008663 2',3'-cyclic-nucleotide 2'-phosphodiesterase activi...         1
    ## 26 GO:0004736                          pyruvate carboxylase activity         1
    ## 27 GO:0004568                                     chitinase activity         1
    ## 28 GO:0003944 N-acetylglucosamine-1-phosphodiester alpha-N-acetyl...         1
    ## 29 GO:0008768                       UDP-sugar diphosphatase activity         1
    ## 30 GO:0004096                                      catalase activity         1
    ## 31 GO:0004114     3',5'-cyclic-nucleotide phosphodiesterase activity         1
    ## 32 GO:0000033                 alpha-1,3-mannosyltransferase activity         1
    ## 33 GO:0005536                                        glucose binding         1
    ##    Significant Expected Fisher               type
    ## 1            2     0.18 0.0110 Biological.Process
    ## 2            2     0.27 0.0260 Biological.Process
    ## 3           11     9.24 0.0440 Biological.Process
    ## 4            1     0.04 0.0450 Biological.Process
    ## 5            1     0.04 0.0450 Biological.Process
    ## 6            1     0.04 0.0450 Biological.Process
    ## 7            1     0.04 0.0450 Biological.Process
    ## 8            1     0.04 0.0450 Biological.Process
    ## 9            1     0.04 0.0450 Biological.Process
    ## 10           1     0.04 0.0450 Biological.Process
    ## 11           1     0.04 0.0450 Biological.Process
    ## 12           1     0.04 0.0450 Biological.Process
    ## 13           1     0.04 0.0450 Biological.Process
    ## 14           1     0.04 0.0450 Biological.Process
    ## 15           1     0.04 0.0450 Biological.Process
    ## 16           1     0.04 0.0450 Biological.Process
    ## 17           4     0.55 0.0016 Molecular.Function
    ## 18           2     0.17 0.0100 Molecular.Function
    ## 19          11     5.14 0.0124 Molecular.Function
    ## 20           1     0.04 0.0421 Molecular.Function
    ## 21           1     0.04 0.0421 Molecular.Function
    ## 22           1     0.04 0.0421 Molecular.Function
    ## 23           1     0.04 0.0421 Molecular.Function
    ## 24           1     0.04 0.0421 Molecular.Function
    ## 25           1     0.04 0.0421 Molecular.Function
    ## 26           1     0.04 0.0421 Molecular.Function
    ## 27           1     0.04 0.0421 Molecular.Function
    ## 28           1     0.04 0.0421 Molecular.Function
    ## 29           1     0.04 0.0421 Molecular.Function
    ## 30           1     0.04 0.0421 Molecular.Function
    ## 31           1     0.04 0.0421 Molecular.Function
    ## 32           1     0.04 0.0421 Molecular.Function
    ## 33           1     0.04 0.0421 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 5 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 5:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 4:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 3:   1 nodes to be scored    (6 eliminated genes)

    ## 
    ##   Level 2:   1 nodes to be scored    (95 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (107 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ##        GO.ID                       Term Annotated Significant Expected Fisher
    ## 1 GO:0006836 neurotransmitter transport         6           1        0 0.0043
    ##                 type
    ## 1 Biological.Process
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 131 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  4 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 10:  5 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 9:   8 nodes to be scored    (143 eliminated genes)

    ## 
    ##   Level 8:   8 nodes to be scored    (186 eliminated genes)

    ## 
    ##   Level 7:   13 nodes to be scored   (216 eliminated genes)

    ## 
    ##   Level 6:   23 nodes to be scored   (259 eliminated genes)

    ## 
    ##   Level 5:   31 nodes to be scored   (369 eliminated genes)

    ## 
    ##   Level 4:   19 nodes to be scored   (480 eliminated genes)

    ## 
    ##   Level 3:   10 nodes to be scored   (813 eliminated genes)

    ## 
    ##   Level 2:   6 nodes to be scored    (986 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1073 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 121 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   5 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   6 nodes to be scored    (107 eliminated genes)

    ## 
    ##   Level 7:   14 nodes to be scored   (452 eliminated genes)

    ## 
    ##   Level 6:   21 nodes to be scored   (475 eliminated genes)

    ## 
    ##   Level 5:   24 nodes to be scored   (738 eliminated genes)

    ## 
    ##   Level 4:   22 nodes to be scored   (892 eliminated genes)

    ## 
    ##   Level 3:   18 nodes to be scored   (1424 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (1747 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2440 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0008218                                        bioluminescence        10
    ## 2  GO:0001731         formation of translation preinitiation complex         1
    ## 3  GO:0002790                                      peptide secretion         2
    ## 4  GO:0000028                       ribosomal small subunit assembly         3
    ## 5  GO:0002933                                    lipid hydroxylation         3
    ## 6  GO:0004368 glycerol-3-phosphate dehydrogenase (quinone) activi...         1
    ## 7  GO:0004142      diacylglycerol cholinephosphotransferase activity         1
    ## 8  GO:0004325                                ferrochelatase activity         1
    ## 9  GO:0004176                       ATP-dependent peptidase activity         2
    ## 10 GO:0000774             adenyl-nucleotide exchange factor activity         2
    ## 11 GO:0004089                         carbonate dehydratase activity         2
    ## 12 GO:0003985                acetyl-CoA C-acetyltransferase activity         2
    ## 13 GO:0004252                     serine-type endopeptidase activity        53
    ## 14 GO:0005319                             lipid transporter activity         3
    ##    Significant Expected   Fisher               type
    ## 1            3     0.08 0.000041 Biological.Process
    ## 2            1     0.01 0.007800 Biological.Process
    ## 3            1     0.02 0.015500 Biological.Process
    ## 4            1     0.02 0.023200 Biological.Process
    ## 5            1     0.02 0.023200 Biological.Process
    ## 6            1     0.01 0.012000 Molecular.Function
    ## 7            1     0.01 0.012000 Molecular.Function
    ## 8            1     0.01 0.012000 Molecular.Function
    ## 9            1     0.02 0.024000 Molecular.Function
    ## 10           1     0.02 0.024000 Molecular.Function
    ## 11           1     0.02 0.024000 Molecular.Function
    ## 12           1     0.02 0.024000 Molecular.Function
    ## 13           3     0.63 0.024000 Molecular.Function
    ## 14           1     0.04 0.035000 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 18 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 7:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 5:   4 nodes to be scored    (105 eliminated genes)

    ## 
    ##   Level 4:   4 nodes to be scored    (123 eliminated genes)

    ## 
    ##   Level 3:   3 nodes to be scored    (240 eliminated genes)

    ## 
    ##   Level 2:   1 nodes to be scored    (285 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (302 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 18 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 7:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 5:   3 nodes to be scored    (122 eliminated genes)

    ## 
    ##   Level 4:   4 nodes to be scored    (188 eliminated genes)

    ## 
    ##   Level 3:   5 nodes to be scored    (443 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (935 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1404 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0000972 transcription-dependent tethering of RNA polymerase...         1
    ## 2 GO:0003341                                        cilium movement        20
    ##   Significant Expected Fisher               type
    ## 1           1     0.00 0.0021 Biological.Process
    ## 2           1     0.04 0.0420 Biological.Process
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 121 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  3 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 9:   5 nodes to be scored    (142 eliminated genes)

    ## 
    ##   Level 8:   5 nodes to be scored    (157 eliminated genes)

    ## 
    ##   Level 7:   9 nodes to be scored    (160 eliminated genes)

    ## 
    ##   Level 6:   20 nodes to be scored   (218 eliminated genes)

    ## 
    ##   Level 5:   27 nodes to be scored   (318 eliminated genes)

    ## 
    ##   Level 4:   21 nodes to be scored   (378 eliminated genes)

    ## 
    ##   Level 3:   14 nodes to be scored   (494 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (880 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1099 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 46 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   4 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 8:   3 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 7:   4 nodes to be scored    (115 eliminated genes)

    ## 
    ##   Level 6:   8 nodes to be scored    (136 eliminated genes)

    ## 
    ##   Level 5:   8 nodes to be scored    (142 eliminated genes)

    ## 
    ##   Level 4:   7 nodes to be scored    (200 eliminated genes)

    ## 
    ##   Level 3:   6 nodes to be scored    (383 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (671 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1626 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0001755                            neural crest cell migration         2
    ## 2 GO:0002218                   activation of innate immune response        15
    ## 3 GO:0001002 RNA polymerase III type 1 promoter sequence-specifi...         1
    ## 4 GO:0004660                   protein farnesyltransferase activity         2
    ## 5 GO:0004719 protein-L-isoaspartate (D-aspartate) O-methyltransf...         3
    ## 6 GO:0004867           serine-type endopeptidase inhibitor activity         4
    ##   Significant Expected Fisher               type
    ## 1           1     0.01 0.0057 Biological.Process
    ## 2           1     0.04 0.0419 Biological.Process
    ## 3           1     0.00 0.0022 Molecular.Function
    ## 4           1     0.00 0.0043 Molecular.Function
    ## 5           1     0.01 0.0065 Molecular.Function
    ## 6           1     0.01 0.0086 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 51 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  2 nodes to be scored    (2 eliminated genes)

    ## 
    ##   Level 9:   4 nodes to be scored    (2 eliminated genes)

    ## 
    ##   Level 8:   5 nodes to be scored    (2 eliminated genes)

    ## 
    ##   Level 7:   7 nodes to be scored    (14 eliminated genes)

    ## 
    ##   Level 6:   6 nodes to be scored    (36 eliminated genes)

    ## 
    ##   Level 5:   8 nodes to be scored    (64 eliminated genes)

    ## 
    ##   Level 4:   7 nodes to be scored    (88 eliminated genes)

    ## 
    ##   Level 3:   6 nodes to be scored    (166 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (362 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (614 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 56 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   4 nodes to be scored    (7 eliminated genes)

    ## 
    ##   Level 7:   5 nodes to be scored    (70 eliminated genes)

    ## 
    ##   Level 6:   7 nodes to be scored    (146 eliminated genes)

    ## 
    ##   Level 5:   8 nodes to be scored    (290 eliminated genes)

    ## 
    ##   Level 4:   12 nodes to be scored   (375 eliminated genes)

    ## 
    ##   Level 3:   12 nodes to be scored   (490 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (935 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2047 eliminated genes)

    ##        GO.ID                                        Term Annotated Significant
    ## 1 GO:0001516          prostaglandin biosynthetic process         2           1
    ## 2 GO:0000028            ribosomal small subunit assembly         3           1
    ## 3 GO:0001523                  retinoid metabolic process         9           1
    ## 4 GO:0005381 iron ion transmembrane transporter activity         1           1
    ## 5 GO:0004252          serine-type endopeptidase activity        53           2
    ## 6 GO:0005245      voltage-gated calcium channel activity         7           1
    ## 7 GO:0005542                          folic acid binding        10           1
    ##   Expected Fisher               type
    ## 1     0.00 0.0042 Biological.Process
    ## 2     0.01 0.0064 Biological.Process
    ## 3     0.02 0.0190 Biological.Process
    ## 4     0.00 0.0029 Molecular.Function
    ## 5     0.15 0.0093 Molecular.Function
    ## 6     0.02 0.0200 Molecular.Function
    ## 7     0.03 0.0285 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 27 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  1 nodes to be scored    (15 eliminated genes)

    ## 
    ##   Level 10:  1 nodes to be scored    (63 eliminated genes)

    ## 
    ##   Level 9:   2 nodes to be scored    (63 eliminated genes)

    ## 
    ##   Level 8:   2 nodes to be scored    (64 eliminated genes)

    ## 
    ##   Level 7:   1 nodes to be scored    (66 eliminated genes)

    ## 
    ##   Level 6:   1 nodes to be scored    (97 eliminated genes)

    ## 
    ##   Level 5:   3 nodes to be scored    (242 eliminated genes)

    ## 
    ##   Level 4:   6 nodes to be scored    (268 eliminated genes)

    ## 
    ##   Level 3:   5 nodes to be scored    (346 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (467 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (595 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 23 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   2 nodes to be scored    (263 eliminated genes)

    ## 
    ##   Level 6:   4 nodes to be scored    (267 eliminated genes)

    ## 
    ##   Level 5:   4 nodes to be scored    (275 eliminated genes)

    ## 
    ##   Level 4:   3 nodes to be scored    (443 eliminated genes)

    ## 
    ##   Level 3:   5 nodes to be scored    (685 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (802 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1010 eliminated genes)

    ##        GO.ID                                       Term Annotated Significant
    ## 1 GO:0000380 alternative mRNA splicing, via spliceosome        15           1
    ## 2 GO:0003341                            cilium movement        20           1
    ##   Expected Fisher               type
    ## 1     0.02  0.021 Biological.Process
    ## 2     0.03  0.028 Biological.Process
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 96 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  2 nodes to be scored    (11 eliminated genes)

    ## 
    ##   Level 11:  4 nodes to be scored    (15 eliminated genes)

    ## 
    ##   Level 10:  6 nodes to be scored    (65 eliminated genes)

    ## 
    ##   Level 9:   6 nodes to be scored    (66 eliminated genes)

    ## 
    ##   Level 8:   7 nodes to be scored    (74 eliminated genes)

    ## 
    ##   Level 7:   7 nodes to be scored    (80 eliminated genes)

    ## 
    ##   Level 6:   10 nodes to be scored   (268 eliminated genes)

    ## 
    ##   Level 5:   16 nodes to be scored   (378 eliminated genes)

    ## 
    ##   Level 4:   16 nodes to be scored   (492 eliminated genes)

    ## 
    ##   Level 3:   11 nodes to be scored   (787 eliminated genes)

    ## 
    ##   Level 2:   7 nodes to be scored    (947 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1069 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 34 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   4 nodes to be scored    (263 eliminated genes)

    ## 
    ##   Level 6:   5 nodes to be scored    (267 eliminated genes)

    ## 
    ##   Level 5:   6 nodes to be scored    (286 eliminated genes)

    ## 
    ##   Level 4:   5 nodes to be scored    (306 eliminated genes)

    ## 
    ##   Level 3:   8 nodes to be scored    (357 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (470 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1649 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0000076                   DNA replication checkpoint signaling         1
    ## 2 GO:0001516                     prostaglandin biosynthetic process         2
    ## 3 GO:0000381 regulation of alternative mRNA splicing, via splice...        11
    ## 4 GO:0004777 succinate-semialdehyde dehydrogenase (NAD+) activit...         1
    ## 5 GO:0004601                                    peroxidase activity         7
    ## 6 GO:0004181                       metallocarboxypeptidase activity        10
    ##   Significant Expected Fisher               type
    ## 1           1     0.00 0.0035 Biological.Process
    ## 2           1     0.01 0.0071 Biological.Process
    ## 3           1     0.04 0.0384 Biological.Process
    ## 4           1     0.00 0.0018 Molecular.Function
    ## 5           1     0.01 0.0126 Molecular.Function
    ## 6           1     0.02 0.0179 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 133 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  4 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 9:   6 nodes to be scored    (143 eliminated genes)

    ## 
    ##   Level 8:   5 nodes to be scored    (149 eliminated genes)

    ## 
    ##   Level 7:   13 nodes to be scored   (167 eliminated genes)

    ## 
    ##   Level 6:   23 nodes to be scored   (187 eliminated genes)

    ## 
    ##   Level 5:   29 nodes to be scored   (299 eliminated genes)

    ## 
    ##   Level 4:   22 nodes to be scored   (395 eliminated genes)

    ## 
    ##   Level 3:   16 nodes to be scored   (644 eliminated genes)

    ## 
    ##   Level 2:   7 nodes to be scored    (758 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (891 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 73 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   2 nodes to be scored    (100 eliminated genes)

    ## 
    ##   Level 7:   4 nodes to be scored    (383 eliminated genes)

    ## 
    ##   Level 6:   13 nodes to be scored   (397 eliminated genes)

    ## 
    ##   Level 5:   16 nodes to be scored   (409 eliminated genes)

    ## 
    ##   Level 4:   16 nodes to be scored   (606 eliminated genes)

    ## 
    ##   Level 3:   14 nodes to be scored   (1054 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (1464 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2009 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0002933                                    lipid hydroxylation         3
    ## 2  GO:0001731         formation of translation preinitiation complex         1
    ## 3  GO:0000028                       ribosomal small subunit assembly         3
    ## 4  GO:0009734                      auxin-activated signaling pathway         4
    ## 5  GO:0004342             glucosamine-6-phosphate deaminase activity         2
    ## 6  GO:0004066 asparagine synthase (glutamine-hydrolyzing) activit...         1
    ## 7  GO:0004142      diacylglycerol cholinephosphotransferase activity         1
    ## 8  GO:0004089                         carbonate dehydratase activity         2
    ## 9  GO:0004252                     serine-type endopeptidase activity        53
    ## 10 GO:0004601                                    peroxidase activity         7
    ##    Significant Expected   Fisher               type
    ## 1            2     0.02 0.000110 Biological.Process
    ## 2            1     0.01 0.006380 Biological.Process
    ## 3            1     0.02 0.019030 Biological.Process
    ## 4            1     0.03 0.025300 Biological.Process
    ## 5            2     0.01 0.000031 Molecular.Function
    ## 6            1     0.01 0.005800 Molecular.Function
    ## 7            1     0.01 0.005800 Molecular.Function
    ## 8            1     0.01 0.011500 Molecular.Function
    ## 9            2     0.31 0.036200 Molecular.Function
    ## 10           1     0.04 0.039700 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 31 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   1 nodes to be scored    (37 eliminated genes)

    ## 
    ##   Level 7:   2 nodes to be scored    (40 eliminated genes)

    ## 
    ##   Level 6:   4 nodes to be scored    (42 eliminated genes)

    ## 
    ##   Level 5:   6 nodes to be scored    (151 eliminated genes)

    ## 
    ##   Level 4:   6 nodes to be scored    (182 eliminated genes)

    ## 
    ##   Level 3:   4 nodes to be scored    (360 eliminated genes)

    ## 
    ##   Level 2:   5 nodes to be scored    (508 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (954 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 9 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 7:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 5:   1 nodes to be scored    (122 eliminated genes)

    ## 
    ##   Level 4:   2 nodes to be scored    (122 eliminated genes)

    ## 
    ##   Level 3:   2 nodes to be scored    (134 eliminated genes)

    ## 
    ##   Level 2:   1 nodes to be scored    (257 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (412 eliminated genes)

    ##        GO.ID                                 Term Annotated Significant
    ## 1 GO:0006836           neurotransmitter transport         6           1
    ## 2 GO:0008218                      bioluminescence        10           1
    ## 3 GO:0003964 RNA-directed DNA polymerase activity       122           1
    ##   Expected Fisher               type
    ## 1     0.02  0.017 Biological.Process
    ## 2     0.03  0.028 Biological.Process
    ## 3     0.04  0.044 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 24 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   2 nodes to be scored    (263 eliminated genes)

    ## 
    ##   Level 6:   3 nodes to be scored    (267 eliminated genes)

    ## 
    ##   Level 5:   3 nodes to be scored    (275 eliminated genes)

    ## 
    ##   Level 4:   3 nodes to be scored    (278 eliminated genes)

    ## 
    ##   Level 3:   7 nodes to be scored    (294 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (306 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1398 eliminated genes)

    ##        GO.ID                                   Term Annotated Significant
    ## 1 GO:0004013        adenosylhomocysteinase activity         1           1
    ## 2 GO:0005200 structural constituent of cytoskeleton         4           1
    ##   Expected Fisher               type
    ## 1        0 0.0011 Molecular.Function
    ## 2        0 0.0043 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 140 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  1 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 11:  2 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 10:  5 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 9:   7 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 8:   10 nodes to be scored   (14 eliminated genes)

    ## 
    ##   Level 7:   16 nodes to be scored   (34 eliminated genes)

    ## 
    ##   Level 6:   27 nodes to be scored   (280 eliminated genes)

    ## 
    ##   Level 5:   30 nodes to be scored   (346 eliminated genes)

    ## 
    ##   Level 4:   20 nodes to be scored   (417 eliminated genes)

    ## 
    ##   Level 3:   12 nodes to be scored   (723 eliminated genes)

    ## 
    ##   Level 2:   7 nodes to be scored    (875 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (989 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 49 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   3 nodes to be scored    (263 eliminated genes)

    ## 
    ##   Level 6:   8 nodes to be scored    (267 eliminated genes)

    ## 
    ##   Level 5:   11 nodes to be scored   (293 eliminated genes)

    ## 
    ##   Level 4:   11 nodes to be scored   (431 eliminated genes)

    ## 
    ##   Level 3:   10 nodes to be scored   (925 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (1283 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1620 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0000290 deadenylation-dependent decapping of nuclear-transc...         1
    ## 2 GO:0001782                                     B cell homeostasis         7
    ## 3 GO:0001889                                      liver development         8
    ## 4 GO:0001523                             retinoid metabolic process         9
    ## 5 GO:0002291 T cell activation via T cell receptor contact with ...        10
    ## 6 GO:0004609              phosphatidylserine decarboxylase activity        13
    ##   Significant Expected Fisher               type
    ## 1           1     0.00 0.0043 Biological.Process
    ## 2           1     0.03 0.0295 Biological.Process
    ## 3           1     0.03 0.0336 Biological.Process
    ## 4           1     0.04 0.0377 Biological.Process
    ## 5           1     0.04 0.0418 Biological.Process
    ## 6           2     0.08 0.0030 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 38 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   3 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 8:   4 nodes to be scored    (43 eliminated genes)

    ## 
    ##   Level 7:   2 nodes to be scored    (58 eliminated genes)

    ## 
    ##   Level 6:   3 nodes to be scored    (127 eliminated genes)

    ## 
    ##   Level 5:   6 nodes to be scored    (289 eliminated genes)

    ## 
    ##   Level 4:   7 nodes to be scored    (338 eliminated genes)

    ## 
    ##   Level 3:   5 nodes to be scored    (612 eliminated genes)

    ## 
    ##   Level 2:   4 nodes to be scored    (679 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (778 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 45 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   3 nodes to be scored    (7 eliminated genes)

    ## 
    ##   Level 7:   5 nodes to be scored    (332 eliminated genes)

    ## 
    ##   Level 6:   6 nodes to be scored    (341 eliminated genes)

    ## 
    ##   Level 5:   7 nodes to be scored    (377 eliminated genes)

    ## 
    ##   Level 4:   7 nodes to be scored    (395 eliminated genes)

    ## 
    ##   Level 3:   9 nodes to be scored    (608 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (633 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1973 eliminated genes)

    ##        GO.ID                                   Term Annotated Significant
    ## 1 GO:0006400                      tRNA modification         4           1
    ## 2 GO:0005245 voltage-gated calcium channel activity         7           1
    ## 3 GO:0004177                aminopeptidase activity         7           1
    ##   Expected Fisher               type
    ## 1     0.01 0.0057 Biological.Process
    ## 2     0.01 0.0130 Molecular.Function
    ## 3     0.01 0.0130 Molecular.Function

``` r
head(results_high0.5_cor_targets)
```

    ##        GO.ID                                               Term Annotated
    ## 1 GO:0006598                        polyamine catabolic process         3
    ## 2 GO:0001745                         compound eye morphogenesis         3
    ## 3 GO:0000014 single-stranded DNA endodeoxyribonuclease activity         3
    ## 4 GO:0005290     L-histidine transmembrane transporter activity         9
    ## 5 GO:0001658   branching involved in ureteric bud morphogenesis         4
    ## 6 GO:0001945                           lymph vessel development         5
    ##   Significant Expected Fisher               type         miRNA
    ## 1           1     0.00 0.0042 Biological.Process Cluster_10452
    ## 2           1     0.00 0.0042 Biological.Process Cluster_10452
    ## 3           1     0.01 0.0065 Molecular.Function Cluster_10452
    ## 4           1     0.02 0.0193 Molecular.Function Cluster_10452
    ## 5           1     0.01 0.0057 Biological.Process Cluster_11565
    ## 6           1     0.01 0.0071 Biological.Process Cluster_11565

Save results

``` r
write.csv(results_high0.5_cor_targets, "../output/27-Apul-mRNA-miRNA-interactions-topGO/miRNA_high0.5_cor_targets_topGO_FE.csv")
```

# 5 FE of specific miRNA’s targets (high 0.6 cor targets)

Loop through all miRNA and run functional enrichment on the miRNA’s
highly correlated targets (PCC magnitude \> 0.6)

``` r
interacting_miRNAs_high0.6 <- unique(high0.6_cor_bind_FA$miRNA) %>% na.omit
results_high0.6_cor_targets <- NULL  # initialize empty df

for(miRNA in interacting_miRNAs_high0.6) {
  
  # Run topGO enrichment function
  miRNA_results <- miRNA_topGO_FE(miRNA, high0.6_cor_bind_FA)
  
  # Only keep results if not empty
  if (nrow(miRNA_results) > 0) {
    
    # Add the miRNA source column
    miRNA_results$miRNA <- miRNA

    # Bind to the accumulating results data frame
    results_high0.6_cor_targets <- rbind(results_high0.6_cor_targets, miRNA_results)
  }
}
```

    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 15 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   2 nodes to be scored    (9 eliminated genes)

    ## 
    ##   Level 5:   3 nodes to be scored    (35 eliminated genes)

    ## 
    ##   Level 4:   3 nodes to be scored    (40 eliminated genes)

    ## 
    ##   Level 3:   1 nodes to be scored    (145 eliminated genes)

    ## 
    ##   Level 2:   1 nodes to be scored    (149 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (186 eliminated genes)

    ##        GO.ID                                           Term Annotated
    ## 1 GO:0005290 L-histidine transmembrane transporter activity         9
    ##   Significant Expected Fisher               type
    ## 1           1        0 0.0032 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 42 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  1 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 12:  2 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 11:  3 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 10:  3 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 9:   5 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 8:   4 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 7:   4 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 6:   4 nodes to be scored    (72 eliminated genes)

    ## 
    ##   Level 5:   6 nodes to be scored    (125 eliminated genes)

    ## 
    ##   Level 4:   4 nodes to be scored    (138 eliminated genes)

    ## 
    ##   Level 3:   1 nodes to be scored    (209 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (239 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (246 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 32 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   5 nodes to be scored    (12 eliminated genes)

    ## 
    ##   Level 5:   7 nodes to be scored    (153 eliminated genes)

    ## 
    ##   Level 4:   7 nodes to be scored    (227 eliminated genes)

    ## 
    ##   Level 3:   6 nodes to be scored    (465 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (716 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1029 eliminated genes)

    ##        GO.ID                                             Term Annotated
    ## 1 GO:0001658 branching involved in ureteric bud morphogenesis         4
    ## 2 GO:0001945                         lymph vessel development         5
    ## 3 GO:0004252               serine-type endopeptidase activity        53
    ## 4 GO:0004771                         sterol esterase activity         5
    ## 5 GO:0004966                        galanin receptor activity        12
    ## 6 GO:0004222                    metalloendopeptidase activity        15
    ## 7 GO:0005201      extracellular matrix structural constituent        17
    ##   Significant Expected Fisher               type
    ## 1           1     0.01 0.0057 Biological.Process
    ## 2           1     0.01 0.0071 Biological.Process
    ## 3           2     0.15 0.0093 Molecular.Function
    ## 4           1     0.01 0.0143 Molecular.Function
    ## 5           1     0.03 0.0341 Molecular.Function
    ## 6           1     0.04 0.0425 Molecular.Function
    ## 7           1     0.05 0.0480 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 47 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   1 nodes to be scored    (15 eliminated genes)

    ## 
    ##   Level 7:   4 nodes to be scored    (15 eliminated genes)

    ## 
    ##   Level 6:   6 nodes to be scored    (15 eliminated genes)

    ## 
    ##   Level 5:   9 nodes to be scored    (64 eliminated genes)

    ## 
    ##   Level 4:   7 nodes to be scored    (78 eliminated genes)

    ## 
    ##   Level 3:   9 nodes to be scored    (136 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (363 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (528 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 17 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 7:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 5:   3 nodes to be scored    (122 eliminated genes)

    ## 
    ##   Level 4:   4 nodes to be scored    (178 eliminated genes)

    ## 
    ##   Level 3:   4 nodes to be scored    (317 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (460 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (934 eliminated genes)

    ##        GO.ID                                 Term Annotated Significant
    ## 1 GO:0002218 activation of innate immune response        15           1
    ## 2 GO:0009034               tryptophanase activity         2           1
    ## 3 GO:0003964 RNA-directed DNA polymerase activity       122           2
    ##   Expected Fisher               type
    ## 1     0.03 0.0320 Biological.Process
    ## 2     0.00 0.0029 Molecular.Function
    ## 3     0.18 0.0108 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 60 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   2 nodes to be scored    (12 eliminated genes)

    ## 
    ##   Level 8:   2 nodes to be scored    (52 eliminated genes)

    ## 
    ##   Level 7:   5 nodes to be scored    (55 eliminated genes)

    ## 
    ##   Level 6:   8 nodes to be scored    (57 eliminated genes)

    ## 
    ##   Level 5:   12 nodes to be scored   (111 eliminated genes)

    ## 
    ##   Level 4:   9 nodes to be scored    (148 eliminated genes)

    ## 
    ##   Level 3:   10 nodes to be scored   (337 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (574 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (743 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 31 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   3 nodes to be scored    (263 eliminated genes)

    ## 
    ##   Level 6:   4 nodes to be scored    (267 eliminated genes)

    ## 
    ##   Level 5:   5 nodes to be scored    (397 eliminated genes)

    ## 
    ##   Level 4:   6 nodes to be scored    (455 eliminated genes)

    ## 
    ##   Level 3:   8 nodes to be scored    (610 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (765 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1479 eliminated genes)

    ##        GO.ID                                 Term Annotated Significant
    ## 1 GO:0002040               sprouting angiogenesis        12           1
    ## 2 GO:0002218 activation of innate immune response        15           1
    ## 3 GO:0003964 RNA-directed DNA polymerase activity       122           3
    ## 4 GO:0009034               tryptophanase activity         2           1
    ##   Expected Fisher               type
    ## 1     0.03 0.0340 Biological.Process
    ## 2     0.04 0.0420 Biological.Process
    ## 3     0.26 0.0015 Molecular.Function
    ## 4     0.00 0.0043 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 56 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   2 nodes to be scored    (12 eliminated genes)

    ## 
    ##   Level 8:   2 nodes to be scored    (52 eliminated genes)

    ## 
    ##   Level 7:   5 nodes to be scored    (55 eliminated genes)

    ## 
    ##   Level 6:   8 nodes to be scored    (57 eliminated genes)

    ## 
    ##   Level 5:   11 nodes to be scored   (111 eliminated genes)

    ## 
    ##   Level 4:   8 nodes to be scored    (148 eliminated genes)

    ## 
    ##   Level 3:   8 nodes to be scored    (311 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (534 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (697 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 17 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 7:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 5:   3 nodes to be scored    (122 eliminated genes)

    ## 
    ##   Level 4:   4 nodes to be scored    (178 eliminated genes)

    ## 
    ##   Level 3:   4 nodes to be scored    (317 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (460 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (934 eliminated genes)

    ##        GO.ID                                 Term Annotated Significant
    ## 1 GO:0002040               sprouting angiogenesis        12           1
    ## 2 GO:0002218 activation of innate immune response        15           1
    ## 3 GO:0003964 RNA-directed DNA polymerase activity       122           3
    ## 4 GO:0009034               tryptophanase activity         2           1
    ##   Expected  Fisher               type
    ## 1     0.03 0.02500 Biological.Process
    ## 2     0.03 0.03200 Biological.Process
    ## 3     0.22 0.00078 Molecular.Function
    ## 4     0.00 0.00360 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 9 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   1 nodes to be scored    (15 eliminated genes)

    ## 
    ##   Level 5:   1 nodes to be scored    (25 eliminated genes)

    ## 
    ##   Level 4:   1 nodes to be scored    (25 eliminated genes)

    ## 
    ##   Level 3:   1 nodes to be scored    (38 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (218 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (246 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 18 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   2 nodes to be scored    (263 eliminated genes)

    ## 
    ##   Level 6:   2 nodes to be scored    (267 eliminated genes)

    ## 
    ##   Level 5:   2 nodes to be scored    (275 eliminated genes)

    ## 
    ##   Level 4:   3 nodes to be scored    (277 eliminated genes)

    ## 
    ##   Level 3:   5 nodes to be scored    (293 eliminated genes)

    ## 
    ##   Level 2:   1 nodes to be scored    (802 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1010 eliminated genes)

    ##        GO.ID                           Term Annotated Significant Expected
    ## 1 GO:0001701 in utero embryonic development        15           1     0.01
    ##   Fisher               type
    ## 1  0.011 Biological.Process
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## [1] GO.ID       Term        Annotated   Significant Expected    Fisher     
    ## <0 rows> (or 0-length row.names)
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 10 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   1 nodes to be scored    (65 eliminated genes)

    ## 
    ##   Level 5:   1 nodes to be scored    (65 eliminated genes)

    ## 
    ##   Level 4:   2 nodes to be scored    (65 eliminated genes)

    ## 
    ##   Level 3:   1 nodes to be scored    (192 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (221 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (246 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 16 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 6:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 5:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 4:   4 nodes to be scored    (40 eliminated genes)

    ## 
    ##   Level 3:   4 nodes to be scored    (61 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (212 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (795 eliminated genes)

    ##        GO.ID                      Term Annotated Significant Expected Fisher
    ## 1 GO:0001822        kidney development        65           1     0.05  0.046
    ## 2 GO:0005096 GTPase activator activity        12           1     0.01  0.013
    ## 3 GO:0004180 carboxypeptidase activity        28           1     0.03  0.030
    ## 4 GO:0001653 peptide receptor activity        40           1     0.04  0.043
    ##                 type
    ## 1 Biological.Process
    ## 2 Molecular.Function
    ## 3 Molecular.Function
    ## 4 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 74 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  3 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 9:   4 nodes to be scored    (142 eliminated genes)

    ## 
    ##   Level 8:   3 nodes to be scored    (143 eliminated genes)

    ## 
    ##   Level 7:   5 nodes to be scored    (145 eliminated genes)

    ## 
    ##   Level 6:   11 nodes to be scored   (172 eliminated genes)

    ## 
    ##   Level 5:   16 nodes to be scored   (268 eliminated genes)

    ## 
    ##   Level 4:   12 nodes to be scored   (323 eliminated genes)

    ## 
    ##   Level 3:   7 nodes to be scored    (455 eliminated genes)

    ## 
    ##   Level 2:   5 nodes to be scored    (768 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (906 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 63 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   3 nodes to be scored    (100 eliminated genes)

    ## 
    ##   Level 7:   7 nodes to be scored    (383 eliminated genes)

    ## 
    ##   Level 6:   11 nodes to be scored   (398 eliminated genes)

    ## 
    ##   Level 5:   12 nodes to be scored   (429 eliminated genes)

    ## 
    ##   Level 4:   11 nodes to be scored   (504 eliminated genes)

    ## 
    ##   Level 3:   10 nodes to be scored   (797 eliminated genes)

    ## 
    ##   Level 2:   4 nodes to be scored    (1255 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1890 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0001825                                   blastocyst formation         1
    ## 2 GO:0004035                          alkaline phosphatase activity         1
    ## 3 GO:0004325                                ferrochelatase activity         1
    ## 4 GO:0004777 succinate-semialdehyde dehydrogenase (NAD+) activit...         1
    ## 5 GO:0008137               NADH dehydrogenase (ubiquinone) activity         1
    ## 6 GO:0004176                       ATP-dependent peptidase activity         2
    ## 7 GO:0003943             N-acetylgalactosamine-4-sulfatase activity         8
    ## 8 GO:0003697                            single-stranded DNA binding         8
    ##   Significant Expected Fisher               type
    ## 1           1     0.00 0.0014 Biological.Process
    ## 2           1     0.00 0.0043 Molecular.Function
    ## 3           1     0.00 0.0043 Molecular.Function
    ## 4           1     0.00 0.0043 Molecular.Function
    ## 5           1     0.00 0.0043 Molecular.Function
    ## 6           1     0.01 0.0086 Molecular.Function
    ## 7           1     0.03 0.0341 Molecular.Function
    ## 8           1     0.03 0.0341 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 45 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   3 nodes to be scored    (28 eliminated genes)

    ## 
    ##   Level 6:   6 nodes to be scored    (37 eliminated genes)

    ## 
    ##   Level 5:   10 nodes to be scored   (58 eliminated genes)

    ## 
    ##   Level 4:   9 nodes to be scored    (151 eliminated genes)

    ## 
    ##   Level 3:   10 nodes to be scored   (497 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (654 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (704 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 59 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   5 nodes to be scored    (375 eliminated genes)

    ## 
    ##   Level 6:   8 nodes to be scored    (397 eliminated genes)

    ## 
    ##   Level 5:   9 nodes to be scored    (560 eliminated genes)

    ## 
    ##   Level 4:   12 nodes to be scored   (580 eliminated genes)

    ## 
    ##   Level 3:   14 nodes to be scored   (790 eliminated genes)

    ## 
    ##   Level 2:   6 nodes to be scored    (1413 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2341 eliminated genes)

    ##        GO.ID                                         Term Annotated Significant
    ## 1 GO:0032259                                  methylation         5           1
    ## 2 GO:0000038 very long-chain fatty acid metabolic process         9           1
    ## 3 GO:0000278                           mitotic cell cycle       160           2
    ## 4 GO:0004568                           chitinase activity         1           1
    ## 5 GO:0005178                             integrin binding         4           1
    ##   Expected Fisher               type
    ## 1     0.02  0.018 Biological.Process
    ## 2     0.03  0.032 Biological.Process
    ## 3     0.57  0.047 Biological.Process
    ## 4     0.01  0.005 Molecular.Function
    ## 5     0.02  0.020 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 5 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 5:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 4:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 3:   1 nodes to be scored    (20 eliminated genes)

    ## 
    ##   Level 2:   1 nodes to be scored    (20 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (63 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 17 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   2 nodes to be scored    (263 eliminated genes)

    ## 
    ##   Level 6:   2 nodes to be scored    (267 eliminated genes)

    ## 
    ##   Level 5:   2 nodes to be scored    (275 eliminated genes)

    ## 
    ##   Level 4:   2 nodes to be scored    (277 eliminated genes)

    ## 
    ##   Level 3:   5 nodes to be scored    (293 eliminated genes)

    ## 
    ##   Level 2:   1 nodes to be scored    (305 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1010 eliminated genes)

    ##        GO.ID            Term Annotated Significant Expected Fisher
    ## 1 GO:0003341 cilium movement        20           1     0.01  0.014
    ##                 type
    ## 1 Biological.Process
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 12 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 7:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 5:   1 nodes to be scored    (35 eliminated genes)

    ## 
    ##   Level 4:   1 nodes to be scored    (44 eliminated genes)

    ## 
    ##   Level 3:   3 nodes to be scored    (71 eliminated genes)

    ## 
    ##   Level 2:   4 nodes to be scored    (246 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (352 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ##        GO.ID         Term Annotated Significant Expected Fisher
    ## 1 GO:0000165 MAPK cascade        35           1     0.02  0.025
    ##                 type
    ## 1 Biological.Process
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 78 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  1 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 12:  2 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 11:  4 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 9:   7 nodes to be scored    (16 eliminated genes)

    ## 
    ##   Level 8:   7 nodes to be scored    (16 eliminated genes)

    ## 
    ##   Level 7:   10 nodes to be scored   (17 eliminated genes)

    ## 
    ##   Level 6:   9 nodes to be scored    (87 eliminated genes)

    ## 
    ##   Level 5:   12 nodes to be scored   (174 eliminated genes)

    ## 
    ##   Level 4:   8 nodes to be scored    (187 eliminated genes)

    ## 
    ##   Level 3:   6 nodes to be scored    (262 eliminated genes)

    ## 
    ##   Level 2:   5 nodes to be scored    (356 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (683 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 23 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 6:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 5:   7 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 4:   5 nodes to be scored    (58 eliminated genes)

    ## 
    ##   Level 3:   5 nodes to be scored    (713 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (956 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1564 eliminated genes)

    ##        GO.ID                                             Term Annotated
    ## 1 GO:0001945                         lymph vessel development         5
    ## 2 GO:0001658 branching involved in ureteric bud morphogenesis         4
    ## 3 GO:0001889                                liver development         8
    ## 4 GO:0004771                         sterol esterase activity         5
    ##   Significant Expected  Fisher               type
    ## 1           2     0.02 0.00015 Biological.Process
    ## 2           1     0.02 0.01692 Biological.Process
    ## 3           1     0.03 0.03360 Biological.Process
    ## 4           1     0.01 0.01100 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## [1] GO.ID       Term        Annotated   Significant Expected    Fisher     
    ## <0 rows> (or 0-length row.names)
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 61 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   7 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 6:   10 nodes to be scored   (87 eliminated genes)

    ## 
    ##   Level 5:   12 nodes to be scored   (246 eliminated genes)

    ## 
    ##   Level 4:   11 nodes to be scored   (439 eliminated genes)

    ## 
    ##   Level 3:   8 nodes to be scored    (731 eliminated genes)

    ## 
    ##   Level 2:   7 nodes to be scored    (933 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1074 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 71 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   4 nodes to be scored    (100 eliminated genes)

    ## 
    ##   Level 7:   10 nodes to be scored   (407 eliminated genes)

    ## 
    ##   Level 6:   11 nodes to be scored   (491 eliminated genes)

    ## 
    ##   Level 5:   13 nodes to be scored   (740 eliminated genes)

    ## 
    ##   Level 4:   13 nodes to be scored   (778 eliminated genes)

    ## 
    ##   Level 3:   10 nodes to be scored   (1160 eliminated genes)

    ## 
    ##   Level 2:   4 nodes to be scored    (1619 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2195 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0000165                                           MAPK cascade        35
    ## 2 GO:0002674     negative regulation of acute inflammatory response         4
    ## 3 GO:0001409 guanine nucleotide transmembrane transporter activi...        30
    ## 4 GO:0004623                              phospholipase A2 activity         4
    ## 5 GO:0005524                                            ATP binding       263
    ## 6 GO:0001640 adenylate cyclase inhibiting G protein-coupled glut...         9
    ## 7 GO:0003924                                        GTPase activity        66
    ##   Significant Expected Fisher               type
    ## 1           2     0.15 0.0084 Biological.Process
    ## 2           1     0.02 0.0169 Biological.Process
    ## 3           2     0.16 0.0110 Molecular.Function
    ## 4           1     0.02 0.0210 Molecular.Function
    ## 5           4     1.42 0.0460 Molecular.Function
    ## 6           1     0.05 0.0480 Molecular.Function
    ## 7           2     0.36 0.0480 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 76 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   3 nodes to be scored    (10 eliminated genes)

    ## 
    ##   Level 8:   5 nodes to be scored    (14 eliminated genes)

    ## 
    ##   Level 7:   6 nodes to be scored    (44 eliminated genes)

    ## 
    ##   Level 6:   10 nodes to be scored   (135 eliminated genes)

    ## 
    ##   Level 5:   18 nodes to be scored   (335 eliminated genes)

    ## 
    ##   Level 4:   16 nodes to be scored   (376 eliminated genes)

    ## 
    ##   Level 3:   10 nodes to be scored   (487 eliminated genes)

    ## 
    ##   Level 2:   5 nodes to be scored    (539 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (691 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 30 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 7:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 5:   5 nodes to be scored    (152 eliminated genes)

    ## 
    ##   Level 4:   7 nodes to be scored    (166 eliminated genes)

    ## 
    ##   Level 3:   7 nodes to be scored    (511 eliminated genes)

    ## 
    ##   Level 2:   5 nodes to be scored    (1033 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1672 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0000038           very long-chain fatty acid metabolic process         9
    ## 2 GO:0000184 nuclear-transcribed mRNA catabolic process, nonsens...        10
    ## 3 GO:0001409 guanine nucleotide transmembrane transporter activi...        30
    ## 4 GO:0005085             guanyl-nucleotide exchange factor activity        12
    ##   Significant Expected Fisher               type
    ## 1           1     0.03 0.0320 Biological.Process
    ## 2           1     0.04 0.0350 Biological.Process
    ## 3           2     0.06 0.0016 Molecular.Function
    ## 4           1     0.03 0.0257 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 161 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  3 nodes to be scored    (2 eliminated genes)

    ## 
    ##   Level 11:  5 nodes to be scored    (143 eliminated genes)

    ## 
    ##   Level 10:  8 nodes to be scored    (143 eliminated genes)

    ## 
    ##   Level 9:   9 nodes to be scored    (145 eliminated genes)

    ## 
    ##   Level 8:   10 nodes to be scored   (156 eliminated genes)

    ## 
    ##   Level 7:   16 nodes to be scored   (176 eliminated genes)

    ## 
    ##   Level 6:   27 nodes to be scored   (306 eliminated genes)

    ## 
    ##   Level 5:   34 nodes to be scored   (353 eliminated genes)

    ## 
    ##   Level 4:   22 nodes to be scored   (432 eliminated genes)

    ## 
    ##   Level 3:   15 nodes to be scored   (756 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (1007 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1102 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 63 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   4 nodes to be scored    (9 eliminated genes)

    ## 
    ##   Level 8:   4 nodes to be scored    (42 eliminated genes)

    ## 
    ##   Level 7:   7 nodes to be scored    (322 eliminated genes)

    ## 
    ##   Level 6:   8 nodes to be scored    (332 eliminated genes)

    ## 
    ##   Level 5:   10 nodes to be scored   (501 eliminated genes)

    ## 
    ##   Level 4:   10 nodes to be scored   (531 eliminated genes)

    ## 
    ##   Level 3:   12 nodes to be scored   (1073 eliminated genes)

    ## 
    ##   Level 2:   5 nodes to be scored    (1355 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2089 eliminated genes)

    ##        GO.ID                                             Term Annotated
    ## 1 GO:0000289 nuclear-transcribed mRNA poly(A) tail shortening         2
    ## 2 GO:0006886                  intracellular protein transport         5
    ## 3 GO:0001782                               B cell homeostasis         7
    ## 4 GO:0001764                                 neuron migration         8
    ## 5 GO:0001222                transcription corepressor binding         3
    ## 6 GO:0005242      inward rectifier potassium channel activity         9
    ## 7 GO:0003777                       microtubule motor activity        16
    ##   Significant Expected Fisher               type
    ## 1           1     0.01 0.0099 Biological.Process
    ## 2           1     0.02 0.0246 Biological.Process
    ## 3           1     0.03 0.0343 Biological.Process
    ## 4           1     0.04 0.0391 Biological.Process
    ## 5           1     0.01 0.0086 Molecular.Function
    ## 6           1     0.03 0.0257 Molecular.Function
    ## 7           1     0.05 0.0452 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 7 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 5:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 4:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 3:   3 nodes to be scored    (293 eliminated genes)

    ## 
    ##   Level 2:   1 nodes to be scored    (293 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (814 eliminated genes)

    ## [1] GO.ID       Term        Annotated   Significant Expected    Fisher     
    ## <0 rows> (or 0-length row.names)
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 39 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   2 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 8:   2 nodes to be scored    (16 eliminated genes)

    ## 
    ##   Level 7:   4 nodes to be scored    (16 eliminated genes)

    ## 
    ##   Level 6:   6 nodes to be scored    (16 eliminated genes)

    ## 
    ##   Level 5:   7 nodes to be scored    (30 eliminated genes)

    ## 
    ##   Level 4:   5 nodes to be scored    (34 eliminated genes)

    ## 
    ##   Level 3:   5 nodes to be scored    (50 eliminated genes)

    ## 
    ##   Level 2:   4 nodes to be scored    (84 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (452 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 6 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 6:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 5:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 4:   1 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 3:   1 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 2:   1 nodes to be scored    (6 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (115 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0002230 positive regulation of defense response to virus by...         1
    ## 2 GO:0002218                   activation of innate immune response        15
    ## 3 GO:0004753                    saccharopine dehydrogenase activity         1
    ##   Significant Expected  Fisher               type
    ## 1           1     0.00 0.00140 Biological.Process
    ## 2           1     0.02 0.02120 Biological.Process
    ## 3           1     0.00 0.00036 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 10 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   1 nodes to be scored    (65 eliminated genes)

    ## 
    ##   Level 5:   1 nodes to be scored    (65 eliminated genes)

    ## 
    ##   Level 4:   2 nodes to be scored    (65 eliminated genes)

    ## 
    ##   Level 3:   1 nodes to be scored    (192 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (221 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (246 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 15 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   2 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 5:   2 nodes to be scored    (125 eliminated genes)

    ## 
    ##   Level 4:   3 nodes to be scored    (201 eliminated genes)

    ## 
    ##   Level 3:   3 nodes to be scored    (213 eliminated genes)

    ## 
    ##   Level 2:   1 nodes to be scored    (336 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (774 eliminated genes)

    ##        GO.ID                                 Term Annotated Significant
    ## 1 GO:0001822                   kidney development        65           1
    ## 2 GO:0004382             GDP phosphatase activity         3           1
    ## 3 GO:0003964 RNA-directed DNA polymerase activity       122           2
    ##   Expected Fisher               type
    ## 1     0.05 0.0460 Biological.Process
    ## 2     0.00 0.0032 Molecular.Function
    ## 3     0.13 0.0056 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 20 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   3 nodes to be scored    (18 eliminated genes)

    ## 
    ##   Level 5:   4 nodes to be scored    (38 eliminated genes)

    ## 
    ##   Level 4:   4 nodes to be scored    (196 eliminated genes)

    ## 
    ##   Level 3:   2 nodes to be scored    (301 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (452 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (626 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 6 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 5:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 4:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 3:   2 nodes to be scored    (264 eliminated genes)

    ## 
    ##   Level 2:   1 nodes to be scored    (497 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (809 eliminated genes)

    ##        GO.ID                                       Term Annotated Significant
    ## 1 GO:0001818 negative regulation of cytokine production        18           1
    ##   Expected Fisher               type
    ## 1     0.01  0.013 Biological.Process
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 26 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   3 nodes to be scored    (263 eliminated genes)

    ## 
    ##   Level 6:   3 nodes to be scored    (267 eliminated genes)

    ## 
    ##   Level 5:   3 nodes to be scored    (397 eliminated genes)

    ## 
    ##   Level 4:   5 nodes to be scored    (399 eliminated genes)

    ## 
    ##   Level 3:   7 nodes to be scored    (427 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (1059 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1422 eliminated genes)

    ##        GO.ID                                 Term Annotated Significant
    ## 1 GO:0003964 RNA-directed DNA polymerase activity       122           2
    ##   Expected Fisher               type
    ## 1     0.18  0.011 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## [1] GO.ID       Term        Annotated   Significant Expected    Fisher     
    ## <0 rows> (or 0-length row.names)
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 5 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 4:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 3:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (10 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (443 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 6 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 6:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 5:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 4:   1 nodes to be scored    (66 eliminated genes)

    ## 
    ##   Level 3:   1 nodes to be scored    (181 eliminated genes)

    ## 
    ##   Level 2:   1 nodes to be scored    (181 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (465 eliminated genes)

    ##        GO.ID                Term Annotated Significant Expected  Fisher
    ## 1 GO:0008218     bioluminescence        10           3     0.02 2.6e-07
    ## 2 GO:0005509 calcium ion binding        66           1     0.02 2.4e-02
    ##                 type
    ## 1 Biological.Process
    ## 2 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 57 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 8:   3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   6 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   8 nodes to be scored    (113 eliminated genes)

    ## 
    ##   Level 5:   8 nodes to be scored    (178 eliminated genes)

    ## 
    ##   Level 4:   11 nodes to be scored   (454 eliminated genes)

    ## 
    ##   Level 3:   11 nodes to be scored   (719 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (922 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1091 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 71 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   8 nodes to be scored    (263 eliminated genes)

    ## 
    ##   Level 6:   9 nodes to be scored    (276 eliminated genes)

    ## 
    ##   Level 5:   15 nodes to be scored   (517 eliminated genes)

    ## 
    ##   Level 4:   15 nodes to be scored   (623 eliminated genes)

    ## 
    ##   Level 3:   14 nodes to be scored   (1124 eliminated genes)

    ## 
    ##   Level 2:   6 nodes to be scored    (1469 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2106 eliminated genes)

    ##        GO.ID                                     Term Annotated Significant
    ## 1 GO:0001822                       kidney development        65           4
    ## 2 GO:0000165                             MAPK cascade        35           3
    ## 3 GO:0000723                     telomere maintenance         4           1
    ## 4 GO:0005178                         integrin binding         4           1
    ## 5 GO:0004674 protein serine/threonine kinase activity        55           2
    ##   Expected Fisher               type
    ## 1     0.51 0.0011 Biological.Process
    ## 2     0.27 0.0020 Biological.Process
    ## 3     0.03 0.0309 Biological.Process
    ## 4     0.02 0.0230 Molecular.Function
    ## 5     0.32 0.0390 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 5 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 4:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 3:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 2:   1 nodes to be scored    (497 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (809 eliminated genes)

    ## [1] GO.ID       Term        Annotated   Significant Expected    Fisher     
    ## <0 rows> (or 0-length row.names)
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## [1] GO.ID       Term        Annotated   Significant Expected    Fisher     
    ## <0 rows> (or 0-length row.names)
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 186 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  5 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 10:  8 nodes to be scored    (210 eliminated genes)

    ## 
    ##   Level 9:   12 nodes to be scored   (212 eliminated genes)

    ## 
    ##   Level 8:   14 nodes to be scored   (293 eliminated genes)

    ## 
    ##   Level 7:   17 nodes to be scored   (330 eliminated genes)

    ## 
    ##   Level 6:   29 nodes to be scored   (428 eliminated genes)

    ## 
    ##   Level 5:   36 nodes to be scored   (545 eliminated genes)

    ## 
    ##   Level 4:   29 nodes to be scored   (794 eliminated genes)

    ## 
    ##   Level 3:   21 nodes to be scored   (996 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (1150 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1348 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 122 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   6 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   6 nodes to be scored    (100 eliminated genes)

    ## 
    ##   Level 7:   14 nodes to be scored   (129 eliminated genes)

    ## 
    ##   Level 6:   21 nodes to be scored   (165 eliminated genes)

    ## 
    ##   Level 5:   28 nodes to be scored   (425 eliminated genes)

    ## 
    ##   Level 4:   21 nodes to be scored   (663 eliminated genes)

    ## 
    ##   Level 3:   17 nodes to be scored   (1596 eliminated genes)

    ## 
    ##   Level 2:   7 nodes to be scored    (1739 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2083 eliminated genes)

    ##         GO.ID                                               Term Annotated
    ## 1  GO:0001843                                neural tube closure         7
    ## 2  GO:0001841                              neural tube formation         8
    ## 3  GO:0001662                           behavioral fear response         1
    ## 4  GO:0002181                            cytoplasmic translation        21
    ## 5  GO:0000054              ribosomal subunit export from nucleus         3
    ## 6  GO:0003746             translation elongation factor activity         7
    ## 7  GO:0004561             alpha-N-acetylglucosaminidase activity         1
    ## 8  GO:0004568                                 chitinase activity         1
    ## 9  GO:0003724                              RNA helicase activity         1
    ## 10 GO:0008641 ubiquitin-like modifier activating enzyme activity         1
    ## 11 GO:0004441   inositol-1,4-bisphosphate 1-phosphatase activity         1
    ## 12 GO:0004305                       ethanolamine kinase activity         1
    ## 13 GO:0004438        phosphatidylinositol-3-phosphatase activity         2
    ## 14 GO:0003727                        single-stranded RNA binding         2
    ## 15 GO:0032450                 maltose alpha-glucosidase activity         2
    ## 16 GO:0004185              serine-type carboxypeptidase activity         2
    ## 17 GO:0005198                       structural molecule activity        54
    ##    Significant Expected Fisher               type
    ## 1            2     0.10 0.0042 Biological.Process
    ## 2            3     0.12 0.0135 Biological.Process
    ## 3            1     0.01 0.0149 Biological.Process
    ## 4            2     0.31 0.0374 Biological.Process
    ## 5            1     0.04 0.0440 Biological.Process
    ## 6            2     0.09 0.0035 Molecular.Function
    ## 7            1     0.01 0.0133 Molecular.Function
    ## 8            1     0.01 0.0133 Molecular.Function
    ## 9            1     0.01 0.0133 Molecular.Function
    ## 10           1     0.01 0.0133 Molecular.Function
    ## 11           1     0.01 0.0133 Molecular.Function
    ## 12           1     0.01 0.0133 Molecular.Function
    ## 13           1     0.03 0.0265 Molecular.Function
    ## 14           1     0.03 0.0265 Molecular.Function
    ## 15           1     0.03 0.0265 Molecular.Function
    ## 16           1     0.03 0.0265 Molecular.Function
    ## 17           4     0.72 0.0366 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 17 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   2 nodes to be scored    (263 eliminated genes)

    ## 
    ##   Level 6:   2 nodes to be scored    (267 eliminated genes)

    ## 
    ##   Level 5:   2 nodes to be scored    (275 eliminated genes)

    ## 
    ##   Level 4:   2 nodes to be scored    (277 eliminated genes)

    ## 
    ##   Level 3:   5 nodes to be scored    (293 eliminated genes)

    ## 
    ##   Level 2:   1 nodes to be scored    (305 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1010 eliminated genes)

    ## [1] GO.ID       Term        Annotated   Significant Expected    Fisher     
    ## <0 rows> (or 0-length row.names)
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 7 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 6:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 5:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 4:   1 nodes to be scored    (28 eliminated genes)

    ## 
    ##   Level 3:   2 nodes to be scored    (35 eliminated genes)

    ## 
    ##   Level 2:   1 nodes to be scored    (146 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (527 eliminated genes)

    ##        GO.ID                      Term Annotated Significant Expected Fisher
    ## 1 GO:0004180 carboxypeptidase activity        28           1     0.01   0.01
    ##                 type
    ## 1 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 36 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   4 nodes to be scored    (9 eliminated genes)

    ## 
    ##   Level 8:   2 nodes to be scored    (9 eliminated genes)

    ## 
    ##   Level 7:   2 nodes to be scored    (9 eliminated genes)

    ## 
    ##   Level 6:   4 nodes to be scored    (9 eliminated genes)

    ## 
    ##   Level 5:   4 nodes to be scored    (21 eliminated genes)

    ## 
    ##   Level 4:   8 nodes to be scored    (22 eliminated genes)

    ## 
    ##   Level 3:   6 nodes to be scored    (307 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (443 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (532 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 30 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   3 nodes to be scored    (263 eliminated genes)

    ## 
    ##   Level 6:   3 nodes to be scored    (267 eliminated genes)

    ## 
    ##   Level 5:   4 nodes to be scored    (397 eliminated genes)

    ## 
    ##   Level 4:   6 nodes to be scored    (399 eliminated genes)

    ## 
    ##   Level 3:   8 nodes to be scored    (462 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (1094 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1588 eliminated genes)

    ##        GO.ID                                       Term Annotated Significant
    ## 1 GO:0015969 guanosine tetraphosphate metabolic process         9           1
    ## 2 GO:0019700      organic phosphonate catabolic process        12           1
    ##   Expected Fisher               type
    ## 1     0.01  0.013 Biological.Process
    ## 2     0.02  0.017 Biological.Process
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 13 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   2 nodes to be scored    (9 eliminated genes)

    ## 
    ##   Level 5:   3 nodes to be scored    (35 eliminated genes)

    ## 
    ##   Level 4:   2 nodes to be scored    (40 eliminated genes)

    ## 
    ##   Level 3:   1 nodes to be scored    (145 eliminated genes)

    ## 
    ##   Level 2:   1 nodes to be scored    (149 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (186 eliminated genes)

    ##        GO.ID                                          Term Annotated
    ## 1 GO:0005302 L-tyrosine transmembrane transporter activity         9
    ##   Significant Expected Fisher               type
    ## 1           1        0 0.0032 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 10 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 7:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 5:   1 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 4:   1 nodes to be scored    (6 eliminated genes)

    ## 
    ##   Level 3:   3 nodes to be scored    (32 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (37 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (532 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 14 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 6:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 5:   3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 4:   3 nodes to be scored    (67 eliminated genes)

    ## 
    ##   Level 3:   3 nodes to be scored    (204 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (278 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (982 eliminated genes)

    ##        GO.ID                            Term Annotated Significant Expected
    ## 1 GO:0002933             lipid hydroxylation         3           1     0.00
    ## 2 GO:0002161 aminoacyl-tRNA editing activity         1           1     0.00
    ## 3 GO:0005509             calcium ion binding        66           1     0.05
    ##    Fisher               type
    ## 1 0.00210 Biological.Process
    ## 2 0.00072 Molecular.Function
    ## 3 0.04699 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 37 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   4 nodes to be scored    (10 eliminated genes)

    ## 
    ##   Level 6:   4 nodes to be scored    (10 eliminated genes)

    ## 
    ##   Level 5:   4 nodes to be scored    (64 eliminated genes)

    ## 
    ##   Level 4:   6 nodes to be scored    (64 eliminated genes)

    ## 
    ##   Level 3:   9 nodes to be scored    (198 eliminated genes)

    ## 
    ##   Level 2:   6 nodes to be scored    (413 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (609 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 21 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 7:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 5:   4 nodes to be scored    (18 eliminated genes)

    ## 
    ##   Level 4:   4 nodes to be scored    (70 eliminated genes)

    ## 
    ##   Level 3:   6 nodes to be scored    (554 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (840 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1612 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0001581 detection of chemical stimulus involved in sensory ...        10
    ## 2 GO:0002221         pattern recognition receptor signaling pathway        11
    ## 3 GO:0001965                        G-protein alpha-subunit binding        16
    ## 4 GO:0005507                                     copper ion binding        18
    ## 5 GO:0004197                   cysteine-type endopeptidase activity        18
    ##   Significant Expected Fisher               type
    ## 1           1     0.02  0.021 Biological.Process
    ## 2           1     0.02  0.023 Biological.Process
    ## 3           1     0.04  0.040 Molecular.Function
    ## 4           1     0.05  0.045 Molecular.Function
    ## 5           1     0.05  0.045 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 70 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  3 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 9:   4 nodes to be scored    (142 eliminated genes)

    ## 
    ##   Level 8:   3 nodes to be scored    (172 eliminated genes)

    ## 
    ##   Level 7:   5 nodes to be scored    (180 eliminated genes)

    ## 
    ##   Level 6:   11 nodes to be scored   (194 eliminated genes)

    ## 
    ##   Level 5:   16 nodes to be scored   (283 eliminated genes)

    ## 
    ##   Level 4:   11 nodes to be scored   (378 eliminated genes)

    ## 
    ##   Level 3:   6 nodes to be scored    (440 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (540 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (660 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ##        GO.ID                       Term Annotated Significant Expected Fisher
    ## 1 GO:0000209 protein polyubiquitination        30           1     0.04  0.042
    ##                 type
    ## 1 Biological.Process
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 11 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 6:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 5:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 4:   2 nodes to be scored    (68 eliminated genes)

    ## 
    ##   Level 3:   2 nodes to be scored    (185 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (186 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (523 eliminated genes)

    ##        GO.ID                           Term Annotated Significant Expected
    ## 1 GO:0004089 carbonate dehydratase activity         2           1     0.00
    ## 2 GO:0005509            calcium ion binding        66           1     0.05
    ##   Fisher               type
    ## 1 0.0014 Molecular.Function
    ## 2 0.0470 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 128 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  7 nodes to be scored    (17 eliminated genes)

    ## 
    ##   Level 10:  7 nodes to be scored    (69 eliminated genes)

    ## 
    ##   Level 9:   10 nodes to be scored   (86 eliminated genes)

    ## 
    ##   Level 8:   8 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 7:   8 nodes to be scored    (158 eliminated genes)

    ## 
    ##   Level 6:   18 nodes to be scored   (165 eliminated genes)

    ## 
    ##   Level 5:   26 nodes to be scored   (338 eliminated genes)

    ## 
    ##   Level 4:   19 nodes to be scored   (460 eliminated genes)

    ## 
    ##   Level 3:   11 nodes to be scored   (735 eliminated genes)

    ## 
    ##   Level 2:   6 nodes to be scored    (876 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1000 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 70 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   7 nodes to be scored    (2 eliminated genes)

    ## 
    ##   Level 6:   14 nodes to be scored   (198 eliminated genes)

    ## 
    ##   Level 5:   13 nodes to be scored   (332 eliminated genes)

    ## 
    ##   Level 4:   14 nodes to be scored   (433 eliminated genes)

    ## 
    ##   Level 3:   12 nodes to be scored   (923 eliminated genes)

    ## 
    ##   Level 2:   5 nodes to be scored    (1330 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1991 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0000454             snoRNA guided rRNA pseudouridine synthesis         1
    ## 2  GO:0000971 tRNA exon ligation utilizing 2',3' cyclic phosphate...         1
    ## 3  GO:0009058                                   biosynthetic process       207
    ## 4  GO:0060271                                        cilium assembly         3
    ## 5  GO:0000463 maturation of LSU-rRNA from tricistronic rRNA trans...         3
    ## 6  GO:0001558                              regulation of cell growth         3
    ## 7  GO:0001508                                       action potential         4
    ## 8  GO:0003735                     structural constituent of ribosome        13
    ## 9  GO:0003944 N-acetylglucosamine-1-phosphodiester alpha-N-acetyl...         1
    ## 10 GO:0000009                 alpha-1,6-mannosyltransferase activity         1
    ## 11 GO:0005536                                        glucose binding         1
    ## 12 GO:0004568                                     chitinase activity         1
    ## 13 GO:0005315 inorganic phosphate transmembrane transporter activ...         1
    ## 14 GO:0000224 peptide-N4-(N-acetyl-beta-glucosaminyl)asparagine a...         2
    ## 15 GO:0032450                     maltose alpha-glucosidase activity         2
    ## 16 GO:0003727                            single-stranded RNA binding         2
    ## 17 GO:0005539                              glycosaminoglycan binding         5
    ##    Significant Expected Fisher               type
    ## 1            1     0.01 0.0099 Biological.Process
    ## 2            1     0.01 0.0099 Biological.Process
    ## 3            2     2.05 0.0108 Biological.Process
    ## 4            1     0.03 0.0295 Biological.Process
    ## 5            1     0.03 0.0295 Biological.Process
    ## 6            1     0.03 0.0295 Biological.Process
    ## 7            1     0.04 0.0391 Biological.Process
    ## 8            2     0.12 0.0057 Molecular.Function
    ## 9            1     0.01 0.0090 Molecular.Function
    ## 10           1     0.01 0.0090 Molecular.Function
    ## 11           1     0.01 0.0090 Molecular.Function
    ## 12           1     0.01 0.0090 Molecular.Function
    ## 13           1     0.01 0.0090 Molecular.Function
    ## 14           1     0.02 0.0179 Molecular.Function
    ## 15           1     0.02 0.0179 Molecular.Function
    ## 16           1     0.02 0.0179 Molecular.Function
    ## 17           1     0.05 0.0443 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 16 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   1 nodes to be scored    (9 eliminated genes)

    ## 
    ##   Level 6:   1 nodes to be scored    (17 eliminated genes)

    ## 
    ##   Level 5:   2 nodes to be scored    (22 eliminated genes)

    ## 
    ##   Level 4:   3 nodes to be scored    (26 eliminated genes)

    ## 
    ##   Level 3:   4 nodes to be scored    (41 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (52 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (532 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 74 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   3 nodes to be scored    (107 eliminated genes)

    ## 
    ##   Level 7:   7 nodes to be scored    (189 eliminated genes)

    ## 
    ##   Level 6:   11 nodes to be scored   (204 eliminated genes)

    ## 
    ##   Level 5:   12 nodes to be scored   (394 eliminated genes)

    ## 
    ##   Level 4:   16 nodes to be scored   (474 eliminated genes)

    ## 
    ##   Level 3:   13 nodes to be scored   (738 eliminated genes)

    ## 
    ##   Level 2:   5 nodes to be scored    (1435 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2157 eliminated genes)

    ##        GO.ID                                         Term Annotated Significant
    ## 1 GO:0008218                              bioluminescence        10           3
    ## 2 GO:0000038 very long-chain fatty acid metabolic process         9           1
    ## 3 GO:0004089               carbonate dehydratase activity         2           1
    ## 4 GO:0003746       translation elongation factor activity         7           1
    ## 5 GO:0005245       voltage-gated calcium channel activity         7           1
    ## 6 GO:0005542                           folic acid binding        10           1
    ##   Expected   Fisher               type
    ## 1     0.03 0.000001 Biological.Process
    ## 2     0.03 0.025000 Biological.Process
    ## 3     0.01 0.008600 Molecular.Function
    ## 4     0.03 0.029900 Molecular.Function
    ## 5     0.03 0.029900 Molecular.Function
    ## 6     0.04 0.042500 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## [1] GO.ID       Term        Annotated   Significant Expected    Fisher     
    ## <0 rows> (or 0-length row.names)
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 20 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   1 nodes to be scored    (2 eliminated genes)

    ## 
    ##   Level 6:   3 nodes to be scored    (2 eliminated genes)

    ## 
    ##   Level 5:   4 nodes to be scored    (2 eliminated genes)

    ## 
    ##   Level 4:   4 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 3:   3 nodes to be scored    (19 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (199 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (277 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 17 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   3 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 8:   2 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 7:   1 nodes to be scored    (112 eliminated genes)

    ## 
    ##   Level 6:   2 nodes to be scored    (133 eliminated genes)

    ## 
    ##   Level 5:   2 nodes to be scored    (133 eliminated genes)

    ## 
    ##   Level 4:   1 nodes to be scored    (136 eliminated genes)

    ## 
    ##   Level 3:   2 nodes to be scored    (264 eliminated genes)

    ## 
    ##   Level 2:   1 nodes to be scored    (497 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (809 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0001755                            neural crest cell migration         2
    ## 2 GO:0001002 RNA polymerase III type 1 promoter sequence-specifi...         1
    ##   Significant Expected  Fisher               type
    ## 1           1        0 0.00140 Biological.Process
    ## 2           1        0 0.00036 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 51 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  2 nodes to be scored    (2 eliminated genes)

    ## 
    ##   Level 9:   4 nodes to be scored    (2 eliminated genes)

    ## 
    ##   Level 8:   5 nodes to be scored    (2 eliminated genes)

    ## 
    ##   Level 7:   7 nodes to be scored    (14 eliminated genes)

    ## 
    ##   Level 6:   6 nodes to be scored    (36 eliminated genes)

    ## 
    ##   Level 5:   8 nodes to be scored    (64 eliminated genes)

    ## 
    ##   Level 4:   7 nodes to be scored    (88 eliminated genes)

    ## 
    ##   Level 3:   6 nodes to be scored    (166 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (362 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (614 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 49 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   2 nodes to be scored    (7 eliminated genes)

    ## 
    ##   Level 7:   4 nodes to be scored    (69 eliminated genes)

    ## 
    ##   Level 6:   6 nodes to be scored    (74 eliminated genes)

    ## 
    ##   Level 5:   7 nodes to be scored    (224 eliminated genes)

    ## 
    ##   Level 4:   11 nodes to be scored   (296 eliminated genes)

    ## 
    ##   Level 3:   12 nodes to be scored   (411 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (856 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2047 eliminated genes)

    ##        GO.ID                                          Term Annotated
    ## 1 GO:0001516            prostaglandin biosynthetic process         2
    ## 2 GO:0000028              ribosomal small subunit assembly         3
    ## 3 GO:0001523                    retinoid metabolic process         9
    ## 4 GO:0004252            serine-type endopeptidase activity        53
    ## 5 GO:0005245        voltage-gated calcium channel activity         7
    ## 6 GO:0005542                            folic acid binding        10
    ## 7 GO:0047869 dimethylpropiothetin dethiomethylase activity        19
    ##   Significant Expected Fisher               type
    ## 1           1     0.00 0.0042 Biological.Process
    ## 2           1     0.01 0.0064 Biological.Process
    ## 3           1     0.02 0.0190 Biological.Process
    ## 4           2     0.11 0.0051 Molecular.Function
    ## 5           1     0.02 0.0150 Molecular.Function
    ## 6           1     0.02 0.0214 Molecular.Function
    ## 7           1     0.04 0.0404 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 2 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 2:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (0 eliminated genes)

    ## [1] GO.ID       Term        Annotated   Significant Expected    Fisher     
    ## <0 rows> (or 0-length row.names)
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## [1] GO.ID       Term        Annotated   Significant Expected    Fisher     
    ## <0 rows> (or 0-length row.names)
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 42 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   2 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 8:   1 nodes to be scored    (7 eliminated genes)

    ## 
    ##   Level 7:   3 nodes to be scored    (21 eliminated genes)

    ## 
    ##   Level 6:   6 nodes to be scored    (24 eliminated genes)

    ## 
    ##   Level 5:   11 nodes to be scored   (47 eliminated genes)

    ## 
    ##   Level 4:   8 nodes to be scored    (76 eliminated genes)

    ## 
    ##   Level 3:   6 nodes to be scored    (468 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (576 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (614 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 46 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   2 nodes to be scored    (100 eliminated genes)

    ## 
    ##   Level 7:   3 nodes to be scored    (383 eliminated genes)

    ## 
    ##   Level 6:   7 nodes to be scored    (397 eliminated genes)

    ## 
    ##   Level 5:   9 nodes to be scored    (408 eliminated genes)

    ## 
    ##   Level 4:   9 nodes to be scored    (534 eliminated genes)

    ## 
    ##   Level 3:   9 nodes to be scored    (866 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (1184 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1706 eliminated genes)

    ##        GO.ID                                           Term Annotated
    ## 1 GO:0002933                            lipid hydroxylation         3
    ## 2 GO:0001731 formation of translation preinitiation complex         1
    ## 3 GO:0004089                 carbonate dehydratase activity         2
    ## 4 GO:0047869  dimethylpropiothetin dethiomethylase activity        19
    ##   Significant Expected  Fisher               type
    ## 1           2     0.01 9.0e-06 Biological.Process
    ## 2           1     0.00 2.1e-03 Biological.Process
    ## 3           1     0.01 5.0e-03 Molecular.Function
    ## 4           1     0.05 4.7e-02 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## [1] GO.ID       Term        Annotated   Significant Expected    Fisher     
    ## <0 rows> (or 0-length row.names)
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 33 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   4 nodes to be scored    (9 eliminated genes)

    ## 
    ##   Level 6:   4 nodes to be scored    (31 eliminated genes)

    ## 
    ##   Level 5:   5 nodes to be scored    (49 eliminated genes)

    ## 
    ##   Level 4:   4 nodes to be scored    (49 eliminated genes)

    ## 
    ##   Level 3:   5 nodes to be scored    (256 eliminated genes)

    ## 
    ##   Level 2:   6 nodes to be scored    (270 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (811 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 36 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   2 nodes to be scored    (263 eliminated genes)

    ## 
    ##   Level 6:   5 nodes to be scored    (267 eliminated genes)

    ## 
    ##   Level 5:   6 nodes to be scored    (275 eliminated genes)

    ## 
    ##   Level 4:   8 nodes to be scored    (353 eliminated genes)

    ## 
    ##   Level 3:   10 nodes to be scored   (426 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (973 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1620 eliminated genes)

    ##        GO.ID                                      Term Annotated Significant
    ## 1 GO:0001782                        B cell homeostasis         7           1
    ## 2 GO:0001889                         liver development         8           1
    ## 3 GO:0001523                retinoid metabolic process         9           1
    ## 4 GO:0001701            in utero embryonic development        15           1
    ## 5 GO:0005542                        folic acid binding        10           1
    ## 6 GO:0004609 phosphatidylserine decarboxylase activity        13           1
    ##   Expected Fisher               type
    ## 1     0.02  0.020 Biological.Process
    ## 2     0.02  0.023 Biological.Process
    ## 3     0.03  0.025 Biological.Process
    ## 4     0.04  0.042 Biological.Process
    ## 5     0.03  0.025 Molecular.Function
    ## 6     0.03  0.032 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 22 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   2 nodes to be scored    (263 eliminated genes)

    ## 
    ##   Level 6:   2 nodes to be scored    (267 eliminated genes)

    ## 
    ##   Level 5:   3 nodes to be scored    (275 eliminated genes)

    ## 
    ##   Level 4:   3 nodes to be scored    (277 eliminated genes)

    ## 
    ##   Level 3:   7 nodes to be scored    (328 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (340 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1559 eliminated genes)

    ##        GO.ID                                   Term Annotated Significant
    ## 1 GO:0004842 ubiquitin-protein transferase activity        35           1
    ##   Expected Fisher               type
    ## 1     0.03  0.025 Molecular.Function

``` r
head(results_high0.6_cor_targets)
```

    ##        GO.ID                                             Term Annotated
    ## 1 GO:0005290   L-histidine transmembrane transporter activity         9
    ## 2 GO:0001658 branching involved in ureteric bud morphogenesis         4
    ## 3 GO:0001945                         lymph vessel development         5
    ## 4 GO:0004252               serine-type endopeptidase activity        53
    ## 5 GO:0004771                         sterol esterase activity         5
    ## 6 GO:0004966                        galanin receptor activity        12
    ##   Significant Expected Fisher               type         miRNA
    ## 1           1     0.00 0.0032 Molecular.Function Cluster_10452
    ## 2           1     0.01 0.0057 Biological.Process Cluster_11565
    ## 3           1     0.01 0.0071 Biological.Process Cluster_11565
    ## 4           2     0.15 0.0093 Molecular.Function Cluster_11565
    ## 5           1     0.01 0.0143 Molecular.Function Cluster_11565
    ## 6           1     0.03 0.0341 Molecular.Function Cluster_11565

Save results

``` r
write.csv(results_high0.6_cor_targets, "../output/27-Apul-mRNA-miRNA-interactions-topGO/miRNA_high0.6_cor_targets_topGO_FE.csv")
```

# 6 FE of specific miRNA’s targets (high 0.7 cor targets)

Loop through all miRNA and run functional enrichment on the miRNA’s
highly correlated targets (PCC magnitude \> 0.7)

``` r
interacting_miRNAs_high0.7 <- unique(high0.7_cor_bind_FA$miRNA) %>% na.omit
results_high0.7_cor_targets <- NULL  # initialize empty df

for(miRNA in interacting_miRNAs_high0.7) {
  
  # Run topGO enrichment function
  miRNA_results <- miRNA_topGO_FE(miRNA, high0.7_cor_bind_FA)
  
  # Only keep results if not empty
  if (nrow(miRNA_results) > 0) {
    
    # Add the miRNA source column
    miRNA_results$miRNA <- miRNA

    # Bind to the accumulating results data frame
    results_high0.7_cor_targets <- rbind(results_high0.7_cor_targets, miRNA_results)
  }
}
```

    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 9 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   1 nodes to be scored    (5 eliminated genes)

    ## 
    ##   Level 5:   1 nodes to be scored    (47 eliminated genes)

    ## 
    ##   Level 4:   1 nodes to be scored    (58 eliminated genes)

    ## 
    ##   Level 3:   1 nodes to be scored    (192 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (218 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (246 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 21 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   3 nodes to be scored    (12 eliminated genes)

    ## 
    ##   Level 5:   5 nodes to be scored    (31 eliminated genes)

    ## 
    ##   Level 4:   5 nodes to be scored    (90 eliminated genes)

    ## 
    ##   Level 3:   3 nodes to be scored    (321 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (467 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (755 eliminated genes)

    ##        GO.ID                      Term Annotated Significant Expected Fisher
    ## 1 GO:0001945  lymph vessel development         5           1     0.00 0.0035
    ## 2 GO:0004771  sterol esterase activity         5           1     0.01 0.0072
    ## 3 GO:0004966 galanin receptor activity        12           1     0.02 0.0172
    ##                 type
    ## 1 Biological.Process
    ## 2 Molecular.Function
    ## 3 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 18 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 7:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 5:   2 nodes to be scored    (35 eliminated genes)

    ## 
    ##   Level 4:   2 nodes to be scored    (44 eliminated genes)

    ## 
    ##   Level 3:   5 nodes to be scored    (97 eliminated genes)

    ## 
    ##   Level 2:   6 nodes to be scored    (307 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (429 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 17 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 7:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 5:   3 nodes to be scored    (122 eliminated genes)

    ## 
    ##   Level 4:   4 nodes to be scored    (178 eliminated genes)

    ## 
    ##   Level 3:   4 nodes to be scored    (317 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (460 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (934 eliminated genes)

    ##        GO.ID                                 Term Annotated Significant
    ## 1 GO:0001649           osteoblast differentiation        26           1
    ## 2 GO:0000165                         MAPK cascade        35           1
    ## 3 GO:0009034               tryptophanase activity         2           1
    ## 4 GO:0003964 RNA-directed DNA polymerase activity       122           2
    ##   Expected Fisher               type
    ## 1     0.04 0.0370 Biological.Process
    ## 2     0.05 0.0490 Biological.Process
    ## 3     0.00 0.0029 Molecular.Function
    ## 4     0.18 0.0108 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 12 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 7:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 5:   1 nodes to be scored    (35 eliminated genes)

    ## 
    ##   Level 4:   1 nodes to be scored    (44 eliminated genes)

    ## 
    ##   Level 3:   3 nodes to be scored    (71 eliminated genes)

    ## 
    ##   Level 2:   4 nodes to be scored    (246 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (352 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 14 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 7:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 5:   2 nodes to be scored    (122 eliminated genes)

    ## 
    ##   Level 4:   3 nodes to be scored    (178 eliminated genes)

    ## 
    ##   Level 3:   3 nodes to be scored    (315 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (438 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (877 eliminated genes)

    ##        GO.ID                                 Term Annotated Significant
    ## 1 GO:0000165                         MAPK cascade        35           1
    ## 2 GO:0003964 RNA-directed DNA polymerase activity       122           2
    ##   Expected Fisher               type
    ## 1     0.02 0.0250 Biological.Process
    ## 2     0.13 0.0056 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 12 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 7:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 5:   1 nodes to be scored    (35 eliminated genes)

    ## 
    ##   Level 4:   1 nodes to be scored    (44 eliminated genes)

    ## 
    ##   Level 3:   3 nodes to be scored    (71 eliminated genes)

    ## 
    ##   Level 2:   4 nodes to be scored    (246 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (352 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 14 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 7:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 5:   2 nodes to be scored    (122 eliminated genes)

    ## 
    ##   Level 4:   3 nodes to be scored    (178 eliminated genes)

    ## 
    ##   Level 3:   3 nodes to be scored    (315 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (438 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (877 eliminated genes)

    ##        GO.ID                                 Term Annotated Significant
    ## 1 GO:0000165                         MAPK cascade        35           1
    ## 2 GO:0003964 RNA-directed DNA polymerase activity       122           2
    ##   Expected Fisher               type
    ## 1     0.02 0.0250 Biological.Process
    ## 2     0.13 0.0056 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## [1] GO.ID       Term        Annotated   Significant Expected    Fisher     
    ## <0 rows> (or 0-length row.names)
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 36 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   3 nodes to be scored    (263 eliminated genes)

    ## 
    ##   Level 6:   4 nodes to be scored    (268 eliminated genes)

    ## 
    ##   Level 5:   6 nodes to be scored    (276 eliminated genes)

    ## 
    ##   Level 4:   6 nodes to be scored    (279 eliminated genes)

    ## 
    ##   Level 3:   9 nodes to be scored    (306 eliminated genes)

    ## 
    ##   Level 2:   4 nodes to be scored    (475 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1833 eliminated genes)

    ##        GO.ID                                     Term Annotated Significant
    ## 1 GO:0008137 NADH dehydrogenase (ubiquinone) activity         1           1
    ## 2 GO:0004176         ATP-dependent peptidase activity         2           1
    ##   Expected Fisher               type
    ## 1        0 0.0011 Molecular.Function
    ## 2        0 0.0022 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## [1] GO.ID       Term        Annotated   Significant Expected    Fisher     
    ## <0 rows> (or 0-length row.names)
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## [1] GO.ID       Term        Annotated   Significant Expected    Fisher     
    ## <0 rows> (or 0-length row.names)
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 17 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   2 nodes to be scored    (5 eliminated genes)

    ## 
    ##   Level 5:   2 nodes to be scored    (71 eliminated genes)

    ## 
    ##   Level 4:   2 nodes to be scored    (82 eliminated genes)

    ## 
    ##   Level 3:   3 nodes to be scored    (221 eliminated genes)

    ## 
    ##   Level 2:   4 nodes to be scored    (247 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (290 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 11 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 5:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 4:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 3:   3 nodes to be scored    (455 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (721 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1042 eliminated genes)

    ##        GO.ID                                            Term Annotated
    ## 1 GO:0001945                        lymph vessel development         5
    ## 2 GO:0002175 protein localization to paranode region of axon        24
    ##   Significant Expected Fisher               type
    ## 1           1     0.01 0.0071 Biological.Process
    ## 2           1     0.03 0.0337 Biological.Process
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 29 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   4 nodes to be scored    (18 eliminated genes)

    ## 
    ##   Level 5:   5 nodes to be scored    (73 eliminated genes)

    ## 
    ##   Level 4:   5 nodes to be scored    (240 eliminated genes)

    ## 
    ##   Level 3:   4 nodes to be scored    (371 eliminated genes)

    ## 
    ##   Level 2:   6 nodes to be scored    (517 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (628 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 34 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   3 nodes to be scored    (100 eliminated genes)

    ## 
    ##   Level 7:   4 nodes to be scored    (383 eliminated genes)

    ## 
    ##   Level 6:   5 nodes to be scored    (463 eliminated genes)

    ## 
    ##   Level 5:   5 nodes to be scored    (474 eliminated genes)

    ## 
    ##   Level 4:   4 nodes to be scored    (492 eliminated genes)

    ## 
    ##   Level 3:   6 nodes to be scored    (636 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (881 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1394 eliminated genes)

    ##        GO.ID                                       Term Annotated Significant
    ## 1 GO:0000165                               MAPK cascade        35           2
    ## 2 GO:0001818 negative regulation of cytokine production        18           1
    ## 3 GO:0003924                            GTPase activity        66           2
    ##   Expected Fisher               type
    ## 1     0.07 0.0018 Biological.Process
    ## 2     0.04 0.0378 Biological.Process
    ## 3     0.10 0.0032 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 7 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 5:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 4:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 3:   3 nodes to be scored    (293 eliminated genes)

    ## 
    ##   Level 2:   1 nodes to be scored    (293 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (814 eliminated genes)

    ## [1] GO.ID       Term        Annotated   Significant Expected    Fisher     
    ## <0 rows> (or 0-length row.names)
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## [1] GO.ID       Term        Annotated   Significant Expected    Fisher     
    ## <0 rows> (or 0-length row.names)
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 21 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   2 nodes to be scored    (65 eliminated genes)

    ## 
    ##   Level 5:   2 nodes to be scored    (100 eliminated genes)

    ## 
    ##   Level 4:   3 nodes to be scored    (109 eliminated genes)

    ## 
    ##   Level 3:   4 nodes to be scored    (263 eliminated genes)

    ## 
    ##   Level 2:   6 nodes to be scored    (467 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (598 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 9 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 6:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 5:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 4:   2 nodes to be scored    (53 eliminated genes)

    ## 
    ##   Level 3:   2 nodes to be scored    (105 eliminated genes)

    ## 
    ##   Level 2:   1 nodes to be scored    (146 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (527 eliminated genes)

    ##        GO.ID                               Term Annotated Significant Expected
    ## 1 GO:0001822                 kidney development        65           2     0.14
    ## 2 GO:0004252 serine-type endopeptidase activity        53           1     0.02
    ##   Fisher               type
    ## 1 0.0061 Biological.Process
    ## 2 0.0190 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 37 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  2 nodes to be scored    (63 eliminated genes)

    ## 
    ##   Level 9:   3 nodes to be scored    (63 eliminated genes)

    ## 
    ##   Level 8:   3 nodes to be scored    (101 eliminated genes)

    ## 
    ##   Level 7:   2 nodes to be scored    (106 eliminated genes)

    ## 
    ##   Level 6:   3 nodes to be scored    (139 eliminated genes)

    ## 
    ##   Level 5:   5 nodes to be scored    (289 eliminated genes)

    ## 
    ##   Level 4:   7 nodes to be scored    (338 eliminated genes)

    ## 
    ##   Level 3:   5 nodes to be scored    (527 eliminated genes)

    ## 
    ##   Level 2:   4 nodes to be scored    (679 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (778 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 21 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   1 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 6:   3 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 5:   5 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 4:   3 nodes to be scored    (33 eliminated genes)

    ## 
    ##   Level 3:   3 nodes to be scored    (196 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (363 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (807 eliminated genes)

    ##        GO.ID                                             Term Annotated
    ## 1 GO:0004305                     ethanolamine kinase activity         1
    ## 2 GO:0004441 inositol-1,4-bisphosphate 1-phosphatase activity         1
    ## 3 GO:0005085       guanyl-nucleotide exchange factor activity        12
    ##   Significant Expected Fisher               type
    ## 1           1     0.00 0.0018 Molecular.Function
    ## 2           1     0.00 0.0018 Molecular.Function
    ## 3           1     0.02 0.0214 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## [1] GO.ID       Term        Annotated   Significant Expected    Fisher     
    ## <0 rows> (or 0-length row.names)
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 12 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 7:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 5:   1 nodes to be scored    (12 eliminated genes)

    ## 
    ##   Level 4:   2 nodes to be scored    (13 eliminated genes)

    ## 
    ##   Level 3:   3 nodes to be scored    (23 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (84 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (532 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 12 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 7:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 5:   2 nodes to be scored    (122 eliminated genes)

    ## 
    ##   Level 4:   3 nodes to be scored    (122 eliminated genes)

    ## 
    ##   Level 3:   3 nodes to be scored    (169 eliminated genes)

    ## 
    ##   Level 2:   1 nodes to be scored    (292 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (578 eliminated genes)

    ##        GO.ID                                   Term Annotated Significant
    ## 1 GO:0019700  organic phosphonate catabolic process        12           1
    ## 2 GO:0004842 ubiquitin-protein transferase activity        35           1
    ##   Expected Fisher               type
    ## 1     0.01 0.0085 Biological.Process
    ## 2     0.03 0.0250 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 7 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 7:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 5:   1 nodes to be scored    (18 eliminated genes)

    ## 
    ##   Level 4:   1 nodes to be scored    (52 eliminated genes)

    ## 
    ##   Level 3:   1 nodes to be scored    (181 eliminated genes)

    ## 
    ##   Level 2:   1 nodes to be scored    (181 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (465 eliminated genes)

    ##        GO.ID               Term Annotated Significant Expected Fisher
    ## 1 GO:0005507 copper ion binding        18           1     0.01 0.0065
    ##                 type
    ## 1 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 27 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  2 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 10:  2 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 9:   3 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 8:   2 nodes to be scored    (6 eliminated genes)

    ## 
    ##   Level 7:   1 nodes to be scored    (82 eliminated genes)

    ## 
    ##   Level 6:   1 nodes to be scored    (83 eliminated genes)

    ## 
    ##   Level 5:   2 nodes to be scored    (242 eliminated genes)

    ## 
    ##   Level 4:   5 nodes to be scored    (268 eliminated genes)

    ## 
    ##   Level 3:   4 nodes to be scored    (326 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (447 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (532 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0000971 tRNA exon ligation utilizing 2',3' cyclic phosphate...         1
    ##   Significant Expected  Fisher               type
    ## 1           1        0 0.00071 Biological.Process
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 5 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 4:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 3:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (10 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (443 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 20 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 7:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 5:   3 nodes to be scored    (152 eliminated genes)

    ## 
    ##   Level 4:   5 nodes to be scored    (156 eliminated genes)

    ## 
    ##   Level 3:   4 nodes to be scored    (170 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (296 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (655 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0008218                                        bioluminescence        10
    ## 2 GO:0004089                         carbonate dehydratase activity         2
    ## 3 GO:0001409 guanine nucleotide transmembrane transporter activi...        30
    ##   Significant Expected   Fisher               type
    ## 1           2     0.01 0.000045 Biological.Process
    ## 2           1     0.00 0.002200 Molecular.Function
    ## 3           1     0.03 0.032100 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 12 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   1 nodes to be scored    (9 eliminated genes)

    ## 
    ##   Level 6:   1 nodes to be scored    (9 eliminated genes)

    ## 
    ##   Level 5:   1 nodes to be scored    (9 eliminated genes)

    ## 
    ##   Level 4:   1 nodes to be scored    (9 eliminated genes)

    ## 
    ##   Level 3:   3 nodes to be scored    (32 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (37 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (532 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 49 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   2 nodes to be scored    (7 eliminated genes)

    ## 
    ##   Level 7:   4 nodes to be scored    (69 eliminated genes)

    ## 
    ##   Level 6:   6 nodes to be scored    (74 eliminated genes)

    ## 
    ##   Level 5:   7 nodes to be scored    (224 eliminated genes)

    ## 
    ##   Level 4:   11 nodes to be scored   (296 eliminated genes)

    ## 
    ##   Level 3:   12 nodes to be scored   (411 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (856 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2047 eliminated genes)

    ##        GO.ID                                          Term Annotated
    ## 1 GO:0001523                    retinoid metabolic process         9
    ## 2 GO:0005245        voltage-gated calcium channel activity         7
    ## 3 GO:0005542                            folic acid binding        10
    ## 4 GO:0047869 dimethylpropiothetin dethiomethylase activity        19
    ##   Significant Expected Fisher               type
    ## 1           1     0.01 0.0064 Biological.Process
    ## 2           1     0.01 0.0130 Molecular.Function
    ## 3           1     0.02 0.0180 Molecular.Function
    ## 4           1     0.03 0.0340 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 6 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 6:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 5:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 4:   1 nodes to be scored    (66 eliminated genes)

    ## 
    ##   Level 3:   1 nodes to be scored    (181 eliminated genes)

    ## 
    ##   Level 2:   1 nodes to be scored    (181 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (465 eliminated genes)

    ##        GO.ID                Term Annotated Significant Expected Fisher
    ## 1 GO:0005509 calcium ion binding        66           1     0.02  0.024
    ##                 type
    ## 1 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 9 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   1 nodes to be scored    (15 eliminated genes)

    ## 
    ##   Level 5:   1 nodes to be scored    (25 eliminated genes)

    ## 
    ##   Level 4:   1 nodes to be scored    (25 eliminated genes)

    ## 
    ##   Level 3:   1 nodes to be scored    (38 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (218 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (246 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 6 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 5:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 4:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 3:   2 nodes to be scored    (103 eliminated genes)

    ## 
    ##   Level 2:   1 nodes to be scored    (146 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (527 eliminated genes)

    ##        GO.ID                           Term Annotated Significant Expected
    ## 1 GO:0001701 in utero embryonic development        15           1     0.01
    ## 2 GO:0004175         endopeptidase activity       103           1     0.04
    ##   Fisher               type
    ## 1  0.011 Biological.Process
    ## 2  0.037 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 6 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 5:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 4:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 3:   2 nodes to be scored    (35 eliminated genes)

    ## 
    ##   Level 2:   1 nodes to be scored    (35 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (549 eliminated genes)

    ##        GO.ID                                   Term Annotated Significant
    ## 1 GO:0004842 ubiquitin-protein transferase activity        35           1
    ##   Expected Fisher               type
    ## 1     0.01  0.013 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 272 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1229 GO terms and 2373 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1411 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## 
    ## Building most specific GOs .....

    ##  ( 460 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 909 GO terms and 1190 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2776 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 0 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## Warning in getSigGroups(object, test.stat): No enrichment can pe performed -
    ## there are no feasible GO terms!

    ## [1] GO.ID       Term        Annotated   Significant Expected    Fisher     
    ## <0 rows> (or 0-length row.names)

``` r
head(results_high0.7_cor_targets)
```

    ##        GO.ID                       Term Annotated Significant Expected Fisher
    ## 1 GO:0001945   lymph vessel development         5           1     0.00 0.0035
    ## 2 GO:0004771   sterol esterase activity         5           1     0.01 0.0072
    ## 3 GO:0004966  galanin receptor activity        12           1     0.02 0.0172
    ## 4 GO:0001649 osteoblast differentiation        26           1     0.04 0.0370
    ## 5 GO:0000165               MAPK cascade        35           1     0.05 0.0490
    ## 6 GO:0009034     tryptophanase activity         2           1     0.00 0.0029
    ##                 type         miRNA
    ## 1 Biological.Process Cluster_11565
    ## 2 Molecular.Function Cluster_11565
    ## 3 Molecular.Function Cluster_11565
    ## 4 Biological.Process Cluster_12081
    ## 5 Biological.Process Cluster_12081
    ## 6 Molecular.Function Cluster_12081

Save results

``` r
write.csv(results_high0.7_cor_targets, "../output/27-Apul-mRNA-miRNA-interactions-topGO/miRNA_high0.7_cor_targets_topGO_FE.csv")
```
