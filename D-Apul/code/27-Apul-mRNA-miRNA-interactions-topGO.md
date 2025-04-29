27-Apul-mRNA-miRNA-interactions-topGO
================
Kathleen Durkin
2025-04-28

- [1 Format topGO files](#1-format-topgo-files)
  - [1.1 Read in and format annotation
    files](#11-read-in-and-format-annotation-files)
  - [1.2 Set up gene2GO object](#12-set-up-gene2go-object)
  - [1.3 Define reference set](#13-define-reference-set)
  - [1.4 Read in PCC/miranda data](#14-read-in-pccmiranda-data)
- [2 FA of all miRNA targets](#2-fa-of-all-mirna-targets)
- [3 FE of specific miRNA’s targets (all
  targets)](#3-fe-of-specific-mirnas-targets-all-targets)
- [4 FE of specific miRNA’s targets (high cor
  targets)](#4-fe-of-specific-mirnas-targets-high-cor-targets)

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

Set function to select genes of interest (ie those that have pvalue \<
0.05)

``` r
topDiffGenes <- function(allScore) {
return(allScore < 0.05)}
```

# 2 FA of all miRNA targets

Functional annotation of all putative miRNA targets

``` r
cor_bind_FA <- left_join(data, annot, by = c("mRNA" = "gene_ID")) %>% distinct()

nrow(cor_bind_FA)
```

    ## [1] 122274

``` r
nrow(cor_bind_FA[!is.na(cor_bind_FA$Gene.Ontology.IDs),])
```

    ## [1] 40154

Of the 122274 putative miRNA targets predicted by miRanda with
significant PCC, 40154 have available annotations

``` r
high_cor_bind_FA <- cor_bind_FA[abs(cor_bind_FA$PCC.cor) > 0.5,]

nrow(high_cor_bind_FA)
```

    ## [1] 16179

``` r
nrow(high_cor_bind_FA[!is.na(high_cor_bind_FA$Gene.Ontology.IDs),])
```

    ## [1] 5307

Of the 16179 putative miRNA targets predicted by miRanda that that have
highly correlated expression (magnitude of correlation \> 0.5), 5307
have available annotations.

Save

``` r
write.csv(cor_bind_FA, "../output/27-Apul-mRNA-miRNA-interactions-topGO/miRNA_targets_FA.csv")
write.csv(high_cor_bind_FA, "../output/27-Apul-mRNA-miRNA-interactions-topGO/miRNA_high_cor_targets_FA.csv")
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

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 153 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  3 nodes to be scored    (10 eliminated genes)

    ## 
    ##   Level 9:   6 nodes to be scored    (15 eliminated genes)

    ## 
    ##   Level 8:   10 nodes to be scored   (33 eliminated genes)

    ## 
    ##   Level 7:   20 nodes to be scored   (55 eliminated genes)

    ## 
    ##   Level 6:   27 nodes to be scored   (259 eliminated genes)

    ## 
    ##   Level 5:   29 nodes to be scored   (564 eliminated genes)

    ## 
    ##   Level 4:   25 nodes to be scored   (657 eliminated genes)

    ## 
    ##   Level 3:   20 nodes to be scored   (904 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (1109 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1282 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 130 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   6 nodes to be scored    (100 eliminated genes)

    ## 
    ##   Level 7:   18 nodes to be scored   (383 eliminated genes)

    ## 
    ##   Level 6:   24 nodes to be scored   (518 eliminated genes)

    ## 
    ##   Level 5:   28 nodes to be scored   (993 eliminated genes)

    ## 
    ##   Level 4:   24 nodes to be scored   (1231 eliminated genes)

    ## 
    ##   Level 3:   19 nodes to be scored   (1776 eliminated genes)

    ## 
    ##   Level 2:   6 nodes to be scored    (2055 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2475 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0000724 double-strand break repair via homologous recombina...        13
    ## 2  GO:0000184 nuclear-transcribed mRNA catabolic process, nonsens...        10
    ## 3  GO:0000712      resolution of meiotic recombination intermediates         1
    ## 4  GO:0006813                                potassium ion transport         1
    ## 5  GO:0000423                                              mitophagy         2
    ## 6  GO:0002933                                    lipid hydroxylation         3
    ## 7  GO:0000028                       ribosomal small subunit assembly         3
    ## 8  GO:0005524                                            ATP binding       263
    ## 9  GO:0004013                        adenosylhomocysteinase activity         1
    ## 10 GO:0004108                         citrate (Si)-synthase activity         1
    ## 11 GO:0004089                         carbonate dehydratase activity         2
    ## 12 GO:0004517                         nitric-oxide synthase activity         2
    ## 13 GO:0004185                  serine-type carboxypeptidase activity         2
    ## 14 GO:0005201            extracellular matrix structural constituent        17
    ##    Significant Expected  Fisher               type
    ## 1            3     0.20 0.00086 Biological.Process
    ## 2            2     0.16 0.00978 Biological.Process
    ## 3            1     0.02 0.01567 Biological.Process
    ## 4            1     0.02 0.01567 Biological.Process
    ## 5            1     0.03 0.03110 Biological.Process
    ## 6            1     0.05 0.04631 Biological.Process
    ## 7            1     0.05 0.04631 Biological.Process
    ## 8           12     5.03 0.00320 Molecular.Function
    ## 9            1     0.02 0.01910 Molecular.Function
    ## 10           1     0.02 0.01910 Molecular.Function
    ## 11           1     0.04 0.03790 Molecular.Function
    ## 12           1     0.04 0.03790 Molecular.Function
    ## 13           1     0.04 0.03790 Molecular.Function
    ## 14           2     0.32 0.04060 Molecular.Function

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

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 276 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  3 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 11:  8 nodes to be scored    (19 eliminated genes)

    ## 
    ##   Level 10:  12 nodes to be scored   (74 eliminated genes)

    ## 
    ##   Level 9:   19 nodes to be scored   (253 eliminated genes)

    ## 
    ##   Level 8:   28 nodes to be scored   (292 eliminated genes)

    ## 
    ##   Level 7:   30 nodes to be scored   (321 eliminated genes)

    ## 
    ##   Level 6:   41 nodes to be scored   (449 eliminated genes)

    ## 
    ##   Level 5:   50 nodes to be scored   (683 eliminated genes)

    ## 
    ##   Level 4:   38 nodes to be scored   (820 eliminated genes)

    ## 
    ##   Level 3:   34 nodes to be scored   (1046 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (1225 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1381 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 184 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   7 nodes to be scored    (9 eliminated genes)

    ## 
    ##   Level 8:   14 nodes to be scored   (144 eliminated genes)

    ## 
    ##   Level 7:   28 nodes to be scored   (438 eliminated genes)

    ## 
    ##   Level 6:   36 nodes to be scored   (549 eliminated genes)

    ## 
    ##   Level 5:   36 nodes to be scored   (926 eliminated genes)

    ## 
    ##   Level 4:   34 nodes to be scored   (1155 eliminated genes)

    ## 
    ##   Level 3:   16 nodes to be scored   (1758 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (2043 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2457 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0000082                  G1/S transition of mitotic cell cycle       104
    ## 2  GO:0001913                           T cell mediated cytotoxicity         1
    ## 3  GO:0002230 positive regulation of defense response to virus by...         1
    ## 4  GO:0001172                            RNA-templated transcription         1
    ## 5  GO:0004176                       ATP-dependent peptidase activity         2
    ## 6  GO:0003964                   RNA-directed DNA polymerase activity       122
    ## 7  GO:0004129                          cytochrome-c oxidase activity         1
    ## 8  GO:0000254                      C-4 methylsterol oxidase activity         1
    ## 9  GO:0004334                           fumarylacetoacetase activity         1
    ## 10 GO:0004862       cAMP-dependent protein kinase inhibitor activity         1
    ## 11 GO:0008137               NADH dehydrogenase (ubiquinone) activity         1
    ## 12 GO:0036381 pyridoxal 5'-phosphate synthase (glutamine hydrolys...         1
    ##    Significant Expected  Fisher               type
    ## 1            9     2.67 0.00083 Biological.Process
    ## 2            1     0.03 0.02564 Biological.Process
    ## 3            1     0.03 0.02564 Biological.Process
    ## 4            1     0.03 0.02564 Biological.Process
    ## 5            2     0.06 0.00100 Molecular.Function
    ## 6           11     3.87 0.00140 Molecular.Function
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

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 275 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  4 nodes to be scored    (2 eliminated genes)

    ## 
    ##   Level 12:  5 nodes to be scored    (18 eliminated genes)

    ## 
    ##   Level 11:  8 nodes to be scored    (22 eliminated genes)

    ## 
    ##   Level 10:  15 nodes to be scored   (71 eliminated genes)

    ## 
    ##   Level 9:   21 nodes to be scored   (222 eliminated genes)

    ## 
    ##   Level 8:   29 nodes to be scored   (258 eliminated genes)

    ## 
    ##   Level 7:   36 nodes to be scored   (318 eliminated genes)

    ## 
    ##   Level 6:   43 nodes to be scored   (389 eliminated genes)

    ## 
    ##   Level 5:   44 nodes to be scored   (598 eliminated genes)

    ## 
    ##   Level 4:   32 nodes to be scored   (685 eliminated genes)

    ## 
    ##   Level 3:   24 nodes to be scored   (979 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (1228 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1344 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 173 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 9:   8 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 8:   14 nodes to be scored   (147 eliminated genes)

    ## 
    ##   Level 7:   25 nodes to be scored   (442 eliminated genes)

    ## 
    ##   Level 6:   30 nodes to be scored   (594 eliminated genes)

    ## 
    ##   Level 5:   33 nodes to be scored   (1020 eliminated genes)

    ## 
    ##   Level 4:   26 nodes to be scored   (1228 eliminated genes)

    ## 
    ##   Level 3:   20 nodes to be scored   (1863 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (2098 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2540 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0006886                        intracellular protein transport         4
    ## 2 GO:0002225 positive regulation of antimicrobial peptide produc...         1
    ## 3 GO:0000289       nuclear-transcribed mRNA poly(A) tail shortening         2
    ## 4 GO:0001578                           microtubule bundle formation         2
    ## 5 GO:0004467              long-chain fatty acid-CoA ligase activity         1
    ## 6 GO:0004997        thyrotropin-releasing hormone receptor activity         1
    ## 7 GO:0004335                                 galactokinase activity         1
    ## 8 GO:0005220 inositol 1,4,5-trisphosphate-gated calcium channel ...         1
    ## 9 GO:0003844             1,4-alpha-glucan branching enzyme activity         1
    ##   Significant Expected Fisher               type
    ## 1           2     0.08 0.0022 Biological.Process
    ## 2           1     0.02 0.0199 Biological.Process
    ## 3           1     0.04 0.0395 Biological.Process
    ## 4           1     0.04 0.0395 Biological.Process
    ## 5           1     0.03 0.0300 Molecular.Function
    ## 6           1     0.03 0.0300 Molecular.Function
    ## 7           1     0.03 0.0300 Molecular.Function
    ## 8           1     0.03 0.0300 Molecular.Function
    ## 9           1     0.03 0.0300 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 200 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  2 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 11:  3 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (70 eliminated genes)

    ## 
    ##   Level 9:   10 nodes to be scored   (101 eliminated genes)

    ## 
    ##   Level 8:   13 nodes to be scored   (145 eliminated genes)

    ## 
    ##   Level 7:   24 nodes to be scored   (200 eliminated genes)

    ## 
    ##   Level 6:   32 nodes to be scored   (437 eliminated genes)

    ## 
    ##   Level 5:   39 nodes to be scored   (667 eliminated genes)

    ## 
    ##   Level 4:   32 nodes to be scored   (744 eliminated genes)

    ## 
    ##   Level 3:   28 nodes to be scored   (978 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (1145 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1329 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 139 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   7 nodes to be scored    (100 eliminated genes)

    ## 
    ##   Level 7:   17 nodes to be scored   (384 eliminated genes)

    ## 
    ##   Level 6:   24 nodes to be scored   (469 eliminated genes)

    ## 
    ##   Level 5:   33 nodes to be scored   (854 eliminated genes)

    ## 
    ##   Level 4:   26 nodes to be scored   (1135 eliminated genes)

    ## 
    ##   Level 3:   16 nodes to be scored   (1717 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (2092 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2462 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0000349 generation of catalytic spliceosome for first trans...         1
    ## 2  GO:0000390                       spliceosomal complex disassembly         1
    ## 3  GO:0002520                              immune system development         1
    ## 4  GO:0001502                                 cartilage condensation         2
    ## 5  GO:0001895                                     retina homeostasis         2
    ## 6  GO:0001503                                           ossification        42
    ## 7  GO:0003964                   RNA-directed DNA polymerase activity       122
    ## 8  GO:0047869          dimethylpropiothetin dethiomethylase activity        19
    ## 9  GO:0000048                           peptidyltransferase activity         1
    ## 10 GO:0005436                    sodium:phosphate symporter activity         1
    ## 11 GO:0004013                        adenosylhomocysteinase activity         1
    ## 12 GO:0004095              carnitine O-palmitoyltransferase activity         1
    ## 13 GO:0004467              long-chain fatty acid-CoA ligase activity         1
    ## 14 GO:0000703 oxidized pyrimidine nucleobase lesion DNA N-glycosy...         1
    ##    Significant Expected  Fisher               type
    ## 1            1     0.02 0.02200 Biological.Process
    ## 2            1     0.02 0.02200 Biological.Process
    ## 3            1     0.02 0.02200 Biological.Process
    ## 4            1     0.04 0.04400 Biological.Process
    ## 5            1     0.04 0.04400 Biological.Process
    ## 6            3     0.93 0.04600 Biological.Process
    ## 7           11     3.65 0.00083 Molecular.Function
    ## 8            4     0.57 0.00205 Molecular.Function
    ## 9            1     0.03 0.02993 Molecular.Function
    ## 10           1     0.03 0.02993 Molecular.Function
    ## 11           1     0.03 0.02993 Molecular.Function
    ## 12           1     0.03 0.02993 Molecular.Function
    ## 13           1     0.03 0.02993 Molecular.Function
    ## 14           1     0.03 0.02993 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 220 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  2 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 11:  4 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 10:  7 nodes to be scored    (70 eliminated genes)

    ## 
    ##   Level 9:   14 nodes to be scored   (241 eliminated genes)

    ## 
    ##   Level 8:   17 nodes to be scored   (286 eliminated genes)

    ## 
    ##   Level 7:   24 nodes to be scored   (339 eliminated genes)

    ## 
    ##   Level 6:   33 nodes to be scored   (457 eliminated genes)

    ## 
    ##   Level 5:   43 nodes to be scored   (653 eliminated genes)

    ## 
    ##   Level 4:   34 nodes to be scored   (702 eliminated genes)

    ## 
    ##   Level 3:   29 nodes to be scored   (1023 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (1154 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1328 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 141 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   6 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   7 nodes to be scored    (100 eliminated genes)

    ## 
    ##   Level 7:   15 nodes to be scored   (386 eliminated genes)

    ## 
    ##   Level 6:   25 nodes to be scored   (469 eliminated genes)

    ## 
    ##   Level 5:   35 nodes to be scored   (831 eliminated genes)

    ## 
    ##   Level 4:   28 nodes to be scored   (1127 eliminated genes)

    ## 
    ##   Level 3:   15 nodes to be scored   (1654 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (2165 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2438 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0000349 generation of catalytic spliceosome for first trans...         1
    ## 2  GO:0002520                              immune system development         1
    ## 3  GO:0007035                                 vacuolar acidification         1
    ## 4  GO:0000390                       spliceosomal complex disassembly         1
    ## 5  GO:0001895                                     retina homeostasis         2
    ## 6  GO:0001503                                           ossification        42
    ## 7  GO:0003964                   RNA-directed DNA polymerase activity       122
    ## 8  GO:0000048                           peptidyltransferase activity         1
    ## 9  GO:0043531                                            ADP binding         1
    ## 10 GO:0004435          phosphatidylinositol phospholipase C activity         1
    ## 11 GO:0005436                    sodium:phosphate symporter activity         1
    ## 12 GO:0004013                        adenosylhomocysteinase activity         1
    ## 13 GO:0004467              long-chain fatty acid-CoA ligase activity         1
    ## 14 GO:0000703 oxidized pyrimidine nucleobase lesion DNA N-glycosy...         1
    ## 15 GO:0003824                                     catalytic activity      1038
    ## 16 GO:0009034                                 tryptophanase activity         2
    ## 17 GO:0016829                                         lyase activity        58
    ##    Significant Expected Fisher               type
    ## 1            1     0.02 0.0230 Biological.Process
    ## 2            1     0.02 0.0230 Biological.Process
    ## 3            1     0.02 0.0230 Biological.Process
    ## 4            1     0.02 0.0230 Biological.Process
    ## 5            1     0.05 0.0450 Biological.Process
    ## 6            3     0.96 0.0480 Biological.Process
    ## 7            9     3.04 0.0028 Molecular.Function
    ## 8            1     0.02 0.0249 Molecular.Function
    ## 9            1     0.02 0.0249 Molecular.Function
    ## 10           1     0.02 0.0249 Molecular.Function
    ## 11           1     0.02 0.0249 Molecular.Function
    ## 12           1     0.02 0.0249 Molecular.Function
    ## 13           1     0.02 0.0249 Molecular.Function
    ## 14           1     0.02 0.0249 Molecular.Function
    ## 15          31    25.83 0.0300 Molecular.Function
    ## 16           1     0.05 0.0492 Molecular.Function
    ## 17           5     1.44 0.0494 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 182 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  4 nodes to be scored    (16 eliminated genes)

    ## 
    ##   Level 10:  6 nodes to be scored    (70 eliminated genes)

    ## 
    ##   Level 9:   12 nodes to be scored   (241 eliminated genes)

    ## 
    ##   Level 8:   15 nodes to be scored   (285 eliminated genes)

    ## 
    ##   Level 7:   22 nodes to be scored   (333 eliminated genes)

    ## 
    ##   Level 6:   29 nodes to be scored   (427 eliminated genes)

    ## 
    ##   Level 5:   32 nodes to be scored   (657 eliminated genes)

    ## 
    ##   Level 4:   25 nodes to be scored   (708 eliminated genes)

    ## 
    ##   Level 3:   23 nodes to be scored   (913 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (1096 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1322 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

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
    ##   Level 8:   6 nodes to be scored    (100 eliminated genes)

    ## 
    ##   Level 7:   15 nodes to be scored   (384 eliminated genes)

    ## 
    ##   Level 6:   20 nodes to be scored   (492 eliminated genes)

    ## 
    ##   Level 5:   30 nodes to be scored   (865 eliminated genes)

    ## 
    ##   Level 4:   25 nodes to be scored   (1000 eliminated genes)

    ## 
    ##   Level 3:   17 nodes to be scored   (1590 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (1988 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2556 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0000390                       spliceosomal complex disassembly         1
    ## 2  GO:0002520                              immune system development         1
    ## 3  GO:0001895                                     retina homeostasis         2
    ## 4  GO:0001503                                           ossification        42
    ## 5  GO:0003964                   RNA-directed DNA polymerase activity       122
    ## 6  GO:0000993                      RNA polymerase II complex binding         1
    ## 7  GO:0043531                                            ADP binding         1
    ## 8  GO:0000703 oxidized pyrimidine nucleobase lesion DNA N-glycosy...         1
    ## 9  GO:0004467              long-chain fatty acid-CoA ligase activity         1
    ## 10 GO:0000048                           peptidyltransferase activity         1
    ## 11 GO:0004013                        adenosylhomocysteinase activity         1
    ##    Significant Expected Fisher               type
    ## 1            1     0.02 0.0210 Biological.Process
    ## 2            1     0.02 0.0210 Biological.Process
    ## 3            1     0.04 0.0420 Biological.Process
    ## 4            3     0.90 0.0430 Biological.Process
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

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 210 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  6 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  9 nodes to be scored    (87 eliminated genes)

    ## 
    ##   Level 9:   13 nodes to be scored   (262 eliminated genes)

    ## 
    ##   Level 8:   17 nodes to be scored   (319 eliminated genes)

    ## 
    ##   Level 7:   22 nodes to be scored   (340 eliminated genes)

    ## 
    ##   Level 6:   30 nodes to be scored   (412 eliminated genes)

    ## 
    ##   Level 5:   41 nodes to be scored   (624 eliminated genes)

    ## 
    ##   Level 4:   33 nodes to be scored   (712 eliminated genes)

    ## 
    ##   Level 3:   25 nodes to be scored   (1014 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (1105 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1373 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 148 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   5 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 8:   11 nodes to be scored   (101 eliminated genes)

    ## 
    ##   Level 7:   22 nodes to be scored   (383 eliminated genes)

    ## 
    ##   Level 6:   27 nodes to be scored   (518 eliminated genes)

    ## 
    ##   Level 5:   31 nodes to be scored   (1016 eliminated genes)

    ## 
    ##   Level 4:   28 nodes to be scored   (1216 eliminated genes)

    ## 
    ##   Level 3:   14 nodes to be scored   (1887 eliminated genes)

    ## 
    ##   Level 2:   6 nodes to be scored    (2152 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2421 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0002042      cell migration involved in sprouting angiogenesis         7
    ## 2  GO:0000103                                   sulfate assimilation         1
    ## 3  GO:0008277 regulation of G protein-coupled receptor signaling ...         2
    ## 4  GO:0003964                   RNA-directed DNA polymerase activity       122
    ## 5  GO:0003676                                   nucleic acid binding       497
    ## 6  GO:0004134                    4-alpha-glucanotransferase activity         1
    ## 7  GO:0001002 RNA polymerase III type 1 promoter sequence-specifi...         1
    ## 8  GO:0003980 UDP-glucose:glycoprotein glucosyltransferase activi...         1
    ## 9  GO:0004829                         threonine-tRNA ligase activity         1
    ## 10 GO:0000703 oxidized pyrimidine nucleobase lesion DNA N-glycosy...         1
    ##    Significant Expected Fisher               type
    ## 1            2     0.16 0.0098 Biological.Process
    ## 2            1     0.02 0.0228 Biological.Process
    ## 3            1     0.05 0.0451 Biological.Process
    ## 4            9     3.17 0.0037 Molecular.Function
    ## 5           14    12.90 0.0224 Molecular.Function
    ## 6            1     0.03 0.0260 Molecular.Function
    ## 7            1     0.03 0.0260 Molecular.Function
    ## 8            1     0.03 0.0260 Molecular.Function
    ## 9            1     0.03 0.0260 Molecular.Function
    ## 10           1     0.03 0.0260 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 125 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   8 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 8:   13 nodes to be scored   (161 eliminated genes)

    ## 
    ##   Level 7:   13 nodes to be scored   (171 eliminated genes)

    ## 
    ##   Level 6:   17 nodes to be scored   (229 eliminated genes)

    ## 
    ##   Level 5:   21 nodes to be scored   (509 eliminated genes)

    ## 
    ##   Level 4:   19 nodes to be scored   (560 eliminated genes)

    ## 
    ##   Level 3:   20 nodes to be scored   (776 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (1025 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1292 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 87 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   6 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 8:   7 nodes to be scored    (101 eliminated genes)

    ## 
    ##   Level 7:   12 nodes to be scored   (396 eliminated genes)

    ## 
    ##   Level 6:   17 nodes to be scored   (445 eliminated genes)

    ## 
    ##   Level 5:   16 nodes to be scored   (664 eliminated genes)

    ## 
    ##   Level 4:   12 nodes to be scored   (901 eliminated genes)

    ## 
    ##   Level 3:   9 nodes to be scored    (1442 eliminated genes)

    ## 
    ##   Level 2:   4 nodes to be scored    (1700 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2111 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0015969             guanosine tetraphosphate metabolic process         9
    ## 2 GO:0000122 negative regulation of transcription by RNA polymer...       140
    ## 3 GO:0001002 RNA polymerase III type 1 promoter sequence-specifi...         1
    ## 4 GO:0003747                    translation release factor activity         1
    ## 5 GO:0004315     3-oxoacyl-[acyl-carrier-protein] synthase activity         2
    ## 6 GO:0003676                                   nucleic acid binding       497
    ## 7 GO:0004622                             lysophospholipase activity         4
    ## 8 GO:0000340                      RNA 7-methylguanosine cap binding         4
    ##   Significant Expected Fisher               type
    ## 1           2     0.08 0.0027 Biological.Process
    ## 2           5     1.30 0.0061 Biological.Process
    ## 3           1     0.01 0.0110 Molecular.Function
    ## 4           1     0.01 0.0110 Molecular.Function
    ## 5           1     0.02 0.0220 Molecular.Function
    ## 6          14     5.56 0.0420 Molecular.Function
    ## 7           1     0.04 0.0440 Molecular.Function
    ## 8           1     0.04 0.0440 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 328 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  5 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  2 nodes to be scored    (11 eliminated genes)

    ## 
    ##   Level 11:  4 nodes to be scored    (31 eliminated genes)

    ## 
    ##   Level 10:  12 nodes to be scored   (74 eliminated genes)

    ## 
    ##   Level 9:   17 nodes to be scored   (230 eliminated genes)

    ## 
    ##   Level 8:   29 nodes to be scored   (317 eliminated genes)

    ## 
    ##   Level 7:   40 nodes to be scored   (345 eliminated genes)

    ## 
    ##   Level 6:   56 nodes to be scored   (526 eliminated genes)

    ## 
    ##   Level 5:   65 nodes to be scored   (821 eliminated genes)

    ## 
    ##   Level 4:   48 nodes to be scored   (887 eliminated genes)

    ## 
    ##   Level 3:   38 nodes to be scored   (1164 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (1308 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1396 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 295 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  6 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 9:   17 nodes to be scored   (13 eliminated genes)

    ## 
    ##   Level 8:   29 nodes to be scored   (151 eliminated genes)

    ## 
    ##   Level 7:   41 nodes to be scored   (479 eliminated genes)

    ## 
    ##   Level 6:   51 nodes to be scored   (647 eliminated genes)

    ## 
    ##   Level 5:   55 nodes to be scored   (1092 eliminated genes)

    ## 
    ##   Level 4:   51 nodes to be scored   (1401 eliminated genes)

    ## 
    ##   Level 3:   28 nodes to be scored   (2033 eliminated genes)

    ## 
    ##   Level 2:   12 nodes to be scored   (2365 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2631 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0009734                      auxin-activated signaling pathway         4
    ## 2  GO:0001731         formation of translation preinitiation complex         1
    ## 3  GO:0000390                       spliceosomal complex disassembly         1
    ## 4  GO:0001922                                 B-1 B cell homeostasis         1
    ## 5  GO:0002230 positive regulation of defense response to virus by...         1
    ## 6  GO:0000103                                   sulfate assimilation         1
    ## 7  GO:0006884                                cell volume homeostasis         1
    ## 8  GO:0005543                                   phospholipid binding         8
    ## 9  GO:0000987    cis-regulatory region sequence-specific DNA binding       112
    ## 10 GO:0003676                                   nucleic acid binding       497
    ## 11 GO:0005044                            scavenger receptor activity        13
    ##    Significant Expected Fisher               type
    ## 1            2     0.18  0.012 Biological.Process
    ## 2            1     0.05  0.046 Biological.Process
    ## 3            1     0.05  0.046 Biological.Process
    ## 4            1     0.05  0.046 Biological.Process
    ## 5            1     0.05  0.046 Biological.Process
    ## 6            1     0.05  0.046 Biological.Process
    ## 7            1     0.05  0.046 Biological.Process
    ## 8            4     0.52  0.001 Molecular.Function
    ## 9            8     7.23  0.030 Molecular.Function
    ## 10          36    32.08  0.037 Molecular.Function
    ## 11           3     0.84  0.047 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 386 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  3 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 12:  4 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 11:  9 nodes to be scored    (6 eliminated genes)

    ## 
    ##   Level 10:  15 nodes to be scored   (67 eliminated genes)

    ## 
    ##   Level 9:   22 nodes to be scored   (246 eliminated genes)

    ## 
    ##   Level 8:   42 nodes to be scored   (298 eliminated genes)

    ## 
    ##   Level 7:   53 nodes to be scored   (352 eliminated genes)

    ## 
    ##   Level 6:   62 nodes to be scored   (606 eliminated genes)

    ## 
    ##   Level 5:   76 nodes to be scored   (910 eliminated genes)

    ## 
    ##   Level 4:   50 nodes to be scored   (1018 eliminated genes)

    ## 
    ##   Level 3:   36 nodes to be scored   (1197 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (1267 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1363 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 300 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   11 nodes to be scored   (10 eliminated genes)

    ## 
    ##   Level 8:   20 nodes to be scored   (145 eliminated genes)

    ## 
    ##   Level 7:   49 nodes to be scored   (457 eliminated genes)

    ## 
    ##   Level 6:   63 nodes to be scored   (616 eliminated genes)

    ## 
    ##   Level 5:   60 nodes to be scored   (1088 eliminated genes)

    ## 
    ##   Level 4:   53 nodes to be scored   (1444 eliminated genes)

    ## 
    ##   Level 3:   25 nodes to be scored   (2044 eliminated genes)

    ## 
    ##   Level 2:   12 nodes to be scored   (2320 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2605 eliminated genes)

    ##        GO.ID                                          Term Annotated
    ## 1 GO:0000050                                    urea cycle         2
    ## 2 GO:0006487                protein N-linked glycosylation        19
    ## 3 GO:0000774    adenyl-nucleotide exchange factor activity         2
    ## 4 GO:0005302 L-tyrosine transmembrane transporter activity         9
    ## 5 GO:0004560                   alpha-L-fucosidase activity         6
    ##   Significant Expected Fisher               type
    ## 1           2     0.12 0.0036 Biological.Process
    ## 2           4     1.15 0.0241 Biological.Process
    ## 3           2     0.12 0.0037 Molecular.Function
    ## 4           3     0.55 0.0142 Molecular.Function
    ## 5           2     0.37 0.0471 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 703 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  8 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 12:  11 nodes to be scored   (16 eliminated genes)

    ## 
    ##   Level 11:  19 nodes to be scored   (41 eliminated genes)

    ## 
    ##   Level 10:  36 nodes to be scored   (102 eliminated genes)

    ## 
    ##   Level 9:   53 nodes to be scored   (296 eliminated genes)

    ## 
    ##   Level 8:   79 nodes to be scored   (427 eliminated genes)

    ## 
    ##   Level 7:   104 nodes to be scored  (537 eliminated genes)

    ## 
    ##   Level 6:   120 nodes to be scored  (740 eliminated genes)

    ## 
    ##   Level 5:   124 nodes to be scored  (1066 eliminated genes)

    ## 
    ##   Level 4:   76 nodes to be scored   (1143 eliminated genes)

    ## 
    ##   Level 3:   55 nodes to be scored   (1290 eliminated genes)

    ## 
    ##   Level 2:   13 nodes to be scored   (1339 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1404 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 477 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  8 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 9:   21 nodes to be scored   (15 eliminated genes)

    ## 
    ##   Level 8:   40 nodes to be scored   (161 eliminated genes)

    ## 
    ##   Level 7:   82 nodes to be scored   (491 eliminated genes)

    ## 
    ##   Level 6:   104 nodes to be scored  (693 eliminated genes)

    ## 
    ##   Level 5:   86 nodes to be scored   (1241 eliminated genes)

    ## 
    ##   Level 4:   78 nodes to be scored   (1623 eliminated genes)

    ## 
    ##   Level 3:   39 nodes to be scored   (2181 eliminated genes)

    ## 
    ##   Level 2:   14 nodes to be scored   (2479 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2678 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0000398                         mRNA splicing, via spliceosome        63
    ## 2  GO:0000165                                           MAPK cascade        35
    ## 3  GO:0002218                   activation of innate immune response        26
    ## 4  GO:0002674     negative regulation of acute inflammatory response         4
    ## 5  GO:0003341                                        cilium movement        20
    ## 6  GO:0001843                                    neural tube closure         7
    ## 7  GO:0005524                                            ATP binding       263
    ## 8  GO:0004364                       glutathione transferase activity         4
    ## 9  GO:0000030                           mannosyltransferase activity        12
    ## 10 GO:0004675 transmembrane receptor protein serine/threonine kin...         3
    ## 11 GO:0003723                                            RNA binding       128
    ## 12 GO:0005198                           structural molecule activity        54
    ## 13 GO:0000146                           microfilament motor activity         9
    ## 14 GO:0005085             guanyl-nucleotide exchange factor activity        12
    ## 15 GO:0000340                      RNA 7-methylguanosine cap binding         4
    ## 16 GO:0005229 intracellularly calcium-gated chloride channel acti...         4
    ##    Significant Expected  Fisher               type
    ## 1           25    14.13 0.00150 Biological.Process
    ## 2           16     7.85 0.00170 Biological.Process
    ## 3           11     5.83 0.03140 Biological.Process
    ## 4            3     0.90 0.03740 Biological.Process
    ## 5            9     4.49 0.04240 Biological.Process
    ## 6            4     1.57 0.04890 Biological.Process
    ## 7           87    62.79 0.00023 Molecular.Function
    ## 8            4     0.95 0.00323 Molecular.Function
    ## 9            9     2.86 0.01018 Molecular.Function
    ## 10           3     0.72 0.01356 Molecular.Function
    ## 11          39    30.56 0.02311 Molecular.Function
    ## 12          21    12.89 0.03889 Molecular.Function
    ## 13           5     2.15 0.04041 Molecular.Function
    ## 14           6     2.86 0.04391 Molecular.Function
    ## 15           3     0.95 0.04456 Molecular.Function
    ## 16           3     0.95 0.04456 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 357 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  6 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  11 nodes to be scored   (27 eliminated genes)

    ## 
    ##   Level 10:  18 nodes to be scored   (85 eliminated genes)

    ## 
    ##   Level 9:   23 nodes to be scored   (261 eliminated genes)

    ## 
    ##   Level 8:   31 nodes to be scored   (352 eliminated genes)

    ## 
    ##   Level 7:   49 nodes to be scored   (403 eliminated genes)

    ## 
    ##   Level 6:   62 nodes to be scored   (576 eliminated genes)

    ## 
    ##   Level 5:   66 nodes to be scored   (880 eliminated genes)

    ## 
    ##   Level 4:   43 nodes to be scored   (988 eliminated genes)

    ## 
    ##   Level 3:   33 nodes to be scored   (1139 eliminated genes)

    ## 
    ##   Level 2:   11 nodes to be scored   (1257 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1365 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 260 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 9:   7 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 8:   18 nodes to be scored   (152 eliminated genes)

    ## 
    ##   Level 7:   38 nodes to be scored   (452 eliminated genes)

    ## 
    ##   Level 6:   50 nodes to be scored   (608 eliminated genes)

    ## 
    ##   Level 5:   50 nodes to be scored   (1122 eliminated genes)

    ## 
    ##   Level 4:   49 nodes to be scored   (1458 eliminated genes)

    ## 
    ##   Level 3:   29 nodes to be scored   (2053 eliminated genes)

    ## 
    ##   Level 2:   12 nodes to be scored   (2360 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2640 eliminated genes)

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
    ## 1            7     2.02 0.0030 Biological.Process
    ## 2            5     1.15 0.0044 Biological.Process
    ## 3            3     0.69 0.0279 Biological.Process
    ## 4            2     0.29 0.0293 Biological.Process
    ## 5            4     0.95 0.0120 Molecular.Function
    ## 6            6     2.20 0.0190 Molecular.Function
    ## 7            4     1.17 0.0250 Molecular.Function
    ## 8            2     0.29 0.0290 Molecular.Function
    ## 9           15     8.93 0.0300 Molecular.Function
    ## 10           8     3.88 0.0360 Molecular.Function
    ## 11          27    19.25 0.0400 Molecular.Function
    ## 12           2     0.37 0.0460 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 440 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  6 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 11:  12 nodes to be scored   (9 eliminated genes)

    ## 
    ##   Level 10:  23 nodes to be scored   (78 eliminated genes)

    ## 
    ##   Level 9:   31 nodes to be scored   (253 eliminated genes)

    ## 
    ##   Level 8:   45 nodes to be scored   (368 eliminated genes)

    ## 
    ##   Level 7:   61 nodes to be scored   (437 eliminated genes)

    ## 
    ##   Level 6:   71 nodes to be scored   (618 eliminated genes)

    ## 
    ##   Level 5:   76 nodes to be scored   (897 eliminated genes)

    ## 
    ##   Level 4:   58 nodes to be scored   (1003 eliminated genes)

    ## 
    ##   Level 3:   42 nodes to be scored   (1206 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (1341 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1391 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 234 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   9 nodes to be scored    (9 eliminated genes)

    ## 
    ##   Level 8:   16 nodes to be scored   (143 eliminated genes)

    ## 
    ##   Level 7:   31 nodes to be scored   (448 eliminated genes)

    ## 
    ##   Level 6:   46 nodes to be scored   (579 eliminated genes)

    ## 
    ##   Level 5:   45 nodes to be scored   (1019 eliminated genes)

    ## 
    ##   Level 4:   43 nodes to be scored   (1307 eliminated genes)

    ## 
    ##   Level 3:   27 nodes to be scored   (1990 eliminated genes)

    ## 
    ##   Level 2:   12 nodes to be scored   (2313 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2649 eliminated genes)

    ##         GO.ID                                         Term Annotated
    ## 1  GO:0000165                                 MAPK cascade        35
    ## 2  GO:0008218                              bioluminescence        10
    ## 3  GO:0003341                              cilium movement        20
    ## 4  GO:0000038 very long-chain fatty acid metabolic process         9
    ## 5  GO:0006886              intracellular protein transport         4
    ## 6  GO:0005524                                  ATP binding       263
    ## 7  GO:0003964         RNA-directed DNA polymerase activity       122
    ## 8  GO:0004342   glucosamine-6-phosphate deaminase activity         2
    ## 9  GO:0005198                 structural molecule activity        54
    ## 10 GO:0000026       alpha-1,2-mannosyltransferase activity         3
    ## 11 GO:0001786                   phosphatidylserine binding         4
    ## 12 GO:0005201  extracellular matrix structural constituent        17
    ## 13 GO:0002020                             protease binding         5
    ##    Significant Expected   Fisher               type
    ## 1           10     2.69 0.000170 Biological.Process
    ## 2            5     0.77 0.000450 Biological.Process
    ## 3            5     1.54 0.014950 Biological.Process
    ## 4            3     0.69 0.026410 Biological.Process
    ## 5            2     0.31 0.031750 Biological.Process
    ## 6           36    18.97 0.000072 Molecular.Function
    ## 7           18     8.80 0.002300 Molecular.Function
    ## 8            2     0.14 0.005200 Molecular.Function
    ## 9           11     3.89 0.005400 Molecular.Function
    ## 10           2     0.22 0.014800 Molecular.Function
    ## 11           2     0.29 0.028200 Molecular.Function
    ## 12           4     1.23 0.029700 Molecular.Function
    ## 13           2     0.36 0.044800 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 259 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  3 nodes to be scored    (8 eliminated genes)

    ## 
    ##   Level 11:  7 nodes to be scored    (13 eliminated genes)

    ## 
    ##   Level 10:  12 nodes to be scored   (68 eliminated genes)

    ## 
    ##   Level 9:   18 nodes to be scored   (225 eliminated genes)

    ## 
    ##   Level 8:   30 nodes to be scored   (281 eliminated genes)

    ## 
    ##   Level 7:   33 nodes to be scored   (306 eliminated genes)

    ## 
    ##   Level 6:   39 nodes to be scored   (426 eliminated genes)

    ## 
    ##   Level 5:   45 nodes to be scored   (679 eliminated genes)

    ## 
    ##   Level 4:   32 nodes to be scored   (743 eliminated genes)

    ## 
    ##   Level 3:   25 nodes to be scored   (965 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (1218 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1338 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 140 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   7 nodes to be scored    (142 eliminated genes)

    ## 
    ##   Level 7:   18 nodes to be scored   (432 eliminated genes)

    ## 
    ##   Level 6:   29 nodes to be scored   (497 eliminated genes)

    ## 
    ##   Level 5:   30 nodes to be scored   (870 eliminated genes)

    ## 
    ##   Level 4:   27 nodes to be scored   (1194 eliminated genes)

    ## 
    ##   Level 3:   15 nodes to be scored   (1751 eliminated genes)

    ## 
    ##   Level 2:   7 nodes to be scored    (2022 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2517 eliminated genes)

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
    ## 11 GO:0004735             pyrroline-5-carboxylate reductase activity         2
    ## 12 GO:0004658                     propionyl-CoA carboxylase activity         2
    ##    Significant Expected Fisher               type
    ## 1            2     0.12 0.0048 Biological.Process
    ## 2            2     0.15 0.0079 Biological.Process
    ## 3            2     0.18 0.0116 Biological.Process
    ## 4            1     0.03 0.0292 Biological.Process
    ## 5            6     1.20 0.0011 Molecular.Function
    ## 6            2     0.20 0.0165 Molecular.Function
    ## 7            1     0.02 0.0227 Molecular.Function
    ## 8            1     0.02 0.0227 Molecular.Function
    ## 9            1     0.02 0.0227 Molecular.Function
    ## 10           6     2.27 0.0244 Molecular.Function
    ## 11           1     0.05 0.0449 Molecular.Function
    ## 12           1     0.05 0.0449 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 2 1 1 2 1 2 1 1 1 2 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 267 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  5 nodes to be scored    (11 eliminated genes)

    ## 
    ##   Level 10:  10 nodes to be scored   (84 eliminated genes)

    ## 
    ##   Level 9:   19 nodes to be scored   (255 eliminated genes)

    ## 
    ##   Level 8:   24 nodes to be scored   (307 eliminated genes)

    ## 
    ##   Level 7:   34 nodes to be scored   (373 eliminated genes)

    ## 
    ##   Level 6:   47 nodes to be scored   (443 eliminated genes)

    ## 
    ##   Level 5:   47 nodes to be scored   (693 eliminated genes)

    ## 
    ##   Level 4:   37 nodes to be scored   (796 eliminated genes)

    ## 
    ##   Level 3:   29 nodes to be scored   (1042 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (1189 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1378 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 188 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  3 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 9:   7 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 8:   9 nodes to be scored    (110 eliminated genes)

    ## 
    ##   Level 7:   22 nodes to be scored   (440 eliminated genes)

    ## 
    ##   Level 6:   32 nodes to be scored   (542 eliminated genes)

    ## 
    ##   Level 5:   40 nodes to be scored   (911 eliminated genes)

    ## 
    ##   Level 4:   40 nodes to be scored   (1290 eliminated genes)

    ## 
    ##   Level 3:   23 nodes to be scored   (1990 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (2304 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2607 eliminated genes)

    ##         GO.ID                                             Term Annotated
    ## 1  GO:0002933                              lipid hydroxylation         3
    ## 2  GO:0006368    transcription elongation by RNA polymerase II         1
    ## 3  GO:0015969       guanosine tetraphosphate metabolic process         9
    ## 4  GO:0004972                 NMDA glutamate receptor activity         3
    ## 5  GO:0005198                     structural molecule activity        54
    ## 6  GO:0005245           voltage-gated calcium channel activity         7
    ## 7  GO:0005509                              calcium ion binding        66
    ## 8  GO:0003909                              DNA ligase activity         1
    ## 9  GO:0004441 inositol-1,4-bisphosphate 1-phosphatase activity         1
    ## 10 GO:0000400                    four-way junction DNA binding         1
    ## 11 GO:0005542                               folic acid binding        10
    ##    Significant Expected Fisher               type
    ## 1            2     0.11 0.0036 Biological.Process
    ## 2            1     0.04 0.0356 Biological.Process
    ## 3            2     0.32 0.0381 Biological.Process
    ## 4            2     0.11 0.0036 Molecular.Function
    ## 5            6     1.91 0.0087 Molecular.Function
    ## 6            2     0.25 0.0231 Molecular.Function
    ## 7            6     2.33 0.0274 Molecular.Function
    ## 8            1     0.04 0.0353 Molecular.Function
    ## 9            1     0.04 0.0353 Molecular.Function
    ## 10           1     0.04 0.0353 Molecular.Function
    ## 11           2     0.35 0.0462 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 465 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  5 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  7 nodes to be scored    (15 eliminated genes)

    ## 
    ##   Level 11:  10 nodes to be scored   (35 eliminated genes)

    ## 
    ##   Level 10:  19 nodes to be scored   (100 eliminated genes)

    ## 
    ##   Level 9:   31 nodes to be scored   (269 eliminated genes)

    ## 
    ##   Level 8:   47 nodes to be scored   (362 eliminated genes)

    ## 
    ##   Level 7:   66 nodes to be scored   (455 eliminated genes)

    ## 
    ##   Level 6:   83 nodes to be scored   (592 eliminated genes)

    ## 
    ##   Level 5:   82 nodes to be scored   (867 eliminated genes)

    ## 
    ##   Level 4:   57 nodes to be scored   (948 eliminated genes)

    ## 
    ##   Level 3:   44 nodes to be scored   (1250 eliminated genes)

    ## 
    ##   Level 2:   11 nodes to be scored   (1321 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1395 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 293 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   12 nodes to be scored   (9 eliminated genes)

    ## 
    ##   Level 8:   22 nodes to be scored   (147 eliminated genes)

    ## 
    ##   Level 7:   55 nodes to be scored   (470 eliminated genes)

    ## 
    ##   Level 6:   62 nodes to be scored   (634 eliminated genes)

    ## 
    ##   Level 5:   57 nodes to be scored   (1113 eliminated genes)

    ## 
    ##   Level 4:   45 nodes to be scored   (1433 eliminated genes)

    ## 
    ##   Level 3:   25 nodes to be scored   (2032 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (2322 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2585 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0002674     negative regulation of acute inflammatory response         4
    ## 2 GO:0001818             negative regulation of cytokine production        18
    ## 3 GO:0001732 formation of cytoplasmic translation initiation com...         6
    ## 4 GO:0005524                                            ATP binding       263
    ## 5 GO:0003964                   RNA-directed DNA polymerase activity       122
    ## 6 GO:0004177                                aminopeptidase activity         7
    ## 7 GO:0003682                                      chromatin binding        13
    ## 8 GO:0000340                      RNA 7-methylguanosine cap binding         4
    ## 9 GO:0004623                              phospholipase A2 activity         4
    ##   Significant Expected Fisher               type
    ## 1           3     0.42 0.0041 Biological.Process
    ## 2           6     1.87 0.0074 Biological.Process
    ## 3           3     0.62 0.0174 Biological.Process
    ## 4          39    24.75 0.0019 Molecular.Function
    ## 5          21    11.48 0.0040 Molecular.Function
    ## 6           3     0.66 0.0217 Molecular.Function
    ## 7           4     1.22 0.0277 Molecular.Function
    ## 8           2     0.38 0.0466 Molecular.Function
    ## 9           2     0.38 0.0466 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 543 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  4 nodes to be scored    (2 eliminated genes)

    ## 
    ##   Level 12:  13 nodes to be scored   (3 eliminated genes)

    ## 
    ##   Level 11:  16 nodes to be scored   (6 eliminated genes)

    ## 
    ##   Level 10:  28 nodes to be scored   (95 eliminated genes)

    ## 
    ##   Level 9:   38 nodes to be scored   (245 eliminated genes)

    ## 
    ##   Level 8:   56 nodes to be scored   (363 eliminated genes)

    ## 
    ##   Level 7:   77 nodes to be scored   (445 eliminated genes)

    ## 
    ##   Level 6:   92 nodes to be scored   (657 eliminated genes)

    ## 
    ##   Level 5:   93 nodes to be scored   (989 eliminated genes)

    ## 
    ##   Level 4:   63 nodes to be scored   (1104 eliminated genes)

    ## 
    ##   Level 3:   48 nodes to be scored   (1246 eliminated genes)

    ## 
    ##   Level 2:   12 nodes to be scored   (1323 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1400 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 258 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 9:   9 nodes to be scored    (13 eliminated genes)

    ## 
    ##   Level 8:   21 nodes to be scored   (146 eliminated genes)

    ## 
    ##   Level 7:   34 nodes to be scored   (447 eliminated genes)

    ## 
    ##   Level 6:   52 nodes to be scored   (567 eliminated genes)

    ## 
    ##   Level 5:   51 nodes to be scored   (1046 eliminated genes)

    ## 
    ##   Level 4:   44 nodes to be scored   (1364 eliminated genes)

    ## 
    ##   Level 3:   26 nodes to be scored   (1960 eliminated genes)

    ## 
    ##   Level 2:   11 nodes to be scored   (2276 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2619 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0000184 nuclear-transcribed mRNA catabolic process, nonsens...        10
    ## 2  GO:0002064                            epithelial cell development         4
    ## 3  GO:0001782                                     B cell homeostasis         7
    ## 4  GO:0003341                                        cilium movement        20
    ## 5  GO:0000723                                   telomere maintenance         4
    ## 6  GO:0002244          hematopoietic progenitor cell differentiation         4
    ## 7  GO:0003676                                   nucleic acid binding       497
    ## 8  GO:0001409 guanine nucleotide transmembrane transporter activi...        30
    ## 9  GO:0003964                   RNA-directed DNA polymerase activity       122
    ## 10 GO:0005509                                    calcium ion binding        66
    ## 11 GO:0003824                                     catalytic activity      1038
    ## 12 GO:0005085             guanyl-nucleotide exchange factor activity        12
    ## 13 GO:0004972                       NMDA glutamate receptor activity         3
    ## 14 GO:0005507                                     copper ion binding        18
    ## 15 GO:0003723                                            RNA binding       128
    ##    Significant Expected  Fisher               type
    ## 1            5     0.96 0.00130 Biological.Process
    ## 2            3     0.38 0.00320 Biological.Process
    ## 3            3     0.67 0.02280 Biological.Process
    ## 4            6     1.92 0.02860 Biological.Process
    ## 5            2     0.38 0.04840 Biological.Process
    ## 6            2     0.38 0.04840 Biological.Process
    ## 7           64    50.90 0.00058 Molecular.Function
    ## 8            8     3.07 0.00859 Molecular.Function
    ## 9           21    12.49 0.01064 Molecular.Function
    ## 10          13     6.76 0.01430 Molecular.Function
    ## 11          86   106.31 0.02042 Molecular.Function
    ## 12           4     1.23 0.02745 Molecular.Function
    ## 13           2     0.31 0.02924 Molecular.Function
    ## 14           5     1.84 0.03043 Molecular.Function
    ## 15          20    13.11 0.03136 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 650 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  5 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  6 nodes to be scored    (2 eliminated genes)

    ## 
    ##   Level 12:  8 nodes to be scored    (23 eliminated genes)

    ## 
    ##   Level 11:  16 nodes to be scored   (32 eliminated genes)

    ## 
    ##   Level 10:  34 nodes to be scored   (95 eliminated genes)

    ## 
    ##   Level 9:   54 nodes to be scored   (271 eliminated genes)

    ## 
    ##   Level 8:   86 nodes to be scored   (400 eliminated genes)

    ## 
    ##   Level 7:   99 nodes to be scored   (494 eliminated genes)

    ## 
    ##   Level 6:   106 nodes to be scored  (681 eliminated genes)

    ## 
    ##   Level 5:   100 nodes to be scored  (987 eliminated genes)

    ## 
    ##   Level 4:   70 nodes to be scored   (1100 eliminated genes)

    ## 
    ##   Level 3:   51 nodes to be scored   (1278 eliminated genes)

    ## 
    ##   Level 2:   13 nodes to be scored   (1352 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1403 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 455 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  10 nodes to be scored   (4 eliminated genes)

    ## 
    ##   Level 9:   23 nodes to be scored   (15 eliminated genes)

    ## 
    ##   Level 8:   43 nodes to be scored   (163 eliminated genes)

    ## 
    ##   Level 7:   73 nodes to be scored   (489 eliminated genes)

    ## 
    ##   Level 6:   94 nodes to be scored   (694 eliminated genes)

    ## 
    ##   Level 5:   81 nodes to be scored   (1212 eliminated genes)

    ## 
    ##   Level 4:   73 nodes to be scored   (1544 eliminated genes)

    ## 
    ##   Level 3:   37 nodes to be scored   (2154 eliminated genes)

    ## 
    ##   Level 2:   14 nodes to be scored   (2492 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2678 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0003341                                        cilium movement        20
    ## 2  GO:0001933         negative regulation of protein phosphorylation         5
    ## 3  GO:0006886                        intracellular protein transport         4
    ## 4  GO:0002244          hematopoietic progenitor cell differentiation         4
    ## 5  GO:0001934         positive regulation of protein phosphorylation         4
    ## 6  GO:0002042      cell migration involved in sprouting angiogenesis         7
    ## 7  GO:0001895                                     retina homeostasis         2
    ## 8  GO:0007224                           smoothened signaling pathway         2
    ## 9  GO:0005524                                            ATP binding       263
    ## 10 GO:0004715 non-membrane spanning protein tyrosine kinase activ...         6
    ## 11 GO:0005388                    P-type calcium transporter activity         4
    ## 12 GO:0003964                   RNA-directed DNA polymerase activity       122
    ## 13 GO:0003676                                   nucleic acid binding       497
    ## 14 GO:0004518                                      nuclease activity        30
    ## 15 GO:0001965                        G-protein alpha-subunit binding        16
    ## 16 GO:0001786                             phosphatidylserine binding         4
    ## 17 GO:0005198                           structural molecule activity        54
    ## 18 GO:0004100                               chitin synthase activity         2
    ## 19 GO:0004579 dolichyl-diphosphooligosaccharide-protein glycotran...         2
    ## 20 GO:0003987                            acetate-CoA ligase activity         2
    ##    Significant Expected   Fisher               type
    ## 1           12     3.62 0.000031 Biological.Process
    ## 2            4     0.90 0.004500 Biological.Process
    ## 3            4     0.72 0.005800 Biological.Process
    ## 4            3     0.72 0.020300 Biological.Process
    ## 5            3     0.72 0.020300 Biological.Process
    ## 6            4     1.27 0.023200 Biological.Process
    ## 7            2     0.36 0.032600 Biological.Process
    ## 8            2     0.36 0.032600 Biological.Process
    ## 9           72    51.12 0.000630 Molecular.Function
    ## 10           5     1.17 0.001380 Molecular.Function
    ## 11           4     0.78 0.001410 Molecular.Function
    ## 12          37    23.71 0.002150 Molecular.Function
    ## 13         104    96.60 0.005260 Molecular.Function
    ## 14           9     5.83 0.015370 Molecular.Function
    ## 15           7     3.11 0.022640 Molecular.Function
    ## 16           3     0.78 0.025000 Molecular.Function
    ## 17          14    10.50 0.030240 Molecular.Function
    ## 18           2     0.39 0.037720 Molecular.Function
    ## 19           2     0.39 0.037720 Molecular.Function
    ## 20           2     0.39 0.037720 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 333 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  6 nodes to be scored    (15 eliminated genes)

    ## 
    ##   Level 11:  10 nodes to be scored   (24 eliminated genes)

    ## 
    ##   Level 10:  21 nodes to be scored   (78 eliminated genes)

    ## 
    ##   Level 9:   28 nodes to be scored   (253 eliminated genes)

    ## 
    ##   Level 8:   33 nodes to be scored   (304 eliminated genes)

    ## 
    ##   Level 7:   45 nodes to be scored   (378 eliminated genes)

    ## 
    ##   Level 6:   54 nodes to be scored   (456 eliminated genes)

    ## 
    ##   Level 5:   55 nodes to be scored   (736 eliminated genes)

    ## 
    ##   Level 4:   37 nodes to be scored   (845 eliminated genes)

    ## 
    ##   Level 3:   28 nodes to be scored   (1023 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (1165 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1306 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 220 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  3 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 9:   9 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 8:   15 nodes to be scored   (143 eliminated genes)

    ## 
    ##   Level 7:   35 nodes to be scored   (452 eliminated genes)

    ## 
    ##   Level 6:   42 nodes to be scored   (607 eliminated genes)

    ## 
    ##   Level 5:   47 nodes to be scored   (1082 eliminated genes)

    ## 
    ##   Level 4:   39 nodes to be scored   (1308 eliminated genes)

    ## 
    ##   Level 3:   18 nodes to be scored   (1944 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (2212 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2559 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0002064                            epithelial cell development         4
    ## 2  GO:0000278                                     mitotic cell cycle       160
    ## 3  GO:0003179                              heart valve morphogenesis         4
    ## 4  GO:0001701                         in utero embryonic development        15
    ## 5  GO:0001732 formation of cytoplasmic translation initiation com...         6
    ## 6  GO:0000712      resolution of meiotic recombination intermediates         1
    ## 7  GO:0000472 endonucleolytic cleavage to generate mature 5'-end ...         1
    ## 8  GO:0002225 positive regulation of antimicrobial peptide produc...         1
    ## 9  GO:0001573                          ganglioside metabolic process         1
    ## 10 GO:0001771                        immunological synapse formation         1
    ## 11 GO:0006513                             protein monoubiquitination         1
    ## 12 GO:0002314                 germinal center B cell differentiation         1
    ## 13 GO:0000987    cis-regulatory region sequence-specific DNA binding       112
    ## 14 GO:0050218                         propionate-CoA ligase activity         1
    ## 15 GO:0004368 glycerol-3-phosphate dehydrogenase (quinone) activi...         1
    ## 16 GO:0004818                         glutamate-tRNA ligase activity         1
    ## 17 GO:0003960                       NADPH:quinone reductase activity         1
    ## 18 GO:0004408                     holocytochrome-c synthase activity         1
    ## 19 GO:0003989                        acetyl-CoA carboxylase activity         1
    ## 20 GO:0004325                                ferrochelatase activity         1
    ## 21 GO:0004802                                 transketolase activity         1
    ## 22 GO:0004458          D-lactate dehydrogenase (cytochrome) activity         1
    ## 23 GO:0005220 inositol 1,4,5-trisphosphate-gated calcium channel ...         1
    ## 24 GO:0010181                                            FMN binding         1
    ## 25 GO:0004332                fructose-bisphosphate aldolase activity         1
    ## 26 GO:0004729 oxygen-dependent protoporphyrinogen oxidase activit...         1
    ## 27 GO:0008957       phenylacetaldehyde dehydrogenase (NAD+) activity         1
    ##    Significant Expected  Fisher               type
    ## 1            3     0.17 0.00029 Biological.Process
    ## 2            9     6.84 0.00311 Biological.Process
    ## 3            2     0.17 0.01020 Biological.Process
    ## 4            3     0.64 0.02340 Biological.Process
    ## 5            2     0.26 0.02412 Biological.Process
    ## 6            1     0.04 0.04274 Biological.Process
    ## 7            1     0.04 0.04274 Biological.Process
    ## 8            1     0.04 0.04274 Biological.Process
    ## 9            1     0.04 0.04274 Biological.Process
    ## 10           1     0.04 0.04274 Biological.Process
    ## 11           1     0.04 0.04274 Biological.Process
    ## 12           1     0.04 0.04274 Biological.Process
    ## 13           7     4.20 0.00075 Molecular.Function
    ## 14           1     0.04 0.03750 Molecular.Function
    ## 15           1     0.04 0.03750 Molecular.Function
    ## 16           1     0.04 0.03750 Molecular.Function
    ## 17           1     0.04 0.03750 Molecular.Function
    ## 18           1     0.04 0.03750 Molecular.Function
    ## 19           1     0.04 0.03750 Molecular.Function
    ## 20           1     0.04 0.03750 Molecular.Function
    ## 21           1     0.04 0.03750 Molecular.Function
    ## 22           1     0.04 0.03750 Molecular.Function
    ## 23           1     0.04 0.03750 Molecular.Function
    ## 24           1     0.04 0.03750 Molecular.Function
    ## 25           1     0.04 0.03750 Molecular.Function
    ## 26           1     0.04 0.03750 Molecular.Function
    ## 27           1     0.04 0.03750 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 345 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  6 nodes to be scored    (15 eliminated genes)

    ## 
    ##   Level 11:  12 nodes to be scored   (24 eliminated genes)

    ## 
    ##   Level 10:  17 nodes to be scored   (81 eliminated genes)

    ## 
    ##   Level 9:   21 nodes to be scored   (263 eliminated genes)

    ## 
    ##   Level 8:   36 nodes to be scored   (309 eliminated genes)

    ## 
    ##   Level 7:   53 nodes to be scored   (381 eliminated genes)

    ## 
    ##   Level 6:   57 nodes to be scored   (509 eliminated genes)

    ## 
    ##   Level 5:   57 nodes to be scored   (786 eliminated genes)

    ## 
    ##   Level 4:   39 nodes to be scored   (850 eliminated genes)

    ## 
    ##   Level 3:   30 nodes to be scored   (1085 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (1248 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1354 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 200 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  5 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   8 nodes to be scored    (2 eliminated genes)

    ## 
    ##   Level 8:   18 nodes to be scored   (155 eliminated genes)

    ## 
    ##   Level 7:   33 nodes to be scored   (446 eliminated genes)

    ## 
    ##   Level 6:   36 nodes to be scored   (601 eliminated genes)

    ## 
    ##   Level 5:   41 nodes to be scored   (1048 eliminated genes)

    ## 
    ##   Level 4:   29 nodes to be scored   (1275 eliminated genes)

    ## 
    ##   Level 3:   18 nodes to be scored   (1958 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (2163 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2428 eliminated genes)

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
    ## 10 GO:0004321                       fatty-acyl-CoA synthase activity         1
    ## 11 GO:0002134                                            UTP binding         1
    ## 12 GO:0004305                           ethanolamine kinase activity         1
    ## 13 GO:0005536                                      D-glucose binding         1
    ## 14 GO:0004492        methyl/ethyl malonyl-CoA decarboxylase activity         1
    ## 15 GO:0005324 long-chain fatty acid transmembrane transporter act...         1
    ##    Significant Expected   Fisher               type
    ## 1            3     0.28 0.002000 Biological.Process
    ## 2            2     0.12 0.004700 Biological.Process
    ## 3            1     0.04 0.040600 Biological.Process
    ## 4            1     0.04 0.040600 Biological.Process
    ## 5            1     0.04 0.040600 Biological.Process
    ## 6           17     5.76 0.000039 Molecular.Function
    ## 7            3     0.61 0.021000 Molecular.Function
    ## 8            7     3.12 0.034000 Molecular.Function
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

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 454 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  6 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 11:  13 nodes to be scored   (9 eliminated genes)

    ## 
    ##   Level 10:  24 nodes to be scored   (81 eliminated genes)

    ## 
    ##   Level 9:   34 nodes to be scored   (260 eliminated genes)

    ## 
    ##   Level 8:   44 nodes to be scored   (380 eliminated genes)

    ## 
    ##   Level 7:   69 nodes to be scored   (457 eliminated genes)

    ## 
    ##   Level 6:   73 nodes to be scored   (582 eliminated genes)

    ## 
    ##   Level 5:   76 nodes to be scored   (903 eliminated genes)

    ## 
    ##   Level 4:   57 nodes to be scored   (975 eliminated genes)

    ## 
    ##   Level 3:   43 nodes to be scored   (1127 eliminated genes)

    ## 
    ##   Level 2:   11 nodes to be scored   (1272 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1399 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 283 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   10 nodes to be scored   (0 eliminated genes)

    ## 
    ##   Level 8:   26 nodes to be scored   (102 eliminated genes)

    ## 
    ##   Level 7:   40 nodes to be scored   (406 eliminated genes)

    ## 
    ##   Level 6:   59 nodes to be scored   (615 eliminated genes)

    ## 
    ##   Level 5:   56 nodes to be scored   (1112 eliminated genes)

    ## 
    ##   Level 4:   48 nodes to be scored   (1466 eliminated genes)

    ## 
    ##   Level 3:   27 nodes to be scored   (2049 eliminated genes)

    ## 
    ##   Level 2:   13 nodes to be scored   (2325 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2623 eliminated genes)

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
    ## 3           2     0.22 0.0152 Biological.Process
    ## 4           4     1.10 0.0199 Biological.Process
    ## 5          17    10.27 0.0216 Biological.Process
    ## 6           3     0.66 0.0233 Biological.Process
    ## 7           4     1.00 0.0140 Molecular.Function
    ## 8           2     0.23 0.0170 Molecular.Function
    ## 9          15     9.42 0.0450 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 436 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  3 nodes to be scored    (2 eliminated genes)

    ## 
    ##   Level 12:  8 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 11:  11 nodes to be scored   (18 eliminated genes)

    ## 
    ##   Level 10:  18 nodes to be scored   (102 eliminated genes)

    ## 
    ##   Level 9:   35 nodes to be scored   (269 eliminated genes)

    ## 
    ##   Level 8:   48 nodes to be scored   (385 eliminated genes)

    ## 
    ##   Level 7:   65 nodes to be scored   (459 eliminated genes)

    ## 
    ##   Level 6:   80 nodes to be scored   (634 eliminated genes)

    ## 
    ##   Level 5:   74 nodes to be scored   (926 eliminated genes)

    ## 
    ##   Level 4:   44 nodes to be scored   (992 eliminated genes)

    ## 
    ##   Level 3:   36 nodes to be scored   (1159 eliminated genes)

    ## 
    ##   Level 2:   11 nodes to be scored   (1264 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1384 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 278 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  5 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 9:   12 nodes to be scored   (12 eliminated genes)

    ## 
    ##   Level 8:   21 nodes to be scored   (153 eliminated genes)

    ## 
    ##   Level 7:   43 nodes to be scored   (459 eliminated genes)

    ## 
    ##   Level 6:   56 nodes to be scored   (597 eliminated genes)

    ## 
    ##   Level 5:   51 nodes to be scored   (1083 eliminated genes)

    ## 
    ##   Level 4:   46 nodes to be scored   (1424 eliminated genes)

    ## 
    ##   Level 3:   27 nodes to be scored   (2064 eliminated genes)

    ## 
    ##   Level 2:   13 nodes to be scored   (2292 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2563 eliminated genes)

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
    ## 1          18     5.42 1.90e-06 Biological.Process
    ## 2           8     2.92 6.30e-03 Biological.Process
    ## 3           3     0.67 2.32e-02 Biological.Process
    ## 4           3     0.83 4.39e-02 Biological.Process
    ## 5          38    20.39 6.80e-05 Molecular.Function
    ## 6           8     2.33 1.50e-03 Molecular.Function
    ## 7          17     9.46 1.15e-02 Molecular.Function
    ## 8           2     0.23 1.70e-02 Molecular.Function
    ## 9           4     1.24 3.05e-02 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 283 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  4 nodes to be scored    (11 eliminated genes)

    ## 
    ##   Level 11:  7 nodes to be scored    (30 eliminated genes)

    ## 
    ##   Level 10:  13 nodes to be scored   (91 eliminated genes)

    ## 
    ##   Level 9:   18 nodes to be scored   (261 eliminated genes)

    ## 
    ##   Level 8:   26 nodes to be scored   (324 eliminated genes)

    ## 
    ##   Level 7:   35 nodes to be scored   (362 eliminated genes)

    ## 
    ##   Level 6:   44 nodes to be scored   (448 eliminated genes)

    ## 
    ##   Level 5:   53 nodes to be scored   (686 eliminated genes)

    ## 
    ##   Level 4:   42 nodes to be scored   (777 eliminated genes)

    ## 
    ##   Level 3:   27 nodes to be scored   (1026 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (1183 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1318 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 204 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 9:   6 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 8:   14 nodes to be scored   (152 eliminated genes)

    ## 
    ##   Level 7:   25 nodes to be scored   (439 eliminated genes)

    ## 
    ##   Level 6:   37 nodes to be scored   (602 eliminated genes)

    ## 
    ##   Level 5:   44 nodes to be scored   (990 eliminated genes)

    ## 
    ##   Level 4:   40 nodes to be scored   (1336 eliminated genes)

    ## 
    ##   Level 3:   21 nodes to be scored   (1994 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (2288 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2585 eliminated genes)

    ##         GO.ID                                               Term Annotated
    ## 1  GO:0002221     pattern recognition receptor signaling pathway        11
    ## 2  GO:0002674 negative regulation of acute inflammatory response         4
    ## 3  GO:0002040                             sprouting angiogenesis        12
    ## 4  GO:0001172                        RNA-templated transcription         1
    ## 5  GO:0000492                            box C/D snoRNP assembly         1
    ## 6  GO:0000076               DNA replication checkpoint signaling         1
    ## 7  GO:0005524                                        ATP binding       263
    ## 8  GO:0005201        extracellular matrix structural constituent        17
    ## 9  GO:0004017                          adenylate kinase activity         1
    ## 10 GO:0004505             phenylalanine 4-monooxygenase activity         1
    ##    Significant Expected   Fisher               type
    ## 1            5     0.47 0.000046 Biological.Process
    ## 2            3     0.17 0.000290 Biological.Process
    ## 3            3     0.51 0.016160 Biological.Process
    ## 4            1     0.04 0.042740 Biological.Process
    ## 5            1     0.04 0.042740 Biological.Process
    ## 6            1     0.04 0.042740 Biological.Process
    ## 7           20    11.29 0.007100 Molecular.Function
    ## 8            3     0.73 0.033700 Molecular.Function
    ## 9            1     0.04 0.042900 Molecular.Function
    ## 10           1     0.04 0.042900 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 366 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  5 nodes to be scored    (11 eliminated genes)

    ## 
    ##   Level 10:  11 nodes to be scored   (74 eliminated genes)

    ## 
    ##   Level 9:   24 nodes to be scored   (232 eliminated genes)

    ## 
    ##   Level 8:   39 nodes to be scored   (316 eliminated genes)

    ## 
    ##   Level 7:   51 nodes to be scored   (402 eliminated genes)

    ## 
    ##   Level 6:   65 nodes to be scored   (593 eliminated genes)

    ## 
    ##   Level 5:   69 nodes to be scored   (880 eliminated genes)

    ## 
    ##   Level 4:   48 nodes to be scored   (968 eliminated genes)

    ## 
    ##   Level 3:   38 nodes to be scored   (1189 eliminated genes)

    ## 
    ##   Level 2:   12 nodes to be scored   (1292 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1366 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 280 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   9 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 8:   23 nodes to be scored   (150 eliminated genes)

    ## 
    ##   Level 7:   41 nodes to be scored   (451 eliminated genes)

    ## 
    ##   Level 6:   52 nodes to be scored   (631 eliminated genes)

    ## 
    ##   Level 5:   56 nodes to be scored   (1089 eliminated genes)

    ## 
    ##   Level 4:   55 nodes to be scored   (1392 eliminated genes)

    ## 
    ##   Level 3:   26 nodes to be scored   (2010 eliminated genes)

    ## 
    ##   Level 2:   12 nodes to be scored   (2380 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2630 eliminated genes)

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
    ## 1            8     1.45 0.000036 Biological.Process
    ## 2            2     0.15 0.005200 Biological.Process
    ## 3           11     4.72 0.005600 Biological.Process
    ## 4            4     0.80 0.005800 Biological.Process
    ## 5            2     0.22 0.014900 Biological.Process
    ## 6            3     0.65 0.022700 Biological.Process
    ## 7            2     0.29 0.028500 Biological.Process
    ## 8           35    19.63 0.000340 Molecular.Function
    ## 9           17     9.11 0.007950 Molecular.Function
    ## 10           2     0.30 0.030090 Molecular.Function
    ## 11           2     0.30 0.030090 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 183 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  2 nodes to be scored    (8 eliminated genes)

    ## 
    ##   Level 11:  4 nodes to be scored    (27 eliminated genes)

    ## 
    ##   Level 10:  9 nodes to be scored    (67 eliminated genes)

    ## 
    ##   Level 9:   13 nodes to be scored   (207 eliminated genes)

    ## 
    ##   Level 8:   18 nodes to be scored   (287 eliminated genes)

    ## 
    ##   Level 7:   18 nodes to be scored   (303 eliminated genes)

    ## 
    ##   Level 6:   23 nodes to be scored   (380 eliminated genes)

    ## 
    ##   Level 5:   33 nodes to be scored   (524 eliminated genes)

    ## 
    ##   Level 4:   29 nodes to be scored   (698 eliminated genes)

    ## 
    ##   Level 3:   20 nodes to be scored   (998 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (1093 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1284 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 173 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   8 nodes to be scored    (10 eliminated genes)

    ## 
    ##   Level 8:   13 nodes to be scored   (143 eliminated genes)

    ## 
    ##   Level 7:   24 nodes to be scored   (449 eliminated genes)

    ## 
    ##   Level 6:   30 nodes to be scored   (566 eliminated genes)

    ## 
    ##   Level 5:   33 nodes to be scored   (985 eliminated genes)

    ## 
    ##   Level 4:   28 nodes to be scored   (1127 eliminated genes)

    ## 
    ##   Level 3:   21 nodes to be scored   (1814 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (2150 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2537 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0003341                                        cilium movement        20
    ## 2  GO:0001676                long-chain fatty acid metabolic process         3
    ## 3  GO:0001736                       establishment of planar polarity         5
    ## 4  GO:0003943             N-acetylgalactosamine-4-sulfatase activity         8
    ## 5  GO:0000048                           peptidyltransferase activity         1
    ## 6  GO:0005337          nucleoside transmembrane transporter activity         1
    ## 7  GO:0003980 UDP-glucose:glycoprotein glucosyltransferase activi...         1
    ## 8  GO:0004587                    ornithine aminotransferase activity         1
    ## 9  GO:0004149 dihydrolipoyllysine-residue succinyltransferase act...         1
    ## 10 GO:0042802                              identical protein binding         1
    ## 11 GO:0001002 RNA polymerase III type 1 promoter sequence-specifi...         1
    ## 12 GO:0004591 oxoglutarate dehydrogenase (succinyl-transferring) ...         1
    ## 13 GO:0004054                               arginine kinase activity         1
    ## 14 GO:0004760                L-serine-pyruvate transaminase activity         1
    ## 15 GO:0004842                 ubiquitin-protein transferase activity        35
    ## 16 GO:0004190                   aspartic-type endopeptidase activity         2
    ## 17 GO:0004714 transmembrane receptor protein tyrosine kinase acti...         2
    ## 18 GO:0005272                                sodium channel activity         2
    ##    Significant Expected  Fisher               type
    ## 1            4     0.40 0.00049 Biological.Process
    ## 2            2     0.06 0.00114 Biological.Process
    ## 3            2     0.10 0.00370 Biological.Process
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
    ## 15           3     0.76 0.03900 Molecular.Function
    ## 16           1     0.04 0.04300 Molecular.Function
    ## 17           1     0.04 0.04300 Molecular.Function
    ## 18           1     0.04 0.04300 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 207 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  2 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 11:  5 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 10:  11 nodes to be scored   (6 eliminated genes)

    ## 
    ##   Level 9:   15 nodes to be scored   (146 eliminated genes)

    ## 
    ##   Level 8:   25 nodes to be scored   (152 eliminated genes)

    ## 
    ##   Level 7:   27 nodes to be scored   (185 eliminated genes)

    ## 
    ##   Level 6:   27 nodes to be scored   (331 eliminated genes)

    ## 
    ##   Level 5:   39 nodes to be scored   (495 eliminated genes)

    ## 
    ##   Level 4:   27 nodes to be scored   (668 eliminated genes)

    ## 
    ##   Level 3:   17 nodes to be scored   (917 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (1186 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1297 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 162 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   10 nodes to be scored   (7 eliminated genes)

    ## 
    ##   Level 7:   23 nodes to be scored   (319 eliminated genes)

    ## 
    ##   Level 6:   32 nodes to be scored   (442 eliminated genes)

    ## 
    ##   Level 5:   33 nodes to be scored   (853 eliminated genes)

    ## 
    ##   Level 4:   33 nodes to be scored   (1060 eliminated genes)

    ## 
    ##   Level 3:   17 nodes to be scored   (1809 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (2153 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2442 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0008218                                        bioluminescence        10
    ## 2  GO:0002025 norepinephrine-epinephrine-mediated vasodilation in...         1
    ## 3  GO:0001731         formation of translation preinitiation complex         1
    ## 4  GO:0001832                                      blastocyst growth         1
    ## 5  GO:0001516                     prostaglandin biosynthetic process         2
    ## 6  GO:0002121                         inter-male aggressive behavior         2
    ## 7  GO:0004089                         carbonate dehydratase activity         2
    ## 8  GO:0005302          L-tyrosine transmembrane transporter activity         9
    ## 9  GO:0004067                                  asparaginase activity         1
    ## 10 GO:0004142      diacylglycerol cholinephosphotransferase activity         1
    ## 11 GO:0002161                        aminoacyl-tRNA editing activity         1
    ## 12 GO:0004017                              adenylate kinase activity         1
    ## 13 GO:0003913                                DNA photolyase activity         1
    ## 14 GO:0004252                     serine-type endopeptidase activity        53
    ## 15 GO:0004176                       ATP-dependent peptidase activity         2
    ## 16 GO:0004566                            beta-glucuronidase activity         2
    ## 17 GO:0004081 bis(5'-nucleosyl)-tetraphosphatase (asymmetrical) a...         2
    ## 18 GO:0003978                       UDP-glucose 4-epimerase activity         2
    ## 19 GO:0004720                      protein-lysine 6-oxidase activity         2
    ##    Significant Expected    Fisher               type
    ## 1            6     0.19 4.400e-09 Biological.Process
    ## 2            1     0.02 1.900e-02 Biological.Process
    ## 3            1     0.02 1.900e-02 Biological.Process
    ## 4            1     0.02 1.900e-02 Biological.Process
    ## 5            1     0.04 3.700e-02 Biological.Process
    ## 6            1     0.04 3.700e-02 Biological.Process
    ## 7            2     0.04 4.300e-04 Molecular.Function
    ## 8            2     0.19 1.409e-02 Molecular.Function
    ## 9            1     0.02 2.092e-02 Molecular.Function
    ## 10           1     0.02 2.092e-02 Molecular.Function
    ## 11           1     0.02 2.092e-02 Molecular.Function
    ## 12           1     0.02 2.092e-02 Molecular.Function
    ## 13           1     0.02 2.092e-02 Molecular.Function
    ## 14           4     1.11 2.357e-02 Molecular.Function
    ## 15           1     0.04 4.140e-02 Molecular.Function
    ## 16           1     0.04 4.140e-02 Molecular.Function
    ## 17           1     0.04 4.140e-02 Molecular.Function
    ## 18           1     0.04 4.140e-02 Molecular.Function
    ## 19           1     0.04 4.140e-02 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 525 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  5 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  6 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  12 nodes to be scored   (32 eliminated genes)

    ## 
    ##   Level 10:  24 nodes to be scored   (99 eliminated genes)

    ## 
    ##   Level 9:   43 nodes to be scored   (290 eliminated genes)

    ## 
    ##   Level 8:   59 nodes to be scored   (406 eliminated genes)

    ## 
    ##   Level 7:   85 nodes to be scored   (499 eliminated genes)

    ## 
    ##   Level 6:   91 nodes to be scored   (660 eliminated genes)

    ## 
    ##   Level 5:   91 nodes to be scored   (999 eliminated genes)

    ## 
    ##   Level 4:   55 nodes to be scored   (1103 eliminated genes)

    ## 
    ##   Level 3:   42 nodes to be scored   (1270 eliminated genes)

    ## 
    ##   Level 2:   11 nodes to be scored   (1317 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1400 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 357 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  5 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   13 nodes to be scored   (12 eliminated genes)

    ## 
    ##   Level 8:   30 nodes to be scored   (153 eliminated genes)

    ## 
    ##   Level 7:   54 nodes to be scored   (459 eliminated genes)

    ## 
    ##   Level 6:   76 nodes to be scored   (669 eliminated genes)

    ## 
    ##   Level 5:   66 nodes to be scored   (1173 eliminated genes)

    ## 
    ##   Level 4:   64 nodes to be scored   (1559 eliminated genes)

    ## 
    ##   Level 3:   33 nodes to be scored   (2124 eliminated genes)

    ## 
    ##   Level 2:   13 nodes to be scored   (2448 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2680 eliminated genes)

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
    ## 10 GO:0004342             glucosamine-6-phosphate deaminase activity         2
    ## 11 GO:0004571 mannosyl-oligosaccharide 1,2-alpha-mannosidase acti...         2
    ## 12 GO:0001965                        G-protein alpha-subunit binding        16
    ## 13 GO:0005085             guanyl-nucleotide exchange factor activity        12
    ##    Significant Expected    Fisher               type
    ## 1           19     8.84 0.0005700 Biological.Process
    ## 2            9     2.72 0.0005800 Biological.Process
    ## 3            5     1.50 0.0102300 Biological.Process
    ## 4           22    14.15 0.0181700 Biological.Process
    ## 5           66    39.17 0.0000034 Molecular.Function
    ## 6            7     1.79 0.0005900 Molecular.Function
    ## 7           29    18.17 0.0053800 Molecular.Function
    ## 8            3     0.60 0.0116800 Molecular.Function
    ## 9            3     0.60 0.0116800 Molecular.Function
    ## 10           2     0.30 0.0221400 Molecular.Function
    ## 11           2     0.30 0.0221400 Molecular.Function
    ## 12           6     2.38 0.0224400 Molecular.Function
    ## 13           5     1.79 0.0229800 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 2 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 295 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  4 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 11:  7 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 10:  12 nodes to be scored   (82 eliminated genes)

    ## 
    ##   Level 9:   19 nodes to be scored   (238 eliminated genes)

    ## 
    ##   Level 8:   28 nodes to be scored   (328 eliminated genes)

    ## 
    ##   Level 7:   40 nodes to be scored   (379 eliminated genes)

    ## 
    ##   Level 6:   43 nodes to be scored   (498 eliminated genes)

    ## 
    ##   Level 5:   52 nodes to be scored   (815 eliminated genes)

    ## 
    ##   Level 4:   44 nodes to be scored   (909 eliminated genes)

    ## 
    ##   Level 3:   32 nodes to be scored   (1090 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (1296 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1376 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 261 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   9 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 8:   22 nodes to be scored   (101 eliminated genes)

    ## 
    ##   Level 7:   41 nodes to be scored   (408 eliminated genes)

    ## 
    ##   Level 6:   63 nodes to be scored   (616 eliminated genes)

    ## 
    ##   Level 5:   46 nodes to be scored   (1115 eliminated genes)

    ## 
    ##   Level 4:   43 nodes to be scored   (1465 eliminated genes)

    ## 
    ##   Level 3:   24 nodes to be scored   (2013 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (2286 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2598 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0000226                  microtubule cytoskeleton organization        43
    ## 2 GO:0000390                       spliceosomal complex disassembly         1
    ## 3 GO:0000052                           citrulline metabolic process         1
    ## 4 GO:0000349 generation of catalytic spliceosome for first trans...         1
    ## 5 GO:0006813                                potassium ion transport         1
    ## 6 GO:0005272                                sodium channel activity         2
    ## 7 GO:0003985                acetyl-CoA C-acetyltransferase activity         2
    ## 8 GO:0005524                                            ATP binding       263
    ##   Significant Expected Fisher               type
    ## 1           6     1.81 0.0280 Biological.Process
    ## 2           1     0.04 0.0420 Biological.Process
    ## 3           1     0.04 0.0420 Biological.Process
    ## 4           1     0.04 0.0420 Biological.Process
    ## 5           1     0.04 0.0420 Biological.Process
    ## 6           2     0.11 0.0033 Molecular.Function
    ## 7           2     0.11 0.0033 Molecular.Function
    ## 8          26    15.08 0.0033 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 153 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  5 nodes to be scored    (63 eliminated genes)

    ## 
    ##   Level 9:   9 nodes to be scored    (203 eliminated genes)

    ## 
    ##   Level 8:   13 nodes to be scored   (246 eliminated genes)

    ## 
    ##   Level 7:   12 nodes to be scored   (269 eliminated genes)

    ## 
    ##   Level 6:   20 nodes to be scored   (368 eliminated genes)

    ## 
    ##   Level 5:   27 nodes to be scored   (404 eliminated genes)

    ## 
    ##   Level 4:   29 nodes to be scored   (474 eliminated genes)

    ## 
    ##   Level 3:   24 nodes to be scored   (748 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (910 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1139 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 144 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   6 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   9 nodes to be scored    (102 eliminated genes)

    ## 
    ##   Level 7:   23 nodes to be scored   (391 eliminated genes)

    ## 
    ##   Level 6:   29 nodes to be scored   (574 eliminated genes)

    ## 
    ##   Level 5:   30 nodes to be scored   (949 eliminated genes)

    ## 
    ##   Level 4:   24 nodes to be scored   (1119 eliminated genes)

    ## 
    ##   Level 3:   13 nodes to be scored   (1790 eliminated genes)

    ## 
    ##   Level 2:   7 nodes to be scored    (2054 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2397 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0001558                              regulation of cell growth         3
    ## 2  GO:0003844             1,4-alpha-glucan branching enzyme activity         1
    ## 3  GO:0004325                                ferrochelatase activity         1
    ## 4  GO:0004095              carnitine O-palmitoyltransferase activity         1
    ## 5  GO:0004777 succinate-semialdehyde dehydrogenase (NAD+) activit...         1
    ## 6  GO:0004618                       phosphoglycerate kinase activity         1
    ## 7  GO:0004066 asparagine synthase (glutamine-hydrolyzing) activit...         1
    ## 8  GO:0005391 P-type sodium:potassium-exchanging transporter acti...         2
    ## 9  GO:0003730                                    mRNA 3'-UTR binding         2
    ## 10 GO:0003987                            acetate-CoA ligase activity         2
    ## 11 GO:0004029                 aldehyde dehydrogenase (NAD+) activity         3
    ## 12 GO:0004031                              aldehyde oxidase activity         3
    ##    Significant Expected Fisher               type
    ## 1            1     0.03  0.028 Biological.Process
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

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 601 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  7 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 12:  9 nodes to be scored    (20 eliminated genes)

    ## 
    ##   Level 11:  14 nodes to be scored   (43 eliminated genes)

    ## 
    ##   Level 10:  32 nodes to be scored   (105 eliminated genes)

    ## 
    ##   Level 9:   50 nodes to be scored   (281 eliminated genes)

    ## 
    ##   Level 8:   71 nodes to be scored   (404 eliminated genes)

    ## 
    ##   Level 7:   91 nodes to be scored   (501 eliminated genes)

    ## 
    ##   Level 6:   103 nodes to be scored  (661 eliminated genes)

    ## 
    ##   Level 5:   96 nodes to be scored   (1001 eliminated genes)

    ## 
    ##   Level 4:   60 nodes to be scored   (1068 eliminated genes)

    ## 
    ##   Level 3:   49 nodes to be scored   (1208 eliminated genes)

    ## 
    ##   Level 2:   13 nodes to be scored   (1326 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1401 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 430 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  9 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 9:   19 nodes to be scored   (4 eliminated genes)

    ## 
    ##   Level 8:   42 nodes to be scored   (157 eliminated genes)

    ## 
    ##   Level 7:   72 nodes to be scored   (473 eliminated genes)

    ## 
    ##   Level 6:   95 nodes to be scored   (707 eliminated genes)

    ## 
    ##   Level 5:   71 nodes to be scored   (1190 eliminated genes)

    ## 
    ##   Level 4:   69 nodes to be scored   (1559 eliminated genes)

    ## 
    ##   Level 3:   35 nodes to be scored   (2137 eliminated genes)

    ## 
    ##   Level 2:   14 nodes to be scored   (2457 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2678 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0003341                                        cilium movement        20
    ## 2  GO:0000184 nuclear-transcribed mRNA catabolic process, nonsens...        10
    ## 3  GO:0002064                            epithelial cell development         4
    ## 4  GO:0000045                                 autophagosome assembly         7
    ## 5  GO:0002790                                      peptide secretion         2
    ## 6  GO:0007224                           smoothened signaling pathway         2
    ## 7  GO:0001933         negative regulation of protein phosphorylation         5
    ## 8  GO:0000038           very long-chain fatty acid metabolic process         9
    ## 9  GO:0005388                    P-type calcium transporter activity         4
    ## 10 GO:0003964                   RNA-directed DNA polymerase activity       122
    ## 11 GO:0003676                                   nucleic acid binding       497
    ## 12 GO:0005044                            scavenger receptor activity        13
    ## 13 GO:0005524                                            ATP binding       263
    ## 14 GO:0004100                               chitin synthase activity         2
    ## 15 GO:0003730                                    mRNA 3'-UTR binding         2
    ## 16 GO:0004573 Glc3Man9GlcNAc2 oligosaccharide glucosidase activit...         2
    ## 17 GO:0004345             glucose-6-phosphate dehydrogenase activity         2
    ## 18 GO:0003756                   protein disulfide isomerase activity         5
    ## 19 GO:0003723                                            RNA binding       128
    ##    Significant Expected Fisher               type
    ## 1            9     3.22 0.0059 Biological.Process
    ## 2            5     1.61 0.0130 Biological.Process
    ## 3            3     0.64 0.0145 Biological.Process
    ## 4            4     1.13 0.0153 Biological.Process
    ## 5            2     0.32 0.0258 Biological.Process
    ## 6            2     0.32 0.0258 Biological.Process
    ## 7            3     0.80 0.0320 Biological.Process
    ## 8            4     1.45 0.0423 Biological.Process
    ## 9            4     0.73 0.0011 Molecular.Function
    ## 10          35    22.39 0.0029 Molecular.Function
    ## 11         108    91.23 0.0033 Molecular.Function
    ## 12           6     2.39 0.0199 Molecular.Function
    ## 13          61    48.28 0.0224 Molecular.Function
    ## 14           2     0.37 0.0336 Molecular.Function
    ## 15           2     0.37 0.0336 Molecular.Function
    ## 16           2     0.37 0.0336 Molecular.Function
    ## 17           2     0.37 0.0336 Molecular.Function
    ## 18           3     0.92 0.0459 Molecular.Function
    ## 19          29    23.50 0.0461 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 2 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 375 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  5 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  7 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 10:  17 nodes to be scored   (88 eliminated genes)

    ## 
    ##   Level 9:   28 nodes to be scored   (240 eliminated genes)

    ## 
    ##   Level 8:   35 nodes to be scored   (339 eliminated genes)

    ## 
    ##   Level 7:   53 nodes to be scored   (415 eliminated genes)

    ## 
    ##   Level 6:   65 nodes to be scored   (558 eliminated genes)

    ## 
    ##   Level 5:   65 nodes to be scored   (907 eliminated genes)

    ## 
    ##   Level 4:   50 nodes to be scored   (1016 eliminated genes)

    ## 
    ##   Level 3:   39 nodes to be scored   (1177 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (1318 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1384 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 280 nontrivial nodes
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
    ##   Level 8:   19 nodes to be scored   (145 eliminated genes)

    ## 
    ##   Level 7:   41 nodes to be scored   (453 eliminated genes)

    ## 
    ##   Level 6:   61 nodes to be scored   (609 eliminated genes)

    ## 
    ##   Level 5:   53 nodes to be scored   (1100 eliminated genes)

    ## 
    ##   Level 4:   50 nodes to be scored   (1410 eliminated genes)

    ## 
    ##   Level 3:   27 nodes to be scored   (1992 eliminated genes)

    ## 
    ##   Level 2:   13 nodes to be scored   (2328 eliminated genes)

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
    ## 1           2     0.27 0.025000 Biological.Process
    ## 2          41    21.15 0.000012 Molecular.Function
    ## 3           9     2.41 0.000380 Molecular.Function
    ## 4          19     9.81 0.003310 Molecular.Function
    ## 5          18    18.01 0.015680 Molecular.Function
    ## 6           2     0.24 0.018290 Molecular.Function
    ## 7          42    39.97 0.019330 Molecular.Function
    ## 8           4     1.29 0.034290 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 293 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  5 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  10 nodes to be scored   (2 eliminated genes)

    ## 
    ##   Level 10:  14 nodes to be scored   (71 eliminated genes)

    ## 
    ##   Level 9:   17 nodes to be scored   (223 eliminated genes)

    ## 
    ##   Level 8:   29 nodes to be scored   (264 eliminated genes)

    ## 
    ##   Level 7:   41 nodes to be scored   (338 eliminated genes)

    ## 
    ##   Level 6:   48 nodes to be scored   (460 eliminated genes)

    ## 
    ##   Level 5:   50 nodes to be scored   (705 eliminated genes)

    ## 
    ##   Level 4:   36 nodes to be scored   (784 eliminated genes)

    ## 
    ##   Level 3:   31 nodes to be scored   (1065 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (1255 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1336 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

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
    ##   Level 9:   5 nodes to be scored    (9 eliminated genes)

    ## 
    ##   Level 8:   10 nodes to be scored   (142 eliminated genes)

    ## 
    ##   Level 7:   23 nodes to be scored   (436 eliminated genes)

    ## 
    ##   Level 6:   32 nodes to be scored   (504 eliminated genes)

    ## 
    ##   Level 5:   34 nodes to be scored   (956 eliminated genes)

    ## 
    ##   Level 4:   29 nodes to be scored   (1253 eliminated genes)

    ## 
    ##   Level 3:   16 nodes to be scored   (1893 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (2135 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2451 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0003341                                        cilium movement        20
    ## 2  GO:0015969             guanosine tetraphosphate metabolic process         9
    ## 3  GO:0001731         formation of translation preinitiation complex         1
    ## 4  GO:0000390                       spliceosomal complex disassembly         1
    ## 5  GO:0000972 transcription-dependent tethering of RNA polymerase...         1
    ## 6  GO:0007175 negative regulation of epidermal growth factor-acti...         1
    ## 7  GO:0006884                                cell volume homeostasis         1
    ## 8  GO:0006513                             protein monoubiquitination         1
    ## 9  GO:0002291 T cell activation via T cell receptor contact with ...        10
    ## 10 GO:0000082                  G1/S transition of mitotic cell cycle       104
    ## 11 GO:0001409 guanine nucleotide transmembrane transporter activi...        30
    ## 12 GO:0004315     3-oxoacyl-[acyl-carrier-protein] synthase activity         2
    ## 13 GO:0003964                   RNA-directed DNA polymerase activity       122
    ## 14 GO:0005524                                            ATP binding       263
    ## 15 GO:0003676                                   nucleic acid binding       497
    ## 16 GO:0005080                               protein kinase C binding         1
    ## 17 GO:0004142      diacylglycerol cholinephosphotransferase activity         1
    ## 18 GO:0004013                        adenosylhomocysteinase activity         1
    ## 19 GO:0005338    nucleotide-sugar transmembrane transporter activity         1
    ## 20 GO:0003863 3-methyl-2-oxobutanoate dehydrogenase (2-methylprop...         1
    ## 21 GO:0004591 oxoglutarate dehydrogenase (succinyl-transferring) ...         1
    ## 22 GO:0001758                         retinal dehydrogenase activity         1
    ##    Significant Expected Fisher               type
    ## 1            3     0.63 0.0230 Biological.Process
    ## 2            2     0.28 0.0300 Biological.Process
    ## 3            1     0.03 0.0310 Biological.Process
    ## 4            1     0.03 0.0310 Biological.Process
    ## 5            1     0.03 0.0310 Biological.Process
    ## 6            1     0.03 0.0310 Biological.Process
    ## 7            1     0.03 0.0310 Biological.Process
    ## 8            1     0.03 0.0310 Biological.Process
    ## 9            2     0.31 0.0370 Biological.Process
    ## 10           7     3.26 0.0390 Biological.Process
    ## 11           6     1.25 0.0012 Molecular.Function
    ## 12           2     0.08 0.0017 Molecular.Function
    ## 13          12     5.10 0.0043 Molecular.Function
    ## 14          18    11.00 0.0229 Molecular.Function
    ## 15          21    20.79 0.0371 Molecular.Function
    ## 16           1     0.04 0.0418 Molecular.Function
    ## 17           1     0.04 0.0418 Molecular.Function
    ## 18           1     0.04 0.0418 Molecular.Function
    ## 19           1     0.04 0.0418 Molecular.Function
    ## 20           1     0.04 0.0418 Molecular.Function
    ## 21           1     0.04 0.0418 Molecular.Function
    ## 22           1     0.04 0.0418 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 256 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  6 nodes to be scored    (23 eliminated genes)

    ## 
    ##   Level 10:  9 nodes to be scored    (88 eliminated genes)

    ## 
    ##   Level 9:   15 nodes to be scored   (253 eliminated genes)

    ## 
    ##   Level 8:   21 nodes to be scored   (307 eliminated genes)

    ## 
    ##   Level 7:   32 nodes to be scored   (347 eliminated genes)

    ## 
    ##   Level 6:   42 nodes to be scored   (448 eliminated genes)

    ## 
    ##   Level 5:   49 nodes to be scored   (729 eliminated genes)

    ## 
    ##   Level 4:   35 nodes to be scored   (833 eliminated genes)

    ## 
    ##   Level 3:   30 nodes to be scored   (1042 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (1162 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1364 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 103 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   12 nodes to be scored   (263 eliminated genes)

    ## 
    ##   Level 6:   18 nodes to be scored   (331 eliminated genes)

    ## 
    ##   Level 5:   21 nodes to be scored   (653 eliminated genes)

    ## 
    ##   Level 4:   21 nodes to be scored   (961 eliminated genes)

    ## 
    ##   Level 3:   18 nodes to be scored   (1652 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (1960 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2504 eliminated genes)

    ##         GO.ID                                       Term Annotated Significant
    ## 1  GO:0001822                         kidney development        65           7
    ## 2  GO:0001818 negative regulation of cytokine production        18           3
    ## 3  GO:0001913               T cell mediated cytotoxicity         1           1
    ## 4  GO:0001696                     gastric acid secretion        10           2
    ## 5  GO:0005524                                ATP binding       263          16
    ## 6  GO:0004017                  adenylate kinase activity         1           1
    ## 7  GO:0003682                          chromatin binding        13           2
    ## 8  GO:0004100                   chitin synthase activity         2           1
    ## 9  GO:0004660       protein farnesyltransferase activity         2           1
    ## 10 GO:0001965            G-protein alpha-subunit binding        16           2
    ##    Expected  Fisher               type
    ## 1      1.90 0.00210 Biological.Process
    ## 2      0.53 0.01390 Biological.Process
    ## 3      0.03 0.02920 Biological.Process
    ## 4      0.29 0.03230 Biological.Process
    ## 5      5.88 0.00013 Molecular.Function
    ## 6      0.02 0.02236 Molecular.Function
    ## 7      0.29 0.03275 Molecular.Function
    ## 8      0.04 0.04422 Molecular.Function
    ## 9      0.04 0.04422 Molecular.Function
    ## 10     0.36 0.04827 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 254 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  4 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 11:  7 nodes to be scored    (2 eliminated genes)

    ## 
    ##   Level 10:  9 nodes to be scored    (79 eliminated genes)

    ## 
    ##   Level 9:   15 nodes to be scored   (112 eliminated genes)

    ## 
    ##   Level 8:   22 nodes to be scored   (169 eliminated genes)

    ## 
    ##   Level 7:   30 nodes to be scored   (212 eliminated genes)

    ## 
    ##   Level 6:   41 nodes to be scored   (428 eliminated genes)

    ## 
    ##   Level 5:   52 nodes to be scored   (643 eliminated genes)

    ## 
    ##   Level 4:   37 nodes to be scored   (724 eliminated genes)

    ## 
    ##   Level 3:   25 nodes to be scored   (1026 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (1208 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1284 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 151 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   7 nodes to be scored    (7 eliminated genes)

    ## 
    ##   Level 7:   21 nodes to be scored   (425 eliminated genes)

    ## 
    ##   Level 6:   34 nodes to be scored   (533 eliminated genes)

    ## 
    ##   Level 5:   30 nodes to be scored   (909 eliminated genes)

    ## 
    ##   Level 4:   28 nodes to be scored   (1257 eliminated genes)

    ## 
    ##   Level 3:   15 nodes to be scored   (1823 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (2106 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2420 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0015969             guanosine tetraphosphate metabolic process         9
    ## 2  GO:0000349 generation of catalytic spliceosome for first trans...         1
    ## 3  GO:0016311                                      dephosphorylation         1
    ## 4  GO:0001731         formation of translation preinitiation complex         1
    ## 5  GO:0006091         generation of precursor metabolites and energy        12
    ## 6  GO:0001569       branching involved in blood vessel morphogenesis         2
    ## 7  GO:0002790                                      peptide secretion         2
    ## 8  GO:0000079 regulation of cyclin-dependent protein serine/threo...         2
    ## 9  GO:0005302          L-tyrosine transmembrane transporter activity         9
    ## 10 GO:0003824                                     catalytic activity      1038
    ## 11 GO:0000822                      inositol hexakisphosphate binding         1
    ## 12 GO:0004140                          dephospho-CoA kinase activity         1
    ## 13 GO:0001054                              RNA polymerase I activity         1
    ## 14 GO:0004334                           fumarylacetoacetase activity         1
    ## 15 GO:0004095              carnitine O-palmitoyltransferase activity         1
    ##    Significant Expected Fisher               type
    ## 1            2     0.19 0.0140 Biological.Process
    ## 2            1     0.02 0.0210 Biological.Process
    ## 3            1     0.02 0.0210 Biological.Process
    ## 4            1     0.02 0.0210 Biological.Process
    ## 5            2     0.26 0.0260 Biological.Process
    ## 6            1     0.04 0.0420 Biological.Process
    ## 7            1     0.04 0.0420 Biological.Process
    ## 8            1     0.04 0.0420 Biological.Process
    ## 9            3     0.24 0.0014 Molecular.Function
    ## 10          30    28.07 0.0122 Molecular.Function
    ## 11           1     0.03 0.0270 Molecular.Function
    ## 12           1     0.03 0.0270 Molecular.Function
    ## 13           1     0.03 0.0270 Molecular.Function
    ## 14           1     0.03 0.0270 Molecular.Function
    ## 15           1     0.03 0.0270 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 333 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  4 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 11:  8 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 10:  17 nodes to be scored   (78 eliminated genes)

    ## 
    ##   Level 9:   26 nodes to be scored   (251 eliminated genes)

    ## 
    ##   Level 8:   39 nodes to be scored   (334 eliminated genes)

    ## 
    ##   Level 7:   47 nodes to be scored   (373 eliminated genes)

    ## 
    ##   Level 6:   49 nodes to be scored   (502 eliminated genes)

    ## 
    ##   Level 5:   56 nodes to be scored   (792 eliminated genes)

    ## 
    ##   Level 4:   42 nodes to be scored   (836 eliminated genes)

    ## 
    ##   Level 3:   32 nodes to be scored   (1044 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (1259 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1367 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 212 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   6 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   8 nodes to be scored    (100 eliminated genes)

    ## 
    ##   Level 7:   30 nodes to be scored   (405 eliminated genes)

    ## 
    ##   Level 6:   47 nodes to be scored   (524 eliminated genes)

    ## 
    ##   Level 5:   48 nodes to be scored   (1047 eliminated genes)

    ## 
    ##   Level 4:   44 nodes to be scored   (1325 eliminated genes)

    ## 
    ##   Level 3:   18 nodes to be scored   (1941 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (2254 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2558 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0002218                   activation of innate immune response        26
    ## 2  GO:0006397                                        mRNA processing        64
    ## 3  GO:0001731         formation of translation preinitiation complex         1
    ## 4  GO:0006601                          creatine biosynthetic process         1
    ## 5  GO:0001113 transcription open complex formation at RNA polymer...         1
    ## 6  GO:0002265       astrocyte activation involved in immune response         1
    ## 7  GO:0009071             serine family amino acid catabolic process         1
    ## 8  GO:0002230 positive regulation of defense response to virus by...         1
    ## 9  GO:0000349 generation of catalytic spliceosome for first trans...         1
    ## 10 GO:0006813                                potassium ion transport         1
    ## 11 GO:0016311                                      dephosphorylation         1
    ## 12 GO:0003730                                    mRNA 3'-UTR binding         2
    ## 13 GO:0003676                                   nucleic acid binding       497
    ## 14 GO:0000016                                       lactase activity         4
    ## 15 GO:0004222                          metalloendopeptidase activity        15
    ## 16 GO:0003997                              acyl-CoA oxidase activity         1
    ## 17 GO:0004325                                ferrochelatase activity         1
    ## 18 GO:0004516          nicotinate phosphoribosyltransferase activity         1
    ## 19 GO:0002161                        aminoacyl-tRNA editing activity         1
    ## 20 GO:0004458          D-lactate dehydrogenase (cytochrome) activity         1
    ## 21 GO:0008810                                     cellulase activity         1
    ## 22 GO:0004777 succinate-semialdehyde dehydrogenase (NAD+) activit...         1
    ## 23 GO:0004334                           fumarylacetoacetase activity         1
    ## 24 GO:0004618                       phosphoglycerate kinase activity         1
    ## 25 GO:0004642    phosphoribosylformylglycinamidine synthase activity         1
    ## 26 GO:0008137               NADH dehydrogenase (ubiquinone) activity         1
    ## 27 GO:0004637            phosphoribosylamine-glycine ligase activity         1
    ## 28 GO:0005524                                            ATP binding       263
    ##    Significant Expected Fisher               type
    ## 1            4     0.91 0.0110 Biological.Process
    ## 2            4     2.23 0.0340 Biological.Process
    ## 3            1     0.03 0.0350 Biological.Process
    ## 4            1     0.03 0.0350 Biological.Process
    ## 5            1     0.03 0.0350 Biological.Process
    ## 6            1     0.03 0.0350 Biological.Process
    ## 7            1     0.03 0.0350 Biological.Process
    ## 8            1     0.03 0.0350 Biological.Process
    ## 9            1     0.03 0.0350 Biological.Process
    ## 10           1     0.03 0.0350 Biological.Process
    ## 11           1     0.03 0.0350 Biological.Process
    ## 12           2     0.08 0.0017 Molecular.Function
    ## 13          28    20.79 0.0097 Molecular.Function
    ## 14           2     0.17 0.0099 Molecular.Function
    ## 15           3     0.63 0.0225 Molecular.Function
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
    ## 28          17    11.00 0.0435 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 483 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  8 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  7 nodes to be scored    (15 eliminated genes)

    ## 
    ##   Level 11:  11 nodes to be scored   (44 eliminated genes)

    ## 
    ##   Level 10:  24 nodes to be scored   (98 eliminated genes)

    ## 
    ##   Level 9:   38 nodes to be scored   (264 eliminated genes)

    ## 
    ##   Level 8:   55 nodes to be scored   (374 eliminated genes)

    ## 
    ##   Level 7:   68 nodes to be scored   (458 eliminated genes)

    ## 
    ##   Level 6:   79 nodes to be scored   (591 eliminated genes)

    ## 
    ##   Level 5:   81 nodes to be scored   (893 eliminated genes)

    ## 
    ##   Level 4:   55 nodes to be scored   (984 eliminated genes)

    ## 
    ##   Level 3:   44 nodes to be scored   (1157 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (1300 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1388 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 292 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  5 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   11 nodes to be scored   (2 eliminated genes)

    ## 
    ##   Level 8:   21 nodes to be scored   (149 eliminated genes)

    ## 
    ##   Level 7:   41 nodes to be scored   (450 eliminated genes)

    ## 
    ##   Level 6:   63 nodes to be scored   (620 eliminated genes)

    ## 
    ##   Level 5:   59 nodes to be scored   (1066 eliminated genes)

    ## 
    ##   Level 4:   49 nodes to be scored   (1406 eliminated genes)

    ## 
    ##   Level 3:   29 nodes to be scored   (2095 eliminated genes)

    ## 
    ##   Level 2:   12 nodes to be scored   (2388 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2630 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0001732 formation of cytoplasmic translation initiation com...         6
    ## 2 GO:0000054                  ribosomal subunit export from nucleus         3
    ## 3 GO:0003676                                   nucleic acid binding       497
    ## 4 GO:0003727                            single-stranded RNA binding         2
    ## 5 GO:0003723                                            RNA binding       128
    ## 6 GO:0003964                   RNA-directed DNA polymerase activity       122
    ## 7 GO:0005198                           structural molecule activity        54
    ##   Significant Expected    Fisher               type
    ## 1           4     0.45 0.0004100 Biological.Process
    ## 2           2     0.23 0.0161100 Biological.Process
    ## 3          59    39.07 0.0000035 Molecular.Function
    ## 4           2     0.16 0.0062000 Molecular.Function
    ## 5          20    10.06 0.0075000 Molecular.Function
    ## 6          17     9.59 0.0131000 Molecular.Function
    ## 7           9     4.25 0.0345000 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 145 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  6 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   8 nodes to be scored    (166 eliminated genes)

    ## 
    ##   Level 8:   11 nodes to be scored   (227 eliminated genes)

    ## 
    ##   Level 7:   16 nodes to be scored   (253 eliminated genes)

    ## 
    ##   Level 6:   22 nodes to be scored   (308 eliminated genes)

    ## 
    ##   Level 5:   30 nodes to be scored   (495 eliminated genes)

    ## 
    ##   Level 4:   23 nodes to be scored   (542 eliminated genes)

    ## 
    ##   Level 3:   17 nodes to be scored   (782 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (1029 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1222 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 76 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   7 nodes to be scored    (2 eliminated genes)

    ## 
    ##   Level 6:   10 nodes to be scored   (10 eliminated genes)

    ## 
    ##   Level 5:   15 nodes to be scored   (124 eliminated genes)

    ## 
    ##   Level 4:   19 nodes to be scored   (349 eliminated genes)

    ## 
    ##   Level 3:   14 nodes to be scored   (1332 eliminated genes)

    ## 
    ##   Level 2:   4 nodes to be scored    (1892 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2430 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0006884                                cell volume homeostasis         1
    ## 2 GO:0001510                                        RNA methylation         4
    ## 3 GO:0001508                                       action potential         4
    ## 4 GO:0004095              carnitine O-palmitoyltransferase activity         1
    ## 5 GO:0001055                             RNA polymerase II activity         1
    ## 6 GO:0004591 oxoglutarate dehydrogenase (succinyl-transferring) ...         1
    ## 7 GO:0004494                      methylmalonyl-CoA mutase activity         1
    ## 8 GO:0004435          phosphatidylinositol phospholipase C activity         1
    ## 9 GO:0016491                                oxidoreductase activity       114
    ##   Significant Expected Fisher               type
    ## 1           1     0.01 0.0085 Biological.Process
    ## 2           1     0.03 0.0338 Biological.Process
    ## 3           1     0.03 0.0338 Biological.Process
    ## 4           1     0.01 0.0090 Molecular.Function
    ## 5           1     0.01 0.0090 Molecular.Function
    ## 6           1     0.01 0.0090 Molecular.Function
    ## 7           1     0.01 0.0090 Molecular.Function
    ## 8           1     0.01 0.0090 Molecular.Function
    ## 9           5     1.03 0.0210 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 220 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  3 nodes to be scored    (11 eliminated genes)

    ## 
    ##   Level 10:  5 nodes to be scored    (18 eliminated genes)

    ## 
    ##   Level 9:   9 nodes to be scored    (173 eliminated genes)

    ## 
    ##   Level 8:   22 nodes to be scored   (204 eliminated genes)

    ## 
    ##   Level 7:   27 nodes to be scored   (229 eliminated genes)

    ## 
    ##   Level 6:   37 nodes to be scored   (432 eliminated genes)

    ## 
    ##   Level 5:   45 nodes to be scored   (609 eliminated genes)

    ## 
    ##   Level 4:   33 nodes to be scored   (727 eliminated genes)

    ## 
    ##   Level 3:   25 nodes to be scored   (974 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (1223 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1358 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 196 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   6 nodes to be scored    (2 eliminated genes)

    ## 
    ##   Level 8:   11 nodes to be scored   (151 eliminated genes)

    ## 
    ##   Level 7:   28 nodes to be scored   (436 eliminated genes)

    ## 
    ##   Level 6:   38 nodes to be scored   (563 eliminated genes)

    ## 
    ##   Level 5:   41 nodes to be scored   (924 eliminated genes)

    ## 
    ##   Level 4:   38 nodes to be scored   (1270 eliminated genes)

    ## 
    ##   Level 3:   20 nodes to be scored   (2010 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (2285 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2575 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0002221         pattern recognition receptor signaling pathway        11
    ## 2  GO:0001666                                    response to hypoxia        14
    ## 3  GO:0006556              S-adenosylmethionine biosynthetic process         1
    ## 4  GO:0001937 negative regulation of endothelial cell proliferati...         1
    ## 5  GO:0005198                           structural molecule activity        54
    ## 6  GO:0003756                   protein disulfide isomerase activity         5
    ## 7  GO:0004321                       fatty-acyl-CoA synthase activity         1
    ## 8  GO:0004561                 alpha-N-acetylglucosaminidase activity         1
    ## 9  GO:0004568                                     chitinase activity         1
    ## 10 GO:0004332                fructose-bisphosphate aldolase activity         1
    ## 11 GO:0001594                          trace-amine receptor activity         1
    ## 12 GO:0004305                           ethanolamine kinase activity         1
    ## 13 GO:0004352                glutamate dehydrogenase (NAD+) activity         1
    ##    Significant Expected   Fisher               type
    ## 1            5     0.32 6.70e-06 Biological.Process
    ## 2            3     0.41 6.70e-03 Biological.Process
    ## 3            1     0.03 2.92e-02 Biological.Process
    ## 4            1     0.03 2.92e-02 Biological.Process
    ## 5            6     1.85 7.60e-03 Molecular.Function
    ## 6            2     0.17 1.09e-02 Molecular.Function
    ## 7            1     0.03 3.43e-02 Molecular.Function
    ## 8            1     0.03 3.43e-02 Molecular.Function
    ## 9            1     0.03 3.43e-02 Molecular.Function
    ## 10           1     0.03 3.43e-02 Molecular.Function
    ## 11           1     0.03 3.43e-02 Molecular.Function
    ## 12           1     0.03 3.43e-02 Molecular.Function
    ## 13           1     0.03 3.43e-02 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 285 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  1 nodes to be scored    (11 eliminated genes)

    ## 
    ##   Level 11:  5 nodes to be scored    (15 eliminated genes)

    ## 
    ##   Level 10:  12 nodes to be scored   (63 eliminated genes)

    ## 
    ##   Level 9:   18 nodes to be scored   (242 eliminated genes)

    ## 
    ##   Level 8:   25 nodes to be scored   (328 eliminated genes)

    ## 
    ##   Level 7:   36 nodes to be scored   (397 eliminated genes)

    ## 
    ##   Level 6:   48 nodes to be scored   (515 eliminated genes)

    ## 
    ##   Level 5:   52 nodes to be scored   (799 eliminated genes)

    ## 
    ##   Level 4:   41 nodes to be scored   (871 eliminated genes)

    ## 
    ##   Level 3:   33 nodes to be scored   (1066 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (1260 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1367 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 233 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   10 nodes to be scored   (9 eliminated genes)

    ## 
    ##   Level 8:   18 nodes to be scored   (151 eliminated genes)

    ## 
    ##   Level 7:   30 nodes to be scored   (459 eliminated genes)

    ## 
    ##   Level 6:   49 nodes to be scored   (583 eliminated genes)

    ## 
    ##   Level 5:   43 nodes to be scored   (937 eliminated genes)

    ## 
    ##   Level 4:   43 nodes to be scored   (1312 eliminated genes)

    ## 
    ##   Level 3:   22 nodes to be scored   (1863 eliminated genes)

    ## 
    ##   Level 2:   12 nodes to be scored   (2193 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2605 eliminated genes)

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
    ##    Significant Expected   Fisher               type
    ## 1            5     0.43 0.000031 Biological.Process
    ## 2            2     0.14 0.007100 Biological.Process
    ## 3            2     0.28 0.030300 Biological.Process
    ## 4            1     0.04 0.035600 Biological.Process
    ## 5            1     0.04 0.035600 Biological.Process
    ## 6            1     0.04 0.035600 Biological.Process
    ## 7            2     0.16 0.008200 Molecular.Function
    ## 8            6     1.28 0.010300 Molecular.Function
    ## 9            5     1.60 0.019500 Molecular.Function
    ## 10           7     2.83 0.020900 Molecular.Function
    ## 11          23    26.53 0.038400 Molecular.Function
    ## 12           2     0.37 0.049800 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 759 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  5 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  9 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 12:  15 nodes to be scored   (23 eliminated genes)

    ## 
    ##   Level 11:  24 nodes to be scored   (47 eliminated genes)

    ## 
    ##   Level 10:  42 nodes to be scored   (114 eliminated genes)

    ## 
    ##   Level 9:   65 nodes to be scored   (293 eliminated genes)

    ## 
    ##   Level 8:   94 nodes to be scored   (425 eliminated genes)

    ## 
    ##   Level 7:   115 nodes to be scored  (542 eliminated genes)

    ## 
    ##   Level 6:   130 nodes to be scored  (750 eliminated genes)

    ## 
    ##   Level 5:   118 nodes to be scored  (1077 eliminated genes)

    ## 
    ##   Level 4:   75 nodes to be scored   (1144 eliminated genes)

    ## 
    ##   Level 3:   52 nodes to be scored   (1277 eliminated genes)

    ## 
    ##   Level 2:   12 nodes to be scored   (1361 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1402 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 496 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  5 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  10 nodes to be scored   (1 eliminated genes)

    ## 
    ##   Level 9:   22 nodes to be scored   (16 eliminated genes)

    ## 
    ##   Level 8:   46 nodes to be scored   (162 eliminated genes)

    ## 
    ##   Level 7:   85 nodes to be scored   (480 eliminated genes)

    ## 
    ##   Level 6:   115 nodes to be scored  (703 eliminated genes)

    ## 
    ##   Level 5:   88 nodes to be scored   (1233 eliminated genes)

    ## 
    ##   Level 4:   74 nodes to be scored   (1628 eliminated genes)

    ## 
    ##   Level 3:   35 nodes to be scored   (2193 eliminated genes)

    ## 
    ##   Level 2:   14 nodes to be scored   (2507 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2678 eliminated genes)

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
    ## 13 GO:0004623                              phospholipase A2 activity         4
    ## 14 GO:0005385            zinc ion transmembrane transporter activity         4
    ## 15 GO:0003727                            single-stranded RNA binding         2
    ## 16 GO:0004365 glyceraldehyde-3-phosphate dehydrogenase (NAD+) (ph...         2
    ## 17 GO:0004517                         nitric-oxide synthase activity         2
    ## 18 GO:0005112                                          Notch binding         2
    ## 19 GO:0003723                                            RNA binding       128
    ##    Significant Expected  Fisher               type
    ## 1            4     0.86 0.00210 Biological.Process
    ## 2            6     2.14 0.00890 Biological.Process
    ## 3            4     1.29 0.02150 Biological.Process
    ## 4            4     1.29 0.02150 Biological.Process
    ## 5           19    13.51 0.02930 Biological.Process
    ## 6            4     1.50 0.04190 Biological.Process
    ## 7            2     0.43 0.04580 Biological.Process
    ## 8           41    26.09 0.00094 Molecular.Function
    ## 9          122   106.28 0.00112 Molecular.Function
    ## 10           4     0.86 0.00207 Molecular.Function
    ## 11           3     0.64 0.00974 Molecular.Function
    ## 12           7     2.78 0.01009 Molecular.Function
    ## 13           3     0.86 0.03274 Molecular.Function
    ## 14           3     0.86 0.03274 Molecular.Function
    ## 15           2     0.43 0.04567 Molecular.Function
    ## 16           2     0.43 0.04567 Molecular.Function
    ## 17           2     0.43 0.04567 Molecular.Function
    ## 18           2     0.43 0.04567 Molecular.Function
    ## 19          34    27.37 0.04845 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 261 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  7 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  11 nodes to be scored   (19 eliminated genes)

    ## 
    ##   Level 9:   16 nodes to be scored   (171 eliminated genes)

    ## 
    ##   Level 8:   22 nodes to be scored   (218 eliminated genes)

    ## 
    ##   Level 7:   28 nodes to be scored   (263 eliminated genes)

    ## 
    ##   Level 6:   34 nodes to be scored   (359 eliminated genes)

    ## 
    ##   Level 5:   48 nodes to be scored   (701 eliminated genes)

    ## 
    ##   Level 4:   44 nodes to be scored   (754 eliminated genes)

    ## 
    ##   Level 3:   36 nodes to be scored   (1061 eliminated genes)

    ## 
    ##   Level 2:   11 nodes to be scored   (1263 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1361 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 187 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   6 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   9 nodes to be scored    (148 eliminated genes)

    ## 
    ##   Level 7:   24 nodes to be scored   (441 eliminated genes)

    ## 
    ##   Level 6:   38 nodes to be scored   (463 eliminated genes)

    ## 
    ##   Level 5:   41 nodes to be scored   (857 eliminated genes)

    ## 
    ##   Level 4:   30 nodes to be scored   (1143 eliminated genes)

    ## 
    ##   Level 3:   23 nodes to be scored   (1773 eliminated genes)

    ## 
    ##   Level 2:   11 nodes to be scored   (2039 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2606 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0006886                        intracellular protein transport         4
    ## 2  GO:0001937 negative regulation of endothelial cell proliferati...         1
    ## 3  GO:0000712      resolution of meiotic recombination intermediates         1
    ## 4  GO:0001731         formation of translation preinitiation complex         1
    ## 5  GO:0008218                                        bioluminescence        10
    ## 6  GO:0001965                        G-protein alpha-subunit binding        16
    ## 7  GO:0004505                 phenylalanine 4-monooxygenase activity         1
    ## 8  GO:0004615                            phosphomannomutase activity         1
    ## 9  GO:0004017                              adenylate kinase activity         1
    ## 10 GO:0004149 dihydrolipoyllysine-residue succinyltransferase act...         1
    ## 11 GO:0003874           6-pyruvoyltetrahydropterin synthase activity         1
    ## 12 GO:0004748 ribonucleoside-diphosphate reductase activity, thio...         1
    ## 13 GO:0004108                         citrate (Si)-synthase activity         1
    ## 14 GO:0004067                                  asparaginase activity         1
    ## 15 GO:0004152                  dihydroorotate dehydrogenase activity         1
    ##    Significant Expected Fisher               type
    ## 1            2     0.11 0.0046 Biological.Process
    ## 2            1     0.03 0.0285 Biological.Process
    ## 3            1     0.03 0.0285 Biological.Process
    ## 4            1     0.03 0.0285 Biological.Process
    ## 5            2     0.28 0.0308 Biological.Process
    ## 6            4     0.50 0.0012 Molecular.Function
    ## 7            1     0.03 0.0314 Molecular.Function
    ## 8            1     0.03 0.0314 Molecular.Function
    ## 9            1     0.03 0.0314 Molecular.Function
    ## 10           1     0.03 0.0314 Molecular.Function
    ## 11           1     0.03 0.0314 Molecular.Function
    ## 12           1     0.03 0.0314 Molecular.Function
    ## 13           1     0.03 0.0314 Molecular.Function
    ## 14           1     0.03 0.0314 Molecular.Function
    ## 15           1     0.03 0.0314 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 447 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  6 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  11 nodes to be scored   (5 eliminated genes)

    ## 
    ##   Level 10:  24 nodes to be scored   (75 eliminated genes)

    ## 
    ##   Level 9:   37 nodes to be scored   (260 eliminated genes)

    ## 
    ##   Level 8:   49 nodes to be scored   (364 eliminated genes)

    ## 
    ##   Level 7:   63 nodes to be scored   (449 eliminated genes)

    ## 
    ##   Level 6:   70 nodes to be scored   (609 eliminated genes)

    ## 
    ##   Level 5:   74 nodes to be scored   (892 eliminated genes)

    ## 
    ##   Level 4:   55 nodes to be scored   (1023 eliminated genes)

    ## 
    ##   Level 3:   43 nodes to be scored   (1162 eliminated genes)

    ## 
    ##   Level 2:   12 nodes to be scored   (1283 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1389 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 341 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  6 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   11 nodes to be scored   (1 eliminated genes)

    ## 
    ##   Level 8:   26 nodes to be scored   (155 eliminated genes)

    ## 
    ##   Level 7:   47 nodes to be scored   (445 eliminated genes)

    ## 
    ##   Level 6:   76 nodes to be scored   (649 eliminated genes)

    ## 
    ##   Level 5:   68 nodes to be scored   (1115 eliminated genes)

    ## 
    ##   Level 4:   62 nodes to be scored   (1473 eliminated genes)

    ## 
    ##   Level 3:   29 nodes to be scored   (2085 eliminated genes)

    ## 
    ##   Level 2:   14 nodes to be scored   (2404 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2634 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0019700                  organic phosphonate catabolic process        12
    ## 2 GO:0002064                            epithelial cell development         4
    ## 3 GO:0001707                                     mesoderm formation         4
    ## 4 GO:0003179                              heart valve morphogenesis         4
    ## 5 GO:0008218                                        bioluminescence        10
    ## 6 GO:0003682                                      chromatin binding        13
    ## 7 GO:0001409 guanine nucleotide transmembrane transporter activi...        30
    ## 8 GO:0005245                 voltage-gated calcium channel activity         7
    ##   Significant Expected  Fisher               type
    ## 1           6     0.97 0.00016 Biological.Process
    ## 2           3     0.32 0.00197 Biological.Process
    ## 3           2     0.32 0.03518 Biological.Process
    ## 4           2     0.32 0.03518 Biological.Process
    ## 5           3     0.81 0.04103 Biological.Process
    ## 6           5     1.28 0.00600 Molecular.Function
    ## 7           7     2.96 0.02400 Molecular.Function
    ## 8           3     0.69 0.02500 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 258 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  2 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 11:  7 nodes to be scored    (26 eliminated genes)

    ## 
    ##   Level 10:  9 nodes to be scored    (74 eliminated genes)

    ## 
    ##   Level 9:   14 nodes to be scored   (255 eliminated genes)

    ## 
    ##   Level 8:   22 nodes to be scored   (301 eliminated genes)

    ## 
    ##   Level 7:   35 nodes to be scored   (324 eliminated genes)

    ## 
    ##   Level 6:   42 nodes to be scored   (441 eliminated genes)

    ## 
    ##   Level 5:   49 nodes to be scored   (680 eliminated genes)

    ## 
    ##   Level 4:   35 nodes to be scored   (825 eliminated genes)

    ## 
    ##   Level 3:   29 nodes to be scored   (1062 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (1189 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1359 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 232 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  6 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   11 nodes to be scored   (10 eliminated genes)

    ## 
    ##   Level 8:   17 nodes to be scored   (155 eliminated genes)

    ## 
    ##   Level 7:   35 nodes to be scored   (448 eliminated genes)

    ## 
    ##   Level 6:   41 nodes to be scored   (576 eliminated genes)

    ## 
    ##   Level 5:   45 nodes to be scored   (1040 eliminated genes)

    ## 
    ##   Level 4:   39 nodes to be scored   (1253 eliminated genes)

    ## 
    ##   Level 3:   23 nodes to be scored   (1930 eliminated genes)

    ## 
    ##   Level 2:   12 nodes to be scored   (2230 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2569 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0003341                                        cilium movement        20
    ## 2  GO:0000082                  G1/S transition of mitotic cell cycle       104
    ## 3  GO:0002224                   toll-like receptor signaling pathway         1
    ## 4  GO:0000972 transcription-dependent tethering of RNA polymerase...         1
    ## 5  GO:0006884                                cell volume homeostasis         1
    ## 6  GO:0000340                      RNA 7-methylguanosine cap binding         4
    ## 7  GO:0004252                     serine-type endopeptidase activity        53
    ## 8  GO:0005245                 voltage-gated calcium channel activity         7
    ## 9  GO:0004714 transmembrane receptor protein tyrosine kinase acti...         2
    ## 10 GO:0003810 protein-glutamine gamma-glutamyltransferase activit...         1
    ## 11 GO:0050218                         propionate-CoA ligase activity         1
    ## 12 GO:0042802                              identical protein binding         1
    ## 13 GO:0004368 glycerol-3-phosphate dehydrogenase (quinone) activi...         1
    ## 14 GO:0001758                         retinal dehydrogenase activity         1
    ## 15 GO:0005005                 transmembrane-ephrin receptor activity         1
    ## 16 GO:0004149 dihydrolipoyllysine-residue succinyltransferase act...         1
    ## 17 GO:0003980 UDP-glucose:glycoprotein glucosyltransferase activi...         1
    ## 18 GO:0005381            iron ion transmembrane transporter activity         1
    ## 19 GO:0002161                        aminoacyl-tRNA editing activity         1
    ## 20 GO:0004458          D-lactate dehydrogenase (cytochrome) activity         1
    ## 21 GO:0004485               methylcrotonoyl-CoA carboxylase activity         1
    ## 22 GO:0001002 RNA polymerase III type 1 promoter sequence-specifi...         1
    ## 23 GO:0004334                           fumarylacetoacetase activity         1
    ## 24 GO:0004427             inorganic diphosphate phosphatase activity         1
    ##    Significant Expected Fisher               type
    ## 1            4     0.60 0.0024 Biological.Process
    ## 2            8     3.11 0.0098 Biological.Process
    ## 3            1     0.03 0.0299 Biological.Process
    ## 4            1     0.03 0.0299 Biological.Process
    ## 5            1     0.03 0.0299 Biological.Process
    ## 6            2     0.17 0.0099 Molecular.Function
    ## 7            6     2.22 0.0217 Molecular.Function
    ## 8            2     0.29 0.0317 Molecular.Function
    ## 9            2     0.08 0.0415 Molecular.Function
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
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 2 1 2 1 2 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 281 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  1 nodes to be scored    (11 eliminated genes)

    ## 
    ##   Level 11:  5 nodes to be scored    (15 eliminated genes)

    ## 
    ##   Level 10:  13 nodes to be scored   (63 eliminated genes)

    ## 
    ##   Level 9:   17 nodes to be scored   (231 eliminated genes)

    ## 
    ##   Level 8:   23 nodes to be scored   (295 eliminated genes)

    ## 
    ##   Level 7:   36 nodes to be scored   (342 eliminated genes)

    ## 
    ##   Level 6:   51 nodes to be scored   (486 eliminated genes)

    ## 
    ##   Level 5:   58 nodes to be scored   (748 eliminated genes)

    ## 
    ##   Level 4:   34 nodes to be scored   (884 eliminated genes)

    ## 
    ##   Level 3:   29 nodes to be scored   (1133 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (1268 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1377 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 241 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   10 nodes to be scored   (10 eliminated genes)

    ## 
    ##   Level 8:   23 nodes to be scored   (50 eliminated genes)

    ## 
    ##   Level 7:   35 nodes to be scored   (447 eliminated genes)

    ## 
    ##   Level 6:   43 nodes to be scored   (610 eliminated genes)

    ## 
    ##   Level 5:   47 nodes to be scored   (999 eliminated genes)

    ## 
    ##   Level 4:   45 nodes to be scored   (1279 eliminated genes)

    ## 
    ##   Level 3:   23 nodes to be scored   (2005 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (2288 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2557 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0000902                                     cell morphogenesis         6
    ## 2  GO:0001889                                      liver development         8
    ## 3  GO:0002230 positive regulation of defense response to virus by...         1
    ## 4  GO:0006488 dolichol-linked oligosaccharide biosynthetic proces...         1
    ## 5  GO:0001832                                      blastocyst growth         1
    ## 6  GO:0000076                   DNA replication checkpoint signaling         1
    ## 7  GO:0005245                 voltage-gated calcium channel activity         7
    ## 8  GO:0003964                   RNA-directed DNA polymerase activity       122
    ## 9  GO:0004435          phosphatidylinositol phospholipase C activity         1
    ## 10 GO:0001002 RNA polymerase III type 1 promoter sequence-specifi...         1
    ## 11 GO:0005246                     calcium channel regulator activity         1
    ## 12 GO:0001758                         retinal dehydrogenase activity         1
    ## 13 GO:0005080                               protein kinase C binding         1
    ## 14 GO:0004683 calcium/calmodulin-dependent protein kinase activit...         1
    ## 15 GO:0004104                                cholinesterase activity         1
    ## 16 GO:0004736                          pyruvate carboxylase activity         1
    ## 17 GO:0004494                      methylmalonyl-CoA mutase activity         1
    ## 18 GO:0004844                      uracil DNA N-glycosylase activity         1
    ## 19 GO:0003863 3-methyl-2-oxobutanoate dehydrogenase (2-methylprop...         1
    ## 20 GO:0005381            iron ion transmembrane transporter activity         1
    ## 21 GO:0005324 long-chain fatty acid transmembrane transporter act...         1
    ## 22 GO:0004760                L-serine-pyruvate transaminase activity         1
    ## 23 GO:0000703 oxidized pyrimidine nucleobase lesion DNA N-glycosy...         1
    ## 24 GO:0005524                                            ATP binding       263
    ##    Significant Expected Fisher               type
    ## 1            2     0.22 0.0180 Biological.Process
    ## 2            2     0.30 0.0330 Biological.Process
    ## 3            1     0.04 0.0370 Biological.Process
    ## 4            1     0.04 0.0370 Biological.Process
    ## 5            1     0.04 0.0370 Biological.Process
    ## 6            1     0.04 0.0370 Biological.Process
    ## 7            3     0.30 0.0023 Molecular.Function
    ## 8           12     5.15 0.0046 Molecular.Function
    ## 9            1     0.04 0.0422 Molecular.Function
    ## 10           1     0.04 0.0422 Molecular.Function
    ## 11           1     0.04 0.0422 Molecular.Function
    ## 12           1     0.04 0.0422 Molecular.Function
    ## 13           1     0.04 0.0422 Molecular.Function
    ## 14           1     0.04 0.0422 Molecular.Function
    ## 15           1     0.04 0.0422 Molecular.Function
    ## 16           1     0.04 0.0422 Molecular.Function
    ## 17           1     0.04 0.0422 Molecular.Function
    ## 18           1     0.04 0.0422 Molecular.Function
    ## 19           1     0.04 0.0422 Molecular.Function
    ## 20           1     0.04 0.0422 Molecular.Function
    ## 21           1     0.04 0.0422 Molecular.Function
    ## 22           1     0.04 0.0422 Molecular.Function
    ## 23           1     0.04 0.0422 Molecular.Function
    ## 24          17    11.10 0.0467 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 192 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  6 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  9 nodes to be scored    (2 eliminated genes)

    ## 
    ##   Level 9:   16 nodes to be scored   (156 eliminated genes)

    ## 
    ##   Level 8:   22 nodes to be scored   (194 eliminated genes)

    ## 
    ##   Level 7:   25 nodes to be scored   (224 eliminated genes)

    ## 
    ##   Level 6:   26 nodes to be scored   (312 eliminated genes)

    ## 
    ##   Level 5:   33 nodes to be scored   (423 eliminated genes)

    ## 
    ##   Level 4:   27 nodes to be scored   (565 eliminated genes)

    ## 
    ##   Level 3:   20 nodes to be scored   (967 eliminated genes)

    ## 
    ##   Level 2:   6 nodes to be scored    (1107 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1242 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 135 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   7 nodes to be scored    (107 eliminated genes)

    ## 
    ##   Level 7:   16 nodes to be scored   (432 eliminated genes)

    ## 
    ##   Level 6:   26 nodes to be scored   (533 eliminated genes)

    ## 
    ##   Level 5:   29 nodes to be scored   (796 eliminated genes)

    ## 
    ##   Level 4:   27 nodes to be scored   (1154 eliminated genes)

    ## 
    ##   Level 3:   17 nodes to be scored   (1721 eliminated genes)

    ## 
    ##   Level 2:   6 nodes to be scored    (1982 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2520 eliminated genes)

    ##         GO.ID                                          Term Annotated
    ## 1  GO:0001825                          blastocyst formation         1
    ## 2  GO:0000492                       box C/D snoRNP assembly         1
    ## 3  GO:0000076          DNA replication checkpoint signaling         1
    ## 4  GO:0006527                    arginine catabolic process         1
    ## 5  GO:0001516            prostaglandin biosynthetic process         2
    ## 6  GO:0000050                                    urea cycle         2
    ## 7  GO:0000028              ribosomal small subunit assembly         3
    ## 8  GO:0001708                       cell fate specification         4
    ## 9  GO:0003179                     heart valve morphogenesis         4
    ## 10 GO:0004129                 cytochrome-c oxidase activity         1
    ## 11 GO:0097602                 cullin family protein binding         1
    ## 12 GO:0004516 nicotinate phosphoribosyltransferase activity         1
    ## 13 GO:0003995               acyl-CoA dehydrogenase activity         1
    ## 14 GO:0004408            holocytochrome-c synthase activity         1
    ## 15 GO:0005381   iron ion transmembrane transporter activity         1
    ## 16 GO:0004222                 metalloendopeptidase activity        15
    ## 17 GO:0004658            propionyl-CoA carboxylase activity         2
    ## 18 GO:0003985       acetyl-CoA C-acetyltransferase activity         2
    ## 19 GO:0004474                      malate synthase activity         2
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

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 218 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  5 nodes to be scored    (19 eliminated genes)

    ## 
    ##   Level 10:  8 nodes to be scored    (70 eliminated genes)

    ## 
    ##   Level 9:   14 nodes to be scored   (218 eliminated genes)

    ## 
    ##   Level 8:   20 nodes to be scored   (261 eliminated genes)

    ## 
    ##   Level 7:   29 nodes to be scored   (275 eliminated genes)

    ## 
    ##   Level 6:   36 nodes to be scored   (397 eliminated genes)

    ## 
    ##   Level 5:   42 nodes to be scored   (690 eliminated genes)

    ## 
    ##   Level 4:   28 nodes to be scored   (802 eliminated genes)

    ## 
    ##   Level 3:   21 nodes to be scored   (1014 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (1090 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1206 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 128 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   8 nodes to be scored    (142 eliminated genes)

    ## 
    ##   Level 7:   16 nodes to be scored   (432 eliminated genes)

    ## 
    ##   Level 6:   24 nodes to be scored   (551 eliminated genes)

    ## 
    ##   Level 5:   24 nodes to be scored   (905 eliminated genes)

    ## 
    ##   Level 4:   26 nodes to be scored   (1133 eliminated genes)

    ## 
    ##   Level 3:   15 nodes to be scored   (1663 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (2048 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2462 eliminated genes)

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
    ## 1            2     0.15 0.00930 Biological.Process
    ## 2            1     0.02 0.02210 Biological.Process
    ## 3            1     0.02 0.02210 Biological.Process
    ## 4            1     0.02 0.02210 Biological.Process
    ## 5           10     2.55 0.00016 Molecular.Function
    ## 6            3     0.27 0.00214 Molecular.Function
    ## 7            1     0.04 0.04140 Molecular.Function
    ## 8            1     0.04 0.04140 Molecular.Function
    ## 9            1     0.04 0.04140 Molecular.Function
    ## 10          17    10.40 0.04940 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 306 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  5 nodes to be scored    (11 eliminated genes)

    ## 
    ##   Level 11:  9 nodes to be scored    (16 eliminated genes)

    ## 
    ##   Level 10:  16 nodes to be scored   (83 eliminated genes)

    ## 
    ##   Level 9:   18 nodes to be scored   (259 eliminated genes)

    ## 
    ##   Level 8:   27 nodes to be scored   (324 eliminated genes)

    ## 
    ##   Level 7:   38 nodes to be scored   (384 eliminated genes)

    ## 
    ##   Level 6:   51 nodes to be scored   (547 eliminated genes)

    ## 
    ##   Level 5:   58 nodes to be scored   (753 eliminated genes)

    ## 
    ##   Level 4:   37 nodes to be scored   (896 eliminated genes)

    ## 
    ##   Level 3:   31 nodes to be scored   (1172 eliminated genes)

    ## 
    ##   Level 2:   11 nodes to be scored   (1247 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1355 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 222 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   7 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 8:   15 nodes to be scored   (150 eliminated genes)

    ## 
    ##   Level 7:   34 nodes to be scored   (438 eliminated genes)

    ## 
    ##   Level 6:   44 nodes to be scored   (571 eliminated genes)

    ## 
    ##   Level 5:   41 nodes to be scored   (996 eliminated genes)

    ## 
    ##   Level 4:   41 nodes to be scored   (1313 eliminated genes)

    ## 
    ##   Level 3:   22 nodes to be scored   (1976 eliminated genes)

    ## 
    ##   Level 2:   12 nodes to be scored   (2278 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2594 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0001951                        intestinal D-glucose absorption         1
    ## 2  GO:0000460                                maturation of 5.8S rRNA         1
    ## 3  GO:0002265       astrocyte activation involved in immune response         1
    ## 4  GO:0001825                                   blastocyst formation         1
    ## 5  GO:0000076                   DNA replication checkpoint signaling         1
    ## 6  GO:0006813                                potassium ion transport         1
    ## 7  GO:0016311                                      dephosphorylation         1
    ## 8  GO:0004560                            alpha-L-fucosidase activity         6
    ## 9  GO:0001758                         retinal dehydrogenase activity         1
    ## 10 GO:0004067                                  asparaginase activity         1
    ## 11 GO:0003960                       NADPH:quinone reductase activity         1
    ## 12 GO:0004134                    4-alpha-glucanotransferase activity         1
    ## 13 GO:0003980 UDP-glucose:glycoprotein glucosyltransferase activi...         1
    ## 14 GO:0004325                                ferrochelatase activity         1
    ## 15 GO:0005381            iron ion transmembrane transporter activity         1
    ## 16 GO:0002161                        aminoacyl-tRNA editing activity         1
    ## 17 GO:0003923                       GPI-anchor transamidase activity         1
    ## 18 GO:0001002 RNA polymerase III type 1 promoter sequence-specifi...         1
    ## 19 GO:0004777 succinate-semialdehyde dehydrogenase (NAD+) activit...         1
    ## 20 GO:0004850                         uridine phosphorylase activity         1
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

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 180 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  5 nodes to be scored    (63 eliminated genes)

    ## 
    ##   Level 9:   11 nodes to be scored   (203 eliminated genes)

    ## 
    ##   Level 8:   16 nodes to be scored   (235 eliminated genes)

    ## 
    ##   Level 7:   23 nodes to be scored   (271 eliminated genes)

    ## 
    ##   Level 6:   30 nodes to be scored   (374 eliminated genes)

    ## 
    ##   Level 5:   30 nodes to be scored   (535 eliminated genes)

    ## 
    ##   Level 4:   30 nodes to be scored   (655 eliminated genes)

    ## 
    ##   Level 3:   23 nodes to be scored   (992 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (1164 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1302 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 154 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   4 nodes to be scored    (100 eliminated genes)

    ## 
    ##   Level 7:   16 nodes to be scored   (389 eliminated genes)

    ## 
    ##   Level 6:   34 nodes to be scored   (403 eliminated genes)

    ## 
    ##   Level 5:   38 nodes to be scored   (788 eliminated genes)

    ## 
    ##   Level 4:   34 nodes to be scored   (1048 eliminated genes)

    ## 
    ##   Level 3:   15 nodes to be scored   (1475 eliminated genes)

    ## 
    ##   Level 2:   7 nodes to be scored    (1820 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2289 eliminated genes)

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
    ## 10 GO:0047830                      D-octopine dehydrogenase activity         1
    ## 11 GO:0004066 asparagine synthase (glutamine-hydrolyzing) activit...         1
    ## 12 GO:0004067                                  asparaginase activity         1
    ## 13 GO:0008137               NADH dehydrogenase (ubiquinone) activity         1
    ## 14 GO:0004325                                ferrochelatase activity         1
    ## 15 GO:0003960                       NADPH:quinone reductase activity         1
    ## 16 GO:0004142      diacylglycerol cholinephosphotransferase activity         1
    ## 17 GO:0046982                    protein heterodimerization activity         1
    ## 18 GO:0002161                        aminoacyl-tRNA editing activity         1
    ## 19 GO:0097602                          cullin family protein binding         1
    ## 20 GO:0004334                           fumarylacetoacetase activity         1
    ## 21 GO:0004152                  dihydroorotate dehydrogenase activity         1
    ## 22 GO:0004458          D-lactate dehydrogenase (cytochrome) activity         1
    ## 23 GO:0004174 electron-transferring-flavoprotein dehydrogenase ac...         1
    ## 24 GO:0000049                                           tRNA binding        30
    ## 25 GO:0004089                         carbonate dehydratase activity         2
    ##    Significant Expected    Fisher               type
    ## 1            5     0.26 0.0000023 Biological.Process
    ## 2            2     0.06 0.0013000 Biological.Process
    ## 3            1     0.02 0.0214000 Biological.Process
    ## 4            1     0.02 0.0214000 Biological.Process
    ## 5            4     1.35 0.0420000 Biological.Process
    ## 6            1     0.04 0.0423000 Biological.Process
    ## 7            1     0.04 0.0423000 Biological.Process
    ## 8            3     0.17 0.0004400 Molecular.Function
    ## 9            2     0.05 0.0005800 Molecular.Function
    ## 10           1     0.02 0.0241600 Molecular.Function
    ## 11           1     0.02 0.0241600 Molecular.Function
    ## 12           1     0.02 0.0241600 Molecular.Function
    ## 13           1     0.02 0.0241600 Molecular.Function
    ## 14           1     0.02 0.0241600 Molecular.Function
    ## 15           1     0.02 0.0241600 Molecular.Function
    ## 16           1     0.02 0.0241600 Molecular.Function
    ## 17           1     0.02 0.0241600 Molecular.Function
    ## 18           1     0.02 0.0241600 Molecular.Function
    ## 19           1     0.02 0.0241600 Molecular.Function
    ## 20           1     0.02 0.0241600 Molecular.Function
    ## 21           1     0.02 0.0241600 Molecular.Function
    ## 22           1     0.02 0.0241600 Molecular.Function
    ## 23           1     0.02 0.0241600 Molecular.Function
    ## 24           3     0.72 0.0343700 Molecular.Function
    ## 25           1     0.05 0.0477500 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 176 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  3 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 11:  7 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 10:  10 nodes to be scored   (69 eliminated genes)

    ## 
    ##   Level 9:   13 nodes to be scored   (221 eliminated genes)

    ## 
    ##   Level 8:   18 nodes to be scored   (250 eliminated genes)

    ## 
    ##   Level 7:   22 nodes to be scored   (256 eliminated genes)

    ## 
    ##   Level 6:   22 nodes to be scored   (341 eliminated genes)

    ## 
    ##   Level 5:   30 nodes to be scored   (577 eliminated genes)

    ## 
    ##   Level 4:   22 nodes to be scored   (626 eliminated genes)

    ## 
    ##   Level 3:   18 nodes to be scored   (828 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (1054 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1287 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 98 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   12 nodes to be scored   (264 eliminated genes)

    ## 
    ##   Level 6:   15 nodes to be scored   (366 eliminated genes)

    ## 
    ##   Level 5:   20 nodes to be scored   (709 eliminated genes)

    ## 
    ##   Level 4:   21 nodes to be scored   (926 eliminated genes)

    ## 
    ##   Level 3:   15 nodes to be scored   (1301 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (1940 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2482 eliminated genes)

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
    ## 1     0.05 0.00083 Biological.Process
    ## 2     0.10 0.00403 Biological.Process
    ## 3     1.09 0.01715 Biological.Process
    ## 4     0.72 0.03255 Biological.Process
    ## 5     0.03 0.03391 Biological.Process
    ## 6     1.72 0.00130 Molecular.Function
    ## 7     0.01 0.01410 Molecular.Function
    ## 8     0.24 0.02290 Molecular.Function
    ## 9     0.04 0.04160 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 153 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  3 nodes to be scored    (10 eliminated genes)

    ## 
    ##   Level 9:   6 nodes to be scored    (15 eliminated genes)

    ## 
    ##   Level 8:   10 nodes to be scored   (33 eliminated genes)

    ## 
    ##   Level 7:   20 nodes to be scored   (55 eliminated genes)

    ## 
    ##   Level 6:   27 nodes to be scored   (259 eliminated genes)

    ## 
    ##   Level 5:   29 nodes to be scored   (564 eliminated genes)

    ## 
    ##   Level 4:   25 nodes to be scored   (657 eliminated genes)

    ## 
    ##   Level 3:   20 nodes to be scored   (904 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (1109 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1282 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 130 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   6 nodes to be scored    (100 eliminated genes)

    ## 
    ##   Level 7:   18 nodes to be scored   (383 eliminated genes)

    ## 
    ##   Level 6:   24 nodes to be scored   (518 eliminated genes)

    ## 
    ##   Level 5:   28 nodes to be scored   (993 eliminated genes)

    ## 
    ##   Level 4:   24 nodes to be scored   (1231 eliminated genes)

    ## 
    ##   Level 3:   19 nodes to be scored   (1776 eliminated genes)

    ## 
    ##   Level 2:   6 nodes to be scored    (2055 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2475 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0000724 double-strand break repair via homologous recombina...        13
    ## 2  GO:0000184 nuclear-transcribed mRNA catabolic process, nonsens...        10
    ## 3  GO:0000712      resolution of meiotic recombination intermediates         1
    ## 4  GO:0006813                                potassium ion transport         1
    ## 5  GO:0000423                                              mitophagy         2
    ## 6  GO:0002933                                    lipid hydroxylation         3
    ## 7  GO:0000028                       ribosomal small subunit assembly         3
    ## 8  GO:0005524                                            ATP binding       263
    ## 9  GO:0004013                        adenosylhomocysteinase activity         1
    ## 10 GO:0004108                         citrate (Si)-synthase activity         1
    ## 11 GO:0004089                         carbonate dehydratase activity         2
    ## 12 GO:0004517                         nitric-oxide synthase activity         2
    ## 13 GO:0004185                  serine-type carboxypeptidase activity         2
    ## 14 GO:0005201            extracellular matrix structural constituent        17
    ##    Significant Expected  Fisher               type
    ## 1            3     0.20 0.00086 Biological.Process
    ## 2            2     0.16 0.00978 Biological.Process
    ## 3            1     0.02 0.01567 Biological.Process
    ## 4            1     0.02 0.01567 Biological.Process
    ## 5            1     0.03 0.03110 Biological.Process
    ## 6            1     0.05 0.04631 Biological.Process
    ## 7            1     0.05 0.04631 Biological.Process
    ## 8           12     5.03 0.00320 Molecular.Function
    ## 9            1     0.02 0.01910 Molecular.Function
    ## 10           1     0.02 0.01910 Molecular.Function
    ## 11           1     0.04 0.03790 Molecular.Function
    ## 12           1     0.04 0.03790 Molecular.Function
    ## 13           1     0.04 0.03790 Molecular.Function
    ## 14           2     0.32 0.04060 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 387 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  5 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 12:  5 nodes to be scored    (7 eliminated genes)

    ## 
    ##   Level 11:  8 nodes to be scored    (38 eliminated genes)

    ## 
    ##   Level 10:  13 nodes to be scored   (84 eliminated genes)

    ## 
    ##   Level 9:   26 nodes to be scored   (257 eliminated genes)

    ## 
    ##   Level 8:   45 nodes to be scored   (346 eliminated genes)

    ## 
    ##   Level 7:   55 nodes to be scored   (431 eliminated genes)

    ## 
    ##   Level 6:   62 nodes to be scored   (610 eliminated genes)

    ## 
    ##   Level 5:   67 nodes to be scored   (897 eliminated genes)

    ## 
    ##   Level 4:   48 nodes to be scored   (958 eliminated genes)

    ## 
    ##   Level 3:   39 nodes to be scored   (1181 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (1317 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1381 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 220 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   14 nodes to be scored   (142 eliminated genes)

    ## 
    ##   Level 7:   31 nodes to be scored   (432 eliminated genes)

    ## 
    ##   Level 6:   45 nodes to be scored   (577 eliminated genes)

    ## 
    ##   Level 5:   48 nodes to be scored   (1072 eliminated genes)

    ## 
    ##   Level 4:   44 nodes to be scored   (1380 eliminated genes)

    ## 
    ##   Level 3:   22 nodes to be scored   (2001 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (2304 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2593 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0002291 T cell activation via T cell receptor contact with ...        10
    ## 2 GO:0002674     negative regulation of acute inflammatory response         4
    ## 3 GO:0002221         pattern recognition receptor signaling pathway        11
    ## 4 GO:0005524                                            ATP binding       263
    ## 5 GO:0003676                                   nucleic acid binding       497
    ## 6 GO:0003824                                     catalytic activity      1038
    ## 7 GO:0004622                             lysophospholipase activity         4
    ## 8 GO:0003964                   RNA-directed DNA polymerase activity       122
    ##   Significant Expected  Fisher               type
    ## 1           3     0.53 0.01300 Biological.Process
    ## 2           2     0.21 0.01500 Biological.Process
    ## 3           3     0.58 0.01700 Biological.Process
    ## 4          29    15.27 0.00038 Molecular.Function
    ## 5          34    28.86 0.00102 Molecular.Function
    ## 6          64    60.27 0.01099 Molecular.Function
    ## 7           2     0.23 0.01860 Molecular.Function
    ## 8          12     7.08 0.04776 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 174 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  1 nodes to be scored    (11 eliminated genes)

    ## 
    ##   Level 11:  2 nodes to be scored    (19 eliminated genes)

    ## 
    ##   Level 10:  7 nodes to be scored    (63 eliminated genes)

    ## 
    ##   Level 9:   8 nodes to be scored    (203 eliminated genes)

    ## 
    ##   Level 8:   17 nodes to be scored   (245 eliminated genes)

    ## 
    ##   Level 7:   20 nodes to be scored   (265 eliminated genes)

    ## 
    ##   Level 6:   25 nodes to be scored   (414 eliminated genes)

    ## 
    ##   Level 5:   32 nodes to be scored   (512 eliminated genes)

    ## 
    ##   Level 4:   27 nodes to be scored   (580 eliminated genes)

    ## 
    ##   Level 3:   21 nodes to be scored   (867 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (1040 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1222 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 147 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   7 nodes to be scored    (7 eliminated genes)

    ## 
    ##   Level 7:   22 nodes to be scored   (325 eliminated genes)

    ## 
    ##   Level 6:   30 nodes to be scored   (353 eliminated genes)

    ## 
    ##   Level 5:   29 nodes to be scored   (673 eliminated genes)

    ## 
    ##   Level 4:   27 nodes to be scored   (1079 eliminated genes)

    ## 
    ##   Level 3:   20 nodes to be scored   (1667 eliminated genes)

    ## 
    ##   Level 2:   7 nodes to be scored    (1900 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2271 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0006598                            polyamine catabolic process         3
    ## 2  GO:0006813                                potassium ion transport         1
    ## 3  GO:0042060                                          wound healing         1
    ## 4  GO:0001867                  complement activation, lectin pathway        13
    ## 5  GO:0001578                           microtubule bundle formation         2
    ## 6  GO:0003985                acetyl-CoA C-acetyltransferase activity         2
    ## 7  GO:0003676                                   nucleic acid binding       497
    ## 8  GO:0003997                              acyl-CoA oxidase activity         1
    ## 9  GO:0004060                 arylamine N-acetyltransferase activity         1
    ## 10 GO:0003944 N-acetylglucosamine-1-phosphodiester alpha-N-acetyl...         1
    ## 11 GO:0004748 ribonucleoside-diphosphate reductase activity, thio...         1
    ## 12 GO:0003923                       GPI-anchor transamidase activity         1
    ## 13 GO:0004777 succinate-semialdehyde dehydrogenase (NAD+) activit...         1
    ## 14 GO:0003913                                DNA photolyase activity         1
    ## 15 GO:0004176                       ATP-dependent peptidase activity         2
    ##    Significant Expected  Fisher               type
    ## 1            2     0.05 0.00090 Biological.Process
    ## 2            1     0.02 0.01780 Biological.Process
    ## 3            1     0.02 0.01780 Biological.Process
    ## 4            2     0.23 0.02110 Biological.Process
    ## 5            1     0.04 0.03530 Biological.Process
    ## 6            2     0.05 0.00058 Molecular.Function
    ## 7           19    12.01 0.00266 Molecular.Function
    ## 8            1     0.02 0.02416 Molecular.Function
    ## 9            1     0.02 0.02416 Molecular.Function
    ## 10           1     0.02 0.02416 Molecular.Function
    ## 11           1     0.02 0.02416 Molecular.Function
    ## 12           1     0.02 0.02416 Molecular.Function
    ## 13           1     0.02 0.02416 Molecular.Function
    ## 14           1     0.02 0.02416 Molecular.Function
    ## 15           1     0.05 0.04775 Molecular.Function

``` r
head(results_all_targets)
```

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0000082                  G1/S transition of mitotic cell cycle       104
    ## 2 GO:0001913                           T cell mediated cytotoxicity         1
    ## 3 GO:0002230 positive regulation of defense response to virus by...         1
    ## 4 GO:0001172                            RNA-templated transcription         1
    ## 5 GO:0004176                       ATP-dependent peptidase activity         2
    ## 6 GO:0003964                   RNA-directed DNA polymerase activity       122
    ##   Significant Expected  Fisher               type         miRNA
    ## 1           9     2.67 0.00083 Biological.Process Cluster_10452
    ## 2           1     0.03 0.02564 Biological.Process Cluster_10452
    ## 3           1     0.03 0.02564 Biological.Process Cluster_10452
    ## 4           1     0.03 0.02564 Biological.Process Cluster_10452
    ## 5           2     0.06 0.00100 Molecular.Function Cluster_10452
    ## 6          11     3.87 0.00140 Molecular.Function Cluster_10452

Save results

``` r
write.csv(results_all_targets, "../output/27-Apul-mRNA-miRNA-interactions-topGO/miRNA_all_targets_topGO_FE.csv")
```

# 4 FE of specific miRNA’s targets (high cor targets)

Loop through all miRNA and run functional enrichment on the miRNA’s
highly correlated targets

``` r
interacting_miRNAs_high <- unique(high_cor_bind_FA$miRNA) %>% na.omit
results_high_cor_targets <- NULL  # initialize empty df

for(miRNA in interacting_miRNAs_high) {
  
  # Run topGO enrichment function
  miRNA_results <- miRNA_topGO_FE(miRNA, high_cor_bind_FA)
  
  # Only keep results if not empty
  if (nrow(miRNA_results) > 0) {
    
    # Add the miRNA source column
    miRNA_results$miRNA <- miRNA

    # Bind to the accumulating results data frame
    results_high_cor_targets <- rbind(results_high_cor_targets, miRNA_results)
  }
}
```

    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 25 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   1 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 7:   1 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 6:   3 nodes to be scored    (11 eliminated genes)

    ## 
    ##   Level 5:   5 nodes to be scored    (11 eliminated genes)

    ## 
    ##   Level 4:   5 nodes to be scored    (14 eliminated genes)

    ## 
    ##   Level 3:   3 nodes to be scored    (191 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (233 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (312 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 37 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   6 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 5:   6 nodes to be scored    (169 eliminated genes)

    ## 
    ##   Level 4:   10 nodes to be scored   (196 eliminated genes)

    ## 
    ##   Level 3:   6 nodes to be scored    (369 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (952 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1775 eliminated genes)

    ##        GO.ID                                               Term Annotated
    ## 1 GO:0006598                        polyamine catabolic process         3
    ## 2 GO:0001745                         compound eye morphogenesis         3
    ## 3 GO:0000014 single-stranded DNA endodeoxyribonuclease activity         3
    ## 4 GO:0005290     L-histidine transmembrane transporter activity         9
    ##   Significant Expected Fisher               type
    ## 1           1     0.00 0.0043 Biological.Process
    ## 2           1     0.00 0.0043 Biological.Process
    ## 3           1     0.01 0.0065 Molecular.Function
    ## 4           1     0.02 0.0193 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 41 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  1 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 11:  2 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 10:  3 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 9:   3 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 8:   7 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 7:   5 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 6:   4 nodes to be scored    (11 eliminated genes)

    ## 
    ##   Level 5:   6 nodes to be scored    (125 eliminated genes)

    ## 
    ##   Level 4:   4 nodes to be scored    (138 eliminated genes)

    ## 
    ##   Level 3:   1 nodes to be scored    (203 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (232 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (245 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

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
    ##   Level 7:   6 nodes to be scored    (371 eliminated genes)

    ## 
    ##   Level 6:   10 nodes to be scored   (406 eliminated genes)

    ## 
    ##   Level 5:   11 nodes to be scored   (561 eliminated genes)

    ## 
    ##   Level 4:   11 nodes to be scored   (652 eliminated genes)

    ## 
    ##   Level 3:   9 nodes to be scored    (1160 eliminated genes)

    ## 
    ##   Level 2:   4 nodes to be scored    (1698 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2034 eliminated genes)

    ##        GO.ID                                             Term Annotated
    ## 1 GO:0001658 branching involved in ureteric bud morphogenesis         4
    ## 2 GO:0001945                         lymph vessel development         5
    ## 3 GO:0004252               serine-type endopeptidase activity        53
    ## 4 GO:0005201      extracellular matrix structural constituent        17
    ## 5 GO:0004771                   sterol ester esterase activity         5
    ##   Significant Expected Fisher               type
    ## 1           1     0.01 0.0057 Biological.Process
    ## 2           1     0.01 0.0071 Biological.Process
    ## 3           3     0.25 0.0016 Molecular.Function
    ## 4           2     0.08 0.0027 Molecular.Function
    ## 5           1     0.02 0.0232 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 70 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   2 nodes to be scored    (38 eliminated genes)

    ## 
    ##   Level 8:   3 nodes to be scored    (63 eliminated genes)

    ## 
    ##   Level 7:   6 nodes to be scored    (66 eliminated genes)

    ## 
    ##   Level 6:   10 nodes to be scored   (116 eliminated genes)

    ## 
    ##   Level 5:   13 nodes to be scored   (188 eliminated genes)

    ## 
    ##   Level 4:   12 nodes to be scored   (222 eliminated genes)

    ## 
    ##   Level 3:   11 nodes to be scored   (331 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (608 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (748 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

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
    ##   Level 7:   5 nodes to be scored    (263 eliminated genes)

    ## 
    ##   Level 6:   5 nodes to be scored    (264 eliminated genes)

    ## 
    ##   Level 5:   6 nodes to be scored    (514 eliminated genes)

    ## 
    ##   Level 4:   8 nodes to be scored    (575 eliminated genes)

    ## 
    ##   Level 3:   7 nodes to be scored    (621 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (1262 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1836 eliminated genes)

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
    ## 4     0.01 0.0094 Molecular.Function
    ## 5     0.57 0.0172 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 73 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   2 nodes to be scored    (38 eliminated genes)

    ## 
    ##   Level 8:   3 nodes to be scored    (63 eliminated genes)

    ## 
    ##   Level 7:   6 nodes to be scored    (66 eliminated genes)

    ## 
    ##   Level 6:   10 nodes to be scored   (116 eliminated genes)

    ## 
    ##   Level 5:   14 nodes to be scored   (188 eliminated genes)

    ## 
    ##   Level 4:   13 nodes to be scored   (222 eliminated genes)

    ## 
    ##   Level 3:   12 nodes to be scored   (351 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (627 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (810 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

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
    ##   Level 7:   5 nodes to be scored    (263 eliminated genes)

    ## 
    ##   Level 6:   6 nodes to be scored    (264 eliminated genes)

    ## 
    ##   Level 5:   6 nodes to be scored    (514 eliminated genes)

    ## 
    ##   Level 4:   7 nodes to be scored    (578 eliminated genes)

    ## 
    ##   Level 3:   7 nodes to be scored    (623 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (767 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1591 eliminated genes)

    ##        GO.ID                                 Term Annotated Significant
    ## 1 GO:0001895                   retina homeostasis         2           1
    ## 2 GO:0005524                          ATP binding       263           6
    ## 3 GO:0003964 RNA-directed DNA polymerase activity       122           4
    ## 4 GO:0009034               tryptophanase activity         2           1
    ## 5 GO:0004322                 ferroxidase activity         3           1
    ## 6 GO:0000287                magnesium ion binding        56           2
    ##   Expected Fisher               type
    ## 1     0.01 0.0100 Biological.Process
    ## 2     1.42 0.0016 Molecular.Function
    ## 3     0.66 0.0033 Molecular.Function
    ## 4     0.01 0.0108 Molecular.Function
    ## 5     0.02 0.0161 Molecular.Function
    ## 6     0.30 0.0355 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 70 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   2 nodes to be scored    (38 eliminated genes)

    ## 
    ##   Level 8:   3 nodes to be scored    (63 eliminated genes)

    ## 
    ##   Level 7:   6 nodes to be scored    (66 eliminated genes)

    ## 
    ##   Level 6:   10 nodes to be scored   (116 eliminated genes)

    ## 
    ##   Level 5:   13 nodes to be scored   (188 eliminated genes)

    ## 
    ##   Level 4:   12 nodes to be scored   (222 eliminated genes)

    ## 
    ##   Level 3:   11 nodes to be scored   (331 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (608 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (748 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

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
    ##   Level 7:   4 nodes to be scored    (263 eliminated genes)

    ## 
    ##   Level 6:   6 nodes to be scored    (264 eliminated genes)

    ## 
    ##   Level 5:   6 nodes to be scored    (448 eliminated genes)

    ## 
    ##   Level 4:   8 nodes to be scored    (578 eliminated genes)

    ## 
    ##   Level 3:   7 nodes to be scored    (623 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (1264 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1591 eliminated genes)

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
    ## 4     0.53 0.0137 Molecular.Function
    ## 5     1.14 0.0212 Molecular.Function
    ## 6     0.24 0.0232 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 23 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   1 nodes to be scored    (12 eliminated genes)

    ## 
    ##   Level 8:   2 nodes to be scored    (37 eliminated genes)

    ## 
    ##   Level 7:   3 nodes to be scored    (40 eliminated genes)

    ## 
    ##   Level 6:   4 nodes to be scored    (57 eliminated genes)

    ## 
    ##   Level 5:   4 nodes to be scored    (137 eliminated genes)

    ## 
    ##   Level 4:   3 nodes to be scored    (148 eliminated genes)

    ## 
    ##   Level 3:   1 nodes to be scored    (211 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (228 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (245 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

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
    ##   Level 7:   5 nodes to be scored    (263 eliminated genes)

    ## 
    ##   Level 6:   6 nodes to be scored    (267 eliminated genes)

    ## 
    ##   Level 5:   8 nodes to be scored    (451 eliminated genes)

    ## 
    ##   Level 4:   8 nodes to be scored    (654 eliminated genes)

    ## 
    ##   Level 3:   7 nodes to be scored    (1192 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (1464 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1789 eliminated genes)

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

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 13 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   1 nodes to be scored    (7 eliminated genes)

    ## 
    ##   Level 5:   3 nodes to be scored    (7 eliminated genes)

    ## 
    ##   Level 4:   2 nodes to be scored    (7 eliminated genes)

    ## 
    ##   Level 3:   3 nodes to be scored    (11 eliminated genes)

    ## 
    ##   Level 2:   1 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (499 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

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

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 95 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   6 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 8:   8 nodes to be scored    (152 eliminated genes)

    ## 
    ##   Level 7:   8 nodes to be scored    (157 eliminated genes)

    ## 
    ##   Level 6:   12 nodes to be scored   (221 eliminated genes)

    ## 
    ##   Level 5:   17 nodes to be scored   (394 eliminated genes)

    ## 
    ##   Level 4:   17 nodes to be scored   (434 eliminated genes)

    ## 
    ##   Level 3:   16 nodes to be scored   (804 eliminated genes)

    ## 
    ##   Level 2:   6 nodes to be scored    (941 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1251 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 92 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   8 nodes to be scored    (263 eliminated genes)

    ## 
    ##   Level 6:   16 nodes to be scored   (294 eliminated genes)

    ## 
    ##   Level 5:   19 nodes to be scored   (620 eliminated genes)

    ## 
    ##   Level 4:   21 nodes to be scored   (836 eliminated genes)

    ## 
    ##   Level 3:   16 nodes to be scored   (1364 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (1698 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2472 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0001676                long-chain fatty acid metabolic process         3
    ## 2  GO:0006598                            polyamine catabolic process         3
    ## 3  GO:0001822                                     kidney development        65
    ## 4  GO:0003777                             microtubule motor activity        16
    ## 5  GO:0005080                               protein kinase C binding         1
    ## 6  GO:0004591 oxoglutarate dehydrogenase (succinyl-transferring) ...         1
    ## 7  GO:0004494                      methylmalonyl-CoA mutase activity         1
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
    ## 8            1     0.02 0.0180 Molecular.Function
    ## 9            1     0.03 0.0268 Molecular.Function
    ## 10           1     0.04 0.0356 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 97 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  1 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 11:  2 nodes to be scored    (2 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (63 eliminated genes)

    ## 
    ##   Level 9:   6 nodes to be scored    (203 eliminated genes)

    ## 
    ##   Level 8:   8 nodes to be scored    (205 eliminated genes)

    ## 
    ##   Level 7:   9 nodes to be scored    (209 eliminated genes)

    ## 
    ##   Level 6:   13 nodes to be scored   (280 eliminated genes)

    ## 
    ##   Level 5:   18 nodes to be scored   (304 eliminated genes)

    ## 
    ##   Level 4:   13 nodes to be scored   (345 eliminated genes)

    ## 
    ##   Level 3:   13 nodes to be scored   (496 eliminated genes)

    ## 
    ##   Level 2:   7 nodes to be scored    (861 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1093 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 117 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   5 nodes to be scored    (100 eliminated genes)

    ## 
    ##   Level 7:   18 nodes to be scored   (396 eliminated genes)

    ## 
    ##   Level 6:   24 nodes to be scored   (426 eliminated genes)

    ## 
    ##   Level 5:   25 nodes to be scored   (583 eliminated genes)

    ## 
    ##   Level 4:   23 nodes to be scored   (823 eliminated genes)

    ## 
    ##   Level 3:   10 nodes to be scored   (1273 eliminated genes)

    ## 
    ##   Level 2:   6 nodes to be scored    (1756 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2141 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0000349 generation of catalytic spliceosome for first trans...         1
    ## 2  GO:0001825                                   blastocyst formation         1
    ## 3  GO:0002790                                      peptide secretion         2
    ## 4  GO:0000050                                             urea cycle         2
    ## 5  GO:0001734                  mRNA m(6)A methyltransferase activity         1
    ## 6  GO:0004325                                ferrochelatase activity         1
    ## 7  GO:0004777 succinate-semialdehyde dehydrogenase (NAD+) activit...         1
    ## 8  GO:0004408                     holocytochrome-c synthase activity         1
    ## 9  GO:0008137               NADH dehydrogenase (ubiquinone) activity         1
    ## 10 GO:0004035                          alkaline phosphatase activity         1
    ## 11 GO:0004067                                  asparaginase activity         1
    ## 12 GO:0004748 ribonucleoside-diphosphate reductase activity, thio...         1
    ## 13 GO:0003730                                    mRNA 3'-UTR binding         2
    ## 14 GO:0004176                       ATP-dependent peptidase activity         2
    ## 15 GO:0003824                                     catalytic activity      1038
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
    ## 15          17    10.48 0.0280 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 301 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  6 nodes to be scored    (27 eliminated genes)

    ## 
    ##   Level 10:  9 nodes to be scored    (85 eliminated genes)

    ## 
    ##   Level 9:   19 nodes to be scored   (247 eliminated genes)

    ## 
    ##   Level 8:   25 nodes to be scored   (278 eliminated genes)

    ## 
    ##   Level 7:   39 nodes to be scored   (332 eliminated genes)

    ## 
    ##   Level 6:   51 nodes to be scored   (456 eliminated genes)

    ## 
    ##   Level 5:   53 nodes to be scored   (743 eliminated genes)

    ## 
    ##   Level 4:   45 nodes to be scored   (839 eliminated genes)

    ## 
    ##   Level 3:   36 nodes to be scored   (1069 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (1253 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1368 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 177 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   5 nodes to be scored    (9 eliminated genes)

    ## 
    ##   Level 8:   10 nodes to be scored   (42 eliminated genes)

    ## 
    ##   Level 7:   24 nodes to be scored   (429 eliminated genes)

    ## 
    ##   Level 6:   34 nodes to be scored   (547 eliminated genes)

    ## 
    ##   Level 5:   38 nodes to be scored   (915 eliminated genes)

    ## 
    ##   Level 4:   32 nodes to be scored   (1266 eliminated genes)

    ## 
    ##   Level 3:   22 nodes to be scored   (1935 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (2207 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2541 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0002674     negative regulation of acute inflammatory response         4
    ## 2  GO:0001960 negative regulation of cytokine-mediated signaling ...         1
    ## 3  GO:0001967                                      suckling behavior         1
    ## 4  GO:0000183                         rDNA heterochromatin formation         1
    ## 5  GO:0001409 guanine nucleotide transmembrane transporter activi...        30
    ## 6  GO:0003682                                      chromatin binding        13
    ## 7  GO:0003723                                            RNA binding       128
    ## 8  GO:0004568                                     chitinase activity         1
    ## 9  GO:0004143           ATP-dependent diacylglycerol kinase activity         1
    ## 10 GO:0005289 high-affinity L-arginine transmembrane transporter ...         1
    ## 11 GO:0004033                   aldo-keto reductase (NADPH) activity         1
    ## 12 GO:0004441       inositol-1,4-bisphosphate 1-phosphatase activity         1
    ## 13 GO:0003676                                   nucleic acid binding       497
    ##    Significant Expected Fisher               type
    ## 1            2     0.15 0.0077 Biological.Process
    ## 2            1     0.04 0.0370 Biological.Process
    ## 3            1     0.04 0.0370 Biological.Process
    ## 4            1     0.04 0.0370 Biological.Process
    ## 5            5     1.07 0.0037 Molecular.Function
    ## 6            3     0.46 0.0097 Molecular.Function
    ## 7           10     4.57 0.0145 Molecular.Function
    ## 8            1     0.04 0.0357 Molecular.Function
    ## 9            1     0.04 0.0357 Molecular.Function
    ## 10           1     0.04 0.0357 Molecular.Function
    ## 11           1     0.04 0.0357 Molecular.Function
    ## 12           1     0.04 0.0357 Molecular.Function
    ## 13          25    17.74 0.0391 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 104 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   4 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 8:   8 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 7:   10 nodes to be scored   (144 eliminated genes)

    ## 
    ##   Level 6:   17 nodes to be scored   (239 eliminated genes)

    ## 
    ##   Level 5:   20 nodes to be scored   (427 eliminated genes)

    ## 
    ##   Level 4:   17 nodes to be scored   (479 eliminated genes)

    ## 
    ##   Level 3:   17 nodes to be scored   (646 eliminated genes)

    ## 
    ##   Level 2:   7 nodes to be scored    (871 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1205 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 53 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   6 nodes to be scored    (263 eliminated genes)

    ## 
    ##   Level 6:   8 nodes to be scored    (265 eliminated genes)

    ## 
    ##   Level 5:   9 nodes to be scored    (432 eliminated genes)

    ## 
    ##   Level 4:   14 nodes to be scored   (463 eliminated genes)

    ## 
    ##   Level 3:   8 nodes to be scored    (749 eliminated genes)

    ## 
    ##   Level 2:   4 nodes to be scored    (1332 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1714 eliminated genes)

    ##        GO.ID                        Term Annotated Significant Expected Fisher
    ## 1 GO:0001755 neural crest cell migration         2           1     0.01 0.0110
    ## 2 GO:0000165                MAPK cascade        35           2     0.20 0.0150
    ## 3 GO:0004017   adenylate kinase activity         1           1     0.01 0.0061
    ## 4 GO:0005524                 ATP binding       263           5     1.61 0.0176
    ##                 type
    ## 1 Biological.Process
    ## 2 Biological.Process
    ## 3 Molecular.Function
    ## 4 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 25 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   4 nodes to be scored    (35 eliminated genes)

    ## 
    ##   Level 5:   4 nodes to be scored    (154 eliminated genes)

    ## 
    ##   Level 4:   4 nodes to be scored    (170 eliminated genes)

    ## 
    ##   Level 3:   5 nodes to be scored    (213 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (442 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (511 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 45 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   5 nodes to be scored    (265 eliminated genes)

    ## 
    ##   Level 6:   7 nodes to be scored    (266 eliminated genes)

    ## 
    ##   Level 5:   8 nodes to be scored    (456 eliminated genes)

    ## 
    ##   Level 4:   7 nodes to be scored    (493 eliminated genes)

    ## 
    ##   Level 3:   8 nodes to be scored    (754 eliminated genes)

    ## 
    ##   Level 2:   4 nodes to be scored    (1458 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1813 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0000165                                           MAPK cascade        35
    ## 2 GO:0005524                                            ATP binding       263
    ## 3 GO:0004698 calcium,diacylglycerol-dependent serine/threonine k...         2
    ## 4 GO:0000900 mRNA regulatory element binding translation repress...         3
    ## 5 GO:0001640 adenylate cyclase inhibiting G protein-coupled glut...         9
    ##   Significant Expected   Fisher               type
    ## 1           4     0.15 4.70e-06 Biological.Process
    ## 2           5     0.95 1.20e-03 Molecular.Function
    ## 3           1     0.01 7.20e-03 Molecular.Function
    ## 4           1     0.01 1.08e-02 Molecular.Function
    ## 5           1     0.03 3.20e-02 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 90 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  1 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 11:  3 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 9:   5 nodes to be scored    (16 eliminated genes)

    ## 
    ##   Level 8:   10 nodes to be scored   (16 eliminated genes)

    ## 
    ##   Level 7:   12 nodes to be scored   (17 eliminated genes)

    ## 
    ##   Level 6:   12 nodes to be scored   (26 eliminated genes)

    ## 
    ##   Level 5:   16 nodes to be scored   (278 eliminated genes)

    ## 
    ##   Level 4:   11 nodes to be scored   (301 eliminated genes)

    ## 
    ##   Level 3:   8 nodes to be scored    (424 eliminated genes)

    ## 
    ##   Level 2:   5 nodes to be scored    (619 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (892 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 47 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 7:   4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   7 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 5:   11 nodes to be scored   (145 eliminated genes)

    ## 
    ##   Level 4:   12 nodes to be scored   (256 eliminated genes)

    ## 
    ##   Level 3:   8 nodes to be scored    (898 eliminated genes)

    ## 
    ##   Level 2:   4 nodes to be scored    (1248 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2094 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0001945                               lymph vessel development         5
    ## 2 GO:0003333                     amino acid transmembrane transport         2
    ## 3 GO:0001658       branching involved in ureteric bud morphogenesis         4
    ## 4 GO:0001889                                      liver development         8
    ## 5 GO:0004252                     serine-type endopeptidase activity        53
    ## 6 GO:0046982                    protein heterodimerization activity         1
    ## 7 GO:0004771                         sterol ester esterase activity         5
    ## 8 GO:0005302          L-tyrosine transmembrane transporter activity         9
    ## 9 GO:0001640 adenylate cyclase inhibiting G protein-coupled glut...         9
    ##   Significant Expected  Fisher               type
    ## 1           2     0.03 0.00028 Biological.Process
    ## 2           1     0.01 0.01137 Biological.Process
    ## 3           1     0.02 0.02262 Biological.Process
    ## 4           1     0.05 0.04479 Biological.Process
    ## 5           3     0.29 0.00260 Molecular.Function
    ## 6           1     0.01 0.00540 Molecular.Function
    ## 7           1     0.03 0.02680 Molecular.Function
    ## 8           1     0.05 0.04770 Molecular.Function
    ## 9           1     0.05 0.04770 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 35 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   4 nodes to be scored    (9 eliminated genes)

    ## 
    ##   Level 7:   3 nodes to be scored    (9 eliminated genes)

    ## 
    ##   Level 6:   4 nodes to be scored    (10 eliminated genes)

    ## 
    ##   Level 5:   5 nodes to be scored    (22 eliminated genes)

    ## 
    ##   Level 4:   4 nodes to be scored    (37 eliminated genes)

    ## 
    ##   Level 3:   8 nodes to be scored    (78 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (345 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (609 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 40 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   4 nodes to be scored    (263 eliminated genes)

    ## 
    ##   Level 6:   5 nodes to be scored    (264 eliminated genes)

    ## 
    ##   Level 5:   6 nodes to be scored    (427 eliminated genes)

    ## 
    ##   Level 4:   10 nodes to be scored   (429 eliminated genes)

    ## 
    ##   Level 3:   9 nodes to be scored    (497 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (1364 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1735 eliminated genes)

    ##        GO.ID                                       Term Annotated Significant
    ## 1 GO:0015969 guanosine tetraphosphate metabolic process         9           1
    ## 2 GO:0019700      organic phosphonate catabolic process        12           1
    ## 3 GO:0001666                        response to hypoxia        14           1
    ## 4 GO:0005178                           integrin binding         4           1
    ## 5 GO:0003964       RNA-directed DNA polymerase activity       122           2
    ##   Expected Fisher               type
    ## 1     0.02  0.019 Biological.Process
    ## 2     0.03  0.025 Biological.Process
    ## 3     0.03  0.030 Biological.Process
    ## 4     0.01  0.011 Molecular.Function
    ## 5     0.35  0.045 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 147 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  2 nodes to be scored    (11 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (11 eliminated genes)

    ## 
    ##   Level 9:   9 nodes to be scored    (166 eliminated genes)

    ## 
    ##   Level 8:   13 nodes to be scored   (184 eliminated genes)

    ## 
    ##   Level 7:   17 nodes to be scored   (229 eliminated genes)

    ## 
    ##   Level 6:   23 nodes to be scored   (297 eliminated genes)

    ## 
    ##   Level 5:   27 nodes to be scored   (544 eliminated genes)

    ## 
    ##   Level 4:   21 nodes to be scored   (600 eliminated genes)

    ## 
    ##   Level 3:   20 nodes to be scored   (879 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (986 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1237 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 103 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   8 nodes to be scored    (142 eliminated genes)

    ## 
    ##   Level 7:   15 nodes to be scored   (432 eliminated genes)

    ## 
    ##   Level 6:   18 nodes to be scored   (541 eliminated genes)

    ## 
    ##   Level 5:   20 nodes to be scored   (793 eliminated genes)

    ## 
    ##   Level 4:   18 nodes to be scored   (1082 eliminated genes)

    ## 
    ##   Level 3:   12 nodes to be scored   (1582 eliminated genes)

    ## 
    ##   Level 2:   5 nodes to be scored    (1983 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2303 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0002674     negative regulation of acute inflammatory response         4
    ## 2 GO:0002221         pattern recognition receptor signaling pathway        11
    ## 3 GO:0000165                                           MAPK cascade        35
    ## 4 GO:0001818             negative regulation of cytokine production        18
    ## 5 GO:0001640 adenylate cyclase inhibiting G protein-coupled glut...         9
    ## 6 GO:0004143           ATP-dependent diacylglycerol kinase activity         1
    ## 7 GO:0005524                                            ATP binding       263
    ## 8 GO:0000026                 alpha-1,2-mannosyltransferase activity         3
    ##   Significant Expected   Fisher               type
    ## 1           3     0.05 8.30e-06 Biological.Process
    ## 2           2     0.15 8.90e-03 Biological.Process
    ## 3           3     0.47 1.05e-02 Biological.Process
    ## 4           2     0.24 2.33e-02 Biological.Process
    ## 5           2     0.12 6.20e-03 Molecular.Function
    ## 6           1     0.01 1.37e-02 Molecular.Function
    ## 7           8     3.60 2.32e-02 Molecular.Function
    ## 8           1     0.04 4.06e-02 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 215 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  4 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 10:  9 nodes to be scored    (74 eliminated genes)

    ## 
    ##   Level 9:   16 nodes to be scored   (218 eliminated genes)

    ## 
    ##   Level 8:   21 nodes to be scored   (268 eliminated genes)

    ## 
    ##   Level 7:   26 nodes to be scored   (327 eliminated genes)

    ## 
    ##   Level 6:   31 nodes to be scored   (406 eliminated genes)

    ## 
    ##   Level 5:   36 nodes to be scored   (492 eliminated genes)

    ## 
    ##   Level 4:   31 nodes to be scored   (583 eliminated genes)

    ## 
    ##   Level 3:   27 nodes to be scored   (971 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (1189 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1336 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 130 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  3 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 9:   6 nodes to be scored    (12 eliminated genes)

    ## 
    ##   Level 8:   10 nodes to be scored   (145 eliminated genes)

    ## 
    ##   Level 7:   15 nodes to be scored   (441 eliminated genes)

    ## 
    ##   Level 6:   20 nodes to be scored   (510 eliminated genes)

    ## 
    ##   Level 5:   24 nodes to be scored   (835 eliminated genes)

    ## 
    ##   Level 4:   23 nodes to be scored   (950 eliminated genes)

    ## 
    ##   Level 3:   17 nodes to be scored   (1698 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (1936 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2386 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0000184 nuclear-transcribed mRNA catabolic process, nonsens...        10
    ## 2  GO:0001771                        immunological synapse formation         1
    ## 3  GO:0000183                         rDNA heterochromatin formation         1
    ## 4  GO:0006513                             protein monoubiquitination         1
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
    ## 1            3     0.20 0.00078 Biological.Process
    ## 2            1     0.02 0.01994 Biological.Process
    ## 3            1     0.02 0.01994 Biological.Process
    ## 4            1     0.02 0.01994 Biological.Process
    ## 5            3     0.42 0.02838 Biological.Process
    ## 6            1     0.04 0.03950 Biological.Process
    ## 7            3     0.56 0.01800 Molecular.Function
    ## 8            1     0.02 0.01900 Molecular.Function
    ## 9            1     0.02 0.01900 Molecular.Function
    ## 10           2     0.23 0.02000 Molecular.Function
    ## 11           6     2.29 0.02500 Molecular.Function
    ## 12           1     0.04 0.03700 Molecular.Function
    ## 13          14     9.32 0.03800 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 337 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  1 nodes to be scored    (2 eliminated genes)

    ## 
    ##   Level 12:  4 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 11:  8 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 10:  17 nodes to be scored   (76 eliminated genes)

    ## 
    ##   Level 9:   28 nodes to be scored   (241 eliminated genes)

    ## 
    ##   Level 8:   34 nodes to be scored   (333 eliminated genes)

    ## 
    ##   Level 7:   46 nodes to be scored   (412 eliminated genes)

    ## 
    ##   Level 6:   57 nodes to be scored   (500 eliminated genes)

    ## 
    ##   Level 5:   53 nodes to be scored   (759 eliminated genes)

    ## 
    ##   Level 4:   42 nodes to be scored   (857 eliminated genes)

    ## 
    ##   Level 3:   34 nodes to be scored   (1100 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (1255 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1390 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 168 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  2 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 9:   8 nodes to be scored    (10 eliminated genes)

    ## 
    ##   Level 8:   14 nodes to be scored   (43 eliminated genes)

    ## 
    ##   Level 7:   21 nodes to be scored   (329 eliminated genes)

    ## 
    ##   Level 6:   27 nodes to be scored   (378 eliminated genes)

    ## 
    ##   Level 5:   34 nodes to be scored   (648 eliminated genes)

    ## 
    ##   Level 4:   28 nodes to be scored   (934 eliminated genes)

    ## 
    ##   Level 3:   20 nodes to be scored   (1727 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (2003 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2566 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0001933         negative regulation of protein phosphorylation         5
    ## 2  GO:0001782                                     B cell homeostasis         7
    ## 3  GO:0001764                                       neuron migration         8
    ## 4  GO:0001731         formation of translation preinitiation complex         1
    ## 5  GO:0002026           regulation of the force of heart contraction         1
    ## 6  GO:0001960 negative regulation of cytokine-mediated signaling ...         1
    ## 7  GO:0001556                                      oocyte maturation         1
    ## 8  GO:0002092        positive regulation of receptor internalization         1
    ## 9  GO:0000038           very long-chain fatty acid metabolic process         9
    ## 10 GO:0005524                                            ATP binding       263
    ## 11 GO:0004197                   cysteine-type endopeptidase activity        18
    ## 12 GO:0004715 non-membrane spanning protein tyrosine kinase activ...         6
    ## 13 GO:0005220 inositol 1,4,5-trisphosphate-gated calcium channel ...         1
    ## 14 GO:0003874           6-pyruvoyltetrahydropterin synthase activity         1
    ## 15 GO:0004745        all-trans-retinol dehydrogenase (NAD+) activity         1
    ## 16 GO:0004423                         iduronate-2-sulfatase activity         1
    ## 17 GO:0004862       cAMP-dependent protein kinase inhibitor activity         1
    ## 18 GO:0015175 neutral L-amino acid transmembrane transporter acti...         1
    ## 19 GO:0005436                    sodium:phosphate symporter activity         1
    ## 20 GO:0003989                        acetyl-CoA carboxylase activity         1
    ## 21 GO:0005096                              GTPase activator activity        12
    ##    Significant Expected Fisher               type
    ## 1            2     0.17 0.0110 Biological.Process
    ## 2            2     0.24 0.0220 Biological.Process
    ## 3            2     0.28 0.0290 Biological.Process
    ## 4            1     0.03 0.0350 Biological.Process
    ## 5            1     0.03 0.0350 Biological.Process
    ## 6            1     0.03 0.0350 Biological.Process
    ## 7            1     0.03 0.0350 Biological.Process
    ## 8            1     0.03 0.0350 Biological.Process
    ## 9            2     0.31 0.0370 Biological.Process
    ## 10          17     7.59 0.0010 Molecular.Function
    ## 11           4     0.52 0.0079 Molecular.Function
    ## 12           2     0.17 0.0114 Molecular.Function
    ## 13           1     0.03 0.0288 Molecular.Function
    ## 14           1     0.03 0.0288 Molecular.Function
    ## 15           1     0.03 0.0288 Molecular.Function
    ## 16           1     0.03 0.0288 Molecular.Function
    ## 17           1     0.03 0.0288 Molecular.Function
    ## 18           1     0.03 0.0288 Molecular.Function
    ## 19           1     0.03 0.0288 Molecular.Function
    ## 20           1     0.03 0.0288 Molecular.Function
    ## 21           2     0.35 0.0450 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 47 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 9:   3 nodes to be scored    (5 eliminated genes)

    ## 
    ##   Level 8:   2 nodes to be scored    (35 eliminated genes)

    ## 
    ##   Level 7:   5 nodes to be scored    (50 eliminated genes)

    ## 
    ##   Level 6:   6 nodes to be scored    (119 eliminated genes)

    ## 
    ##   Level 5:   7 nodes to be scored    (276 eliminated genes)

    ## 
    ##   Level 4:   7 nodes to be scored    (348 eliminated genes)

    ## 
    ##   Level 3:   6 nodes to be scored    (593 eliminated genes)

    ## 
    ##   Level 2:   4 nodes to be scored    (662 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (849 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

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
    ##   Level 5:   2 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 4:   3 nodes to be scored    (272 eliminated genes)

    ## 
    ##   Level 3:   3 nodes to be scored    (303 eliminated genes)

    ## 
    ##   Level 2:   1 nodes to be scored    (485 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1008 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0000472 endonucleolytic cleavage to generate mature 5'-end ...         1
    ## 2 GO:0001889                                      liver development         8
    ## 3 GO:0010181                                            FMN binding         1
    ##   Significant Expected  Fisher               type
    ## 1           1     0.00 0.00210 Biological.Process
    ## 2           1     0.02 0.01700 Biological.Process
    ## 3           1     0.00 0.00036 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 82 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  1 nodes to be scored    (10 eliminated genes)

    ## 
    ##   Level 9:   2 nodes to be scored    (14 eliminated genes)

    ## 
    ##   Level 8:   3 nodes to be scored    (14 eliminated genes)

    ## 
    ##   Level 7:   9 nodes to be scored    (34 eliminated genes)

    ## 
    ##   Level 6:   15 nodes to be scored   (191 eliminated genes)

    ## 
    ##   Level 5:   16 nodes to be scored   (353 eliminated genes)

    ## 
    ##   Level 4:   14 nodes to be scored   (406 eliminated genes)

    ## 
    ##   Level 3:   13 nodes to be scored   (688 eliminated genes)

    ## 
    ##   Level 2:   6 nodes to be scored    (799 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1036 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 64 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   7 nodes to be scored    (375 eliminated genes)

    ## 
    ##   Level 6:   10 nodes to be scored   (418 eliminated genes)

    ## 
    ##   Level 5:   14 nodes to be scored   (663 eliminated genes)

    ## 
    ##   Level 4:   12 nodes to be scored   (728 eliminated genes)

    ## 
    ##   Level 3:   11 nodes to be scored   (1262 eliminated genes)

    ## 
    ##   Level 2:   4 nodes to be scored    (1592 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2066 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0001960 negative regulation of cytokine-mediated signaling ...         1
    ## 2 GO:0000302                    response to reactive oxygen species         6
    ## 3 GO:0000184 nuclear-transcribed mRNA catabolic process, nonsens...        10
    ## 4 GO:0005536                                      D-glucose binding         1
    ## 5 GO:0004492        methyl/ethyl malonyl-CoA decarboxylase activity         1
    ## 6 GO:0003964                   RNA-directed DNA polymerase activity       122
    ##   Significant Expected Fisher               type
    ## 1           1     0.00 0.0036 Biological.Process
    ## 2           1     0.02 0.0212 Biological.Process
    ## 3           1     0.04 0.0352 Biological.Process
    ## 4           1     0.01 0.0061 Molecular.Function
    ## 5           1     0.01 0.0061 Molecular.Function
    ## 6           3     0.75 0.0360 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 142 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  8 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 10:  9 nodes to be scored    (11 eliminated genes)

    ## 
    ##   Level 9:   7 nodes to be scored    (42 eliminated genes)

    ## 
    ##   Level 8:   10 nodes to be scored   (48 eliminated genes)

    ## 
    ##   Level 7:   17 nodes to be scored   (72 eliminated genes)

    ## 
    ##   Level 6:   23 nodes to be scored   (242 eliminated genes)

    ## 
    ##   Level 5:   25 nodes to be scored   (300 eliminated genes)

    ## 
    ##   Level 4:   18 nodes to be scored   (365 eliminated genes)

    ## 
    ##   Level 3:   13 nodes to be scored   (549 eliminated genes)

    ## 
    ##   Level 2:   6 nodes to be scored    (615 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (757 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

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
    ##   Level 6:   8 nodes to be scored    (264 eliminated genes)

    ## 
    ##   Level 5:   11 nodes to be scored   (392 eliminated genes)

    ## 
    ##   Level 4:   9 nodes to be scored    (500 eliminated genes)

    ## 
    ##   Level 3:   8 nodes to be scored    (1017 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (1325 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1900 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0007175 negative regulation of epidermal growth factor-acti...         1
    ## 2 GO:0002230 positive regulation of defense response to virus by...         1
    ## 3 GO:0002933                                    lipid hydroxylation         3
    ## 4 GO:0000184 nuclear-transcribed mRNA catabolic process, nonsens...        10
    ## 5 GO:0004753                    saccharopine dehydrogenase activity         1
    ## 6 GO:0000016                                       lactase activity         4
    ## 7 GO:0008376               acetylgalactosaminyltransferase activity        14
    ##   Significant Expected Fisher               type
    ## 1           1     0.00 0.0036 Biological.Process
    ## 2           1     0.00 0.0036 Biological.Process
    ## 3           1     0.01 0.0107 Biological.Process
    ## 4           1     0.04 0.0352 Biological.Process
    ## 5           1     0.00 0.0032 Molecular.Function
    ## 6           1     0.01 0.0129 Molecular.Function
    ## 7           1     0.05 0.0446 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 105 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (63 eliminated genes)

    ## 
    ##   Level 9:   7 nodes to be scored    (203 eliminated genes)

    ## 
    ##   Level 8:   12 nodes to be scored   (213 eliminated genes)

    ## 
    ##   Level 7:   12 nodes to be scored   (216 eliminated genes)

    ## 
    ##   Level 6:   15 nodes to be scored   (352 eliminated genes)

    ## 
    ##   Level 5:   17 nodes to be scored   (430 eliminated genes)

    ## 
    ##   Level 4:   12 nodes to be scored   (473 eliminated genes)

    ## 
    ##   Level 3:   14 nodes to be scored   (699 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (789 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1024 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 65 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   7 nodes to be scored    (263 eliminated genes)

    ## 
    ##   Level 6:   10 nodes to be scored   (267 eliminated genes)

    ## 
    ##   Level 5:   12 nodes to be scored   (464 eliminated genes)

    ## 
    ##   Level 4:   15 nodes to be scored   (554 eliminated genes)

    ## 
    ##   Level 3:   12 nodes to be scored   (1227 eliminated genes)

    ## 
    ##   Level 2:   5 nodes to be scored    (1778 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2312 eliminated genes)

    ##        GO.ID                              Term Annotated Significant Expected
    ## 1 GO:0001822                kidney development        65           7     0.79
    ## 2 GO:0000165                      MAPK cascade        35           4     0.42
    ## 3 GO:0004104           cholinesterase activity         1           1     0.01
    ## 4 GO:0004382          GDP phosphatase activity         3           1     0.03
    ## 5 GO:0000340 RNA 7-methylguanosine cap binding         4           1     0.04
    ## 6 GO:0005178                  integrin binding         4           1     0.04
    ## 7 GO:0005524                       ATP binding       263           6     2.66
    ##    Fisher               type
    ## 1 4.4e-06 Biological.Process
    ## 2 6.1e-04 Biological.Process
    ## 3 1.0e-02 Molecular.Function
    ## 4 3.0e-02 Molecular.Function
    ## 5 4.0e-02 Molecular.Function
    ## 6 4.0e-02 Molecular.Function
    ## 7 4.3e-02 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 78 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   7 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 8:   9 nodes to be scored    (167 eliminated genes)

    ## 
    ##   Level 7:   9 nodes to be scored    (189 eliminated genes)

    ## 
    ##   Level 6:   11 nodes to be scored   (202 eliminated genes)

    ## 
    ##   Level 5:   14 nodes to be scored   (286 eliminated genes)

    ## 
    ##   Level 4:   9 nodes to be scored    (320 eliminated genes)

    ## 
    ##   Level 3:   9 nodes to be scored    (409 eliminated genes)

    ## 
    ##   Level 2:   4 nodes to be scored    (507 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (656 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 41 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   4 nodes to be scored    (263 eliminated genes)

    ## 
    ##   Level 6:   5 nodes to be scored    (264 eliminated genes)

    ## 
    ##   Level 5:   7 nodes to be scored    (427 eliminated genes)

    ## 
    ##   Level 4:   10 nodes to be scored   (429 eliminated genes)

    ## 
    ##   Level 3:   9 nodes to be scored    (933 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (1563 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1907 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0015969             guanosine tetraphosphate metabolic process         9
    ## 2 GO:0000122 negative regulation of transcription by RNA polymer...       140
    ## 3 GO:0003964                   RNA-directed DNA polymerase activity       122
    ## 4 GO:0005178                                       integrin binding         4
    ##   Significant Expected Fisher               type
    ## 1           2     0.05  0.001 Biological.Process
    ## 2           3     0.80  0.037 Biological.Process
    ## 3           3     0.53  0.014 Molecular.Function
    ## 4           1     0.02  0.017 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 90 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   3 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 8:   8 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 7:   9 nodes to be scored    (142 eliminated genes)

    ## 
    ##   Level 6:   13 nodes to be scored   (262 eliminated genes)

    ## 
    ##   Level 5:   19 nodes to be scored   (400 eliminated genes)

    ## 
    ##   Level 4:   14 nodes to be scored   (434 eliminated genes)

    ## 
    ##   Level 3:   12 nodes to be scored   (669 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (844 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1072 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

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
    ##   Level 8:   6 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 7:   11 nodes to be scored   (375 eliminated genes)

    ## 
    ##   Level 6:   13 nodes to be scored   (417 eliminated genes)

    ## 
    ##   Level 5:   15 nodes to be scored   (715 eliminated genes)

    ## 
    ##   Level 4:   12 nodes to be scored   (760 eliminated genes)

    ## 
    ##   Level 3:   11 nodes to be scored   (1137 eliminated genes)

    ## 
    ##   Level 2:   5 nodes to be scored    (1503 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2009 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0006511          ubiquitin-dependent protein catabolic process         7
    ## 2 GO:0000703 oxidized pyrimidine nucleobase lesion DNA N-glycosy...         1
    ## 3 GO:0003980 UDP-glucose:glycoprotein glucosyltransferase activi...         1
    ## 4 GO:0001002 RNA polymerase III type 1 promoter sequence-specifi...         1
    ## 5 GO:0016836                                   hydro-lyase activity         4
    ## 6 GO:0000340                      RNA 7-methylguanosine cap binding         4
    ## 7 GO:0005524                                            ATP binding       263
    ##   Significant Expected Fisher               type
    ## 1           1     0.03 0.0300 Biological.Process
    ## 2           1     0.01 0.0076 Molecular.Function
    ## 3           1     0.01 0.0076 Molecular.Function
    ## 4           1     0.01 0.0076 Molecular.Function
    ## 5           1     0.03 0.0300 Molecular.Function
    ## 6           1     0.03 0.0300 Molecular.Function
    ## 7           5     1.99 0.0423 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

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

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 41 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   4 nodes to be scored    (375 eliminated genes)

    ## 
    ##   Level 6:   7 nodes to be scored    (394 eliminated genes)

    ## 
    ##   Level 5:   7 nodes to be scored    (438 eliminated genes)

    ## 
    ##   Level 4:   8 nodes to be scored    (444 eliminated genes)

    ## 
    ##   Level 3:   7 nodes to be scored    (609 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (1107 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1772 eliminated genes)

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

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

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
    ##   Level 1:   1 nodes to be scored    (422 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 54 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   3 nodes to be scored    (7 eliminated genes)

    ## 
    ##   Level 7:   7 nodes to be scored    (312 eliminated genes)

    ## 
    ##   Level 6:   9 nodes to be scored    (334 eliminated genes)

    ## 
    ##   Level 5:   9 nodes to be scored    (545 eliminated genes)

    ## 
    ##   Level 4:   10 nodes to be scored   (687 eliminated genes)

    ## 
    ##   Level 3:   8 nodes to be scored    (925 eliminated genes)

    ## 
    ##   Level 2:   4 nodes to be scored    (1077 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1889 eliminated genes)

    ##        GO.ID                                   Term Annotated Significant
    ## 1 GO:0008218                        bioluminescence        10           6
    ## 2 GO:0004089         carbonate dehydratase activity         2           1
    ## 3 GO:0005245 voltage-gated calcium channel activity         7           1
    ##   Expected  Fisher               type
    ## 1     0.04 2.0e-14 Biological.Process
    ## 2     0.01 5.0e-03 Molecular.Function
    ## 3     0.02 1.8e-02 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 252 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  5 nodes to be scored    (11 eliminated genes)

    ## 
    ##   Level 10:  9 nodes to be scored    (88 eliminated genes)

    ## 
    ##   Level 9:   17 nodes to be scored   (248 eliminated genes)

    ## 
    ##   Level 8:   25 nodes to be scored   (314 eliminated genes)

    ## 
    ##   Level 7:   30 nodes to be scored   (349 eliminated genes)

    ## 
    ##   Level 6:   41 nodes to be scored   (480 eliminated genes)

    ## 
    ##   Level 5:   46 nodes to be scored   (757 eliminated genes)

    ## 
    ##   Level 4:   33 nodes to be scored   (896 eliminated genes)

    ## 
    ##   Level 3:   30 nodes to be scored   (1042 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (1135 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1374 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 144 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   5 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   9 nodes to be scored    (142 eliminated genes)

    ## 
    ##   Level 7:   17 nodes to be scored   (441 eliminated genes)

    ## 
    ##   Level 6:   26 nodes to be scored   (482 eliminated genes)

    ## 
    ##   Level 5:   28 nodes to be scored   (779 eliminated genes)

    ## 
    ##   Level 4:   29 nodes to be scored   (1006 eliminated genes)

    ## 
    ##   Level 3:   19 nodes to be scored   (1617 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (2063 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2393 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0001822                                     kidney development        65
    ## 2 GO:0000165                                           MAPK cascade        35
    ## 3 GO:0001409 guanine nucleotide transmembrane transporter activi...        30
    ## 4 GO:0005524                                            ATP binding       263
    ## 5 GO:0005381            iron ion transmembrane transporter activity         1
    ## 6 GO:0004140                          dephospho-CoA kinase activity         1
    ## 7 GO:0005324 long-chain fatty acid transmembrane transporter act...         1
    ## 8 GO:0003777                             microtubule motor activity        16
    ##   Significant Expected   Fisher               type
    ## 1          12     2.27 8.40e-07 Biological.Process
    ## 2           4     1.22 3.10e-02 Biological.Process
    ## 3           4     0.65 3.60e-03 Molecular.Function
    ## 4          12     5.69 9.10e-03 Molecular.Function
    ## 5           1     0.02 2.16e-02 Molecular.Function
    ## 6           1     0.02 2.16e-02 Molecular.Function
    ## 7           1     0.02 2.16e-02 Molecular.Function
    ## 8           2     0.35 4.55e-02 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 22 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   2 nodes to be scored    (30 eliminated genes)

    ## 
    ##   Level 7:   3 nodes to be scored    (37 eliminated genes)

    ## 
    ##   Level 6:   2 nodes to be scored    (38 eliminated genes)

    ## 
    ##   Level 5:   2 nodes to be scored    (44 eliminated genes)

    ## 
    ##   Level 4:   4 nodes to be scored    (66 eliminated genes)

    ## 
    ##   Level 3:   4 nodes to be scored    (124 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (164 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (522 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 35 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   5 nodes to be scored    (263 eliminated genes)

    ## 
    ##   Level 6:   6 nodes to be scored    (283 eliminated genes)

    ## 
    ##   Level 5:   6 nodes to be scored    (330 eliminated genes)

    ## 
    ##   Level 4:   6 nodes to be scored    (461 eliminated genes)

    ## 
    ##   Level 3:   5 nodes to be scored    (556 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (1065 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1772 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0000052                           citrulline metabolic process         1
    ## 2 GO:0000209                             protein polyubiquitination        30
    ## 3 GO:0003980 UDP-glucose:glycoprotein glucosyltransferase activi...         1
    ## 4 GO:0004559                             alpha-mannosidase activity         7
    ## 5 GO:0005507                                     copper ion binding        18
    ##   Significant Expected Fisher               type
    ## 1           1     0.00 0.0014 Biological.Process
    ## 2           1     0.04 0.0423 Biological.Process
    ## 3           1     0.00 0.0025 Molecular.Function
    ## 4           1     0.02 0.0176 Molecular.Function
    ## 5           1     0.05 0.0446 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

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

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

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
    ##   Level 8:   1 nodes to be scored    (2 eliminated genes)

    ## 
    ##   Level 7:   4 nodes to be scored    (2 eliminated genes)

    ## 
    ##   Level 6:   2 nodes to be scored    (6 eliminated genes)

    ## 
    ##   Level 5:   6 nodes to be scored    (60 eliminated genes)

    ## 
    ##   Level 4:   6 nodes to be scored    (90 eliminated genes)

    ## 
    ##   Level 3:   3 nodes to be scored    (519 eliminated genes)

    ## 
    ##   Level 2:   4 nodes to be scored    (914 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1180 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0005391 P-type sodium:potassium-exchanging transporter acti...         2
    ##   Significant Expected Fisher               type
    ## 1           1        0 0.0036 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 352 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  6 nodes to be scored    (11 eliminated genes)

    ## 
    ##   Level 11:  8 nodes to be scored    (27 eliminated genes)

    ## 
    ##   Level 10:  17 nodes to be scored   (93 eliminated genes)

    ## 
    ##   Level 9:   25 nodes to be scored   (255 eliminated genes)

    ## 
    ##   Level 8:   33 nodes to be scored   (348 eliminated genes)

    ## 
    ##   Level 7:   48 nodes to be scored   (424 eliminated genes)

    ## 
    ##   Level 6:   58 nodes to be scored   (532 eliminated genes)

    ## 
    ##   Level 5:   57 nodes to be scored   (855 eliminated genes)

    ## 
    ##   Level 4:   45 nodes to be scored   (927 eliminated genes)

    ## 
    ##   Level 3:   38 nodes to be scored   (1129 eliminated genes)

    ## 
    ##   Level 2:   11 nodes to be scored   (1275 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1375 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 267 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  5 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   14 nodes to be scored   (3 eliminated genes)

    ## 
    ##   Level 8:   18 nodes to be scored   (153 eliminated genes)

    ## 
    ##   Level 7:   40 nodes to be scored   (465 eliminated genes)

    ## 
    ##   Level 6:   53 nodes to be scored   (615 eliminated genes)

    ## 
    ##   Level 5:   51 nodes to be scored   (1063 eliminated genes)

    ## 
    ##   Level 4:   46 nodes to be scored   (1374 eliminated genes)

    ## 
    ##   Level 3:   27 nodes to be scored   (2047 eliminated genes)

    ## 
    ##   Level 2:   11 nodes to be scored   (2282 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2568 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0002064                            epithelial cell development         4
    ## 2  GO:0000122 negative regulation of transcription by RNA polymer...       140
    ## 3  GO:0001732 formation of cytoplasmic translation initiation com...         6
    ## 4  GO:0005388                    P-type calcium transporter activity         4
    ## 5  GO:0005198                           structural molecule activity        54
    ## 6  GO:0003756                   protein disulfide isomerase activity         5
    ## 7  GO:0001786                             phosphatidylserine binding         4
    ## 8  GO:0002020                                       protease binding         5
    ## 9  GO:0005539                              glycosaminoglycan binding         5
    ## 10 GO:0003682                                      chromatin binding        13
    ## 11 GO:0004843                  cysteine-type deubiquitinase activity         6
    ##    Significant Expected  Fisher               type
    ## 1            3     0.24 0.00076 Biological.Process
    ## 2           14     8.28 0.03037 Biological.Process
    ## 3            2     0.35 0.04434 Biological.Process
    ## 4            4     0.23 0.00001 Molecular.Function
    ## 5           10     3.08 0.00025 Molecular.Function
    ## 6            3     0.28 0.00167 Molecular.Function
    ## 7            2     0.23 0.01794 Molecular.Function
    ## 8            2     0.28 0.02879 Molecular.Function
    ## 9            2     0.28 0.02879 Molecular.Function
    ## 10           3     0.74 0.03403 Molecular.Function
    ## 11           2     0.34 0.04158 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

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
    ##   Level 8:   3 nodes to be scored    (4 eliminated genes)

    ## 
    ##   Level 7:   7 nodes to be scored    (5 eliminated genes)

    ## 
    ##   Level 6:   10 nodes to be scored   (8 eliminated genes)

    ## 
    ##   Level 5:   10 nodes to be scored   (99 eliminated genes)

    ## 
    ##   Level 4:   7 nodes to be scored    (175 eliminated genes)

    ## 
    ##   Level 3:   6 nodes to be scored    (282 eliminated genes)

    ## 
    ##   Level 2:   5 nodes to be scored    (310 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (793 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

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
    ##   Level 6:   8 nodes to be scored    (270 eliminated genes)

    ## 
    ##   Level 5:   8 nodes to be scored    (425 eliminated genes)

    ## 
    ##   Level 4:   12 nodes to be scored   (512 eliminated genes)

    ## 
    ##   Level 3:   9 nodes to be scored    (737 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (1470 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2070 eliminated genes)

    ##        GO.ID                                  Term Annotated Significant
    ## 1 GO:0001947                         heart looping         4           1
    ## 2 GO:0000902                    cell morphogenesis         6           1
    ## 3 GO:0019700 organic phosphonate catabolic process        12           1
    ## 4 GO:0001666                   response to hypoxia        14           1
    ## 5 GO:0000049                          tRNA binding        30           2
    ## 6 GO:0004565           beta-galactosidase activity         3           1
    ## 7 GO:0005525                           GTP binding         6           1
    ##   Expected Fisher               type
    ## 1     0.01  0.014 Biological.Process
    ## 2     0.02  0.021 Biological.Process
    ## 3     0.04  0.042 Biological.Process
    ## 4     0.05  0.049 Biological.Process
    ## 5     0.16  0.011 Molecular.Function
    ## 6     0.02  0.016 Molecular.Function
    ## 7     0.03  0.032 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 38 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   4 nodes to be scored    (9 eliminated genes)

    ## 
    ##   Level 7:   3 nodes to be scored    (9 eliminated genes)

    ## 
    ##   Level 6:   4 nodes to be scored    (10 eliminated genes)

    ## 
    ##   Level 5:   6 nodes to be scored    (114 eliminated genes)

    ## 
    ##   Level 4:   5 nodes to be scored    (122 eliminated genes)

    ## 
    ##   Level 3:   8 nodes to be scored    (204 eliminated genes)

    ## 
    ##   Level 2:   4 nodes to be scored    (535 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (766 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

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
    ##   Level 6:   5 nodes to be scored    (264 eliminated genes)

    ## 
    ##   Level 5:   4 nodes to be scored    (300 eliminated genes)

    ## 
    ##   Level 4:   6 nodes to be scored    (332 eliminated genes)

    ## 
    ##   Level 3:   6 nodes to be scored    (370 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (665 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1724 eliminated genes)

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

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 48 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   5 nodes to be scored    (9 eliminated genes)

    ## 
    ##   Level 7:   5 nodes to be scored    (9 eliminated genes)

    ## 
    ##   Level 6:   5 nodes to be scored    (45 eliminated genes)

    ## 
    ##   Level 5:   6 nodes to be scored    (122 eliminated genes)

    ## 
    ##   Level 4:   7 nodes to be scored    (132 eliminated genes)

    ## 
    ##   Level 3:   10 nodes to be scored   (318 eliminated genes)

    ## 
    ##   Level 2:   6 nodes to be scored    (649 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (900 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 32 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   4 nodes to be scored    (263 eliminated genes)

    ## 
    ##   Level 6:   5 nodes to be scored    (264 eliminated genes)

    ## 
    ##   Level 5:   4 nodes to be scored    (427 eliminated genes)

    ## 
    ##   Level 4:   7 nodes to be scored    (429 eliminated genes)

    ## 
    ##   Level 3:   6 nodes to be scored    (474 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (1306 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1585 eliminated genes)

    ##        GO.ID                                       Term Annotated Significant
    ## 1 GO:0001822                         kidney development        65           3
    ## 2 GO:0000165                               MAPK cascade        35           2
    ## 3 GO:0015969 guanosine tetraphosphate metabolic process         9           1
    ##   Expected Fisher               type
    ## 1     0.32 0.0029 Biological.Process
    ## 2     0.17 0.0117 Biological.Process
    ## 3     0.04 0.0441 Biological.Process
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

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
    ##   Level 3:   1 nodes to be scored    (222 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (232 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (245 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 32 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 7:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   6 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 5:   7 nodes to be scored    (75 eliminated genes)

    ## 
    ##   Level 4:   8 nodes to be scored    (246 eliminated genes)

    ## 
    ##   Level 3:   5 nodes to be scored    (733 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (1245 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1709 eliminated genes)

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
    ## 4           1     0.02  0.023 Molecular.Function
    ## 5           1     0.05  0.045 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 57 nontrivial nodes
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
    ##   Level 5:   9 nodes to be scored    (62 eliminated genes)

    ## 
    ##   Level 4:   10 nodes to be scored   (82 eliminated genes)

    ## 
    ##   Level 3:   9 nodes to be scored    (180 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (274 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (723 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 42 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 7:   3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   7 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 5:   9 nodes to be scored    (106 eliminated genes)

    ## 
    ##   Level 4:   11 nodes to be scored   (301 eliminated genes)

    ## 
    ##   Level 3:   8 nodes to be scored    (908 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (1262 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1868 eliminated genes)

    ##        GO.ID                                  Term Annotated Significant
    ## 1 GO:0000028      ribosomal small subunit assembly         3           1
    ## 2 GO:0000054 ribosomal subunit export from nucleus         3           1
    ## 3 GO:0002933                   lipid hydroxylation         3           1
    ## 4 GO:0001561            fatty acid alpha-oxidation         3           1
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

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 147 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 13:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 11:  2 nodes to be scored    (11 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (11 eliminated genes)

    ## 
    ##   Level 9:   7 nodes to be scored    (166 eliminated genes)

    ## 
    ##   Level 8:   12 nodes to be scored   (184 eliminated genes)

    ## 
    ##   Level 7:   17 nodes to be scored   (216 eliminated genes)

    ## 
    ##   Level 6:   22 nodes to be scored   (304 eliminated genes)

    ## 
    ##   Level 5:   27 nodes to be scored   (548 eliminated genes)

    ## 
    ##   Level 4:   22 nodes to be scored   (580 eliminated genes)

    ## 
    ##   Level 3:   23 nodes to be scored   (710 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (846 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1059 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 95 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   5 nodes to be scored    (100 eliminated genes)

    ## 
    ##   Level 7:   11 nodes to be scored   (385 eliminated genes)

    ## 
    ##   Level 6:   18 nodes to be scored   (438 eliminated genes)

    ## 
    ##   Level 5:   21 nodes to be scored   (741 eliminated genes)

    ## 
    ##   Level 4:   16 nodes to be scored   (966 eliminated genes)

    ## 
    ##   Level 3:   11 nodes to be scored   (1645 eliminated genes)

    ## 
    ##   Level 2:   7 nodes to be scored    (1837 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2335 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0000050                                             urea cycle         2
    ## 2 GO:0003676                                   nucleic acid binding       497
    ## 3 GO:0005085             guanyl-nucleotide exchange factor activity        12
    ## 4 GO:0004143           ATP-dependent diacylglycerol kinase activity         1
    ## 5 GO:0004568                                     chitinase activity         1
    ## 6 GO:0004571 mannosyl-oligosaccharide 1,2-alpha-mannosidase acti...         2
    ##   Significant Expected Fisher               type
    ## 1           1     0.02 0.0170 Biological.Process
    ## 2          12     6.99 0.0076 Molecular.Function
    ## 3           2     0.17 0.0116 Molecular.Function
    ## 4           1     0.01 0.0141 Molecular.Function
    ## 5           1     0.01 0.0141 Molecular.Function
    ## 6           1     0.03 0.0279 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 69 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   4 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 8:   7 nodes to be scored    (170 eliminated genes)

    ## 
    ##   Level 7:   7 nodes to be scored    (179 eliminated genes)

    ## 
    ##   Level 6:   11 nodes to be scored   (248 eliminated genes)

    ## 
    ##   Level 5:   13 nodes to be scored   (324 eliminated genes)

    ## 
    ##   Level 4:   9 nodes to be scored    (372 eliminated genes)

    ## 
    ##   Level 3:   9 nodes to be scored    (424 eliminated genes)

    ## 
    ##   Level 2:   4 nodes to be scored    (523 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (680 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

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

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

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
    ##   Level 5:   1 nodes to be scored    (12 eliminated genes)

    ## 
    ##   Level 4:   3 nodes to be scored    (13 eliminated genes)

    ## 
    ##   Level 3:   3 nodes to be scored    (39 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (74 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (452 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 69 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   4 nodes to be scored    (7 eliminated genes)

    ## 
    ##   Level 7:   9 nodes to be scored    (312 eliminated genes)

    ## 
    ##   Level 6:   13 nodes to be scored   (339 eliminated genes)

    ## 
    ##   Level 5:   13 nodes to be scored   (580 eliminated genes)

    ## 
    ##   Level 4:   13 nodes to be scored   (816 eliminated genes)

    ## 
    ##   Level 3:   9 nodes to be scored    (1095 eliminated genes)

    ## 
    ##   Level 2:   4 nodes to be scored    (1178 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2053 eliminated genes)

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
    ## 4           1     0.01 0.0094 Molecular.Function
    ## 5           1     0.02 0.0232 Molecular.Function
    ## 6           1     0.03 0.0324 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 366 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 12:  6 nodes to be scored    (11 eliminated genes)

    ## 
    ##   Level 11:  9 nodes to be scored    (27 eliminated genes)

    ## 
    ##   Level 10:  18 nodes to be scored   (87 eliminated genes)

    ## 
    ##   Level 9:   24 nodes to be scored   (251 eliminated genes)

    ## 
    ##   Level 8:   36 nodes to be scored   (355 eliminated genes)

    ## 
    ##   Level 7:   53 nodes to be scored   (411 eliminated genes)

    ## 
    ##   Level 6:   59 nodes to be scored   (493 eliminated genes)

    ## 
    ##   Level 5:   64 nodes to be scored   (619 eliminated genes)

    ## 
    ##   Level 4:   46 nodes to be scored   (816 eliminated genes)

    ## 
    ##   Level 3:   35 nodes to be scored   (1144 eliminated genes)

    ## 
    ##   Level 2:   10 nodes to be scored   (1294 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1387 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 222 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   16 nodes to be scored   (42 eliminated genes)

    ## 
    ##   Level 7:   35 nodes to be scored   (327 eliminated genes)

    ## 
    ##   Level 6:   47 nodes to be scored   (569 eliminated genes)

    ## 
    ##   Level 5:   44 nodes to be scored   (982 eliminated genes)

    ## 
    ##   Level 4:   37 nodes to be scored   (1307 eliminated genes)

    ## 
    ##   Level 3:   25 nodes to be scored   (1979 eliminated genes)

    ## 
    ##   Level 2:   12 nodes to be scored   (2278 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2616 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0001508                                       action potential         4
    ## 2  GO:0001732 formation of cytoplasmic translation initiation com...         6
    ## 3  GO:0000183                         rDNA heterochromatin formation         1
    ## 4  GO:0006777         Mo-molybdopterin cofactor biosynthetic process         1
    ## 5  GO:0000454             snoRNA guided rRNA pseudouridine synthesis         1
    ## 6  GO:0001937 negative regulation of endothelial cell proliferati...         1
    ## 7  GO:0002943                          tRNA dihydrouridine synthesis         1
    ## 8  GO:0000105                       L-histidine biosynthetic process         1
    ## 9  GO:0002100               tRNA wobble adenosine to inosine editing         1
    ## 10 GO:0001771                        immunological synapse formation         1
    ## 11 GO:0007601                                      visual perception         1
    ## 12 GO:0001706                                     endoderm formation         1
    ## 13 GO:0001967                                      suckling behavior         1
    ## 14 GO:0000398                         mRNA splicing, via spliceosome        63
    ## 15 GO:0003735                     structural constituent of ribosome        13
    ## 16 GO:0001786                             phosphatidylserine binding         4
    ## 17 GO:0003964                   RNA-directed DNA polymerase activity       122
    ## 18 GO:0008768                       UDP-sugar diphosphatase activity         1
    ## 19 GO:0004736                          pyruvate carboxylase activity         1
    ## 20 GO:0004568                                     chitinase activity         1
    ## 21 GO:0004114     3',5'-cyclic-nucleotide phosphodiesterase activity         1
    ## 22 GO:0000033                 alpha-1,3-mannosyltransferase activity         1
    ## 23 GO:0000009                 alpha-1,6-mannosyltransferase activity         1
    ## 24 GO:0004108                         citrate (Si)-synthase activity         1
    ## 25 GO:0003944 N-acetylglucosamine-1-phosphodiester alpha-N-acetyl...         1
    ## 26 GO:0004856                                D-xylulokinase activity         1
    ## 27 GO:0005536                                      D-glucose binding         1
    ## 28 GO:0008663 2',3'-cyclic-nucleotide 2'-phosphodiesterase activi...         1
    ## 29 GO:0004096                                      catalase activity         1
    ## 30 GO:0004638      phosphoribosylaminoimidazole carboxylase activity         1
    ##    Significant Expected Fisher               type
    ## 1            2     0.17 0.0110 Biological.Process
    ## 2            2     0.26 0.0250 Biological.Process
    ## 3            1     0.04 0.0430 Biological.Process
    ## 4            1     0.04 0.0430 Biological.Process
    ## 5            1     0.04 0.0430 Biological.Process
    ## 6            1     0.04 0.0430 Biological.Process
    ## 7            1     0.04 0.0430 Biological.Process
    ## 8            1     0.04 0.0430 Biological.Process
    ## 9            1     0.04 0.0430 Biological.Process
    ## 10           1     0.04 0.0430 Biological.Process
    ## 11           1     0.04 0.0430 Biological.Process
    ## 12           1     0.04 0.0430 Biological.Process
    ## 13           1     0.04 0.0430 Biological.Process
    ## 14           7     2.74 0.0490 Biological.Process
    ## 15           4     0.54 0.0016 Molecular.Function
    ## 16           2     0.17 0.0099 Molecular.Function
    ## 17          11     5.10 0.0117 Molecular.Function
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
    ## 29           1     0.04 0.0418 Molecular.Function
    ## 30           1     0.04 0.0418 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

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
    ##   Level 2:   1 nodes to be scored    (93 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (106 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

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

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 107 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   7 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 8:   12 nodes to be scored   (178 eliminated genes)

    ## 
    ##   Level 7:   14 nodes to be scored   (201 eliminated genes)

    ## 
    ##   Level 6:   18 nodes to be scored   (294 eliminated genes)

    ## 
    ##   Level 5:   21 nodes to be scored   (390 eliminated genes)

    ## 
    ##   Level 4:   14 nodes to be scored   (434 eliminated genes)

    ## 
    ##   Level 3:   9 nodes to be scored    (731 eliminated genes)

    ## 
    ##   Level 2:   6 nodes to be scored    (961 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1060 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 122 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   4 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   5 nodes to be scored    (107 eliminated genes)

    ## 
    ##   Level 7:   16 nodes to be scored   (432 eliminated genes)

    ## 
    ##   Level 6:   23 nodes to be scored   (466 eliminated genes)

    ## 
    ##   Level 5:   24 nodes to be scored   (763 eliminated genes)

    ## 
    ##   Level 4:   23 nodes to be scored   (922 eliminated genes)

    ## 
    ##   Level 3:   16 nodes to be scored   (1395 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (1898 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2438 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0008218                                        bioluminescence        10
    ## 2  GO:0001731         formation of translation preinitiation complex         1
    ## 3  GO:0002790                                      peptide secretion         2
    ## 4  GO:0000028                       ribosomal small subunit assembly         3
    ## 5  GO:0002933                                    lipid hydroxylation         3
    ## 6  GO:0004142      diacylglycerol cholinephosphotransferase activity         1
    ## 7  GO:0004325                                ferrochelatase activity         1
    ## 8  GO:0004368 glycerol-3-phosphate dehydrogenase (quinone) activi...         1
    ## 9  GO:0000774             adenyl-nucleotide exchange factor activity         2
    ## 10 GO:0004089                         carbonate dehydratase activity         2
    ## 11 GO:0003985                acetyl-CoA C-acetyltransferase activity         2
    ## 12 GO:0004176                       ATP-dependent peptidase activity         2
    ## 13 GO:0004252                     serine-type endopeptidase activity        53
    ## 14 GO:0005319                             lipid transporter activity         3
    ##    Significant Expected   Fisher               type
    ## 1            3     0.08 0.000042 Biological.Process
    ## 2            1     0.01 0.007800 Biological.Process
    ## 3            1     0.02 0.015600 Biological.Process
    ## 4            1     0.02 0.023300 Biological.Process
    ## 5            1     0.02 0.023300 Biological.Process
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

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

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
    ##   Level 3:   3 nodes to be scored    (238 eliminated genes)

    ## 
    ##   Level 2:   1 nodes to be scored    (283 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (300 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 18 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 7:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 6:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 5:   3 nodes to be scored    (188 eliminated genes)

    ## 
    ##   Level 4:   4 nodes to be scored    (303 eliminated genes)

    ## 
    ##   Level 3:   4 nodes to be scored    (443 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (1222 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1405 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0000972 transcription-dependent tethering of RNA polymerase...         1
    ## 2 GO:0003341                                        cilium movement        20
    ##   Significant Expected Fisher               type
    ## 1           1     0.00 0.0021 Biological.Process
    ## 2           1     0.04 0.0422 Biological.Process
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 107 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   5 nodes to be scored    (166 eliminated genes)

    ## 
    ##   Level 8:   8 nodes to be scored    (166 eliminated genes)

    ## 
    ##   Level 7:   10 nodes to be scored   (170 eliminated genes)

    ## 
    ##   Level 6:   18 nodes to be scored   (277 eliminated genes)

    ## 
    ##   Level 5:   22 nodes to be scored   (350 eliminated genes)

    ## 
    ##   Level 4:   16 nodes to be scored   (388 eliminated genes)

    ## 
    ##   Level 3:   13 nodes to be scored   (469 eliminated genes)

    ## 
    ##   Level 2:   9 nodes to be scored    (789 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1085 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

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
    ##   Level 3:   6 nodes to be scored    (382 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (670 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1621 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0001755                            neural crest cell migration         2
    ## 2 GO:0001002 RNA polymerase III type 1 promoter sequence-specifi...         1
    ## 3 GO:0004660                   protein farnesyltransferase activity         2
    ## 4 GO:0004719 protein-L-isoaspartate (D-aspartate) O-methyltransf...         3
    ## 5 GO:0004867           serine-type endopeptidase inhibitor activity         4
    ##   Significant Expected Fisher               type
    ## 1           1     0.01 0.0057 Biological.Process
    ## 2           1     0.00 0.0022 Molecular.Function
    ## 3           1     0.00 0.0043 Molecular.Function
    ## 4           1     0.01 0.0065 Molecular.Function
    ## 5           1     0.01 0.0086 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 49 nontrivial nodes
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
    ##   Level 4:   6 nodes to be scored    (88 eliminated genes)

    ## 
    ##   Level 3:   5 nodes to be scored    (176 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (460 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (602 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 56 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   4 nodes to be scored    (7 eliminated genes)

    ## 
    ##   Level 7:   7 nodes to be scored    (49 eliminated genes)

    ## 
    ##   Level 6:   7 nodes to be scored    (137 eliminated genes)

    ## 
    ##   Level 5:   9 nodes to be scored    (291 eliminated genes)

    ## 
    ##   Level 4:   13 nodes to be scored   (375 eliminated genes)

    ## 
    ##   Level 3:   10 nodes to be scored   (736 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (1103 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (2043 eliminated genes)

    ##        GO.ID                                        Term Annotated Significant
    ## 1 GO:0001516          prostaglandin biosynthetic process         2           1
    ## 2 GO:0000028            ribosomal small subunit assembly         3           1
    ## 3 GO:0001523                  retinoid metabolic process         9           1
    ## 4 GO:0005381 iron ion transmembrane transporter activity         1           1
    ## 5 GO:0004252          serine-type endopeptidase activity        53           2
    ## 6 GO:0005245      voltage-gated calcium channel activity         7           1
    ## 7 GO:0005542                          folic acid binding        10           1
    ##   Expected Fisher               type
    ## 1     0.00 0.0043 Biological.Process
    ## 2     0.01 0.0064 Biological.Process
    ## 3     0.02 0.0191 Biological.Process
    ## 4     0.00 0.0029 Molecular.Function
    ## 5     0.15 0.0093 Molecular.Function
    ## 6     0.02 0.0200 Molecular.Function
    ## 7     0.03 0.0285 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

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
    ##   Level 8:   1 nodes to be scored    (64 eliminated genes)

    ## 
    ##   Level 7:   2 nodes to be scored    (65 eliminated genes)

    ## 
    ##   Level 6:   3 nodes to be scored    (82 eliminated genes)

    ## 
    ##   Level 5:   4 nodes to be scored    (239 eliminated genes)

    ## 
    ##   Level 4:   3 nodes to be scored    (292 eliminated genes)

    ## 
    ##   Level 3:   5 nodes to be scored    (353 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (375 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (586 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

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
    ##   Level 6:   5 nodes to be scored    (264 eliminated genes)

    ## 
    ##   Level 5:   4 nodes to be scored    (270 eliminated genes)

    ## 
    ##   Level 4:   4 nodes to be scored    (438 eliminated genes)

    ## 
    ##   Level 3:   3 nodes to be scored    (695 eliminated genes)

    ## 
    ##   Level 2:   2 nodes to be scored    (982 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1008 eliminated genes)

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

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 95 nontrivial nodes
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
    ##   Level 9:   5 nodes to be scored    (66 eliminated genes)

    ## 
    ##   Level 8:   7 nodes to be scored    (74 eliminated genes)

    ## 
    ##   Level 7:   10 nodes to be scored   (76 eliminated genes)

    ## 
    ##   Level 6:   12 nodes to be scored   (292 eliminated genes)

    ## 
    ##   Level 5:   17 nodes to be scored   (400 eliminated genes)

    ## 
    ##   Level 4:   12 nodes to be scored   (491 eliminated genes)

    ## 
    ##   Level 3:   10 nodes to be scored   (798 eliminated genes)

    ## 
    ##   Level 2:   6 nodes to be scored    (868 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1056 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

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
    ##   Level 6:   6 nodes to be scored    (264 eliminated genes)

    ## 
    ##   Level 5:   6 nodes to be scored    (281 eliminated genes)

    ## 
    ##   Level 4:   6 nodes to be scored    (301 eliminated genes)

    ## 
    ##   Level 3:   6 nodes to be scored    (367 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (650 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1648 eliminated genes)

    ##        GO.ID                                                   Term Annotated
    ## 1 GO:0000076                   DNA replication checkpoint signaling         1
    ## 2 GO:0001516                     prostaglandin biosynthetic process         2
    ## 3 GO:0000381 regulation of alternative mRNA splicing, via splice...        11
    ## 4 GO:0004777 succinate-semialdehyde dehydrogenase (NAD+) activit...         1
    ## 5 GO:0004601                                    peroxidase activity         7
    ## 6 GO:0004181                       metallocarboxypeptidase activity        10
    ##   Significant Expected Fisher               type
    ## 1           1     0.00 0.0036 Biological.Process
    ## 2           1     0.01 0.0071 Biological.Process
    ## 3           1     0.04 0.0386 Biological.Process
    ## 4           1     0.00 0.0018 Molecular.Function
    ## 5           1     0.01 0.0126 Molecular.Function
    ## 6           1     0.02 0.0179 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 108 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 11:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 10:  3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   5 nodes to be scored    (140 eliminated genes)

    ## 
    ##   Level 8:   8 nodes to be scored    (141 eliminated genes)

    ## 
    ##   Level 7:   13 nodes to be scored   (152 eliminated genes)

    ## 
    ##   Level 6:   17 nodes to be scored   (218 eliminated genes)

    ## 
    ##   Level 5:   20 nodes to be scored   (320 eliminated genes)

    ## 
    ##   Level 4:   18 nodes to be scored   (349 eliminated genes)

    ## 
    ##   Level 3:   16 nodes to be scored   (562 eliminated genes)

    ## 
    ##   Level 2:   6 nodes to be scored    (731 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (886 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 71 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   3 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   2 nodes to be scored    (100 eliminated genes)

    ## 
    ##   Level 7:   6 nodes to be scored    (383 eliminated genes)

    ## 
    ##   Level 6:   14 nodes to be scored   (394 eliminated genes)

    ## 
    ##   Level 5:   15 nodes to be scored   (480 eliminated genes)

    ## 
    ##   Level 4:   15 nodes to be scored   (717 eliminated genes)

    ## 
    ##   Level 3:   11 nodes to be scored   (1059 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (1468 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1995 eliminated genes)

    ##         GO.ID                                                   Term Annotated
    ## 1  GO:0002933                                    lipid hydroxylation         3
    ## 2  GO:0001731         formation of translation preinitiation complex         1
    ## 3  GO:0000028                       ribosomal small subunit assembly         3
    ## 4  GO:0009734                      auxin-activated signaling pathway         4
    ## 5  GO:0004342             glucosamine-6-phosphate deaminase activity         2
    ## 6  GO:0004142      diacylglycerol cholinephosphotransferase activity         1
    ## 7  GO:0004066 asparagine synthase (glutamine-hydrolyzing) activit...         1
    ## 8  GO:0004089                         carbonate dehydratase activity         2
    ## 9  GO:0004252                     serine-type endopeptidase activity        53
    ## 10 GO:0004601                                    peroxidase activity         7
    ##    Significant Expected   Fisher               type
    ## 1            2     0.02 0.000110 Biological.Process
    ## 2            1     0.01 0.006410 Biological.Process
    ## 3            1     0.02 0.019120 Biological.Process
    ## 4            1     0.03 0.025420 Biological.Process
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

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

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
    ##   Level 3:   4 nodes to be scored    (354 eliminated genes)

    ## 
    ##   Level 2:   5 nodes to be scored    (499 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (932 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

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

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

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

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

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
    ##   Level 6:   3 nodes to be scored    (264 eliminated genes)

    ## 
    ##   Level 5:   3 nodes to be scored    (270 eliminated genes)

    ## 
    ##   Level 4:   4 nodes to be scored    (272 eliminated genes)

    ## 
    ##   Level 3:   5 nodes to be scored    (304 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (486 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1394 eliminated genes)

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

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 102 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 15:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 14:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 13:  1 nodes to be scored    (1 eliminated genes)

    ## 
    ##   Level 12:  2 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 11:  3 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 10:  2 nodes to be scored    (3 eliminated genes)

    ## 
    ##   Level 9:   4 nodes to be scored    (14 eliminated genes)

    ## 
    ##   Level 8:   8 nodes to be scored    (14 eliminated genes)

    ## 
    ##   Level 7:   12 nodes to be scored   (42 eliminated genes)

    ## 
    ##   Level 6:   16 nodes to be scored   (224 eliminated genes)

    ## 
    ##   Level 5:   18 nodes to be scored   (302 eliminated genes)

    ## 
    ##   Level 4:   13 nodes to be scored   (357 eliminated genes)

    ## 
    ##   Level 3:   12 nodes to be scored   (619 eliminated genes)

    ## 
    ##   Level 2:   8 nodes to be scored    (775 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (983 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 49 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 9:   1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 7:   4 nodes to be scored    (263 eliminated genes)

    ## 
    ##   Level 6:   9 nodes to be scored    (282 eliminated genes)

    ## 
    ##   Level 5:   10 nodes to be scored   (332 eliminated genes)

    ## 
    ##   Level 4:   11 nodes to be scored   (556 eliminated genes)

    ## 
    ##   Level 3:   8 nodes to be scored    (922 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (1280 eliminated genes)

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
    ## 2           1     0.03 0.0296 Biological.Process
    ## 3           1     0.03 0.0338 Biological.Process
    ## 4           1     0.04 0.0379 Biological.Process
    ## 5           1     0.04 0.0421 Biological.Process
    ## 6           2     0.08 0.0030 Molecular.Function
    ##  Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
    ##  - attr(*, "names")= chr [1:30094] "FUN_002326" "FUN_002315" "FUN_002316" "FUN_002303" ...

    ## 
    ## Building most specific GOs .....

    ##  ( 266 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 1149 GO terms and 2126 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 1404 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 36 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   2 nodes to be scored    (41 eliminated genes)

    ## 
    ##   Level 7:   4 nodes to be scored    (45 eliminated genes)

    ## 
    ##   Level 6:   5 nodes to be scored    (124 eliminated genes)

    ## 
    ##   Level 5:   6 nodes to be scored    (274 eliminated genes)

    ## 
    ##   Level 4:   5 nodes to be scored    (362 eliminated genes)

    ## 
    ##   Level 3:   5 nodes to be scored    (528 eliminated genes)

    ## 
    ##   Level 2:   4 nodes to be scored    (645 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (768 eliminated genes)

    ## 
    ## Building most specific GOs .....

    ##  ( 458 GO terms found. )

    ## 
    ## Build GO DAG topology ..........

    ##  ( 912 GO terms and 1199 relations. )

    ## 
    ## Annotating nodes ...............

    ##  ( 2773 genes annotated to the GO terms. )

    ## 
    ##           -- Weight01 Algorithm -- 
    ## 
    ##       the algorithm is scoring 47 nontrivial nodes
    ##       parameters: 
    ##           test statistic: fisher

    ## 
    ##   Level 10:  1 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 9:   2 nodes to be scored    (0 eliminated genes)

    ## 
    ##   Level 8:   3 nodes to be scored    (7 eliminated genes)

    ## 
    ##   Level 7:   6 nodes to be scored    (312 eliminated genes)

    ## 
    ##   Level 6:   8 nodes to be scored    (334 eliminated genes)

    ## 
    ##   Level 5:   8 nodes to be scored    (392 eliminated genes)

    ## 
    ##   Level 4:   8 nodes to be scored    (424 eliminated genes)

    ## 
    ##   Level 3:   7 nodes to be scored    (590 eliminated genes)

    ## 
    ##   Level 2:   3 nodes to be scored    (812 eliminated genes)

    ## 
    ##   Level 1:   1 nodes to be scored    (1969 eliminated genes)

    ##        GO.ID                                   Term Annotated Significant
    ## 1 GO:0006400                      tRNA modification         4           1
    ## 2 GO:0005245 voltage-gated calcium channel activity         7           1
    ## 3 GO:0004177                aminopeptidase activity         7           1
    ##   Expected Fisher               type
    ## 1     0.01 0.0057 Biological.Process
    ## 2     0.01 0.0130 Molecular.Function
    ## 3     0.01 0.0130 Molecular.Function

``` r
head(results_high_cor_targets)
```

    ##        GO.ID                                               Term Annotated
    ## 1 GO:0006598                        polyamine catabolic process         3
    ## 2 GO:0001745                         compound eye morphogenesis         3
    ## 3 GO:0000014 single-stranded DNA endodeoxyribonuclease activity         3
    ## 4 GO:0005290     L-histidine transmembrane transporter activity         9
    ## 5 GO:0001658   branching involved in ureteric bud morphogenesis         4
    ## 6 GO:0001945                           lymph vessel development         5
    ##   Significant Expected Fisher               type         miRNA
    ## 1           1     0.00 0.0043 Biological.Process Cluster_10452
    ## 2           1     0.00 0.0043 Biological.Process Cluster_10452
    ## 3           1     0.01 0.0065 Molecular.Function Cluster_10452
    ## 4           1     0.02 0.0193 Molecular.Function Cluster_10452
    ## 5           1     0.01 0.0057 Biological.Process Cluster_11565
    ## 6           1     0.01 0.0071 Biological.Process Cluster_11565

Save results

``` r
write.csv(results_high_cor_targets, "../output/27-Apul-mRNA-miRNA-interactions-topGO/miRNA_high_cor_targets_topGO_FE.csv")
```
