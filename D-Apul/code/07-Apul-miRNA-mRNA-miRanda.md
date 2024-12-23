07-Apul-miRNA-mRNA-miRanda
================
Kathleen Durkin
2024-12-19

miRanda is a target prediction software, used to identify likely
miRNA-mRNA interactions.

Inputs:

- FASTA of A.pulchra 3’UTRs `Apul_3UTR_1kb.fasta`, generated in
  `05-Apul-annotate-UTRs`

- FASTA of A.pulchra mature miRNAs `miRNA_mature-Apul.fasta`. miRNAs
  identified in `04-Apul-sRNA-discovery-ShortStack`, matures isolated
  for use in `06-Apul-miRNA-mRNA-RNAhybrid`

Outputs:

``` bash

# score cutoff >100
# energy cutoff <-10
# strict binding

/home/shared/miRanda-3.3a/src/miranda \
../data/06-Apul-miRNA-mRNA-RNAhybrid/miRNA_mature-Apul.fasta \
../output/05-Apul-annotate-UTRs/Apul_3UTR_1kb.fasta \
-sc 100 \
-en -10 \
-strict \
-out ../output/07-Apul-miRNA-mRNA-miRanda/Apul-miRanda-3UTR-strict_all.tab
```

Let’s look at the output

``` bash

echo "miranda run finished!"
echo "counting number of putative interactions predicted"

zgrep -c "Performing Scan" ../output/07-Apul-miRNA-mRNA-miRanda/Apul-miRanda-3UTR-strict_all.tab

echo "Parsing output"
grep -A 1 "Scores for this hit:" ../output/07-Apul-miRNA-mRNA-miRanda/Apul-miRanda-3UTR-strict_all.tab | sort | grep '>' > ../output/07-Apul-miRNA-mRNA-miRanda/Apul-miRanda-3UTR-strict-parsed.txt

echo "counting number of putative interactions predicted"
wc -l ../output/07-Apul-miRNA-mRNA-miRanda/Apul-miRanda-3UTR-strict_all.tab
```

    ## miranda run finished!
    ## counting number of putative interactions predicted
    ## 1905309
    ## Parsing output
    ## counting number of putative interactions predicted
    ## 16835258 ../output/07-Apul-miRNA-mRNA-miRanda/Apul-miRanda-3UTR-strict_all.tab
