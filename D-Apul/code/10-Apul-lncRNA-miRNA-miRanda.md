10-Apul-lncRNA-miRNA-miRanda
================
Kathleen Durkin
2024-12-19

- [1 Obtain lncRNA fasta](#1-obtain-lncrna-fasta)
- [2 Run miRanda](#2-run-miranda)

miRanda is a target prediction software, used to identify likely
miRNA-mRNA interactions.

Inputs:

- FASTA of A.pulchra lncRNAs, generated in `08-Apul-lncRNA`

- FASTA of A.pulchra mature miRNAs `miRNA_mature-Apul.fasta`. miRNAs
  identified in `04-Apul-sRNA-discovery-ShortStack`, matures isolated
  for use in `06-Apul-miRNA-mRNA-RNAhybrid`

# 1 Obtain lncRNA fasta

``` bash

# from lncRNA gtf
/home/shared/bedtools2/bin/bedtools getfasta -fi ../data/Apulchra-genome.fa -bed ../output/08-Apul-lncRNA/lncRNAs.gtf -fo ../output/10-Apul-lncRNA-miRNA-miRanda/lncRNAs.fa
```

# 2 Run miRanda

``` bash

# score cutoff >100
# energy cutoff <-10
# strict binding

/home/shared/miRanda-3.3a/src/miranda \
../data/06-Apul-miRNA-mRNA-RNAhybrid/miRNA_mature-Apul.fasta \
../output/10-Apul-lncRNA-miRNA-miRanda/lncRNAs.fa \
-sc 100 \
-en -10 \
-strict \
-out ../output/10-Apul-lncRNA-miRNA-miRanda/Apul-miRanda-lncRNA-strict_all.tab
```

Let’s look at the output

``` bash

echo "miranda run finished!"
echo "Counting number of interacting miRNA-lncRNA pairs"

zgrep -c "Performing Scan" ../output/10-Apul-lncRNA-miRNA-miRanda/Apul-miRanda-lncRNA-strict_all.tab

echo "Parsing output"
grep -A 1 "Scores for this hit:" ../output/10-Apul-lncRNA-miRNA-miRanda/Apul-miRanda-lncRNA-strict_all.tab | sort | grep '>' > ../output/10-Apul-lncRNA-miRNA-miRanda/Apul-miRanda-lncRNA-strict-parsed.txt

echo "counting number of putative interactions predicted (can include multiple interactions between single miRNA-lncRNA pair)"
wc -l ../output/10-Apul-lncRNA-miRNA-miRanda/Apul-miRanda-lncRNA-strict_all.tab
```

    ## miranda run finished!
    ## Counting number of interacting miRNA-lncRNA pairs
    ## 1233231
    ## Parsing output
    ## counting number of putative interactions predicted (can include multiple interactions between single miRNA-lncRNA pair)
    ## 12105116 ../output/10-Apul-lncRNA-miRNA-miRanda/Apul-miRanda-lncRNA-strict_all.tab

This is a lot of putative interactions! We can probably narrow it down
though. In vertebrates, miRNA-mRNA binding only requires complementarity
of an miRNA seed region of ~8 nucleotides. This requirement is built in
to miRanda target prediction. In cnidarians, however, miRNA-mRNA binding
is believed to require near-complete complementarity of the full mature
miRNA, similarly to plants. While I couldn’t find any information on
expected requirements for miRNA-lncRNA sponges, its possible the binding
will function similarly to miRNA-mRNa binding. Let’s look at how many
putative interactions are predicted for a binding length of at least 21
nucleotides (the length of our smallest mature miRNA).

``` bash
echo "number of putative interactions of at least 21 nucleotides"
awk -F'\t' '$7 >= 21' ../output/10-Apul-lncRNA-miRNA-miRanda/Apul-miRanda-lncRNA-strict-parsed.txt | wc -l
echo ""
echo "check some:"
awk -F'\t' '$7 >= 21' ../output/10-Apul-lncRNA-miRNA-miRanda/Apul-miRanda-lncRNA-strict-parsed.txt | head -5
```

    ## number of putative interactions of at least 21 nucleotides
    ## 24018
    ## 
    ## check some:
    ## >Cluster_10452.mature::ptg000020l:10483758-10483779(-)   ntLink_3:56720-62176    150.00  -18.45  2 21    351 375 22  63.64%  77.27%
    ## >Cluster_10452.mature::ptg000020l:10483758-10483779(-)   ntLink_3:76351-81676    150.00  -18.45  2 21    220 244 22  63.64%  77.27%
    ## >Cluster_10452.mature::ptg000020l:10483758-10483779(-)   ntLink_3:95704-100623   150.00  -18.45  2 21    378 402 22  63.64%  77.27%
    ## >Cluster_10452.mature::ptg000020l:10483758-10483779(-)   ntLink_4:152868-157343  159.00  -17.21  2 21    2472 2495   21  66.67%  80.95%
    ## >Cluster_10452.mature::ptg000020l:10483758-10483779(-)   ntLink_4:205718-210236  159.00  -17.21  2 21    2451 2474   21  66.67%  80.95%

We can also see from the alignment percentages (last 2 entries) that
this number includes alignments with multiple mismatches. Let’s filter
again to reduce the number of permissible mismatches. Let’s say we want
no more than 3 mismatches. For an alignment of 21 nucleotides, this
would be an alignment rate of (21-3)/21 = 85.7%.

``` bash
echo "number of putative interactions of at least 21 nucleotides, with at most 3 mismatches"
awk -F'\t' '$7 >= 21' ../output/10-Apul-lncRNA-miRNA-miRanda/Apul-miRanda-lncRNA-strict-parsed.txt | awk -F'\t' '$8 >= 85' | wc -l
echo ""
echo "check some:"
awk -F'\t' '$7 >= 21' ../output/10-Apul-lncRNA-miRNA-miRanda/Apul-miRanda-lncRNA-strict-parsed.txt | awk -F'\t' '$8 >= 85' | head -5
```

    ## number of putative interactions of at least 21 nucleotides, with at most 3 mismatches
    ## 22
    ## 
    ## check some:
    ## >Cluster_14532.mature::ptg000025l:7472581-7472603(-) ptg000021l:14212189-14212865    174.00  -23.17  2 22    384 406 21  85.71%  85.71%
    ## >Cluster_14532.mature::ptg000025l:7472581-7472603(-) ptg000023l:17245108-17245424    174.00  -18.76  2 22    49 71   21  85.71%  85.71%
    ## >Cluster_14532.mature::ptg000025l:7472581-7472603(-) ptg000024l:8680062-8686568  178.00  -21.90  2 22    3861 3883   21  85.71%  90.48%
    ## >Cluster_14532.mature::ptg000025l:7472581-7472603(-) ptg000024l:8680105-8686841  178.00  -21.90  2 22    3818 3840   21  85.71%  90.48%
    ## >Cluster_14532.mature::ptg000025l:7472581-7472603(-) ptg000026l:4000983-4003199  180.00  -21.97  2 22    185 208 21  85.71%  85.71%
