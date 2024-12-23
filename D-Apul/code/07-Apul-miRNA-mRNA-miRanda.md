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

This is a lot of putative interactions! We can probably narrow it down
though. In vertebrates, miRNA-mRNA binding only requires complementarity
of an miRNA seed region of ~8 nucleotides. This requirement is built in
to miRanda target prediction. In cnidarians, however, miRNA-mRNA binding
is believed to require near-complete complementarity of the full mature
miRNA, similarly to plants. Let’s look at how many putative interactions
are predicted for a binding length of at least 21 nucleotides (the
length of our smallest mature miRNA).

``` bash
echo "number of putative interactions of at least 21 nucleotides"
awk -F'\t' '$7 >= 21' ../output/07-Apul-miRNA-mRNA-miRanda/Apul-miRanda-3UTR-strict-parsed.txt | wc -l
echo ""
echo "check some:"
awk -F'\t' '$7 >= 21' ../output/07-Apul-miRNA-mRNA-miRanda/Apul-miRanda-3UTR-strict-parsed.txt | head -5
```

    ## number of putative interactions of at least 21 nucleotides
    ## 18420
    ## 
    ## check some:
    ## >Cluster_10452.mature::ptg000020l:10483758-10483779(-)   ntLink_4:203294-204294  159.00  -17.21  2 21    826 849 21  66.67%  80.95%
    ## >Cluster_10452.mature::ptg000020l:10483758-10483779(-)   ntLink_4:241021-242021  159.00  -17.21  2 21    396 419 21  66.67%  80.95%
    ## >Cluster_10452.mature::ptg000020l:10483758-10483779(-)   ntLink_6:11395524-11396524  151.00  -13.49  2 21    890 914 22  63.64%  72.73%
    ## >Cluster_10452.mature::ptg000020l:10483758-10483779(-)   ntLink_6:12015318-12016318  150.00  -12.28  2 20    975 999 21  61.90%  76.19%
    ## >Cluster_10452.mature::ptg000020l:10483758-10483779(-)   ntLink_6:12240644-12241644  154.00  -18.86  2 21    361 387 24  70.83%  70.83%

We can also see from the alignment percentages (last 2 entries) that
this number includes alignments with multiple mismatches. Let’s filter
again to reduce the number of permissible mismatches. Let’s say we want
no more than 3 mismatches. For an alignment of 21 nucleotides, this
would be an alignment rate of (21-3)/21 = 85.7%.

``` bash
echo "number of putative interactions of at least 21 nucleotides, with at most 3 mismatches"
awk -F'\t' '$7 >= 21' ../output/07-Apul-miRNA-mRNA-miRanda/Apul-miRanda-3UTR-strict-parsed.txt | awk -F'\t' '$8 >= 85' | wc -l
echo ""
echo "check some:"
awk -F'\t' '$7 >= 21' ../output/07-Apul-miRNA-mRNA-miRanda/Apul-miRanda-3UTR-strict-parsed.txt | awk -F'\t' '$8 >= 85' | head -5
```

    ## number of putative interactions of at least 21 nucleotides, with at most 3 mismatches
    ## 33
    ## 
    ## check some:
    ## >Cluster_14532.mature::ptg000025l:7472581-7472603(-) ptg000007l:4113238-4114238  174.00  -19.59  2 22    352 374 21  85.71%  85.71%
    ## >Cluster_14532.mature::ptg000025l:7472581-7472603(-) ptg000016l:7511874-7512874  180.00  -21.41  2 22    573 596 21  85.71%  85.71%
    ## >Cluster_14610.mature::ptg000025l:10668923-10668945(-)   ptg000016l:1601190-1602190  180.00  -19.49  2 22    785 808 21  85.71%  85.71%
    ## >Cluster_14610.mature::ptg000025l:10668923-10668945(-)   ptg000021l:2346486-2347486  180.00  -16.96  2 22    838 861 21  85.71%  85.71%
    ## >Cluster_14633.mature::ptg000025l:11346116-11346137(+)   ptg000004l:10107854-10108854    178.00  -17.36  2 21    54 77   21  85.71%  90.48%
