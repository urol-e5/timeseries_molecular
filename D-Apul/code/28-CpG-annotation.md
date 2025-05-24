28-Apul CpG Motifs
================
Steven Roberts
14 May, 2025

- <a href="#1-genome" id="toc-1-genome">1 genome</a>
- <a href="#2-full-run" id="toc-2-full-run">2 Full run</a>
- <a href="#3-intersectbed" id="toc-3-intersectbed">3 Intersectbed</a>

Annotating CpGs ….

methylation calls derived at
<https://sr320.github.io/tumbling-oysters/posts/37-Apul-meth/>

Now there are files for each of the Apul samples

@
<https://gannet.fish.washington.edu/seashell/bu-github/timeseries_molecular/D-Apul/output/15.5-Apul-bismark/>

Maybe do all CpG in genome..

# 1 genome

``` bash
head ../data/Apulcra-genome.fa
```

    >ntLink_7
    TGTATTTCTAGAGATCTAAAGTACTCAGGATAGCTGAAAAAAAAACCGTC
    GTTCTTTCCCTCTGAATTTCTCACTTAGGTTTTCTGGGATATGTGACGCG
    CTGCAATCCTAAACAGCCACTAACGCCGGATGTCTTCAACCGAGCGCCTC
    TGCTCGAAAGCTATCTGGACCATTTGAAGGTATGGTTCAGCTTGATTTAC
    TGTTGCTTCGTAGTACTTTCATCCTGATGGTTCTAATGTGCTAAGTTTTC
    GCACATATTCTAGGTCATTAGGTTCTAATAGCAATGAGCTTTCAGTGATT
    CCCGAGAATACGCTTTCTGTAGCAGCAACTAAATTATGGAATGGTCTACC
    CTAAAATCTTAGAAAGTCTACCTCGTATATCAGCAGGCTAACTGCTACGT
    TCATTTTTCAAGCATTTGAAAGACGAAAGAAACAATAATAGTTAACCATA

``` r
library(seqinr)

# Replace 'input.fasta' with the name of your multi-sequence fasta file
input_file <- "../data/Apulcra-genome.fa"
sequences <- read.fasta(input_file)
```

``` r
# Set the seed for reproducibility (optional)
set.seed(42)

number_of_sequences_to_select <- 10

if (length(sequences) < number_of_sequences_to_select) {
  warning("There are fewer than 10 sequences in the fasta file. All sequences will be selected.")
  number_of_sequences_to_select <- length(sequences)
}

selected_indices <- sample(length(sequences), number_of_sequences_to_select)
selected_sequences <- sequences[selected_indices]
```

``` r
# Replace 'output.fasta' with your desired output file name
output_file <- "../output/28-Apul-CpG-Annotation/output.fasta"
write.fasta(selected_sequences, names(selected_sequences), output_file, open = "w")
```

``` bash
#likely will not need; fix issue where gff and fa name did not match
# sed -i 's/>lcl|/>/g' ../output/10_seqs.fa
```

``` bash
#needed downstream for IGV
/home/shared/samtools-1.12/samtools faidx \
../output/28-Apul-CpG-Annotation/output.fasta
```

``` bash
fuzznuc -sequence ../output/28-Apul-CpG-Annotation/output.fasta -pattern CG -rformat gff -outfile ../output/28-Apul-CpG-Annotation/CGoutput-10seq.gff
```

``` bash
tail ../output/28-Apul-CpG-Annotation/CGoutput-10seq.gff
```

# 2 Full run

``` bash
fuzznuc -sequence ../data/Apulcra-genome.fa -pattern CG -rformat gff -outfile ../output/28-Apul-CpG-Annotation/Apul-CG-motifs.gff
```

# 3 Intersectbed

``` bash
bedtools intersect \
-a ../output/28-Apul-CpG-Annotation/Apul-CG-motifs.gff \
-b ../data/Apulchra-genome.gff \
-wb \
> ../output/28-Apul-CpG-Annotation/intersect_both.gff
```

``` bash
head ../output/28-Apul-CpG-Annotation/intersect_both.gff
```

``` bash
tail -50 ../data/Apulchra-genome.gff
```

    ptg000184l  funannotate exon    27291   27416   .   +   .   ID=FUN_044364-T1.exon2;Parent=FUN_044364-T1;
    ptg000184l  funannotate exon    27559   27784   .   +   .   ID=FUN_044364-T1.exon3;Parent=FUN_044364-T1;
    ptg000184l  funannotate exon    28179   28245   .   +   .   ID=FUN_044364-T1.exon4;Parent=FUN_044364-T1;
    ptg000184l  funannotate exon    33400   35016   .   +   .   ID=FUN_044364-T1.exon5;Parent=FUN_044364-T1;
    ptg000184l  funannotate CDS 26899   27163   .   +   0   ID=FUN_044364-T1.cds;Parent=FUN_044364-T1;
    ptg000184l  funannotate CDS 27291   27416   .   +   2   ID=FUN_044364-T1.cds;Parent=FUN_044364-T1;
    ptg000184l  funannotate CDS 27559   27784   .   +   2   ID=FUN_044364-T1.cds;Parent=FUN_044364-T1;
    ptg000184l  funannotate CDS 28179   28245   .   +   1   ID=FUN_044364-T1.cds;Parent=FUN_044364-T1;
    ptg000184l  funannotate CDS 33400   35016   .   +   0   ID=FUN_044364-T1.cds;Parent=FUN_044364-T1;
    ptg000185l  funannotate gene    8952    9296    .   -   .   ID=FUN_044365;
    ptg000185l  funannotate mRNA    8952    9296    .   -   .   ID=FUN_044365-T1;Parent=FUN_044365;product=hypothetical protein;
    ptg000185l  funannotate exon    8952    9296    .   -   .   ID=FUN_044365-T1.exon1;Parent=FUN_044365-T1;
    ptg000185l  funannotate CDS 8952    9296    .   -   0   ID=FUN_044365-T1.cds;Parent=FUN_044365-T1;
    ptg000185l  funannotate gene    30137   34494   .   +   .   ID=FUN_044366;
    ptg000185l  funannotate mRNA    30137   34494   .   +   .   ID=FUN_044366-T1;Parent=FUN_044366;product=hypothetical protein;
    ptg000185l  funannotate exon    30137   30432   .   +   .   ID=FUN_044366-T1.exon1;Parent=FUN_044366-T1;
    ptg000185l  funannotate exon    30536   30737   .   +   .   ID=FUN_044366-T1.exon2;Parent=FUN_044366-T1;
    ptg000185l  funannotate exon    30865   31414   .   +   .   ID=FUN_044366-T1.exon3;Parent=FUN_044366-T1;
    ptg000185l  funannotate exon    33568   34163   .   +   .   ID=FUN_044366-T1.exon4;Parent=FUN_044366-T1;
    ptg000185l  funannotate exon    34251   34494   .   +   .   ID=FUN_044366-T1.exon5;Parent=FUN_044366-T1;
    ptg000185l  funannotate CDS 30137   30432   .   +   0   ID=FUN_044366-T1.cds;Parent=FUN_044366-T1;
    ptg000185l  funannotate CDS 30536   30737   .   +   1   ID=FUN_044366-T1.cds;Parent=FUN_044366-T1;
    ptg000185l  funannotate CDS 30865   31414   .   +   0   ID=FUN_044366-T1.cds;Parent=FUN_044366-T1;
    ptg000185l  funannotate CDS 33568   34163   .   +   2   ID=FUN_044366-T1.cds;Parent=FUN_044366-T1;
    ptg000185l  funannotate CDS 34251   34494   .   +   0   ID=FUN_044366-T1.cds;Parent=FUN_044366-T1;
    ptg000186l  funannotate gene    587 1375    .   +   .   ID=FUN_044367;
    ptg000186l  funannotate mRNA    587 1375    .   +   .   ID=FUN_044367-T1;Parent=FUN_044367;product=hypothetical protein;
    ptg000186l  funannotate exon    587 835 .   +   .   ID=FUN_044367-T1.exon1;Parent=FUN_044367-T1;
    ptg000186l  funannotate exon    1070    1375    .   +   .   ID=FUN_044367-T1.exon2;Parent=FUN_044367-T1;
    ptg000186l  funannotate CDS 587 835 .   +   0   ID=FUN_044367-T1.cds;Parent=FUN_044367-T1;
    ptg000186l  funannotate CDS 1070    1375    .   +   0   ID=FUN_044367-T1.cds;Parent=FUN_044367-T1;
    ptg000186l  funannotate gene    2304    2977    .   +   .   ID=FUN_044368;
    ptg000186l  funannotate mRNA    2304    2977    .   +   .   ID=FUN_044368-T1;Parent=FUN_044368;product=hypothetical protein;
    ptg000186l  funannotate exon    2304    2504    .   +   .   ID=FUN_044368-T1.exon1;Parent=FUN_044368-T1;
    ptg000186l  funannotate exon    2576    2977    .   +   .   ID=FUN_044368-T1.exon2;Parent=FUN_044368-T1;
    ptg000186l  funannotate CDS 2304    2504    .   +   0   ID=FUN_044368-T1.cds;Parent=FUN_044368-T1;
    ptg000186l  funannotate CDS 2576    2977    .   +   0   ID=FUN_044368-T1.cds;Parent=FUN_044368-T1;
    ptg000186l  funannotate gene    8523    11277   .   +   .   ID=FUN_044369;
    ptg000186l  funannotate mRNA    8523    11277   .   +   .   ID=FUN_044369-T1;Parent=FUN_044369;product=hypothetical protein;
    ptg000186l  funannotate exon    8523    8549    .   +   .   ID=FUN_044369-T1.exon1;Parent=FUN_044369-T1;
    ptg000186l  funannotate exon    10894   11277   .   +   .   ID=FUN_044369-T1.exon2;Parent=FUN_044369-T1;
    ptg000186l  funannotate CDS 8523    8549    .   +   0   ID=FUN_044369-T1.cds;Parent=FUN_044369-T1;
    ptg000186l  funannotate CDS 10894   11277   .   +   0   ID=FUN_044369-T1.cds;Parent=FUN_044369-T1;
    ptg000186l  funannotate gene    18980   19053   .   -   .   ID=FUN_044370;
    ptg000186l  funannotate tRNA    18980   19053   .   -   .   ID=FUN_044370-T1;Parent=FUN_044370;product=tRNA-Thr;
    ptg000186l  funannotate exon    18980   19053   .   -   .   ID=FUN_044370-T1.exon1;Parent=FUN_044370-T1;
    ptg000187l  funannotate gene    16329   17045   .   -   .   ID=FUN_044371;
    ptg000187l  funannotate mRNA    16329   17045   .   -   .   ID=FUN_044371-T1;Parent=FUN_044371;product=hypothetical protein;
    ptg000187l  funannotate exon    16329   17045   .   -   .   ID=FUN_044371-T1.exon1;Parent=FUN_044371-T1;
    ptg000187l  funannotate CDS 16329   17045   .   -   0   ID=FUN_044371-T1.cds;Parent=FUN_044371-T1;

``` bash
awk -F'\t' '
BEGIN {
  header = "##gff-version 3"
}
$0 ~ /^#/ { next }
{
  outfile = "../output/28-Apul-CpG-Annotation/" $3 ".gff"
  if (!(outfile in written)) {
    print header > outfile
    written[outfile] = 1
  }
  print >> outfile
}' ../data/Apulchra-genome.gff
```

``` bash
ls ../output/28-Apul-CpG-Annotation/*gff
```

    ../output/28-Apul-CpG-Annotation/Apul-CG-motifs.gff
    ../output/28-Apul-CpG-Annotation/CDS.gff
    ../output/28-Apul-CpG-Annotation/CG_intersect_mRNA.gff
    ../output/28-Apul-CpG-Annotation/CGoutput-10seq.gff
    ../output/28-Apul-CpG-Annotation/exon.gff
    ../output/28-Apul-CpG-Annotation/gene.gff
    ../output/28-Apul-CpG-Annotation/intersect_both.gff
    ../output/28-Apul-CpG-Annotation/mRNA.gff
    ../output/28-Apul-CpG-Annotation/tRNA.gff

``` bash
bedtools intersect \
-a ../output/28-Apul-CpG-Annotation/Apul-CG-motifs.gff \
-b ../output/28-Apul-CpG-Annotation/mRNA.gff \
-wb \
> ../output/28-Apul-CpG-Annotation/CG_intersect_mRNA.gff
```

``` bash
head ../output/28-Apul-CpG-Annotation/CG_intersect_mRNA.gff
```

    ntLink_7    fuzznuc nucleotide_motif    97  98  2   +   .   ID=ntLink_7.3;note=*pat pattern:CG  ntLink_7    funannotate mRNA    79  4679    .   +   .   ID=FUN_002303-T1;Parent=FUN_002303;product=hypothetical protein;
    ntLink_7    fuzznuc nucleotide_motif    99  100 2   +   .   ID=ntLink_7.4;note=*pat pattern:CG  ntLink_7    funannotate mRNA    79  4679    .   +   .   ID=FUN_002303-T1;Parent=FUN_002303;product=hypothetical protein;
    ntLink_7    fuzznuc nucleotide_motif    124 125 2   +   .   ID=ntLink_7.5;note=*pat pattern:CG  ntLink_7    funannotate mRNA    79  4679    .   +   .   ID=FUN_002303-T1;Parent=FUN_002303;product=hypothetical protein;
    ntLink_7    fuzznuc nucleotide_motif    127 128 2   +   .   ID=ntLink_7.6;note=*pat pattern:CG  ntLink_7    funannotate mRNA    79  4679    .   +   .   ID=FUN_002303-T1;Parent=FUN_002303;product=hypothetical protein;
    ntLink_7    fuzznuc nucleotide_motif    141 142 2   +   .   ID=ntLink_7.7;note=*pat pattern:CG  ntLink_7    funannotate mRNA    79  4679    .   +   .   ID=FUN_002303-T1;Parent=FUN_002303;product=hypothetical protein;
    ntLink_7    fuzznuc nucleotide_motif    145 146 2   +   .   ID=ntLink_7.8;note=*pat pattern:CG  ntLink_7    funannotate mRNA    79  4679    .   +   .   ID=FUN_002303-T1;Parent=FUN_002303;product=hypothetical protein;
    ntLink_7    fuzznuc nucleotide_motif    155 156 2   +   .   ID=ntLink_7.9;note=*pat pattern:CG  ntLink_7    funannotate mRNA    79  4679    .   +   .   ID=FUN_002303-T1;Parent=FUN_002303;product=hypothetical protein;
    ntLink_7    fuzznuc nucleotide_motif    209 210 2   +   .   ID=ntLink_7.10;note=*pat pattern:CG ntLink_7    funannotate mRNA    79  4679    .   +   .   ID=FUN_002303-T1;Parent=FUN_002303;product=hypothetical protein;
    ntLink_7    fuzznuc nucleotide_motif    250 251 2   +   .   ID=ntLink_7.11;note=*pat pattern:CG ntLink_7    funannotate mRNA    79  4679    .   +   .   ID=FUN_002303-T1;Parent=FUN_002303;product=hypothetical protein;
    ntLink_7    fuzznuc nucleotide_motif    303 304 2   +   .   ID=ntLink_7.12;note=*pat pattern:CG ntLink_7    funannotate mRNA    79  4679    .   +   .   ID=FUN_002303-T1;Parent=FUN_002303;product=hypothetical protein;
