#!/usr/bin/env bash
# 58.0-colony-genotypes.sh
#
# Call SNP genotypes for every RNA-seq sample of one species, reading only the
# reads in a set of moderately expressed genes straight from the BAMs on
# gannet (htslib range requests), so no BAM is downloaded in full.
#
# Each sample is its own VCF sample, so the four timepoints of a colony can be
# checked for genetic identity before they are pooled into colony genotypes
# (done in 58-genotype-methylation.Rmd).
#
# usage: 58.0-colony-genotypes.sh <species: Apul|Peve|Ptua> <reference.fa> <out_dir> [parallel_jobs]
#   reference.fa must be the genome the BAMs were aligned to, with a .fai index.
#   Needs: bcftools, samtools, curl, Rscript (data.table), run from M-multi-species/scripts.
set -euo pipefail

species=$1; ref=$2; out=$3; chunks=${4:-6}
mkdir -p "$out"
G=https://gannet.fish.washington.edu

case $species in
  Apul) prefix=ACR; gff=../../D-Apul/data/Apulchra-genome.gff
        counts=../../D-Apul/output/02.20-D-Apul-RNAseq-alignment-HiSat2/apul-gene_count_matrix.csv
        bam_url() { echo "$G/v1_web/owlshell/bu-github/timeseries_molecular/D-Apul/output/31.5-Apul-lncRNA-discovery/$1.sorted.bam"; } ;;
  Peve) prefix=POR; gff=../../E-Peve/data/Porites_evermanni_validated.gff3
        counts=../../E-Peve/output/02.20-E-Peve-RNAseq-alignment-HiSat2/peve-gene_count_matrix.csv
        bam_url() { echo "$G/gitrepos/urol-e5/timeseries_molecular/E-Peve/output/02.20-E-Peve-RNAseq-alignment-HiSat2/$1/$1.sorted.bam"; } ;;
  Ptua) prefix=POC; gff=../../F-Ptua/data/Pocillopora_meandrina_HIv1.genes-validated.gff3
        counts=../../F-Ptua/output/02.20-F-Ptua-RNAseq-alignment-HiSat2/ptua-gene_count_matrix.csv
        bam_url() { echo "$G/gitrepos/urol-e5/timeseries_molecular/F-Ptua/output/02.20-F-Ptua-RNAseq-alignment-HiSat2/$1/$1.sorted.bam"; } ;;
  *) echo "unknown species $species" >&2; exit 1 ;;
esac

# 1. Target genes: ranks 501-3500 by mean expression. The most highly
#    expressed genes are skipped because they dominate read volume (and so
#    transfer time) without adding many more variable sites.
Rscript -e '
  suppressMessages(library(data.table))
  a <- commandArgs(TRUE)
  m <- fread(a[1]); x <- as.matrix(m[, -1]); x <- t(t(x) / colSums(x)) * 1e6
  ids <- m[[1]][order(rowMeans(x), decreasing = TRUE)][501:3500]
  g <- fread(cmd = paste("grep -v \"^#\"", a[2]), header = FALSE, sep = "\t", quote = "")[V3 == "gene"]
  g[, id := sub(";.*", "", sub(".*ID=", "", V9))]
  g <- g[id %in% ids]
  fwrite(g[order(V1, V4), .(V1, V4 - 1, V5)], a[3], sep = "\t", col.names = FALSE)
  cat(nrow(g), "of", length(ids), "target genes found in GFF\n")
' "$counts" "$gff" "$out/$species-targets.bed"

# 2. BAM list from the master metadata (one line per sample with RNA-seq).
awk -F, -v p="$prefix" 'NR > 1 && $1 ~ "^"p"-" {print $1"-"$2}' ../../master-molecular-metadata.csv |
  sort -u > "$out/$species-samples.txt"
: > "$out/$species-bams.txt"; : > "$out/$species-sample-names.txt"
while read -r s; do
  u=$(bam_url "$s")
  if curl -sfI --max-time 60 "$u" > /dev/null && curl -sfI --max-time 60 "$u.bai" > /dev/null; then
    echo "$u" >> "$out/$species-bams.txt"
    printf '*\t%s\t%s\n' "$u" "$s" >> "$out/$species-sample-names.txt"
  else
    echo "missing BAM or index: $s" >&2
  fi
done < "$out/$species-samples.txt"
echo "$(wc -l < "$out/$species-bams.txt") BAMs found"

# 3. Call variants in chunks of 100 target genes, `chunks` at a time, then
#    concatenate. A remote read can fail mid-stream and leave a truncated VCF,
#    so each chunk is retried up to 3 times and marked .done only on success.
#    Rerunning the script skips chunks already done.
# Region files without a .bed suffix are read as 1-based, so convert first.
awk -v OFS='\t' '{print $1, $2 + 1, $3}' "$out/$species-targets.bed" > "$out/$species-targets.regions"
mkdir -p "$out/$species-chunks"
split -a 3 -l 100 "$out/$species-targets.regions" "$out/$species-chunks/c-"
call_chunk() {
  local c=$1 try
  [[ -e $c.done ]] && return 0
  for try in 1 2 3; do
    if bcftools mpileup -f "$ref" -R "$c" -b "$out/$species-bams.txt" -G "$out/$species-sample-names.txt" \
         -a FORMAT/AD,FORMAT/DP -q 20 -Q 20 -d 500 --skip-indels -Ou 2> "$c.mpileup.log" |
       bcftools call -mv -Oz -o "$c.vcf.gz" 2> "$c.call.log"; then
      touch "$c.done"; echo "done $(basename "$c")"; return 0
    fi
    echo "retry $try failed: $(basename "$c")" >&2; sleep 30
  done
  return 1
}
export -f call_chunk; export out species ref
ls "$out/$species-chunks/"c-??? | xargs -P "$chunks" -I{} bash -c 'set -o pipefail; call_chunk {}' || true
missing=$(for c in "$out/$species-chunks/"c-???; do [[ -e $c.done ]] || echo "$c"; done | wc -l)
if (( missing > 0 )); then
  echo "$missing chunk(s) failed after retries; rerun this script to resume" >&2
  exit 1
fi
ls "$out/$species-chunks/"c-???.vcf.gz > "$out/$species-chunks.txt"
bcftools concat -f "$out/$species-chunks.txt" -Ou |
  bcftools sort -Oz -o "$out/$species-raw.vcf.gz"
bcftools index -t "$out/$species-raw.vcf.gz"

# 4. Biallelic SNPs, QUAL >= 30.
bcftools view -m2 -M2 -v snps -i 'QUAL>=30' "$out/$species-raw.vcf.gz" -Oz -o "$out/$species-snps.vcf.gz"
bcftools index -t "$out/$species-snps.vcf.gz"
echo "$species: $(bcftools view -H "$out/$species-snps.vcf.gz" | wc -l) biallelic SNPs (QUAL >= 30)"
rm -rf "$out/$species-chunks" "$out/$species-chunks.txt"
