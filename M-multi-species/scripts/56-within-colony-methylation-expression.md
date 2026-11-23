56-within-colony-methylation-expression
================
Steven Roberts
2026-11-23

- [1 Goal](#1-goal)
- [2 Inputs](#2-inputs)
- [3 Analysis](#3-analysis)
  - [3.1 1. Across genes](#31-1-across-genes)
  - [3.2 2. Between colonies](#32-2-between-colonies)
  - [3.3 3. Within colonies](#33-3-within-colonies)
  - [3.4 Strongest coupled genes](#34-strongest-coupled-genes)

# 1 Goal

`55-variance-partitioning` showed that gene-body methylation varies
mostly among colonies and very little with shared time, while expression
tracks time much more. That leaves open whether methylation and
expression **change together within a colony** over time. This script
tests that directly for every gene that has both gene-body methylation
and expression.

Three levels are examined:

1.  **Across genes** - are more methylated genes more highly expressed,
    and less variable over time? (the classic gene-body methylation
    pattern)
2.  **Between colonies** - for a given gene, do colonies with higher
    methylation have higher expression? (10 colonies per gene, so low
    power per gene; the distribution across genes is what matters)
3.  **Within colonies** - when a colony’s methylation of a gene goes up
    from one timepoint to another, does its expression of that gene move
    with it? Tested two ways:
    - **within colony**: both layers centred on each colony’s mean, so
      this includes change shared by all colonies (season)
    - **within colony, net of shared time**: site-by-timepoint means
      also removed, leaving only colony-specific change. This is the
      stricter test of gene-specific coupling.

The null for the within-colony tests comes from shuffling methylation
samples among timepoints within each colony, the same shuffle for every
gene, which keeps the correlation among genes intact.

``` r
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE,
                      fig.width = 10, fig.height = 6, dpi = 150)
library(data.table)
library(DESeq2)
library(ggplot2)
library(matrixStats)

out_dir  <- "../output/56-within-colony-methylation-expression"
meth_dir <- "../output/55-variance-partitioning"   # gene-body methylation from 55.1
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
n_perm <- 1000
set.seed(56)
```

# 2 Inputs

Gene-body methylation (mean % over at least 10 CpGs per gene and sample)
comes from `55.1-gene-body-methylation.py`, run by
`55-variance-partitioning.Rmd`. Only genes with a methylation value in
every sample shared by both layers are used, so every gene is tested on
the same samples.

``` r
meta <- fread("../../master-molecular-metadata.csv")
meta[, sample := paste(colony, timepoint, sep = "-")]
meta[, st := paste(site, timepoint, sep = ":")]
setkey(meta, sample)

canon <- function(x) {
  hit <- grepl("(ACR|POR|POC)-[0-9]+[-_]TP[1-4]", x)
  out <- rep(NA_character_, length(x))
  out[hit] <- sub("_", "-", regmatches(x, regexpr("(ACR|POR|POC)-[0-9]+[-_]TP[1-4]", x)))
  out
}
as_matrix <- function(dt) {
  keep <- names(dt)[-1]; s <- canon(keep)
  keep <- keep[!is.na(s)]; s <- s[!is.na(s)]
  m <- as.matrix(dt[, ..keep]); storage.mode(m) <- "double"
  if (anyDuplicated(s)) m <- t(rowsum(t(m), s)) else colnames(m) <- s
  rownames(m) <- dt[[1]]
  m[, sort(colnames(m)), drop = FALSE]
}
vst_counts <- function(m, min_count = 10) {
  m <- round(m)
  m <- m[rowSums(m >= min_count) >= ceiling(ncol(m) / 4), , drop = FALSE]
  dds <- DESeqDataSetFromMatrix(m, data.frame(row.names = colnames(m)), ~ 1)
  dds <- estimateSizeFactors(dds, type = "poscounts")
  assay(varianceStabilizingTransformation(dds, blind = TRUE))
}

sp <- list(
  Apul = "../../D-Apul/output/02.20-D-Apul-RNAseq-alignment-HiSat2/apul-gene_count_matrix.csv",
  Peve = "../../E-Peve/output/02.20-E-Peve-RNAseq-alignment-HiSat2/peve-gene_count_matrix.csv",
  Ptua = "../../F-Ptua/output/02.20-F-Ptua-RNAseq-alignment-HiSat2/ptua-gene_count_matrix.csv"
)
```

# 3 Analysis

``` r
# residuals of every row of Y after regressing on design X (rank-safe)
resid_rows <- function(Y, X) {
  q <- qr(X); Q <- qr.Q(q)[, seq_len(q$rank), drop = FALSE]
  Y - (Y %*% Q) %*% t(Q)
}
rowcor <- function(A, B) rowSums(A * B) / sqrt(rowSums(A^2) * rowSums(B^2))
p_from_r <- function(r, df) 2 * pt(-abs(r * sqrt(df / (1 - r^2))), df)
shuffle_within <- function(groups) {
  idx <- seq_along(groups)
  for (g in split(idx, groups)) if (length(g) > 1) idx[g] <- g[sample.int(length(g))]
  idx
}
meth_classes <- function(x) cut(x, c(-Inf, 5, 20, Inf), labels = c("<5%", "5-20%", ">20%"))
```

``` r
per_gene <- list(); perm_stats <- list(); inputs <- list()

for (s in names(sp)) {
  E <- vst_counts(as_matrix(fread(sp[[s]])))
  M <- as_matrix(fread(file.path(meth_dir, paste0(s, "-gene-body-methylation.csv.gz"))))
  smp <- intersect(colnames(E), colnames(M))
  genes <- intersect(rownames(E), rownames(M))
  M <- M[genes, smp, drop = FALSE]
  M <- M[rowSums(is.na(M)) == 0 & rowVars(M) > 0, , drop = FALSE]
  E <- E[rownames(M), smp, drop = FALSE]
  md <- meta[smp]
  colony <- factor(md$colony)
  n <- length(smp)
  inputs[[s]] <- data.table(species = s, samples = n, colonies = nlevels(colony),
                            genes = nrow(M), genes_methylated_5pct = sum(rowMeans(M) >= 5))

  designs <- list(within_colony          = model.matrix(~ colony),
                  within_colony_net_time = model.matrix(~ colony + st, md))
  cls <- meth_classes(rowMeans(M))
  res <- data.table(species = s, gene = rownames(M), mean_meth = rowMeans(M),
                    meth_class = cls, mean_expr = rowMeans(E))

  # 1. within-colony variability of expression (plasticity)
  Ewc <- resid_rows(E, designs$within_colony)
  res[, expr_within_colony_sd := sqrt(rowSums(Ewc^2) / (n - nlevels(colony)))]

  # 2. between colonies: correlation of colony means
  cm <- function(Y) t(rowsum(t(Y), colony) / as.vector(table(colony)))
  res[, r_between := rowcor(cm(M) - rowMeans(cm(M)), cm(E) - rowMeans(cm(E)))]

  # 3. within colonies, with permutation null
  for (d in names(designs)) {
    X <- designs[[d]]
    df <- n - qr(X)$rank - 1
    Er <- resid_rows(E, X)
    r <- rowcor(resid_rows(M, X), Er)
    res[, (paste0("r_", d)) := r]
    res[, (paste0("p_", d)) := p_from_r(r, df)]

    perm <- rbindlist(lapply(seq_len(n_perm), function(b) {
      rp <- rowcor(resid_rows(M[, shuffle_within(colony), drop = FALSE], X), Er)
      pp <- p_from_r(rp, df)
      rbind(data.table(meth_class = "all", mean_r = mean(rp), n_p01 = sum(pp < 0.01)),
            data.table(meth_class = as.character(cls), rp, pp)[
              , .(mean_r = mean(rp), n_p01 = sum(pp < 0.01)), by = meth_class])
    }), idcol = "perm")
    perm_stats[[paste(s, d)]] <- cbind(species = s, test = d, perm)
  }
  per_gene[[s]] <- res
}

per_gene <- rbindlist(per_gene)
for (d in c("within_colony", "within_colony_net_time"))
  per_gene[, (paste0("fdr_", d)) := p.adjust(get(paste0("p_", d)), "BH"), by = species]
perm_stats <- rbindlist(perm_stats)
fwrite(per_gene, file.path(out_dir, "per-gene-methylation-expression.csv.gz"))
knitr::kable(rbindlist(inputs), caption = "Genes and samples analysed")
```

| species | samples | colonies | genes | genes_methylated_5pct |
|:--------|--------:|---------:|------:|----------------------:|
| Apul    |      39 |       10 | 13193 |                  6954 |
| Peve    |      37 |       10 |  7414 |                  3579 |
| Ptua    |      32 |       10 | 17173 |                  6946 |

Genes and samples analysed

## 3.1 1. Across genes

Spearman correlations across genes between mean gene-body methylation
and (a) mean expression, (b) within-colony standard deviation of
expression.

``` r
across <- per_gene[, .(
  genes = .N,
  rho_meth_vs_mean_expr = cor(mean_meth, mean_expr, method = "spearman"),
  rho_meth_vs_expr_sd = cor(mean_meth, expr_within_colony_sd, method = "spearman"),
  rho_meth_vs_expr_sd_given_mean = {
    # partial: expression variability adjusted for expression level
    a <- resid(lm(rank(expr_within_colony_sd) ~ rank(mean_expr)))
    b <- resid(lm(rank(mean_meth) ~ rank(mean_expr)))
    cor(a, b)
  }), by = species]
fwrite(across, file.path(out_dir, "across-genes.csv"))
knitr::kable(across, digits = 3)
```

| species | genes | rho_meth_vs_mean_expr | rho_meth_vs_expr_sd | rho_meth_vs_expr_sd_given_mean |
|:---|---:|---:|---:|---:|
| Apul | 13193 | 0.457 | -0.167 | -0.189 |
| Peve | 7414 | 0.398 | -0.412 | -0.358 |
| Ptua | 17173 | 0.388 | -0.146 | -0.190 |

``` r
per_gene[, bin := cut(frank(mean_meth) / .N, seq(0, 1, 0.05)), by = species]
binned <- per_gene[, .(mean_expr = mean(mean_expr), expr_sd = mean(expr_within_colony_sd),
                       mean_meth = mean(mean_meth), genes = .N), by = .(species, bin)]
per_gene[, bin := NULL]
bl <- melt(binned, id.vars = c("species", "bin", "mean_meth", "genes"),
           measure.vars = c("mean_expr", "expr_sd"))
bl[, variable := factor(variable, c("mean_expr", "expr_sd"),
                        c("Mean expression (VST)", "Within-colony SD of expression"))]
p1 <- ggplot(bl, aes(mean_meth, value, colour = species)) +
  geom_line() + geom_point(size = 1.5) +
  facet_wrap(~ variable, scales = "free_y") +
  labs(x = "Gene-body methylation (%), 5% bins of genes", y = NULL, colour = NULL) +
  theme_bw()
ggsave(file.path(out_dir, "across-genes.png"), p1, width = 10, height = 4.5, dpi = 200)
p1
```

![](56-within-colony-methylation-expression_files/figure-gfm/across-plot-1.png)<!-- -->

## 3.2 2. Between colonies

``` r
between <- per_gene[, .(genes = .N, mean_r = mean(r_between),
                        frac_positive = mean(r_between > 0)), by = .(species, meth_class)]
setorder(between, species, meth_class)
knitr::kable(between, digits = 3)
```

| species | meth_class | genes | mean_r | frac_positive |
|:--------|:-----------|------:|-------:|--------------:|
| Apul    | \<5%       |  6239 |  0.011 |         0.506 |
| Apul    | 5-20%      |  2890 |  0.072 |         0.578 |
| Apul    | \>20%      |  4064 |  0.050 |         0.544 |
| Peve    | \<5%       |  3835 |  0.015 |         0.515 |
| Peve    | 5-20%      |  1364 |  0.053 |         0.557 |
| Peve    | \>20%      |  2215 |  0.027 |         0.532 |
| Ptua    | \<5%       | 10227 |  0.010 |         0.510 |
| Ptua    | 5-20%      |  4324 |  0.020 |         0.518 |
| Ptua    | \>20%      |  2622 |  0.009 |         0.503 |

## 3.3 3. Within colonies

For each test: the mean within-colony correlation across genes, the
number of genes with nominal p \< 0.01, the number expected under the
permutation null, and permutation p-values for both statistics.

``` r
obs <- rbindlist(lapply(c("within_colony", "within_colony_net_time"), function(d) {
  x <- per_gene[, .(species, meth_class = as.character(meth_class),
                    r = get(paste0("r_", d)), p = get(paste0("p_", d)),
                    fdr = get(paste0("fdr_", d)))]
  rbind(x[, .(meth_class = "all", genes = .N, mean_r = mean(r), n_p01 = sum(p < 0.01),
              n_fdr10 = sum(fdr < 0.1)), by = species],
        x[, .(genes = .N, mean_r = mean(r), n_p01 = sum(p < 0.01),
              n_fdr10 = sum(fdr < 0.1)), by = .(species, meth_class)])[, test := d]
}))
keys <- c("species", "test", "meth_class")
pp <- merge(perm_stats, obs[, .(species, test, meth_class, obs_r = mean_r, obs_n = n_p01)],
            by = keys)
null <- pp[, .(null_mean_r = mean(mean_r),
               null_n_p01 = mean(n_p01),
               # two-sided for the mean correlation, one-sided (excess) for the count
               perm_p_mean_r = (1 + sum(abs(mean_r - mean(mean_r)) >=
                                          abs(obs_r[1] - mean(mean_r)))) / (n_perm + 1),
               perm_p_n_p01 = (1 + sum(n_p01 >= obs_n[1])) / (n_perm + 1)), by = keys]
within <- merge(obs, null, by = keys)
within[, meth_class := factor(meth_class, c("all", "<5%", "5-20%", ">20%"))]
setorder(within, test, species, meth_class)
fwrite(within, file.path(out_dir, "within-colony-summary.csv"))
knitr::kable(within[, .(test, species, meth_class, genes, mean_r, null_mean_r, perm_p_mean_r,
                        n_p01, null_n_p01, perm_p_n_p01, n_fdr10)], digits = 3)
```

| test | species | meth_class | genes | mean_r | null_mean_r | perm_p_mean_r | n_p01 | null_n_p01 | perm_p_n_p01 | n_fdr10 |
|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| within_colony | Apul | all | 13193 | 0.007 | 0.000 | 0.395 | 1174 | 140.081 | 0.001 | 1039 |
| within_colony | Apul | \<5% | 6239 | -0.012 | 0.000 | 0.367 | 137 | 65.377 | 0.001 | 112 |
| within_colony | Apul | 5-20% | 2890 | 0.014 | 0.000 | 0.430 | 359 | 29.162 | 0.001 | 317 |
| within_colony | Apul | \>20% | 4064 | 0.032 | 0.001 | 0.116 | 678 | 45.542 | 0.001 | 610 |
| within_colony | Peve | all | 7414 | 0.017 | 0.000 | 0.013 | 118 | 79.425 | 0.010 | 4 |
| within_colony | Peve | \<5% | 3835 | 0.000 | 0.000 | 0.952 | 43 | 41.094 | 0.423 | 2 |
| within_colony | Peve | 5-20% | 1364 | 0.051 | 0.000 | 0.001 | 23 | 15.112 | 0.093 | 1 |
| within_colony | Peve | \>20% | 2215 | 0.026 | 0.000 | 0.034 | 52 | 23.219 | 0.002 | 1 |
| within_colony | Ptua | all | 17173 | 0.016 | 0.000 | 0.111 | 357 | 223.154 | 0.005 | 23 |
| within_colony | Ptua | \<5% | 10227 | 0.026 | 0.000 | 0.010 | 174 | 121.128 | 0.004 | 11 |
| within_colony | Ptua | 5-20% | 4324 | 0.002 | 0.001 | 0.965 | 110 | 62.308 | 0.023 | 4 |
| within_colony | Ptua | \>20% | 2622 | 0.001 | 0.001 | 1.000 | 73 | 39.718 | 0.023 | 8 |
| within_colony_net_time | Apul | all | 13193 | -0.005 | 0.000 | 0.474 | 422 | 142.589 | 0.005 | 26 |
| within_colony_net_time | Apul | \<5% | 6239 | -0.010 | 0.000 | 0.434 | 93 | 64.987 | 0.018 | 0 |
| within_colony_net_time | Apul | 5-20% | 2890 | 0.006 | 0.000 | 0.724 | 101 | 29.989 | 0.008 | 4 |
| within_colony_net_time | Apul | \>20% | 4064 | -0.006 | 0.000 | 0.748 | 228 | 47.613 | 0.005 | 22 |
| within_colony_net_time | Peve | all | 7414 | 0.006 | 0.000 | 0.423 | 85 | 76.071 | 0.250 | 4 |
| within_colony_net_time | Peve | \<5% | 3835 | -0.006 | 0.000 | 0.496 | 47 | 40.566 | 0.219 | 1 |
| within_colony_net_time | Peve | 5-20% | 1364 | 0.033 | 0.001 | 0.013 | 15 | 13.654 | 0.412 | 1 |
| within_colony_net_time | Peve | \>20% | 2215 | 0.008 | -0.001 | 0.381 | 23 | 21.851 | 0.420 | 2 |
| within_colony_net_time | Ptua | all | 17173 | 0.011 | 0.000 | 0.346 | 266 | 206.175 | 0.052 | 7 |
| within_colony_net_time | Ptua | \<5% | 10227 | 0.027 | 0.000 | 0.045 | 129 | 113.856 | 0.154 | 2 |
| within_colony_net_time | Ptua | 5-20% | 4324 | -0.013 | 0.000 | 0.602 | 86 | 56.302 | 0.046 | 2 |
| within_colony_net_time | Ptua | \>20% | 2622 | -0.012 | 0.000 | 0.717 | 51 | 36.017 | 0.114 | 3 |

``` r
hl <- melt(per_gene, id.vars = c("species", "gene", "meth_class"),
           measure.vars = c("r_within_colony", "r_within_colony_net_time"),
           variable.name = "test", value.name = "r")
hl[, test := sub("^r_", "", test)]
p2 <- ggplot(hl[meth_class != "<5%"], aes(r, fill = species)) +
  geom_histogram(binwidth = 0.05, boundary = 0, alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = 2) +
  facet_grid(test ~ species, scales = "free_y") +
  labs(x = "Within-colony correlation of methylation and expression (genes >= 5% methylated)",
       y = "Genes") +
  theme_bw() + theme(legend.position = "none")
ggsave(file.path(out_dir, "within-colony-r.png"), p2, width = 10, height = 6, dpi = 200)
p2
```

![](56-within-colony-methylation-expression_files/figure-gfm/within-plot-1.png)<!-- -->

## 3.4 Strongest coupled genes

Genes at 10% FDR in the stricter test (colony-specific change), if any.

``` r
top <- per_gene[fdr_within_colony_net_time < 0.1][order(species, p_within_colony_net_time)]
fwrite(top, file.path(out_dir, "coupled-genes-fdr10.csv"))
if (nrow(top)) knitr::kable(top[, .(species, gene, mean_meth, r_within_colony_net_time,
                                    p_within_colony_net_time, fdr_within_colony_net_time)],
                            digits = 4) else cat("No genes at FDR < 0.1.\n")
```

| species | gene | mean_meth | r_within_colony_net_time | p_within_colony_net_time | fdr_within_colony_net_time |
|:---|:---|---:|---:|---:|---:|
| Apul | FUN_000537 | 23.1149 | 0.7981 | 0e+00 | 0.0209 |
| Apul | FUN_000352 | 29.4800 | -0.7961 | 0e+00 | 0.0209 |
| Apul | FUN_026369 | 14.8308 | -0.7863 | 0e+00 | 0.0209 |
| Apul | FUN_041965 | 56.7120 | 0.7821 | 0e+00 | 0.0209 |
| Apul | FUN_038439 | 9.8883 | -0.7708 | 0e+00 | 0.0276 |
| Apul | FUN_009661 | 74.8637 | -0.7647 | 0e+00 | 0.0298 |
| Apul | FUN_004400 | 22.9832 | 0.7566 | 0e+00 | 0.0335 |
| Apul | FUN_016780 | 68.2891 | 0.7547 | 0e+00 | 0.0335 |
| Apul | FUN_023795 | 48.0181 | -0.7327 | 0e+00 | 0.0583 |
| Apul | FUN_022869 | 8.2658 | 0.7320 | 0e+00 | 0.0583 |
| Apul | FUN_008936 | 64.7663 | 0.7316 | 0e+00 | 0.0583 |
| Apul | FUN_035276 | 65.5111 | -0.7269 | 1e-04 | 0.0632 |
| Apul | FUN_005725 | 67.1937 | 0.7213 | 1e-04 | 0.0687 |
| Apul | FUN_007759 | 63.1600 | -0.7200 | 1e-04 | 0.0687 |
| Apul | FUN_031709 | 65.2075 | -0.7114 | 1e-04 | 0.0818 |
| Apul | FUN_015256 | 51.4770 | -0.7104 | 1e-04 | 0.0818 |
| Apul | FUN_038220 | 25.1344 | -0.7074 | 1e-04 | 0.0818 |
| Apul | FUN_023189 | 42.4361 | 0.7066 | 1e-04 | 0.0818 |
| Apul | FUN_037461 | 34.7056 | -0.7052 | 1e-04 | 0.0818 |
| Apul | FUN_032288 | 71.8579 | 0.7027 | 1e-04 | 0.0818 |
| Apul | FUN_006778 | 24.2319 | 0.7019 | 1e-04 | 0.0818 |
| Apul | FUN_009304 | 54.5193 | -0.7009 | 1e-04 | 0.0818 |
| Apul | FUN_037208 | 32.2610 | -0.6930 | 2e-04 | 0.0958 |
| Apul | FUN_040402 | 14.9296 | 0.6927 | 2e-04 | 0.0958 |
| Apul | FUN_024778 | 57.8506 | -0.6916 | 2e-04 | 0.0958 |
| Apul | FUN_015293 | 74.1174 | 0.6890 | 2e-04 | 0.0999 |
| Peve | gene-Peve_00043373 | 61.8427 | -0.7992 | 0e+00 | 0.0604 |
| Peve | gene-Peve_00027826 | 0.3835 | -0.7699 | 0e+00 | 0.0755 |
| Peve | gene-Peve_00004045 | 35.5701 | 0.7676 | 0e+00 | 0.0755 |
| Peve | gene-Peve_00031737 | 19.8930 | -0.7566 | 0e+00 | 0.0855 |
| Ptua | gene-Pocillopora_meandrina_HIv1\_\_\_RNAseq.g1587.t1 | 30.3281 | -0.9049 | 0e+00 | 0.0103 |
| Ptua | gene-Pocillopora_meandrina_HIv1\_\_\_RNAseq.g3784.t1 | 0.8249 | 0.8938 | 0e+00 | 0.0115 |
| Ptua | gene-Pocillopora_meandrina_HIv1\_\_\_RNAseq.g17956.t1 | 21.3060 | -0.8787 | 0e+00 | 0.0197 |
| Ptua | gene-Pocillopora_meandrina_HIv1\_\_\_RNAseq.g10783.t1a | 16.0926 | -0.8553 | 0e+00 | 0.0428 |
| Ptua | gene-Pocillopora_meandrina_HIv1\_\_\_RNAseq.3753_t | 0.6029 | 0.8547 | 0e+00 | 0.0428 |
| Ptua | gene-Pocillopora_meandrina_HIv1\_\_\_RNAseq.g3991.t1a | 7.6566 | 0.8376 | 0e+00 | 0.0703 |
| Ptua | gene-Pocillopora_meandrina_HIv1\_\_\_RNAseq.g1582.t1 | 25.6533 | -0.8364 | 0e+00 | 0.0703 |

``` r
sessionInfo()
```

    ## R version 4.5.1 (2025-06-13)
    ## Platform: aarch64-apple-darwin20
    ## Running under: macOS Sequoia 15.6.1
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
    ## 
    ## locale:
    ## [1] en_US/en_US/en_US/C/en_US/en_US
    ## 
    ## time zone: America/Los_Angeles
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] ggplot2_4.0.1               DESeq2_1.50.2              
    ##  [3] SummarizedExperiment_1.40.0 Biobase_2.70.0             
    ##  [5] MatrixGenerics_1.22.0       matrixStats_1.5.0          
    ##  [7] GenomicRanges_1.62.1        Seqinfo_1.0.0              
    ##  [9] IRanges_2.44.0              S4Vectors_0.48.0           
    ## [11] BiocGenerics_0.56.0         generics_0.1.4             
    ## [13] data.table_1.17.8          
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] rappdirs_0.3.3      SparseArray_1.10.10 lattice_0.22-7     
    ##  [4] digest_0.6.39       magrittr_2.0.4      evaluate_1.0.5     
    ##  [7] grid_4.5.1          RColorBrewer_1.1-3  fastmap_1.2.0      
    ## [10] R.oo_1.27.1         jsonlite_2.0.0      Matrix_1.7-4       
    ## [13] R.utils_2.13.0      scales_1.4.0        coro_1.1.0         
    ## [16] textshaping_1.0.4   codetools_0.2-20    httr2_1.2.1        
    ## [19] abind_1.4-8         cli_3.6.5           rlang_1.3.0        
    ## [22] XVector_0.50.0      R.methodsS3_1.8.2   withr_3.0.2        
    ## [25] DelayedArray_0.36.1 yaml_2.3.10         S4Arrays_1.10.1    
    ## [28] tools_4.5.1         parallel_4.5.1      BiocParallel_1.44.0
    ## [31] dplyr_1.2.1         locfit_1.5-9.12     curl_7.0.0         
    ## [34] vctrs_0.7.3         R6_2.6.1            lifecycle_1.0.5    
    ## [37] ragg_1.5.0          pkgconfig_2.0.3     pillar_1.11.1      
    ## [40] gtable_0.3.6        glue_1.8.0          Rcpp_1.1.2         
    ## [43] systemfonts_1.3.1   tidyselect_1.2.1    xfun_0.54          
    ## [46] tibble_3.3.0        knitr_1.50          farver_2.1.2       
    ## [49] htmltools_0.5.8.1   labeling_0.4.3      rmarkdown_2.30     
    ## [52] ellmer_0.4.0        compiler_4.5.1      S7_0.2.1
