I have run barnacle as a means for gene expresssion pattern tensor decomposition on a three species corals; *Acropora pulchra*, *Pocillopora tuahiniensis*, and *Porites evermanni* 

The general implementation can be found in

M-multi-species/scripts/13.00-multiomics-barnacle.Rmd

with output at

M-multi-species/output/26-rank35-optimization/lambda_gene_0.2/barnacle_factors

1.  **Reads all CSV files** from `M-multi-species/output/26-rank35-optimization/lambda_gene_0.2/top_genes_per_component/`
2.  **Joins each file** with `M-multi-species/output/12-ortho-annot/ortholog_groups_annotated.csv` based on OG ID (ortholog group identifier)
3.  **Creates annotated files** with suffix `_annotation.csv` containing:
    -   Original data (OG ID and component values)
    -   Ortholog information (apul, peve, ptua gene IDs)
    -   Protein annotations (protein name, organism, UniProt IDs)
    -   GO annotations (GO biological process, cellular component, molecular function)
    -   GOSlim annotations

Then based on the tensor decomposition files in

M-multi-species/output/26-rank35-optimization/lambda_gene_0.2/barnacle_factors

component_weights.csv metadata.json sample_mapping.csv gene_factors.csv sample_factors.csv time_factors.csv

and

Physiological data to correspoding to the samples avialbale and to be read in from

<http://gannet.fish.washington.edu/seashell/snaps/MOSAiC_Physiology_Data__MOSAiC.csv>

Develop a narrative story about what is physiologically occurring in the 3 species of corals across the 4 time points.

Generate compelling visuals to support the narrative.