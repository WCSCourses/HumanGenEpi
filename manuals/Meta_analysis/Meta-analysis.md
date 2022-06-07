
In this module you will learn how to perform GWAS meta-analysis using the software **Metal** (https://genome.sph.umich.edu/wiki/METAL_Documentation)

We will start by familiarising with the config file that we will indicate the studies which whose summary statistics we want to aggregate.
``` bash

cat metal.txt

# Meta-analysis weighted by standard error does not work well
# because different studies used very different transformations
# SCHEME   STDERR

# Not sure if genomic control is a good idea, given the large
# number of true associations in these three regions ...
# GENOMICCONTROL ON

# To help identify allele flips, it can be useful to track
# allele frequencies in the meta-analysis.
# AVERAGEFREQ ON
# MINMAXFREQ ON

MARKER   SNP
WEIGHT   N
ALLELE   EFFECT_ALLELE NON_EFFECT_ALLELE
FREQ     EFFECT_ALLELE_FREQ
EFFECT   BETA
STDERR   SE
PVAL     P_VAL

PROCESS DGI_three_regions.txt

MARKER   SNP
ALLELE   EFFECT_ALLELE NON_EFFECT_ALLELE
FREQ     FREQ_EFFECT
WEIGHT   N
EFFECT   BETA
STDERR   SE
PVAL     PVALUE

PROCESS MAGIC_FUSION_Results.txt.gz

MARKER   SNP
DEFAULT  4108
ALLELE   AL1 AL2
FREQ     FREQ1
EFFECT   EFFECT
STDERR   SE
PVAL     PVALUE

PROCESS magic_SARDINIA.tbl

ANALYZE
```

We have two GWAS summary statistics files indicate below and we should edit the *metal.txt* config file above to match with these files

```bash
EAS.txt
SAS.txt
```
Before we run our meta-analysis we need to check the following

*1.* Are the SNPs in chr:bp format on the same genome build
*2.* Is the GWAS trait transformation similar between the studies if not we have to do a p-value based meta-analysis
*2.* If the quantitative trait transformation is similar we will run the meta-analysis based on the beta and standard error
*3.* There is no neef to worry about trait transformation for the case-control meta-analysis

By specifying the scheme in the *metal.txt* config file , we are now going to run

1. P value based meta-analyis

2. Beta and standard error based meta-analysis

What differences do you see in the output of these analysis ?

Lastly we would like to assess the heterogeneity of accross the study results. This we will do by modifying the metal.txt file .

Can you state the SNP which has the most significant heterogeneity across these studies?



