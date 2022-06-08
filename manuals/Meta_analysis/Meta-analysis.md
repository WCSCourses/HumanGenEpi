
In this module you will learn how to perform GWAS meta-analysis using the software **Metal** (https://genome.sph.umich.edu/wiki/METAL_Documentation)

We will start by familiarising with the config file that indicates the studies for the metaanalysis.
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

We have two GWAS summary statistics files indicated below and we should edit the *metal.txt* config file above to match with the headers of these files

```bash
EAS.txt
SAS.txt
```
Before we run our meta-analysis we need to check the following:

*1.* If the SNPs in chr:bp format on the same genome build
*2.* If the GWAS quantitative trait transformation similar is among the studies; if not we have to do a p-value based meta-analysis
*2.* If the quantitative trait transformation is similar we will run the meta-analysis based on the beta and standard error
*3.* There is no neef to worry about trait transformation for the case-control meta-analysis

By specifying the scheme in the *metal.txt* config file , we are able to run these types of analysis indicated below:

1. P value based meta-analyis
```bash
SCHEME SAMPLESIZE
```
You will need to specify N in the individual studies

2. Beta and standard error based meta-analysis

```bash
SCHEME   STDERR
```

You will need to specify BETA and SE in the individual studies

What differences do you see in the output of these analysis ?

Lastly we would like to assess the heterogeneity accross the study results. This we will do by modifying the metal.txt file .

```bash
ANALYZE HETEROGENEITY
```

Can you state the SNP which has the most significant heterogeneity across these studies?



