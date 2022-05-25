# Sample-based and Variant-based Quality Control (QC)

## Objectives
In this practical, you will learn the basic steps of data quality assessment and perform QC on genotype data that are typically carried out in genome wide association studies (GWAS).

## Why do we need the QC on GWAS data?
Suboptimal study design, poor sample quality, and errors in genotype calling can introduce systematic biases into GWAS, leading to spurious associations. A thorough QC can help us identify samples and variants that should be removed prior to association analysis to minimize false-positive and false-negative associations.

Here we assume that the study design has been conducted appropriately and basic QC to genotypes (e.g. QC based on GenTrain score in Illumina GenomeStudio or rare SNPs calling using zCall) has been performed after clustering and calling from probe intensity data. 

## The QC protocol: from Sample QC to Variant QC 
The QC protocol for GWAS typically split into two broad categories, `Sample QC` and `Variant QC`. Sample QC is usually done prior to Variant QC in order to maximise the number of variants remained for association analysis.

### Sample QC
- Missingness. High missingness indicates poor clustering -> Examples of genomeStudio clustering figure
- 

### Variant QC
- Missingness
- Violation of Hardy-Weinberg equilibrium
- Minor allele frequency (MAF): Rare variants are difficult to call and usually have higher genotyping errors

## Prerequisites
For sample and variant QC, we need the following software that is pre-installed in this VM
- `PLINK`[1.9 beta](https://www.cog-genomics.org/plink/) for genotype data management 
- `R` for visualization of different QC measures

## Step 0: Download the genotype data
- First, create a directory named `QC` for this practical section

```bash
mkdir ~/QC
cd ~/QC
```

- Download the raw genotype data in binary format and save it to the current directory
```bash
FILEURL="http://"

wget $FILEURL
```


## Dataset
- Genotyping data simulated from 1000 Genomes East Asian population
- SNP extracted from Illumina Asian Screening array

