# Sample-based and Variant-based Quality Control (QC)

## Objectives
In this practical, you will learn the basic steps of data quality assessment and perform QC on genotype data that are typically carried out in genome wide association studies (GWAS).
* [Sample QC](#Sample-QC)
* [Variant QC](#Variant-QC)

## Prerequisites
For sample and variant QC, we need the following software that is pre-installed in the VM
- `PLINK`[1.9 beta](https://www.cog-genomics.org/plink/) for data management 
- `R` for visualization of different QC measures

## Why do we need the QC on GWAS data?
Suboptimal study design, poor sample quality, and errors in genotype calling can introduce systematic biases into GWAS, leading to spurious associations. A thorough QC can help us identify samples and variants that should be removed prior to association analysis to minimize false-positive and false-negative associations.

Here we assume that the study design has been conducted appropriately and basic QC to genotypes (e.g. QC based on GenTrain score in Illumina GenomeStudio or rare SNPs calling using zCall) has been performed after clustering and calling from probe intensity data. 

## The QC protocol: from Sample QC to Variant QC 
The QC protocol for GWAS typically split into two broad categories, `Sample QC` and `Variant QC`. Sample QC is usually done prior to Variant QC in order to maximise the number of variants remained for association analysis.

<img src="https://user-images.githubusercontent.com/8644480/170681398-e29f945e-fc94-4876-b695-9a8f2250968e.png"  width="1024" height="468">

In this practical, we will first generate the summary profiles of 
+ missingness (Sample and variant)
+ sex chromosome heterogeneity and missingness
+ 

## Step 0: Download the genotype data
- First, create a directory named `QC` for this practical section

```bash
mkdir ~/practical2_QC
cd ~/practical2_QC
```

- Download the raw genotype data in binary format and save it to the current directory
```bash
wget https://github.com/WCSCourses/HumanGenEpi/raw/main/course_data/SNP_array_QC/practical2.tar.gz
```

### Dataset
The dataset used in this practical was simuated from haplotypes of East Asian samples of the 1000 Genomes Project ([Phase 3](https://www.internationalgenome.org/category/phase-3/)). SNPs included in the dataset reflect to those included in the [Illumina Asian Screening array](https://www.illumina.com/products/by-type/microarray-kits/infinium-asian-screening.html) designed to maximize the genomic coverage for East Asian population.

## Sample-QC
### Step_1: Individuals with excessive missing genotypes
- Obtain profile of missingness per individual and per SNP
```bash
plink --bfile chrAll.ASA --missing --out chrAll.ASA.beforeQC
```
This command generates two files `.imiss` for sample-based and `.lmiss` for variant-based
> chrAll.ASA.beforeQC.imiss    # per individual
> chrAll.ASA.beforeQC.lmiss    # per variant

For both files, the last three columns measure the missingness for each individual or variant 
> N_MISS: Number of missing genotype call(s)
> N_GENO: Number of genotype call(s)
> F_MISS: Missing call rate

- Plot the distribution of missingness 
We can then use R script to generate histogram of missingness
<pre>
# ========================== R code ==========================
# ----------------------------------------------------------#
#         (1)  SAMPLE CALL RATE    - threshold = 99%        #
# ----------------------------------------------------------#
imiss<-read.table("chrAll.ASA.beforeQC.imiss",h=T)
head(imiss)
summary(imiss$F_MISS)
</pre>


<pre>
# ========================== R code ==========================
hist(imiss$F_MISS, freq=TRUE, col="blue", border ="black", main = "Sample Call Rate", sub = Cohort, xlab="F_MISS", ylab="Frequency")

pdf(imiss,"missingness.bySample.pdf")
plot(imiss$)
dev.off()
</pre>

### Step_2: Individuals with sex discrepancy
- Sex check
Extract xchr SNPs
plink --bfile chrAll.ASA --chr 23 --make-bed --out $DIR/$FILE-xchr

Run missingness on xchr SNPs
plink --bfile $DIR/$FILE-xchr --missing --out $DIR/$FILE-xchr-missin
g

### Step_3: Individuals with outlying heterozygosity rate
To avoid bias by genotyping error of rare variants and SNPs in strong LD, we usually perform the heterogeneity check using only common variants (MAF>=5%), excluding complex regions and SNPs in strong LD
```bash
plink --bfile chrAll.ASA --autosome --maf 0.05 --make-bed --out chrAll.ASA.autosome.maf05
plink --bfile chrAll.ASA --autosome --extract chrAll.ASA.autosome.maf05.bim --missing --out chrAll.ASA.autosome.maf05
plink --bfile chrAll.ASA --autosome --exclude chrAll.ASA.autosome.maf05.bim --missing --out chrAll.ASA.autosome.maf-lt-05
```
- LD pruning using R-squared 0.1
```bash
plink --bfile chrAll.ASA.autosome.maf05 --indep 200 50 0.1 --out chrAll.ASA.autosome.maf05.pruning
plink --bfile chrAll.ASA.autosome.maf05 --extract chrAll.ASA.autosome.maf05.pruning.prune.in --het --out chrAll.ASA.autosome.maf05.pruned
```

plink --bfile chrAll.ASA.autosome.maf05.pruned --het --out chrAll.ASA.autosome.maf05.pruned

```
### Step_4: Duplicated or related individuals
- Obtain pair-wise IBD for relatedness checking
```bash
plink --bfile chrAll.ASA.autosome.maf05.pruned --genome
```

### Step_5: Ethnicity outliers
