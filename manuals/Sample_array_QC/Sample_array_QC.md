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

<p align="center">
<img src="https://user-images.githubusercontent.com/8644480/170681398-e29f945e-fc94-4876-b695-9a8f2250968e.png"  width="500" height="300">
</p>

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
> chrAll.ASA.beforeQC.imiss    # per individual<br>
> chrAll.ASA.beforeQC.lmiss    # per variant<br>

For both files, the last three columns measure the missingness for each individual or variant 
> N_MISS: Number of missing genotype call(s)<br>
> N_GENO: Number of genotype call(s)<br>
> F_MISS: Missing call rate<br>

- Plot the distribution of missingness 
We can then use R script to generate histogram of missingness
<pre>
# ========================== R code ==========================
# ----------------------------------------------------------#
#         (1)  SAMPLE CALL RATE    - threshold = 98%        #
# ----------------------------------------------------------#
imiss<-read.table("chrAll.ASA.beforeQC.imiss",h=T)
head(imiss)
summary(imiss$F_MISS)

# Plot missingness across samples
hist(imiss$F_MISS, freq=T, col="darkred", border ="black", main="Sample Call Rate", xlab="F_MISS", ylab="Number of samples")

# Plot missingness with altered y-axis for a zoom in view
hist(imiss$F_MISS, breaks=seq(0,0.2,0.01), freq=T, col="darkred", border ="black", main="Sample Call Rate", xlab="F_MISS", ylab="Number of samples",ylim=c(0,20))
# ============================================================
</pre>
![practical2 missing-hist1](https://user-images.githubusercontent.com/8644480/170730926-95e94bab-26a7-487b-beed-92cb352237bc.png)
![practical2 missing-hist2](https://user-images.githubusercontent.com/8644480/170731155-cad32ec4-a5a9-48a3-bec6-3492cf4e3471.png)

:closed_book: **Q:** Can you try to plot the right number of excluded samples?
<details>
  <summary>You can try some basic R codes by yourself first</summary>

- Answer 1
<pre><code>
#==== R =====
hist(imiss$F_MISS, breaks=50, freq=T, col="darkred", border="black", main="Sample Call Rate", xlab="F_MISS", ylab="Number of samples", ylim=c(0,100), xlim=c(0,0.2))
abline(v=0.02, lwd=2, lty=2, col="darkblue")
#abline(v=0.01, lwd=2, lty=2, col="darkgreen")
</code></pre>
![practical2 missing-hist3](https://user-images.githubusercontent.com/8644480/170732092-20f91ff2-1aa2-4d70-9d5b-7943f4e1954e.png)

- Answer 2
<pre><code>
#==== R =====
plot(sort(imiss$F_MISS), pch=20, col="darkred", main="Sample Call Rate", xlab="ASA samples", ylab="F_MISS")
abline(h=0.02, lwd=2, lty=2, col="darkblue")
#abline(h=0.01, lwd=2, lty=2, col="darkgreen")
</code></pre>
![practical2 missing-hist4](https://user-images.githubusercontent.com/8644480/170732149-7791c2b8-48e9-4f73-a8ea-83b2b6025dda.png)

- Answer 3
<pre><code>
#==== R =====
plot(sort(imiss$F_MISS), pch=20, col="darkred", main="Sample Call Rate", xlab="ASA samples", ylab="F_MISS")
abline(v=0.02, lwd=2, lty=2, col="darkblue")
#abline(v=0.01, lwd=2, lty=2, col="darkgreen")
</code></pre>
![practical2 missing-hist5](https://user-images.githubusercontent.com/8644480/170732217-70b6a656-e270-4b1f-a09d-be908371f860.png)

</details>

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
