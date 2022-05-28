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
## Step_1: Individuals with excessive missing genotypes
- Obtain profile of missingness per individual and per SNP
```bash
plink --bfile chrAll.ASA --missing --out chrAll.ASA.beforeQC
```
This command generates two files, including `.imiss` for sample-based and `.lmiss` for variant-based missingness
> chrAll.ASA.beforeQC.imiss    # per individual<br>
> chrAll.ASA.beforeQC.lmiss    # per variant<br>

For both files, the last three columns measure the missingness for each individual or variant 
> N_MISS: Number of missing genotype call(s)<br>
> N_GENO: Number of genotype call(s)<br>
> F_MISS: Missing call rate<br>

- Plot the distribution of missingness 
We can then use R script to generate histogram of missingness
```R
# ========================== R code ==========================
# ----------------------------------------------------------#
#         (1)  SAMPLE CALL RATE    - threshold = 98%        #
# ----------------------------------------------------------#
imiss <- read.table("chrAll.ASA.beforeQC.imiss",h=T)
head(imiss)
summary(imiss$F_MISS)

# Plot missingness across samples
hist(imiss$F_MISS, freq=T, col="darkred", border ="black", main="Sample Call Rate", 
xlab="F_MISS", ylab="Number of samples")

# Plot missingness with altered y-axis for a zoom in view
hist(imiss$F_MISS, breaks=seq(0,0.6,0.05), freq=T, col="darkred", border ="black", 
main="Sample Call Rate", xlab="F_MISS", ylab="Number of samples",ylim=c(0,10))
# ============================================================
```
![practical2 missing-hist1](https://user-images.githubusercontent.com/8644480/170830897-4eb7a450-bc51-42ea-b68b-bf5ea98197ab.png)
![practical2 missing-hist2](https://user-images.githubusercontent.com/8644480/170830900-ef75c54a-ba80-4393-9c77-7502ac7c319a.png)

:closed_book: **Q:** Can you try to plot the right number of excluded samples?
<details>
  <summary> Try your own R codes </summary>
<p></p>

- **Answer 1**
```R
#==== R =====
hist(imiss$F_MISS, breaks=seq(0,0.6,0.01), freq=T, col="darkred", border="black", main="Sample Call Rate", xlab="F_MISS", ylab="Number of samples", ylim=c(0,100))
abline(v=0.02, lwd=2, lty=2, col="darkblue")
abline(v=0.01, lwd=2, lty=2, col="darkgreen")
```
![practical2 missing-hist3](https://user-images.githubusercontent.com/8644480/170831053-f941b2df-090a-4880-8711-0dd7a4b3fb77.png)

- **Answer 2**
```R
 #==== R =====
plot(sort(imiss$F_MISS), pch=20, col="darkred", main="Sample Call Rate", xlab="ASA samples", ylab="F_MISS")
abline(h=0.02, lwd=2, lty=2, col="darkblue")
abline(h=0.01, lwd=2, lty=2, col="darkgreen")
```
![practical2 missing-hist4](https://user-images.githubusercontent.com/8644480/170831061-7d3811a0-8502-4a6e-bc22-e4d5a45eb52e.png)

- **Answer 3**
```R
#==== R =====
plot(imiss$F_MISS,pch=20,col="darkred", main="Sample Call Rate", xlab="ASA samples", ylab="F_MISS")
abline(h=0.02, lwd=2, lty=2, col="darkblue")
abline(h=0.01, lwd=2, lty=2, col="darkgreen")
```
![practical2 missing-hist5](https://user-images.githubusercontent.com/8644480/170831069-58dd61a8-b4db-4ee2-a6a7-e63968b35108.png)

</details>

## Step_2: Individuals with sex discrepancy
- Obtain missingness of chr X
```bash
plink --bfile chrAll.ASA --chr 23 --missing --out chrX.ASA
```
- Check sex
```bash
plink --bfile chrAll.ASA --check-sex --out chrX.ASA
```
The function of `--check-sex` normally compares sex assignments in the input pedigree data with inbreeding coefficients (F) imputed from SNPs on X chromosome. By default, the F estimates smaller than 0.2 yield female calls, and values larger than 0.8 yield male calls. 

```R
# ========================== R code ==========================
sexcheck <- read.table("chrX.ASA.sexcheck",h=T)
head(sexcheck)

#       FID      IID PEDSEX SNPSEX  STATUS       F
# 1 id1_0211 id2_0211      2      2      OK -1.0180
# 2 id1_1202 id2_1202      2      0 PROBLEM  0.6545
# 3 id1_1130 id2_1130      2      2      OK -0.7449
# 4 id1_1064 id2_1064      2      2      OK -0.2048
# 5 id1_1109 id2_1109      2      0 PROBLEM  0.7954
# 6 id1_0108 id2_0108      2      0 PROBLEM  0.3521

mismatch <- sexcheck[sexcheck$STATUS=="PROBLEM",]

colsex <- c("darkblue","darkred")
plot(sexcheck$F, col=colsex[sexcheck$PEDSEX], main="Sex check", pch=20, xlab="ASA samples", ylab="chrX Inbreeding coefficient (F)")
points(row.names(mismatch),mismatch$F,pch=22,bg="green",col="black",lwd=2,cex=2)
abline(h=0.2, lwd=2, lty=2, col="red")
abline(h=0.8, lwd=2, lty=2, col="blue")
legend("bottomright",c("Male PEDSEX","Female PEDSEX","sample with PROBELM"), col=c(colsex,"black"),pt.bg="green", pch=c(20,20,22))
# =============================================================
```
![practical2 sexcheck](https://user-images.githubusercontent.com/8644480/170831468-747d54ae-1b5a-4057-aa00-319f9a2a4864.png)

- What is the relationship between missingness and inbreeeding coefficient for chrX?
```R
# ========================== R code ==========================
imiss.X <- read.table("chrX.ASA.imiss",h=T)
sexcheck.imiss <- merge(imiss.X, sexcheck, by="IID")
mismatch.imiss <- sexcheck.imiss[sexcheck.imiss$STATUS=="PROBLEM",]
plot(sexcheck.imiss$F_MISS, sexcheck.imiss$F, pch=20, col=colsex[sexcheck.imiss$PEDSEX], xlab="Number of missing genotypes on chrX", ylab="chrX Inbreeding coefficient (F)")
points(mismatch.imiss$F_MISS, mismatch.imiss$F, pch=22, bg="green", col="black", lwd=2, cex=1.5)
abline(h=0.2, lwd=2, lty=2, col="red")
abline(h=0.8, lwd=2, lty=2, col="blue")
legend("bottomleft",c("Male PEDSEX","Female PEDSEX","sample with PROBELM"), col=c(colsex,"black"),pt.bg="green", pch=c(20,20,22))
# =============================================================
```
![practical2 sexcheck-fmiss](https://user-images.githubusercontent.com/8644480/170832150-fb1471ed-be27-4e1e-80eb-cdf78989d21f.png)

In addition to poor sample or genotyping quality, X chromosome aneuploidy, such as Turner syndrome (e.g. 45,X0) and Klinefelter syndrome (47, XXY), may lead to abnormal heterogenity. Missingness for SNPs on chr Y can be used to impute sex using `--check-sex ycount [female max F] [male min F] [female max Y obs] [male min Y obs]`. Before determining the minimum and maximum threshold of observed chrY variants, we can first set [female max Y obs] to maximum number of chrY variants and [male min Y obs] to 0  
```bash
plink --bfile chrAll.ASA --check-sex ycount 0.2 0.8 805 0 --out chrXY.ASA
```
:closed_book: **Q:** Do samples with sex discrepancy or abnormal heterogenity look like having X chromosome aneuploidy?
<details>
  <summary> Try your own PLINK / R codes </summary>
<p></p>
  
- ** 1 **
```bash
egrep PROBlEM chrXY.ASA.sexcheck
#         FID       IID PEDSEX SNPSEX  STATUS       F YCOUNT
#    id1_1202  id2_1202      2      0 PROBLEM 0.65450    607
#    id1_1109  id2_1109      2      0 PROBLEM 0.79540    636
#    id1_0108  id2_0108      2      0 PROBLEM 0.35210    569
#   CHSHet002 CHSHet002      1      0 PROBLEM 0.67790    804
#     id1_300   id2_300      1      2 PROBLEM 0.03826      2
#     id1_301   id2_301      2      1 PROBLEM 1.00000    799
#     id1_500   id2_500      1      2 PROBLEM 0.01425      2 
#     id1_501   id2_501      2      1 PROBLEM 1.00000    805
```
- ** 2 **
```R
# ========================== R code ==========================
sexcheck.XY <- read.table("chrXY.ASA.sexcheck",h=T)
mismatch.XY <- sexcheck.XY[sexcheck.XY$STATUS=="PROBLEM",]
colsex <- c("darkblue","darkred")
plot(sexcheck.XY$YCOUNT, sexcheck.XY$F, pch=20, col=colsex[sexcheck.XY$PEDSEX], xlab="Number of non-missing genotypes on chrY", ylab="chrX Inbreeding coefficient (F)")
points(mismatch.XY$YCOUNT, mismatch.XY$F, pch=22, bg=colsex[mismatch.XY$PEDSEX], col="white", lwd=1, cex=1.5)
abline(h=0.2, lwd=2, lty=2, col="red")
abline(h=0.8, lwd=2, lty=2, col="blue")
# =============================================================
```  
![practical2 sexcheck-XY](https://user-images.githubusercontent.com/8644480/170838631-9463a53c-2f21-4151-bf44-b91d5d33111e.png)
```bash
plink --bfile chrAll.ASA --check-sex ycount 0.2 0.8 100 700 --out chrXY.ASA.2
egrep -h PROBLEM chrXY.ASA.sexcheck chrXY.ASA.2.sexcheck | sort | uniq -u
```
  
</details>  

### Step_3: Individuals with outlying heterozygosity rate
To avoid bias by genotyping error of rare variants and linkage disequilibrium, we usually perform the heterogeneity check using only common variants (MAF>=5%), excluding complex regions and SNPs in strong LD
```bash
plink --bfile chrAll.ASA --autosome --maf 0.05 --make-bed --out chr1-22.ASA.maf05
plink --bfile chr1-22.ASA.maf05 --missing --out chr1-22.ASA.maf05
```
- LD pruning using R-squared 0.1
```bash
plink --bfile chr1-22.ASA.maf05 --indep-pairwise 200 50 0.1 --out chr1-22.ASA.maf05.pruning
plink --bfile chr1-22.ASA.maf05 --extract chr1-22.ASA.maf05.pruning.prune.in --make-bed --out chr1-22.ASA.maf05.pruned
plink --bfile chr1-22.ASA.maf05.pruned --het --out chr1-22.ASA.maf05.pruned
```
```R
#imiss.autosome <- read.table("chr1-22.ASA.maf05.imiss",h=T)
#het.autosome <- read.table("chr1-22.ASA.maf05.pruned.het",h=T)
imiss.autosome <- read.table("chrAll.ASA.autosome.maf05.imiss",h=T)
het.autosome <- read.table("chrAll.ASA.autosome.maf05.pruned.het",h=T)
imiss.het <- merge(imiss.autosome, het.autosome, by="IID")
imiss.het$PCT_HET <- (imiss.het$N.NM. - imiss.het$O.HOM.)/imiss.het$N.NM.
plot(imiss.het$F_MISS, imiss.het$PCT_HET, pch=20, col="darkred", xlab="Number of missing genotypes on chrX", ylab="chrX Inbreeding coefficient (F)")

```


:closed_book: **Q:** Can you plot the distribution of missingness of SNPs on chrY?
<details>
  <summary> Try your own R codes </summary>
<p></p>
- **Answer 3**
<pre><code>#==== R =====
plot(imiss$F_MISS,1:nrow(imiss),pch=20,col="darkred", main="Sample Call Rate", xlab="F_MISS", ylab="ASA samples")
abline(v=0.02, lwd=2, lty=2, col="darkblue")
#abline(v=0.01, lwd=2, lty=2, col="darkgreen")
</code></pre>
</details>

<img width="466" alt="BA ASA FMISS_hetF" src="https://user-images.githubusercontent.com/8644480/170811103-cb0d77fe-4290-4b9f-9cd6-a270730b9c8a.png">


### Step_4: Duplicated or related individuals
- Obtain pair-wise IBD for relatedness checking
```bash
plink --bfile chrAll.ASA.autosome.maf05.pruned --genome
```

### Step_5: Ethnicity outliers
