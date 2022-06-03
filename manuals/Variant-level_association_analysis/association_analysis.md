# Variant-level Association Analysis

## Objectives
In this practical, you will learn how to perform association analysis in genome-wide scale while adjusting for confounding factors.

## Prequisites
We will start with the dataset of 480,780 variants for 992 samples passing QC from the previous **Sample array QC** practical. To ensure that we are using the same dataset, let's check the md5sum of these PLINK files

<pre><code>md5sum ~/Practicals/<b>[yourname]</b>/practical2_QC/chrAll.ASA.afterSampleQC.afterVariantQC.*
</code></pre>

| Md5sum                           | File                                        |
| -------------------------------- |---------------------------------------------|
| 7c6619616f0391886b41bd54eb340ee9 | chrAll.ASA.afterSampleQC.afterVariantQC.bed |
| dd8f8e4c40207a8b5aff1ec60c25fe83 | chrAll.ASA.afterSampleQC.afterVariantQC.bim |
| 798746e8bcc2c18110dac62ee78ff602 | chrAll.ASA.afterSampleQC.afterVariantQC.fam |
 
If the md5sum values are consistent, you can copy the PLINK binary files to a new folder and named it as `practical3_associationtest`
<pre><code>mkdir ~/Practicals/<b>[yourname]</b>/practical3_associationtest
cd ~/Practicals/<b>[yourname]</b>/practical3_associationtest
cp ~/Practicals/<b>[yourname]</b>/practical2_QC/chrAll.ASA.afterSampleQC.afterVariantQC.* .
</code></pre>
Otherwise, please copy the whole dataset passing QC from the `~/Day2_association_analysis/` directory
<details>
 <summary>Code to copy</summary>
 
```bash
cp ~/Day2_association_analysis/chrAll.ASA.afterSampleQC.afterVariantQC.* .
```
</details>

## Association analysis
## Step 1: Examining the phenotype file
```bash
cp ~/Day2_association_analysis/CAD_LDL.pheno .
```
Besides the 6th column in PLINK .ped and .fam file, you can use another phenotype file via the command of `--pheno [filename]`.

First, let's use R to examine the content of the `CAD_LDL.pheno` phenotype file<br>
:closed_book: **Q1.** How many cornoary artery disease (CAD) patients and controls are there in this dataset? 
<details>
<summary> Answer </summary>

+ 496 cases and 496 controls
```R
# ========================== R code ==========================
pheno <- read.table("CAD_LDL.pheno",h=T)
summary(pheno)
table(pheno$CAD)
# ============================================================
```
</details>
 
:closed_book: **Q2.** Is the low density lipoprotein (LDL) level normally distributed?
<details>
<summary> Answer </summary>

+ Yes. The LDL level is normally distributed and no further transformation is needed
```R
# ========================== R code ==========================
# Plot the histogram to assess if the LDL level is normally distributed
hist(pheno$LDL, freq=T, col="darkred", border ="black", main="Distribution of LDL", xlab="LDL level", ylab="Number of samples")
```
```R
# Plot the QQ plot to assess if the LDL level is normally distributed
qqnorm(pheno$LDL)
qqline(pheno$LDL)
# ============================================================
```
</details>

:closed_book: **Q3.** Is there any relationship between age, LDL level and CAD?
<details>
<summary> Answer </summary>

 + CAD patients are generally older than controls (_P_=2.6x10<sup>-43</sup>) and LDL level increases with age (_P_=0.005).
```R
# ========================== R code ==========================
# Test if age is associated with LDL and CAD
summary(glm(LDL ~ AGE, data=pheno))
summary(glm(CAD==2 ~ AGE, family="binomial", data=pheno))
by(pheno$AGE, pheno$CAD, summary)
# ============================================================
```
</details>

## Step 2: Performing association analysis
Next, we will use PLINK to perform association test
 
### 1. Perform linear regression WITHOUT adjusting for age to test for association with LDL 
```bash  
plink --bfile chrAll.ASA.afterSampleQC.afterVariantQC \
  --chr 1-22,X,XY \
  --pheno CAD_LDL.pheno \
  --pheno-name LDL \
  --linear --adjust \
  --out LDL
```
By default, sex is automatically added as a covariate for all X chromosome SNPs but not anywhere else. The `--adjust` option generates an .adjusted file with basic multiple testing corrections (e.g. FDR) for the raw p-values. It also estimates the genomic control factor (lambda) from the data.
 
:closed_book: **Q1.** Which SNP is the top SNP? Can you give the effect size and the confidence interval (CI)?
<details>
  <summary> Answers </summary>
  
+ **Answer 1:** The top SNP is **rs2075650** on chr19 with association _p_-value of 3.8x10<sup>-32</sup>
```bash
head -n 2 LDL.assoc.linear.adjusted 
```
 
+ **Answer 2**: SNP effect beta (95% CI) = -0.9111 (-1.057, -0.7652).<br>
You can add `--ci 0.95` while running the regression test in PLINK or you can compute the 95% CI in R
<pre><code>plink --bfile chrAll.ASA.afterSampleQC.afterVariantQC \
  <b>--snp rs2075650</b> \
  --pheno CAD_LDL.pheno --pheno-name LDL \
  --linear <b>--ci 0.95</b> \
  --out LDL.topsnp 
</code></pre>
 
```R
# ========================== R code ==========================
L95 <- -0.9111-1.96*0.07448
L95
U95 <- -0.9111+1.96*0.07448
U95
# ============================================================
```
</details>  

:closed_book: **Q2.** What is the lamda? Does it indicate inflation of the association statistics?
<details>
  <summary> Answer </summary>

+ lamda=1.00236. No inflation in summary statistics.
</details>
 
- Plot the manhattan plot to visualize the association results  
```bash 
# Rscript practical3.manhattanPlot.R <assoc_file> <outpu_prefix>
Rscript practical3.manhattanPlot.R LDL.assoc.linear manhattanPlot.LDL
```
 
### 2.  Perform linear regression while adjusted for age
To adjust for confounding factor(s), you can use `--covar <filename>` to specify the covariate file, i.e. `CAD_LDL.pheno` in this example. The covariate(s) used for adjustment is specified through `--covar-name`.
 
<pre><code>plink --bfile chrAll.ASA.afterSampleQC.afterVariantQC \
  --chr 1-22,X,XY \
  --pheno CAD_LDL.pheno \
  --pheno-name LDL \
  <b>--covar CAD_LDL.pheno \</b>
  <b>--covar-name AGE \</b>
  <b>--hide-covar \</b>
  --linear --adjust \
  --out LDL.adj-AGE
</code></pre>

- Plot the association results without adjustment of age and compare it with the previous one without adjustment
```bash 
Rscript practical3.manhattanPlot.R LDL.adj-AGE.assoc.linear manhattanPlot.LDL.adj-AGE
```

### 3. Perform logistic regression with and without adjusting for age to test for association with CAD
- Without adjusting for age
```bash 
plink --bfile chrAll.ASA.afterSampleQC.afterVariantQC \
 --chr 1-22,X,XY \
 --pheno CAD_LDL.pheno --pheno-name CAD \
 --logistic \
 --out CAD
```
- Adjusting for age
```bash 
plink --bfile chrAll.ASA.afterSampleQC.afterVariantQC \
 --chr 1-22,X,XY \
 --pheno CAD_LDL.pheno --pheno-name CAD \
 --covar CAD_LDL.pheno --covar-name AGE --hide-covar \
 --logistic \
 --out CAD.adj-AGE
```
- Plot the manhattan plots to compare the difference
```bash 
Rscript practical3.manhattanPlot.R CAD.assoc.logistic manhattanPlot.CAD
Rscript practical3.manhattanPlot.R CAD.adj-AGE.assoc.logistic manhattanPlot.CAD.adj-AGE
```

### Step 3: Annotating the association findings
#### LocusZoom (https://statgen.github.io/localzoom/)<br>
LocalZoom is a tool for generating regional association plot via web browser or offline. For the web-based js version without uploading the summary statistics, a tabix-indexed tab-delimited text file is needed. If you are comfortable uploading your data to a server, you can consider using the my.locuszoom.org upload service with more custom function.
```bash
# (i) convert "fixed width" PLINK association result format to "tab-delimited"
# (ii) add the A2 allele (specifying the right Ref and Alt alleles are needed for LD plot)
awk 'BEGIN { OFS="\t"; a2["SNP"]="A2" } NR==FNR { a1[$2]=$5; a2[$2]=$6 } NR!=FNR { $4=$4"\t"a2[$2]; if (FNR!=1) $2=$1":"$3; print $0 }' \
 chrAll.ASA.afterSampleQC.afterVariantQC.bim \
 LDL.adj-AGE.assoc.linear | \
 bgzip > LDL.adj-AGE.assoc.linear.forLocuszoom.gz
 
# (iii) index the association file
tabix -s 1 -b 3 -e 3 --skip-lines 1 -f LDL.adj-AGE.assoc.linear.forLocuszoom.gz
```
+ First, make sure that `GWAS` and the right genome assembly are chosen, i.e. `GRCh37` for this example dataset
<img src="https://user-images.githubusercontent.com/8644480/171670740-86b20b26-cc20-4512-81cc-d43940be970c.png" width=800>
 
+ Click `Add tab-indexed datafile` > `Browse` to add the LDL.adj-AGE.assoc.linear.forLocuszoom **.gz** and **.gz.tbi** files
<img src="https://user-images.githubusercontent.com/8644480/171672531-00848874-2530-446a-aad4-102475c2ea62.png" width=800>
+ Choose `GWAS scattered plot` as `Type`
<img src="https://user-images.githubusercontent.com/8644480/171672261-3cf5120d-a24b-455a-96e6-959e262b2ef6.png" width=250>
 
+ Under `Variant from columns` > Choose the major allele `A2` as `Ref` allele and the minor allele `A1` as the `Alt` allele
<img src="https://user-images.githubusercontent.com/8644480/171672887-bc1c4203-0880-44a5-9f66-64b376a155a9.png" width=500>

+ Next > Accept option > Type `rs2075650` in the box > Go to the region
+ Choose `LD population: EAS`
+ You can also change the range of region to `19:45300000-455000000` for a zoom in view
<img src="https://user-images.githubusercontent.com/8644480/171678660-f12f0ff9-c4d1-4472-bf9d-089ad58d0164.png" width=800>

#### FUMA (https://fuma.ctglab.nl/)
