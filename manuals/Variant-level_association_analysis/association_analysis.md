# Variant-level Association Analysis

## Objectives
In this practical, you will learn how to perform association analysis in genome-wide scale while adjusting for confounding factors.

## Prequisites
We will start with the dataset of 480,780 variants for 992 samples passing QC from the previous **Sample array QC** practical. To ensure that we are on the same dataset,

<pre><code>md5sum ~/Practicals/<b>[yourname]</b>/practical2_QC/chrAll.ASA.afterSampleQC.afterVariantQC.*
</code></pre>

| Md5sum                           | File                                        |
| -------------------------------- |---------------------------------------------|
| 7c6619616f0391886b41bd54eb340ee9 | chrAll.ASA.afterSampleQC.afterVariantQC.bed |
| dd8f8e4c40207a8b5aff1ec60c25fe83 | chrAll.ASA.afterSampleQC.afterVariantQC.bim |
| 798746e8bcc2c18110dac62ee78ff602 | chrAll.ASA.afterSampleQC.afterVariantQC.fam |
 
If the md5sum values are consistent, you can copy the PLINK binary file to a new folder named as `practical3_associationtest`
<pre><code>mkdir ~/Practicals/<b>[yourname]</b>/practical3_associationtest
cd ~/Practicals/<b>[yourname]</b>/practical3_associationtest
cp ~/Practicals/<b>[yourname]</b>/practical2_QC/chrAll.ASA.afterSampleQC.afterVariantQC.* .
</code></pre>
Otherwise, please copy the dataset passing QC from the `~/Day2_association_analysis/` directory
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
Besides the 6th column in PLINK .ped and .fam file, you can supply another phenotype file via the command of `--pheno [filename]`.

First, let's use R to examine the content of the `CAD_LDL.pheno` phenotype file<br>
:closed_book: **Q1.** How many cornoary artery disease (CAD) cases are there? 
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

:closed_book: **Q3.** Is there any relationship between age, LDL and CAD?
<details>
<summary> Answer </summary>

 + CAD patients are generally older than controls (_P_=2.6x10<sup>-43</sup>) and LDL level increases by age (_P_=0.005) in general.
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
The `--adjust` option generates an .adjusted file containing several basic multiple testing corrections (e.g. FDR) for the raw p-values. By default, it will also estimate the genomic control factor (lambda) from the data.
 
By default, sex is automatically added as a covariate for all X chromosome SNPs but not anywhere else.

:closed_book: **Q1.** Which SNP is the top SNP? Can you give the effect size in beta and the confidence intervals (CI)?
<details>
  <summary> Answers </summary>
  
 + **Answer 1:** The top SNP is **rs2075650** on chr19 with association _p_-value of 3.8x10<sup>-32</sup>
```bash
head -n 2 LDL.assoc.linear.adjusted 
```
 
+ **Answer 2**: SNP effect (95% CI) = -0.9111 (-1.057, -0.7652).<br>
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
 
- Plot the manhattan plot to visualize the association results  
```bash 
# Rscript practical3.manhattanPlot.R <assoc_file> <outpu_prefix>
Rscript practical3.manhattanPlot.R LDL.assoc.linear QQPlot.LDL
```
 
- Plot the quantile-quantile plot to examine for inflation 
```bash 
Rscript practical3.QQPlot.R LDL.assoc.linear QQPlot.LDL
```
 
:closed_book: **Q2.** What is the lamda? Do the lamda and QQ plot indicate inflation of the association statistics?
<details>
  <summary> Answer </summary>

+ lamda=1.00219. No inflation in summary statistics.
</details>

### 2.  Perform linear regression while adjusted for age
You can use `--covar <filename>` to specify the file with covariate(s), i.e. `CAD_LDL.pheno` in this example. The covariate(s) used for adjustment is specified through `--covar-name`.
<pre><code>
plink --bfile chrAll.ASA.afterSampleQC.afterVariantQC \
  --chr 1-22,X,XY \
  --pheno CAD_LDL.pheno \
  --pheno-name LDL \
  <b>--covar CAD_LDL.pheno \</b>
  <b>--covar-name AGE \</b>
  <b>--hide-covar \</b>
  --linear --adjust \
  --out LDL.adj-AGE
```

- Plot the association results without adjustment of age and compare with the one without adjustment
```bash 
Rscript practical3.manhattanPlot.R LDL.adj-AGE.assoc.linear manhattanPlot.LDL.adj-AGE
```

3.  Perform logistic regression adjusted for age to test for association with CAD
```bash 
plink --bfile chrAll.ASA.afterSampleQC.afterVariantQC \
 --chr 1-22,X,XY \
 --pheno CAD_LDL.pheno --pheno-name CAD \
 --covar CAD_LDL.pheno --covar-name AGE --hide-covar \
 --logistic \
 --out CAD.adj-AGE
```

- Plot the manhattan plot and QQ plot to evaluate the association
```bash 
Rscript practical3.QQPlot.R CAD.adj-AGE.assoc.logistic QQPlot.CAD.adj-AGE
Rscript practical3.manhattanPlot.R CAD.adj-AGE.assoc.logistic manhattanPlot.CAD.adj-AGE
```
### Step 3: Visualizing the association results

### Step 4: Annotating the association findings
