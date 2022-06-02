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
<pre><code>
mkdir ~/Practicals/<b>[yourname]</b>/practical3_associationtest
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
### Step 1: Examining the phenotype file
```bash
cp ~/Day2_association_analysis/CAD_LDL.pheno .
```
Besides the 6th column in PLINK .ped and .fam file, you can supply another phenotype file via the command of `--pheno [filename]`.

Use R to examine the phenotypes
1. How many cornoary artery disease (CAD) cases are there?
2. Is the low density lipoprotein (LDL) level normally distributed?
3. Is there any relationship between age, LDL and CAD?

<details>
<summary></summary>
```R
# ========================== R code ==========================
pheno <- read.table("CAD_LDL.pheno",h=T)
summary(pheno)
table(pheno$CAD)
```
```R
# Plot the histogram to assess if the LDL level is normally distributed
hist(pheno$LDL, freq=T, col="darkred", border ="black", xlab="LDL level", ylab="Number of samples")

# Plot the QQ plot to assess if the LDL level is normally distributed
qqnorm(pheno$LDL)
qqline(pheno$LDL)
```
```R
# Test if age is associated with LDL and CAD
summary(glm(LDL ~ AGE, data=pheno))
summary(glm(CAD==2 ~ AGE, family="binomial", data=pheno))
by(pheno$AGE, pheno$CAD, summary)
```
# ============================================================
```
</details>

### Step 2: Performing association analysis
Next, we will use PLINK to perform association test
1. Perform linear regression WITHOUT adjusting for age to test for association with LDL 
```bash  
plink --bfile chrAll.ASA.afterSampleQC.afterVariantQC --set-hh-missing --pheno CAD_LDL.pheno --pheno-name LDL --linear --adjust --out LDL
```
The `--adjust` option generates an .adjusted file containing several basic multiple testing corrections (e.g. FDR) for the raw p-values. By default, it will also estimate the genomic control factor (lambda) from the data.

:closed_book: **Q:** Which SNP is the top SNP? Can you give the effect size in beta (95% CI)?
<details>
  <summary> Answer </summary>
  
** Answer 1 **
** Answer 2 **
Add `--ci 0.95` while running the regression test in PLINK
</details>  

- Plot the quantile-quantile plot to examine for inflation 
```bash 
Rscript practical3.QQPlot.R LDL.assoc.linear QQPlot.LDL
```
:closed_book: **Q:** What is the lamda? Do the lamda and QQ plot indicate inflation of the association statistics?
<details>
  <summary> Answer </summary>

- lamda=1.00219. No inflation of summary statistics.
<details>

2.  Perform linear regression adjusted for age 
```bash  
plink --bfile chrAll.ASA.afterSampleQC.afterVariantQC --set-hh-missing --pheno CAD_LDL.pheno --pheno-name LDL --covar CAD_LDL.pheno --covar-name AGE --hide-covar --linear --adjust --out LDL.adj-AGE
```
- Plot and compare the association results with and without adjustment of age
```bash 
Rscript practical3.manhattanPlot.R LDL.assoc.linear manhattanPlot.LDL
Rscript practical3.manhattanPlot.R LDL.adj-AGE.assoc.linear manhattanPlot.LDL.adj-AGE
```

3.  Perform logistic regression adjusted for age to test for association with CAD
```bash 
plink --bfile chrAll.ASA.afterSampleQC.afterVariantQC --set-hh-missing --pheno practical3.pheno --pheno-name CAD --covar CAD_LDL.pheno --covar-name AGE --logistic --out CAD.adj-AGE
```

### Step 3: Visualizing the association results

### Step 4: Annotating the association findings
