# Variant-level Association Analysis

## Objectives
In this practical, you will learn how to perform association analysis in genome-wide scale while adjusting for confounding factors.

## Prequisites
We will start with the dataset of 480,780 variants for 992 samples passing QC from the previous **Sample array QC** practical. To ensure that we are on the same dataset,

<pre><code>
md5sum ~/Practicals/<b>[yourname]</b>/practical2_QC/chrAll.ASA.afterSampleQC.afterVariantQC.*
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
- How many cornoary artery disease (CAD) cases are there?
- Is the low density lipoprotein (LDL) level normally distributed?
- Is there any relationship between age, LDL and CAD?
