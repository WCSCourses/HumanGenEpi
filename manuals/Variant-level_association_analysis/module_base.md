# Variant-level Association Analysis

## Objectives
In this practical, you will learn how to perform association analysis in genome-wide scale while adjusting for confounding factors.

## Prequisites
We will start with the dataset of 480,780 variants for 992 samples passing QC from the previous **Sample array QC** practical. To ensure that we are on the same dataset,

```bash
md5sum ~/practical2_QC/chrAll.ASA.afterSampleQC.afterVariantQC.*
```

| Md5sum                           | File                                        |
| -------------------------------- |---------------------------------------------|
| 7c6619616f0391886b41bd54eb340ee9 | chrAll.ASA.afterSampleQC.afterVariantQC.bed |
| dd8f8e4c40207a8b5aff1ec60c25fe83 | chrAll.ASA.afterSampleQC.afterVariantQC.bim |
| 798746e8bcc2c18110dac62ee78ff602 | chrAll.ASA.afterSampleQC.afterVariantQC.fam |
 
If the md5sum values are consistent, you can copy the PLINK binary file to a new folder named as `practical3_associationtest`
```bash
mkdir ~/practical3_associationtest
cd ~/practical3_associationtest
cp ~/practical2_QC/chrAll.ASA.afterSampleQC.afterVariantQC.* ~/practical3_associationtest/
```
Otherwise, please download the dataset passing QC using the following command
```bash
wget xxxx
```
  
  
