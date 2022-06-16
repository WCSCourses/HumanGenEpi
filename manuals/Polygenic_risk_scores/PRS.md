#  Polygenic risk score manual
In this manual you will learn how to compute polygenic risck scores using **PRSice** (https://www.prsice.info/) 

**You will need to download the following files for this practical** 

```bash
# In your home directory create a folder called PRS and download the files into it

mkdir PRS
cd PRS

wget https://github.com/WCSCourses/HumanGenEpi/raw/main/manuals/Polygenic_risk_scores/ASN.bed.xz

wget https://github.com/WCSCourses/HumanGenEpi/raw/main/manuals/Polygenic_risk_scores/ASN.pheno

wget https://github.com/WCSCourses/HumanGenEpi/raw/main/manuals/Polygenic_risk_scores/ASN.bim

wget https://github.com/WCSCourses/HumanGenEpi/raw/main/manuals/Polygenic_risk_scores/test

wget https://github.com/WCSCourses/HumanGenEpi/raw/main/manuals/Polygenic_risk_scores/validate

wget https://github.com/WCSCourses/HumanGenEpi/raw/main/manuals/Polygenic_risk_scores/ASN.gwas.txt.xz

wget https://github.com/WCSCourses/HumanGenEpi/raw/main/manuals/Polygenic_risk_scores/ASN.cov

wget https://github.com/WCSCourses/HumanGenEpi/raw/main/manuals/Polygenic_risk_scores/ASN.fam

wget https://github.com/WCSCourses/HumanGenEpi/raw/main/manuals/Polygenic_risk_scores/PRSice_linux

wget https://github.com/WCSCourses/HumanGenEpi/raw/main/manuals/Polygenic_risk_scores/PRSice.R

#Unzip the following files as shown below

xz -dv ASN.gwas.txt.xz
xz -dv ASN.bed.xz
```


## 1. Training the PRS to find the best predictive one using the validation target data set


***Its key to do some QC on your discovery data set***

``` bash
Rscript PRSice.R --prsice PRSice_linux --base ASN.gwas.txt --target ASN --binary F --keep validate --pheno ASN.pheno --cov ASN.cov --out trial1
```
Check out the trial1.log file which contains the errors shown above. So, we see that initially
the all the SNPs in the base file are read. Ambiguous variants are removed to avoid strand
errors, however, other strand flips are automatically detected and accounted for. You can
filter the variants to include based on the quality of imputation and if duplicate variants are
available the latest version of PRSice will stop as shown above. You can use the –extract
trial1.valid however this option will not show you the detailed breakdown of SNPs included
so we will remove these duplicates using the code below.

``` bash
uniq ASN.gwas.txt > ASN.gwasqc.txt
```

***Preliminary analysis using default settings***

``` bash
Rscript PRSice.R --prsice PRSice_linux --base ASN.gwasqc.txt --target ASN --keep validate --pheno ASN.pheno --binary F --cov ASN.cov --out Prelim
```

Let’s look at the .summary file and the plots and ensure you understand them. What is the
PRS R2 and how many SNPs are in the best preforming PRS ?

***Optimise computation to get the most predictive PRS***

*   clump-kb 500 clump-r2 0.1

``` bash 
Rscript PRSice.R --prsice PRSice_linux --base ASN.gwasqc.txt --target ASN --binary F --keep validate --pheno ASN.pheno --cov ASN.cov –-clump-kb 500 –-clump-r2 0.1 –-base-info INFO:0.4 --out Opt500_0.1
```
*  clump-kb 250 clump-r2 0.3

``` bash
Rscript PRSice.R --prsice PRSice_linux --base 
ASN.gwasqc.txt --binary F --target ASN --keep validate --pheno ASN.pheno --cov ASN.cov –-clump-kb 
250 –-clump-r2 0.3 –-base-info INFO:0.4 --out Opt250_0.5
```
***Which parametes give the best predictive PRS ? Lets print out the SNPs of the best predictive PRS***

``` bash
Rscript PRSice.R --prsice PRSice_linux --base 
ASN.gwasqc.txt --target ASN --keep validate --pheno ASN.pheno --binary F --cov ASN.cov –-clump-kb 500 –-clump-r2 0.1 –-base-info INFO:0.4 –-print-snp –-out validation

```

## 2. Apply the best PRS in the testing data set


Look for the validation.snp file and filter SNPs that are working best at the optimal p value threshold of **0.46875** indicated in the **validation.summary** file

``` bash 
awk '$4 <= 0.46875' validation.snp | awk '{print $2}' > PRS_snps
```
**Lets make a discovery data set of the SNPs for the best predictive PRS***

```bash
grep -wf PRS_snps ASN.gwasqc.txt | sed '1i\CHR BP SNP A1 A2	N SE P OR INFO MAF' > Bestprs_disc

```
***Finally lets run this best PRS in the test dataset and create a decile plot***

``` bash
Rscript PRSice.R --prsice PRSice_linux --base 
Bestprs_disc --target ASN --keep test --pheno ASN.pheno --binary F --cov ASN.cov –-no-clump –-keep-ambig --fastscore --bar-levels 1 –-base-info INFO:0.4 --quantile 10 --quant-break 1,2,3,4,5,6,7,8,9,10 --quant-ref 1 –-out test

```
