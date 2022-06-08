#  Polygenic risk score manual
In this manual you will learn how to compute polygenic risck scores using **PRSice** (https://www.prsice.info/) and **PRScsx** (https://www.nature.com/articles/s41588-022-01054-7)

To begin with we will start with splitting our target data set into the validation and testing file using the Rscript called splitting.

**Rscript splitting.R**

This will generate two files that we will use in the down stream analysis as  test and validation data sets

 We will now start with the PRS software using the following data
```bash
1.ASN.fam
2.ASN.bim
3.ASN.bed
````
## PRSice

* Visualise the files
``` bash
head heightgwas.txt
head test
head valid

```
* QC parameters for the base file

``` bash
Rscript /usr/local/bin/PRSice.R --prsice usr/local/bin/PRSice_linux --base
Height.gwas.txt --target ASN --pheno ASN.height --cov ASN.cov --out trial1
```
Check out the trial1.log file which contains the errors shown above. So, we see that initially
the all the SNPs in the base file are read. Ambiguous variants are removed to avoid strand
errors, however, other strand flips are automatically detected and accounted for. You can
filter the variants to include based on the quality of imputation and if duplicate variants are
available the latest version of PRSice will stop as shown above. You can use the –extract
trial1.valid however this option will not show you the detailed breakdown of SNPs included
so we will remove these duplicates using the code below.

``` bash
uniq Height.gwas.txt > Height.gwasqc.txt
```

* Preliminary analysis using default settings

``` bash
Rscript /usr/local/bin/PRSice.R --prsice usr/local/bin/PRSice_linux --base
Height.gwasqc.txt --target ASN --pheno ASN.height --cov ASN.cov –-base-info
INFO:0.4 --out Prelim
```

Let’s look at the .summary file and the plots and ensure you understand them. What is the
PRS R2 and how many SNPs are in the best preforming PRS ?

* Optimise computation to get the most predictive PRS

    * clump-kb 500 clump-r2 0.1

``` bash 
Rscript /usr/local/bin/PRSice.R --prsice usr/local/bin/PRSice_linux --base
Height.gwasqc.txt --target ASN --pheno ASN.height --cov ASN.cov–-clump-kb
500 –clump-r2 0.1 –-base-info INFO:0.4 --out Opt500_0.1
```
*  clump-kb 250 clump-r2 0.3

``` bash
Rscript /usr/local/bin/PRSice.R --prsice usr/local/bin/PRSice_linux --base
Height.gwasqc.txt --target ASN --pheno ASN.height --cov ASN.cov –-clump-kb
250 –clump-r2 0.3 –-base-info INFO:0.4 --out Opt250_0.5
```

## Replicate the best PRS in the testing data set

``` bash
Rscript /usr/local/bin/PRSice.R --prsice usr/local/bin/PRSice_linux --base
Height.gwasqc.txt --target ASN --pheno ASN.height --cov ASN.cov –-print-snp
–-out Validation

```
Look for the Validation.snp file and filter SNPs that are working best at the optimal p value threshold.

``` bash 
awk '$4 <= 0.13995' Validation.snp | awk '{print $2}' >
Replication_snps
```

