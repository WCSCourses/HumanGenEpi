# Introduction to PLINK

## Objectives
In this practical, you will learn how to use PLINK to manage genotyping data.

## Step 1: Understand the input data formats in PLINK
First, create a working directory named `practical1` and change to the directory
```bash
mkdir practical1
cd practical1
```
Download the genotype data in regular PLINK text file format (`PED / MAP`) into the directory
```bash
wget 'http://github..../practical1.zip'
```
```bash
unzip practical1.zip'
```
#### - PLINK text variant file (MAP)
Let's have a look at the variant `MAP` file
```bash
head practical1.map
```
MAP file stores the variant information and must contain as many markers as are in the PED file.<br>
It contains no header line and each line describes a single marker with 4 columns:

> 1. Chromosome code (1-22, 23/X, 24/Y, 25/XY, 26/MT; '0' indicates unknown)<br>
> 2. Variant identifier<br>
> 3. Position in morgans or centimorgans (safe to use dummy value of '0')<br>
> 4. Base-pair coordinate (1-based)

:green_book: **Q:** How many SNPs are there?
```bash
wc -l practical1.map
```
:green_book: **Q:** Where are these SNPs located in?
```bash
cut -f 1 practical1.map | sort | uniq -c
```

#### - PLINK text pedigree and genotype file (PED)
Let's have a look at the `PED` file
```bash
less -S practical1.ped   # type 'q' to quit
```
The PED file is a white-space (space or tab) delimited file<br>
It has no header line, and one line per sample with 2**V**+6 fields where **V** is the number of variants. The first six columns are mandatory and are the same as  in a `FAM` file. 

> 1. Family ID ('FID')<br>
> 2. Individual ID('IID')<br>
> 3. Within-family ID of father ('0' if father isn't in dataset)<br>
> 4. Within-family ID of mother ('0' if mother isn't in dataset)<br>
> 5. Sex (1=male; 2=female; 0=unknown)<br>
> 6. Phenotype (1=control, 2=case, -9 / 0=missing for case-control; numeric for quantitative trait)<br>

The 7th and 8th fields are the alleles for the first variant in the MAP file. The 9th and 10th are allele calls for the second variant and so on.

:green_book: **Q:** How many samples are there?
```bash
wc practical1.ped
```
:green_book: **Q:** Which samples are related?
```bash
awk '$3!=0 || $4!=0' practical1.ped | cut -d' ' -f 1-6     # space as delimiter
```

## Step 2: Data conversion in PLINK
Read the `practical1` PLINK text fileset and convert to the PLINK binary fileset (`BED / BIM / FAM`)
<pre><code>plink <b>--file practical1</b> --make-bed --out practical1_1

## Equivalent to 
# plink <b>--file practical1</b> --out practical1_1
# plink <b>--ped practical1.ped --map practical1.map</b> --make-bed --out practical1_1
</code></pre>
The command `--file practical1` reads in **practical1**.ped and **practical1**.map and outputs the PLINK binary files with prefix of **practical1_1**

> practical1_1.bed (binary file; not human readable)<br>
> practical1_1.bim <br>
> practical1_1.fam <br>

--> Paste and explain the log from the new VM here <---

The FAM file stores the pedigree information of the PED file (i.e. the first 6 columns).<br>
The BIM file is the extended MAP file with first four columns same as the MAP file. The 5th and 6th columns record the A1 and A2 alleles:<br>

> 5. Allele 1 (**A1**; corresponding to the **minor allele** by default; 0 is monomorphic)<br>
> 6. Allele 2 (**A2**; corresponding to the **major allele** by default)<br>

:closed_book: **Q:** What is the minor allele of the second SNP?
<details>
  <summary>You can try some basic unix commands by yourself first</summary>
  
  <pre>head -n 2 practical1_1.bim     # output the first 2 rows </pre>
  <pre>head -n 2 practical1_1.bim | tail -n 1 | cut -f 5   # output the first 2 rows -> output the last row -> cut out the 5th column </pre>
  <pre>sed -n 2p practical1_1.bim | cut -f 5 # output the 2nd row and cut out the 5th column </pre>
  <pre>awk 'NR==2 { print $5 }' practical1_1.bim   # output the 2nd row/record (NR==2) and print the 5th column ($5) </pre>
</details>

## Step 3: Data management in PLINK
In fact, the PLINK file records the genotype information for SNPs in three chromosomal regions
```bash
cat LDLgenes.bed
```
> 1 55505149 55530526 PCSK9<br>
> 2 21224301 21266945 APOB<br>
> 19 11200038 11244505 LDLR<br>

Next, we can try some basic data management functions in PLINK
#### - SNP management (Extract / Exclude)
##### -- Extract variants by by SNP ID(s)
E.g. A missense variant _APOB_:NM_000384.3:c.293C>T:p.Thr98Ile ([rs1367117](https://www.ebi.ac.uk/gwas/variants/rs1367117)) was previously reported to be associated with LDL level.<br>
:closed_book: **Q:** Who and how many of the samples carry at least one of the risk allele **A**?
```bash
plink --bfile practical1_1 --snp rs1367117 --recode --out practical1_1.rs1367117
```
<details>
  <summary>You can try some basic unix commands by yourself first</summary>
  
  <pre> awk '$7=="A" || $8=="A"' practical1_1.rs1367117.ped </pre>  
  <pre> awk '$7=="A" || $8=="A"' practical1_1.rs1367117.ped | wc </pre>
</details>

In addition to extract a single SNP, you can also extract multiple SNPs simultaneously using `--snps`<br>
E.g. `--snps rs1042034-rs1042031,rs693,exm175886`  for SNPs from rs1042034-rs1042031 as well as rs693 and exm175886. Similarly, you can exclude multiple SNPs simultaneously using `--exclude-snps`

##### --  Extract variants by chromosomal position
You may specify a chromosomal region and extract all genotypes fall within the region<br>
E.g. You can extract genotypes of all SNPs in _PCSK9_ (chr1:55505149-55530526)
```bash
plink --bfile practical1_1 --chr 1 --from-bp 55505149 --to-bp 55530526 --recode --out practical1_1.PCSK9_byChrPos
```
Instead of base pair, you are allowed to specific the position in kb using `--from-kb <kb pos> --to-kb <kb pos>` or in mb using `--from-mb <mb pos> --to-mb <mb pos>`<br>
To extract multiple regions, you can also use `--extract range <set file>` and specify a set file (in UCSC 1-based coordinate bed file format)
```bash
head -n 2 LDLgenes.bed > non-LDLR.bed
plink --bfile practical1_1 --extract range non-LDLR.bed --make-bed --out practical1_1.nonLDLR
```

##### -- Extracting or excluding multiple variants
You may specify a list of variants to be extracted or excluded<br>
```bash
awk '$1==1 { print $2 }' practical1_1.bim > PCSK9.snp
plink --bfile practical1_1 --extract PCSK9.snp --make-bed --out practical1_1.PCSK9_byExtract
```
```bash
awk '$1!=1 { print $2 }' practical1_1.bim > non-PCSK9.snp
plink --bfile practical1_1 --exclude non-PCSK9.snp --make-bed --out practical1_1.PCSK9_byExclude
```

#### - Sample management (Keep / Remove)
##### --  Keeping or removing samples
To obtain genotypes of a particular list of samples, you can use the `--keep` or `--remove` commands to include or exclude a list of samples with specified family and individual IDs, respectively.
```bash
awk '$1=="HG00103" || $1=="" { print $1,$2 }' practical1_1.fam > related.indiv
plink --bfile practical1_1 --remove related.indiv --make-bed --out practical1_1.unrelated
```

##### -- Basic summary statistics
Other functions of PLINK are implemented in a similar manner. For example, if you would like to obtain allele frequencies for SNPs, you can use the `--freq` function.
```bash
plink --bfile practical1_1 --freq --out practical1_1
```
You can combine all commands in one line. For example, to obtain frequency of SNPs in PCSK9 across all unrelated samples, you can use 
```bash
plink --bfile practical1_1 --remove related.indiv --extract PCSK9.snp --freq --out practical1_1.unrelated.PCSK9
```
:closed_book: **Q:** What are the minor allele frequencies of for all unrelated samples?
