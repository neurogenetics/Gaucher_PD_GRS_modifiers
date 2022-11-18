## GD modifiers

`LNG ❤️ Open Science 😍`

 - **Project:** Penetrance modifiers of PD in Gaucher patients
 - **Author(s):** Cornelis
 - **Others involved:** Dan V. Mike N.
 - **Externals:** Ellen Sidransky group, Manisha Balwani group, Tony Schapira group
 - **Date Last Updated:** September 2022

---
### Quick Description: 
GBA is an important gene for PD and both damaging missense variants (including N370S) and common non-coding variants have been associated with PD risk. Rare highly damaging missense variants often results in an increased of risk for PD >2-5 OR, while the common non-coding variant have the usual moderate risk increase of ~1.15. 
Homozyous or bi-allelic mutations in GBA cause Gaucher's disease which is an autosomal recessive lysosomal storage disease disorder. Gaucher's disease patients are also at higher risk of developing. The goal here is to identify modifiers of PD penetrance using a Gaucher's disease cohort.

### Motivation/Goals:
1) Understanding the underlying data and creating an overview of the data...
2) Generate covariate files
3) Generate genetic risk score files
4) Performing genetic risk score analysis

### Link to Manuscript:
TBD

## Structure of README:

### [1. Understanding data and subset data](#1-understanding-data-and-subset-data)
This section goes through:
- understanding the data and where is the data located

### [2. Generate covariate files](#2-Generate-covariate-files)
This section goes through:
- Remove carriers of other PD variants
- Subset for carriers of GBA variants
- Calculate PCs
- Merge PCs with phenotype info

### [3. Generate genetic risk score files](#3-Generate-genetic-risk-score-files)
This section goes through:
- Generating the genetic risk score

### [4. Performing genetic risk score analysis](#4-Performing-genetic-risk-score-analysis)
This section goes through:
- Performing the analysis


--------------------------------------


## 1. Understanding data and subset data
This section goes through:
- understanding the data and where is the data located

```
General working folder for end products:
cd /data/XXX/projects/GD_modifiers/

Location of genotype data:
/data/XXX/PD/GP2/genotypes/GD/
 clean 
 imputed

Location of raw data:
/data/XXX/PD/GP2/raw_genotypes/GD/QC AND /data/XXX/projects/GD_modifiers/RAW_DATA/

Location of phenotype data:
/data/XXX/projects/GD_modifiers/Clinical_data/
 sample_sheet_newyorkv2
 sample_sheet_nhgriv2.csv
 GD Sample Phenotypes - Cornelis.xlsx 

```

```
Sample overlap/relatedness and sex check

cd /data/XXX/projects/GD_modifiers/clean_plink_files_final/

module load plink

ls | grep fam | grep -v REMOVE | grep -v GD_EUR_AJ | sed -e 's/.fam//g' > merge_list.txt

plink --merge-list merge_list.txt --make-bed --out GD_clean_merged

grep -e MBA -e NHGR -e UCL GD_clean_merged.fam > GD_samples.txt

plink --bfile GD_clean_merged --keep GD_samples.txt --make-bed --out GD_clean_merged_cases_only

plink --bfile GD_clean_merged_cases_only --genome --min 0.2 --maf 0.05 --geno 0.05
## randomly removing one sample is PIHAT >0.2. Unless one of the samples is a GD-PD then the GD-control was excluded

plink --bfile GD --genome --min 0.2 --maf 0.05 --geno 0.05 --keep QC_PASSING.txt

module load plink
plink --bfile GD_clean_merged_cases_only --check-sex --maf 0.05 --geno 0.05 --out SEX_check

## all good!

```

```
Covariate inventory overview of QC passed samples (excluding relatedness ones above):

module load plink
plink --bfile GD_clean_merged_cases_only --remove REMOVE_related.txt --out analysis_ready_cohortv1 --make-bed

329 individuals

```
GD1 (in brackets is number of PD cases included)

| ANCESTRY | MTSINAI | NHGRI   | UCL    | TOTAL    |
|----------|---------|---------|--------|----------|
| AJ       | 83 (4)  | 66 (14) | 3      | 152      |
| EUR      | 9       | 42 (6)  | 23 (2) | 76       |
| AMR      | 1       | 5 (1)   | 0      | 6        |
| AAC      | 0       | 4       | 0      | 4        |
| EAS      | 0       | 1       | 0      | 1        |
| TOTAL    | 93      | 118     | 26     | 238      |


```
Create ancestry plot of QC passed samples

- make cleaned QC passed sample sheet
cd /data/XXX/projects/GD_modifiers/clean_plink_files_final

module load plink
plink --bfile analysis_ready_cohortv1 --maf 0.05 --geno 0.05 --hwe 1E-5 --make-bed --out GD_pre_PC

674310 variants and 329 people pass filters and QC.

plink --bfile GD_pre_PC --indep-pairwise 500 10 0.05 --autosome --out pruned_data
plink --bfile GD_pre_PC --extract pruned_data.prune.in --make-bed --out GD_pre_PC_prune 
plink --bfile GD_pre_PC_prune --pca --out PCA

module load R
R
library(ggplot2)
PC <- read.table("PCA.eigenvec",header=F)
ancestry <- read.table("../all_pheno.txt",header=T)
MM <- merge(PC,ancestry,by.x="V2",by.y="IID")

ggplot(MM, aes(x=V3, y=V4, color=ancestry))  +
    geom_point(size=2) + xlab("PC1") + ylab("PC2")

ggsave("ANCESTRY_breakdown_GD_QC_passed_march_2022.jpg")
# only GD1
newdata <- subset(MM, GD_Type == 1 )
ggplot(newdata, aes(x=V3, y=V4, color=ancestry))  +
    geom_point(size=2) + xlab("PC1") + ylab("PC2")

ggsave("ANCESTRY_breakdown_GD1_ONLY_QC_passed_march_2022.jpg")

```

## 2. Generate covariate files
This section goes through:
- Harmonize phenotype files
- Calculate PCs
- Merge PCs with phenotype info

Create covariate files:

```
all GD1's AJ only
all GD1's EUR only
all GD1's AJ-EUR combined

location GD (plink2):
cd /data/XXX/projects/GD_modifiers/clean_plink_files_final

awk '$4 == "AJ"' ../all_pheno.txt | awk '$8 == "1"' | cut -f 1,2 > AJ_Yes_analysis_GD1.txt
awk '$4 == "EUR"' ../all_pheno.txt | awk '$8 == "1"' | cut -f 1,2 > EUR_Yes_analysis_GD1.txt
cat AJ_Yes_analysis_GD1.txt EUR_Yes_analysis_GD1.txt > AJEUR_Yes_analysis_GD1.txt

module load plink

plink --bfile analysis_ready_cohortv1 --keep AJ_Yes_analysis_GD1.txt --make-bed --out AJ_Yes_analysis_GD1
plink --bfile analysis_ready_cohortv1 --keep EUR_Yes_analysis_GD1.txt --make-bed --out EUR_Yes_analysis_GD1
plink --bfile analysis_ready_cohortv1 --keep AJEUR_Yes_analysis_GD1.txt --make-bed --out AJEUR_Yes_analysis_GD1

plink --bfile AJ_Yes_analysis_GD1 --indep-pairwise 500 10 0.05 --autosome --out pruned_data --maf 0.05 --geno 0.05
plink --bfile AJ_Yes_analysis_GD1 --extract pruned_data.prune.in --make-bed --out GD_pre_PC_prune 
plink --bfile GD_pre_PC_prune --pca --out AJ_Yes_analysis_GD1_PCA

plink --bfile EUR_Yes_analysis_GD1 --indep-pairwise 500 10 0.05 --autosome --out pruned_data --maf 0.05 --geno 0.05
plink --bfile EUR_Yes_analysis_GD1 --extract pruned_data.prune.in --make-bed --out GD_pre_PC_prune 
plink --bfile GD_pre_PC_prune --pca --out EUR_Yes_analysis_GD1_PCA

plink --bfile AJEUR_Yes_analysis_GD1 --indep-pairwise 500 10 0.05 --autosome --out pruned_data --maf 0.05 --geno 0.05
plink --bfile AJEUR_Yes_analysis_GD1 --extract pruned_data.prune.in --make-bed --out GD_pre_PC_prune 
plink --bfile GD_pre_PC_prune --pca --out AJEUR_Yes_analysis_GD1_PCA

all GD1's with Coriell data AJ only
all GD1's with Coriell data EUR only
all GD1's with Coriell AJ-EUR combined

location Coriell (plink2):
/data/XXX/PD/GP2/releases/gp2_tier2/release1_29112021/raw_genotypes/

cd /data/XXX/PD/GP2/releases/gp2_tier2/release1_29112021/raw_genotypes/
module load plink/2.3-alpha

grep CORIELL_ GP2_round1_AJ_release1.samples > /data/XXX/projects/GD_modifiers/clean_plink_files_final/CORIELL_AJ.txt
grep CORIELL_ GP2_round1_EUR_release1.samples > /data/XXX/projects/GD_modifiers/clean_plink_files_final/CORIELL_EUR.txt

plink2 --pfile GP2_round1_AJ_release1 --keep /data/XXX/projects/GD_modifiers/clean_plink_files_final/CORIELL_AJ.txt --make-bed --out /data/XXX/projects/GD_modifiers/clean_plink_files_final/CORIELL_AJ
plink2 --pfile GP2_round1_EUR_release1 --keep /data/XXX/projects/GD_modifiers/clean_plink_files_final/CORIELL_EUR.txt --make-bed --out /data/XXX/projects/GD_modifiers/clean_plink_files_final/CORIELL_EUR

cd /data/XXX/projects/GD_modifiers/clean_plink_files_final/
module load plink
# AJ process
plink --bfile AJ_Yes_analysis_GD1 --bmerge CORIELL_AJ --out pass1 --make-bed
plink --bfile AJ_Yes_analysis_GD1 --exclude pass1-merge.missnp --make-bed --out temp
plink --bfile CORIELL_AJ --exclude pass1-merge.missnp --make-bed --out temp2
plink --bfile temp --bmerge temp2 --make-bed --out AJ_Yes_analysis_GD1_plus_coriell

plink --bfile AJ_Yes_analysis_GD1_plus_coriell --indep-pairwise 500 10 0.05 --autosome --out pruned_data --maf 0.05 --geno 0.05
plink --bfile AJ_Yes_analysis_GD1_plus_coriell --extract pruned_data.prune.in --make-bed --out GD_pre_PC_prune 
plink --bfile GD_pre_PC_prune --pca --out AJ_Yes_analysis_GD1_plus_coriell_PCA

# EUR process
plink --bfile EUR_Yes_analysis_GD1 --bmerge CORIELL_EUR --out pass1 --make-bed
plink --bfile EUR_Yes_analysis_GD1 --exclude pass1-merge.missnp --make-bed --out temp
plink --bfile CORIELL_EUR --exclude pass1-merge.missnp --make-bed --out temp2
plink --bfile temp --bmerge temp2 --make-bed --out EUR_Yes_analysis_GD1_plus_coriell

plink --bfile EUR_Yes_analysis_GD1_plus_coriell --indep-pairwise 500 10 0.05 --autosome --out pruned_data --maf 0.05 --geno 0.05
plink --bfile EUR_Yes_analysis_GD1_plus_coriell --extract pruned_data.prune.in --make-bed --out GD_pre_PC_prune 
plink --bfile GD_pre_PC_prune --pca --out EUR_Yes_analysis_GD1_plus_coriell_PCA

# AJ-EUR process
plink --bfile AJ_Yes_analysis_GD1_plus_coriell --bmerge EUR_Yes_analysis_GD1_plus_coriell --out AJEUR_Yes_analysis_GD1_plus_coriell --make-bed

plink --bfile AJEUR_Yes_analysis_GD1_plus_coriell --indep-pairwise 500 10 0.05 --autosome --out pruned_data --maf 0.05 --geno 0.05
plink --bfile AJEUR_Yes_analysis_GD1_plus_coriell --extract pruned_data.prune.in --make-bed --out GD_pre_PC_prune 
plink --bfile GD_pre_PC_prune --pca --out AJEUR_Yes_analysis_GD1_plus_coriell_PCA

### PC files:
module load R/3.6.1
R
require("plyr")
samples <- read.table("../all_pheno.txt",header=T)
PCAJ <- read.table("AJ_Yes_analysis_GD1_PCA.eigenvec",header=F)
PCAJ <- setNames(PCAJ, c("FID","IID","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"))
PCAJ <- PCAJ[ c(1:12) ]

PCEUR <- read.table("EUR_Yes_analysis_GD1_PCA.eigenvec",header=F)
PCEUR <- setNames(PCEUR, c("FID","IID","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"))
PCEUR <- PCEUR[ c(1:12) ]

PCAJEUR <- read.table("AJEUR_Yes_analysis_GD1_PCA.eigenvec",header=F)
PCAJEUR <- setNames(PCAJEUR, c("FID","IID","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"))
PCAJEUR <- PCAJEUR[ c(1:12) ]

PCAJC <- read.table("AJ_Yes_analysis_GD1_plus_coriell_PCA.eigenvec",header=F)
PCAJC <- setNames(PCAJC, c("FID","IID","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"))
PCAJC <- PCAJC[ c(1:12) ]

PCEURC <- read.table("EUR_Yes_analysis_GD1_plus_coriell_PCA.eigenvec",header=F)
PCEURC <- setNames(PCEURC, c("FID","IID","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"))
PCEURC <- PCEURC[ c(1:12) ]

PCAJEURC <- read.table("AJEUR_Yes_analysis_GD1_plus_coriell_PCA.eigenvec",header=F)
PCAJEURC <- setNames(PCAJEURC, c("FID","IID","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"))
PCAJEURC <- PCAJEURC[ c(1:12) ]

PCAJ2 <- merge(samples,PCAJ,by="IID")
PCAJ2$FID.y <- NULL
names(PCAJ2)[2] <- "FID"
names(PCAJ2)[1] <- "IID"
names(PCAJ2)[9] <- "PD_Status"
names(PCAJ2)[10] <- "PD_age_of_onset"
PCAJ22 <- PCAJ2[, c(2, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23)]

PCEUR2 <- merge(samples,PCEUR,by="IID")
PCEUR2$FID.y <- NULL
names(PCEUR2)[2] <- "FID"
names(PCEUR2)[1] <- "IID"
names(PCAJ2)[9] <- "PD_Status"
names(PCAJ2)[10] <- "PD_age_of_onset"
PCEUR22 <- PCEUR2[, c(2, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23)]

PCAJEUR2 <- merge(samples,PCAJEUR,by="IID")
PCAJEUR2$FID.y <- NULL
names(PCAJEUR2)[2] <- "FID"
names(PCAJEUR2)[1] <- "IID"
names(PCAJ2)[9] <- "PD_Status"
names(PCAJ2)[10] <- "PD_age_of_onset"
PCAJEUR22 <- PCAJEUR2[, c(2, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23)]

write.table(PCAJ22,file = "AJ_Yes_analysis_GD1_PCA_COV.txt", quote = FALSE,row.names=F,sep="\t")
write.table(PCEUR22,file = "EUR_Yes_analysis_GD1_PCA_COV.txt", quote = FALSE,row.names=F,sep="\t")
write.table(PCAJEUR22,file = "AJEUR_Yes_analysis_GD1_PCA_COV.txt", quote = FALSE,row.names=F,sep="\t")

# next wrap up adding in Coriell
CORIELL_AJ <- read.table("CORIELL_AJ.fam",header=F)
CORIELL_AJ <- setNames(CORIELL_AJ, c("FID","IID","trash1","trash2","sex","PD_PHENO"))
CORIELL_AJ$trash1 <- CORIELL_AJ$trash2 <- NULL
CORIELL_AJ$ANCESTRY <- "AJ"
CORIELL_AJ$COHORT <- "CORIELL"
CORIELL_AJ$PD_PHENO <- CORIELL_AJ$PD_PHENO-1
CORIELL_AJ[CORIELL_AJ == -10] <- NA

CORIELL_EUR <- read.table("CORIELL_EUR.fam",header=F)
CORIELL_EUR <- setNames(CORIELL_EUR, c("FID","IID","trash1","trash2","sex","PD_PHENO"))
CORIELL_EUR$trash1 <- CORIELL_EUR$trash2 <- NULL
CORIELL_EUR$ANCESTRY <- "EUR"
CORIELL_EUR$COHORT <- "CORIELL"
CORIELL_EUR$PD_PHENO <- CORIELL_EUR$PD_PHENO-1
CORIELL_EUR[CORIELL_EUR == -10] <- NA

CORIELL_AJ2 <- merge(CORIELL_AJ,PCAJC,by="IID")
CORIELL_AJ2 <- CORIELL_AJ2[c("FID.x", "IID", "sex", "PD_PHENO", "ANCESTRY", "COHORT", "PC1", "PC2", "PC3", "PC4", "PC5")]
names(CORIELL_AJ2)[1] <- "FID"

CORIELL_EUR2 <- merge(CORIELL_EUR,PCEURC,by="IID")
CORIELL_EUR2 <- CORIELL_EUR2[c("FID.x", "IID", "sex", "PD_PHENO", "ANCESTRY", "COHORT", "PC1", "PC2", "PC3", "PC4", "PC5")]
names(CORIELL_EUR2)[1] <- "FID"

CORIELL_AJEUR <- rbind(CORIELL_AJ,CORIELL_EUR)
CORIELL_AJEUR2 <- merge(CORIELL_AJEUR,PCAJEURC,by="IID")
CORIELL_AJEUR2 <- CORIELL_AJEUR2[c("FID.x", "IID", "sex", "PD_PHENO", "ANCESTRY", "COHORT", "PC1", "PC2", "PC3", "PC4", "PC5")]
names(CORIELL_AJEUR2)[1] <- "FID"

# combine with GD data...
PCAJ222 <- PCAJ22[, c(1, 2, 3, 9, 4, 12, 14, 15, 16, 17, 18)]
names(PCAJ222)[3] <- "sex"
names(PCAJ222)[4] <- "PD_PHENO"
names(PCAJ222)[5] <- "ANCESTRY"
names(PCAJ222)[6] <- "COHORT"

PCEUR222 <- PCEUR22[, c(1, 2, 3, 9, 4, 12, 14, 15, 16, 17, 18)]
names(PCEUR222)[3] <- "sex"
names(PCEUR222)[4] <- "PD_PHENO"
names(PCEUR222)[5] <- "ANCESTRY"
names(PCEUR222)[6] <- "COHORT"

PCAJEUR222 <- PCAJEUR22[, c(1, 2, 3, 9, 4, 12, 14, 15, 16, 17, 18)]
names(PCAJEUR222)[3] <- "sex"
names(PCAJEUR222)[4] <- "PD_PHENO"
names(PCAJEUR222)[5] <- "ANCESTRY"
names(PCAJEUR222)[6] <- "COHORT"

CORIELL_GD_AJ <- rbind(CORIELL_AJ2,PCAJ222)
CORIELL_GD_EUR <- rbind(CORIELL_EUR2,PCEUR222)
CORIELL_GD_AJEUR <- rbind(CORIELL_AJEUR2,PCAJEUR222)

write.table(CORIELL_GD_AJ,file = "AJ_Yes_analysis_GD1_PCA_COV_with_CORIELL.txt", quote = FALSE,row.names=F,sep="\t")
write.table(CORIELL_GD_EUR,file = "EUR_Yes_analysis_GD1_PCA_COV_with_CORIELL.txt", quote = FALSE,row.names=F,sep="\t")
write.table(CORIELL_GD_AJEUR,file = "AJEUR_Yes_analysis_GD1_PCA_COV_with_CORIELL.txt", quote = FALSE,row.names=F,sep="\t")

```

## 3 Generate genetic risk score files
This section goes through:
- Generating the genetic risk score

```
Using imputed clean data
### CORIELL
cd /data/XXX/PD/GP2/releases/gp2_tier2/release1_29112021/imputed_genotypes/EUR/
module load plink/2.3-alpha
ls | grep pvar | sed -e 's/.pvar//g' > merge_list.txt
# merge plink files
plink2 --pmerge-list merge_list.txt --make-pgen --out /data/XXX/projects/GD_modifiers/clean_plink_files_final/CORIEL_IMPUTE_EUR

cd /data/XXX/PD/GP2/releases/gp2_tier2/release1_29112021/imputed_genotypes/AJ/
module load plink/2.3-alpha
ls | grep pvar | sed -e 's/.pvar//g' > merge_list.txt
# merge plink files
plink2 --pmerge-list merge_list.txt --make-pgen --out /data/XXX/projects/GD_modifiers/clean_plink_files_final/CORIEL_IMPUTE_AJ

### GD EUR AJ data
cd /data/XXX/projects/GD_modifiers/
module load plink/2.3-alpha
/data/XXX/PD/GP2/genotypes/GD/imputed
module load plink/2.3-alpha
# unzip imputed files
for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
  do
  unzip -P imputer chr_"$chnum".zip
done
for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
  do
  plink2 --vcf chr"$chnum".dose.vcf.gz --make-pgen --out PLINK_CHR"$chnum"
done

ls | grep pvar | sed -e 's/.pvar//g' > merge_list.txt

plink2 --pmerge-list merge_list.txt --make-pgen --out GD_MERGED

scp GD_MERGED* /data/XXX/projects/GD_modifiers/clean_plink_files_final/

# generate HG38 TopMed ready risk score...

Take /data/XXX/GENERAL/META5_GRS_chr_bp.txt

Covert chr:bp to hg38 using liftover UCSC (hg19 => hg38)

add new hg38 coordinates to file and extract from PVAR

grep -f RISK_SCORE_VAR.txt MERGED.pvar > RISK_SCORE_VAR_PVAR.txt 

Then merge data, assume the same risk allele...

saved as /data/XXX/GENERAL/META5_GRS_TOPMED_style.txt

# generate PD genetic risk score
cd /data/XXX/projects/GD_modifiers/clean_plink_files_final/
module load plink/2.3-alpha

plink2 --pfile CORIEL_IMPUTE_EUR --extract /data/XXX/GENERAL/META5_GRS_TOPMED_style.txt --out CORIEL_IMPUTE_EUR_GRS --make-bed
plink2 --pfile CORIEL_IMPUTE_AJ --extract /data/XXX/GENERAL/META5_GRS_TOPMED_style.txt --out CORIEL_IMPUTE_AJ_GRS --make-bed

# manually update psam add UCL_ to UCL samples
sed -i 's/0_MBA_/MBA_/g' GD_MERGED.psam
sed -i 's/0_NHGRI_/NHGRI_/g' GD_MERGED.psam
sed -i 's/0_UCL_/UCL_/g' GD_MERGED.psam

plink2 --pfile GD_MERGED --extract /data/XXX/GENERAL/META5_GRS_TOPMED_style.txt --out IMPUTE_GD_MERGED_GRS --make-bed --keep analysis_ready_cohortv1.fam

ls | grep _GRS | grep -v var | grep bim | sed -e 's/.bim//g' > GRS_merge_list.txt
module load plink
plink --merge-list GRS_merge_list.txt --make-bed --out GRS_file_GD_CORIEL_AJ_EUR --geno 0.05

plink --bfile GRS_file_GD_CORIEL_AJ_EUR --score /data/XXX/GENERAL/META5_GRS_TOPMED_style.txt --out GRS_file_GD_CORIEL_AJ_EUR_PD_SCORE

plink --bfile GRS_file_GD_CORIEL_AJ_EUR --score /data/XXX/GENERAL/META5_GRS_TOPMED_style_no_GBA.txt --out GRS_file_GD_CORIEL_AJ_EUR_PD_SCORE_NO_GBA

plink --bfile GRS_file_GD_CORIEL_AJ_EUR --score /data/XXX/GENERAL/META5_GRS_TOPMED_style_no_GBA_no_LRRK2.txt --out GRS_file_GD_CORIEL_AJ_EUR_PD_SCORE_NO_GBA_NO_LRRK2

```

## 4 Performing genetic risk score analysis
This section goes through:
- Analysis groups
- Performing the analysis

#### Analysis groups

```
AJ GD1 non PD vs PD + covariates 5 PCs + sex
Euro GD1 non PD vs PD + covariates 5 PCs + sex
Euro + AJ GD1 non PD vs PD + covariates 5 PCs + sex

Then also adding in PD and controls...

```

#### Performing the analysis

```

cd /data/XXX/projects/GD_modifiers/clean_plink_files_final/

# score files
GRS_file_GD_CORIEL_AJ_EUR_PD_SCORE.profile
GRS_file_GD_CORIEL_AJ_EUR_PD_SCORE_NO_GBA.profile
GRS_file_GD_CORIEL_AJ_EUR_PD_SCORE_NO_GBA_NO_LRRK2.profile

# covariate files
AJ_Yes_analysis_GD1_PCA_COV.txt
EUR_Yes_analysis_GD1_PCA_COV.txt
AJEUR_Yes_analysis_GD1_PCA_COV.txt
AJ_Yes_analysis_GD1_PCA_COV_with_CORIELL.txt
EUR_Yes_analysis_GD1_PCA_COV_with_CORIELL.txt
AJEUR_Yes_analysis_GD1_PCA_COV_with_CORIELL.txt

# quickly fixing dementia 
sed -i 's/dementia/0/g' AJ_Yes_analysis_GD1_PCA_COV.txt
sed -i 's/dementia/0/g' EUR_Yes_analysis_GD1_PCA_COV.txt
sed -i 's/dementia/0/g' AJEUR_Yes_analysis_GD1_PCA_COV.txt
sed -i 's/dementia/0/g' AJ_Yes_analysis_GD1_PCA_COV_with_CORIELL.txt
sed -i 's/dementia/0/g' EUR_Yes_analysis_GD1_PCA_COV_with_CORIELL.txt
sed -i 's/dementia/0/g' AJEUR_Yes_analysis_GD1_PCA_COV_with_CORIELL.txt

# quicly fix header mistake
sed -i 's/PD_Status_/PD_Status/g' EUR_Yes_analysis_GD1_PCA_COV.txt
sed -i 's/PD_age_of_onset_/PD_age_of_onset/g' EUR_Yes_analysis_GD1_PCA_COV.txt
sed -i 's/PD_Status_/PD_Status/g' AJEUR_Yes_analysis_GD1_PCA_COV.txt
sed -i 's/PD_age_of_onset_/PD_age_of_onset/g' AJEUR_Yes_analysis_GD1_PCA_COV.txt

module load R/3.6.1
R
library(ggplot2)
# data load in

#####
# multiple ancestry...
#####

dataAJ <- read.table("AJ_Yes_analysis_GD1_PCA_COV_with_CORIELL.txt",header=T)
dataEUR <- read.table("EUR_Yes_analysis_GD1_PCA_COV_with_CORIELL.txt",header=T)
data <- read.table("AJEUR_Yes_analysis_GD1_PCA_COV_with_CORIELL.txt",header=T)
# GRS load in
GRS <- read.table("GRS_file_GD_CORIEL_AJ_EUR_PD_SCORE_NO_GBA.profile",header=T)
#GRS <- read.table("GRS_file_GD_CORIEL_AJ_EUR_PD_SCORE_NO_GBA_NO_LRRK2.profile",header=T)


MM_AJ <- merge(dataAJ,GRS,by="IID")
MM_AJcontrol <- subset(MM_AJ, COHORT == "CORIELL" & PHENO == "1" )
meanGRS <- mean(MM_AJcontrol$SCORE)
sdGRS <- sd(MM_AJcontrol$SCORE)
MM_AJ$SCOREZ <- (MM_AJ$SCORE - meanGRS)/sdGRS

MM_EUR <- merge(dataEUR,GRS,by="IID")
MM_EURcontrol <- subset(MM_EUR, COHORT == "CORIELL" & PHENO == "1" )
meanGRS <- mean(MM_EURcontrol$SCORE)
sdGRS <- sd(MM_EURcontrol$SCORE)
MM_EUR$SCOREZ <- (MM_EUR$SCORE - meanGRS)/sdGRS

MM_AJ_EUR <- merge(data,GRS,by="IID")
MM_AJEURcontrol <- subset(MM_AJ_EUR, COHORT == "CORIELL" & PHENO == "1" )
meanGRS <- mean(MM_AJEURcontrol$SCORE)
sdGRS <- sd(MM_AJEURcontrol$SCORE)
MM_AJ_EUR$SCOREZ <- (MM_AJ_EUR$SCORE - meanGRS)/sdGRS


# combine data
dataZ <- rbind(MM_AJ, MM_EUR)

dataZ$combo <- paste(dataZ$ANCESTRY, dataZ$COHORT,dataZ$PD_PHENO, sep = "_")
dataZ$combo[dataZ$combo=="AJ_CORIELL_0"] <- "AJ_control"  
dataZ$combo[dataZ$combo=="AJ_CORIELL_1"] <- "AJ_PD"  
dataZ$combo[dataZ$combo=="AJ_NY_0"] <- "AJ_GD"  
dataZ$combo[dataZ$combo=="AJ_NY_1"] <- "AJ_GD_PD"    
dataZ$combo[dataZ$combo=="AJ_NHGRI_0"] <- "AJ_GD" 
dataZ$combo[dataZ$combo=="AJ_NHGRI_1"] <- "AJ_GD_PD" 
dataZ$combo[dataZ$combo=="EUR_CORIELL_0"] <- "EUR_control" 
dataZ$combo[dataZ$combo=="EUR_CORIELL_1"] <- "EUR_PD" 
dataZ$combo[dataZ$combo=="EUR_NY_0"] <- "EUR_GD"
dataZ$combo[dataZ$combo=="EUR_NY_1"] <- "EUR_GD_PD"   
dataZ$combo[dataZ$combo=="EUR_NHGRI_0"] <- "EUR_GD" 
dataZ$combo[dataZ$combo=="EUR_NHGRI_1"] <- "EUR_GD_PD" 
dataZ$combo[dataZ$combo=="EUR_UCL_1"] <- "EUR_GD_PD"
dataZ$combo[dataZ$combo=="EUR_UCL_0"] <- "EUR_GD"
dataZ$combo[dataZ$combo=="AJ_UCL_1"] <- "AJ_GD_PD"
dataZ$combo[dataZ$combo=="AJ_UCL_0"] <- "AJ_GD"


## AJ_Yes_analysis_GD1_PCA_COV_with_CORIELL
table(dataZ$combo)
```
| group       | count |
|-------------|-------|
| AJ_control  | 109   |
| AJ_GD       | 134   |
| AJ_GD_PD    | 18    |
| AJ_PD       | 335   |
| EUR_control | 933   |
| EUR_GD      | 65    |
| EUR_GD_PD   | 8    |
| EUR_PD      | 2050  |

```
pdf("AJ_Yes_analysis_GD1_PCA_COV_with_CORIELL_normalized_control_V2.pdf",width=12)
p <- ggplot(dataZ, aes(x=as.factor(combo), y=SCOREZ, fill=as.factor(combo))) + 
  geom_violin(trim=FALSE)
p2 <- p+geom_boxplot(width=0.4, fill="white" ) + theme_minimal()
p2 + scale_fill_manual(values=c("azure1", "azure3","azure1", "azure3","azure1", "azure3","azure1", "azure3")) + theme_bw() + 
labs(title="PD GRS",x="control-case groups", y = "GRS Z-score") + theme(legend.position="none")
dev.off()

# save file for later use
write.table(dataZ, file="FINAL_FILE_GRS_INCLUDED_normalized_control.txt", quote=FALSE,row.names=F,sep="\t")


####
#AJ_Yes_analysis_GD1_PCA_COV_with_CORIELL.txt
####

AJ_only <- subset(dataZ, ANCESTRY =="AJ")

pdf("AJ_GRS_ONLY_PLOT.pdf",width=12)
p <- ggplot(AJ_only, aes(x=as.factor(combo), y=SCOREZ, fill=as.factor(combo))) + 
  geom_violin(trim=FALSE)
p2 <- p+geom_boxplot(width=0.4, fill="white" ) + theme_minimal()
p2 + scale_fill_manual(values=c("azure1", "azure3","azure1", "azure3")) + theme_bw() + 
labs(title="PD GRS",x="control-case groups", y = "GRS Z-score") + theme(legend.position="none")
dev.off()


pdf("AJ_GRS_ONLY_PLOT_grey.pdf",width=12)
p <- ggplot(AJ_only, aes(x=as.factor(combo), y=SCOREZ, fill=as.factor(combo))) + 
  geom_violin(trim=FALSE)
p2 <- p+geom_boxplot(width=0.4, fill="white" ) + theme_minimal()
p2 + scale_fill_manual(values=c("gray100", "gray75","gray50", "gray25")) + theme_bw() + 
labs(title="PD GRS",x="control-case groups", y = "GRS Z-score") + theme(legend.position="none")
dev.off()

GD_only <- subset(dataZ, COHORT == "NHGRI" | COHORT == "NY" | COHORT == "UCL")
GD_only_AJ <- subset(GD_only, ANCESTRY =="AJ")
table(GD_only_AJ$PD_PHENO)
# 18 cases and 134 controls

PD_only <- subset(MM_AJ, COHORT == "CORIELL")
table(PD_only$PD_PHENO)
# 335 cases and 109 controls

# NO AGE
thisFormula1 <- formula(paste("PD_PHENO ~ SCOREZ + sex + PC1 + PC2 + PC3 + PC4 + PC5"))
model1 <- glm(thisFormula1, data = GD_only_AJ, ,family=binomial)
print(summary(model1))

            Estimate Std. Error z value Pr(>|z|)  
(Intercept)  -2.1331     0.8309  -2.567   0.0102 *
SCOREZ        0.4104     0.2289   1.793   0.0730 .
sex          -0.1602     0.5383  -0.298   0.7660  
PC1          -3.9359     2.7539  -1.429   0.1529  
PC2           1.4538     3.8317   0.379   0.7044  
PC3           0.8794     3.6583   0.240   0.8100  
PC4          -7.6105     3.4258  -2.222   0.0263 *
PC5          -3.7701     3.2834  -1.148   0.2509  

# WITH AGE
thisFormula1 <- formula(paste("PD_PHENO ~ SCOREZ + sex + PC1 + PC2 + PC3 + PC4 + PC5 + AGE"))
model1 <- glm(thisFormula1, data = GD_only_AJ, ,family=binomial)
print(summary(model1))

            Estimate Std. Error z value Pr(>|z|)  
(Intercept) -0.61228    1.41405  -0.433   0.6650  
SCOREZ       0.44583    0.23737   1.878   0.0603 .
sex         -0.26602    0.54827  -0.485   0.6275  
PC1         -4.84006    2.87100  -1.686   0.0918 .
PC2          2.08739    4.05602   0.515   0.6068  
PC3          0.47716    3.68586   0.129   0.8970  
PC4         -7.34858    3.46109  -2.123   0.0337 *
PC5         -3.61220    3.35557  -1.076   0.2817  
AGE         -0.02548    0.01981  -1.286   0.1983  

# NO AGE
thisFormula1 <- formula(paste("PD_PHENO ~ SCOREZ + sex + PC1 + PC2 + PC3 + PC4 + PC5"))
model1 <- glm(thisFormula1, data = PD_only, ,family=binomial)
print(summary(model1))

            Estimate Std. Error z value Pr(>|z|)    
(Intercept)  2.45482    0.37960   6.467 1.00e-10 ***
SCOREZ       0.56451    0.10452   5.401 6.62e-08 ***
sex         -1.11337    0.24187  -4.603 4.16e-06 ***
PC1         -2.17485    3.12796  -0.695    0.487    
PC2         -0.08652    4.10115  -0.021    0.983    
PC3         -5.30393    4.65617  -1.139    0.255    
PC4          8.95688    6.86249   1.305    0.192    
PC5         -4.04094    6.34813  -0.637    0.524    

# WITH AGE
thisFormula1 <- formula(paste("PD_PHENO ~ SCOREZ + sex + PC1 + PC2 + PC3 + PC4 + PC5 + AGE"))
model1 <- glm(thisFormula1, data = PD_only, ,family=binomial)
print(summary(model1))

            Estimate Std. Error z value Pr(>|z|)    
(Intercept)   9.86460    1.11601   8.839  < 2e-16 ***
SCOREZ        0.54290    0.11897   4.563 5.04e-06 ***
sex          -1.34096    0.28135  -4.766 1.88e-06 ***
PC1          -0.01470    4.24542  -0.003    0.997    
PC2          -3.22680    8.15324  -0.396    0.692    
PC3          -7.04231    6.14172  -1.147    0.252    
PC4           7.73286    5.67885   1.362    0.173    
PC5         -10.26529    8.10112  -1.267    0.205    
AGE          -0.11150    0.01435  -7.771 7.80e-15 ***


####
#EUR_Yes_analysis_GD1_PCA_COV_with_CORIELL.txt
####

EUR_only <- subset(dataZ, ANCESTRY =="EUR")

pdf("EUR_GRS_ONLY_PLOT.pdf",width=12)
p <- ggplot(EUR_only, aes(x=as.factor(combo), y=SCOREZ, fill=as.factor(combo))) + 
  geom_violin(trim=FALSE)
p2 <- p+geom_boxplot(width=0.4, fill="white" ) + theme_minimal()
p2 + scale_fill_manual(values=c("azure1", "azure3","azure1", "azure3")) + theme_bw() + 
labs(title="PD GRS",x="control-case groups", y = "GRS Z-score") + theme(legend.position="none")
dev.off()

pdf("EUR_GRS_ONLY_PLOT_gray.pdf",width=12)
p <- ggplot(EUR_only, aes(x=as.factor(combo), y=SCOREZ, fill=as.factor(combo))) + 
  geom_violin(trim=FALSE)
p2 <- p+geom_boxplot(width=0.4, fill="white" ) + theme_minimal()
p2 + scale_fill_manual(values=c("gray100", "gray75","gray50", "gray25")) + theme_bw() + 
labs(title="PD GRS",x="control-case groups", y = "GRS Z-score") + theme(legend.position="none")
dev.off()

GD_only <- subset(EUR_only, COHORT == "NHGRI" | COHORT == "NY" | COHORT == "UCL")
table(GD_only$PD_PHENO)
# 8 cases and 65 controls

pdf("TESTINGv2_EUR_GD_onlyv2.pdf",width=12)
p <- ggplot(GD_only, aes(x=as.factor(combo), y=SCOREZ, fill=as.factor(combo))) + 
  geom_violin(trim=FALSE)
p2 <- p+geom_boxplot(width=0.4, fill="white" ) + theme_minimal()
p2 + scale_fill_manual(values=c("azure1", "azure3")) + theme_bw() + 
labs(title="PD GRS",x="control-case groups", y = "GRS Z-score") + theme(legend.position="none")
dev.off()

GD_only <- subset(MM_EUR, COHORT == "NHGRI" | COHORT == "NY" | COHORT == "UCL")
table(GD_only$PD_PHENO)
# 8 cases and 65 controls

PD_only <- subset(MM_EUR, COHORT == "CORIELL")
table(PD_only$PD_PHENO)
# 2050 cases and 933 controls

# NO AGE
thisFormula1 <- formula(paste("PD_PHENO ~ SCOREZ + sex + PC1 + PC2 + PC3 + PC4 + PC5"))
model1 <- glm(thisFormula1, data = GD_only, ,family=binomial)
print(summary(model1))

            Estimate Std. Error z value Pr(>|z|)  
(Intercept)  -3.4257     1.6557  -2.069   0.0385 *
SCOREZ        0.9450     0.5570   1.696   0.0898 .
sex           0.2945     0.9404   0.313   0.7541  
PC1          -1.1983     6.5535  -0.183   0.8549  
PC2           6.8849     4.9662   1.386   0.1656  
PC3          -3.1165     4.2205  -0.738   0.4603  
PC4          -7.0512     3.8512  -1.831   0.0671 .
PC5           3.1989     5.3502   0.598   0.5499  

# WITH AGE
thisFormula1 <- formula(paste("PD_PHENO ~ SCOREZ + sex + PC1 + PC2 + PC3 + PC4 + PC5 + AGE"))
model1 <- glm(thisFormula1, data = GD_only, ,family=binomial)
print(summary(model1))

            Estimate Std. Error z value Pr(>|z|)  
(Intercept) -2.997239   2.383835  -1.257   0.2086  
SCOREZ       0.952018   0.560632   1.698   0.0895 .
sex          0.209423   1.006443   0.208   0.8352  
PC1         -1.292685   6.646555  -0.194   0.8458  
PC2          7.245435   5.230261   1.385   0.1660  
PC3         -3.329378   4.321230  -0.770   0.4410  
PC4         -7.322375   4.019884  -1.822   0.0685 .
PC5          3.523594   5.561762   0.634   0.5264  
AGE         -0.006377   0.025672  -0.248   0.8038  

# NO AGE
thisFormula1 <- formula(paste("PD_PHENO ~ SCOREZ + sex + PC1 + PC2 + PC3 + PC4 + PC5"))
model1 <- glm(thisFormula1,data = PD_only, ,family=binomial)
print(summary(model1))

            Estimate Std. Error z value Pr(>|z|)    
(Intercept)  1.65640    0.12889  12.851   <2e-16 ***
SCOREZ       0.52315    0.04224  12.385   <2e-16 ***
sex         -0.69014    0.08268  -8.347   <2e-16 ***
PC1          2.44333    2.31224   1.057    0.291    
PC2         -5.41943    2.22697  -2.434    0.015 *  
PC3         -3.71880    3.91527  -0.950    0.342    
PC4          1.96848    2.76194   0.713    0.476    
PC5          0.55881    2.33084   0.240    0.811    

coef(summary(model1))[,4]
# real P-value => 3.154674e-35

# WITH AGE
thisFormula1 <- formula(paste("PD_PHENO ~ SCOREZ + sex + PC1 + PC2 + PC3 + PC4 + PC5 + AGE"))
model1 <- glm(thisFormula1,data = PD_only, ,family=binomial)
print(summary(model1))

            Estimate Std. Error z value Pr(>|z|)    
(Intercept)  8.192232   0.354405  23.115  < 2e-16 ***
SCOREZ       0.439785   0.047736   9.213  < 2e-16 ***
sex         -0.888227   0.095186  -9.331  < 2e-16 ***
PC1          7.691461   2.678266   2.872 0.004081 ** 
PC2         -9.603514   2.567720  -3.740 0.000184 ***
PC3         -4.470774   3.449473  -1.296 0.194950    
PC4          1.974919   3.294412   0.599 0.548856    
PC5         -4.617924   2.898108  -1.593 0.111064    
AGE         -0.102437   0.004765 -21.498  < 2e-16 ***

coef(summary(model1))[,4]
# real P-value => 3.174363e-20


#### next up is forest plot

### AJ - GD
thisFormula1 <- formula(paste("PD_PHENO ~ SCOREZ + sex + PC1 + PC2 + PC3 + PC4 + PC5"))
model1 <- glm(thisFormula1, data = AJ_GD_only, ,family=binomial)
print(summary(model1))
            Estimate Std. Error z value Pr(>|z|)    
SCOREZ        0.4104     0.2289   1.793   0.0730 .

### EUR - GD
thisFormula1 <- formula(paste("PD_PHENO ~ SCOREZ + sex + PC1 + PC2 + PC3 + PC4 + PC5"))
model1 <- glm(thisFormula1, data = EUR_GD_only, ,family=binomial)
print(summary(model1))
            Estimate Std. Error z value Pr(>|z|)    
SCOREZ        0.9450     0.5570   1.696   0.0898 .

GD (per ancestry separate, GRS normalized using Coriell data)
            Estimate Std. Error z value Pr(>|z|)    
AJ         0.4104     0.2289   1.793   0.0730  
EUR        0.9450     0.5570   1.696   0.0898

# input GD-PD NO AGE
Ancestry	Beta	SE	Z	P	
AJ	0.4104	0.2289	1.793	0.073	
EUR	0.9450	0.5570	1.696	0.0898	

# input GD-PD WITH AGE
Ancestry	Beta	SE	Z	P
AJ	0.44583	0.23737	1.878	0.0603
EUR	0.952018	0.560632	1.698	0.0895


PD (per ancestry separate, GRS normalized using Coriell data)
            Estimate Std. Error z value Pr(>|z|)    
EUR       0.53072    0.03903  13.598 3.154674e-35
AJ        0.56451    0.10452   5.401 6.62e-08

# input PD NO AGE
Ancestry	Beta	SE	Z	P
AJ	0.56451	0.10452	5.401	6.62E-08
EUR	0.53072	0.03903	13.598	3.15E-35

# input PD WITH AGE
Ancestry	Beta	SE	Z	P
AJ	0.54290	0.11897	4.563	5.04E-06
EUR	0.439785	0.047736	9.213	3.17E-20


### forest plotting

R
require("rmeta")
library(metafor)
library(data.table)
# pick one:
data <- read.table("input_forest_GD.txt", header = T)
data <- read.table("input_forest_PD.txt", header = T)
data <- read.table("input_forest_GD_AGE.txt", header = T)
data <- read.table("input_forest_PD_AGE.txt", header = T)

labs <- data$Ancestry
yi   <- data$Beta
sei  <- data$SE
resFe  <- rma(yi=yi, sei=sei, method="FE")
resRe  <- rma(yi=yi, sei=sei)
print(summary(resFe))
print(summary(resRe))

##### GD NO AGE
Model Results:
Fixed-Effects Model
estimate      se    zval    pval   ci.lb   ci.ub    
  0.4876  0.2117  2.3032  0.0213  0.0727  0.9026  * 

Random-Effects Model
estimate      se    zval    pval   ci.lb   ci.ub   ​ 
  0.4876  0.2117  2.3032  0.0213  0.0727  0.9026  * 

##### GD WITH AGE
Model Results:
Fixed-Effects Model
estimate      se    zval    pval   ci.lb   ci.ub    
  0.5228  0.2186  2.3916  0.0168  0.0944  0.9512  * 

Random-Effects Model
estimate      se    zval    pval   ci.lb   ci.ub   ​ 
  0.5228  0.2186  2.3916  0.0168  0.0944  0.9512  * 

# random and fixed is the same, which is great :)

pdf(file = "GD_only_forest_FINALv2_AGE.pdf", width = 8, height = 6)
forest(resRe, xlim=c(resFe$beta-0.9931472,resFe$beta+0.6931472), main="Meta-analysis of PD GRS",atransf=exp, xlab=paste("Odds Ratio (95%CI)",sep=""), slab=labs, mlab="Fixed Effects", col = "red", border = "red", cex=.9)
dev.off()

pdf(file = "GD_only_forest_FINAL_arrowv2_AGE.pdf", width = 8, height = 6)
forest(resRe, xlim=c(resFe$beta-0.9931472,resFe$beta+0.9931472), main="Meta-analysis of PD GRS", at=log(c(0.5,0.75, 1, 2, 14)), atransf=exp, xlab=paste("Odds Ratio (95%CI)",sep=""), slab=labs, mlab="Random Effects", col = "red", border = "red", cex=.9)
dev.off()

##### PD

Model Results:
Fixed-Effects Model
estimate      se     zval    pval   ci.lb   ci.ub
 0.4541  0.0443  10.2496  <.0001  0.3673  0.5409  *** 

Random-Effects Model
estimate      se     zval    pval   ci.lb   ci.ub 
  0.4541  0.0443  10.2496  <.0001  0.3673  0.5409  *** 

# random and fixed is the same, which is great :)

resRe$pval
1.188681e-24

pdf(file = "PD_only_forest_FINAL_order_OK_v2_AGE.pdf", width = 8, height = 6)
forest(resRe, xlim=c(resFe$beta-0.9931472,resFe$beta+0.6931472), main="Meta-analysis of PD GRS",atransf=exp, xlab=paste("Odds Ratio (95%CI)",sep=""), slab=labs, mlab="Fixed Effects", col = "red", border = "red", cex=.9)
dev.off()

pdf(file = "PD_only_forest_FINAL_arrow_order_OK_v2_AGE.pdf", width = 8, height = 6)
forest(resRe, xlim=c(resFe$beta-0.9931472,resFe$beta+0.9931472), main="Meta-analysis of PD GRS", at=log(c(0.5,0.75, 1, 2, 12)), atransf=exp, xlab=paste("Odds Ratio (95%CI)",sep=""), slab=labs, mlab="Random Effects", col = "red", border = "red", cex=.9)
dev.off()

## Done....
```
![myImage](https://media.giphy.com/media/Wj7lNjMNDxSmc/giphy.gif)

