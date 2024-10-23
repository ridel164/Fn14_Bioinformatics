# TNFRSF12A_Bioinformatics-
This repository is for storing bioinformatic data as part of my MSc Genetics Thesis. Included within is commands for downloading Whole Genome Sequence (WGS) SRA files, through to variant calling and annotation.
All work outlined in this repository belongs to The Cunliffe Laboratory at The University of Otago. 

## Outline of Project
The premise for the successive bioinformatic work follows Aim 1 of my MSc Genetics Research, which is to characterize regualtory factors within the TNFRSF12A promoter region. Further, the work focuses on identification of promoter sequence abberations consistent with _TNFRSF12A_ transcription and expression. 

Objective: Analyse breast cancer cell line sequences and compare variants between cell lines with high and low expression of _TNFRSF12A_, to analyse whether sequence abberations contribute to differential expression of _TNFRSF12A_ in breast cancer. 

Breast cacner cell lines MDAMB231 and BT549 have been chosen as high expressing cell lines where as, MCF7 and T47D have been chosen as low expressing. These decisions have been made in line with published literature, and are outlined in (0.0) Preliminary Results.

### Project Setup and Directories
The Cunliffe laboratory high capacity storage (HCS) space has been used as the main directory  for outpu tdata and storage of commands/scripts. Additional HCS spaces within the University of Otago have been used for access to modules, the reference genome, and designed annotation script. 

All scripts and commands were written in BBEdit and are stored in the `dsm/molecular-oncology` HCS space. The scripts run using an `sbatch` command had parameters outlined as described below where relvenat, and were finalized using `exit 0`. 

## 1. SRA Download

### 1.1. Load Modules
`module load SRA-Toolkit`

### 1.2 Set Working Directories 
```
SRAID=$1                                        
wrkdr=/mnt/hcs/dsm-molecularoncology/Elyse_Bioinformatics
```
Change directory to the working directory described above using `cd $wrkdr`

### 1.3. Downloading SRA files using prefetch
```
if [[ ! -s ${wrkdr}/${SRAID}/${SRAID}.sra ]]
then
	prefetch --output-directory $wrkdr ${SRAID}	--max-size 1000000000
fi
```

## 2. Creating FASTQ Files

Use `mkdir -p  $wrkdr/${SRAID}/FASTQ_Files` to create fastq file directory 

### 2.1 Splitting SRA into FASTQ Files 

```
if [[ ! -s ${wrkdr}/${SRAID}/${SRAID}.sra ]]
then 
	fastq-dump --outdir $wrkdr/${SRAID}/FASTQ_Files --split-files ${wrkdr}/${SRAID}/${SRAID}.sra
fi
```
## 3. Align samples to Reference sequence (hg38)

Use `mkdir -p  $wrkdr/${SRAID}/Alignments` to create output alignment directory 

### 3.1 Load Modules
```
module load BWA
module load SAMtools
```

### 3.2. Set Working Directories
The whole human genome sequence hg38 is used as the refernece sequence. This is stored within a shared high capacity storage space at The University of Otago.

`ref="/resource/bundles/hg38/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"`

### 3.3. Create bam file to be aligned to refernce.

This command uses the fastq files created in 2.1 to develop a bam file, which is aligned to the reference sequence described in 3.2. The aligned sequences are then piped sirectly into samtools for sorting. 
```
if [[ ! -s $wrkdr/${SRAID}/Alignments/$SRAID.bam ]]
then
	bwa mem -t 12 -k 21 -Y $ref $wrkdr/${SRAID}/FASTQ_Files/${SRAID}_1.fastq $wrkdr/${SRAID}/FASTQ_Files/${SRAID}_2.fastq | \
	samtools sort -@ 4 -o $wrkdr/${SRAID}/Alignments/$SRAID.bam 
fi
```

## 4. Variant call using DeepVariant
The following script for section 4 was processed as an entire script using `sbatch` with the outlined parameters below
```
#!/bin/bash

#SBATCH --job-name=DeepVariant
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=24:00:00
```

### 4.1 Set Inputs
The below inputs are required to specify directory locations 
```
SRAID=$1
wrkdr=/mnt/hcs/dsm-molecularoncology/Elyse_Bioinformatics
ref="/resource/bundles/hg38/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
bam_file=$wrkdr/${SRAID}/Alignments/${SRAID}.bam
target_regions=chr16:3017000-3022000
sample_id=${SRAID}
```

The target region spans upstream of the TNFRSF12A promoter (within upstream gene CLDN6) and beyond exon 1.


### 4.2 Run DeepVariant 
```
/resource/pipelines/Fastq2VCFplus/DeepVariant.sh -i $bam_file -s $sample_id -a $target_regions -k config_DV_hg38.sh
```
This command requires manual changing of the SRAID dependant on which bam file is being processed. 

## 5. Annotation of vcf

### 5.1 Load Modules 
```
module load BCFtools/1.17-GCC-12.2.0
```

### 5.2. Set Working Directories 

Path to the annotation script (opens help menu)
`/resource/pipelines/Variant_Annotation/Annotate_my_VCF.sh`

Changing working directory to match individual vcf being annotated in 5.3
`cd /mnt/hcs/dsm-molecularoncology/Elyse_Bioinformatics/SRR8670674/Variant_Call/`

### 5.3. Running Annotation Script 

To run annotation script for each individual vcf the command below is used. Note that, the working directory and vcf location must be adjusted each time the scrpt is submitted. 
```
/resource/pipelines/Variant_Annotation/Annotate_my_VCF.sh -i /mnt/hcs/dsm-molecularoncology/Elyse_Bioinformatics/SRR8670674/Variant_Call/SRR8670674.vcf.gz -d /resource/pipelines/Variant_Annotation/defaults_hg38_gnomADv4.txt
```

### 5.4 Exporting data from annotated vcf 

For both commands below the correct file name must be adjusted in the script each time. Here, `${SRR8670674}.vcf.gz` is used, and must be adjusted for differnt files names (incl vcf siffix).  

To generate human readable data from the annotated files, the command below is used to convert vcf to tsv. 
```
    vcf_human_readable () {
          input=$1
           out=$(echo ${input} | sed "s/\.vcf.gz/_out\.tsv/g")
info=; for var in `zcat ${input} | grep "##INFO" | cut -d\, -f1 | sed "s/##INFO=<ID=//g" | grep -v "-"`; do if [ -z "$info" ]; then info="%INFO/${var}"; else info="${info}\t%INFO/${var}"; fi; done;
bcftools query --print-header -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t'$info'\n' ${input} > ${out}
          wait < <(jobs -p)
}
 
vcf_human_readable ${SRR8670674}.vcf.gz
```

To exract further sample information the script below is used.
```
     vcf_human_readable () {
          input=$1
           out=$(echo ${input} | sed "s/\.vcf.gz/_out\.tsv/g")
format=; for var in `less ${input} | grep "^##FORMAT" | cut -d\, -f1 | sed "s/##FORMAT=<ID=//g"`; do if [ -z "$format" ]; then format="%${var}"; else format="${format}\t%${var}"; fi; done;
info=; for var in `zcat ${input} | grep "##INFO" | cut -d\, -f1 | sed "s/##INFO=<ID=//g" | grep -v "-"`; do if [ -z "$info" ]; then info="%INFO/${var}"; else info="${info}\t%INFO/${var}"; fi; done;
bcftools query --print-header -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t'$info'\t['$format'\t]\n' ${input} > ${out}  
          wait < <(jobs -p)
}
 
vcf_human_readable ${SRR8670674}.vcf.gz
```
------------------------------------------------------------------------------------------------------------------









