# TNFRSF12A_Bioinformatics-
This repository is for storing bioinformatic data as part of my MSc Genetics Thesis. Included within is commands for downloading WGS SRA files through to Variant calling and Annotating them.  
All work outlined in this repository is part of the Cunliffe Laboratory at The University of Otago. 

## Outline of Project
The purpose of this aim is to analyse breast cancer cell line sequences, and compare variants between cell lines with high and low expression of TNFRSF12A. Breast cacner cell lines MDAMB231 and BT549 have been chosen as high expressing cell lines where as, MCF7 and T47D have been chosen as low expressing. These decision have been made in line with published literature. 

### Project Setup and Directories
The Cunliffe laboratory high capacity storage (HCS) space has been used as the main directory  for outputdata and storage of scripts. Additional HCS spaces within the University of Otago have been used for access to modules, the reference genome, and designed annotation script. 

## SRA Download

### Setting Working Directories 
`SRAID=$1`                                         
`wrkdr=/mnt/hcs/dsm-molecularoncology/Elyse_Bioinformatics`
`cd $wrkdr`

### Downloading SRA files using prefetch
`module load SRA-Toolkit`
```
if [[ ! -s ${wrkdr}/${SRAID}/${SRAID}.sra ]]
then
	prefetch --output-directory $wrkdr ${SRAID}	--max-size 1000000000
fi
```

## Creating FASTQ Files
To create fastq file directory 
`mkdir -p  $wrkdr/${SRAID}/FASTQ_Files`

### Splitting SRA into FASTQ Files 
```
if [[ ! -s ${wrkdr}/${SRAID}/${SRAID}.sra ]]
then 
	fastq-dump --outdir $wrkdr/${SRAID}/FASTQ_Files --split-files ${wrkdr}/${SRAID}/${SRAID}.sra
fi
```
### Alingmnet of Sampls to Reference Sequence (hg38)
To create alignment directory 
`mkdir -p  $wrkdr/${SRAID}/Alignments`

### Creating bam file to be aligned to reference sequence and sorted
`module load BWA`
`module load SAMtools`
`ref="/resource/bundles/hg38/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"`
```
if [[ ! -s $wrkdr/${SRAID}/Alignments/$SRAID.bam ]]
then
	bwa mem -t 12 -k 21 -Y $ref $wrkdr/${SRAID}/FASTQ_Files/${SRAID}_1.fastq $wrkdr/${SRAID}/FASTQ_Files/${SRAID}_2.fastq | \
	samtools sort -@ 4 -o $wrkdr/${SRAID}/Alignments/$SRAID.bam 
fi
```

