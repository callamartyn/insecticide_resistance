#!/bin/bash
# testing how to connect with github

optspec=":r:i:h:s:t"

MYDIR="$(dirname "$(realpath "$0")")"

snpfinder="${MYDIR}/snpfinder.py"
echo $snpfinder
SNPtable="${MYDIR}/all_vgsc_snps.csv"



trimmomatic="/home/calla/bin/Trimmomatic-0.38/trimmomatic-0.38.jar"
adapters="/home/calla/bin/Trimmomatic-0.38/adapters/TruSeq3-PE.fa"
freebayes="/home/calla/bin/freebayes/bin/freebayes"
minimap="/home/calla/bin/minimap2-2.9_x64-linux/minimap2"
samtools="/home/calla/bin/samtools"
snpfinder="snpfinder.py"
scriptname=${0##*/}

while getopts "$optspec" option; do
	case "${option}" in
		r) ref=${OPTARG};; # choose a reference sequence
		i) input=${OPTARG};; # choose input files
    s) suffix=${OPTARG};; # specify fastq suffix (to be stripped for file name)
		t) SNPtable=${OPTARG};; # specify table containing snps of interest
		h) HELP=1;;
	esac
done

if [[ $HELP -eq 1  ||  $# -lt 1 ]]
then
        cat <<USAGE

${bold}${scriptname}${normal}

This program will identify SNPs in raw sequencing data

It requires paired end data

If you have more than one species/reference, you will need to run each reference separately

${bold}Command Line Switches:${normal}

        -h      Show this help & ignore all other switches

        -r      Specify reference sequence in fasta format, required

        -i      Specify input fastq file; specify R1 file only, program will identify R2 file

This switch is optional and will use all fastq files in directory if unspecified

        -s      Specify suffix for R1 files

This switch is optional and will be set to "_R1_001.fastq.gz" if unspecified

	-t Specify a SNP table containing SNPs of interest

This switch is optional and will be set to "all_vgsc_snps.csv" if unspecified

${bold}Usage:${normal}

        Run example 1
                $scriptname -r JN695777.1.fasta

        Run example 2
                $scriptname -r AM422833.1.fasta -i Sample_23_RNA_lib_R1.fastq -s "_R1.fastq"

USAGE
	exit
fi

if [[ -z $suffix ]]
then
	suffix="_R1_001.fastq.gz"
fi

if [[ -z $input ]]
then
	input=$(find . -maxdepth 1 -name "*${suffix}")
fi
echo "fastq files are $input"

if [[ -z $SNPtable ]]
then
	SNPtable="all_vgsc_snps.csv"
fi
echo "using table $SNPtable"

for fastq in $input
do
        name=$(basename $fastq $suffix)
        #check if file has already been trimmed and paired
        if [ ! -f "${name}_R2_paired.fastq.gz" ]; then
          # if there is no paired file, continue to trimmomatic
          echo "trimming $name ..."
          # cut adapters and quality trim the file
          java -jar $trimmomatic PE \
      		  -threads 16 \
      		    -phred33 -trimlog "${name}.trim.log" -basein $fastq "${name}_R1_paired.fastq.gz" "${name}_R1_unpaired.fastq.gz" "${name}_R2_paired.fastq.gz" "${name}_R2_unpaired.fastq.gz" ILLUMINACLIP:$adapters:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75
        fi
        # continue if file exists or once trimmomatic has run
        # check if sam file already exists
        if [ ! -f "${name}.sam" ]; then
          echo "mapping reads for $name ..."
          # map the trimmed reads to the reference using minimap2
          $minimap -ax sr $ref "${name}_R1_paired.fastq.gz" "${name}_R2_paired.fastq.gz" > "${name}.sam"
        fi
        # convert bam to sam and remove unmapped readlines
        if [ ! -f "${name}_sort.bam" ]; then
          echo "editing bam file for $name"
          samtools view -F4 -b "${name}.sam" > "${name}_mapped.bam"
          samtools rmdup "${name}_mapped.bam" "${name}_dedup.bam"
          samtools fixmate -r -O bam "${name}_dedup.bam" "${name}_clean.bam"
          # sort the sam file
          samtools sort "${name}_clean.bam" -o "${name}_sort.bam" -O bam
        fi
        if [ ! -f "${name}.vcf" ]; then
          # generate a vcf using freebayes
          echo "generating vcf for $name"
          $freebayes -f $ref "${name}_sort.bam" > "${name}.vcf"
        fi
done

files=$(find . -maxdepth 1 -name "*.vcf")

python "${MYDIR}/snpfinder.py" $ref $SNPtable $files
