# CAPP_SEQ_pipeline

Documentation of pipeline for processing of hybrid capture end reads of WGS data. Pipeline accepts raw FASTQ forward and reverse reads and performs duplex UMI consensus and/or single strand UMI consensus variant calling.

![image](https://user-images.githubusercontent.com/92883998/157540760-79000d43-a81b-4e36-8d7f-ac1d4b62aeb3.png)

# Required Packages

fastqc https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

picard https://github.com/broadinstitute/picard

bwa http://bio-bwa.sourceforge.net/bwa.shtml

trimmomatic http://www.usadellab.org/cms/?page=trimmomatic

samtools http://samtools.sourceforge.net/

GATK4: https://gatk.broadinstitute.org/hc/en-us

fgbio:http://fulcrumgenomics.github.io/fgbio/


# Usage

The pipeline perform duplex read consensus calling and single-strand consensus calling using the '--call_type' flag.
