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

Inputs:

``` Sample_ID``` 

``` Panel_ID``` 

```call_type```  Options = ```Duplex``` or ```Single```

```min_read```   Options = ```0``` or ```1```

```READ_STRUCTURE```

Required Files:

``` ${Sample_ID}_R1.fastq.gz```  raw R1 FASTQ saved under "data/FASTQ"
``` ${Sample_ID}_R2.fastq.gz```  raw R2 FASTQ saved under "data/FASTQ"

``` ${Panel_ID}.baits.hs38DH.interval_list```   BAIT interval file saved under "data/panel"

``` ${Panel_ID}.targets.hs38DH.interval_list```   Target interval file saved under "data/panel"

```REFGENOME``` path to ucsc Homo_sapiens_assembly38.fasta

```BWA_INDEX``` path to hg38 bwa index

The pipeline perform duplex read consensus calling and single-strand consensus calling using the ```--call_type``` flag.

The # of single strand consensus reads required for dupelx calling can also be specificed using the ```--min_read``` flag. Default = 1

Read strucutre must be specified for extracting and grouping UMIs. For more on read structure, see https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures

CAPP SEQ Metrics are computed using Picard ```CollectHsMetrics```

Variants are called using GATK ```HaplotypeCaller``` and ```Mutect2```

