# CAPP_SEQ_pipeline

Documentation of pipeline for processing of hybrid capture end reads of WGS data. Pipeline accepts raw FASTQ forward and reverse reads and performs duplex UMI consensus and/or single strand UMI consensus variant calling.

![image](https://user-images.githubusercontent.com/92883998/157908176-ae0d0d86-72b1-4748-ac74-dd3d2ea2d08f.png)


# Required Packages

fastqc https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

picard https://github.com/broadinstitute/picard

bwa http://bio-bwa.sourceforge.net/bwa.shtml

trimmomatic http://www.usadellab.org/cms/?page=trimmomatic

samtools http://samtools.sourceforge.net/

GATK4: https://gatk.broadinstitute.org/hc/en-us

fgbio:http://fulcrumgenomics.github.io/fgbio/


# Usage

Inputs: (modify inside shell script)

``` Sample_ID``` 

``` Panel_ID``` 

```call_type```  Options = ```Duplex``` or ```Single``` Default = ```Duplex```

```min_read``` Default = ```1```

```READ_STRUCTURE``` Default = ```9M+T 9M+T```




Required Files:

``` ${Sample_ID}_R1.fastq.gz```  raw R1 FASTQ saved under "data/FASTQ"
``` ${Sample_ID}_R2.fastq.gz```  raw R2 FASTQ saved under "data/FASTQ"

``` ${Panel_ID}.baits.hs38DH.interval_list```   BAIT interval file saved under "data/panel"

``` ${Panel_ID}.targets.hs38DH.interval_list```   Target interval file saved under "data/panel"

```REFGENOME``` path to ucsc Homo_sapiens_assembly38.fasta

```BWA_INDEX``` path to hg38 bwa index

Optional:  ```PoN.csv``` pipeline output for PoNs saved under "data/PoN" 

The pipeline perform duplex read consensus calling and single-strand consensus calling using the ```call_type``` flag.

The # of single strand consensus reads required for dupelx calling can also be specificed using the ```min_read``` flag. Default = 1 see http://fulcrumgenomics.github.io/fgbio/tools/latest/CallDuplexConsensusReads.html

Read structure must be specified for extracting and grouping UMIs with ```READ_STRUCTURE``` flag. For more on read structure, see https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures


CAPP SEQ Metrics are computed using Picard ```CollectHsMetrics``` output to directory in csv format.

Variants are called using GATK ```HaplotypeCaller``` and ```Mutect2```

Variants combined and post-processed filtered in Python script with an optional PoNs filter and output in csv format.

