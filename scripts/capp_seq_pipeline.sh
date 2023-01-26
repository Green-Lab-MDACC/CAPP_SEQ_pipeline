#!/bin/bash

#USER INPUTS
#SAVE FASTQ FILES UNDER FASTQ FOLDER WITH FILE NAME MATCHING SAMPLE_ID_R1 AND SAMPLE_ID_R2
#CONSENSUS: "DUPLEX" FOR DUPLEX CONSENSUS CALLING. "SINGLE" FOR SINGLE STRAND CONSENSUS CALLING. (CREATES A NEW FOLDER NAMED THIS UNDER SAMPLE NAME)


#SPECIFY PROJECT FOLDER
PROJECT="INSERT_FOLDER_NAME"
#SPECIFY WORKING DIRECTORY
WORKING_DIR=$PWD/bnsugg/duplex_pipeline
#SPECIFY VARIANT CALLING TYPE
CONSENSUS="DUPLEX" #"DUPLEX" or "SINGLE"
SINGLE_CONSENSUS=3 #NUMBER OF CONSENSUS READS REQUIRED FOR SINGLE STRAND CALL #http://fulcrumgenomics.github.io/fgbio/tools/latest/CallMolecularConsensusReads.html
DOUBLE_CONSENSUS=6 #NUMBER OF CONSENSUS READS REQUIRED FOR DUPLEX CALL #http://fulcrumgenomics.github.io/fgbio/tools/latest/CallDuplexConsensusReads.html
MIN_MAP_QUAL=30



#SPECIFY NAMES OF SAMPLES (USED FOR FASTQ/BAM/VCF LABELING AND FINDING FILE PATHS) 
for i in 'INSERT_SAMPLE_NAME' \
'INSERT_SAMPLE_NAME' \
'INSERT_SAMPLE_NAME' \
'INSERT_SAMPLE_NAME'

do

    SAMPLE_NAME=$i

    #SAMPLES FILE NAMES
    SAMPLE_ID_R1="${SAMPLE_NAME}_R1"
    SAMPLE_ID_R2="${SAMPLE_NAME}_R2"

    #ASSIGN INPUT FILE PATHS FOR FASTQS
    FIRST_FASTQ=$WORKING_DIR/data/$PROJECT/FASTQ/${SAMPLE_ID_R1}.fastq.gz
    SECOND_FASTQ=$WORKING_DIR/data/$PROJECT/FASTQ/${SAMPLE_ID_R2}.fastq.gz

    #CREATE DIRECTORY FOR SAVING OUTPUTS
    mkdir $WORKING_DIR/data/$PROJECT/$SAMPLE_NAME
    mkdir $WORKING_DIR/data/$PROJECT/$SAMPLE_NAME/$CONSENSUS

    #ASSIGN OUTPUT FILE PATHS AND NAMES
    FGBIO_FASTQTOBAM_OUT=$WORKING_DIR/data/$PROJECT/$SAMPLE_NAME/$CONSENSUS/${SAMPLE_NAME}.ubam.bam
    FGBIO_FASTQTOBAM_EXTRACTED_OUT=$WORKING_DIR/data/$PROJECT/$SAMPLE_NAME/$CONSENSUS/${SAMPLE_NAME}.ubam.umi.extracted.bam
    PICARD_SAMTOFASTQ_OUT=$WORKING_DIR/data/$PROJECT/$SAMPLE_NAME/$CONSENSUS/${SAMPLE_NAME}.ubam.umi.extracted.fastq
    BWA_TEMP_SAM=$WORKING_DIR/data/$PROJECT/$SAMPLE_NAME/$CONSENSUS/${SAMPLE_NAME}.ubam.aligned.sam
    BWA_TEMP_SAM_2=$WORKING_DIR/data/$PROJECT/$SAMPLE_NAME/$CONSENSUS/${SAMPLE_NAME}.ubam.aligned_2.sam
    PICARD_MERGEBAMALIGNMENT_OUT=$WORKING_DIR/data/$PROJECT/$SAMPLE_NAME/$CONSENSUS/${SAMPLE_NAME}.ubam.umi.extracted.aligned.bam
    FGBIO_GROUPED_OUT=$WORKING_DIR/data/$PROJECT/$SAMPLE_NAME/$CONSENSUS/${SAMPLE_NAME}.ubam.umi.extracted.aligned.grouped.bam
    FGBIO_GROUPED_OUT_METRICS=$WORKING_DIR/data/$PROJECT/$SAMPLE_NAME/$CONSENSUS/${SAMPLE_NAME}.ubam.umi.extracted.aligned.grouped.txt
    UMI_CALLS_FAMILY_HISTOGRAM=$WORKING_DIR/data/$PROJECT/$SAMPLE_NAME/$CONSENSUS/${SAMPLE_NAME}.umi_calls.histogram.csv
    FGBIO_CALLDUPLEXCONSENSUS_GROUPED_OUT=$WORKING_DIR/data/$PROJECT/$SAMPLE_NAME/$CONSENSUS/${SAMPLE_NAME}.ubam.umi.extracted.aligned.grouped.umicalled.bam
    FGBIO_CALLDUPLEXCONSENSUS_GROUPED_OUT_METRICS=$WORKING_DIR/data/$PROJECT/$SAMPLE_NAME/$CONSENSUS/${SAMPLE_NAME}.ubam.umi.extracted.aligned.grouped.umicalled.metrics.txt
    BAITS=$WORKING_DIR/data/Baits_Merged.interval_list
    TARGETS=$WORKING_DIR/data/Targets_Merged.interval_list
    PICARD_SAMTOFASTQ_CALLDUPLEXCONSENSUS_OUT=$WORKING_DIR/data/$PROJECT/$SAMPLE_NAME/$CONSENSUS/${SAMPLE_NAME}.ubam.umi.extracted.aligned.grouped.umicalled.fastq
    PICARD_SAMTOFASTQ_CALLDUPLEXCONSENSUS_MAPPED_OUT=$WORKING_DIR/data/$PROJECT/$SAMPLE_NAME/$CONSENSUS/${SAMPLE_NAME}.ubam.umi.extracted.aligned.grouped.umicalled.mapped.bam
    GATK_REALIGNTARGETS=$WORKING_DIR/data/$PROJECT/$SAMPLE_NAME/$CONSENSUS/${SAMPLE_NAME}.ubam.umi.extracted.aligned.grouped.umicalled.mapped.realigntargets.intervals
    GATK_REALIGNEDBAM=$WORKING_DIR/data/$PROJECT/$SAMPLE_NAME/$CONSENSUS/${SAMPLE_NAME}.ubam.umi.extracted.aligned.grouped.umicalled.mapped.realigned.bam
    GATK_REALIGNTARGETSFILT=$WORKING_DIR/data/$PROJECT/$SAMPLE_NAME/$CONSENSUS/${SAMPLE_NAME}.ubam.umi.extracted.aligned.grouped.umicalled.mapped.realigned.filtered.bam
    QUALITY_DIST=$WORKING_DIR/data/$PROJECT/$SAMPLE_NAME/$CONSENSUS/${SAMPLE_NAME}.ubam.umi.extracted.aligned.grouped.umicalled.mapped.realigned.BQSR.quality.txt
    QUALITY_DISTPDF=$WORKING_DIR/data/$PROJECT/$SAMPLE_NAME/$CONSENSUS/${SAMPLE_NAME}.ubam.umi.extracted.aligned.grouped.umicalled.mapped.realigned.BQSR.quality.pdf
    QUALITY_DIST2=$WORKING_DIR/data/$PROJECT/$SAMPLE_NAME/$CONSENSUS/${SAMPLE_NAME}.ubam.umi.extracted.aligned.grouped.umicalled.mapped.realigned.BQSR.quality2.txt
    QUALITY_DISTPDF2=$WORKING_DIR/data/$PROJECT/$SAMPLE_NAME/$CONSENSUS/${SAMPLE_NAME}.ubam.umi.extracted.aligned.grouped.umicalled.mapped.realigned.BQSR.quality2.pdf
    FGBIO_CALLDUPLEXCONSENSUS_GROUPED_OUT_FILTERED_METRICS=$WORKING_DIR/data/$PROJECT/$SAMPLE_NAME/$CONSENSUS/${SAMPLE_NAME}.ubam.umi.extracted.aligned.grouped.umicalled.mapped.realigned.filtered.txt
    GATK_REALIGNTARGETSFILTCLIPPED=$WORKING_DIR/data/$PROJECT/$SAMPLE_NAME/$CONSENSUS/${SAMPLE_NAME}.ubam.umi.extracted.aligned.grouped.umicalled.mapped.realigned.filtered.clipped.bam
    FGBIO_CALLDUPLEXCONSENSUS_GROUPED_OUT_FILTERED_CLIPPED_METRICS=$WORKING_DIR/data/$PROJECT/$SAMPLE_NAME/$CONSENSUS/${SAMPLE_NAME}.ubam.umi.extracted.aligned.grouped.umicalled.mapped.realigned.filtered.clipped.txt
    VARIANTS_CALLED=$WORKING_DIR/data/$PROJECT/$SAMPLE_NAME/$CONSENSUS/${SAMPLE_NAME}.variants.vcf.gz
    VARIANTS_CALLED_TABLE=$WORKING_DIR/data/$PROJECT/$SAMPLE_NAME/$CONSENSUS/${SAMPLE_NAME}.variants.table
    READS_PRINTED=$WORKING_DIR/data/$PROJECT/$SAMPLE_NAME/$CONSENSUS/${SAMPLE_NAME}.READS_PRINTED.txt
    FINAL_CALLS=$WORKING_DIR/data/$PROJECT/$SAMPLE_NAME/$CONSENSUS/${SAMPLE_NAME}.haplo.vcf.gz
    FINAL_CALLS_TABLE=$WORKING_DIR/data/$PROJECT/$SAMPLE_NAME/$CONSENSUS/${SAMPLE_NAME}.haplo.table
    PONPATH=$WORKING_DIR/data/PON.xlsx
    ENUM_VCF=$WORKING_DIR/data/$PROJECT/$SAMPLE_NAME/$CONSENSUS/${SAMPLE_NAME}.enum.${CONSENSUS}.vcf
    FINAL_ENUM_TABLE=$WORKING_DIR/data/$PROJECT/$SAMPLE_NAME/$CONSENSUS/${SAMPLE_NAME}.enum.filt.${CONSENSUS}.table
    ENUM_FILT_VCF=$WORKING_DIR/data/$PROJECT/$SAMPLE_NAME/$CONSENSUS/${SAMPLE_NAME}.enum.filt.${CONSENSUS}.vcf
    FASTQC_OUT=$WORKING_DIR/data/$PROJECT/$SAMPLE_NAME
    VARDICT_CALLS=$WORKING_DIR/data/$PROJECT/$SAMPLE_NAME/$CONSENSUS/${SAMPLE_NAME}.vardict.table
    BQSR_TABLE=$WORKING_DIR/data/$PROJECT/$SAMPLE_NAME/$CONSENSUS/${SAMPLE_NAME}.ubam.umi.extracted.aligned.grouped.umicalled.mapped.realigned.filtered.clipped.bqsr.table
    BQSR_BAM=$WORKING_DIR/data/$PROJECT/$SAMPLE_NAME/$CONSENSUS/${SAMPLE_NAME}.ubam.umi.extracted.aligned.grouped.umicalled.mapped.realigned.filtered.clipped.bqsr.bam
    VARDICT_CALLS_TABLE=$WORKING_DIR/data/$PROJECT/$SAMPLE_NAME/$CONSENSUS/${SAMPLE_NAME}_vardict.table

    #ASSIGN REFERENCE PATHS
    BWA_INDEX=$WORKING_DIR/data/hg38/bwa/hg38
    REFGENOME=$WORKING_DIR/data/hg38/Homo_sapiens_assembly38.fasta
    REFGENOMEDICT=$WORKING_DIR/data/hg38/Homo_sapiens_assembly38.dict
    MILLS_INDEL=$WORKING_DIR/data/hg38/hg38_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
    ONEKG_INDEL=$WORKING_DIR/data/hg38/beta/Homo_sapiens_assembly38.known_indels.vcf.gz
    DBSNP=$WORKING_DIR/data/hg38/hg38_bundle/dbsnp_138.hg38.vcf.gz
    BAITS_INTERVAL=$WORKING_DIR/data/$PROJECT/Probes_Merged_hg38.interval_list
    TARGETS_INTERVAL=$WORKING_DIR/data/$PROJECT/Targets_Merged_hg38.interval_list
    BAITS=$WORKING_DIR/data/$PROJECT/Probes_Merged_hg38.bed
    TARGETS=$WORKING_DIR/data/$PROJECT/Targets_Merged_hg38.bed

    #MODULE PATHS
    FGBIO=/risapps/rhel7/fgbio/1.3.0/fgbio-1.3.0.jar
    PICARD=$PWD/bnsugg/programs/picard.jar
    GATK=/risapps/rhel7/gatk/4.1.0.0/gatk-package-4.1.0.0-local.jar
    VARDICT=$PWD/bnsugg/programs/VarDictJava/build/install/VarDict/bin
    #SET READGROUPS AND STRUCTURE (USED FOR UMI DETERMINATION)
    RG_PU=`zcat $FIRST_FASTQ | head -n 1 | awk -F ":" '{ OFS="."; print $3, $4, $10; }' `
    READGROUP="@RG\tPL:ILLUMINA\tID:${RG_PU}\tPU:${RG_PU}\tSM:${SAMPLE_NAME}"

    #SPECICY READ STRUCTURE BASED ON UMI STRUCTURE
    READ_STRUCTURE='8M1S+T 8M1S+T' #https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures

    #RUN FASTQC
    fastqc $FIRST_FASTQ -o $FASTQC_OUT
    fastqc $SECOND_FASTQ -o $FASTQC_OUT

    # #RUN FASTQ TO UBAM
    if test -f "$SECOND_FASTQ"; then
        java -jar $FGBIO FastqToBam --input $FIRST_FASTQ $SECOND_FASTQ \
        --output $FGBIO_FASTQTOBAM_OUT --read-structures $READ_STRUCTURE --sort true --read-group-id $RG_PU --sample $SAMPLE_NAME --platform-unit $RG_PU --library $SAMPLE_NAME
    fi

    # #EXTRACT UMIs
    if test -f "$FGBIO_FASTQTOBAM_OUT"; then
        java -jar $FGBIO ExtractUmisFromBam -i $FGBIO_FASTQTOBAM_OUT -o $FGBIO_FASTQTOBAM_EXTRACTED_OUT -r $READ_STRUCTURE -t ZA ZB -s RX
    fi

    # #CONVERT UBAM TO INTERLEAVE FASTQ
    if test -f "$FGBIO_FASTQTOBAM_EXTRACTED_OUT"; then
    java -jar $PICARD SamToFastq -I $FGBIO_FASTQTOBAM_EXTRACTED_OUT -F $PICARD_SAMTOFASTQ_OUT -INTER true
    fi

    #ALIGN TO HG38 AND MERGE
    if test -f "$PICARD_SAMTOFASTQ_OUT"; then
        bwa mem -t 4 -M -w 2 -p -R $READGROUP $BWA_INDEX $PICARD_SAMTOFASTQ_OUT > $BWA_TEMP_SAM
    fi

    if test -f "$BWA_TEMP_SAM"; then
        java -jar $PICARD MergeBamAlignment -UNMAPPED $FGBIO_FASTQTOBAM_OUT -ALIGNED $BWA_TEMP_SAM -O $PICARD_MERGEBAMALIGNMENT_OUT -R $REFGENOME -SO coordinate --ALIGNER_PROPER_PAIR_FLAGS true -MAX_GAPS -1 -ORIENTATIONS FR -VALIDATION_STRINGENCY SILENT -CREATE_INDEX true
    fi


    #GROUPED PAIRED READS BY UMI
    if test -f "$PICARD_MERGEBAMALIGNMENT_OUT"; then
        if [ $CONSENSUS = "DUPLEX" ]; then
            java -jar $FGBIO GroupReadsByUmi -i $PICARD_MERGEBAMALIGNMENT_OUT -o $FGBIO_GROUPED_OUT -s paired --edits 0 --min-map-q $MIN_MAP_QUAL -f $UMI_CALLS_FAMILY_HISTOGRAM
        else
            java -jar $FGBIO GroupReadsByUmi -i $PICARD_MERGEBAMALIGNMENT_OUT -o $FGBIO_GROUPED_OUT -s adjacency --edits 0 --min-map-q $MIN_MAP_QUAL -f $UMI_CALLS_FAMILY_HISTOGRAM
        fi
    fi


    #CALL CONSENSUS READS
    if test -f "$FGBIO_GROUPED_OUT"; then
        if [ $CONSENSUS = "DUPLEX" ]; then
            java -jar $FGBIO CallDuplexConsensusReads -i $FGBIO_GROUPED_OUT -o $FGBIO_CALLDUPLEXCONSENSUS_GROUPED_OUT --error-rate-pre-umi 45 --error-rate-post-umi 40 --min-input-base-quality $MIN_MAP_QUAL --min-reads $DOUBLE_CONSENSUS $SINGLE_CONSENSUS $SINGLE_CONSENSUS
        else
            java -jar $FGBIO CallMolecularConsensusReads -i $FGBIO_GROUPED_OUT -o $FGBIO_CALLDUPLEXCONSENSUS_GROUPED_OUT --error-rate-pre-umi 30 --min-reads $SINGLE_CONSENSUS --error-rate-post-umi 30 --min-input-base-quality $MIN_MAP_QUAL
        fi
    fi

    #RE-ALIGN
    if test -f "$FGBIO_CALLDUPLEXCONSENSUS_GROUPED_OUT"; then
        java -jar $PICARD SamToFastq -I $FGBIO_CALLDUPLEXCONSENSUS_GROUPED_OUT -F $PICARD_SAMTOFASTQ_CALLDUPLEXCONSENSUS_OUT -INTER true
    fi

    if test -f "$PICARD_SAMTOFASTQ_CALLDUPLEXCONSENSUS_OUT"; then
        bwa mem -t 4 -p $BWA_INDEX $PICARD_SAMTOFASTQ_CALLDUPLEXCONSENSUS_OUT > $BWA_TEMP_SAM_2
    fi

    if test -f "$BWA_TEMP_SAM_2"; then
        java -jar $PICARD MergeBamAlignment -UNMAPPED $FGBIO_CALLDUPLEXCONSENSUS_GROUPED_OUT -ALIGNED $BWA_TEMP_SAM_2 -O $PICARD_SAMTOFASTQ_CALLDUPLEXCONSENSUS_MAPPED_OUT -R $REFGENOME -SO coordinate --ALIGNER_PROPER_PAIR_FLAGS true -MAX_GAPS -1 -ORIENTATIONS FR -VALIDATION_STRINGENCY SILENT -CREATE_INDEX true

    fi

    #CREATE QUALITY SCORE DISTRIBUTION
    if test -f "$PICARD_SAMTOFASTQ_CALLDUPLEXCONSENSUS_MAPPED_OUT"; then
        java -jar $PICARD QualityScoreDistribution -I $PICARD_SAMTOFASTQ_CALLDUPLEXCONSENSUS_MAPPED_OUT -O $QUALITY_DIST -CHART $QUALITY_DISTPDF
    fi

    #FILTER READS
    if test -f "$PICARD_SAMTOFASTQ_CALLDUPLEXCONSENSUS_MAPPED_OUT"; then
        if [ $CONSENSUS = "DUPLEX" ]; then
            java -jar $FGBIO FilterConsensusReads --input $PICARD_SAMTOFASTQ_CALLDUPLEXCONSENSUS_MAPPED_OUT --output $GATK_REALIGNTARGETSFILT --ref $REFGENOME --min-reads $DOUBLE_CONSENSUS $SINGLE_CONSENSUS $SINGLE_CONSENSUS --max-read-error-rate 0.025 --max-base-error-rate .1 --min-base-quality $MIN_MAP_QUAL	 
        else
            java -jar $FGBIO FilterConsensusReads --input $PICARD_SAMTOFASTQ_CALLDUPLEXCONSENSUS_MAPPED_OUT --output $GATK_REALIGNTARGETSFILT --ref $REFGENOME -M $SINGLE_CONSENSUS --max-read-error-rate 0.025 --max-base-error-rate .1 --min-base-quality $MIN_MAP_QUAL
        fi
    fi

    #CLIP BAM
    if test -f "$GATK_REALIGNTARGETSFILT"; then
        java -jar $FGBIO ClipBam --input $GATK_REALIGNTARGETSFILT --output $GATK_REALIGNTARGETSFILTCLIPPED --ref $REFGENOME --clip-overlapping-reads true
        
        #REMOVE TEMP FILES
        rm -f "$FGBIO_CALLDUPLEXCONSENSUS_GROUPED_OUT"
        rm -f "$BWA_TEMP_SAM"
        rm -f "$BWA_TEMP_SAM_2"
        rm -f "$PICARD_SAMTOFASTQ_CALLDUPLEXCONSENSUS_OUT"
        rm -f "$FGBIO_GROUPED_OUT"
        rm -f "$PICARD_MERGEBAMALIGNMENT_OUT"
        rm -f "$FGBIO_FASTQTOBAM_EXTRACTED_OUT"
        rm -f "$PICARD_SAMTOFASTQ_OUT"
        rm -f "$BWA_TEMP_SAM"
        rm -f "$FGBIO_FASTQTOBAM_OUT"
        rm -f "$FGBIO_FASTQTOBAM_EXTRACTED_OUT"
        rm -f "$PICARD_SAMTOFASTQ_OUT"

    fi

    #BQSR
    if test -f "$GATK_REALIGNTARGETSFILTCLIPPED"; then
        java -jar $GATK BaseRecalibrator -I $GATK_REALIGNTARGETSFILTCLIPPED -R $REFGENOME --known-sites $MILLS_INDEL --known-sites $DBSNP --known-sites $ONEKG_INDEL -O $BQSR_TABLE
    fi

    if test -f "$BQSR_TABLE"; then
        java -jar $GATK ApplyBQSR -R $REFGENOME -I $GATK_REALIGNTARGETSFILTCLIPPED --bqsr-recal-file $BQSR_TABLE  -O $BQSR_BAM
    fi


    #COLLECT HYBRID CAPTURE METRICS
    if test -f "$BQSR_BAM"; then
        java -jar $PICARD CollectHsMetrics I=$BQSR_BAM O=$FGBIO_CALLDUPLEXCONSENSUS_GROUPED_OUT_FILTERED_CLIPPED_METRICS BAIT_INTERVALS=$BAITS_INTERVAL TARGET_INTERVALS=$TARGETS_INTERVAL R=$REFGENOME
    fi

    #CALL VARIANTS USING HAPLOTYPE CALLER AND MUTECT
    if test -f "$BQSR_BAM"; then
        java -jar $GATK Mutect2 -I $BQSR_BAM -R $REFGENOME -O $VARIANTS_CALLED --max-reads-per-alignment-start 0 
        java -jar $GATK HaplotypeCaller  -I $BQSR_BAM -R $REFGENOME -O $FINAL_CALLS
        AF_THR="0"
        $VARDICT/VarDict -b $BQSR_BAM -G $REFGENOME -N $SAMPLE_NAME -h -f $AF_THR -c 1 -S 2 -E 3 $BAITS > $VARDICT_CALLS
    fi


    if test -f "$VARIANTS_CALLED"; then
        java -jar $GATK VariantsToTable \
        -V $VARIANTS_CALLED \
        -F CHROM -F POS -F TYPE -F ALT -F REF -GF AD -GF DP \
        -O $VARIANTS_CALLED_TABLE
    fi

    if test -f "$FINAL_CALLS"; then
        java -jar $GATK VariantsToTable \
        -V $FINAL_CALLS \
        -F CHROM -F POS -F TYPE -F ALT -F REF -GF AD -GF DP \
        -O $FINAL_CALLS_TABLE
    fi


    if test -f "$BQSR_BAM"; then
        bcftools mpileup -f $REFGENOME --max-depth 0 $BQSR_BAM -a "FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP" | bcftools call -A -m > $ENUM_VCF
    fi

    if test -f "$ENUM_VCF"; then
        java -jar $GATK VariantsToTable -V $ENUM_VCF -F CHROM -F POS -F TYPE -F ALT -F REF -F QUAL -F DP -F DP4 -GF AD -O $FINAL_ENUM_TABLE
    fi

    # #FILTER AND POST PROCESS CALL VARIANTS
    if [ -f "$VARIANTS_CALLED_TABLE" ] && [ -f "$FINAL_CALLS_TABLE" ] ; then
        python $WORKING_DIR/scripts/seq_post_process.py ${SAMPLE_NAME} ${VARIANTS_CALLED_TABLE} ${FINAL_CALLS_TABLE} ${PONPATH} ${PROJECT} ${CONSENSUS}
    fi


done
