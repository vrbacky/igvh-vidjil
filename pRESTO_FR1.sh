#!/bin/bash

FILES=../data/*L001_R1*.fastq
PRIMERS1=./Primers_FR1.fasta
PRIMERS2=./Primers_JH.fasta
QUALITY=20
NPROC=8 

for r1_file in $FILES
do
    no_path=$(basename $r1_file)
    dir_name=${no_path%_S2?_L001_R1_001.fastq}
    r2_file=${r1_file%1_001.fastq}"2_001.fastq"

    echo "***** Processing sample: $dir_name *****"

    # Create final directory
    mkdir $dir_name
    cd $dir_name

    mkdir qc

    # Assemble pairs
    AssemblePairs.py align -1 ../$r1_file -2 ../$r2_file --coord illumina --rc tail --outdir ./ --outname 01-${dir_name}_AP --log 01-${dir_name}_AP.log | tee 00-log.txt
    ParseLog.py -l 01-${dir_name}_AP.log -f ID LENGTH OVERLAP ERROR PVALUE | tee -a 00-log.txt
    fastqc -t $NPROC -q -o ./qc 01-${dir_name}_AP_assemble-pass.fastq

    # Filter quality
    FilterSeq.py quality -s 01-${dir_name}_AP_assemble-pass.fastq -q $QUALITY --outname 02-${dir_name}_AP_FQ --log 02-${dir_name}_FQ.log | tee -a 00-log.txt 
    ParseLog.py -l 02-${dir_name}_FQ.log -f ID LENGTH OVERLAP ERROR PVALUE | tee -a 00-log.txt
    fastqc -t $NPROC -q -o ./qc 02-${dir_name}_AP_FQ_quality-pass.fastq

    # Mask primers
    MaskPrimers.py score -s 02-${dir_name}_AP_FQ_quality-pass.fastq -p $PRIMERS1 --start 0 --mode cut --outname 03-${dir_name}_AP_FQ_MPF --log 03-${dir_name}_MPF.log | tee -a 00-log.txt
    MaskPrimers.py score -s 03-${dir_name}_AP_FQ_MPF_primers-pass.fastq -p $PRIMERS2 --revpr --start 0 --mode cut --outname 03-${dir_name}_AP_FQ_MPFR --log 03-${dir_name}_MPF.log | tee -a 00-log.txt
    ParseHeaders.py expand -s 03-${dir_name}_AP_FQ_MPFR_primers-pass.fastq -f PRIMER | tee -a 00-log.txt
    ParseHeaders.py rename -s 03-${dir_name}_AP_FQ_MPFR_primers-pass_reheader.fastq -f PRIMER1 PRIMER2 -k VPRIMER CPRIMER --outname 03-${dir_name}_AP_FQ_MPFR_for_Vidjil | tee -a 00-log.txt
    fastq_to_fasta -n -v -i 03-${dir_name}_AP_FQ_MPFR_for_Vidjil_reheader.fastq -o 03-${dir_name}_AP_FQ_MPFR_for_Vidjil_reheader.fasta
    fastqc -t $NPROC -q -o ./qc 03-${dir_name}_AP_FQ_MPF_primers-pass.fastq
    fastqc -t $NPROC -q -o ./qc 03-${dir_name}_AP_FQ_MPFR_primers-pass.fastq
    fastqc -t $NPROC -q -o ./qc 03-${dir_name}_AP_FQ_MPFR_for_Vidjil_reheader.fastq

    # Remove duplicate sequences
    CollapseSeq.py -s 03-${dir_name}_AP_FQ_MPFR_for_Vidjil_reheader.fastq -n 20 --inner --uf CPRIMER --cf VPRIMER --act set --outname 04-${dir_name}_AP_FQ_MPFR_CS | tee -a 00-log.txt
    fastqc -t $NPROC -q -o ./qc 04-${dir_name}_AP_FQ_MPFR_CS_collapse-unique.fastq

    # Extract repeating sequences
    SplitSeq.py group -s 04-${dir_name}_AP_FQ_MPFR_CS_collapse-unique.fastq -f DUPCOUNT --num 2 --outname 05-${dir_name}_AP_FQ_MPFR_CS_ER | tee -a 00-log.txt
    fastqc -t $NPROC -q -o ./qc 05-${dir_name}_AP_FQ_MPFR_CS_ER_atleast-2.fastq

    # Create an annotation table
    ParseHeaders.py table -s 05-${dir_name}_AP_FQ_MPFR_CS_ER_atleast-2.fastq -f ID DUPCOUNT CPRIMER VPRIMER | tee -a 00-log.txt

    rm ./qc/*.zip
    cd ..
    echo "***** Processing of sample $dir_name finished *****"
    echo ""

done
