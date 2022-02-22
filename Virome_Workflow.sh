#!/bin/bash

## Adapter and quality trimming, removal of phiX and human reads, merging and deduplication. ##
mkdir /path/to/logFiles;
mkdir /path/to/qc_fasta;
mkdir /path/to/qc_fastq;
gunzip -k /path/to/SampleName_S1_L001_R1_001.fastq.gz;
gunzip -k /path/to/SampleName_S1_L001_R2_001.fastq.gz;
mkdir /path/to/SampleName;
bbduk.sh in=/path/to/SampleName_S1_L001_R1_001.fastq in2=/path/to/SampleName_S1_L001_R2_001.fastq ref=/path/to/kappa_M029.fa out=/path/to/SampleName/SampleName-adaptTrimQC_R1.fastq out2=/path/to/SampleName/SampleName-adaptTrimQC_R2.fastq k=25 hdist=1 ktrim=r qtrim=rl mink=11 trimq=30 minlength=75 minavgquality=20 removeifeitherbad=f otm=t tpe=t overwrite=t 1> /path/to/SampleName/SampleName-adaptTrimQC.log.txt 2>&1;
bbduk.sh in=/path/to/SampleName/SampleName-adaptTrimQC_R1.fastq in2=/path/to/SampleName/SampleName-adaptTrimQC_R2.fastq ref=/path/to/phix174_ill.ref.fa.gz out=/path/to/SampleName/SampleName-R1-phixRemoved.fastq out2=/path/to/SampleName/SampleName-R2-phixRemoved.fastq k=31 hdist=1 overwrite=t 1> /path/to/SampleName/SampleNamephixRemoved.log.txt 2>&1;
bbmap.sh minid=.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 -Xmx64g path=/path/to/Human_GRCh38 in=/path/to/SampleName/SampleName-R1-phixRemoved.fastq in2=/path/to/SampleName/SampleName-R2-phixRemoved.fastq outu=/path/to/SampleName/SampleNamehostRemoved.fastq outm=/path/to/SampleName/SampleNamehostMatched.fastq 1>/path/to/SampleName/SampleNamehostRemoval.log.txt 2>&1;
dedupe.sh in=/path/to/SampleName/SampleNamehostRemoved.fastq out=/path/to/SampleName/SampleNamefirstDeduplication.fastq outd=/path/to/SampleName/SampleNamefirstDuplication.fastq csf=dedupe.cluster.stats overwrite=t minidentity=99 1> /path/to/SampleName/SampleNamefirstDeduplication.log.txt 2>&1;
bbmerge.sh in=/path/to/SampleName/SampleNamefirstDeduplication.fastq out=/path/to/SampleName/SampleNamefirstDeduplicationMerged.fastq outu=/path/to/SampleName/SampleNamefirstDeduplicationUnMerged.fastq 1>/path/to/SampleName/SampleNamefirstDeduplicationMerged.log.txt 2>&1;
cat /path/to/SampleName/SampleNamefirstDeduplicationMerged.fastq /path/to/SampleName/SampleNamefirstDeduplicationUnMerged.fastq > /path/to/SampleName/SampleNamefirstDeduplicationMerged_UnMerged.fastq;
dedupe.sh in=/path/to/SampleName/SampleNamefirstDeduplicationMerged_UnMerged.fastq out=/path/to/SampleName/SampleNamesecondDeduplication.fastq outd=/path/to/SampleName/SampleNamesecondtDuplication.fastq csf=dedupe.cluster.stats overwrite=t minidentity=100 ac=f 1> /path/to/SampleName/SampleNamesecondDeduplication.log.txt 2>&1;
bbduk.sh in=/path/to/SampleName/SampleNamesecondDeduplication.fastq out=/path/to/SampleName/SampleNamesecondDeduplication_filtered.fastq minlength=75 overwrite=t 1>/path/to/SampleName/SampleNamesecondDeduplication_filtered.log.txt 2>&1; 
sed -n '1~4s/^@/>/p;2~4p' /path/to/SampleName/SampleNamesecondDeduplication_filtered.fastq > /path/to/SampleName/SampleNamesecondDeduplication_filtered.fasta;
cp /path/to/SampleName/SampleName-adaptTrimQC.log.txt /path/to/logFiles;
cp /path/to/SampleName/SampleNamephixRemoved.log.txt /path/to/logFiles;
cp /path/to/SampleName/SampleNamehostRemoval.log.txt /path/to/logFiles;
cp /path/to/SampleName/SampleNamefirstDeduplication.log.txt /path/to/logFiles;
cp /path/to/SampleName/SampleNamefirstDeduplicationMerged.log.txt /path/to/logFiles;
cp /path/to/SampleName/SampleNamesecondDeduplication.log.txt /path/to/logFiles;
cp /path/to/SampleName/SampleNamesecondDeduplication_filtered.log.txt /path/to/logFiles;
cp /path/to/SampleName/SampleNamesecondDeduplication_filtered.fastq /path/to/qc_fastq;
cp /path/to/SampleName/SampleNamesecondDeduplication_filtered.fasta /path/to/qc_fasta;

## Contig building, length filtering and deduplication. ##
mkdir /path/to/qc_fasta/Contigs;
mkdir /path/to/qc_fasta/Contigs/finalContigs;
idba_ud -r SampleNamesecondDeduplication_filtered.fasta --num_threads 120 -o /path/to/qc_fasta/Contigs/SampleName_Contigs;
mv /path/to/qc_fasta/Contigs/SampleName_Contigs/contig.fa /path/to/qc_fasta/Contigs/SampleName_Contigs/SampleNamecontig.fa;
bbduk.sh in=/path/to/qc_fasta/Contigs/SampleName_Contigs/SampleNamecontig.fa out=/path/to/qc_fasta/Contigs/SampleName_Contigs/SampleNamecontig_filtered.fa minlen=500;
dedupe.sh in=/path/to/qc_fasta/Contigs/SampleName_Contigs/SampleNamecontig_filtered.fa out=/path/to/qc_fasta/Contigs/SampleName_Contigs/SampleNamecontig_filtered_deduped.fa minidentity=99;
toAmos -s /path/to/qc_fasta/Contigs/SampleName_Contigs/SampleNamecontig_filtered_deduped.fa -o /path/to/qc_fasta/Contigs/SampleName_Contigs/SampleNamecontig_amos.afg;
toAmos -s /path/to/qc_fasta/Contigs/SampleName_Contigs/SampleNamecontig_filtered_deduped.fa -o /path/to/qc_fasta/Contigs/SampleName_Contigs/SampleNamecontig_amos;
minimus2 /path/to/qc_fasta/Contigs/SampleName_Contigs/SampleNamecontig_amos -D REFCOUNT=0;
cat /path/to/qc_fasta/Contigs/SampleName_Contigs/SampleNamecontig_amos.fasta /path/to/qc_fasta/Contigs/SampleName_Contigs/SampleNamecontig_amos.singletons.seq > /path/to/qc_fasta/Contigs/SampleName_Contigs/SampleName_minimus.fasta;
awk '/^>/{print ">SampleName_Contig_1" ++i; next}{print}' < /path/to/qc_fasta/Contigs/SampleName_Contigs/SampleName_minimus.fasta > /path/to/qc_fasta/Contigs/SampleName_Contigs/SampleName_ReNamed_contig.fa;
cp /path/to/qc_fasta/Contigs/SampleName_Contigs/SampleName_ReNamed_contig.fa /path/to/qc_fasta/Contigs/finalContigs;

## Map reads back to contigs. ##
bwa index /path/to/allContigs500deduped.fa
bwa mem -M -L 97,97 /path/to/allContigs500deduped.fa /path/to/SampleNamesecondDeduplication_filtered_2018.fasta > /path/to/SampleNamereadsMapped.sam;
samtools view -h -F 0x900 /path/to/SampleNamereadsMapped.sam > /path/to/SampleNamesecondaryRemoved.sam;
samtools view -h -F 0x4 /path/to/SampleNamesecondaryRemoved.sam > /path/to/SampleNamesecondaryUnMappedRemoved.sam;
samtools view -S -b /path/to/SampleNamesecondaryUnMappedRemoved.sam > /path/to/SampleNamesecondaryUnMappedRemoved.bam;
samtools sort /path/to/SampleNamesecondaryUnMappedRemoved.bam > /path/to/SampleNamesecondaryUnMappedRemoved_sorted.bam;
samtools index /path/to/SampleNamesecondaryUnMappedRemoved_sorted.bam;
samtools idxstats /path/to/SampleNamesecondaryUnMappedRemoved_sorted.bam > /path/to/SampleNamecounts.txt;

## BLASTx parameters ##
blastx -db /path/to/RefSeqViralDB -query /path/to/allContigs.fa -evalue 1e-3 -num_threads 28 -max_target_seqs 1 -outfmt "6 qseqid sseqid evalue bitscore pident nident qcovs length mismatch qlen slen" -out /path/to/concatenatedContig.blastx.out

## VirSorter parameters ##
./wrapper_phage_contigs_sorter_iPlant.pl -f /path/to/dedupAllContigs.fa --virome --wdir /path/to/dedupAllContigs --ncpu 60 --data-dir virsorter-data/

## Megablast parameters ##
blastn -task megablast -db /path/to/ntDB -query SampleName_unmapped_reads.fasta -evalue 1e-10 -num_threads 60 -max_target_seqs 1 -max_hsps 1 -outfmt "6 qseqid sseqid staxids evalue bitscore pident nident qcovs length mismatch qlen slen" -out SampleName_megablastNT.out
