#!/bin/bash

path1='/data/hemberg/nullomers/IEDB/HCC2/'
cd $path1
cores=2
overhang=149  #99
path_prefix="/data/hemberg/shared_resources/genomes/human/GRCh38.p13.107."
path_to_indices=$path_prefix"STARindices"$overhang
genome_fasta=$path_prefix"genome.fa"
dbSNP_vcf=/data/hemberg/shared_resources/genomes/human/"1000GENOMES-phase_3.vcf.gz"

for x in $(seq $1 $2 $3);do n="SRR"$x
/data/hemberg/shared_resources/sratoolkit.2.11.2-centos_linux64/bin/fasterq-dump --split-3 $1
sleep 5
module load STAR samtools bcftools htslib bzip2 gatk
STAR --runThreadN $cores --genomeDir $path_to_indices --readFilesIn _1.fastq _2.fastq --readFilesPrefix $n --outStd BAM_SortedByCoordinate \
--outSAMtype BAM SortedByCoordinate --outSAMmultNmax 3 \
--outFilterMultimapNmax 20 --outFilterMismatchNmax 15 --alignSJDBoverhangMin 2 --alignIntronMax 1000000 \
--outFileNamePrefix ./${n}_mapped/ --twopassMode None  --peOverlapNbasesMin 10 > $n.bam

rm ${n}_1.fastq ${n}_2.fastq

##java gatk --java-options "-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails -Xloggc:gc_log.log -Xms10000m"
gatk MarkDuplicates --INPUT ${n}.bam --OUTPUT ${n}.dedup.bam --METRICS_FILE ${n}.dedup.metrics --REMOVE_DUPLICATES true --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT
rm $n.bam

#--------- bcftools ----------
bcftools mpileup -q 20 -Q 30 -g 0 -f $genome_fasta $n.dedup.bam | bcftools call -mv -Ov -o $n.bcftools.vcf &

#--------- both haplotypeCaller and Lofreq ----------
gatk AddOrReplaceReadGroups -I ${n}.dedup.bam -O ${n}.dedupRG.bam -RGID 4 -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM 20
gatk SplitNCigarReads -R ${genome_fasta} -I ${n}.dedupRG.bam -O ${n}.split.bam
rm ${n}.dedupRG.bam 
#rm ${n}.dedup.metrics ${n}.dedup.bam ${n}.dedup.bai

gatk BaseRecalibrator --use-original-qualities -known-sites ${dbSNP_vcf} -R ${genome_fasta} -I ${n}.split.bam -O ${n}.recal.out
gatk ApplyBQSR --add-output-sam-program-record --use-original-qualities --bqsr-recal-file ${n}.recal.out -R ${genome_fasta} -I ${n}.split.bam -O ${n}.bsqr.bam 
rm ${n}.recal.out ${n}.split.bam ${n}.split.bai

#--------- haplotypeCaller ----------

gatk HaplotypeCaller -dont-use-soft-clipped-bases -R ${genome_fasta} --dbsnp ${dbSNP_vcf} -I ${n}.bsqr.bam -O ${n}.hc.vcf.gz
gatk VariantFiltration --window 35 --cluster 3 --filter-name "FS" --filter "FS > 30.0" --filter-name "QD" --filter "QD < 2.0" --R ${genome_fasta} --V ${n}.hc.vcf.gz -O ${n}.hc.filtered.vcf
awk '$7=="PASS"' ${n}.hc.filtered.vcf > ${n}.hapCall.vcf
rm  ${n}.hc.filtered.vcf ${n}.hc.vcf.gz

#--------- Lofreq ----------
/PHShome/mm1169/usr/bin/lofreq_star-2.1.5/bin/lofreq call -f ${genome_fasta} -o ${n}.lofreq.vcf ${n}.bsqr.bam 
#/PHShome/mm1169/usr/bin/lofreq_star-2.1.5/lofreq call-parallel --pp-threads 8 -f ref.fa -o vars.vcf aln.bam
rm ${n}.bsqr.bam ${n}.bsqr.bai

wait

module purge
module load perl gcc htslib python bcftools samtools
#awk '{if($0!~/^#.*/){if($4~/[ATGC]/){if(substr($5,length($5)-3)==",<*>"){$5=substr($5,1,length($5)-4)};print "chr"$1"\t"$2"\t.\t"$4"\t"$5}}}' $1.vcf 
~/ensembl-vep/vep --filter_common -af_gnomade --no_stats --biotype --buffer_size 20000 --distance 5000 --polyphen b --regulatory --sift b \
--species homo_sapiens --symbol --transcript_version --plugin NMD --cache --offline --use_given_ref --merged --fork $cores \
-i ${n}.hapCall.vcf -o STDOUT --format vcf --force_overwrite --tab| awk '{\
if($1~/^#.*/){next};if($9~/-/){split($9,u,"-")}else{u[1]=$9;u[2]=$9;};\
split($5,v,".");split($4,z,".");split($25,w,"[()]");if(w[1]=="-")w[2]=0;if($27=="-")$27=0;\
if(w[2]<0.06 && $27<1e-5)print $2"\t"z[1]"\t"v[1]"\t"$7"\t"u[1]"\t"u[2]"\t"$12"\t"w[2]"\t"$26"\t"$27"\t"$45}' >  ${n}.hapCall.vep

~/ensembl-vep/vep --filter_common -af_gnomade --no_stats --biotype --buffer_size 20000 --distance 5000 --polyphen b --regulatory --sift b \
--species homo_sapiens --symbol --transcript_version --plugin NMD --cache --offline --use_given_ref --merged --fork $cores \
-i ${n}.lofreq.vcf -o STDOUT --format vcf --force_overwrite --tab| awk '{\
if($1~/^#.*/){next};if($9~/-/){split($9,u,"-")}else{u[1]=$9;u[2]=$9;};\
split($5,v,".");split($4,z,".");split($25,w,"[()]");if(w[1]=="-")w[2]=0;if($27=="-")$27=0;\
if(w[2]<0.06 && $27<1e-5)print $2"\t"z[1]"\t"v[1]"\t"$7"\t"u[1]"\t"u[2]"\t"$12"\t"w[2]"\t"$26"\t"$27"\t"$45}' >  ${n}.lofreq.vep

~/ensembl-vep/vep --filter_common -af_gnomade --no_stats --biotype --buffer_size 20000 --distance 5000 --polyphen b --regulatory --sift b \
--species homo_sapiens --symbol --transcript_version --plugin NMD --cache --offline --use_given_ref --merged --fork $cores \
-i ${n}.bcftools.vcf -o STDOUT --format vcf --force_overwrite --tab| awk '{\
if($1~/^#.*/){next};if($9~/-/){split($9,u,"-")}else{u[1]=$9;u[2]=$9;};\
split($5,v,".");split($4,z,".");split($25,w,"[()]");if(w[1]=="-")w[2]=0;if($27=="-")$27=0;\
if(w[2]<0.06 && $27<1e-5)print $2"\t"z[1]"\t"v[1]"\t"$7"\t"u[1]"\t"u[2]"\t"$12"\t"w[2]"\t"$26"\t"$27"\t"$45}' >  ${n}.bcftools.vep

#create a fasta with upto 50 codons from both upstream and downstream regions of each frameshift or missense variation
awk '(FNR==NR){if($3~/^ENS/ && ($4~/frameshift/ || $4~/missense/)){trCnt[$3]++;conseq[$3"_"trCnt[$3]]=$4;start[$3"_"trCnt[$3]]=$5;end[$3"_"trCnt[$3]]=$6;gene[$3]=$2;gsub(/[a-z]/,"",$7);split($7,u,"/");org[$3"_"trCnt[$3]]=u[1];var[$3"_"trCnt[$3]]=u[2];next}}{if($0~/^>.*/){if(p!=""){for(i=1;i<=p;i++){neoSeq="";e1=end[pid[2]"_"i];if(conseq[pid[2]"_"i]~/frameshift_variant/){c1="f"}else{c1="m"};if(c1=="m"){c1=c1 "_" substr(cdsSeq,start[pid[2]"_"i],1)}printf(">%s:%s:%d:%d:%s",pid[2],gene[pid[2]],i,start[pid[2]"_"i],c1);if(start[pid[2]"_"i]>1){neoSeqStart=start[pid[2]"_"i]-148;if(neoSeqStart<1){neoSeqStart=1};neoSeqStart=int(neoSeqStart/3)*3+1;neoSeq = substr(cdsSeq,neoSeqStart,start[pid[2]"_"i]-neoSeqStart);}else{start[pid[2]"_"i]=1;neoSeqStart=0};printf(":%s\t%s\n",start[pid[2]"_"i]-neoSeqStart+1,cdsDesc);if(org[pid[2]"_"i]=="-" || org[pid[2]"_"i]==""){neoSeq=neoSeq substr(cdsSeq,start[pid[2]"_"i],1) var[pid[2]"_"i] substr(cdsSeq,e1,1);}else{neoSeq=neoSeq var[pid[2]"_"i]};if(e1<length(cdsSeq)){if(c1=="f"){neoSeqEnd=length(cdsSeq);}else{neoSeqEnd=e1+149;if(neoSeqEnd>length(cdsSeq)){neoSeqEnd=length(cdsSeq)}};neoSeq = neoSeq substr(cdsSeq,e1+1,neoSeqEnd-e1)};print neoSeq}};split($1,pid,"[>.]");p=trCnt[pid[2]];if(p!=""){$1="";cdsDesc=$0;cdsSeq=""};}else{if(p!=""){gsub(/\s/,"",$0);cdsSeq=cdsSeq $0}}}' ${n}.hapCall.vep ../Homo_sapiens.GRCh38.cds.all.fa > ${n}.hapCall.neoCDS.fasta
awk '(FNR==NR){if($3~/^ENS/ && ($4~/frameshift/ || $4~/missense/)){trCnt[$3]++;conseq[$3"_"trCnt[$3]]=$4;start[$3"_"trCnt[$3]]=$5;end[$3"_"trCnt[$3]]=$6;gene[$3]=$2;gsub(/[a-z]/,"",$7);split($7,u,"/");org[$3"_"trCnt[$3]]=u[1];var[$3"_"trCnt[$3]]=u[2];next}}{if($0~/^>.*/){if(p!=""){for(i=1;i<=p;i++){neoSeq="";e1=end[pid[2]"_"i];if(conseq[pid[2]"_"i]~/frameshift_variant/){c1="f"}else{c1="m"};if(c1=="m"){c1=c1 "_" substr(cdsSeq,start[pid[2]"_"i],1)}printf(">%s:%s:%d:%d:%s",pid[2],gene[pid[2]],i,start[pid[2]"_"i],c1);if(start[pid[2]"_"i]>1){neoSeqStart=start[pid[2]"_"i]-148;if(neoSeqStart<1){neoSeqStart=1};neoSeqStart=int(neoSeqStart/3)*3+1;neoSeq = substr(cdsSeq,neoSeqStart,start[pid[2]"_"i]-neoSeqStart);}else{start[pid[2]"_"i]=1;neoSeqStart=0};printf(":%s\t%s\n",start[pid[2]"_"i]-neoSeqStart+1,cdsDesc);if(org[pid[2]"_"i]=="-" || org[pid[2]"_"i]==""){neoSeq=neoSeq substr(cdsSeq,start[pid[2]"_"i],1) var[pid[2]"_"i] substr(cdsSeq,e1,1);}else{neoSeq=neoSeq var[pid[2]"_"i]};if(e1<length(cdsSeq)){if(c1=="f"){neoSeqEnd=length(cdsSeq);}else{neoSeqEnd=e1+149;if(neoSeqEnd>length(cdsSeq)){neoSeqEnd=length(cdsSeq)}};neoSeq = neoSeq substr(cdsSeq,e1+1,neoSeqEnd-e1)};print neoSeq}};split($1,pid,"[>.]");p=trCnt[pid[2]];if(p!=""){$1="";cdsDesc=$0;cdsSeq=""};}else{if(p!=""){gsub(/\s/,"",$0);cdsSeq=cdsSeq $0}}}' ${n}.bcftools.vep ../Homo_sapiens.GRCh38.cds.all.fa > ${n}.bcftools.neoCDS.fasta
awk '(FNR==NR){if($3~/^ENS/ && ($4~/frameshift/ || $4~/missense/)){trCnt[$3]++;conseq[$3"_"trCnt[$3]]=$4;start[$3"_"trCnt[$3]]=$5;end[$3"_"trCnt[$3]]=$6;gene[$3]=$2;gsub(/[a-z]/,"",$7);split($7,u,"/");org[$3"_"trCnt[$3]]=u[1];var[$3"_"trCnt[$3]]=u[2];next}}{if($0~/^>.*/){if(p!=""){for(i=1;i<=p;i++){neoSeq="";e1=end[pid[2]"_"i];if(conseq[pid[2]"_"i]~/frameshift_variant/){c1="f"}else{c1="m"};if(c1=="m"){c1=c1 "_" substr(cdsSeq,start[pid[2]"_"i],1)}printf(">%s:%s:%d:%d:%s",pid[2],gene[pid[2]],i,start[pid[2]"_"i],c1);if(start[pid[2]"_"i]>1){neoSeqStart=start[pid[2]"_"i]-148;if(neoSeqStart<1){neoSeqStart=1};neoSeqStart=int(neoSeqStart/3)*3+1;neoSeq = substr(cdsSeq,neoSeqStart,start[pid[2]"_"i]-neoSeqStart);}else{start[pid[2]"_"i]=1;neoSeqStart=0};printf(":%s\t%s\n",start[pid[2]"_"i]-neoSeqStart+1,cdsDesc);if(org[pid[2]"_"i]=="-" || org[pid[2]"_"i]==""){neoSeq=neoSeq substr(cdsSeq,start[pid[2]"_"i],1) var[pid[2]"_"i] substr(cdsSeq,e1,1);}else{neoSeq=neoSeq var[pid[2]"_"i]};if(e1<length(cdsSeq)){if(c1=="f"){neoSeqEnd=length(cdsSeq);}else{neoSeqEnd=e1+149;if(neoSeqEnd>length(cdsSeq)){neoSeqEnd=length(cdsSeq)}};neoSeq = neoSeq substr(cdsSeq,e1+1,neoSeqEnd-e1)};print neoSeq}};split($1,pid,"[>.]");p=trCnt[pid[2]];if(p!=""){$1="";cdsDesc=$0;cdsSeq=""};}else{if(p!=""){gsub(/\s/,"",$0);cdsSeq=cdsSeq $0}}}' ${n}.lofreq.vep ../Homo_sapiens.GRCh38.cds.all.fa > ${n}.lofreq.neoCDS.fasta

#search neoepitopes from various epitopeDBs in each frameshift or missense variation
python ../cds2neopitopes.py ../epitope_DBs/IEDB_neoepitopes_per_ENST2.tab ../epitope_DBs/TSNAdb_frequent_ICGC_per_ENST3.tab ../epitope_DBs/TSNAdb_frequent_TCGA_per_ENST3.tab ${n}.hapCall.neoCDS.fasta ${n}_hapCall
python ../cds2neopitopes.py ../epitope_DBs/IEDB_neoepitopes_per_ENST2.tab ../epitope_DBs/TSNAdb_frequent_ICGC_per_ENST3.tab ../epitope_DBs/TSNAdb_frequent_TCGA_per_ENST3.tab ${n}.lofreq.neoCDS.fasta ${n}_lofreq
python ../cds2neopitopes.py ../epitope_DBs/IEDB_neoepitopes_per_ENST2.tab ../epitope_DBs/TSNAdb_frequent_ICGC_per_ENST3.tab ../epitope_DBs/TSNAdb_frequent_TCGA_per_ENST3.tab ${n}.bcftools.neoCDS.fasta ${n}_bcftools

done
