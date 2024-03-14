#!/bin/bash
path1='$(pwd)/'

cores=2
overhang=99 #149
path_prefix="/data/hemberg/shared_resources/genomes/human/GRCh38.p13.107."
path_to_indices=$path_prefix"STARindices"$overhang
genome_fasta=$path_prefix"genome.fa"
dbSNP_vcf=/data/hemberg/shared_resources/genomes/human/"1000GENOMES-phase_3.vcf.gz"

for x in $(seq $1 $2 $3);do n="SRR"$x
/data/hemberg/shared_resources/sratoolkit.2.11.2-centos_linux64/bin/fasterq-dump --split-3 $n
sleep 5

stSTAR=$(date +%s)
cspeed=$(lscpu|grep 'Model name'|awk -F ":" '{gsub(/^ */,"",$2);print $2}')

module load STAR samtools bcftools htslib bzip2 gatk
rm -r ${n}_mapped
STAR --runThreadN $cores --genomeDir $path_to_indices --readFilesIn _1.fastq _2.fastq --readFilesPrefix $n --outStd BAM_SortedByCoordinate \
--outSAMtype BAM SortedByCoordinate --outSAMmultNmax 3 \
--outFilterMultimapNmax 20 --outFilterMismatchNmax 15 --alignSJDBoverhangMin 2 --alignIntronMax 1000000 \
--outFileNamePrefix ./${n}_mapped/ --twopassMode None  --peOverlapNbasesMin 10 > $n.bam

rm -r ${n}_mapped ${n}_1.fastq ${n}_2.fastq
##java gatk --java-options "-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails -Xloggc:gc_log.log -Xms10000m"
gatk MarkDuplicates --INPUT ${n}.bam --OUTPUT ${n}.dedup.bam --METRICS_FILE ${n}.dedup.metrics --REMOVE_DUPLICATES true --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT
rm $n.bam
spSTAR=$(date +%s)
tidSTAR=$(( $spSTAR - $stSTAR ))

#--------- bcftools ----------
stbcf=$(date +%s)
bcftools mpileup -q 10 -Q 30 -g 0 -f $genome_fasta $n.dedup.bam | bcftools call -mv -Ov -o $n.bcftools.vcf
#DP4 altFr+altRv >= 3
awk -F "\t" '{split($8,v,";");for(i=1;i<=length(v);i++){if(v[i]~/^DP4/){split(v[i],u,"[,=]");if(u[4]+u[5]>=3){print};next}}}' ${n}.bcftools.vcf > ${n}.bcftools.vcf1
mv ${n}.bcftools.vcf1 ${n}.bcftools.vcf
spbcf=$(date +%s)

#--------- both haplotypeCaller and Lofreq ----------
sthap=$(date +%s)
stlo=$(date +%s)
gatk AddOrReplaceReadGroups -I ${n}.dedup.bam -O ${n}.dedupRG.bam -RGID 4 -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM 20
gatk SplitNCigarReads -R ${genome_fasta} -I ${n}.dedupRG.bam -O ${n}.split.bam
rm ${n}.dedupRG.bam
rm ${n}.dedup.metrics # ${n}.dedup.bam ${n}.dedup.bai

gatk BaseRecalibrator --use-original-qualities -known-sites ${dbSNP_vcf} -R ${genome_fasta} -I ${n}.split.bam -O ${n}.recal.out
gatk ApplyBQSR --add-output-sam-program-record --use-original-qualities --bqsr-recal-file ${n}.recal.out -R ${genome_fasta} -I ${n}.split.bam -O ${n}.bsqr.bam
rm ${n}.recal.out ${n}.split.bam ${n}.split.bai
splo=$(date +%s)

#--------- haplotypeCaller ----------
gatk HaplotypeCaller -dont-use-soft-clipped-bases -R ${genome_fasta} --dbsnp ${dbSNP_vcf} -I ${n}.bsqr.bam -O ${n}.hc.vcf.gz
gatk VariantFiltration --window 35 --cluster 3 --filter-name "FS" --filter "FS > 30.0" --filter-name "QD" --filter "QD < 3.0" --R ${genome_fasta} --V ${n}.hc.vcf.gz -O ${n}.hc.filtered.vcf
#AD >= 3
awk -F "\t" '($7=="PASS"){split($10,v,":");split(v[2],u,",");for(i=2;i<=length(u);i++){if(u[i]>=3){print;next}}}' ${n}.hc.filtered.vcf > ${n}.hapCall.vcf
rm  ${n}.hc.filtered.vcf ${n}.hc.vcf.gz
sphap=$(date +%s)

#--------- Lofreq ----------
stlo1=$(date +%s)
rm ${n}.lofreq.vcf
/PHShome/mm1169/usr/bin/lofreq_star-2.1.5/bin/lofreq call -f ${genome_fasta} -o ${n}.lofreq.vcf ${n}.bsqr.bam
#/PHShome/mm1169/usr/bin/lofreq_star-2.1.5/lofreq call-parallel --pp-threads 8 -f ref.fa -o vars.vcf aln.bam
rm ${n}.bsqr.bam ${n}.bsqr.bai
splo1=$(date +%s)
wait

module purge
#--------- all tools: VEP and cds2neopitopes ----------
module load perl gcc htslib python bcftools samtools
sthap1=$(date +%s)
#awk '{if($0!~/^#.*/){if($4~/[ATGC]/){if(substr($5,length($5)-3)==",<*>"){$5=substr($5,1,length($5)-4)};print "chr"$1"\t"$2"\t.\t"$4"\t"$5}}}' $1.vcf
~/ensembl-vep/vep -af_gnomade --no_stats --biotype --buffer_size 20000 --distance 5000 --polyphen b --regulatory --sift b \
--species homo_sapiens --symbol --transcript_version --plugin NMD --cache --offline --use_given_ref --merged --fork $cores \
-i ${n}.hapCall.vcf -o STDOUT --format vcf --force_overwrite --tab| awk '{\
if($1~/^#.*/){next};if($9~/-/){split($9,u,"-")}else{u[1]=$9;u[2]=$9;};\
split($5,v,".");split($4,z,".");split($25,w,"[()]");if(w[1]=="-")w[2]=0;if($27=="-")$27=0;\
if($7~/frameshift/ || $7~/missense/)print $2"\t"z[1]"\t"v[1]"\t"$7"\t"u[1]"\t"u[2]"\t"$12"\t"w[2]"\t"$26"\t"$27"\t"$45}' >  ${n}.hapCall.vep

awk '(FNR==NR){if($3~/^ENS/ && ($4~/frameshift/ || $4~/missense/)){trCnt[$3]++;conseq[$3"_"trCnt[$3]]=$4;start[$3"_"trCnt[$3]]=$5;end[$3"_"trCnt[$3]]=$6;gene[$3]=$2;gsub(/[a-z]/,"",$7);split($7,u,"/");org[$3"_"trCnt[$3]]=u[1];var[$3"_"trCnt[$3]]=u[2];next}}{if($0~/^>.*/){if(p!=""){for(i=1;i<=p;i++){neoSeq="";e1=end[pid[2]"_"i];if(conseq[pid[2]"_"i]~/frameshift_variant/){c1="f"}else{c1="m"};if(c1=="m"){c1=c1 "_" substr(cdsSeq,start[pid[2]"_"i],1)}printf(">%s:%s:%d:%d:%s",pid[2],gene[pid[2]],i,start[pid[2]"_"i],c1);if(start[pid[2]"_"i]>1){neoSeqStart=start[pid[2]"_"i]-148;if(neoSeqStart<1){neoSeqStart=1};neoSeqStart=int(neoSeqStart/3)*3+1;neoSeq = substr(cdsSeq,neoSeqStart,start[pid[2]"_"i]-neoSeqStart);}else{start[pid[2]"_"i]=1;neoSeqStart=0};printf(":%s\t%s\n",start[pid[2]"_"i]-neoSeqStart+1,cdsDesc);if(org[pid[2]"_"i]=="-" || org[pid[2]"_"i]==""){neoSeq=neoSeq substr(cdsSeq,start[pid[2]"_"i],1) var[pid[2]"_"i] substr(cdsSeq,e1,1);}else{neoSeq=neoSeq var[pid[2]"_"i]};if(e1<length(cdsSeq)){if(c1=="f"){neoSeqEnd=length(cdsSeq);}else{neoSeqEnd=e1+149;if(neoSeqEnd>length(cdsSeq)){neoSeqEnd=length(cdsSeq)}};neoSeq = neoSeq substr(cdsSeq,e1+1,neoSeqEnd-e1)};print neoSeq}};split($1,pid,"[>.]");p=trCnt[pid[2]];if(p!=""){$1="";cdsDesc=$0;cdsSeq=""};}else{if(p!=""){gsub(/\s/,"",$0);cdsSeq=cdsSeq $0}}}' ${n}.hapCall.vep ../Homo_sapiens.GRCh38.cds.all.fa > ${n}.hapCall.neoCDS.fasta

module purge
python ../cds2neopitopes.py ../epitope_DBs/IEDB_neoepitopes_per_ENST2.tab ../epitope_DBs/TSNAdb_frequent_ICGC_per_ENST3.tab ../epitope_DBs/TSNAdb_frequent_TCGA_per_ENST3.tab ${n}.hapCall.neoCDS.fasta ${n}_hapCall
sphap1=$(date +%s)

module load perl gcc htslib python bcftools samtools
stlo2=$(date +%s)
~/ensembl-vep/vep  -af_gnomade --no_stats --biotype --buffer_size 20000 --distance 5000 --polyphen b --regulatory --sift b \
--species homo_sapiens --symbol --transcript_version --plugin NMD --cache --offline --use_given_ref --merged --fork $cores \
-i ${n}.lofreq.vcf -o STDOUT --format vcf --force_overwrite --tab| awk '{\
if($1~/^#.*/){next};if($9~/-/){split($9,u,"-")}else{u[1]=$9;u[2]=$9;};\
split($5,v,".");split($4,z,".");split($25,w,"[()]");if(w[1]=="-")w[2]=0;if($27=="-")$27=0;\
if($7~/frameshift/ || $7~/missense/)print $2"\t"z[1]"\t"v[1]"\t"$7"\t"u[1]"\t"u[2]"\t"$12"\t"w[2]"\t"$26"\t"$27"\t"$45}' >  ${n}.lofreq.vep

awk '(FNR==NR){if($3~/^ENS/ && ($4~/frameshift/ || $4~/missense/)){trCnt[$3]++;conseq[$3"_"trCnt[$3]]=$4;start[$3"_"trCnt[$3]]=$5;end[$3"_"trCnt[$3]]=$6;gene[$3]=$2;gsub(/[a-z]/,"",$7);split($7,u,"/");org[$3"_"trCnt[$3]]=u[1];var[$3"_"trCnt[$3]]=u[2];next}}{if($0~/^>.*/){if(p!=""){for(i=1;i<=p;i++){neoSeq="";e1=end[pid[2]"_"i];if(conseq[pid[2]"_"i]~/frameshift_variant/){c1="f"}else{c1="m"};if(c1=="m"){c1=c1 "_" substr(cdsSeq,start[pid[2]"_"i],1)}printf(">%s:%s:%d:%d:%s",pid[2],gene[pid[2]],i,start[pid[2]"_"i],c1);if(start[pid[2]"_"i]>1){neoSeqStart=start[pid[2]"_"i]-148;if(neoSeqStart<1){neoSeqStart=1};neoSeqStart=int(neoSeqStart/3)*3+1;neoSeq = substr(cdsSeq,neoSeqStart,start[pid[2]"_"i]-neoSeqStart);}else{start[pid[2]"_"i]=1;neoSeqStart=0};printf(":%s\t%s\n",start[pid[2]"_"i]-neoSeqStart+1,cdsDesc);if(org[pid[2]"_"i]=="-" || org[pid[2]"_"i]==""){neoSeq=neoSeq substr(cdsSeq,start[pid[2]"_"i],1) var[pid[2]"_"i] substr(cdsSeq,e1,1);}else{neoSeq=neoSeq var[pid[2]"_"i]};if(e1<length(cdsSeq)){if(c1=="f"){neoSeqEnd=length(cdsSeq);}else{neoSeqEnd=e1+149;if(neoSeqEnd>length(cdsSeq)){neoSeqEnd=length(cdsSeq)}};neoSeq = neoSeq substr(cdsSeq,e1+1,neoSeqEnd-e1)};print neoSeq}};split($1,pid,"[>.]");p=trCnt[pid[2]];if(p!=""){$1="";cdsDesc=$0;cdsSeq=""};}else{if(p!=""){gsub(/\s/,"",$0);cdsSeq=cdsSeq $0}}}' ${n}.lofreq.vep ../Homo_sapiens.GRCh38.cds.all.fa > ${n}.lofreq.neoCDS.fasta

module purge
python ../cds2neopitopes.py ../epitope_DBs/IEDB_neoepitopes_per_ENST2.tab ../epitope_DBs/TSNAdb_frequent_ICGC_per_ENST3.tab ../epitope_DBs/TSNAdb_frequent_TCGA_per_ENST3.tab ${n}.lofreq.neoCDS.fasta ${n}_lofreq
splo2=$(date +%s)

module load perl gcc htslib python bcftools samtools
stbcf1=$(date +%s)
~/ensembl-vep/vep -af_gnomade --no_stats --biotype --buffer_size 20000 --distance 5000 --polyphen b --regulatory --sift b \
--species homo_sapiens --symbol --transcript_version --plugin NMD --cache --offline --use_given_ref --merged --fork $cores \
-i ${n}.bcftools.vcf -o STDOUT --format vcf --force_overwrite --tab| awk '{\
if($1~/^#.*/){next};if($9~/-/){split($9,u,"-")}else{u[1]=$9;u[2]=$9;};\
split($5,v,".");split($4,z,".");split($25,w,"[()]");if(w[1]=="-")w[2]=0;if($27=="-")$27=0;\
if($7~/frameshift/ || $7~/missense/)print $2"\t"z[1]"\t"v[1]"\t"$7"\t"u[1]"\t"u[2]"\t"$12"\t"w[2]"\t"$26"\t"$27"\t"$45}' >  ${n}.bcftools.vep

awk '(FNR==NR){if($3~/^ENS/ && ($4~/frameshift/ || $4~/missense/)){trCnt[$3]++;conseq[$3"_"trCnt[$3]]=$4;start[$3"_"trCnt[$3]]=$5;end[$3"_"trCnt[$3]]=$6;gene[$3]=$2;gsub(/[a-z]/,"",$7);split($7,u,"/");org[$3"_"trCnt[$3]]=u[1];var[$3"_"trCnt[$3]]=u[2];next}}{if($0~/^>.*/){if(p!=""){for(i=1;i<=p;i++){neoSeq="";e1=end[pid[2]"_"i];if(conseq[pid[2]"_"i]~/frameshift_variant/){c1="f"}else{c1="m"};if(c1=="m"){c1=c1 "_" substr(cdsSeq,start[pid[2]"_"i],1)}printf(">%s:%s:%d:%d:%s",pid[2],gene[pid[2]],i,start[pid[2]"_"i],c1);if(start[pid[2]"_"i]>1){neoSeqStart=start[pid[2]"_"i]-148;if(neoSeqStart<1){neoSeqStart=1};neoSeqStart=int(neoSeqStart/3)*3+1;neoSeq = substr(cdsSeq,neoSeqStart,start[pid[2]"_"i]-neoSeqStart);}else{start[pid[2]"_"i]=1;neoSeqStart=0};printf(":%s\t%s\n",start[pid[2]"_"i]-neoSeqStart+1,cdsDesc);if(org[pid[2]"_"i]=="-" || org[pid[2]"_"i]==""){neoSeq=neoSeq substr(cdsSeq,start[pid[2]"_"i],1) var[pid[2]"_"i] substr(cdsSeq,e1,1);}else{neoSeq=neoSeq var[pid[2]"_"i]};if(e1<length(cdsSeq)){if(c1=="f"){neoSeqEnd=length(cdsSeq);}else{neoSeqEnd=e1+149;if(neoSeqEnd>length(cdsSeq)){neoSeqEnd=length(cdsSeq)}};neoSeq = neoSeq substr(cdsSeq,e1+1,neoSeqEnd-e1)};print neoSeq}};split($1,pid,"[>.]");p=trCnt[pid[2]];if(p!=""){$1="";cdsDesc=$0;cdsSeq=""};}else{if(p!=""){gsub(/\s/,"",$0);cdsSeq=cdsSeq $0}}}' ${n}.bcftools.vep ../Homo_sapiens.GRCh38.cds.all.fa > ${n}.bcftools.neoCDS.fasta

module purge
python ../cds2neopitopes.py ../epitope_DBs/IEDB_neoepitopes_per_ENST2.tab ../epitope_DBs/TSNAdb_frequent_ICGC_per_ENST3.tab ../epitope_DBs/TSNAdb_frequent_TCGA_per_ENST3.tab ${n}.bcftools.neoCDS.fasta ${n}_bcftools
spbcf1=$(date +%s)

tidhap=$(( $sphap1 - $sthap1 + $sphap - $sthap + $spSTAR - $stSTAR))
tidlo=$(( $splo2 - $stlo2 + $splo1 - $stlo1 + $splo - $stlo + $spSTAR - $stSTAR))
tidbcf=$(( $spbcf1 - $stbcf1 + $spbcf - $stbcf + $spSTAR - $stSTAR))
echo "${n},${tidhap},${tidlo},${tidbcf},${cspeed}" >> runtime.log

done
