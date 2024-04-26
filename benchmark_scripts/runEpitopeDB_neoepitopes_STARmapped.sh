pathDS=$(pwd)/
path1="" #specify full path of the install directory eg. /usr/bin/cfRNA-neoepitopes
pathDB=${path1}/nullomer_lists/
pathMI=${path1}/mapping_indices/
cd ${pathDS}
cores=2

sample=$1

#/data/hemberg/shared_resources/sratoolkit.2.11.2-centos_linux64/bin/fasterq-dump --split-3 ${sample}
julia ${path1}/cfRNAnullomers.v0.2.jl -f ${pathDS}${sample}_1.fastq -n ${pathDB}epitopeDB_nullomers.tsv -l 16 -S ${pathDS}${sample}_1.nullomers -p 0.0 -q 0 --logfile ${pathDS}${sample}_1.nullomers.json -N / &
julia ${path1}/cfRNAnullomers.v0.2.jl -f ${pathDS}${sample}_2.fastq -n ${pathDB}epitopeDB_nullomersRevComp.tsv -l 16 -S ${pathDS}${sample}_2.nullomers -p 0.0 -q 0 --logfile ${pathDS}${sample}_2.nullomers.json -N /
wait
awk -v pid=${sample}_1 '{if(FNR==1){f++};if(f<=2){if(FNR%4==1){a[$1]=1}}else{if(FNR%4==1){if(a[$1]==1){x=1}else{x=0}};if(x==1){print $0 > pid ".nullomers.union.epitopeDB.fastq"}}}' ${sample}_*.nullomers_epitopeDB.fastq ${sample}_1.fastq &
awk -v pid=${sample}_2 '{if(FNR==1){f++};if(f<=2){if(FNR%4==1){a[$1]=1}}else{if(FNR%4==1){if(a[$1]==1){x=1}else{x=0}};if(x==1){print $0 > pid ".nullomers.union.epitopeDB.fastq"}}}' ${sample}_*.nullomers_epitopeDB.fastq ${sample}_2.fastq
wait
#rm ${sample}_1.fastq ${sample}_2.fastq
rm ${sample}_*.nullomers_epitopeDB.fastq
rm ${sample}_1.nullomers.json ${sample}_2.nullomers.json

module load STAR samtools bcftools htslib bzip2 gatk
rm -r ${sample}_mapped
STAR --runThreadN $cores --genomeDir $path_to_indices --readFilesIn _1.nullomers.union.epitopeDB.fastq _2.nullomers.union.epitopeDB.fastq --readFilesPrefix $sample --outStd BAM_SortedByCoordinate \
--outSAMtype BAM SortedByCoordinate --outSAMmultNmax 3 \
--outFilterMultimapNmax 20 --outFilterMismatchNmax 15 --alignSJDBoverhangMin 2 --alignIntronMax 1000000 \
--outFileNamePrefix ./${sample}_mapped/ --twopassMode None  --peOverlapNbasesMin 10 > $sample.ns.bam
rm -r ${sample}_mapped 
rm ${sample}_*.nullomers.union.epitopeDB.fastq

##java gatk --java-options "-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails -Xloggc:gc_log.log -Xms10000m"
gatk MarkDuplicates --INPUT ${sample}.ns.bam --OUTPUT ${sample}.dedup.null.bam --METRICS_FILE ${sample}.dedup.metrics --REMOVE_DUPLICATES true --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT

rm test1 ${sample}.ns.bam ${sample}.dedup.metrics
mudskipper bulk --alignment ${sample}.dedup.null.bam --index /data/hemberg/shared_resources/genomes/human/GRCh38.p13.107.mskdb -l --out ${sample}.STAR2CDS.bam

#find coverage for each nullomer correspoding to specific neoepitope ; filter reads with (MAPQ < 10 AND (read_length <= (clipped_length*3) OR alignment_score <= expected_score)
#expected score calculated using linear equation (218 * aligned_length + 75)/115) allows 1 high quality mismatch (penalty=3) for mapped_length=35 and ~4 lower quality mismatches (penalty=15) for mapped_length=150
samtools view ${sample}.STAR2CDS.bam | awk -v "sid=${sample}" -F "\t" '{if(FNR==1){fn++; split(FILENAME,dbtmp,"/");split(dbtmp[length(dbtmp)],dbName,"_n")};if(fn==1){gEnsg[$1]=$2;gSname[$1]=$3}else{if(fn<=4){if($1!=""){split($1,u,";");split($3,uv,"->");if(null_list[$4]==""){ann[$4]=$5;;null_list[$4]=$1;for(i=1;i<=length(u);i++){null_neo[$4"_"u[i]]=$3;null_db[$4"_"u[i]]=dbName[1];null_wtMHC[$4"_"u[i]]="Epitope:"uv[1]","$6;null_neoMHC[$4"_"u[i]]="Epitope:"uv[2]","$7;null_pop[$4"_"u[i]]=$8;}}else{for(i=1;i<=length(u);i++){if(index(null_list[$4],u[i])==0){null_list[$4]=null_list[$4]";"u[i]};if(null_neo[$4"_"u[i]]==""){null_neo[$4"_"u[i]]=$3}else{if(index(null_neo[$4"_"u[i]],$3)==0){null_neo[$4"_"u[i]]=null_neo[$4"_"u[i]]";"$3}};if(null_db[$4"_"u[i]]==""){null_db[$4"_"u[i]]=dbName[1]}else{if(index(null_db[$4"_"u[i]],dbName[1])==0){null_db[$4"_"u[i]]=null_db[$4"_"u[i]]";"dbName[1]}};wtMHC_p="Epitope:"uv[1]",";neoMHC_p="Epitope:"uv[2]",";if(null_wtMHC[$4"_"u[i]]==""){null_wtMHC[$4"_"u[i]]=wtMHC_p""$6}else{if(index(null_wtMHC[$4"_"u[i]],wtMHC_p)==0){null_wtMHC[$4"_"u[i]]=null_wtMHC[$4"_"u[i]]","wtMHC_p""$6}};if(null_neoMHC[$4"_"u[i]]==""){null_neoMHC[$4"_"u[i]]=neoMHC_p""$7}else{if(index(null_neoMHC[$4"_"u[i]],neoMHC_p)==0){null_neoMHC[$4"_"u[i]]=null_neoMHC[$4"_"u[i]]","neoMHC_p""$7}};if(null_pop[$4"_"u[i]]==""){null_pop[$4"_"u[i]]=$8}else{if(index(null_pop[$4"_"u[i]],$8)==0){null_pop[$4"_"u[i]]=null_pop[$4"_"u[i]]""$8}};}}}}else{if($5<10){if($12!~/^AS/){next};s=0;sc=0;split($6,u,"[A-Z]",v);for(i=1;i<=length(u);i++){if(v[i]=="S"){s=s+(u[i]);sc++}};if(length($10)<(s*3)){next};trulen=length($10)-s;expS=(218*trulen+75)/115; split($12,alnS,":");if(expS>alnS[3]){next}};split($3,v,".");if(et!=v[1]){et=v[1];nt=split(null_list[et],nlu,";")};for(i=1;i<=length(nlu);i++){if(index($10,nlu[i])>0){null_count[et"_"nlu[i]]++}}}}}END{for(k in null_count){split(k,w,"_");print sid"\t"w[2]"\t"null_neo[k]"\t"null_count[k]"\t"w[1]"\t"gEnsg[w[1]]"\t"gSname[w[1]]"\t"null_db[k]"\t"ann[w[1]]"\t"null_wtMHC[k]"\t"null_neoMHC[k]"\t"null_pop[k]}}' ${pathDB}ENST_ENSG_HGNC.map ${pathDB}IEDB_neoepitopes-nullomers.filtered.tsv ${pathDB}TSNAdb_TCGA_neoepitopes-nullomers.filtered.tsv ${pathDB}TSNAdb_ICGC_neoepitopes-nullomers.filtered.tsv - > ${sample}_epitopeDB_neoepitopes_STAR_counts_raw.tsv
#Calculate coverage of the each neoepitope as the coverage of most covered corresponding nullomer; read coverage and nullomers discovered on the read are reported
sort -k1,1 -k7,7 -k3,3 -k4n -t$'\t' ${sample}_epitopeDB_neoepitopes_STAR_counts_raw.tsv | awk -F "\t" 'BEGIN{print "sample_id\tgene_id\tHGNC_symbol\ttop_nullomer\tneoepitopes\t#reads\t#nullomers\tdb_name\tannotation\twildTypeHLA\tneoEpitopeHLA\tGeneticAncestry"}{if(neo_list!=$3 || symbol!=$7){if(FNR!=1 && rcount>1){print pid"\t"et"\t"symbol"\t"null"\t"neo_list"\t"rcount"\t"ncount"\t"dbName"\t"ann"\t"wtHLA"\t"neoHLA"\t"popVar}; pid=$1; neo_list=$3; et=$6; symbol=$7; dbName=$8; ann=$9; wtHLA=$10; neoHLA=$11; popVar=$12; ncount=0}; ncount++;rcount=$4; null=$2}END{if(rcount>1){print pid"\t"et"\t"symbol"\t"null"\t"neo_list"\t"rcount"\t"ncount"\t"dbName"\t"ann"\t"wtHLA"\t"neoHLA"\t"popVar}}' > ${sample}_epitopeDB_neoepitopes_STAR_readCounts.tsv

rm ${sample}.nullomers.union.epitopeDB.bam
