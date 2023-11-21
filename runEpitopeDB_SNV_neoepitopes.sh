pathDS=$(pwd)/
path1="" #specify full path of the install directory eg. /usr/bin/cfRNA-neoepitopes
pathDB=${path1}/nullomer_lists/
pathMI=${path1}/mapping_indices/
cd ${pathDS}
cores=2

sample=$1
#julia
#/data/hemberg/shared_resources/sratoolkit.2.11.2-centos_linux64/bin/fasterq-dump --split-3 ${sample}
julia ${path1}/cfRNAnullomers.v0.2.jl -f ${pathDS}${sample}_1.fastq -n ${pathDB}epitopeDB_nullomers.tsv -l 16 -S ${pathDS}${sample}_1.nullomers -p 0.0 -q 0 --logfile ${pathDS}${sample}_1.nullomers.json -N / &
julia ${path1}/cfRNAnullomers.v0.2.jl -f ${pathDS}${sample}_2.fastq -n ${pathDB}epitopeDB_nullomersRevComp.tsv -l 16 -S ${pathDS}${sample}_2.nullomers -p 0.0 -q 0 --logfile ${pathDS}${sample}_2.nullomers.json -N /
wait
awk -v pid=${sample}_1 '{if(FNR==1){f++};if(f<=2){if(FNR%4==1){a[$1]=1}}else{if(FNR%4==1){if(a[$1]==1){x=1}else{x=0}};if(x==1){print $0 > pid ".nullomers.union.epitopeDB.fastq"}}}' ${sample}_*.nullomers_epitopeDB.fastq ${sample}_1.fastq &
awk -v pid=${sample}_2 '{if(FNR==1){f++};if(f<=2){if(FNR%4==1){a[$1]=1}}else{if(FNR%4==1){if(a[$1]==1){x=1}else{x=0}};if(x==1){print $0 > pid ".nullomers.union.epitopeDB.fastq"}}}' ${sample}_*.nullomers_epitopeDB.fastq ${sample}_2.fastq

#rm ${sample}_1.fastq ${sample}_2.fastq
rm ${sample}_*.nullomers_epitopeDB.fastq ${sample}_*.nullomers_ChimerDB.fastq ${sample}_*.nullomers_TE.fastq
rm ${sample}_1.nullomers.json ${sample}_1.nullomers.json
wait

#bowtie2
bowtie2 -x ${pathMI}Homo_sapiens.GRCh38.cds.all -p $cores --very-sensitive-local -1 ${sample}_1.nullomers.union.epitopeDB.fastq -2 ${sample}_2.nullomers.union.epitopeDB.fastq -S ${sample}.nullomers.union.epitopeDB.sam
rm ${sample}_*.nullomers.union.epitopeDB.fastq

#samtools
samtools sort ${sample}.nullomers.union.epitopeDB.sam | samtools rmdup - ${sample}.nullomers.union.epitopeDB.bam 2>test1
rm test test1 test2 ${sample}.nullomers.union.epitopeDB.sam

##find coverage for each nullomer correspoding to specific neoepitope ; filter read length >60 bp
#samtools view ${sample}.nullomers.union.epitopeDB.bam | awk -v "sid=${sample}" -F "\t" '{if(FNR==1){fn++; split(FILENAME,dbtmp,"/");split(dbtmp[length(dbtmp)],dbName,"_n")};if(fn<=3){if($1!="" || NF==5){split($1,u,";");if(null_list[$4]==""){ann[$4]=$5;null_list[$4]=$1;for(i=1;i<=length(u);i++){null_neo[$4"_"u[i]]=$3;null_db[$4"_"u[i]]=dbName[1]}}else{for(i=1;i<=length(u);i++){if(index(null_list[$4],u[i])==0){null_list[$4]=null_list[$4]";"u[i]};if(null_neo[$4"_"u[i]]==""){null_neo[$4"_"u[i]]=$3}else{if(index(null_neo[$4"_"u[i]],$3)==0){null_neo[$4"_"u[i]]=null_neo[$4"_"u[i]]";"$3}};if(null_db[$4"_"u[i]]==""){null_db[$4"_"u[i]]=dbName[1]}else{if(index(null_db[$4"_"u[i]],dbName[1])==0){null_db[$4"_"u[i]]=null_db[$4"_"u[i]]";"dbName[1]}}}}}}else{if(length($10)>60){split($3,v,".");if(et!=v[1]){et=v[1];nt=split(null_list[et],nlu,";")};for(i=1;i<=length(nlu);i++){if(index($10,nlu[i])>0){null_count[et"_"nlu[i]]++}}}}}END{for(k in null_count){split(k,w,"_");print sid"\t"w[2]"\t"null_neo[k]"\t"null_count[k]"\t"w[1]"\t"null_db[k]"\t"ann[w[1]]}}' ${pathDB}IEDB_neoepitopes-nullomers.filtered.tsv ${pathDB}TSNAdb_TCGA_neoepitopes-nullomers.filtered.tsv ${pathDB}TSNAdb_ICGC_neoepitopes-nullomers.filtered.tsv  - > ${sample}_epitopeDB_neoepitopes_counts_raw.tsv
#find coverage for each nullomer correspoding to specific neoepitope ; filter reads with (MAPQ < 10 AND (read_length <= (clipped_length*3) OR alignment_score <= expected_score)
#expected score calculated using linear equation (218 * aligned_length + 75)/115) allows 1 high quality mismatch (penalty=3) for mapped_length=35 and ~4 lower quality mismatches (penalty=15) for mapped_length=150
samtools view ${sample}.nullomers.union.epitopeDB.bam | awk -v "sid=${sample}" -F "\t" '{if(FNR==1){fn++; split(FILENAME,dbtmp,"/");split(dbtmp[length(dbtmp)],dbName,"_n")};if(fn<=3){if($1!="" || NF==5){split($1,u,";");if(null_list[$4]==""){ann[$4]=$5;null_list[$4]=$1;for(i=1;i<=length(u);i++){null_neo[$4"_"u[i]]=$3;null_db[$4"_"u[i]]=dbName[1]}}else{for(i=1;i<=length(u);i++){if(index(null_list[$4],u[i])==0){null_list[$4]=null_list[$4]";"u[i]};if(null_neo[$4"_"u[i]]==""){null_neo[$4"_"u[i]]=$3}else{if(index(null_neo[$4"_"u[i]],$3)==0){null_neo[$4"_"u[i]]=null_neo[$4"_"u[i]]";"$3}};if(null_db[$4"_"u[i]]==""){null_db[$4"_"u[i]]=dbName[1]}else{if(index(null_db[$4"_"u[i]],dbName[1])==0){null_db[$4"_"u[i]]=null_db[$4"_"u[i]]";"dbName[1]}}}}}}else{if($5<10){if($12~/^AS/){s=0;sc=0;split($6,u,"[A-Z]",v);for(i=1;i<=length(u);i++){if(v[i]=="S"){s=s+(u[i]);sc++}};if(length($10)<(s*3)){next};trulen=length($10)-s;expS=(218*trulen+75)/115; split($12,alnS,":");if(expS>alnS[3]){next}}};split($3,v,".");if(et!=v[1]){et=v[1];nt=split(null_list[et],nlu,";")};for(i=1;i<=length(nlu);i++){if(index($10,nlu[i])>0){null_count[et"_"nlu[i]]++}}}}END{for(k in null_count){split(k,w,"_");print sid"\t"w[2]"\t"null_neo[k]"\t"null_count[k]"\t"w[1]"\t"null_db[k]"\t"ann[w[1]]}}' ${pathDB}IEDB_neoepitopes-nullomers.filtered.tsv ${pathDB}TSNAdb_TCGA_neoepitopes-nullomers.filtered.tsv ${pathDB}TSNAdb_ICGC_neoepitopes-nullomers.filtered.tsv  - > ${sample}_epitopeDB_neoepitopes_counts_raw.tsv

#Calculate coverage of the each neoepitope as the coverage of most covered corresponding nullomer; read coverage and nullomers discovered on the read are reported
sort -k1,1 -k5,5 -k3,3 -k4n -t$'\t' ${sample}_epitopeDB_neoepitopes_counts_raw.tsv | awk -F "\t" 'BEGIN{print "sample_id\tENST_id\ttop_nullomer\tneoepitopes\t#reads\t#nullomers\tdb_name\tannotation"}{if(neo_list!=$3 || et!=$5){if(FNR!=1 && rcount>1){print pid"\t"et"\t"null"\t"neo_list"\t"rcount"\t"ncount"\t"dbName"\t"ann}; pid=$1; neo_list=$3; et=$5; dbName=$6; ann=$7; ncount=0}; ncount++;rcount=$4; null=$2}END{if(rcount>1){print pid"\t"et"\t"null"\t"neo_list"\t"rcount"\t"ncount"\t"dbName"\t"ann}}' > ${sample}_epitopeDB_neoepitopes_readCounts.tsv
rm ${sample}.nullomers.union.epitopeDB.bam
