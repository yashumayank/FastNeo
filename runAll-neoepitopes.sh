pathDS=$(pwd)/
path1="" #specify full path of the install directory eg. /usr/bin/cfRNA-neoepitopes
pathDB=${path1}/nullomer_lists/
pathMI=${path1}/mapping_indices/
cd ${pathDS}
cores=2

sample=$1
#julia
#/data/hemberg/shared_resources/sratoolkit.2.11.2-centos_linux64/bin/fasterq-dump --split-3 ${sample}
julia ${path1}/cfRNAnullomers.v0.2.jl -f ${pathDS}${sample}_1.fastq -n ${pathDB}epitopeDB_nullomers.tsv,${pathDB}TE_nullomers.tsv,${pathDB}ChimerDB_nullomers.tsv -l 16 -S ${pathDS}${sample}_1.nullomers -p 0.0 -q 0 --logfile ${pathDS}${sample}_1.nullomers.json -N / &
julia ${path1}/cfRNAnullomers.v0.2.jl -f ${pathDS}${sample}_2.fastq -n ${pathDB}epitopeDB_nullomersRevComp.tsv,${pathDB}TE_nullomersRevComp.tsv,${pathDB}ChimerDB_nullomersRevComp.tsv -l 16 -S ${pathDS}${sample}_2.nullomers -p 0.0 -q 0 --logfile ${pathDS}${sample}_2.nullomers.json -N /
wait
awk -v pid=${sample}_1 '{if(FNR==1){f++};if(f<=2){if(FNR%4==1){a[$1]=1}}else{if(f==3 || f==4){if(FNR%4==1){b[$1]=1}}else{if(f==5 || f==6){if(FNR%4==1){c[$1]=1}}else{if(FNR%4==1){if(a[$1]==1){x=1}else{x=0};if(b[$1]==1){y=1}else{y=0};if(c[$1]==1){z=1}else{z=0}};if(x==1){print $0 > pid ".nullomers.union.epitopeDB.fastq"};if(y==1){print $0 > pid ".nullomers.union.ChimerDB.fastq"};if(z==1){print $0 > pid ".nullomers.union.TE.fastq"}}}}}' ${sample}_*.nullomers_epitopeDB.fastq ${sample}_*.nullomers_ChimerDB.fastq ${sample}_*.nullomers_TE.fastq ${sample}_1.fastq &
awk -v pid=${sample}_2 '{if(FNR==1){f++};if(f<=2){if(FNR%4==1){a[$1]=1}}else{if(f==3 || f==4){if(FNR%4==1){b[$1]=1}}else{if(f==5 || f==6){if(FNR%4==1){c[$1]=1}}else{if(FNR%4==1){if(a[$1]==1){x=1}else{x=0};if(b[$1]==1){y=1}else{y=0};if(c[$1]==1){z=1}else{z=0}};if(x==1){print $0 > pid ".nullomers.union.epitopeDB.fastq"};if(y==1){print $0 > pid ".nullomers.union.ChimerDB.fastq"};if(z==1){print $0 > pid ".nullomers.union.TE.fastq"}}}}}' ${sample}_*.nullomers_epitopeDB.fastq ${sample}_*.nullomers_ChimerDB.fastq ${sample}_*.nullomers_TE.fastq ${sample}_2.fastq
wait
rm ${sample}_1.fastq ${sample}_2.fastq
rm ${sample}_*.nullomers_epitopeDB.fastq ${sample}_*.nullomers_ChimerDB.fastq ${sample}_*.nullomers_TE.fastq
rm ${sample}_1.nullomers.json ${sample}_1.nullomers.json

#bowtie2
bowtie2 -x ${pathMI}Homo_sapiens.GRCh38.cds.all -p $cores --very-sensitive-local -1 ${sample}_1.nullomers.union.epitopeDB.fastq -2 ${sample}_2.nullomers.union.epitopeDB.fastq -S ${sample}.nullomers.union.epitopeDB.sam
rm ${sample}_*.nullomers.union.epitopeDB.fastq
bowtie2 -x ${pathMI}ChimerDB_bowtie -p $cores --very-sensitive --mp 3 -1 ${sample}_1.nullomers.union.ChimerDB.fastq -2 ${sample}_2.nullomers.union.ChimerDB.fastq -S ${sample}.nullomers.union.ChimerDB.sam
rm ${sample}_*.nullomers.union.ChimerDB.fastq
bowtie2 -x ${pathMI}GRCh38.p7.gencode25 -p $cores -1 ${sample}_1.nullomers.union.TE.fastq -2 ${sample}_2.nullomers.union.TE.fastq -S ${sample}.nullomers.union.TE.sam
rm ${sample}_*.nullomers.union.TE.fastq

#samtools
samtools sort ${sample}.nullomers.union.epitopeDB.sam | samtools rmdup - ${sample}.nullomers.union.epitopeDB.bam 2>test1 &
samtools sort ${sample}.nullomers.union.TE.sam | samtools rmdup - ${sample}.nullomers.union.TE.bam 2>test2
samtools sort ${sample}.nullomers.union.ChimerDB.sam | samtools rmdup - ${sample}.nullomers.union.ChimerDB.bam 2>test
wait
rm test test1 test2 ${sample}.nullomers.union.epitopeDB.sam ${sample}.nullomers.union.TE.sam ${sample}.nullomers.union.ChimerDB.sam

#find coverage for each nullomer correspoding to specific neoepitope or gene fusion junction; apply read length filter >60 bp
samtools view ${sample}.nullomers.union.epitopeDB.bam | awk -v "sid=${sample}" -F "\t" '{if(FNR==1){fn++; split(FILENAME,dbtmp,"/");split(dbtmp[length(dbtmp)],dbName,"_n")};if(fn<=3){if($1!="" || NF==5){split($1,u,";");if(null_list[$4]==""){ann[$4]=$5;null_list[$4]=$1;for(i=1;i<=length(u);i++){null_neo[$4"_"u[i]]=$3;null_db[$4"_"u[i]]=dbName[1]}}else{for(i=1;i<=length(u);i++){if(index(null_list[$4],u[i])==0){null_list[$4]=null_list[$4]";"u[i]};if(null_neo[$4"_"u[i]]==""){null_neo[$4"_"u[i]]=$3}else{if(index(null_neo[$4"_"u[i]],$3)==0){null_neo[$4"_"u[i]]=null_neo[$4"_"u[i]]";"$3}};if(null_db[$4"_"u[i]]==""){null_db[$4"_"u[i]]=dbName[1]}else{if(index(null_db[$4"_"u[i]],dbName[1])==0){null_db[$4"_"u[i]]=null_db[$4"_"u[i]]";"dbName[1]}}}}}}else{if(length($10)>60){split($3,v,".");if(et!=v[1]){et=v[1];nt=split(null_list[et],nlu,";")};for(i=1;i<=length(nlu);i++){if(index($10,nlu[i])>0){null_count[et"_"nlu[i]]++}}}}}END{for(k in null_count){split(k,w,"_");print sid"\t"w[2]"\t"null_neo[k]"\t"null_count[k]"\t"w[1]"\t"null_db[k]"\t"ann[w[1]]}}' ${pathDB}IEDB_neoepitopes-nullomers.filtered.tsv ${pathDB}TSNAdb_TCGA_neoepitopes-nullomers.filtered.tsv ${pathDB}TSNAdb_ICGC_neoepitopes-nullomers.filtered.tsv  - > ${sample}_epitopeDB_neoepitopes_counts_raw.tsv &
samtools view ${sample}.nullomers.union.ChimerDB.bam | awk -v "sid=${sample}" -F "\t" '{if(FNR==NR){if($1!=""||NF==7){split($1,u,";");null_list[$3]=$1;junc5[$3]=$4;junc3[$3]=$6;for(i=1;i<length(u);i++){null_neo[$3"_"u[i]]=$2;null_genes[$3"_"u[i]]=$5}}}else{if(length($10)>60){split($3,v,".");if(et!=v[1]){et=v[1];nt=split(null_list[et],nlu,";")};for(i=1;i<length(nlu);i++){if(index($10,nlu[i])>0){null_count[et"_"nlu[i]]++}}}}}END{for(k in null_count){split(k,w,"_");print sid"\t"w[2]"\t"null_neo[k]"\t"w[1]"\t"null_genes[k]"\t"null_count[k]"\t"junc5[w[1]]"\t"junc3[w[1]]}}' ${pathDB}ChimerDB-nullomersEdge6.tsv - > ${sample}_fusion_neoepitopes_counts_raw.tsv
samtools view -h ${sample}.nullomers.union.TE.bam | awk 'length($10)>60 || $1~/^@/' | samtools view -bS - > ${sample}.nullomers.union2.TE.bam
wait

#bedtools
#Calculate coverage of the each neoepitope or gene fusion junction as the coverage of most covered corresponding nullomer; Coverage of the TE is calculated as the total reads mapped
sort -k1,1 -k5,5 -k3,3 -k4n -t$'\t' ${sample}_epitopeDB_neoepitopes_counts_raw.tsv | awk -F "\t" 'BEGIN{print "sample_id\tENST_id\ttop_nullomer\tneoepitopes\t#reads\t#nullomers\tdb_name\tannotation"}{if(neo_list!=$3 || et!=$5){if(FNR!=1 && rcount>1){print pid"\t"et"\t"null"\t"neo_list"\t"rcount"\t"ncount"\t"dbName"\t"ann}; pid=$1; neo_list=$3; et=$5; dbName=$6; ann=$7; ncount=0}; ncount++;rcount=$4; null=$2}END{if(rcount>1){print pid"\t"et"\t"null"\t"neo_list"\t"rcount"\t"ncount"\t"dbName"\t"ann}}' > ${sample}_epitopeDB_neoepitopes_readCounts.tsv &
sort -k1,1 -k4,4 -k6n -t$'\t' ${sample}_fusion_neoepitopes_counts_raw.tsv | awk -F "\t" 'BEGIN{print "sample_id\tChimerKB_id\ttop_nullomer\tneoepitopes\t#reads\t#nullomers\tannotation\tjunction5\tjunction3"}{if(et!=$4){if(FNR!=1 && rcount>1){print pid"\t"et"\t"null"\t"neo_list"\t"rcount"\t"ncount"\t"ann"\t"junc5"\t"junc3}; pid=$1; neo_list=$3; et=$4; ann=$5; junc5=$7; junc3=$8; ncount=0}; ncount++;rcount=$6; null=$2}END{if(rcount>1){print pid"\t"et"\t"null"\t"neo_list"\t"rcount"\t"ncount"\t"ann"\t"junc5"\t"junc3}}' > ${sample}_fusion_neoepitopes_readCounts.tsv
bedtools intersect -C -a ${pathDB}TransposElements/TE_cancerRegions.bed -b ${sample}.nullomers.union2.TE.bam | awk -v "sid=${sample}" 'BEGIN{print "sample_id\tTE_ID\tposition\t#reads"}{if($7>1){print sid"\t"$4"\t"$1":"$2":"$3":"$6"\t"$7}}' >  ${sample}_TE_neoepitopes_readCounts.tsv
wait

rm ${sample}.nullomers.union.epitopeDB.bam ${sample}.nullomers.union.TE.bam ${sample}.nullomers.union.ChimerDB.bam ${sample}.nullomers.union2.TE.bam
