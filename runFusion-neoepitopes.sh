pathDS=$(pwd)/
path1="" #specify full path of the install directory eg. /usr/bin/cfRNA-neoepitopes
pathDB=${path1}/nullomer_lists/
pathMI=${path1}/mapping_indices/
cd ${pathDS}
cores=2

sample=$1
#julia
#/data/hemberg/shared_resources/sratoolkit.2.11.2-centos_linux64/bin/fasterq-dump --split-3 ${sample}
julia ${path1}/cfRNAnullomers.v0.2.jl -f ${pathDS}${sample}_1.fastq -n ${pathDB}ChimerDB_nullomers.tsv -l 16 -S ${pathDS}${sample}_1.nullomers -p 0.0 -q 0 --logfile ${pathDS}${sample}_1.nullomers.json -N / &
julia ${path1}/cfRNAnullomers.v0.2.jl -f ${pathDS}${sample}_2.fastq -n ${pathDB}ChimerDB_nullomersRevComp.tsv -l 16 -S ${pathDS}${sample}_2.nullomers -p 0.0 -q 0 --logfile ${pathDS}${sample}_2.nullomers.json -N /
wait

awk -v pid=${sample}_1 '{if(FNR==1){f++};if(f<=2){if(FNR%4==1){a[$1]=1}}else{if(FNR%4==1){if(a[$1]==1){x=1}else{x=0}};if(x==1){print $0 > pid ".nullomers.union.ChimerDB.fastq"}}}' ${sample}_*.nullomers_ChimerDB.fastq ${sample}_1.fastq &
awk -v pid=${sample}_2 '{if(FNR==1){f++};if(f<=2){if(FNR%4==1){a[$1]=1}}else{if(FNR%4==1){if(a[$1]==1){x=1}else{x=0}};if(x==1){print $0 > pid ".nullomers.union.ChimerDB.fastq"}}}' ${sample}_*.nullomers_ChimerDB.fastq ${sample}_2.fastq
wait

rm ${sample}_1.fastq ${sample}_2.fastq
rm ${sample}_*.nullomers_ChimerDB.fastq
rm ${sample}_1.nullomers.json ${sample}_2.nullomers.json

#bowtie2
bowtie2 -x ${pathMI}ChimerDB_bowtie -p $cores --no-unal --omit-sec-seq --very-sensitive-local -1 ${sample}_1.nullomers.union.ChimerDB.fastq -2 ${sample}_2.nullomers.union.ChimerDB.fastq -S ${sample}.nullomers.union.ChimerDB.sam
rm ${sample}_*.nullomers.union.ChimerDB.fastq

#samtools
samtools sort ${sample}.nullomers.union.ChimerDB.sam | samtools rmdup - ${sample}.nullomers.union.ChimerDB.bam 2 > test
rm test ${sample}.nullomers.union.ChimerDB.sam

#find coverage for each nullomer correspoding to specific fusion ; filter reads with (MAPQ < 30 AND (read_length <= (clipped_length*3) OR alignment_score <= expected_score)
#expected score calculated using linear equation (218 * aligned_length + 75) / 115 , which allows 1 mismatch (penalty=2) for mapped_length=35 and ~4 lower quality mismatches (penalty=15) for mapped_length=150
#alternatively a stringent expected score could be used (218 * aligned_length + 420) / 115 , which allows no mismatch (penalty=0) for mapped_length=35 and ~3 lower quality mismatches (penalty=12) for mapped_length=150
samtools view ${sample}.nullomers.union.ChimerDB.bam | awk -v "sid=${sample}" -F "\t" '{if(FNR==NR){if($1!=""||NF==7){split($1,u,";");null_list[$3]=$1;junc5[$3]=$4;junc3[$3]=$6;for(i=1;i<length(u);i++){null_neo[$3"_"u[i]]=$2;null_genes[$3"_"u[i]]=$5}}}else{s=0;sc=0;split($6,u,"[A-Z]",v);for(i=1;i<=length(u);i++){if(v[i]=="S"){s=s+(u[i]);sc++}};if(length($10)<(s*3)){next};trulen=length($10)-s;if($4<505-trulen || $4>496){next};if($5<30){if($12!~/^AS/){next};expS=(218*trulen+75)/115; split($12,alnS,":");if(expS>alnS[3]){next}};split($3,v,".");if(et!=v[1]){et=v[1];nt=split(null_list[et],nlu,";")};for(i=1;i<length(nlu);i++){if(index($10,nlu[i])>0){null_count[et"_"nlu[i]]++}}}}END{for(k in null_count){split(k,w,"_");print sid"\t"w[2]"\t"null_neo[k]"\t"w[1]"\t"null_genes[k]"\t"null_count[k]"\t"junc5[w[1]]"\t"junc3[w[1]]}}' ${pathDB}ChimerDB-nullomersEdge6.tsv - > ${sample}_fusion_neoepitopes_counts_raw.tsv

#Calculate coverage of the each neoepitope or gene fusion junction as the coverage of most covered corresponding nullomer; Coverage of the TE is calculated as the total reads mapped
sort -k1,1 -k4,4 -k6n -t$'\t' ${sample}_fusion_neoepitopes_counts_raw.tsv | awk -F "\t" 'BEGIN{print "sample_id\tChimerKB_id\ttop_nullomer\tneoepitopes\t#reads\t#nullomers\tannotation\tjunction5\tjunction3"}{if(et!=$4){if(FNR!=1 && rcount>1){print pid"\t"et"\t"null"\t"neo_list"\t"rcount"\t"ncount"\t"ann"\t"junc5"\t"junc3}; pid=$1; neo_list=$3; et=$4; ann=$5; junc5=$7; junc3=$8; ncount=0}; ncount++;rcount=$6; null=$2}END{if(rcount>1){print pid"\t"et"\t"null"\t"neo_list"\t"rcount"\t"ncount"\t"ann"\t"junc5"\t"junc3}}' > ${sample}_fusion_neoepitopes_readCounts.tsv
#rm ${sample}.nullomers.union.ChimerDB.bam

rm ${sample}.nullomers.union.epitopeDB.bam ${sample}.nullomers.union.TE.bam ${sample}.nullomers.union.ChimerDB.bam ${sample}.nullomers.union2.TE.bam
