dataset=${2}
path="/data/hemberg/nullomers/IEDB/"
pathDS=${path}${dataset}"/epitopeDB_neoepitopes/"
pathDB=${path}epitope_DBs/
cd ${pathDS}

cores=2
genome_fasta=$path_prefix"genome.fa"

sample=${1}"
module load julia Bowtie2
/data/hemberg/shared_resources/sratoolkit.2.11.2-centos_linux64/bin/fasterq-dump --split-3 ${sample}
julia ${path}cfDNAnullomers5.jl -f ${pathDS}${sample}_1.fastq -n ${pathDB}epitopeDB_nullomers.tsv -l 16 -S ${pathDS}${sample}_1.nullomers.fastq -p 0.0 -q 0 --logfile ${pathDS}${sample}_1.nullomers.json -N / &
julia ${path}cfDNAnullomers5.jl -f ${pathDS}${sample}_2.fastq -n ${pathDB}epitopeDB_nullomers_revComp.tsv -l 16 -S ${pathDS}${sample}_2.nullomers.fastq -p 0.0 -q 0 --logfile ${pathDS}${sample}_2.nullomers.json -N /
wait

grep -h -e '^@SRR' ${sample}_1.nullomers.fastq | cut -f 1 --delim=' ' > ${sample}_1.nullomers.sorted &
grep -h -e '^@SRR' ${sample}_2.nullomers.fastq | cut -f 1 --delim=' ' > ${sample}_2.nullomers.sorted
wait
awk '(FNR==1){f++;}{if(f<=2){a[$1]=1}else{if(FNR%4==1){if(a[$1]==1){x=1}else{x=0}};if(x==1)print $0 > FILENAME".nullomers.union"}}' ${sample}_1.nullomers.sorted ${sample}_2.nullomers.sorted ${sample}_1.fastq ${sample}_2.fastq

export LD_LIBRARY_PATH=/apps/software-compiled/tbb/2018_U5-GCCcore-7.3.0/build/linux_intel64_gcc_cc7.3.0_libc2.12_kernel2.6.32_release/:$LD_LIBRARY_PATH
bowtie2 -x ${path}Homo_sapiens.GRCh38.cds.all -p $cores --very-sensitive-local --no-discordant -1 ${sample}_1.fastq.nullomers.union -2 ${sample}_2.fastq.nullomers.union -S ${sample}.nullomers.union.cds.sam

module purge
module load samtools
samtools sort ${sample}.nullomers.union.cds.sam |samtools rmdup - ${sample}.nullomers.union.cds.bam 2>test

#convert this to python
samtools view ${sample}.nullomers.union.cds.bam| awk -v "sid=${sample}" -F "\t" '{if(FNR==1){fn++; split(FILENAME,dbName,"_n")};if(fn<=3){if($1!="" || NF==5){split($1,u,";");if(null_list[$4]==""){null_list[$4]=$1;for(i=1;i<=length(u);i++){null_neo[$4"_"u[i]]=$3;null_db[$4"_"u[i]]=dbName[1]}}else{for(i=1;i<=length(u);i++){if(index(null_list[$4],u[i])==0){null_list[$4]=null_list[$4]";"u[i]};if(null_neo[$4"_"u[i]]==""){null_neo[$4"_"u[i]]=$3}else{if(index(null_neo[$4"_"u[i]],$3)==0){null_neo[$4"_"u[i]]=null_neo[$4"_"u[i]]";"$3}};if(null_db[$4"_"u[i]]==""){null_db[$4"_"u[i]]=dbName[1]}else{if(index(null_db[$4"_"u[i]],dbName[1])==0){null_db[$4"_"u[i]]=null_db[$4"_"u[i]]";"dbName[1]}}}}}}else{split($3,v,".");if(et!=v[1]){et=v[1];nt=split(null_list[et],nlu,";")};for(i=1;i<=length(nlu);i++){if(index($10,nlu[i])>0){null_count[et"_"nlu[i]]++}}}}END{for(k in null_count){split(k,w,"_");print sid"\t"w[2]"\t"null_neo[k]"\t"null_count[k]"\t"w[1]"\t"neo_db[k]}}' ../../epitope_DBs/IEDB_neoepitopes-nullomersTop20.tsv ../../epitope_DBs/TSNAdb_TCGA_neoepitopes-nullomersTop20.tsv ../../epitope_DBs/TSNAdb_ICGC_neoepitopes-nullomersTop20.tsv  - > ${sample}_nullomers_neoepitopes_counts_raw.tsv
sort -k1,1 -k5,5 -k3,3 -k4n ${sample}_nullomers_neoepitopes_counts_raw.tsv |awk '{if(neo_list!=$3 || et!=$5){if(FNR!=1 && rcount>1){print pid"\t"null"\t"neo_list"\t"rcount"\t"ncount"\t"et}; pid=$1; neo_list=$3; et=$5; ncount=0}; ncount++;rcount=$4; null=$2}END{print pid"\t"null"\t"neo_list"\t"rcount}' > ${sample}_nullomers_neoepitopes_readCounts.tsv

rm ${sample}*.union ${sample}_1.nullomers.fastq ${sample}_2.nullomers.fastq ${sample}_1.fastq ${sample}_2.fastq ${sample}_1.nullomers.sorted ${sample}_2.nullomers.sorted ${sample}.nullomers.union.cds.sam
