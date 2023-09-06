dataset=${2}
path="/data/hemberg/nullomers/IEDB/"
pathDS=${path}${dataset}"/fusion_neoepitopes/"
pathDB=${path}epitope_DBs/fusion/
cd ${pathDS}

cores=2
sample=${1}"

alignment_file="ChimerKB_bowtie"

module load julia Bowtie2
/data/hemberg/shared_resources/sratoolkit.2.11.2-centos_linux64/bin/fasterq-dump --split-3 ${sample}
julia ${path}cfDNAnullomers5.jl -f ${pathDS}${sample}_1.fastq -n ${pathDB}ChimerKB_nullomers.tsv -l 16 -S ${pathDS}${sample}_1.nullomers.fastq -p 0.0 -q 0 --logfile ${pathDS}${sample}_1.nullomers.json -N / &
julia ${path}cfDNAnullomers5.jl -f ${pathDS}${sample}_2.fastq -n ${pathDB}ChimerKB_nullomers_revComp.tsv -l 16 -S ${pathDS}${sample}_2.nullomers.fastq -p 0.0 -q 0 --logfile ${pathDS}${sample}_2.nullomers.json -N /
wait

#grep -h -e '^@SRR' ${sample}_1.nullomers.fastq | cut -f 1 --delim=' ' > ${sample}_1.nullomers.sorted &
#grep -h -e '^@SRR' ${sample}_2.nullomers.fastq | cut -f 1 --delim=' ' > ${sample}_2.nullomers.sorted
#wait
#awk '(FNR==1){f++;}{if(f<=2){a[$1]=1}else{if(FNR%4==1){if(a[$1]==1){x=1}else{x=0}};if(x==1)print $0 > FILENAME".nullomers.union"}}' ${sample}_1.nullomers.sorted ${sample}_2.nullomers.sorted ${sample}_1.fastq ${sample}_2.fastq

awk '(FNR==1){f++;}{if(f<=2){if(FNR%4==1){a[$1]=1}}else{if(FNR%4==1){if(a[$1]==1){x=1}else{x=0}};if(x==1)print $0 > FILENAME".nullomers.union"}}' ${sample}_1.nullomers.fastq ${sample}_2.nullomers.fastq ${sample}_1.fastq ${sample}_2.fastq

export LD_LIBRARY_PATH=/apps/software-compiled/tbb/2018_U5-GCCcore-7.3.0/build/linux_intel64_gcc_cc7.3.0_libc2.12_kernel2.6.32_release/:$LD_LIBRARY_PATH
bowtie2 -x ${pathDB}ChimerKB_bowtie -p $cores --very-sensitive --mp 3 -1 ${sample}_1.fastq.nullomers.union -2 ${sample}_2.fastq.nullomers.union -S ${sample}.nullomers.union.sam

module purge
module load htslib bcftools samtools bedtools
samtools sort ${sample}.nullomers.union.sam |samtools rmdup - ${sample}.nullomers.union.bam
samtools view ${sample}.nullomers.union.bam| awk -v "sid=${sample}" -F "\t" '{if(FNR==NR){if($1!="" || NF==7){split($1,u,";");null_list[$3]=$1;for(i=1;i<length(u);i++){null_neo[$3"_"u[i]]=$2;null_genes[$3"_"u[i]]=$5}}}else{split($3,v,".");if(et!=v[1]){et=v[1];nt=split(null_list[et],nlu,";")};for(i=1;i<length(nlu);i++){if(index($10,nlu[i])>0){null_count[et"_"nlu[i]]++}}}}END{for(k in null_count){split(k,w,"_");print sid"\t"w[2]"\t"null_neo[k]"\t"w[1]"\t"null_genes[k]"\t"null_count[k]}}' ${pathDB}ChimerKB_n_Pub_junction_nullomers.tsv - > ${sample}_nullomers_fusions_counts_raw.tsv

#samtools idxstats ${sample}.nullomers.union.bam| awk ' {print $1" "$3}' > ${sample}.fusionJunction_Reads.tab

rm ${sample}*.union ${sample}*.union.sam ${sample}_1.nullomers.fastq ${sample}_2.nullomers.fastq ${sample}_1.fastq ${sample}_2.fastq
