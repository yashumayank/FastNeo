dataset=${2}
path="/data/hemberg/nullomers/IEDB/"
pathDS=${path}${dataset}"/TE_neoepitopes/"
pathDB=${path}epitope_DBs/TransposElements/
cd ${pathDS}

cores=2
#alignment_fasta=$path_prefix"genome.fa"

sample=${1}"

module load julia Bowtie2

/data/hemberg/shared_resources/sratoolkit.2.11.2-centos_linux64/bin/fasterq-dump --split-3 ${sample}
julia ${path}cfDNAnullomers5.jl -f ${pathDS}${sample}_1.fastq -n ${pathDB}/TE_nullomers.tsv -l 16 -S ${pathDS}${sample}_1.nullomers.fastq -p 0.0 -q 0 --logfile ${pathDS}${sample}_1.nullomers.json -N / &
julia ${path}cfDNAnullomers5.jl -f ${pathDS}${sample}_2.fastq -n ${pathDB}/TE_nullomers_revComp.tsv -l 16 -S ${pathDS}${sample}_2.nullomers.fastq -p 0.0 -q 0 --logfile ${pathDS}${sample}_2.nullomers.json -N /
wait

#grep -h -e '^@SRR' ${sample}_1.nullomers.fastq | cut -f 1 --delim=' ' > ${sample}_1.nullomers.sorted &
#grep -h -e '^@SRR' ${sample}_2.nullomers.fastq | cut -f 1 --delim=' ' > ${sample}_2.nullomers.sorted
#wait
#awk '(FNR==1){f++;}{if(f<=2){a[$1]=1}else{if(FNR%4==1){if(a[$1]==1){x=1}else{x=0}};if(x==1)print $0 > FILENAME".nullomers.union"}}' ${sample}_1.nullomers.sorted ${sample}_2.nullomers.sorted ${sample}_1.fastq ${sample}_2.fastq

awk '(FNR==1){f++;}{if(f<=2){if(FNR%4==1){a[$1]=1}}else{if(FNR%4==1){if(a[$1]==1){x=1}else{x=0}};if(x==1)print $0 > FILENAME".nullomers.union"}}' ${sample}_1.nullomers.fastq ${sample}_2.nullomers.fastq ${sample}_1.fastq ${sample}_2.fastq

rm  ${sample}_1.nullomers.sorted ${sample}_2.nullomers.sorted ${sample}_1.fastq ${sample}_2.fastq

export LD_LIBRARY_PATH=/apps/software-compiled/tbb/2018_U5-GCCcore-7.3.0/build/linux_intel64_gcc_cc7.3.0_libc2.12_kernel2.6.32_release:$LD_LIBRARY_PATH
bowtie2 -x /data/hemberg/shared_resources/genomes/human/GRCh38.p7.gencode25 -p $cores --very-sensitive-local -1 ${sample}_1.fastq.nullomers.union -2 ${sample}_2.fastq.nullomers.union -S ${sample}.nullomers.union.sam

module purge
module load perl gcc htslib python bcftools samtools bedtools
samtools sort ${sample}.nullomers.union.sam |samtools rmdup - ${sample}.nullomers.union.bam
bedtools intersect -C -a ${pathDB}TE_antigenRegions.bed -b ${sample}.nullomers.union.bam >  ${sample}.TE_antigenReads.tab

rm ${sample}*.union ${sample}*.union.sam ${sample}_1.nullomers.fastq ${sample}_2.nullomers.fastq ${sample}_1.fastq ${sample}_2.fastq
