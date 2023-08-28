dataset=${2}
path="/data/hemberg/nullomers/IEDB/"
pathDS=${path}${dataset}"/fusion_neoepitopes/"
pathDB=${path}epitope_DBs/fusion/
cd ${pathDS}

cores=2
genome_fasta=$path_prefix"genome.fa"

sample=${1}"

alignment_file=ChimerKB_bowtie

module load julia Bowtie2
/data/hemberg/shared_resources/sratoolkit.2.11.2-centos_linux64/bin/fasterq-dump --split-3 ${sample}
julia ${path}cfDNAnullomers5.jl -f ${pathDS}${sample}_1.fastq -n ${pathDB}ChimerKB_nullomers.tsv -l 16 -S ${pathDS}${sample}_1.neomers.fastq -p 0.0 -q 0 --logfile ${pathDS}${sample}_1.neomers.json -N / &
julia ${path}cfDNAnullomers5.jl -f ${pathDS}${sample}_2.fastq -n ${pathDB}ChimerKB_nullomers_revComp.tsv -l 16 -S ${pathDS}${sample}_2.neomers.fastq -p 0.0 -q 0 --logfile ${pathDS}${sample}_2.neomers.json -N /
wait

grep -h -e '^@SRR' ${sample}_1.neomers.fastq | cut -f 1 --delim=' ' > ${sample}_1.neomers.sorted &
grep -h -e '^@SRR' ${sample}_2.neomers.fastq | cut -f 1 --delim=' ' > ${sample}_2.neomers.sorted
wait
awk '(FNR==1){f++;}{if(f<=2){a[$1]=1}else{if(FNR%4==1){if(a[$1]==1){x=1}else{x=0}};if(x==1)print $0 > FILENAME".neomers.union"}}' ${sample}_1.neomers.sorted ${sample}_2.neomers.sorted ${sample}_1.fastq ${sample}_2.fastq

export LD_LIBRARY_PATH=/apps/software-compiled/tbb/2018_U5-GCCcore-7.3.0/build/linux_intel64_gcc_cc7.3.0_libc2.12_kernel2.6.32_release/:$LD_LIBRARY_PATH
bowtie2 -x ${pathDB}ChimerKB_bowtie -p $cores --very-sensitive --mp 3 -1 ${sample}_1.fastq.neomers.union -2 ${sample}_2.fastq.neomers.union -S ${sample}.neomers.union.sam

module purge
module load perl gcc htslib python bcftools samtools bedtools
samtools sort ${sample}.neomers.union.sam |samtools rmdup - ${sample}.neomers.union.bam
samtools idxstats ${sample}.neomers.union.bam | awk ' {print $1" "$3}' > ${sample}.fusionJunction_Reads.tab

rm ${sample}*.union ${sample}*.union.sam ${sample}_1.neomers.fastq ${sample}_2.neomers.fastq ${sample}_1.fastq ${sample}_2.fastq
