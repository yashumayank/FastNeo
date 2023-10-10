#!/bin/bash
#BSUB -J epiDB-Null
#BSUB -o test-%J.out
#BSUB -e test-%J.err
#BSUB -q bigmem
#BSUB -n 1

#eval "$(/PHShome/mm1169/miniconda3/bin/conda shell.bash hook)"
#conda activate scanpy
cd /data/hemberg/nullomers/IEDB/epitope_DBs/fusion
#--extract Human cds 16mers
if [ ! -e "Homo_sapiens_cds_16mers.tab" ] || [ ! -s "Homo_sapiens_cds_16mers.tab" ] ; then
  awk '{if($0~/^>.*/){ln=length(x);for(i=1;i+15<=ln;i++){a[substr(x,i,16)]+=1};x=""}else{x=x $1}}END{for(j in a){print j "\t" a[j]}}' ../Homo_sapiens.GRCh38.cds.all.fa > Homo_sapiens_cds_16mers.tab
fi
#1342 have valid junctions out of 2243 fusions that are part of both ChimerKB & Pub
#1096 unique junctions and 182 with predicted neoepitopes
#9(17)-15(577) cds nullomers per junction
#0(16)-15(157) genome nullomers per junction 2(28), 3(35)

#create protein kmer to get the null peptides (not used)
#python3 ../createProteinKmers.py  ../../Homo_sapiens.GRCh38.pep.all.fa.gz > ../Homo_sapiens_protein_9mers.tab
#create genome nullomers from hg19 genome fasta (not used)
#zcat /data/hemberg/shared_resources/genomes/human/hg19.fa.gz|awk '{if($0~/^>.*/){ln=length(x);for(i=1;i+15<=ln;i++){a[substr(x,i,16)]=1};x=""}else{x=x $1}}END{for(j in a){print j}}' - > ../Homo_sapiens_genome_16mers.tab

python3 extractChimerDbNullomers.py ChimerKB_n_Pub.tab /data/hemberg/shared_resources/genomes/human/hg19.fa.gz /data/hemberg/shared_resources/genomes/human/hg19.ensGene.gtf ChimerKB_n_Pub
#Fusion junctions ouside CDS: 680
python3 ChimerKB_5_3_seqs2nullomers.py /data/hemberg/nullomers/IEDB/epitope_DBs/Homo_sapiens_cds_16mers.tab ../../Homo_sapiens.GRCh38.pep.all.fa.gz ChimerKB_n_Pub

awk '{split($1,u,";");for(i=1;i<length(u);i++)print toupper(u[i])}' ChimerKB_n_Pub_junction_nullomers.tsv > ChimerDB_nullomers.all.tsv

#Calculate nullomer frequencies in the genome
awk '(NR==FNR){a[$1]=1;next}{if($1~/^>.*/){if(seqid!=""){for(i=1;i<=(length(seq)-15);i++){xtr=substr(seq,i,16);if(xtr in a){a[xtr]++}}};seqid=$1;seq=""}else{seq=seq""$1}}END{for(j=1;j<=(length(seq)-15);j++){xtr=substr(seq,j,16);if(xtr in a){a[xtr]++}};for(k in a){print k"\t"a[k]-1}}' ChimerDB_nullomers.all.tsv /data/hemberg/shared_resources/genomes/human/GRCh38.p13.107.genome.fa > ChimerDB_nullomers_freq.tsv
#---take 20 least frequent per neoepitope; filter out nullomers with genome-frequency > 100
#awk -F "\t" 'BEGIN{ PROCINFO["sorted_in"] = "@ind_num_asc"}{if(NR==FNR){f[$1]=$2}else{split($1,u,";");delete(v);delete(vX);for(i=1;i<length(u);i++){v[i] = f[u[i]]};n = asorti(v,vX,"@val_num_asc");c=1;nlist="";for(x=1;x<=100;x++){while(v[vX[c]]<x && c<=n){nlist=nlist";"u[vX[c]];c++};if(c>20 || c>n){break}};print substr(nlist,2)"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}}' ChimerDB_nullomers_freq.tsv ChimerKB_n_Pub_junction_nullomers.tsv > ChimerDB-nullomersTop20.tsv
#---take first 3 and last 3 nullomers with genome-frequency > 100
awk -F "\t" '{if(NR==FNR){f[$1]=$2}else{split($1,u,";");c=0;nlist="";olist="";for(x=1;x<length(u);x++){if(f[u[x]]<300){nlist=nlist";"u[x];c++};if(c>6){split(substr(nlist,2),v,";");olist=v[1]";"v[2]";"v[3]";"v[length(v)-2]";"v[length(v)-1]";"v[length(v)]";"}else{olist=substr(nlist,2)}};print olist"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}}' ChimerDB_nullomers_freq.tsv ChimerKB_n_Pub_junction_nullomers.tsv > ChimerDB-nullomersEdge6.tsv

awk '{split($1,u,";");for(i=1;i<length(u);i++)print toupper(u[i])}' ChimerDB-nullomersEdge6.tsv > ChimerDB_nullomers.tsv
#revComp transcriptome_nullomers to search on the read 2
awk 'BEGIN{c["A"]="T";c["T"]="A";c["G"]="C";c["C"]="G";}{y="";for(j=16;j>=1;j--){y=y c[substr($1,j,1)]};print y;}' ChimerDB_nullomers.tsv > ChimerDB_nullomersRevComp.tsv

module load bowtie2
#export LD_LIBRARY_PATH=/apps/software-compiled/tbb/2018_U5-GCCcore-7.3.0/build/linux_intel64_gcc_cc7.3.0_libc2.12_kernel2.6.32_release/:$LD_LIBRARY_PATH
bowtie2-build ChimerKB_n_Pub_for_mapping.fasta ChimerDB_bowtie
