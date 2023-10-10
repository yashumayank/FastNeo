#!/bin/bash
#BSUB -J epiDB-Null
#BSUB -o test-%J.out
#BSUB -e test-%J.err
#BSUB -q bigmem
#BSUB -n 1

eval "$(/PHShome/mm1169/miniconda3/bin/conda shell.bash hook)"
conda activate scanpy
cd /data/hemberg/nullomers/IEDB/epitope_DBs/
#Run R script to extract neoepitopes from the IEDB database
Rscript extractIEDBneoepitopes.R

#create maps from protein id (uniprot) to ENST id (ensembleUgencode)
awk -F "\t" '{gene="";tran="";for(x=1;x<=NF;x++){if($x~/ENSG/){gene=$x}else if($x~/ENST/){tran=$x}}; print $1"\t"gene"\t"tran}' HUMAN_9606_idmapping_selected.tab > uniprot2ensembl_map.tab
#awk '{a[$2]=a[$2]"; "$1}END{for(x in a)print x"\t\t"substr(a[x],3)}' gencode.v40.metadata.TrEMBL gencode.v40.metadata.SwissProt > uniprot2gencode.tab
#merge maps from protein id (uniprot) to ENST id (ensembleUgencode)
awk -F "\t" '{split($2,genes,"; ");split($3,trans,"; ");for(i in genes)if(g[$1]!~genes[i])g[$1]=g[$1]";"genes[i];for(i in trans)if(t[$1]!~trans[i])t[$1]=t[$1]";"trans[i]}END{for(x in t){print x"\t"substr(g[x],2)"\t"substr(t[x],2)}}'  uniprot2ensembl_map.tab  > uniprot2ensembl_map2.tab

#convert IEDB neoepitopes to the format required for downstream analysis
awk -F "\t" '(NR==FNR){g[$1]=$2;t[$1]=$3;next}{if(t[$6]!=""){split(g[$6],genes,";");\
split(t[$6],trans,";");for(i in trans){split(trans[i],j,".");plen=length($1);\
trpep[j[1]]=trpep[j[1]]";"$1;trlen[j[1]]=trlen[j[1]]";"plen;trWT[j[1]]=trWT[j[1]]";"$5}\
}else{print $6 " not found ERROR." > "convertID.err"}}\
END{for(x in trpep)print x"\t"substr(trpep[x],2)"\t"substr(trlen[x],2)"\t"substr(trWT[x],2)}' \
uniprot2ensembl_map2.tab IEDB_neoepitopes_sapien.tab > IEDB_neoepitopes_per_ENST2.tab
#convert TSNAdb neoepitopes to the format required for downstream analysis
awk '(FNR!=1 ){trpep[$4]=trpep[$4]";"$14;trlen[$4]=trlen[$4]";"length($14);trWT[$4]=trWT[$4]";"$10;}\
END{for(x in trpep)print x"\t"substr(trpep[x],2)"\t"substr(trlen[x],2)"\t"substr(trWT[x],2)}' \
TSNAdb_frequent_neoantigen_TCGA_4.0_adj.txt > TSNAdb_frequent_TCGA_per_ENST3.tab

#--extract Human cds 16mers
if [ ! -e "Homo_sapiens_cds_16mers.tab" ] || [ ! -s "Homo_sapiens_cds_16mers.tab" ] ; then
  awk '{if($0~/^>.*/){ln=length(x);for(i=1;i+15<=ln;i++){a[substr(x,i,16)]+=1};x=""}else{x=x $1}}END{for(j in a){print j "\t" a[j]}}' ../Homo_sapiens.GRCh38.cds.all.fa > Homo_sapiens_cds_16mers.tab
fi
#--extract extract nullomers associated to neoepitopes in IEDB and TSNAdb
python IEDB_TSNAdb2nullomer.py TSNAdb_frequent_ICGC_per_ENST3.tab Homo_sapiens_cds_16mers.tab ../Homo_sapiens.GRCh38.cds.all.fa TSNAdb_ICGC
python IEDB_TSNAdb2nullomer.py TSNAdb_frequent_TCGA_per_ENST3.tab Homo_sapiens_cds_16mers.tab ../Homo_sapiens.GRCh38.cds.all.fa TSNAdb_TCGA
python IEDB_TSNAdb2nullomer.py IEDB_neoepitopes_per_ENST2.tab Homo_sapiens_cds_16mers.tab ../Homo_sapiens.GRCh38.cds.all.fa IEDB

#---Count the WT-neoepitope pairs extracted from each database
echo "TSNAdb_ICGC"
awk 'NF==4' TSNAdb_ICGC_neoepitopes-nullomers.tsv |cut -f3 |sort|uniq -c|sort|wc -l
echo "TSNAdb_TCGA"
awk 'NF==4' TSNAdb_TCGA_neoepitopes-nullomers.tsv |cut -f3 |sort|uniq -c|sort|wc -l
echo "IEDB"
awk 'NF==4' IEDB_neoepitopes-nullomers.tsv |cut -f3 |sort|uniq -c|sort|wc -l

#--merge the nullomers associated to all the above neoepitopes
awk -F "\t" '{split($1,u,";");for(i in u){a[u[i]]=1}}END{for(i in a)print i}' IEDB_neoepitopes-nullomers.tsv TSNAdb_TCGA_neoepitopes-nullomers.tsv TSNAdb_ICGC_neoepitopes-nullomers.tsv > epitopeDB_nullomers.all.tsv

#--Calculate genome frequency of each nullomer
awk -F "\t" '(NR==FNR){a[$1]=1;next}{if($1~/^>.*/){if(seqid!=""){for(i=1;i<=(length(seq)-15);i++){xtr=substr(seq,i,16);if(xtr in a){a[xtr]++}}};seqid=$1;seq=""}else{seq=seq""$1}}END{for(j=1;j<=(length(seq)-15);j++){xtr=substr(seq,j,16);if(xtr in a){a[xtr]++}};for(k in a){print k"\t"a[k]-1}}' epitopeDB_nullomers.all.tsv /data/hemberg/shared_resources/genomes/human/GRCh38.p13.107.genome.fa > epitopeDB_nullomers_freq.tsv

#---take 20 least frequent per neoepitope; filter out nullomers with genome-frequency > 100
#awk -F "\t" 'BEGIN{PROCINFO["sorted_in"] = "@ind_num_asc"}{if(NR==FNR){f[$1]=$2}else{split($1,u,";");delete(v);delete(vX);for(i=1;i<=length(u);i++){v[i] = f[u[i]]};n = asorti(v,vX,"@val_num_asc");c=1;nlist="";for(x=1;x<=100;x++){while(v[vX[c]]<x && c<=n){nlist=nlist";"u[vX[c]];c++};if(c>20 || c>n){break}};print substr(nlist,2)"\t"$2"\t"$3"\t"$4"\t"$5}}' epitopeDB_nullomers_freq.tsv TSNAdb_ICGC_neoepitopes-nullomers.tsv > TSNAdb_ICGC_neoepitopes-nullomersTop20.tsv
#awk -F "\t" 'BEGIN{PROCINFO["sorted_in"] = "@ind_num_asc"}{if(NR==FNR){f[$1]=$2}else{split($1,u,";");delete(v);delete(vX);for(i=1;i<=length(u);i++){v[i] = f[u[i]]};n = asorti(v,vX,"@val_num_asc");c=1;nlist="";for(x=1;x<=100;x++){while(v[vX[c]]<x && c<=n){nlist=nlist";"u[vX[c]];c++};if(c>20 || c>n){break}};print substr(nlist,2)"\t"$2"\t"$3"\t"$4"\t"$5}}' epitopeDB_nullomers_freq.tsv TSNAdb_TCGA_neoepitopes-nullomers.tsv > TSNAdb_TCGA_neoepitopes-nullomersTop20.tsv
#awk -F "\t" 'BEGIN{PROCINFO["sorted_in"] = "@ind_num_asc"}{if(NR==FNR){f[$1]=$2}else{split($1,u,";");delete(v);delete(vX);for(i=1;i<=length(u);i++){v[i] = f[u[i]]};n = asorti(v,vX,"@val_num_asc");c=1;nlist="";for(x=1;x<=100;x++){while(v[vX[c]]<x && c<=n){nlist=nlist";"u[vX[c]];c++};if(c>20 || c>n){break}};print substr(nlist,2)"\t"$2"\t"$3"\t"$4"\t"$5}}' epitopeDB_nullomers_freq.tsv IEDB_neoepitopes-nullomers.tsv > IEDB_neoepitopes-nullomersTop20.tsv

#---For nullomers with genome-frequency > 100, keep first 3 and last 3 nullomers per DNA confirmation of each neoepitope
awk -F "\t" '{if(NR==FNR){f[$1]=$2}else{split($1,u,";");c=0;nlist="";olist="";for(x=1;x<=length(u);x++){if(f[u[x]]<300){nlist=nlist";"u[x];c++};if(c>6){split(substr(nlist,2),v,";");olist=v[1]";"v[2]";"v[3]";"v[length(v)-2]";"v[length(v)-1]";"v[length(v)]}else{olist=substr(nlist,2)}};print olist"\t"$2"\t"$3"\t"$4"\t"$5}}' epitopeDB_nullomers_freq.tsv TSNAdb_ICGC_neoepitopes-nullomers.tsv > TSNAdb_ICGC_neoepitopes-nullomersEdge6.tsv
awk -F "\t" '{if(NR==FNR){f[$1]=$2}else{split($1,u,";");c=0;nlist="";olist="";for(x=1;x<=length(u);x++){if(f[u[x]]<300){nlist=nlist";"u[x];c++};if(c>6){split(substr(nlist,2),v,";");olist=v[1]";"v[2]";"v[3]";"v[length(v)-2]";"v[length(v)-1]";"v[length(v)]}else{olist=substr(nlist,2)}};print olist"\t"$2"\t"$3"\t"$4"\t"$5}}' epitopeDB_nullomers_freq.tsv TSNAdb_TCGA_neoepitopes-nullomers.tsv > TSNAdb_TCGA_neoepitopes-nullomersEdge6.tsv
awk -F "\t" '{if(NR==FNR){f[$1]=$2}else{split($1,u,";");c=0;nlist="";olist="";for(x=1;x<=length(u);x++){if(f[u[x]]<300){nlist=nlist";"u[x];c++};if(c>6){split(substr(nlist,2),v,";");olist=v[1]";"v[2]";"v[3]";"v[length(v)-2]";"v[length(v)-1]";"v[length(v)]}else{olist=substr(nlist,2)}};print olist"\t"$2"\t"$3"\t"$4"\t"$5}}' epitopeDB_nullomers_freq.tsv IEDB_neoepitopes-nullomers.tsv > IEDB_neoepitopes-nullomersEdge6.tsv

#--merge the filtered nullomers associated to IEDB and TSNAdb neoepitopes
awk -F "\t" '($1!=""){split($1,u,";");for(i in u){a[u[i]]=1}}END{for(i in a)print i}' IEDB_neoepitopes-nullomersEdge6.tsv TSNAdb_TCGA_neoepitopes-nullomersEdge6.tsv TSNAdb_ICGC_neoepitopes-nullomersEdge6.tsv > epitopeDB_nullomers.tsv
#--revComp cds_nullomers to search on the read 2 of the RNAseq data
awk 'BEGIN{c["A"]="T";c["T"]="A";c["G"]="C";c["C"]="G";}{y="";for(j=16;j>=1;j--){y=y c[substr($1,j,1)]};print y;}' epitopeDB_nullomers.tsv > epitopeDB_nullomersRevComp.tsv

#nulomers  DNA  epitope  ENST_id  annotation  HLA  wild_Affinity  mut_Affinity
sort TSNAdb_ICGC_neoepitopes-nullomersEdge6.tsv|uniq | awk -F "\t" '(NR==FNR){hla[$4"_"$10"->"$14]=$9; wt[$4"_"$10"->"$14]=$11; neo[$4"_"$10"->"$14]=$15;next}{OFS="\t";print $0"\t"hla[$4"_"$3] "\t" wt[$4"_"$3] "\t" neo[$4"_"$3]}' TSNAdb_frequent_neoantigen_ICGC_4.0_adj.txt - > TSNAdb_ICGC_neoepitopes-nullomers.filtered.tsv
sort TSNAdb_TCGA_neoepitopes-nullomersEdge6.tsv|uniq | awk -F "\t" '(NR==FNR){hla[$4"_"$10"->"$14]=$9; wt[$4"_"$10"->"$14]=$11; neo[$4"_"$10"->"$14]=$15;next}{OFS="\t";print $0"\t"hla[$4"_"$3] "\t" wt[$4"_"$3] "\t" neo[$4"_"$3]}' TSNAdb_frequent_neoantigen_TCGA_4.0_adj.txt - > TSNAdb_TCGA_neoepitopes-nullomers.filtered.tsv
sort IEDB_neoepitopes-nullomersEdge6.tsv|uniq > IEDB_neoepitopes-nullomers.filtered.tsv

