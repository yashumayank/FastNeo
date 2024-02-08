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
awk -F "\t" '(FNR!=1 ){trpep[$4]=trpep[$4]";"$14;trlen[$4]=trlen[$4]";"length($14);trWT[$4]=trWT[$4]";"$10;}\
END{for(x in trpep)print x"\t"substr(trpep[x],2)"\t"substr(trlen[x],2)"\t"substr(trWT[x],2)}' TSNAdb_frequent_neoantigen_TCGA_4.0_adj.txt > TSNAdb_frequent_TCGA_per_ENST3.tab
awk -F "\t" '(FNR!=1 ){trpep[$4]=trpep[$4]";"$14;trlen[$4]=trlen[$4]";"length($14);trWT[$4]=trWT[$4]";"$10;}\
END{for(x in trpep)print x"\t"substr(trpep[x],2)"\t"substr(trlen[x],2)"\t"substr(trWT[x],2)}' TSNAdb_frequent_neoantigen_ICGC_4.0_adj.txt > TSNAdb_frequent_ICGC_per_ENST3.tab

# crete a list of all epitopes and run netMHC I (<=11 aa) & II (>=12 aa) on all epitopes 
awk '{split($2,u,";");split($4,v,";"); for(i=1;i<=length(u);i++){if(length(u[i])!=length(v[i])){print u[i]"\n"v[i] > "IEDB_neoepitopes_per_ENST2.unequal.tab"}else{\
if(length(u[i])<=11){a[u[i]]=length(u[i]);a[v[i]]=length(u[i])}else{a[u[i]]="12+";a[v[i]]="12+"}}}}END{for(x in a) print x > "IEDB_neoepitopes_per_ENST2."a[x]"mers.tab"}' IEDB_neoepitopes_per_ENST2.tab
awk '{split($2,u,";");split($4,v,";"); for(i=1;i<=length(u);i++){if(length(u[i])!=length(v[i])){print u[i]"\n"v[i] > "TSNAdb_neoepitopes_per_ENST2.unequal.tab"}else{\
if(length(u[i])<=11){a[u[i]]=length(u[i]);a[v[i]]=length(u[i])}else{a[u[i]]="12+";a[v[i]]="12+"}}}}END{for(x in a) print x > "TSNAdb_neoepitopes_per_ENST2."a[x]"mers.tab"}' TSNAdb_frequent_ICGC_per_ENST3.tab TSNAdb_frequent_TCGA_per_ENST3.tab

netMHCpan -BA -s -l 8 \
-a HLA-A01:01,HLA-A02:01,HLA-A03:01,HLA-A11:01,HLA-A24:02,HLA-B07:02,HLA-B35:01,HLA-B40:01,HLA-B51:01,HLA-C01:02,HLA-C03:03,HLA-C03:04,HLA-C04:01,HLA-C06:02,HLA-C07:01,HLA-C07:02,HLA-C08:01,HLA-A33:03,HLA-B08:01 \
-p  IEDB_neoepitopes_per_ENST2.8mers.tab | awk '($18=="WB"||$18=="SB"){print $3"\t"$2"\t"$16}' > IEDB_neoantigen_candidate_8mers_netMHC.out
netMHCpan -BA -s -l 9 \
-a HLA-A01:01,HLA-A02:01,HLA-A03:01,HLA-A11:01,HLA-A24:02,HLA-B07:02,HLA-B35:01,HLA-B40:01,HLA-B51:01,HLA-C01:02,HLA-C03:03,HLA-C03:04,HLA-C04:01,HLA-C06:02,HLA-C07:01,HLA-C07:02,HLA-C08:01,HLA-A33:03,HLA-B08:01 \
-p  IEDB_neoepitopes_per_ENST2.9mers.tab | awk '($18=="WB"||$18=="SB"){print $3"\t"$2"\t"$16}' > IEDB_neoantigen_candidate_9mers_netMHC.out
netMHCpan -BA -s -l 10 \
-a HLA-A01:01,HLA-A02:01,HLA-A03:01,HLA-A11:01,HLA-A24:02,HLA-B07:02,HLA-B35:01,HLA-B40:01,HLA-B51:01,HLA-C01:02,HLA-C03:03,HLA-C03:04,HLA-C04:01,HLA-C06:02,HLA-C07:01,HLA-C07:02,HLA-C08:01,HLA-A33:03,HLA-B08:01 \
-p  IEDB_neoepitopes_per_ENST2.10mers.tab | awk '($18=="WB"||$18=="SB"){print $3"\t"$2"\t"$16}' > IEDB_neoantigen_candidate_10mers_netMHC.out
netMHCpan -BA -s -l 11 \
-a HLA-A01:01,HLA-A02:01,HLA-A03:01,HLA-A11:01,HLA-A24:02,HLA-B07:02,HLA-B35:01,HLA-B40:01,HLA-B51:01,HLA-C01:02,HLA-C03:03,HLA-C03:04,HLA-C04:01,HLA-C06:02,HLA-C07:01,HLA-C07:02,HLA-C08:01,HLA-A33:03,HLA-B08:01 \
-p  IEDB_neoepitopes_per_ENST2.11mers.tab | awk '($18=="WB"||$18=="SB"){print $3"\t"$2"\t"$16}' > IEDB_neoantigen_candidate_11mers_netMHC.out
netMHCIIpan -BA -s -inptype 1 \
-a DRB1_0701,DRB1_1501,DRB1_0301,DRB1_1101,DRB1_0101,DRB1_1302,DRB1_1301,DRB1_1502,DRB1_0401,DRB1_1201,DRB1_0403,DRB4_0101,DRB3_0202,DRB3_0101,DRB5_0101,DRB3_0301,HLA-DQA10102-DQB10501,HLA-DQA10103-DQB10501,HLA-DQA10501-DQB10201,HLA-DPA10103-DPB10201,HLA-DPA10103-DPB10402 \
-f IEDB_neoepitopes_per_ENST2.12+mers.tab | awk '($10<=10 && $1>0 && $1<61){print $3"\t"$2"\t"$14}' > IEDB_neoantigen_candidate_12+mers_netMHC.out

netMHCpan -BA -s -l 8 \
-a HLA-A01:01,HLA-A02:01,HLA-A03:01,HLA-A11:01,HLA-A24:02,HLA-B07:02,HLA-B35:01,HLA-B40:01,HLA-B51:01,HLA-C01:02,HLA-C03:03,HLA-C03:04,HLA-C04:01,HLA-C06:02,HLA-C07:01,HLA-C07:02,HLA-C08:01,HLA-A33:03,HLA-B08:01 \
-p  TSNAdb_neoepitopes_per_ENST2.8mers.tab | awk '($18=="WB"||$18=="SB"){print $3"\t"$2"\t"$16}' > TSNAdb_neoantigen_candidate_8mers_netMHC.out
netMHCpan -BA -s -l 9 \
-a HLA-A01:01,HLA-A02:01,HLA-A03:01,HLA-A11:01,HLA-A24:02,HLA-B07:02,HLA-B35:01,HLA-B40:01,HLA-B51:01,HLA-C01:02,HLA-C03:03,HLA-C03:04,HLA-C04:01,HLA-C06:02,HLA-C07:01,HLA-C07:02,HLA-C08:01,HLA-A33:03,HLA-B08:01 \
-p  TSNAdb_neoepitopes_per_ENST2.9mers.tab | awk '($18=="WB"||$18=="SB"){print $3"\t"$2"\t"$16}' > TSNAdb_neoantigen_candidate_9mers_netMHC.out
netMHCpan -BA -s -l 10 \
-a HLA-A01:01,HLA-A02:01,HLA-A03:01,HLA-A11:01,HLA-A24:02,HLA-B07:02,HLA-B35:01,HLA-B40:01,HLA-B51:01,HLA-C01:02,HLA-C03:03,HLA-C03:04,HLA-C04:01,HLA-C06:02,HLA-C07:01,HLA-C07:02,HLA-C08:01,HLA-A33:03,HLA-B08:01 \
-p  TSNAdb_neoepitopes_per_ENST2.10mers.tab | awk '($18=="WB"||$18=="SB"){print $3"\t"$2"\t"$16}' > TSNAdb_neoantigen_candidate_10mers_netMHC.out
netMHCpan -BA -s -l 11 \
-a HLA-A01:01,HLA-A02:01,HLA-A03:01,HLA-A11:01,HLA-A24:02,HLA-B07:02,HLA-B35:01,HLA-B40:01,HLA-B51:01,HLA-C01:02,HLA-C03:03,HLA-C03:04,HLA-C04:01,HLA-C06:02,HLA-C07:01,HLA-C07:02,HLA-C08:01,HLA-A33:03,HLA-B08:01 \
-p  TSNAdb_neoepitopes_per_ENST2.11mers.tab | awk '($18=="WB"||$18=="SB"){print $3"\t"$2"\t"$16}' > TSNAdb_neoantigen_candidate_11mers_netMHC.out

awk '{if(b[$1"_"$2"_"$3]!=1){a[$1]=a[$1]","$2";"$3;b[$1"_"$2"_"$3]=1}}END{for(i in a){print i"\t"substr(a[i],2)}}' IEDB_neoantigen_candidate_*mers_netMHC.out > IEDB_neoepitope_bindingAffinity_netMHC.out
awk '{if(b[$1"_"$2"_"$3]!=1){a[$1]=a[$1]","$2";"$3;b[$1"_"$2"_"$3]=1}}END{for(i in a){print i"\t"substr(a[i],2)}}' TSNAdb_neoantigen_candidate_*mers_netMHC.out > TSNAdb_neoepitope_bindingAffinity_netMHC.out

#--extract Human cds 16mers
if [ ! -e "Homo_sapiens_cds_16mers.tab" ] || [ ! -s "Homo_sapiens_cds_16mers.tab" ] ; then
  awk '{if($0~/^>.*/){ln=length(x);for(i=1;i+15<=ln;i++){a[substr(x,i,16)]+=1};x=""}else{x=x $1}}END{for(j in a){print j "\t" a[j]}}' ../Homo_sapiens.GRCh38.cds.all.fa > Homo_sapiens_cds_16mers.tab
fi

#get HGNC gene_symbols and Ensembl gene ids
#grep ">" ../Homo_sapiens.GRCh38.cds.all.fa|awk '{gn=" ";en=" ";split($0,u,"[ :>]");for(i=1;i<=length(u);i++){if(u[i]~/ENST/){split(u[i],v,".");en=v[1]};if(u[i]~/ENSG/){split(u[i],w,".");gn=w[1]}};print en"\t"gn}' > ENST2ENSG.map
#awk -F "\t" '(FNR==NR){gn[substr($8,1,15)]=$2; etn[substr($7,1,15)]=$2;next}{if(gn[$2]!=""){print $1"\t"$2"\t"gn[$2]}else{print $1"\t"$2"\t"etn[$2]}}' HGNCnames.tab ENST2ENSG.map > ENST_ENSG_HGNC.map
grep ">" ../Homo_sapiens.GRCh38.cds.all.fa|awk '{gn="";en="";sn="";split($0,u,"[ :>]");for(i=1;i<=length(u);i++){if(u[i]~/ENST/){split(u[i],v,".");en=v[1]};if(u[i]~/ENSG/){split(u[i],w,".");gn=w[1]};if(u[i]~/gene_symbol/){sn=u[i+1]}};if(sn==""){sn=gn};print en"\t"gn"\t"sn}' > ENST_ENSG_HGNC.map

#extract population variants from gnomad database ('AF' field > 1e-07 and 'vep' field must have Consequence == missense_variant, BIOTYPE == protein_coding, Feature ~ /ENST/)
module load bcftools
for i in $(seq 1 22) X Y ;do
  curl -O https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/exomes/gnomad.exomes.v4.0.sites.chr${1}.vcf.bgz
  bcftools view gnomad.exomes.v4.0.sites.chr${1}.vcf.bgz | grep "missense" |  awk -F "\t" '{split($8,u,"[;]");split(u[3],v,"=");if(v[2] > 1e-07){split(u[499],w,",");for(i=1;i<=length(w);i++){split(w[i],x,"|");if(x[2]=="missense_variant" && x[8]=="protein_coding"){if(x[7]~/ENST/){pid=x[7]}else{next};print pid"\t"x[14]"\t"x[15]"\t"x[16]"\t"x[17]"\t"$1";"$2";"$4";"$5"\t"u[3]";"u[4]";"u[8]";"u[12]";"u[25]";"u[37]";"u[49]";"u[61]";"u[73]";"u[85]";"u[97]";"u[249]}}}}' > gnomad.exomes.v4.0.sites.chr${1}.vep2
  rm gnomad.exomes.v4.0.sites.chr${1}.vcf.bgz
done
cat gnomad.exomes.v4.0.sites.chr*.vep2 > gnomad.exomes_all_missense.vep

#--extract extract nullomers associated to neoepitopes in IEDB and TSNAdb
python IEDB_TSNAdb2nullomer.py IEDB_neoepitopes_per_ENST2.tab Homo_sapiens_cds_16mers.tab IEDB_neoepitope_bindingAffinity_netMHC.out gnomad.exomes_all_missense.vep ../Homo_sapiens.GRCh38.cds.all.fa IEDB
python IEDB_TSNAdb2nullomer.py TSNAdb_frequent_ICGC_per_ENST3.tab Homo_sapiens_cds_16mers.tab TSNAdb_neoepitope_bindingAffinity_netMHC.out gnomad.exomes_all_missense.vep ../Homo_sapiens.GRCh38.cds.all.fa TSNAdb_ICGC
python IEDB_TSNAdb2nullomer.py TSNAdb_frequent_TCGA_per_ENST3.tab Homo_sapiens_cds_16mers.tab TSNAdb_neoepitope_bindingAffinity_netMHC.out gnomad.exomes_all_missense.vep ../Homo_sapiens.GRCh38.cds.all.fa TSNAdb_TCGA

#---Count the WT-neoepitope pairs extracted from each database
echo "TSNAdb_ICGC"
awk '{split($2,u,";");split($4,v,";"); for(i=1;i<=length(u);i++){a[u[i]"_"v[i]]=1}}END{for(x in a){print x} }'  TSNAdb_frequent_ICGC_per_ENST3.tab|wc -l
cut -f3 TSNAdb_ICGC_neoepitopes-nullomers.tsv |sort|uniq -c|sort -gk1,1 |wc -l
echo "TSNAdb_TCGA"
awk '{split($2,u,";");split($4,v,";"); for(i=1;i<=length(u);i++){a[u[i]"_"v[i]]=1}}END{for(x in a){print x} }'  TSNAdb_frequent_TCGA_per_ENST3.tab|wc -l
cut -f3 TSNAdb_TCGA_neoepitopes-nullomers.tsv |sort|uniq -c|sort -gk1,1 |wc -l
echo "IEDB"
awk '{split($2,u,";");split($4,v,";"); for(i=1;i<=length(u);i++){a[u[i]"_"v[i]]=1}}END{for(x in a){print x} }' IEDB_neoepitopes_per_ENST2.tab|wc -l
cut -f3 IEDB_neoepitopes-nullomers.tsv |sort|uniq -c|sort -gk1,1 |wc -l

#--merge the nullomers associated to all the above neoepitopes
awk -F "\t" '{split($1,u,";");for(i in u){a[u[i]]=1}}END{for(i in a)print i}' IEDB_neoepitopes-nullomers.tsv TSNAdb_TCGA_neoepitopes-nullomers.tsv TSNAdb_ICGC_neoepitopes-nullomers.tsv > epitopeDB_nullomers.all.tsv

#--Calculate genome frequency of each nullomer
zcat /data/hemberg/shared_resources/genomes/human/GRCh38.p13.107.genome.fa.gz| awk -F "\t" 'BEGIN{c["A"]="T";c["T"]="A";c["G"]="C";c["C"]="G";}{if(NR==FNR){a[$1]=1;y="";for(j=16;j>=1;j--){y=y c[substr($1,j,1)]};rv[y]=$1;next};if($1~/^>.*/){if(seqid!=""){for(i=1;i<=(length(seq)-15);i++){xtr=substr(seq,i,16);if(xtr in a){a[xtr]++};if(xtr in rv){a[rv[xtr]]++}}};seqid=$1;seq=""}else{seq=seq""$1}}END{for(j=1;j<=(length(seq)-15);j++){xtr=substr(seq,j,16);if(xtr in a){a[xtr]++};if(xtr in rv){a[rv[xtr]]++}};for(k in a){print k"\t"a[k]-1}}' epitopeDB_nullomers.all.tsv - > epitopeDB_nullomers_freq.tsv

#---take 10 least frequent per neoepitope; filter out nullomers with genome-frequency > 300
#awk -F "\t" 'BEGIN{PROCINFO["sorted_in"] = "@ind_num_asc"}{if(NR==FNR){f[$1]=$2}else{split($1,u,";");delete(v);delete(vX);for(i=1;i<=length(u);i++){v[i] = f[u[i]]};n = asorti(v,vX,"@val_num_asc");c=1;nlist="";for(x=1;x<=300;x++){while(v[vX[c]]<x && c<=n){nlist=nlist";"u[vX[c]];c++};if(c>10 || c>n){break}};print substr(nlist,2)"\t"$2"\t"$3"\t"$4"\t"$5}}' epitopeDB_nullomers_freq.tsv TSNAdb_ICGC_neoepitopes-nullomers.tsv > TSNAdb_ICGC_neoepitopes-nullomersTop20.tsv
#awk -F "\t" 'BEGIN{PROCINFO["sorted_in"] = "@ind_num_asc"}{if(NR==FNR){f[$1]=$2}else{split($1,u,";");delete(v);delete(vX);for(i=1;i<=length(u);i++){v[i] = f[u[i]]};n = asorti(v,vX,"@val_num_asc");c=1;nlist="";for(x=1;x<=300;x++){while(v[vX[c]]<x && c<=n){nlist=nlist";"u[vX[c]];c++};if(c>10 || c>n){break}};print substr(nlist,2)"\t"$2"\t"$3"\t"$4"\t"$5}}' epitopeDB_nullomers_freq.tsv TSNAdb_TCGA_neoepitopes-nullomers.tsv > TSNAdb_TCGA_neoepitopes-nullomersTop20.tsv
#awk -F "\t" 'BEGIN{PROCINFO["sorted_in"] = "@ind_num_asc"}{if(NR==FNR){f[$1]=$2}else{split($1,u,";");delete(v);delete(vX);for(i=1;i<=length(u);i++){v[i] = f[u[i]]};n = asorti(v,vX,"@val_num_asc");c=1;nlist="";for(x=1;x<=300;x++){while(v[vX[c]]<x && c<=n){nlist=nlist";"u[vX[c]];c++};if(c>10 || c>n){break}};print substr(nlist,2)"\t"$2"\t"$3"\t"$4"\t"$5}}' epitopeDB_nullomers_freq.tsv IEDB_neoepitopes-nullomers.tsv > IEDB_neoepitopes-nullomersTop20.tsv

#---For nullomers with genome-frequency < 300, keep first 3 and last 3 nullomers per DNA confirmation of each neoepitope
awk -F "\t" '{if(NR==FNR){f[$1]=$2}else{split($1,u,";");c=0;nlist="";olist="";for(x=1;x<=length(u);x++){if(f[u[x]]<300){nlist=nlist";"u[x];c++};if(c>6){split(substr(nlist,2),v,";");olist=v[1]";"v[2]";"v[3]";"v[length(v)-2]";"v[length(v)-1]";"v[length(v)]}else{olist=substr(nlist,2)}};print olist"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}}' epitopeDB_nullomers_freq.tsv TSNAdb_ICGC_neoepitopes-nullomers.tsv > TSNAdb_ICGC_neoepitopes-nullomersEdge6.tsv
awk -F "\t" '{if(NR==FNR){f[$1]=$2}else{split($1,u,";");c=0;nlist="";olist="";for(x=1;x<=length(u);x++){if(f[u[x]]<300){nlist=nlist";"u[x];c++};if(c>6){split(substr(nlist,2),v,";");olist=v[1]";"v[2]";"v[3]";"v[length(v)-2]";"v[length(v)-1]";"v[length(v)]}else{olist=substr(nlist,2)}};print olist"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}}' epitopeDB_nullomers_freq.tsv TSNAdb_TCGA_neoepitopes-nullomers.tsv > TSNAdb_TCGA_neoepitopes-nullomersEdge6.tsv
awk -F "\t" '{if(NR==FNR){f[$1]=$2}else{split($1,u,";");c=0;nlist="";olist="";for(x=1;x<=length(u);x++){if(f[u[x]]<300){nlist=nlist";"u[x];c++};if(c>6){split(substr(nlist,2),v,";");olist=v[1]";"v[2]";"v[3]";"v[length(v)-2]";"v[length(v)-1]";"v[length(v)]}else{olist=substr(nlist,2)}};print olist"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}}' epitopeDB_nullomers_freq.tsv IEDB_neoepitopes-nullomers.tsv > IEDB_neoepitopes-nullomersEdge6.tsv

#--merge the filtered nullomers associated to IEDB and TSNAdb neoepitopes
awk -F "\t" '($1!=""){split($1,u,";");for(i in u){a[u[i]]=1}}END{for(i in a)print i}' IEDB_neoepitopes-nullomersEdge6.tsv TSNAdb_TCGA_neoepitopes-nullomersEdge6.tsv TSNAdb_ICGC_neoepitopes-nullomersEdge6.tsv > epitopeDB_nullomers.tsv
#--revComp cds_nullomers to search on the read 2 of the RNAseq data
awk 'BEGIN{c["A"]="T";c["T"]="A";c["G"]="C";c["C"]="G";}{y="";for(j=16;j>=1;j--){y=y c[substr($1,j,1)]};print y;}' epitopeDB_nullomers.tsv > epitopeDB_nullomersRevComp.tsv

#nulomers  DNA  epitope  ENST_id  annotation  HLA  wild_Affinity  mut_Affinity
#sort TSNAdb_ICGC_neoepitopes-nullomers.tsv|uniq | awk -F "\t" '(NR==FNR){hla[$4"_"$10"->"$14]=$9; wt[$4"_"$10"->"$14]=$11; neo[$4"_"$10"->"$14]=$15;next}{OFS="\t";print $0"\t"hla[$4"_"$3] "\t" wt[$4"_"$3] "\t" neo[$4"_"$3]}' TSNAdb_frequent_neoantigen_ICGC_4.0_adj.txt - > TSNAdb_ICGC_neoepitopes-nullomers.filtered.tsv
sort TSNAdb_ICGC_neoepitopes-nullomersEdge6.tsv|uniq > TSNAdb_ICGC_neoepitopes-nullomers.filtered.tsv 
sort TSNAdb_TCGA_neoepitopes-nullomersEdge6.tsv|uniq > TSNAdb_TCGA_neoepitopes-nullomers.filtered.tsv
sort IEDB_neoepitopes-nullomersEdge6.tsv|uniq > IEDB_neoepitopes-nullomers.filtered.tsv


