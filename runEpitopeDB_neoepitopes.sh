pathDS=$(pwd)/
pathDB=${path1}/nullomer_lists/
pathMI=${path1}/mapping_indices/
#path1 must have full path of the install directory eg. /usr/bin/cfRNA-neoepitopes

cd ${pathDS}
cores=2

MINQUAL=20
NULLOMERLEN=16

MINMQ=10
MINMQf=30
MAXCLIP=3
#alignments scores and aligned bases to calculate expected_score (2 * #matches - [3-6] * #mismatches)
AS1=64
AS2=277
AB1=35
AB2=150
OVERLAP=5
OUTFILE="-"
#parse options
POSITIONAL_ARGS=()

while [[ $# -gt 0 ]]; do
  case $1 in
    -m|--mapq)
      MINMQ="$2"
      shift
      shift 
      ;;
    -c|--clippedbases)
      MAXCLIP="$2"
      shift
      shift 
      ;;
    -x|--alignscore35)
      AS1="$2"
      shift
      shift 
      ;;
    -y|--alignscore150)
      AS2="$2"
      shift
      shift 
      ;;
    -q|--basequality)
      MINQUAL="$2"
      shift 
      shift 
      ;;
    -o|--outprefix)
      OUTFILE="$2"
      shift
      shift
      ;;
    -*|--*|-h)
      echo "FsatNeo v0.1 by Hemberg Lab (https://hemberg-lab.github.io)"
      echo "Usage:"
      echo "  run_neoepitopes.sh [options] {input_filename_prefix}"
      echo "  paired end data must be in 2 fastq files named as {input_filename_prefix}_1.fastq and {input_filename_prefix}_2.fastq"
      echo "  other options (default value):"
      echo "    -m|--mapq [INT] Expected alignment score filter (MES) is used if MAPQ is less than this value (10)"
      echo "    -x|--alignscore35 [INT] Minumum expected alignment score if read length = 35 nucleotides (64)"
      echo "    -y|--alignscore150 [INT] Minumum expected alignment score if read length = 150 nucleotides(277)"
      echo "    -c|--clippedbases [NUM] Minimum value of (read length) / (clipped length) (3)"
      echo "    -q|--basequality [INT] Minimum squencing quality of all the bases in the nullomer (20)"
      echo "    -o|--outprefix [STR] Prefix for the output files (input_filename_prefix])"
      exit 1
      ;;
    *)
      POSITIONAL_ARGS+=("$1") # save positional arg
      shift # past argument
      ;;
  esac
done

set -- "${POSITIONAL_ARGS[@]}" # restore positional parameters

if [[ -n $1 ]]; then
  sample=$1
  if [[ "$OUTFILE" == "-" ]]; then
    OUTFILE=$1
  fi
else
  echo "Input file required. Check help to run the command"
  echo "run_neoepitopes.sh -h"
  exit 1
fi

#SLOPE=$(echo "scale=2;(${AS2}-${AS1})/(${AB2}-${AB1})"|bc) #bc not found error
SLOPE=$(awk -v as1="${AS1}" -v as2="${AS2}" -v ab1="${AB1}" -v ab2="${AB2}" 'BEGIN { printf "%.2f", (as2-as1)/(ab2-ab1) }' </dev/null)

#julia
#/data/hemberg/shared_resources/sratoolkit.2.11.2-centos_linux64/bin/fasterq-dump --split-3 ${sample}
julia ${path1}/cfRNAnullomers.v0.2.jl -f ${pathDS}${sample}_1.fastq -n ${pathDB}epitopeDB_nullomers.tsv -l ${NULLOMERLEN} -S ${pathDS}${sample}_1.nullomers -q ${MINQUAL} --logfile ${pathDS}${sample}_1.nullomers.json -N / &
julia ${path1}/cfRNAnullomers.v0.2.jl -f ${pathDS}${sample}_2.fastq -n ${pathDB}epitopeDB_nullomersRevComp.tsv -l ${NULLOMERLEN} -S ${pathDS}${sample}_2.nullomers -q ${MINQUAL} --logfile ${pathDS}${sample}_2.nullomers.json -N /
wait
awk -v pid=${sample}_1 '{if(FNR==1){f++};if(f<=2){if(FNR%4==1){a[$1]=1}}else{if(FNR%4==1){if(a[$1]==1){x=1}else{x=0}};if(x==1){print $0 > pid ".nullomers.union.epitopeDB.fastq"}}}' ${sample}_*.nullomers_epitopeDB.fastq ${sample}_1.fastq &
awk -v pid=${sample}_2 '{if(FNR==1){f++};if(f<=2){if(FNR%4==1){a[$1]=1}}else{if(FNR%4==1){if(a[$1]==1){x=1}else{x=0}};if(x==1){print $0 > pid ".nullomers.union.epitopeDB.fastq"}}}' ${sample}_*.nullomers_epitopeDB.fastq ${sample}_2.fastq
wait
#rm ${sample}_1.fastq ${sample}_2.fastq
rm ${sample}_*.nullomers_epitopeDB.fastq
rm ${sample}_1.nullomers.json ${sample}_2.nullomers.json

#bowtie2
bowtie2 -x ${pathMI}Homo_sapiens.GRCh38.cds.all -p $cores  --no-unal --omit-sec-seq --very-sensitive-local -1 ${sample}_1.nullomers.union.epitopeDB.fastq -2 ${sample}_2.nullomers.union.epitopeDB.fastq -S ${sample}.nullomers.union.epitopeDB.sam
rm ${sample}_*.nullomers.union.epitopeDB.fastq

#samtools
samtools sort ${sample}.nullomers.union.epitopeDB.sam | samtools rmdup - ${sample}.nullomers.union.epitopeDB.bam 2>test1
rm test1 ${sample}.nullomers.union.epitopeDB.sam

#find coverage for each nullomer correspoding to specific neoepitope ; filter reads with (MAPQ < 10 AND (read_length <= (clipped_length*3) OR alignment_score <= expected_score)
#expected score calculated using linear equation (218 * aligned_length + 75)/115) allows 1 high quality mismatch (penalty=3) for mapped_length=35 and ~4 lower quality mismatches (penalty=15) for mapped_length=150
samtools view ${sample}.nullomers.union.epitopeDB.bam | awk -v "sid=${sample}" -v "minMQ=${MINMQ}" -v "maxClip=${MAXCLIP}" -v "S1=${AS1}" -v "B1=${AB1}" -v "slope=${SLOPE}" -F "\t" '{if(FNR==1){fn++; split(FILENAME,dbtmp,"/");split(dbtmp[length(dbtmp)],dbName,"_n")};if(fn==1){gEnsg[$1]=$2;gSname[$1]=$3}else{if(fn<=4){if($1!=""){split($1,u,";");split($3,uv,"->");if(null_list[$4]==""){ann[$4]=$5;;null_list[$4]=$1;for(i=1;i<=length(u);i++){null_neo[$4"_"u[i]]=$3;null_db[$4"_"u[i]]=dbName[1];null_wtMHC[$4"_"u[i]]="Epitope:"uv[1]","$6;null_neoMHC[$4"_"u[i]]="Epitope:"uv[2]","$7;null_pop[$4"_"u[i]]=$8;}}else{for(i=1;i<=length(u);i++){if(index(null_list[$4],u[i])==0){null_list[$4]=null_list[$4]";"u[i]};if(null_neo[$4"_"u[i]]==""){null_neo[$4"_"u[i]]=$3}else{if(index(null_neo[$4"_"u[i]],$3)==0){null_neo[$4"_"u[i]]=null_neo[$4"_"u[i]]";"$3}};if(null_db[$4"_"u[i]]==""){null_db[$4"_"u[i]]=dbName[1]}else{if(index(null_db[$4"_"u[i]],dbName[1])==0){null_db[$4"_"u[i]]=null_db[$4"_"u[i]]";"dbName[1]}};wtMHC_p="Epitope:"uv[1]",";neoMHC_p="Epitope:"uv[2]",";if(null_wtMHC[$4"_"u[i]]==""){null_wtMHC[$4"_"u[i]]=wtMHC_p""$6}else{if(index(null_wtMHC[$4"_"u[i]],wtMHC_p)==0){null_wtMHC[$4"_"u[i]]=null_wtMHC[$4"_"u[i]]","wtMHC_p""$6}};if(null_neoMHC[$4"_"u[i]]==""){null_neoMHC[$4"_"u[i]]=neoMHC_p""$7}else{if(index(null_neoMHC[$4"_"u[i]],neoMHC_p)==0){null_neoMHC[$4"_"u[i]]=null_neoMHC[$4"_"u[i]]","neoMHC_p""$7}};if(null_pop[$4"_"u[i]]==""){null_pop[$4"_"u[i]]=$8}else{if(index(null_pop[$4"_"u[i]],$8)==0){null_pop[$4"_"u[i]]=null_pop[$4"_"u[i]]""$8}};}}}}else{if($5<minMQ){if($12!~/^AS/){next};s=0;sc=0;split($6,u,"[A-Z]",v);for(i=1;i<=length(u);i++){if(v[i]=="S"){s=s+(u[i]);sc++}};if(length($10)<(s*maxClip)){next};trulen=length($10)-s;expS=slope*(trulen - B1) + S1; split($12,alnS,":");if(expS>alnS[3]){next}};split($3,v,".");if(et!=v[1]){et=v[1];nt=split(null_list[et],nlu,";")};for(i=1;i<=length(nlu);i++){if(index($10,nlu[i])>0){null_count[et"_"nlu[i]]++}}}}}END{for(k in null_count){split(k,w,"_");print sid"\t"w[2]"\t"null_neo[k]"\t"null_count[k]"\t"w[1]"\t"gEnsg[w[1]]"\t"gSname[w[1]]"\t"null_db[k]"\t"ann[w[1]]"\t"null_wtMHC[k]"\t"null_neoMHC[k]"\t"null_pop[k]}}' ${pathDB}ENST_ENSG_HGNC.map ${pathDB}IEDB_neoepitopes-nullomers.filtered.tsv ${pathDB}TSNAdb_TCGA_neoepitopes-nullomers.filtered.tsv ${pathDB}TSNAdb_ICGC_neoepitopes-nullomers.filtered.tsv - > ${sample}_epitopeDB_neoepitopes_counts_raw.tsv &
#Calculate coverage of the each neoepitope as the coverage of most covered corresponding nullomer; read coverage and nullomers discovered on the read are reported
sort -k1,1 -k7,7 -k3,3 -k4n -t$'\t' ${sample}_epitopeDB_neoepitopes_counts_raw.tsv | awk -F "\t" 'BEGIN{print "sample_id\tgene_id\tHGNC_symbol\ttop_nullomer\tneoepitopes\tread_count\tnullomer_count\tdb_name\tannotation\twildTypeHLA\tneoEpitopeHLA\tGeneticAncestry"}{if(neo_list!=$3 || symbol!=$7){if(FNR!=1 && rcount>1){print pid"\t"et"\t"symbol"\t"null"\t"neo_list"\t"rcount"\t"ncount"\t"dbName"\t"ann"\t"wtHLA"\t"neoHLA"\t"popVar}; pid=$1; neo_list=$3; et=$6; symbol=$7; dbName=$8; ann=$9; wtHLA=$10; neoHLA=$11; popVar=$12; ncount=0}; ncount++;rcount=$4; null=$2}END{if(rcount>1){print pid"\t"et"\t"symbol"\t"null"\t"neo_list"\t"rcount"\t"ncount"\t"dbName"\t"ann"\t"wtHLA"\t"neoHLA"\t"popVar}}' > ${OUTFILE}_neoepitopes.tsv &

rm ${sample}.nullomers.union.epitopeDB.bam ${sample}_epitopeDB_neoepitopes_counts_raw.tsv
