# FastNeo
This tool has been developed to detect known human neoepitopes and gene fusions in the cell-free RNA. It supports detection of neoepitopes from IEDB, TSNAdb, and neoepitopes produced via gene fusion events described in ChimerKB and ChimerPub (YE Jang et al. 2020: https://doi.org/10.1093/nar/gkz1013).

# RELEASEnotes 
Release 0.1: First release as documented in the manuscript

# SYSTEM REQUIREMENTS
x86-64 compatible processors, 4GB RAM and 64-bit Linux. The whole workflow runs on 2 cores and processes 10^10 bases (two 10GB fastq files) in ~15 minutes

# INSTALL PRE-REQUISITES

This tool requires bowtie2 (https://bowtie-bio.sourceforge.net/bowtie2/), samtools (http://www.htslib.org) and julia (https://julialang.org/downloads/). Use the following command to download and install through conda
```
conda install -c bioconda bowtie2 samtools julia
```
Or if this tools are already installed on your high performance cluster, use the following command to load them
```
module load bowtie2 samtools julia
```
Use the following code to install the required packages in Julia
```
julia pkginstall.jl
```

# INSTALL AND RUN

Download and install the tool
```
git clone https://github.com/yashumayank/FastNeo.git
cd FastNeo
bash INSTALL
export PATH=`pwd`:$PATH
```

FastNeo uses paired-end RNAseq data. The read 1 and read 2 must be is seperate fastq files must be named as {sample_name}_1.fastq and {sample_name}_2.fastq

Detect IEDB/TSNAdb neoepitopes in your sample using the following command
```
runEpitopeDB_neoepitopes.sh {sample_name}
```

Detect both IEDB/TSNAdb neoepitopes and gene fusions in your sample using the following command
```
run_neoepitopes.sh {sample_name}
```

The output is saved in the folder from where the command was run and is written to files with name {sample_name}_*_readCounts.tsv

# TEST RUN 

Search for neoepitopes in the sample dataset
```
cd test_run
run_neoepitopes.sh test_sample
```
Crosscheck the output with the expected output files that are in the test_run folder. 
Below is an example that downloads published cfRNA data using SRA toolkit and runs the tool using custom options:

```
fasterq-dump --split-3 SRR25143498
run_neoepitopes.sh -o outpie -q 10 -m 5 -c 2.7 -x 60 -y 260 SRR25143498
```

# OUTPUT format

For IEDB/TSNEdb neoepitopes, the output file shows:-

1) sample_id: fastq file name
2) gene_id: Gene ID with ENSG prefix
3) HGNC_symbol: Gene symbol
4) top_nullomer: nullomer with most read covergage among the nullomers that are associated to the mutation/mutations that produce these neoepitopes
5) neoepitopes: all wildtype-neoepitopes pairs associated to the top nullomer for the mutation:
   
   [{Wildtype epitope}->{Noepitope};] ...
6) read_count: number of mapped reads 
7) nullomer_count: number of nullomers on the read with most nullomers
8) db_name: Name of the database/databases that contain these neoepitopes.
9) annotation: known function of the protein that contains the epitopes
10) wildTypeHLA: binding affinities of all HLA alleles predicted to bind to each wildtype epitope. column output format:
 
    [Epitope:{wildtype epitope};[{HLA allele};{binding affinity},] ...] ... 
11) neoEpitopeHLA: binding affinities of all HLA alleles predicted to bind to each neoepitope column output format:

    [Epitope:{neoepitope};[{HLA allele};{binding affinity},] ...] ... 
12) GeneticAncestry: frequency of the neoepitope producing mutations in healthy individuals of various genetic ancestry groups. column format:
    
    [SNV[1|2],{chromosome};{position};{old base};{new base},AF={AF};grpmax={group with highest frequency};AF_XX={AF_XX};AF_XY={AF_XY};AF_afr={AF_afr};AF_amr={AF_amr};AF_asj={asj};AF_eas={eas};AF_fin={fin};AF_mid={mid};AF_nfe={nfe};AF_sas={sas};] ...


Output file format for the gene fusions:-

1) sample_id: fastq file name
2) ChimerKB_id: ChimerKB ID
3) top_nullomer: most covered nullomer per gene fusion
4) neoepitopes: neoepitopes associated with the gene fusion
5) read_count: number of mapped reads
6) nullomer_count: number of nullomers on the read with most nullomers
7) annotation: gene names
8) junction5: genomic loci of 5’ junction
9) junction3: genomic loci of 3’ junction.

