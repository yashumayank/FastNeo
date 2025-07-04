# FastNeo
FastNeo detects known human neoepitopes and gene fusions in the stranded bulk RNA-seq data. It is especially suited for the cell-free RNA data with short fragmented sequences and low allele coverage. FastNeo detects known neoepitopes that are characterised in IEDB (https://www.iedb.org) and TSNAdb (http://biopharm.zju.edu.cn/). It also detects neoepitopes produced via the gene fusion events described in ChimerKB and ChimerPub (https://www.kobic.re.kr/chimerdb/chimerkb). The neoepitopes, gene fusions, and corrresponding metadata are listed in the files provided in the nullomer_lists folder.

#### RELEASE NOTES 
Release 0.1: First release as documented in the manuscript, https://doi.org/10.1093/bioinformatics/btaf138

#### SYSTEM REQUIREMENTS
x86-64 compatible processors, 4GB RAM and 64-bit Linux. The whole workflow runs on 2 cores and processes 10^10 bases (two 10GB fastq files) in ~15 minutes

### PREREQUISITES

bowtie2 (https://bowtie-bio.sourceforge.net/bowtie2/)
samtools (http://www.htslib.org)
julia (https://julialang.org/downloads/) 
The following command can be used to download and install the prerequisites via conda :-
```
conda install -c bioconda bowtie2 samtools julia
```

In case of a high performance cluster (HPC), some or all of these prerequisites might already be installed on the (HPC). Please check if the modules are already installed using the following command: `module spider bowtie2 samtools julia` and load the ones that are installed using the following command `module load {module names}`. The prerequisites that are not installed on the HPC must be installed separately. Please use conda, as shown above, or the install instructions on the tools' webpage. 

### INSTALLATION

Please make sure you are in the folder where you want to install. Then download and install FastNeo using the following commands:-
```
git clone https://github.com/yashumayank/FastNeo.git
cd FastNeo
bash INSTALL
export PATH=`pwd`:$PATH
```

### TEST RUN 

Search for neoepitopes and gene fusions in the sample dataset
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

### USAGE AND OPTIONS

Command:

`run_neoepitopes.sh [optional arguments] {input prefix}`

The `optional arguments` provide an interface to customize some of the heuristics that are used to maximize the signal to noise ratio in the stranded RNAseq data. The only required parameter to run this command is the `input prefix`, which is the prefix of the paired end RNA-seq data. Two fastq files named `input prefix_1.fastq` and `input prefix_2.fastq` must be present in the working directory.

Command to detect only IEDB/TSNAdb neoepitopes and skip the gene fusions:

`runEpitopeDB_neoepitopes.sh [optional arguments] {input prefix}`

#### Optional arguments

` -m|--mapq [INTEGER] `
 
The minimum MAPQ score for the reads mapped to human consensus coding sequence can be specified here. As MAPQ score is biased against short reads. Hence, a different quality filter is used for all the reads with MAPQ below this value. If the MAPQ score is below this value then a filter specified by options `--alignscore35` and `--alignscore150` is used (default: 10)

` -x|--alignscore35 [INTEGER] `

Minumum expected alignment score if mapped read length = 35 nucleotides. This filter is used together with `--alignscore150` for the reads with MAPQ score below the value specified in `--mapq` (default: 64; allows 1 medium quality mismatch with penalty=4)

` -y|--alignscore150 [INTEGER] `

Minumum expected alignment score if mapped read length = 150 nucleotides. This filter is used together with `--alignscore35` for the reads with MAPQ score below the value specified in `--mapq` (default: 277; allows ~4 medium quality mismatch with penalty=15)

` -c|--clippedbases [NUMBER] `

This filter can be used to exclude the reads with mapped length less than the allowed fraction of the total length. Minimum allowed value of the (read length) / (clipped length) (default: 3)

` -q|--basequality [INTEGER] `

Minimum sequencing quality of all the nucleotides in the nullomer, which usually includes a few nucleotides around the mutated bases (default: 20)

` -f|--mapqf [INTEGER] `

The minimum MAPQ score for the reads mapped to gene fusion junctions can be specified here. As MAPQ score is irrelevant when mapping to a small subset of sequences, hence the value specified here is much larger than used for the `--mapq`. If the MAPQ score is below this value then a filter specified by options `--alignscore35` and `--alignscore150` is used (default: 30)

` -v|--overlap [INTEGER] `

Minimum number of nucleotides from the read that should be mapped to both sides of the gene fusion junction (default: 5)

` -o|--outprefix [STRING] `

Prefix for the output files (default: {input prefix})

### OUTPUT

The output file/files are written to the same folder from where the command was run. Output files use the `input prefix` by default unless the different `--outprefix` is specified.

The descriptions of columns in the `{prefix}_neoepitopes.tsv`:-

1) sample_id: fastq file name
2) gene_id: Gene ID with ENSG prefix
3) HGNC_symbol: Gene symbol
4) top_nullomer: nullomer with most read coverage among the nullomers that are associated to the mutation/mutations that produce these neoepitopes
5) neoepitopes: all wildtype-neoepitopes pairs associated to the top nullomer for the mutation:

   column format:
   [{wildtype epitope}->{neoepitope};] ...
7) read_count: Number of mapped reads 
8) nullomer_count: Number of nullomers on the read with most nullomers
9) db_name: Name of the database/databases that contain these neoepitopes
10) annotation: Known function of the protein that contains the epitopes
11) wildTypeHLA: binding affinities of all HLA alleles predicted to bind to each wildtype epitope

    column format:
    [Epitope:{wildtype epitope};[{HLA allele};{binding affinity},] ...] ... 
13) neoEpitopeHLA: Binding affinities of all HLA alleles predicted to bind to each neoepitope

    column format:
    [Epitope:{neoepitope};[{HLA allele};{binding affinity},] ...] ...
15) GeneticAncestry: Allelic frequency of the neoepitope in germline variants from the genetic ancestry groups

    column format:
    [SNV[1|2],{chromosome};{position};{old base};{new base},AF={AF};grpmax={group with highest frequency};AF_XX={AF_XX};AF_XY={AF_XY};AF_afr={AF_afr};AF_amr={AF_amr};AF_asj={asj};AF_eas={eas};AF_fin={fin};AF_mid={mid};AF_nfe={nfe};AF_sas={sas};] ...


The descriptions of columns in the `{prefix}_fusions.tsv`:-

1) sample_id: fastq file name
2) ChimerKB_id: ChimerKB ID
3) top_nullomer: most covered nullomer per gene fusion
4) neoepitopes: neoepitopes associated with the gene fusion
5) read_count: number of mapped reads
6) nullomer_count: number of nullomers on the read with most nullomers
7) annotation: gene names
8) junction5: genomic loci of 5’ junction
9) junction3: genomic loci of 3’ junction.

### AUXILIARY PIPELINES

Scripts in the benchmarking_scripts folder can be used to run the pipelines that were used to benchmarking this tool. These scripts are written to work on a local cluster and each of them has its specific software requirements. Please check the commands used in the pipeline and install all the tools before running them.

1) run_gatk_bcftools_lofreq.sh: Pipeline to search for IEDB and TSNEdb neoepitopes using GATK HaplotypeCaller (McKenna et al., 2010), samtools/bcftools mpileup (Li, 2011), and Lofreq (Wilm et al., 2012). The python script cds2neopitopes.py is required to run this pipeline
2) runEpitopeDB_neoepitopes_STARmapped.sh: Pipeline for neoepitopes search using FastNeo with STAR mapping instead of Bowtie2 (Langmead and Salzberg, 2012). This pipeline requires mudskipper (https://github.com/OceanGenomics/mudskipper)


### CITATION

All the pipelines and tools in this repository were released with the following manuscript :-

Mahajan M. & Hemberg M. (2024). Detecting known neoepitopes, gene fusions, transposable elements, and circular RNAs in cell-free RNA. medRxiv 2024.06.07.24308622; doi: https://doi.org/10.1101/2024.06.07.24308622
