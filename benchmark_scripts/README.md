### NEOEPITOPE DETECTION: BENCHMARKING PIPELINE 

These scripts were used to benchmarking FastNeo against the most famous variant calling tools. They detect the known neoepitopes that are characterised in IEDB (https://www.iedb.org) and TSNAdb (http://biopharm.zju.edu.cn/). These scripts are written to work on a local cluster and each of them has its specific software requirements. Please check the commands used in the pipeline and install all the tools before running them.

1) run_gatk_bcftools_lofreq.sh: Pipeline to search for IEDB and TSNEdb neoepitopes using GATK HaplotypeCaller (McKenna et al., 2010), samtools/bcftools mpileup (Li, 2011), and Lofreq (Wilm et al., 2012). The python script cds2neopitopes.py is required to run this pipeline
2) runEpitopeDB_neoepitopes_STARmapped.sh: Pipeline for neoepitopes search using FastNeo with STAR mapping instead of Bowtie2 (Langmead and Salzberg, 2012). This pipeline requires mudskipper (https://github.com/OceanGenomics/mudskipper)

