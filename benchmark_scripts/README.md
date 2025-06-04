### NEOEPITOPE DETECTION: BENCHMARKING PIPELINES

The pipeline in `run_gatk_bcftools_lofreq.sh` was used to benchmark FastNeo against some of the most famous variant calling methods. This pipeline detects the known neoepitopes that are characterised in IEDB (https://www.iedb.org) and TSNAdb (http://biopharm.zju.edu.cn/) using GATK HaplotypeCaller (McKenna et al., 2010), samtools/bcftools mpileup (Li, 2011), and Lofreq (Wilm et al., 2012). The python script `cds2neopitopes.py` is required to run this pipeline. 

The pipeline in `runEpitopeDB_neoepitopes_STARmapped.sh` was used to detect neoepitopes using FastNeo with STAR mapping instead of Bowtie2 (Langmead and Salzberg, 2012). This pipeline requires mudskipper (https://github.com/OceanGenomics/mudskipper) to convert Genome mapped co-ordinates to CDS mapped co-ordinates .

Both these pipelines are written to work on a local cluster and each of them has its specific software requirements. Please check the commands used in the pipeline and install the prerequisites before running.


