# neoepitope-nullomers
Workflow to produce nullomers corresponding to neoepitopes from various sources. These scripts also create the alignment file for gene-fusion reads and the nullomer to neoepitope map files

## Prerequisites
- Python (Biopython, numpy)
- R (dplyr, data.table, stringr)
- bedtools
- samtools
- genome fasta and cds fasta, and protein fasta (https://ftp.ensembl.org/pub/release-107/fasta/homo_sapiens/)
- TSNAdb_frequent_neoantigen_TCGA_4.0_adj.txt, TSNAdb_frequent_neoantigen_ICGC_4.0_adj.txt from TSNAdb (http://biopharm.zju.edu.cn/tsnadb/download/)
- epitope_full_v3.csv, mhc_ligand_full.csv, tcell_full_v3.csv, bcell_full_v3.csv from IEDB database (downloaded on Sept-19-2022)
- ChimerKB4.xlsx from ChimerDB. Keep rows common with 'Pub' and export as a tsv file to ChimerKB_n_Pub.tab  (downloaded on May-19-2023)
- Genome fasta and gtf file from hg19 (https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/)
- TE_neoantigens.tsv and TE_chimeric_2297.tsv from Shah et al. 2023 repository 
- protein ID to transcript id mapiing files from uniprot and gencode HUMAN_9606_idmapping_selected.tab, gencode.v40.metadata.TrEMBL, gencode.v40.metadata.SwissProt  (downloaded on June-27-2022)

## Create the files required to run cfRNA-neoepitopes pipeline
 
```
createEpitopeDB-nullomers.sh
createFusion-nullomers.sh
createTE-nullomers.sh
```
