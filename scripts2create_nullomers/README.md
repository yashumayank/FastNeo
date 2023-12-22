# neoepitope-to-nullomers
Workflow to produce nullomers corresponding to neoepitopes from various sources. These scripts also create the alignment file for gene-fusion reads and the maps to convert nullomer to neoepitope

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
- protein ID to transcript id mapiing files from uniprot and gencode HUMAN_9606_idmapping_selected.tab, gencode.v40.metadata.TrEMBL, gencode.v40.metadata.SwissProt  (downloaded on June-27-2022)

## Create nullomer lists, mapping indices and associated files

The following two shell scripts contain the workflow to create nullomers and mapping files that are required to run this pipeline.

```
createEpitopeDB-nullomers.sh
createFusion-nullomers.sh
```
