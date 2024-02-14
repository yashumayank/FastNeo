# RNA-neoepitopes
This tool has been developed to detect neoepitopes from published resources in the cell free RNA. It supports detection of neoepitopes from IEDB, TSNAdb, and neoepitopes produced via the gene fusion events described in ChimerKB and ChimerPub (YE Jang et al. 2020: https://doi.org/10.1093/nar/gkz1013). Cancer specific neoepitopes: TSNAdb

# REFERENCES/AUTHOR

# SUPPORT/REQUESTS

# HARDWARE/SOFTWARE REQUIREMENTS
x86-64 compatible processors
64 bit Linux or Mac OS X

# RELEASEnotes 
contains detailed information about the latest major release CHANGES contains detailed information about all the changes in all releases

# INSTALL AND SETUP (PRE-REQUISITES)

Download the tool
```
git clone https://github.com/yashumayank/cfRNA-neoepitopes.git
```

Download and install bowtie2, samtools and julia, which are required to run this tool. The following command can be used to install through conda :-
```
conda install -c bioconda bowtie2 samtools julia
```

Use the following code to install the required packages in Julia.
```
julia pkginstall.jl
```

# RUN the pipelines 

The paired-end RNAseq data must be in the 2-file fastq format and the files must be named as {sample_name}_1.fastq and {sample_name}_2.fastq

The IEDB/TSNAdb neoepitope detection can be run together using the following command :- 
```
sh runEpitopeDB_neoepitopes.sh {sample_name}
```

The gene fusion detection can be run together using the following command:- 
```
sh runFusion_neoepitopes.sh {sample_name}
```
Both IEDB/TSNAdb neoepitope and gene fusion detection can be run together using the following command
```
sh run_neoepitopes.sh {sample_name}
```


# OUTPUT format

For IEDB/TSNEdb neoepitopes, the output file shows:

1) sample_id: fastq file name
2) ENST_id: transcript ensembl ID 
3) top_nullomer: nullomer with most read covergage among the nullomers that are associated to the mutation/mutations that produce these neoepitopes
4) neoepitopes: all wildtype-neoepitopes pairs associated to the top nullomer for the mutation:
   
   [{Wildtype epitope}->{Noepitope};] ...
5) #reads: number of mapped reads 
6) #nullomers: number of nullomers on the read with most nullomers
7) db_name: Name of the database/databases that contain these neoepitopes.
8) annotation: known function of the protein that contains the epitopes
9) wildTypeHLA: binding affinities of all HLA alleles predicted to bind to each wildtype epitope. column output format:
 
    [Epitope:{wildtype epitope};[{HLA allele};{binding affinity},] ...] ... 
10) neoEpitopeHLA: binding affinities of all HLA alleles predicted to bind to each neoepitope column output format:

    [Epitope:{neoepitope};[{HLA allele};{binding affinity},] ...] ... 
11) GeneticAncestry: frequency of the neoepitope producing mutations in healthy individuals of various genetic ancestry groups. column format:
    
    [SNV[1|2],{chromosome};{position};{old base};{new base},AF={AF};grpmax={group with highest frequency};AF_XX={AF_XX};AF_XY={AF_XY};AF_afr={AF_afr};AF_amr={AF_amr};AF_asj={asj};AF_eas={eas};AF_fin={fin};AF_mid={mid};AF_nfe={nfe};AF_sas={sas};] ...


Output file format for the gene fusions:

1) sample_id: fastq file name
2) ChimerKB_id: ChimerKB ID
3) top_nullomer: most covered nullomer per gene fusion
4) neoepitopes: neoepitopes associated with the gene fusion
5) #reads: number of mapped reads
6) #nullomers: number of nullomers on the read with most nullomers
7) annotation: gene names
8) junction5: genomic loci of 5’ junction
9) junction3: genomic loci of 3’ junction.
10) (to be added) bases_mapped: number of bases mapped upstream and downstream of the fusion junction

# LIMITATIONS
The neoepitopes found in plasma-cell-free RNA represents cells from the whole body and only a miniscule proportion of the RNA is from the tumour. Thus, even the most discriminative neoepitopes might not have originated in the tumour itself.

# FUNDING

