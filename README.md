# RNA-neoepitopes
This tool has been developed to detect neoepitopes from published resources in the cell free RNA. It supports detection of neoepitopes from IEDB, TSNAdb, neoepitopes produced via the gene fusion events described in ChimerKB and ChimerPub (YE Jang et al. 2020: https://doi.org/10.1093/nar/gkz1013), and neoepitopes produced transcription of Transposable elements (NM Shah et al. 2023: https://doi.org/10.1038/s41588-023-01349-3). Cancer specific neoepitopes: TSNAdb and Transposable elements

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

Open 'runAll-neoepitopes.sh' and 'runEpitopeDB_SNV_neoepitopes.sh' in a text editor and update path1 variable with the full path of the 'cfRNA-neoepitopes' directory  


Download mapping indices


# LIMITATIONS
The neoepitopes found in plasma-cell-free RNA represents cells from the whole body and only a miniscule proportion of the RNA is from the tumour. Thus, even the most discriminative neoepitopes might not have originated in the tumour itself.

# FUNDING

