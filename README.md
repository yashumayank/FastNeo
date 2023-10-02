# RNA-neoepitopes
This tool has been developed to detect neoepitopes from published resources in the cell free RNA. It supports detection of neoepitopes from IEDB, TSNAdb, neoepitopes produced via the gene fusion events described in ChimerKB and ChimerPub (YE Jang et al. 2020: https://doi.org/10.1093/nar/gkz1013), and neoepitopes produced transcription of Transposable elements (NM Shah et al. 2023: https://doi.org/10.1038/s41588-023-01349-3). Cancer specific neoepitopes: TSNAdb and Transposable elements

# REFERENCES/AUTHOR

# SUPPORT/REQUESTS

# HARDWARE/SOFTWARE REQUIREMENTS
x86-64 compatible processors
64 bit Linux or Mac OS X

# PRE-REQUISITES
- For Mac OSx:
Install Xcode or to install only git: "brew install git"
- For Windows:
git for windows (https://gitforwindows.org/),
- Common requirements fro linux, windows and OSx
Julia (modules required: --- ), Python (Biopython), Bowtie2, samtools, bedtools

RELEASEnotes 
contains detailed information about the latest major release CHANGES contains detailed information about all the changes in all releases

# INSTALL AND SETUP
git clone https://github.com/yashumayank/cfRNA-neoepitopes.git
Download mapping indices
-   for fusions (https://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/)
-   proteins (https://ftp.ensembl.org/pub/release-107/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz)

# LIMITATIONS
The neoepitopes found in plasma-cell-free RNA represents cells from the whole body and only a miniscule proportion of the RNA is from the tumour. Thus, even the most discriminative neoepitopes might not have originated in the tumour itself.

# FUNDING

