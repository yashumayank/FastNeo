#Usage:
#python createGenomeKmers.py <genome_gzip_file> <kmer_size>

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
import sys, csv, re
import gzip


genKmer={}
kmer_len = int(sys.argv[2])
print("Reading genome k-mers")
with gzip.open(sys.argv[1], "rt") as handle:
    for recScaf in SeqIO.parse(handle, "fasta"):
        scafSeq=recScaf.seq
        scafSeqRev=recScaf.seq.reverse_complement()
        if len(scafSeq) >= kmer_len:
            slen=len(scafSeq)-kmer_len+1
            for kstt in range(slen):
                kstp=kstt+kmer_len
                if not "*" in scafSeq[kstt:kstp]:
                    if(scafSeq[kstt:kstp] in genKmer):
                        genKmer[scafSeq[kstt:kstp]]+=1
                    else:
                        genKmer[scafSeq[kstt:kstp]]=1
                if not "*" in scafSeqRev[kstt:kstp]:
                    if(scafSeqRev[kstt:kstp] in genKmer):
                        genKmer[scafSeqRev[kstt:kstp]]+=1
                    else:
                        genKmer[scafSeqRev[kstt:kstp]]=1

for i in genKmer:
    print(i + "\t" + genKmer[i])
