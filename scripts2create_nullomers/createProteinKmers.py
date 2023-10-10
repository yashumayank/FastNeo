from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
import sys, csv, re
import gzip


pepKmer=set()
nullPep_len = 9
print("Reading protein k-mers")
with gzip.open(sys.argv[1], "rt") as handle:
    for recProt in SeqIO.parse(handle, "fasta"):
        protSeq=recProt.seq
        if len(protSeq) >= 9:
            slen=len(protSeq)-nullPep_len+1
            for kstt in range(slen):
                kstp=kstt+nullPep_len
                if not "*" in protSeq[kstt:kstp]:
                    pepKmer.add(protSeq[kstt:kstp])

for i in pepKmer:
    print(i)
