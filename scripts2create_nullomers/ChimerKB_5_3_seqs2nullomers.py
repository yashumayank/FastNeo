#!/usr/bin/python
# encoding: utf-8
#Sample input header: >ENST00000390237:1:SRR10822565:HCC:104:m:35
#Produces a single file with peptides of various lengths The header consists of transcript_ID, peptide_ID, #Cancer/#control samples
#peptide IDs generated for each transcript as follows: first,last,mutated peptide wrt. original inpuut peptide
#Usage:
#python3 ChimerKB_5_3_seqs2nullomers.py <Homo_sapiens_cds_16mers> <Homo_sapiens_GRCh38.pep.all.fa.gz> <output_prefix>
#Example:
#python3 ChimerKB_5_3_seqs2nullomers.py Homo_sapiens_cds_16mers.tab Homo_sapiens.GRCh38.pep.all.fa.gz ChimerKB_n_Pub

'''
@author:     Mayank Mahajan
@copyright:  2022 Brigham and women's hospital. All rights reserved.
@license:    GPL3
@contact:    yashumayank@gmail.com
@deffield    updated: Updated
'''

fusPep=set()
pepKmer=set()
nullPep_len = 9
null_len = 0
print("Reading cds k-mers")
with open(sys.argv[1]) as cdsKmer:
    for lines in cdsKmer:
        Chimertok = lines.split()
        fusPep.add(Chimertok[0])
        null_len = len(Chimertok[0])

print("Reading protein k-mers")
with gzip.open(sys.argv[2], "rt") as handle:
    for recProt in SeqIO.parse(handle, "fasta"):
        protSeq=recProt.seq
        if len(protSeq) >=9:
            slen=len(protSeq)-nullPep_len+1
            for kstt in range(slen):
                kstp=kstt+nullPep_len
                if not "*" in protSeq[kstt:kstp]:
                    pepKmer.add(protSeq[kstt:kstp])

flen=null_len-1
nf = open(sys.argv[3] + "_junction_nullomers.tsv", "w")
jf = open(sys.argv[3] + "_15bp_around_junction.fasta", "w")
mfasta = open(sys.argv[3] + "_for_mapping.fasta", "w")
print("Reading sequence around junctions, extracting junction-cds_nullomers and mapable fasta")
for rec5, rec3 in zip(SeqIO.parse(sys.argv[3] + "_5primeJunction.fasta", "fasta"),
                      SeqIO.parse(sys.argv[3] + "_3primeJunction.fasta", "fasta")):

    juncFrame5 = rec5.description.split()[1].split(':')[3]
    juncFrame3 = rec3.description.split()[1].split(':')[3]
    neopepStr=""
    if not juncFrame5 == "x":
        plen5=int(juncFrame5)
        plen3=int(juncFrame3)
        junc_pseq = Seq(rec5.seq[-plen5:] + rec3.seq[:plen3])
        fusion_pep= junc_pseq.translate(to_stop=True)
        if len(fusion_pep) >= 9:
#            nf.write(str(rec5.seq[-plen5:] + " " + rec3.seq[:plen3]) + "\n")
#            nf.write(str(fusion_pep) + "\n")
            plen = len(fusion_pep)-nullPep_len + 1
            for pstt in range(plen):
                pstp=pstt+nullPep_len
                if not(fusion_pep[pstt:pstp] in pepKmer):
                    neopepStr=neopepStr + ";" + str(fusion_pep[pstt:pstp])
    #nf.write("\t"+ str(rec5.description) +" <-> "+ str(rec3.description)+"\n")

    junc_seq = (rec5.seq[-flen:] + rec3.seq[:flen]).upper()
    for stt in range(flen):
        stp=stt+null_len
        if not(junc_seq[stt:stp] in fusPep):
             nf.write(str(junc_seq[stt:stp]) +";")
    nf.write("\t" + juncFrame5 + neopepStr + "\t"+ str(rec5.description) +"\t"+ str(rec3.description.split()[1])+"\n")

    null_record = SeqRecord(Seq(junc_seq), id=rec5.id,
                       description = rec5.description + "\t" + rec3.description.split()[1])
    SeqIO.write(null_record, jf, "fasta")
    map_record = SeqRecord(Seq(rec5.seq + rec3.seq), id=rec5.id,
                       description = rec5.description + "\t" + rec3.description.split()[1])
    SeqIO.write(map_record, mfasta, "fasta")

nf.close()
jf.close()
mfasta.close()
