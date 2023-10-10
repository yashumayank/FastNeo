#!/usr/bin/python
# encoding: utf-8
#Sample input header: >ENST00000390237:1:SRR10822565:HCC:104:m:35
#Produces a fasta file each for 1. bases upstream of 5' junction, and 2. bases downstream of 3' junction for each row in ChimerDB gene fusion file.
#
#Usage:
#python3 extractChimerDbNullomers.py <ChimerKB_n_Pub.tab> <hg19.fa.gz> <hg19.ensGene.gtf> <output_prefix>

'''
@author:     Mayank Mahajan

@copyright:  2022 Brigham and women's hospital. All rights reserved.

@license:    GPL3

@contact:    yashumayank@gmail.com
@deffield    updated: Updated

'''

from Bio import SeqIO
from Bio.Seq import MutableSeq
import sys, csv, re
import gzip

#from optparse import OptionParser

__all__ = []
__version__ = 0.1
__date__ = '2022-07-06'
__updated__ = '2023-08-31'

strand5={}
strand3={}
juncName={}
juncChr5={}
juncChr3={}
juncPos5={}
juncPos3={}
juncSeq5={}
juncSeq3={}
juncHead5={}
juncHead3={}
juncFrame5={}
juncFrame3={}
startcChr={}
cdsChr={}
startcStart={}
cdsStart={}
cdsStop={}
cdsStrand={}
cdsFrame={}
xCount=0
nullPep_len=9
includeRes3end=20
#fick = sys.stdout
#Read the geneIDs, neopeptides/gene & their lengths from IEDB
#Read the geneIDs, neopeptides/gene & their lengths from IEDB
print("Reading cds frame and strand from hg19 gtf")
cdsCount=0
startcCount=0
with open(sys.argv[3]) as hg19_cds:
    for ensts in hg19_cds:
        ensttok = ensts.split("\t")
        if(ensttok[2]=="CDS"):
            cdsChr[cdsCount]=ensttok[0][3:]
            cdsStart[cdsCount]=int(ensttok[3])
            cdsStop[cdsCount]=int(ensttok[4])
            cdsStrand[cdsCount]=ensttok[6]
            cdsFrame[cdsCount]=int(ensttok[7])
#            print(cdsChr[cdsCount] + " " + cdsStart[cdsCount] + " " + cdsStop[cdsCount] + " " + cdsStrand[cdsCount] + " "+ cdsFrame[cdsCount])
            cdsCount+=1
        elif(ensttok[2]=="start_codon"):
            startcChr[cdsCount]=ensttok[0][3:]
            startcStart[cdsCount]=int(ensttok[3])
            startcCount+=1

print("Reading fusion junctions from ChimeraKB_Pub database")
with open(sys.argv[1]) as Chimerf:
    next(Chimerf)
    for lines in Chimerf:
        Chimertok = lines.split("\t")
        if(len(Chimertok)<15):
            continue
#        valid = re.compile(r"^[Cc]hr\d+:\d.*")
#        valid.match(Chimertok[5])
#        valid.match(Chimertok[6])
        junc5 = re.split(r"[r :]",Chimertok[5])
        junc3 = re.split(r"[r :]",Chimertok[6])
        if(len(junc5)<3 or len(junc3)<3):
            continue
        try:
            if(int(junc5[1]) < 1 or int(junc5[2]) < 1 or int(junc3[1]) < 1 or int(junc3[2]) < 1 or (Chimertok[10]!="+" andChimertok[10]!="-") or (Chimertok[14]!="+" and Chimertok[14]!="-")):
                raise ValueError
            juncChr5[Chimertok[0]]= junc5[1]
            juncChr3[Chimertok[0]]= junc3[1]
            juncPos5[Chimertok[0]]= int(junc5[2])
            juncPos3[Chimertok[0]]= int(junc3[2])
            strand5[Chimertok[0]] = Chimertok[10]
            strand3[Chimertok[0]] = Chimertok[14]
            juncName[Chimertok[0]]= Chimertok[4]
        except ValueError:
            continue
#        print(f'{juncName[Chimertok[0]]}\t{juncChr5[Chimertok[0]]}\t{juncPos5[Chimertok[0]]}\t{strand5[Chimertok[0]]}\t{juncChr3[Chimertok[0]]}\t{juncPos3[Chimertok[0]]}\t{strand3[Chimertok[0]]}')

        dist5=""
        frame5=""
        for cdsi in cdsChr:           #search 5' cds
            if(cdsChr[cdsi] == junc5[1]):
                if(cdsStart[cdsi] <= int(junc5[2]) and cdsStop[cdsi] >= int(junc5[2]) and Chimertok[10] == cdsStrand[cdsi]):
#                    print(Chimertok[0] + " " + Chimertok[10] + " " + cdsStrand[cdsi])
#---if junction inside the cds: distance from exon start to 5' junction and then substract frame in gtf
                    if(Chimertok[10]=="+"):
                        dist5 = int(junc5[2]) - cdsStart[cdsi] +1
                    elif(Chimertok[10]=="-"):
                        dist5 = cdsStop[cdsi] - int(junc5[2]) +1
                    frame5=cdsFrame[cdsi]
                    continue
        if(not(dist5=="")):
            dist3=""
            left3=""
            frame3=""
            for cdsi in cdsChr:       #search 3' cds
                if(cdsChr[cdsi] == junc3[1]):
                    if(cdsStart[cdsi] <= int(junc3[2]) and cdsStop[cdsi] >= int(junc3[2]) and Chimertok[14] == cdsStrand[cdsi]):
                        if(Chimertok[14]=="+"):
                            dist3 = int(junc3[2]) - cdsStart[cdsi]
                            left3 = cdsStop[cdsi] - int(junc3[2]) +1
                        elif(Chimertok[14]=="-"):
                            dist3 = cdsStop[cdsi] - int(junc3[2])
                            left3 = int(junc3[2]) - cdsStart[cdsi] +1
                        frame3=cdsFrame[cdsi]
                        continue
            if(not(dist3=="")):        #check start codon is intact
                for startci in startcChr:
                    if(startcChr[startci] == junc5[1]):
                        tmpsDist = int(junc5[2]) - int(startcStart[startci])
                        if(tmpsDist>=0 and tmpsDist<=2):
                            dist5==""
                            continue
        if(dist5=="" or dist3==""):
            juncFrame3[Chimertok[0]]="x"
            juncFrame5[Chimertok[0]]="x"
            xCount+=1
        else:
#-------- if original 3' frames == new 3' frame  --------
            frame3j = (3-((dist5 -frame5) %3)) %3
            frame3j_original = (3-((dist3 -frame3) %3)) %3
            if(frame3j == frame3j_original):
                juncFrame5[Chimertok[0]] =str((nullPep_len *3) +(3-frame3j) +3)
                juncFrame3[Chimertok[0]] =str((nullPep_len *3) + frame3j +3)
                if(int(juncFrame5[Chimertok[0]])>dist5):
                    juncFrame5[Chimertok[0]] = str(dist5-frame5)
                if(int(juncFrame3[Chimertok[0]])>left3):
                    juncFrame3[Chimertok[0]] = str(left3)
                if(int(juncFrame5[Chimertok[0]])<1 or int(juncFrame3[Chimertok[0]])<1): #remove edge junctions
                    juncFrame5[Chimertok[0]] = "x"
                    juncFrame3[Chimertok[0]] = "x"
                    xCount+=1
            else:
                juncFrame5[Chimertok[0]] = "x"
                juncFrame3[Chimertok[0]] = "x"
                xCount+=1

#           print(Chimertok[0] + " " + junc5[2] + " " + Chimertok[10] + " " + cdsStrand[cdsi])
print(" - Fusion junctions ouside CDS: " + str(xCount))

print("Reading genome sequence and extracting 500bp near the junction")

with gzip.open(sys.argv[2], "rt") as handle:
    for seq_record in SeqIO.parse(handle, "fasta"):     #read one fasta record at a time
        token=str(seq_record.id).split('hr')            #split header into tokens

        for jid in juncName:
            if(juncChr5[jid]==token[1]):
                juncHead5[jid] = ">" + jid + "\tchr" + juncChr5[jid] + ":" + str(juncPos5[jid]) + ":" + strand5[jid] + ":"+ juncFrame5[jid] + "\t" + juncName[jid]
                if strand5[jid]=="+":                   #5prime + strand
                    stt = juncPos5[jid]-500
                    stp = juncPos5[jid]
                    juncSeq5[jid] = str(seq_record.seq[stt:stp])
                else:                                   #5prime - strand
                    stt = juncPos5[jid]-1
                    stp = juncPos5[jid]+499
                    juncSeq5[jid] = str(seq_record.seq[stt:stp].reverse_complement())
            if(juncChr3[jid]==token[1]):
                juncHead3[jid] = ">" + jid + "\tchr" + juncChr3[jid] + ":" + str(juncPos3[jid]) + ":" + strand3[jid] + ":"+ juncFrame3[jid] + "\t" + juncName[jid]
                if strand3[jid]=="+":                   #3prime + strand
                    stt = juncPos3[jid]-1
                    stp = juncPos3[jid]+499
                    juncSeq3[jid]=str(seq_record.seq[stt:stp])
                else:                                   #3prime - strand
                    stt = juncPos3[jid]-500
                    stp = juncPos3[jid]
                    juncSeq3[jid]=str(seq_record.seq[stt:stp].reverse_complement())

fasta5 = open(sys.argv[4] + "_5primeJunction.fasta", "w")
fasta3 = open(sys.argv[4] + "_3primeJunction.fasta", "w")
for jid in juncName:
    fasta5.write(juncHead5[jid] +"\n")
    fasta5.write(juncSeq5[jid] +"\n")
    fasta3.write(juncHead3[jid] +"\n")
    fasta3.write(juncSeq3[jid] +"\n")
fasta5.close()
fasta3.close()
