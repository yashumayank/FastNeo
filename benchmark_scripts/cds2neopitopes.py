#!/usr/bin/python
# encoding: utf-8
#Usage:
#python3 cds2peptides.py <IEDB_neoepitopes.tab> <TSNAdb_TCGA_neoepitopes.tab> <TSNAdb_TCGA_neoepitopes.tab> <codon_sequences_with mutations.fasta>  > <neoepitopes_found.tsv>
#<codon_sequences_with mutations.fasta> is a single file with neoepitopes from IEDB and TSNAdb that appear in a list of mutated codon sequences.
#Sample header for input fasta: >ENST00000390237:1:ENSG00000110338:104:m:35
#<neoepitope_found.tsv> is a list of: neoepitope     wildtype_epitope     transcript_ID     gene_ID     mutation_ID     epitope_DBs
#transcript_ID is ENSTxx,position_of_first_mutation,[m_WTbase|f] m and f symbolize missene mutation and frameshift mutation respectively 
#mutation_ID is generated for each transcript as follows: first,last,mutated peptide wrt. original input

'''
@author:     Mayank Mahajan

@copyright:  2022 Brigham and women's hospital. All rights reserved.

@license:    GPL3

@contact:    yashumayank@gmail.com
@deffield    updated: Updated

'''

from Bio import SeqIO
from Bio.Seq import MutableSeq
import sys, csv
#import os
#from optparse import OptionParser

__all__ = []
__version__ = 0.1
__date__ = '2022-07-06'
__updated__ = '2024-02-06'

trscptID={}
pepID={}
pepDB={}
IEDBpep={}
TICGCpep={}
TTCGApep={}
IEDBlen={}
TICGClen={}
TTCGAlen={}
WTpeps={}
gid={}
#fick = sys.stdout
# 2. use ENST_pep as id to account for same neoepitopes in multiple genes
def scanEpi(pList,pLens,dbName,pepI,mutR,pepWT, mFlag, mPos):
     # print(pList)
     # print(pLens)
     # create peptide of lengths that occur in the list of IEDB peptides for this gene
    for i in sorted(set(pLens)):                     #loop to get different petide lengths 
        start=0                                      #first res to start the loop
        if mutR > i:                                 #check if mutated res. is farther that size of output peptide 
            start=mutR-i                             #if yes shift the start farther to inculde the mutation
        end=len(pepI)-i                              #last res to finish the loop
        if end < start:                              #check output peptide is longer that input 
            break
        if mFlag != "f" and end >= mutR:             #if not a frameshift_variant check that the last peptide start is <= the mutated peptide
            end=mutR-1
        for j in range(start, end+1):
            pep = pepI[j:j+i]
            #peptideID is count within trscpt, first, last and muitated peptide wrt. input peptide
            tmpPepID = mPos + "," + str(j+1) + "," + str(j+i)  + "," + str(mutR)

        #Check if the sample neoepitope is among the IEDB neoepitopes
            if(pep in pList):
                if pep in trscptID:                     #Is pep already in dict
                    if not (tmpTrID in trscptID[pep]):  #exclude if duplicate occurances
                        trscptID[pep] = trscptID[pep] + ":" + tmpTrID
                        pepID[pep] = pepID[pep] + ":" + tmpPepID
                    if not (dbName in pepDB[pep]):
                        pepDB[pep] = pepDB[pep] + ":" + dbName
                else:
                    trscptID[pep] = tmpTrID
                    pepID[pep] = tmpPepID
                    pepDB[pep] = dbName
                    if mFlag != "f":
                        WTpeps[pep]=pepWT[j:j+i]        #create an dictionary of WT peptides
                    else:
                        WTpeps[pep]=pepWT

#Read the geneIDs, neopeptides/gene & their lengths from IEDB
print("Reading IEDB neoepitopes")
with open(sys.argv[1]) as IEDBf:
    for lines in IEDBf:
        IEDBtok = lines.split()
#        print(IEDBtok[0] + " " + IEDBtok[1] + " " + IEDBtok[2])
        IEDBpep[IEDBtok[0]] = IEDBtok[1]
        IEDBlen[IEDBtok[0]] = IEDBtok[2]
print("Reading TSNAdbTCGA neoepitopes")
with open(sys.argv[2]) as TICGCf:
    for lines in TICGCf:
        itok = lines.split()
#        print(IEDBtok[0] + " " + IEDBtok[1] + " " + IEDBtok[2])
        TICGCpep[itok[0]] = itok[1]
        TICGClen[itok[0]] = itok[2]
print("Reading TSNAdbICGC neoepitopes")
with open(sys.argv[3]) as TTCGAf:
    for lines in TTCGAf:
        ttok = lines.split()
#        print(IEDBtok[0] + " " + IEDBtok[1] + " " + IEDBtok[2])
        TTCGApep[ttok[0]] = ttok[1]
        TTCGAlen[ttok[0]] = ttok[2]

print("Reading neoepitopes in samples and comparing to IEDB ")
for seq_record in SeqIO.parse(sys.argv[4], "fasta"):     #read one fasta record at a time
    token=str(seq_record.id).split(':')                  #split header into tokens
    #  print(token[0] + "," + token[1] + "," + token[2] + "," + token[3] + "," + token[4] + "," + token[5])
    pepI=str(seq_record.seq.translate(to_stop=True))     #get the full input peptide
    mutR=int((int(token[5])+2)/3)                        #get the position of mutated residue 
    gid[token[0]]=token[1]
    if token[4] != "f":
        seq_WT=MutableSeq(seq_record.seq)
        try:
            seq_WT[int(token[5])-1] = token[4][2]        #create the WT cds with token[5][2] at token[6]
        except:
            print("Token4: " + token[4])
        pepWT=str(seq_WT.translate(to_stop=True))        #get the full WT peptide
    else:
        pepWT="-"
     
    #Transcript ID, start, type_of_variant for each input peptide cds
    tmpTrID = token[0] + "," + token[3] + "," + token[4]
    ##Check if this gene is in IEDB and length of the neo-peptides in this gene
   # print(token[0])
    if token[0] in IEDBpep:
        pepList=IEDBpep[token[0]].split(';')
        pepLens=list(map(int,IEDBlen[token[0]].split(';')))
        scanEpi(pepList, pepLens, "IEDB", pepI, mutR, pepWT, token[4], token[2])
    if token[0] in TICGCpep:
        pepList=TICGCpep[token[0]].split(';')
        pepLens=list(map(int,TICGClen[token[0]].split(';')))
        scanEpi(pepList, pepLens, "ICGC", pepI, mutR, pepWT, token[4], token[2])
    if token[0] in TTCGApep:
        pepList=TTCGApep[token[0]].split(';')
        pepLens=list(map(int,TTCGAlen[token[0]].split(';')))
        scanEpi(pepList, pepLens, "TCGA", pepI, mutR, pepWT, token[4], token[2])

neoF = open(sys.argv[5] + "_neoepitopes.tsv", "w")
for p in pepID:
  #The header consists of transcript_ID - peptide_ID - #Cancer/#control samples eg. "ENST00000362079,700,m-618,42,51,50-4/1"
  #transcript_ID: EnsemblID, mutation start position, type_of_variant for each input peptide cds
  #Peptide_ID: mutated residue in protein, first, last and mutated residue in the full peptide
  tks = trscptID[p].split(",")
  try:
      neoF.write(p + "\t" + WTpeps[p] + "\t" + trscptID[p]  + "\t" + gid[tks[0]] + "\t" + pepID[p] + "\t" + pepDB[p] +"\n")
  except:
      print("Failed:\t" + trscptID[p]  + "\t" + tks[0])
neoF.close()
