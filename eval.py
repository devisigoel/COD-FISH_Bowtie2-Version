#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys, re
from Bio.SeqUtils import MeltingTemp
import primer3
import glob
from Bio import pairwise2

ValueFile = "Values.csv"

def loadENSTmap(filename):
    EtoN = {}
    NtoE = {}
    for ii in open(filename):
        ii = ii.strip()
        iiL = ii.split(" ")
        EtoN[iiL[2]] = iiL[0]
        NtoE[iiL[0]] = iiL[2]

    return(EtoN, NtoE)

def complement(seq):
    compl = {'A':'T','T':'A','C':'G','G':'C'}
    compList = []
    for nt in seq:
        compList.append(compl[nt])
    compSeq = ''.join(compList)
    return(compSeq)


def load_TranscriptomeDict(transcriptomeFile):
    transcriptomeDict = {}
    GeneSymbolDict = {}
    for transcriptomemRNA in open(transcriptomeFile):
        transcriptomemRNA = transcriptomemRNA.strip()
        if transcriptomemRNA[0] == ">":
            seqname = transcriptomemRNA.replace(">", "")
            seqNL = seqname.split(" ") 
            ENSTid = re.sub("\..*", "", seqNL[0])
            if len(seqNL) >= 4:
                Gene = seqNL[3].split(":")
                if len(Gene) == 2:
                    GeneSymbolDict[ENSTid] = Gene[1]
                else:
                    GeneSymbolDict[ENSTid] = ENSTid
            transcriptomeDict[ENSTid] = ""
        else:
            transcriptomeDict[ENSTid] = transcriptomeDict[ENSTid] + transcriptomemRNA
    return(transcriptomeDict, GeneSymbolDict)

#def alignmentScore(mRNA_seq, probe_seq, ProbeTm):
#    nt_Dict = {('A','T'): 1, ('A','A'): -1, ('A','C'): -1, ('A','G'): -1, 
#           ('C','A'): -1, ('C','T'): -1, ('C','C'): -1, ('C','G'): 2,
#           ('G','A'): -1, ('G','C'): 2, ('G','G'): -1, ('G','T'): -1,
#           ('T','A'): 1, ('T','C'): -1, ('T','T'): -1, ('T','G'): -1}
#    mRNA_Align = pairwise2.align.globalds(probe_seq, mRNA_seq, nt_Dict, -1, -1, penalize_end_gaps=(False, False), score_only = True)  
#    PerfectAlign = int(pairwise2.align.globalds(probe_seq, complement(probe_seq), nt_Dict, -1, -1, penalize_end_gaps=(False, False), score_only = True))    
#    score = (mRNA_Align/PerfectAlign) * ProbeTm
#    return(score)

#inPrefix= sys.argv[1]
BLASTfile = sys.argv[1]
TranscriptomeFile = sys.argv[2]
probeFastaFile = sys.argv[3]
#startIDX = int(sys.argv[4])
mRNAGene = sys.argv[4]
targetmRNA = sys.argv[5]
mapFile = sys.argv[6]

transcriptomeDict, GeneSymbolDict = load_TranscriptomeDict(TranscriptomeFile)

def load_Probes(probeFastaFile):
    probeDict = {}
    probe, pName = '', ''
    for line in open(probeFastaFile):
        line = line.strip()
        if line.find(">") == 0: 
            if probe != '':
                probeDict[pName] = probe.upper()
            probe = ''
            pName = line.replace(">", "")
        else:
            probe = probe + line

    if probe != '':
        probeDict[pName] = probe.upper()
    
    return(probeDict)

EtoN, NtoE = loadENSTmap(mapFile)
probeDict = load_Probes(probeFastaFile)
ScoreList = []
numProbes = len(probeDict)
#for aa in range(numProbes):
#    infoDict = {}
#    fileIDX = startIDX + aa
#    filePath = inPrefix + "-" + str(fileIDX)
#    print("Processing ", filePath)

allProbes = {}
for ll in open(BLASTfile):
    ll = ll.strip()
    
    if ll.find('probe') == 7:
        lList = ll.split(" ")
        probeName = lList[1]
        allProbes[lList[1]] = {}
    if ll.find('>') == 0:
        name = ll
        name = name.split(' ')
        name = name[0][1:]
        allProbes[probeName][name] = []
    elif ll.find('Query ') == 0:
        qList = []
        x = re.split('\s+',ll)
        allProbes[probeName][name].append((x[1],x[3]))
    elif ll.find('Sbjct') == 0:
        sList = []
        y = re.split('\s+',ll)
        allProbes[probeName][name].append((y[1], y[3]))

for aa in allProbes:
    infoDict = allProbes[aa]
    probeSeq = probeDict[aa]
    probeLen = len(probeSeq)
    for ii in infoDict:
        ENSTid = re.sub("\..*", "", ii) 
        Longseq = transcriptomeDict[ENSTid]        
        if int(infoDict[ii][1][1]) < int(infoDict[ii][1][0]):
            end = (int(infoDict[ii][0][0]) - 1) + int(infoDict[ii][1][0])
            start = int(infoDict[ii][1][1]) - (probeLen - int(infoDict[ii][0][1]))
            if start < 1:
                start = 1
            if end > len(Longseq):
                end = len(Longseq) 
            seq = Longseq[start-1:end]
            seq = seq[::-1]
        else: continue
        Tm = MeltingTemp.Tm_NN(probeSeq, check=True, strict=True, shift=0, 
                                          nn_table=None, tmm_table=None, imm_table=None, 
                                          de_table=None, dnac1=25, dnac2=25, selfcomp=False, 
                                          Na=50, K=0, Tris=0, Mg=0, dNTPs=0, saltcorr=5)
        score = primer3.calcHeterodimer(seq, probeSeq).tm
        ScoreList.append((score, aa, ii))

PtoComp = {}
for ii in ScoreList:
    if ii[2] in PtoComp:
        PtoComp[ii[2]] += 1
    else:
        PtoComp[ii[2]] = 1


PtoCompList = list(PtoComp.items())
PtoCompList.sort(key=lambda p: p[1], reverse=True)
    
Scores = []

if probeFastaFile.find("ENST") >= 0:
    FastaFileShort= probeFastaFile.split("/")[-1]
    FastaFileName = FastaFileShort.split("_")[0]
    type = EtoN[FastaFileName]
    print(FastaFileName, type, GeneSymbolDict[FastaFileName])
else:
    FastaFileShort= probeFastaFile.split("/")[-1]
    FastaFileName = FastaFileShort.split(".")[0]
    ENST = NtoE[FastaFileName]
    print(ENST, FastaFileName, GeneSymbolDict[ENST])

for temp in range(20,40):
    numBoundtoTargetRNA = 0
    numBoundtoCompRNA = {}
    numBoundtoCompRNAFreq = {}
    for pp in ScoreList:
        if pp[2] == targetmRNA:
            if pp[0] >= temp:
                numBoundtoTargetRNA += 1
        else:
            ENSTid = pp[2].split('.')
            ENSTid = ENSTid[0] 
            if GeneSymbolDict[ENSTid] == mRNAGene:
                continue
            if pp[0] >= temp:
                if pp[2] in numBoundtoCompRNA:
                    numBoundtoCompRNA[pp[2]] += 1
                else: 
                    numBoundtoCompRNA[pp[2]] = 1
    numBoundtoCompRNAValues = list(numBoundtoCompRNA.values())
    numBoundtoCompRNAValues.sort(reverse = True)
#    idx = 0
#    for ii in numBoundtoCompRNAValues:
#        if ii == 1:
#            idx += 1
#            numBoundtoCompRNAFreq[ii] = idx
#        else:
#            numBoundtoCompRNAFreq[ii] = 1

    for ii in numBoundtoCompRNAValues:
        if ii in numBoundtoCompRNAFreq:
            numBoundtoCompRNAFreq[ii] += 1
        else:
            numBoundtoCompRNAFreq[ii] = 1

    value = 0
    for ii in numBoundtoCompRNAFreq:
        value =(((ii/50)**2) * numBoundtoCompRNAFreq[ii]) + value

    Dvalue = ((numBoundtoTargetRNA/50)**2) - value
   
#    print("T = ", str(temp) , "N = ", numBoundtoTargetRNA, "C = ", numBoundtoCompRNAFreq)
    print("T = ", str(temp) , "N = ", numBoundtoTargetRNA, "C = ", numBoundtoCompRNAFreq, "Evaluation Value = ", Dvalue, "CompetingRNAs = ", numBoundtoCompRNA)
    


#    maxComp = 0
#    maxCompRNA = ''
#    for cc in numBoundtoCompRNA:
#        if numBoundtoCompRNA[cc] > maxComp:
#            maxComp = numBoundtoCompRNA[cc]
#            maxCompRNA = cc
#    Scores.append((numBoundtoTargetRNA, maxComp, maxCompRNA,temp))
      
#print(probeFastaFile,  Scores)    
print('\n')
