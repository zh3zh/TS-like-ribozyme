import sys
import os
import argparse
from collections import defaultdict

#barcode id to seq
sBcMap = {}
#barcode id + num of reads
dBcMap = {}


revBase = {}
revBase['A'] = 'T'
revBase['T'] = 'A'
revBase['G'] = 'C'
revBase['C'] = 'G'

varSet = {}


varFullMap = defaultdict(int)
varPartMap = defaultdict(int)
varDNAMap = defaultdict(int)
runPath = ""

def reverse_seq(dnaSeq):
    rseq = ""
    for i in range(len(dnaSeq)):
        c = dnaSeq[i]
        if c == 'A':
            rseq = 'T'+rseq
        elif c == 'T':
            rseq = 'A'+rseq
        elif c == 'G':
            rseq = 'C'+rseq
        elif c == 'C':
            rseq = 'G' + rseq
        else:
            rseq = 'N' + rseq
    return rseq


def readFastaFile(fastaFile):
    seq = ""
    with open(fastaFile) as f:
        for line in f:
            if(line.startswith('>')):
                continue
            seq += line.strip()
    return seq


def mutSeq(refSeq,readSeq):
# 'num of mutations' + ' ' + WT + pos + Mut
    if(len(refSeq) != len(readSeq)):
        return 'length not equal'
    mutPosNum = 0
    varList = []
    for i in range(len(refSeq)):
        a = refSeq[i]
        b = readSeq[i]
        if a != b:
            mutPosNum += 1
            varList.append(a+str(i)+b)
    return str(mutPosNum)+' '+','.join(varList)


def numMismatch(mutinfo1, mutinfo2):
    mutinfo1_1 = mutinfo1.split(' ')
    mutinfo1_2 = mutinfo1_1[1]
    mutinfo1_3 = mutinfo1_2.split(',')
    mutinfo2_1 = mutinfo2.split(' ')
    mutinfo2_2 = mutinfo2_1[1]
    mutinfo2_3 = mutinfo2_2.split(',')
    count = 0
    for i in mutinfo1_3:
        if i in mutinfo2_3:
            count += 1
    numismatch = len(mutinfo1_3) + len(mutinfo2_3) - 2*count
    return numismatch



def generateNewBarcodeLib(bcFile, refSeq, varOut):
# bcfile: bc2Variants.txt
# varOut: var.stat
    varCount = 0
    f = open(bcFile, 'r')
    of = open(varOut, 'w')
    for line in f:
        if line.startswith('bc'):
            continue
        spt = line.split('\t')
        bcID = spt[0]
        varNum = int(spt[1])
        if varNum == 0:
            continue
        elif varNum == 1:
            seqSpt = spt[2].split('_')
            mutInfo = mutSeq(refSeq, seqSpt[0])
            # mutInfo = mutSeq: 'num of mutations' + ' ' + Nuc(WT) + pos + Nuc(Mut)
            sBcMap[bcID] = mutInfo
            dBcMap[bcID] = int(seqSpt[2])
            infoSpt = mutInfo.split(' ')
            varSet[mutInfo] = int(infoSpt[0])
            varCount += 1
        
    sBcNum = len(sBcMap)
   
    varInfoCount = [0]*100


    of.write('var: ' + str(varCount))

    for key in varSet:
        n = varSet[key]
        varInfoCount[n] += 1

    for i in range(10):
        of.write(str(i) + " " + str(varInfoCount[i]) + "\n")
    of.close()




def readFullShortCounts(RefSeq, output, discard):
 
    ext5 = 'GGCAGAGCTC'
    df = open(discard, 'w')
    for mutInfo in varSet:
        varFullMap[mutInfo] = 0.0
        varPartMap[mutInfo] = 0.0
        varDNAMap[mutInfo] = 0.0

  
    matchedRead = 0
    shortLenRead = 0
    shortLenMatchBcRead = 0
    shortLenNotMatchedBcRead = 0

    f = open(runPath + "chop_LR.txt", 'r')
    for line in f:
        spt = line.split('\t')
        diff = 0
        num_mut = 0
        seqLen = len(spt[1])

        matchedRead += 1
        shortLenRead += 1
        seq = spt[1]
        bcID = spt[0]
        rev_seq = ext5 + seq
        mutInfo = mutSeq(fasta, rev_seq)
        if bcID in sBcMap:
            mutInfo_map = sBcMap[bcID]
            mutseq = varToMutSeq(fasta, mutInfo_map)
            mutInfo_short = mutSeq(seq, mutseq[-36:])
            spt2 = mutInfo_short.split(' ')
            spt3 = spt2[1].split(',')
            mismatch = len(spt3)
            if mismatch <= 2 :
                shortLenMatchBcRead += 1
                varPartMap[mutInfo_map] += 1

                    
    for bcID in sBcMap:
        mutInfo = sBcMap[bcID]
        varDNAMap[mutInfo] += dBcMap[bcID] 
 
               

    df.write("matchedRead: " + str(matchedRead) + "\n")
    df.write("short length: " + str(shortLenRead) + " short length matched: " + str(shortLenMatchBcRead) + "\n")
    df.close()

    of = open(output,'w')
    for mutInfo in varSet:

        pCount = varPartMap[mutInfo]
        of.write(mutInfo+' '+str(varDNAMap[mutInfo])+' '+str(pCount)+'\n')
    of.close()
    f.close()



def argParseInit():
    global parser
    parser = argparse.ArgumentParser(description='count variants')
    parser.add_argument('--runPath', default='./', help='path to working directory')
    parser.add_argument('--refSeq', required=True, help='reference sequence in fasta format')

if __name__ == "__main__":

    argParseInit()
    args = parser.parse_args()
    runPath = args.runPath
    if not os.path.exists(runPath):
        os.makedirs(runPath)

    fasta = readFastaFile(args.refSeq)
    bcFileName = runPath + '/bc2Variants.txt'
    varOutput = runPath + '/var.stat'

    output1 = runPath + '/var.count'
    discard1 = runPath + '/var.discard'

    generateNewBarcodeLib(bcFileName, fasta, varOutput)
    readFullShortCounts(fasta, output1, discard1)

