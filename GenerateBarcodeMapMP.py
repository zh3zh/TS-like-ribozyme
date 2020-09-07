import os
import argparse
from collections import defaultdict



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


runPath = ""
refSeq = ""
umi_cluster = {}
cluster_umi = {}
consensusCutoff = 0.9

def readUmisFromMergedFile():
    umiCounts = defaultdict(int)

    fq_read = open(runPath + "LD.assembled.fastq")
    tmpOut = open(runPath + "tmp_LD.fq",'w')

    while True:
        try:
            read_lines = [next(fq_read).decode("ascii") for i in range(4)]
        except:
            break
        else:
            m2 = read_lines[1].rstrip()
            m1 = m2[19:-17]
            umi = m1[-18:]
            seq = m1[:-20]
            if (len(seq) != len(refSeq)):
                continue
            diffCount = 0
            for a,b in zip(seq, refSeq):
                if a!=b:
                    diffCount += 1
            if diffCount > 9:
                continue
            umiCounts[umi] += 1
            m4 = read_lines[3].rstrip()
            m3 = m4[19:-17]
            qseq = m3[:-20]
            tmpOut.write(umi + " " + seq + " " + qseq + "\n")

    bcOut = open(runPath + "bc.index", 'w')
    index = 0
    effMatchedReadNum = 0
    for bc in umiCounts:
        n = umiCounts[bc]
        if n > 7:
            umi_cluster[bc] = index
            cluster_umi[index] = bc
            bcOut.write(bc+'\t'+str(index)+'\t'+str(n)+'\n')
            index += 1
            effMatchedReadNum += n
    bcOut.close()

   


def writePseudoSamFile():

    samOut = open(runPath+"pseudo.sam",'w')
    matchCigar = str(len(refSeq))+'M'

    head_1 = ['@HD', 'VN:1.0', 'SO:unsorted']
    head_2 = ['@SQ', 'SN:pseudo_genome', 'LN:'+str(len(umi_cluster)+len(refSeq))]
    head_3 = ['@PG', 'ID:bowtie2', 'PN:bowtie2', 'VN:2.0.2']
    samOut.write('\t'.join(head_1)+'\n')
    samOut.write('\t'.join(head_2)+'\n')
    samOut.write('\t'.join(head_3)+'\n')

    f = open(runPath + "tmp_LD.fq", 'r')
    for line in f:
        spt = line.split(' ')
        umi = spt[0]
        if umi in umi_cluster:
            umiIndex = umi_cluster[umi]
            samLine = [umi + '_' + str(umiIndex), '0', '*', str(umiIndex), '50', matchCigar, '*', '0', '0', \
                           spt[1], spt[2]]
            samOut.write('\t'.join(samLine) + '\n')
    f.close()

    samOut.close()


def analysisMergedBarcodes(seqs, refSeq):
    atgc = ['A', 'T', 'G', 'C']
    seqNum = len(seqs)
    if(seqNum == 0):
        return []
    seqLen = len(seqs[0])
    baseCount = [[0 for col in range(4)] for row in range(seqLen)]
    for seq in seqs:
        for i in range(seqLen):
            s = seq[i]
            if s == 'A':
                baseCount[i][0] += 1
            elif s == 'T':
                baseCount[i][1] += 1
            elif s == 'G':
                baseCount[i][2] += 1
            elif s == 'C':
                baseCount[i][3] += 1
    badPosList = []
    consSeq = ''
    mutSeq = ''
    minP = 1.0
    minPCount = 0
    for i in range(seqLen):
        maxP = 0.0
        maxPCount = 0
        conBase = 'A'
        for j in range(4):
            p = baseCount[i][j]*1.0/seqNum
            if p > maxP:
                maxP = p
                maxPCount = baseCount[i][j]
                conBase = atgc[j]

        consSeq += conBase
        if (conBase != refSeq[i]):
            mutSeq += str(i)+"-"+conBase if len(mutSeq) == 0 else '-'+str(i)+'-'+conBase
        if maxP < consensusCutoff:
            badPosList.append(i)
        if maxP < minP:
            minP = maxP
            minPCount = maxPCount


    varList = []
    if len(badPosList) == 0:
        varList.append(consSeq+'_'+str(round(minP,2))+'_'+str(minPCount))
    else:

        keyCounts = defaultdict(int)

        for seq in seqs:
            key = ''
            for badPos in badPosList:
                key += seq[badPos]
            keyCounts[key] += 1
        for key in keyCounts:
            nKey = keyCounts[key]
            pKey = nKey*1.0/len(seqs)
            if (pKey > 0.5 and nKey > 9):

                varTag = ''
                vList = list(consSeq)
                for x in range(len(badPosList)):
                    vList[badPosList[x]] = key[x]
                    varTag += key[x]+str(badPosList[x])+','
                varTag = varTag[0:len(varTag)-1]
                varList.append("".join(vList)+'_'+str(round(pKey,2))+'_'+str(nKey)+'_'+varTag)
    return varList


def readSortedSamFile(inputSam, outputFile, refSeq):
    rf = open(inputSam, 'r')
    of = open(outputFile, 'w')
    head_2 = ['bcID', 'variantNum', 'variantSeq_pVar_nVar']
    of.write('\t'.join(head_2)+'\n')
    preClusterIndex = '0'

    seqs = []
    rawUmis = []
    for line in rf:
        line_split = line.split('\t')
        head_split = line_split[0].split('_')
        umi_index = head_split[1]
        rawUmi = head_split[0]
        if(umi_index != preClusterIndex and umi_index != '0'):
            vList = analysisMergedBarcodes(seqs,refSeq)

            if(len(vList) > 0):
                of.write(str(preClusterIndex)+'\t'+str(len(vList)))
                for line in vList:
                    of.write('\t'+line)
                of.write('\n')
            preClusterIndex = umi_index

            del seqs[:]
            del rawUmis[:]
        seq = line_split[9]
        if len(seq) == len(refSeq):
            seqs.append(seq)
            rawUmis.append(rawUmi)
    vList = analysisMergedBarcodes(seqs,refSeq)
    if(len(vList) > 0):
        of.write(str(preClusterIndex)+'\t'+str(len(vList)))
        for line in vList:
            of.write('\t'+line)
        of.write('\n')

    rf.close()
    of.close()

parser = None
def argParseInit():
    """
    :rtype: object
    """
    global parser
    parser = argparse.ArgumentParser(description='generate barcode map')
    parser.add_argument('--runPath', required=True, help='path to working directory')
    parser.add_argument('--refSeq', required=True, help='reference sequence in fasta format')



if __name__ == "__main__":

    argParseInit()
    args = parser.parse_args()
    runPath = args.runPath
    if not os.path.exists(runPath):
        os.makedirs(runPath)
    refSeq = readFastaFile(args.refSeq)
    readUmisFromMergedFile()
    writePseudoSamFile()
    os.system('samtools view -bS ' + runPath + 'pseudo.sam > ' + runPath + 'pseudo.bam')
    os.system('samtools sort ' + runPath + 'pseudo.bam ' + runPath + 'pseudo.sort')
    os.system('samtools view ' + runPath + 'pseudo.sort.bam >' + runPath + 'pseudo.sort.sam')
    readSortedSamFile(runPath + 'pseudo.sort.sam', runPath + 'bc2Variants.txt', refSeq)
