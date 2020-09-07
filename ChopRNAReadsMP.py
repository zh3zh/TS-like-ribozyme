import os
import regex
import argparse
from collections import defaultdict
import gzip
from multiprocessing import Pool, cpu_count

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



umi_cluster = {}
cluster_umi = {}

fullLengthCount = []
partLengthCount = []
runPath = ""
refSeq = ""
shortlength = 36

def readUmiClusterFromFile(fileName):
    f = open(fileName, 'r')
    for line in f:
        spt = line.split('\t')
        umi = spt[0]
        id = int(spt[1])
        umi_cluster[umi] = spt[1]

        cluster_umi[id] = umi
        fullLengthCount = []*len(umi_cluster)
        partLengthCount = []*len(umi_cluster)


def readMergedFile():
    fq_read = open(runPath + "LR.assembled.fastq")
    of = open(runPath + "chop_LR.txt", 'w')

    while True:
        try:
            read_lines = [next(fq_read).decode("ascii") for i in range(4)]
        except:
            break
        else:
            m2 = read_lines[1].rstrip()
            m1 = m2[18:-17]
            umi = m1[-18:]
            if umi in umi_cluster:
                id = int(umi_cluster[umi])
                seq = m1[:-20]
                if len(seq) == shortlength:
                    of.write(str(id)+'\t'+seq+"\n")
    of.close()
    fq_read.close()





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
    readUmiClusterFromFile(runPath+'/bc.index')

    readMergedFile()
