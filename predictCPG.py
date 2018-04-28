import sys
import numpy as np
from Bio import SeqIO
import pandas as pd

STATE = ['A+', 'A-', 'G+', 'G-', 'C+', 'C-', 'T+', 'T-', 'N+', 'N-'] 
OBS = ['A', 'G', 'C', 'T', 'N']


def load_cpg(fpath):
    """
    trn type: pandas DataFrame
    """
    columns = ['bin','chrom','chromStart','chromEnd','name','length','cpgNum','gcNum','perCpg','perGc','obsExp']
    cpg = pd.read_csv(fpath,sep='\t', names = columns)
    return cpg


def getFreq(seq, cpg_df):
    cpg_df = cpg_df.sort_values(by=['chromStart'])
    cpg_starts = cpg_df['chromStart']
    cpg_ends = cpg_df['chromEnd']
    d = {c:{s:c+s for s in ['+','-']} for c in OBS}
    # initialize counter
    pseudo_count = 0
    transition = {n_prev:{n_next:pseudo_count for n_next in STATE} for n_prev in STATE}
    # count
    state = '-'   # indicator of whether i is in CpG or not
    idx = 0     # index of cpg info
    state_prev = seq[0] + state
    for i in range(len(seq)):
        if i % 1000000 == 0:
            print(i)
        if i == cpg_starts[idx]:
            state = '+'
        if i ==  cpg_ends[idx]:
            state = '-'
            idx += 1
        state_next = d[seq[i]][state]
        transition[state_prev][state_next] += 1
        state_prev = state_next
    return pd.DataFrame(transition)


def getTransitionProb(trans_freq):
    """
    trans_freq: pandas DataFrame, transition frequency from count
    """

def getObs(data):
    global OBSSEQ
    dataFile = open(data,'r')
    for line in dataFile:
        OBSSEQ += line.rstrip()
    dataFile.close()
    return OBSSEQ


state_dict = {"A" : {"A+", "A-"}, "T" : {"T+", "T-"}, "G" : {"G+", "G-"}, "C" : {"C+", "C-"}}

def viterbi(seq, log_trans_prob):

    path = []
    prob_mem = {j:{i:0 for i in state} for j in range(len(obsSeq))}
    prev_state_mem = {j:{i:0 for i in state} for j in range(len(obsSeq))}

    for i in xrange(len(seq)):

        for c in seq:
            cur_state = state_dict[c]

            for cur in cur_state:
                max_prob = -1
                max_prev = 0

                for prev in STATE:
                    p = prob_mem[i - 1][prev] + log_trans_prob[prev][cur] + 0

                    if p > max_prob:
                        max_prob = p
                        max_prev = prev

            prob_mem[i][cur] = max_prob
            prev_state_mem[i][cur] = max_prev


    cur_prob = max(prob_mem[len(seq) - 1].values())
    best_score = max_prob

    for s in cur_prob.keys():
        if cur_prob[s]==cur_prob:
            max_state = s

    path.append(max_state)

    for i in range(len(seq) - 1, 0, -1):
        prev_max_state = prev_state_mem[i][max_state]
        path.append(prev_max_state)
        max_state = prev_max_state

    result_path = []

    for i in range(len(seq) - 1, -1, -1):
        result_path.append(path[i])

    return result_path, best_score


if __name__ == '__main__':
    if len(sys.argv) == 3:
        seqPath = sys.argv[1]
        cpgPath = sys.argv[2]
    else:
        seqPath = 'data/chr1.fa'
        cpgPath = 'data/hg38-cpgIslandExt.txt'
    
    chr = SeqIO.read(seqPath,'fasta').seq.upper()
    cpg_df = load_cpg(cpgPath)
    cpg_df_chr1 = cpg_df[cpg_df['chrom'] == 'chr1']
