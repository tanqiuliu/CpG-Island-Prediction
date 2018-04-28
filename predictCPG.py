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


def sample_seq(chr_seq, cpg_df, start, end):
    """
    sample a subseq from chromosome and build cpg_df for it
    """
    seq = chr_seq[start:end]
    cpg_sample = []
    for i in range(len(cpg_df)):
        if start < cpg_df.iloc[i]['chromStart'] and cpg_df.iloc[i]['chromEnd'] < end:
            cpg_sample.append([cpg_df.iloc[i]['chromStart']-start, cpg_df.iloc[i]['chromEnd']-start])
        elif cpg_df.iloc[i]['chromStart'] < start < cpg_df.iloc[i]['chromEnd']:
            cpg_sample.append([0, cpg_df.iloc[i]['chromEnd']-start])
        elif cpg_df.iloc[i]['chromStart'] < end < cpg_df.iloc[i]['chromEnd']:
            cpg_sample.append([cpg_df.iloc[i]['chromStart']-start, end-start])
    cpg_df_sample = pd.DataFrame(cpg_sample, columns=['chromStart','chromEnd'])
    return seq, cpg_df_sample


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


def getLogTransitionProb(freq):
    """
    freq: pandas DataFrame, transition frequency from count
    """
    neg_inf = 1e-10
    freq = freq.drop(columns=['N+','N-'],index=['N+','N-'])
    prob = freq / np.sum(freq, axis = 0) + neg_inf
    return np.log(prob)


def score(cpg_gt, cpg_pred):
    pass



if __name__ == '__main__':
    if len(sys.argv) == 3:
        seqPath = sys.argv[1]
        cpgPath = sys.argv[2]
    else:
        seqPath = 'data/chr1.fa'
        cpgPath = 'data/hg38-cpgIslandExt.txt'
    
    chr = SeqIO.read(seqPath,'fasta')
    chr_id = chr.id
    chr_seq = chr.seq.upper()
    cpg_df = load_cpg(cpgPath)
    cpg_df_chr1 = cpg_df[cpg_df['chrom'] == chr_id]

    train_seq, train_cpg = sample_seq(chr_seq, cpg_df_chr1, start=0, end=1000000) 
    test_seq, test_cpg = sample_seq(chr_seq, cpg_df_chr1, start=1000000, end=2000000) 
    freq = getFreq(train_seq, train_cpg)
    log_trans_prob = getLogTransitionProb(freq)


