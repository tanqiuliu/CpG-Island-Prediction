import sys
import numpy as np
from Bio import SeqIO
import pandas as pd
from util import *

STATE = ['A+', 'A-', 'G+', 'G-', 'C+', 'C-', 'T+', 'T-', 'N+', 'N-'] 
OBS = ['A', 'G', 'C', 'T', 'N']


def load_cpg(fpath):
    """
    trn type: pandas DataFrame
    """
    columns = ['bin','chrom','chromStart','chromEnd','name','length','cpgNum','gcNum','perCpg','perGc','obsExp']
    cpg = pd.read_csv(fpath,sep='\t', names = columns)
    return cpg


def count(seq, cpg_df):
    cpg_df = cpg_df.sort_values(by=['chromStart'])
    cpg_starts = cpg_df['chromStart']
    cpg_ends = cpg_df['chromEnd']
    # initialize counter
    pseudo_count = 1
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
        state_next = seq[i] + state
        transition[state_prev][state_next] += 1
        state_prev = state_next
    return transition
        



if __name__ == '__main__':
    if len(sys.argv) == 3:
        seqPath = sys.argv[1]
        cpgPath = sys.argv[2]
    else:
        seqPath = 'data/chr1.fa'
        cpgPath = 'data/hg38-cpgIslandExt.txt'
    
    chr = SeqIO.read(seqPath,'fasta').seq.upper()
    cpg_df = load_cpg(cpgPath)
    cpg_df_chr1 = cpg_df[cpg_df['chrom'] == chr.id]
