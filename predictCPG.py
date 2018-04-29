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
    #freq = freq.drop(columns=['N+','N-'],index=['N+','N-'])
    prob = freq / np.sum(freq, axis = 0) + neg_inf
    return np.log(prob)




def viterbi(seq, log_trans_prob, log_prior_prob):
    path = []
    prob_mem = {j : { i : 0 for i in STATE} for j in range(len(seq))}
    prev_state_mem = {j : { i : None for i in STATE} for j in range(len(seq))}
    state_dict = {"A" : {"A+", "A-"}, "T" : {"T+", "T-"}, "G" : {"G+", "G-"}, "C" : {"C+", "C-"}}
    # prior
    for s in STATE:
        prob_mem[0][s] = prior[s]
    #
    for i in range(1,len(seq)):
        print(i)
        cur_state = state_dict[seq[i]]
        for cur in cur_state:
            max_prob = -np.inf
            max_prev = -1
            for prev in STATE:
                p = prob_mem[i - 1][prev] + log_trans_prob[prev][cur] + 0
                if p > max_prob:
                    max_prob = p
                    max_prev = prev         
            prob_mem[i][cur] = max_prob
            prev_state_mem[i][cur] = max_prev
        for s in STATE:
            if s not in cur_state:
                prob_mem[i][s] = -np.inf
    #
    cur_prob = prob_mem[len(seq) - 1]
    max_prob = max(cur_prob.values())
    best_score = max_prob
    for s in cur_prob.keys():
        if cur_prob[s] == max_prob:
            max_state = s
    path.append(max_state)
    #
    for i in range(len(seq) - 1, 0, -1):
        prev_max_state = prev_state_mem[i][max_state]
        path.append(prev_max_state)
        max_state = prev_max_state
    #
    result_path = []
    for i in range(len(seq) - 1, -1, -1):
        result_path.append(path[i])
    return result_path, best_score


def iou(start1, end1, start2, end2):
    if start1 < end1 < start2 < end2:
        return 0
    elif start2 < end2 < start1 < end1: 
        return 0
    else:
        intersect = min(end1, end2) - max(start1, start2)
        union = max(end1, end2) - min(start1, start2)
        return intersect / union


def score(cpg_gt, cpg_pred, thresholds=[0.5]):
    """
    cpg_gt = train_cpg.iloc[0:5]
    aa = np.array([[29000,29500],[134800,135300],[200000,200400],[350000,360000],[380000,383000]])
    cpg_pred = pd.DataFrame(aa,columns=['chromStart','chromEnd'])
    """
    iou_matrix = np.zeros((len(cpg_gt), len(cpg_pred)))
    for i in range(len(cpg_gt)):
        for j in range(len(cpg_pred)):
            iou_matrix[i,j] = iou(cpg_gt.iloc[i]['chromStart'], cpg_gt.iloc[i]['chromEnd'], cpg_pred.iloc[j]['chromStart'], cpg_pred.iloc[j]['chromEnd'])
    scores = []
    for threshold in thresholds:
        match = (iou_matrix > threshold) * 1
        TP = np.sum(np.sum(match))
        FP = len(cpg_pred) - TP
        FN = len(cpg_gt) - TP
        scores.append(TP / (TP + FP + FN))
    return np.sum(scores) / len(thresholds)


# TO DO:
# post processing
# further improvement

    




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
    test_seq, test_cpg = sample_seq(chr_seq, cpg_df_chr1, start=1000000, end=1001000) 
    freq = getFreq(train_seq, train_cpg)
    log_trans_prob = getLogTransitionProb(freq)
    pred_path, best_score = viterbi(test_seq, log_trans_prob)


