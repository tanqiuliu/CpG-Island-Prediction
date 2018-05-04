import sys
import numpy as np
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import pylab as pl
from matplotlib import collections as mc
from matplotlib import colors as mcolors

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
    prior = {s:pseudo_count for s in STATE}
    # count
    state = '-'   # indicator of whether i is in CpG or not
    idx = 0     # index of cpg info
    state_prev = seq[0] + state
    for i in range(len(seq)):
        # if i % 1000000 == 0:
        #     print(i)
        if idx < cpg_df.shape[0]:
            if i == cpg_starts[idx]:
                state = '+'
            if i ==  cpg_ends[idx]:
                state = '-'
                idx += 1
        # count transition
        state_next = d[seq[i]][state]
        transition[state_prev][state_next] += 1
        state_prev = state_next
        # count prior
        prior[d[seq[i]][state]] += 1
    return pd.DataFrame(transition), pd.Series(prior)


def getLogTransitionProb(freq):
    """
    freq: pandas DataFrame, transition frequency from count
    """
    neg_inf = 1e-30
    #freq = freq.drop(columns=['N+','N-'],index=['N+','N-'])
    prob = freq / np.sum(freq, axis = 0) + neg_inf
    return np.log(prob)


def getLogPriorProb(prior_count):
    neg_inf = 1e-30
    prob = prior_count / np.sum(prior_count) + neg_inf
    return np.log(prob)


def viterbi(seq, log_trans_prob, log_prior_prob):
    path = []
    prob_mem = {j : { i : 0 for i in STATE} for j in range(len(seq))}
    prev_state_mem = {j : { i : None for i in STATE} for j in range(len(seq))}
    state_dict = {"A" : {"A+", "A-"}, "T" : {"T+", "T-"}, "G" : {"G+", "G-"}, "C" : {"C+", "C-"}}
    # prior
    for s in STATE:
        prob_mem[0][s] = log_prior_prob[s]
    #
    for i in range(1,len(seq)):
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
    result_path = path[::-1]
    return result_path, best_score


def path2cpg(result_path):
    cpg_df = []
    cpg_state = ''
    for idx in range(len(result_path)):
        if result_path[idx][1] == '+' and cpg_state != '+':
            cpg_state = '+'
            start = idx
        if result_path[idx][1] == '-' and cpg_state == '+':
            cpg_state = '-'
            end = idx
            cpg_df.append([start, end])
    if cpg_state == '+':
        cpg_df.append([start, len(result_path)])
    return pd.DataFrame(np.array(cpg_df), columns=['chromStart', 'chromEnd'])



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


def visualize(gt_cpg, pred_cpg):
    #line1 = []
    gt_line\
        = []
    pred_line = []

    for i in range(len(pred_cpg)):
        s = pred_cpg.iloc[i]['chromStart']
        e = pred_cpg.iloc[i]['chromEnd']
        pred_line.append(s)
        pred_line.append(e)

    pred_y = [44] * len(pred_line)
    pred_line = list(zip(pred_line, pred_y))
    pred_line = [pred_line[x:x + 2] for x in range(0, len(pred_line), 2)]

    for i in range(len(gt_cpg)):
        s = gt_cpg.iloc[i]['chromStart']
        e = gt_cpg.iloc[i]['chromEnd']
        gt_line.append(s)
        gt_line.append(e)
        #print start, end

    gt_y = [45] * len(gt_line)
    gt_line = list(zip(gt_line, gt_y))
    gt_line = [gt_line[x:x + 2] for x in range(0, len(gt_line), 2)]
    # zip into tuple

    # merge into two tuple list
    #ground = [[(0, 50), (2030, 50)], [(13290, 50), (13514, 50)], [(13949, 50), (14471, 50)]]
    #lines = [[(0, 40), (2093, 40)], [(5759, 40), (6003, 40)], [(18348, 40), (18681, 40)]]
    # c = np.array([(1, 1, 0, 0), (0, 0, 1, 0)])
    colors = [mcolors.to_rgba(c) for c in plt.rcParams['axes.prop_cycle'].by_key()['color']]


    lc = mc.LineCollection(pred_line, colors=colors, linewidths=3, label ='Predicted CPG' )
    gc = mc.LineCollection(gt_line, colors=colors, linewidths=3, label = 'Grounded CPG')

    #ax.legend((pred_line, gt_line), ('Predicted CPG', 'Grouded CPG'))
    fig, ax = pl.subplots()
    ax.set_xlim(0, 95550)
    ax.add_collection(lc)
    ax.add_collection(gc)
    #ax.autoscale()
    ax.margins(0.1)
    pl.show()


# TO DO:
# post processing
# further improvement

def load2():
    train = SeqIO.read('./data/gene_data/training.txt','fasta')
    test = SeqIO.read('./data/gene_data/testing.txt','fasta')
    train_seq = train.seq.upper()
    test_seq = test.seq.upper()
    columns = ['chromStart', 'chromEnd']
    train_cpg = pd.read_csv('./data/gene_data/cpg_island_training.txt',sep=' ', names = columns)
    test_cpg = pd.read_csv('./data/gene_data/cpg_island_testing.txt',sep=' ', names = columns)
    return train_seq, test_seq, train_cpg, test_cpg





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

    train_seq, train_cpg = sample_seq(chr_seq, cpg_df_chr1, start=2000000, end=3000000) 
    test_seq, test_cpg = sample_seq(chr_seq, cpg_df_chr1, start=1000000, end=1100000) 

    # train_seq, test_seq, train_cpg, test_cpg = load2()

    transitionFreq, priorFreq = getFreq(train_seq, train_cpg)
    log_trans_prob = getLogTransitionProb(transitionFreq)
    log_prior_prob = getLogPriorProb(priorFreq)
    pred_path, best_score = viterbi(test_seq, log_trans_prob, log_prior_prob)
    pred_cpg = path2cpg(pred_path)
    # print("GT:")
    # print(test_cpg)
    # print("PRED:")
    # print(pred_cpg)
    print(score(test_cpg, pred_cpg, thresholds = [0.1]))
    visualize(test_cpg, pred_cpg)


