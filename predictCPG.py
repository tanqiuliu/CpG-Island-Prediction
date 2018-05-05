import sys
import numpy as np
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import pylab as pl
from matplotlib import collections as mc
from matplotlib.colors import ListedColormap

STATE = ['A+', 'A-', 'G+', 'G-', 'C+', 'C-', 'T+', 'T-', 'N+', 'N-']
OBS = ['A', 'G', 'C', 'T', 'N']


def load_cpg(fpath):
    """
    trn type: pandas DataFrame
    """
    columns = ['bin', 'chrom', 'chromStart', 'chromEnd', 'name', 'length', 'cpgNum', 'gcNum', 'perCpg', 'perGc',
               'obsExp']
    cpg = pd.read_csv(fpath, sep='\t', names=columns)
    return cpg

def load_cpg2(fpath):
    """
    trn type: pandas DataFrame
    """
    columns = ['chrom', 'chromStart', 'chromEnd']
    cpg = pd.read_csv(fpath, sep='\t', names=columns)
    return cpg


def sample_seq(chr_seq, cpg_df, start, end):
    """
    sample a subseq from chromosome and build cpg_df for it
    """
    seq = chr_seq[start:end]
    cpg_sample = []
    for i in range(len(cpg_df)):
        if start < cpg_df.iloc[i]['chromStart'] and cpg_df.iloc[i]['chromEnd'] < end:
            cpg_sample.append([cpg_df.iloc[i]['chromStart'] - start, cpg_df.iloc[i]['chromEnd'] - start])
        elif cpg_df.iloc[i]['chromStart'] < start < cpg_df.iloc[i]['chromEnd']:
            cpg_sample.append([0, cpg_df.iloc[i]['chromEnd'] - start])
        elif cpg_df.iloc[i]['chromStart'] < end < cpg_df.iloc[i]['chromEnd']:
            cpg_sample.append([cpg_df.iloc[i]['chromStart'] - start, end - start])
    cpg_df_sample = pd.DataFrame(cpg_sample, columns=['chromStart', 'chromEnd'])
    return seq, cpg_df_sample


def getFreq(seq, cpg_df):
    cpg_df = cpg_df.sort_values(by=['chromStart'])
    cpg_starts = cpg_df['chromStart']
    cpg_ends = cpg_df['chromEnd']
    d = {c: {s: c + s for s in ['+', '-']} for c in OBS}
    # initialize counter
    pseudo_count = 0
    transition = {n_prev: {n_next: pseudo_count for n_next in STATE} for n_prev in STATE}
    prior = {s: pseudo_count for s in STATE}
    # count
    state = '-'  # indicator of whether i is in CpG or not
    idx = 0  # index of cpg info
    state_prev = seq[0] + state
    for i in range(len(seq)):
        # if i % 1000000 == 0:
        #     print(i)
        if idx < cpg_df.shape[0]:
            if i == cpg_starts[idx]:
                state = '+'
            if i == cpg_ends[idx]:
                state = '-'
                idx += 1
        # count transition
        state_next = d[seq[i]][state]
        if state == "-":
            transition[state_prev][state_next] += 1
        if state == "+":
            transition[state_prev][state_next] += 1

        state_prev = state_next
        # count prior
        if state == "-":
            prior[d[seq[i]][state]] += 1
        if state == "+":
            prior[d[seq[i]][state]] += 1
    print prior
    return pd.DataFrame(transition), pd.Series(prior)


def getLogTransitionProb(freq):
    """
    freq: pandas DataFrame, transition frequency from count
    """
    neg_inf = 1e-30
    # freq = freq.drop(columns=['N+','N-'],index=['N+','N-'])
    prob = freq / np.sum(freq, axis=0) + neg_inf
    return np.log(prob)


def getLogPriorProb(prior_count):
    neg_inf = 1e-30
    prob = prior_count / np.sum(prior_count) + neg_inf
    print prob
    return np.log(prob)


def viterbi(seq, log_trans_prob, log_prior_prob):
    path = []
    prob_mem = {j: {i: 0 for i in STATE} for j in range(len(seq))}
    prev_state_mem = {j: {i: None for i in STATE} for j in range(len(seq))}
    state_dict = {"A": {"A+", "A-"}, "T": {"T+", "T-"}, "G": {"G+", "G-"}, "C": {"C+", "C-"}, "N": {"N+", "N-"}}
    # prior
    for s in STATE:
        prob_mem[0][s] = log_prior_prob[s]
    #
    for i in range(1, len(seq)):
        cur_state = state_dict[seq[i]]
        for cur in cur_state:
            max_prob = -5000000000
            max_prev = -1
            for prev in STATE:
                p = prob_mem[i - 1][prev] + log_trans_prob[prev][cur]  + 0
                if p >= max_prob:
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
    #print prob_mem, prev_state_mem
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
        TP = float(TP)
        FP = len(cpg_pred) - TP
        FP = float(FP)
        FN = len(cpg_gt) - TP
        FN = float(FN)
        scores.append(TP / (TP + FP + FN))
    return np.sum(scores) / len(thresholds)

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
    cpg_df = pd.DataFrame(np.array(cpg_df), columns=['chromStart', 'chromEnd'])
    return cpg_df

def getCpgInfo(cpg_df, seq):
    # extract CpG information from predicted CpG islands
    cpg_df['length'] = np.nan
    cpg_df['cpgNum'] = np.nan
    cpg_df['gcNum'] = np.nan
    cpg_df['perCpg'] = np.nan
    cpg_df['perGc'] = np.nan
    cpg_df['obsExp'] = np.nan
    info = cpg_df.values
    for idx in range(len(cpg_df)):
        start = info[idx,0]
        end = info[idx,1]
        subseq = seq[int(start):int(end)]
        length = end - start
        cpgNum = len(subseq.split('CG')) - 1
        gNum = np.sum(np.array(list(subseq)) == 'G')
        cNum = np.sum(np.array(list(subseq)) == 'C')
        obsExp = cpgNum * length / (gNum * cNum)
        info[idx,2] = length
        info[idx,3] = cpgNum
        info[idx,4] = gNum + cNum
        info[idx,5] = (2 * cpgNum / length) * 100
        info[idx,6] = ((gNum + cNum) / length) * 100
        info[idx,7] = obsExp
        # cpg_df['cpgNum'][idx] = cpgNum
        # cpg_df['gcNum'][idx] =  gNum + cNum
        # cpg_df['perCpg'][idx] = (2 * cpgNum / length) * 100
        # cpg_df['perGc'][idx] = ((gNum + cNum) / length) * 100
        # cpg_df['obsExp'][idx] = obsExp
    return pd.DataFrame(info, columns=cpg_df.columns)

def visualize(gt_cpg, pred_cpg, perGc_list, win_size):    #line1 = []
    gt_line= []
    pred_line = []

    for i in range(len(pred_cpg)):
        s = pred_cpg.iloc[i]['chromStart']
        e = pred_cpg.iloc[i]['chromEnd']
        pred_line.append(s)
        pred_line.append(e)

    pred_y = [1] * len(pred_line)
    pred_line = list(zip(pred_line, pred_y))
    pred_line = [pred_line[x:x + 2] for x in range(0, len(pred_line), 2)]

    for i in range(len(gt_cpg)):
        s = gt_cpg.iloc[i]['chromStart']
        e = gt_cpg.iloc[i]['chromEnd']
        gt_line.append(s)
        gt_line.append(e)
        #print start, end

    gt_y = [0] * len(gt_line)
    gt_line = list(zip(gt_line, gt_y))
    gt_line = [gt_line[x:x + 2] for x in range(0, len(gt_line), 2)]
    # zip into tuple

    x = np.linspace(0, win_size*len(perGc_list), len(perGc_list))
    perGc_list = list(zip(x, perGc_list/100))
    perGc_list = [perGc_list[x:x + 2] for x in range(0, len(perGc_list), 1)]

    # colors = [mcolors.to_rgba(c) for c in plt.rcParams['axes.prop_cycle'].by_key()['color']]
    gc = mc.LineCollection(gt_line,  colors = 'c', linewidths=3, label = 'Grounded CPG Islands')
    lc = mc.LineCollection(pred_line,  colors= 'm', linewidths=3, label ='Predicted CPG Islands' )
    gcp =  mc.LineCollection(perGc_list,  colors= 'k', linewidths=1, label ='CG Percent in Sequence' )
    cMap = ListedColormap(['white', 'green', 'blue', 'red'])

    #ax.legend((pred_line, gt_line), ('Predicted CPG', 'Grouded CPG'))
    fig, ax1 = pl.subplots(2, sharex=True)
    ax1[1].set_xlim(0, 100000)
    ax1[1].set_yticklabels(list(["0", "0.2", "0.4", "0.6", "0.8","1.0"]), minor=False)
    ax1[1].set_ylim(-0.5, 2.5)
    ax1[1].add_collection(lc)
    ax1[1].add_collection(gc)
    ax1[0].add_collection(gcp)
#    ax1[0].set_ylabels("Cg Percent", minor=False)
    ax1[1].set_ylim(-0.5, 2.5)

    ax1[1].legend(frameon=False, loc='upper left', shadow = 'True', title="Predicted results v.s Grounded Truth")
    #ax1.autoscale()
    ax1[0].set_title('CG Percent in Squence')
    #ax.xlabel('Sequence Index')
#    heatmap = ax1[0].pcolor(perGc_list, cmap=cMap)
   # cbar = plt.colorbar(heatmap)
  #  cbar.ax.set_yticklabels(['0', '0.25', '0.5', '> 0.75'])
   # cbar.set_label('Percentage of CG', rotation=270)

    pl.show()


def seqwiseGcPer(seq, win_size = 200):
    num_win = int(np.ceil(len(seq) / win_size))
    perGc_list = np.zeros(num_win)
    perCpg_list = np.zeros(num_win)
    for idx in range(num_win):
        subseq = seq[idx * win_size:(idx + 1) * win_size]
        perGc = float((np.sum(np.array(list(subseq) ) == 'G') + float(np.sum(np.array(list(subseq)) == 'C')) / win_size * 100))
        perCpg = float((len(subseq.split('CG')) - 1) * 2 / win_size * 100)
        perGc_list[idx] = perGc
        perCpg_list[idx] = perCpg
    return perCpg_list, perGc_list, len(seq)

def visual_info(perGc_list):    #line1 = []
    x = np.linspace(0, 100000, len(perGc_list))
    # percent_line = []
    # population = np.array([2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0])
    plt.plot(x, perGc_list, linestyle='-')
    # cbaxes = plt.axes([10, 20, 30, 40, 50, 60, 70, 80, 90,100, 110, 120])  # This is the position for the colorbar
    #    plt.colorbar()
    plt.legend(loc='upper center', shadow='True', title="GC Percentage%")
    plt.show()

   ##| x = np.linspace(0, 1000000, len(perGc_list))
    #pt = list(zip(x, list(perGc_list)))
    #gcp = mc.LineCollection(pt, linestyles='solid')
    #fig, ax = pl.subplots()
    #ax.add_collection(gcp)
    #fig = plt.gcf()
    #axcb = fig.colorbar(gcp)
    #axcb.set_label('CG percent')
    #a.sci(gcp)  # This allows interactive changing of the colormap.
    #ax.margins(0.1)
   #


    # colors = [mcolors.to_rgba(c) for c in plt.rcParams['axes.prop_cycle'].by_key()['color']]
    #gc = mc.LineCollection(gt_line,  colors = 'c', linewidths=3, label = 'Grounded CPG Islands')
    #lc = mc.LineCollection(pred_line,  colors= 'm', linewidths=3, label ='Predicted CPG Islands' )

    #ax.legend((pred_line, gt_line), ('Predicted CPG', 'Grouded CPG'))
    fig, ax = pl.subplots()
    #ax.set_xlim(0, 80000)
    #ax.set_ylim(-0.5, 1.8)
    #ax.add_collection(lc)
    #ax.add_collection(gc)
    #ax.legend(loc='upper center', shadow = 'True', title="Predicted results v.s Grounded Truth")
    #ax.autoscale()
    #ax.xlabel('Sequence Index')
    #ax.margins(0.1)




if __name__ == '__main__':
    seqPath = 'chr1.fa'
    cpgPath = 'cpgIslandExt.txt'
    #seqPath = 'training.fa'
    #cpgPath = 'cpg_island_training.txt'

    if len(sys.argv) == 3:
        seqPath = sys.argv[1]
        cpgPath = sys.argv[2]

    chr = SeqIO.read(seqPath, 'fasta')
    chr_id = chr.id
    chr_seq = chr.seq.upper()
    cpg_df = load_cpg(cpgPath)
    cpg_df_chr1 = cpg_df[cpg_df['chrom'] == chr_id]

    train_seq, train_cpg = sample_seq(chr_seq, cpg_df_chr1, start=2500000, end=3500000)
    test_seq, test_cpg = sample_seq(chr_seq, cpg_df_chr1, start=1600000, end=1700000)
    transitionFreq, priorFreq = getFreq(train_seq, train_cpg)
    log_trans_prob = getLogTransitionProb(transitionFreq)
    log_prior_prob = getLogPriorProb(priorFreq)
    pred_path, best_score = viterbi(test_seq, log_trans_prob, log_prior_prob)
    print pred_path
    pred_cpg = path2cpg(pred_path)
    #print visual(test_cpg, pred_cpg)

    print("GT:")
    print(test_cpg)
    print("PRED:")
    print(pred_cpg)
    print ("Score: ")
    print score(test_cpg, pred_cpg)
    perCpg_list, perGc_list,num = seqwiseGcPer(test_seq)
    print len(perGc_list)
    win_size = 200

    print visualize(test_cpg, pred_cpg, perGc_list, win_size)
    print visual_info(perGc_list)