# CpG-Island-Prediction

#### Team members: Tanqiu Liu, Wenjie Zhu

## Background: 
#### (1) CpG island Background
CpG sites are sites in DNA sequence where a cytosine (“C”) nucleotide is followed by a guanine (“G”) nucleotide. Cytosine in CpG sites are usually methylated by DNA-methyltransferases. Methylated CpG sites affect the expression level of genes they related to and play an important role in gene regulation network in mammalian cells. CpG islands are regions (usually > 200 bp in length) of DNA sequence with high frequency of CpG sites. CpG island is usually associated with gene promoters and therefore becomes an important feature in gene/promoter prediction and epigenetic analysis. 

#### (2) HMM & Viterbi
The Hidden Markov Model (HMM) is a statistical model for sequences of discrete symbols. HMM is defined by Alphabets, Sequence, i-th letter in sequence, and a set of sequence. In HMM the states in the machine are not directly visible but some output at certain states are observable. Each state has a probability distribution over possible output states.
Viterbi algorithm is a dynamic programming algorithm for finding the most likely sequence of hidden states. for each intermediate state, until it reaches the end state. At each time only the most likely path leading to each state survives.

## Goal:
Implement Viterbi algorithm to predict CpG islands within a DNA sequence. 

## Data: 
We use human genome assembly hg38 with annotation to build the training sequence data and CpG island label for our project. The dataset can be downloaded from:
http://hgdownload.soe.ucsc.edu/downloads.html   
