# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 12:02:00 2024

@author: kmazo
"""

import matplotlib.pyplot as plt
import ast

import numpy as np
import numpy.matlib
import numpy.linalg as linalg
import random
import math
import scipy
import scipy.stats



file1 = open('scores-AMINO-B100-K9.txt', 'r')
#file1 = open('scores-16S3-B1000-K15.txt', 'r')
data1 = file1.read()
both_scores1 = ast.literal_eval(data1)
both_scores = both_scores1



HMM_to_scores = {}
for set1 in both_scores[list(both_scores.keys())[0]]:
    HMM_to_scores[set1] = []

for qseq in both_scores:
    for set1 in both_scores[qseq]:
        HMM_to_scores[set1].append(both_scores[qseq][set1])
num_seqs = len(both_scores)      
HMMs = list(HMM_to_scores.keys())        

## plot for HMMs
# for j in range(len(HMMs)):
#     plt.scatter([HMM_to_scores[HMMs[j]][i][0] for i in range(num_seqs)], [HMM_to_scores[HMMs[j]][i][1] for i in range(num_seqs)])
#     plt.show()

seq1 = 'SEQ1' 
soi = ['SEQ68', 'SEQ283']
dpi1 = 300

for seq1 in soi: #list(both_scores.keys())[-50:-1]:
    fig = plt.figure( dpi=dpi1)
    plt.scatter([both_scores[seq1][i][0] for i in range(len(HMMs))], [both_scores[seq1][i][1] for i in range(len(HMMs))])
    plt.title(seq1)
    plt.xlabel("bitscore")
    plt.ylabel("J-score")
    plt.xlim([0, 120])
    plt.ylim([-0.001, 0.03])
    plt.show()
    
#plt.scatter([both_scores_vector[i][0] for i in range(len(both_scores_vector))], [both_scores_vector[i][1] for i in range(len(both_scores_vector))])