#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 21:42:58 2023

@author: kmazooji
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


#file1 = open('scores-ROSE1000M1-B500-K15.txt', 'r')
#file1 = open('scores-16S3-B1000-K15.txt', 'r')
#file1 = open('scores-homfam-adh-cd-B1000-K10.txt', 'r')
#file1 = open('scores-homfam-adh-cd-B1000-K5.txt', 'r')
#file1 = open('scores-homfam-adh-cd-B1000-K15.txt', 'r')
#file1 = open('scores-homfam-adh-cd-B1000-K20.txt', 'r')
#file1 = open('scores-16ST-B1000-K20.txt', 'r')
#file1 = open('scores-16ST-B1000-K15.txt', 'r')
file1 = open('scores-16ST-B1000-K10.txt', 'r')
#file1 = open('scores-indelible10000M2-B1000-K20.txt', 'r')

data1 = file1.read()
both_scores = ast.literal_eval(data1)
both_scores1 = both_scores

K = 10

HMM_to_scores = {}
# for set1 in both_scores[list(both_scores.keys())[0]]:
#     HMM_to_scores[set1] = []

# for qseq in both_scores:
#     for set1 in both_scores[qseq]:
#         HMM_to_scores[set1].append(both_scores[qseq][set1])

for qseq in both_scores:
    for set1 in both_scores[qseq]:
        if set1 in HMM_to_scores:
            HMM_to_scores[set1].append(both_scores[qseq][set1])
        else:
            HMM_to_scores[set1] = [both_scores[qseq][set1]]
num_seqs = len(both_scores)      
HMMs = list(HMM_to_scores.keys())   
num_HMMs = len(HMMs)    
c = 0.1
R = int(np.round(num_HMMs * c))

# for j in range(len(HMMs)):
#     plt.scatter([HMM_to_scores[HMMs[j]][i][0] for i in range(num_seqs)], [HMM_to_scores[HMMs[j]][i][1] for i in range(num_seqs)])
#     plt.show()

#seq1 = 'SEQ1' 

if False:
    for seq1 in list(both_scores.keys())[-50:-1]:
        plt.scatter([both_scores[seq1][i][0] for i in range(len(HMMs))], [both_scores[seq1][i][1] for i in range(len(HMMs))])
        plt.title(seq1)
        plt.xlabel("bitscore")
        plt.ylabel("J-score")
        plt.show()
#plt.scatter([both_scores_vector[i][0] for i in range(len(both_scores_vector))], [both_scores_vector[i][1] for i in range(len(both_scores_vector))])



query_seqs = {}
for qseq in both_scores:
    query_seqs[qseq] = True

qseq_to_score = {} # holds biscores
query_set_comparison = {} # holds J-scores

for qseq in query_seqs:
    query_set_comparison[qseq] = {}
    qseq_to_score[qseq] = {}
    for hmm in both_scores[qseq].keys():        
        query_set_comparison[qseq][hmm] = [None, None, both_scores[qseq][hmm][1]]
        qseq_to_score[qseq][hmm]  = both_scores[qseq][hmm][0]


    
    
query_assignment = {}
set_to_queries = {}
for qseq in query_set_comparison:
#for qseq in list(qseq_to_kmers.keys())[0:2]:
    best_seen = None
    best_val = 0
    for set_num in query_set_comparison[qseq]:
        if query_set_comparison[qseq][set_num][2] > best_val:
            best_seen = set_num
            best_val = query_set_comparison[qseq][set_num][2]
            #print(best_seen, best_val)
    
    query_assignment[qseq] = best_seen
    
    # if best_val > 0:
    #     query_assignment[qseq] = best_seen
    # else:
    #     query_assignment[qseq] = 0
    
    #print("")
    #print(qseq)
    #print(best_seen)
    #print(qseq_cands[qseq])
    
    if best_seen in set_to_queries:
        set_to_queries[best_seen].append(qseq)
    else:
        set_to_queries[best_seen] = [qseq]
     
        
     
        
     







query_assignment_bitscore = {}
set_to_queries_bitscore = {}
for qseq in query_seqs:
#for qseq in list(qseq_to_kmers.keys())[0:2]:
    best_seen = None
    best_val = -np.infty
    for set_num in qseq_to_score[qseq]:
        if qseq_to_score[qseq][set_num] > best_val:
            best_seen = set_num
            best_val = qseq_to_score[qseq][set_num] 
            #print(best_seen, best_val)
    
    query_assignment_bitscore[qseq] = best_seen
    
    # if best_val > 0:
    #     query_assignment_bitscore[qseq] = best_seen
    # else:
    #     query_assignment[qseq] = 0
    
    #print("")
    #print(qseq)
    #print(best_seen)
    #print(qseq_cands[qseq])
    
    if best_seen in set_to_queries_bitscore:
        set_to_queries_bitscore[best_seen].append(qseq)
    else:
        set_to_queries_bitscore[best_seen] = [qseq]

same_assignment_qseqs = {}
# for qseq in query_assignment:
#     #print(query_assignment[qseq] , query_assignment_bitscore[qseq])
#     if query_assignment[qseq]  == query_assignment_bitscore[qseq]:
#         same_assignment_qseqs[qseq] = True
for qseq in query_seqs:
    #print(query_assignment[qseq] , query_assignment_bitscore[qseq])
    if query_assignment[qseq]  == query_assignment_bitscore[qseq]:
        same_assignment_qseqs[qseq] = True


#R = 5
same_assignment_qseqs = {}
in_top_R_Mscore_qseqs = {}
in_top_R_bitscore_qseqs = {}
qseqs_to_top_jaccard = {}
qseqs_to_int_score = {}
qseqs_to_spearman = {}
qseqs_to_spearman_top = {}
for qseq in query_assignment:
    #print(query_assignment[qseq] , query_assignment_bitscore[qseq])
    cands_bit = qseq_to_score[qseq]
    cands_sorted_bit = sorted(cands_bit.items(), key = lambda kv: kv[1], reverse=True)  
    top_sets_bit = set()
    
    cands = query_set_comparison[qseq]
    cands_sorted = sorted(cands.items(), key = lambda kv: kv[1][2], reverse=True)  
    top_sets = set()
    
    for ts in range(R):
        set_num = cands_sorted[ts][0]  
        top_sets.add(set_num)

    for ts in range(R):
        set_num = cands_sorted_bit[ts][0]  
        top_sets_bit.add(set_num)

    if query_assignment[qseq]  == query_assignment_bitscore[qseq]:
        same_assignment_qseqs[qseq] = True
    # if query_assignment[qseq] in top_sets:
    #     in_top_R_qseqs[qseq] = True
    if query_assignment_bitscore[qseq] in top_sets:
        in_top_R_Mscore_qseqs[qseq] = True
    
    top_union = top_sets.union(top_sets_bit)
    qseqs_to_top_jaccard[qseq] = len(top_sets.intersection(top_sets_bit)) / len(top_union)
    
    bitscores = [qseq_to_score[qseq][hmm] for hmm in qseq_to_score[qseq]]
    Jscores = [query_set_comparison[qseq][hmm][2] for hmm in qseq_to_score[qseq]]
    qseqs_to_spearman[qseq] = scipy.stats.spearmanr(bitscores, Jscores)
    
    '''
    bitscores_top = [qseq_to_score[qseq][hmm] for hmm in top_union]
    Jscores_top = [query_set_comparison[qseq][hmm][2] for hmm in top_union]
    qseqs_to_spearman_top[qseq] = scipy.stats.spearmanr(bitscores_top, Jscores_top)
    '''
    
    
    
    
    
        
        
print("")

print(len(same_assignment_qseqs), len(query_assignment))
print(len(same_assignment_qseqs) / len(query_assignment))

print("")
print(len(in_top_R_Mscore_qseqs ), len(query_assignment))        
print(len(in_top_R_Mscore_qseqs ) / len(query_assignment))    

avg_jaccard = 0
for qseq in qseqs_to_top_jaccard:
    avg_jaccard += qseqs_to_top_jaccard[qseq]
avg_jaccard = avg_jaccard / len(qseqs_to_top_jaccard)
print("")
print(avg_jaccard)

print("")
avg_spearman = 0
avg_spearman_pvalue = 0
# for qseq in qseqs_to_spearman:
#     avg_spearman += qseqs_to_spearman[qseq].correlation
#     avg_spearman_pvalue += qseqs_to_spearman[qseq].pvalue
for qseq in qseqs_to_spearman:
    if not np.isnan(qseqs_to_spearman[qseq].correlation):
        avg_spearman += qseqs_to_spearman[qseq].correlation
        avg_spearman_pvalue += qseqs_to_spearman[qseq].pvalue
    else:
        avg_spearman += 0
        avg_spearman_pvalue += 1      
avg_spearman = avg_spearman / len(qseqs_to_spearman)
avg_spearman_pvalue = avg_spearman_pvalue / len(qseqs_to_spearman)
print(avg_spearman)
print(avg_spearman_pvalue)

'''
print("")
avg_spearman_top = 0
avg_spearman_pvalue_top = 0
for qseq in qseqs_to_spearman_top:
    avg_spearman_top += qseqs_to_spearman_top[qseq].correlation
    avg_spearman_pvalue_top += qseqs_to_spearman_top[qseq].pvalue
avg_spearman_top = avg_spearman_top / len(qseqs_to_spearman_top)
avg_spearman_pvalue_top = avg_spearman_pvalue_top / len(qseqs_to_spearman_top)
print(avg_spearman_top)
print(avg_spearman_pvalue_top)        
'''

for qseq in list(query_seqs.keys())[0:30]:
#for qseq in list(query_seqs.keys()):
    bitscores = [qseq_to_score[qseq][hmm] for hmm in qseq_to_score[qseq]]
    Jscores = [query_set_comparison[qseq][hmm][2] for hmm in qseq_to_score[qseq]]
    plt.scatter(bitscores, Jscores)
    #plt.title(qseq + " (length = " + str(len(query_seqs[qseq])) + ")")
    plt.title(qseq)
    plt.ylabel("J-score (K = " + str(K) + ")" )
    plt.xlabel("bitscore")
    plt.show()  


