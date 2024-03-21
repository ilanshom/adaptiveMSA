#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 22 17:02:12 2023

@author: kmazooji
"""
import os
import numpy as np
#import subprocess
# import Bio
# from Bio import Phylo
# from Bio import SeqIO
# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord

final_output_directory = '/home/mazooji2/UPP-datasets-runs/final_experiments_amino/'   
output_directory = '/home/mazooji2/UPP-datasets-runs/homfam/'                     
input_directory =  '/home/mazooji2/UPP-datasets/homfam/'

FastSP_folder = '/home/mazooji2/FastSP-code/FastSP/'

first_run = False

#datas = ["adh", "PDZ", "p450"] #, "p450"]
datas = os.listdir(input_directory)

#K = 25
#Ks = [i for i in range(5, 16)]
Ks = [5, 10, 15, 20]
frac_B = 0.2
rounds = 3
stc = 10    
backbone = 1000
processors = 24


info = {}
times = {}
ref_comp = {}
num_seqs = {}
seq_len_stats = {}
sp_errors = {}
tc_scores = {}
sp_scores = {}
modeler_scores = {}
#datas = ["zf-CCHH", "tRNA-synt_2b", "aat_test-only", "p450_test-only"]
#datas = os.listdir(input_directory)
#datas  = ["10000M2"] #, "16S.T", "16S.3"]
#datas  = ["16S.T", "16S.3"]
#Ks = [5, 10, 15, 20, 25]


for data in datas:

    reference_file = input_directory + data + '/model/true.reduced.fasta'
    
    reference_record_ids = {}
    f4 = open(reference_file)
    for line in f4.readlines():
        if len(line) > 0:
            if line[0] == ">":
                reference_record_ids[line[1:].strip()] = True
 
                
        
        
    # reference_records = Bio.SeqIO.parse(reference_file, "fasta")
    # reference_record_ids = {}
    # for rec in reference_records:
    #     reference_record_ids[rec.id] = True

    seqs_file = input_directory + data + '/model/initial.fas'
    
    

    # f0 = open(seqs_file)
    # num_lines = 0
    # for line in f0.readlines():
    #     if len(line) > 0:
    #         if line[0] == ">":
    #             num_lines += 1
    # num_seqs[data] = num_lines    


    f0 = open(seqs_file)
    num_lines = 0
    i = 0
    lens = []
    lines = f0.readlines()
    for line in lines:
        if len(line) > 0:
            if line[0] == ">":
                num_lines += 1
                j = i +1
                len1 = 0
                while j < len(lines):
                    if lines[j][0] == ">":
                        break
                    len1 += len(lines[j].strip().replace("-", ""))
                    j += 1
                lens.append(len1)
        i += 1
    num_seqs[data] = num_lines   
    seq_len_stats[data] = [min(lens), max(lens), np.mean(lens), lens]

    # compare upp without pasta run
    
    output_subdirectory = output_directory + data + "UPP_bb-B" + str(backbone) + "-x" + str(processors) + "/"

    output_file = output_subdirectory + "output_alignment.fasta"
    
    output_subset_file = output_subdirectory + "output_subset_alignment.fasta"
    

    if first_run:
    #if 1:    
        f0 = open(output_file, 'r')
        f1 = open(output_subset_file, 'w')
        lines = f0.readlines()
        for i in range(len(lines)):
            line = lines[i]
            if len(line) > 0:
                if line[0] == ">":
                    if line[1:].strip() in reference_record_ids:
                        f1.write(lines[i])
                        f1.write(lines[i+1])
        f1.close()

        os.system("script -c \"java -jar " + FastSP_folder + "FastSP.jar -r " + reference_file + " -e " + output_subset_file + " \" " + output_subdirectory + "reference_comparison_scores.txt")

    terminal_file = output_subdirectory + "terminal_output.txt"

    times[data] = {}
    ref_comp[data] = {}
    sp_errors[data] = {}
    tc_scores[data] = {}
    sp_scores[data] = {}
    modeler_scores[data] = {}
    
    alg_name = "UPP_bb"

	
    f2 = open(terminal_file)
    for line in f2.readlines(): 
        if len(line.split("line 197")) > 1:
            time_taken = float(line.split("line 197")[1].split(" ")[-2])
    times[data][alg_name] = time_taken   

    ref_comp[data] = {}	
    comp_file = output_subdirectory + "reference_comparison_scores.txt"
    f3 = open(comp_file)
    spfn = 0
    spfp = 0
    tc_score = 0
    for line in f3.readlines(): 
        if len(line.split()) > 0:
            if line.split()[0] == "SP-Score":
                sp_ref = round(float(line.split()[1]), 3)
            if line.split()[0] == "SPFN":
                spfn = round(float(line.split()[1]), 3)
            if line.split()[0] == "SPFP":
                spfp = round(float(line.split()[1]), 3)
            if line.split()[0] == "TC":
                tc_score = round(float(line.split()[1]), 3)
    ref_comp[data][alg_name] = sp_ref
    sp_errors[data][alg_name] = (spfn + spfp)/2
    tc_scores[data][alg_name] = tc_score
    sp_scores[data][alg_name] = 1 - spfn
    modeler_scores[data][alg_name]  = 1 - spfp

    
    

    # compare upp2 same backbone
    
    output_subdirectory = output_directory + data + "-UPP2-common_backbone-early_stop-x" + str(processors) + "/"

    output_file = output_subdirectory + "output_alignment.fasta"

    output_subset_file = output_subdirectory + "output_subset_alignment.fasta"
    alg_name = "UPP2_bb"
    
    if first_run:
    #if 1:
        f0 = open(output_file, 'r')
        f1 = open(output_subset_file, 'w')
        lines = f0.readlines()
        for i in range(len(lines)):
            line = lines[i]
            if len(line) > 0:
                if line[0] == ">":
                    if line[1:].strip() in reference_record_ids:
                        f1.write(lines[i])
                        f1.write(lines[i+1])
        f1.close()

        os.system("script -c \"java -jar " + FastSP_folder + "FastSP.jar -r " + reference_file + " -e " + output_subset_file + " \" " + output_subdirectory + "reference_comparison_scores.txt")

    terminal_file = output_subdirectory + "terminal_output.txt"
    f2 = open(terminal_file)
    for line in f2.readlines(): 
        if len(line.split("line 197")) > 1:
            time_taken = float(line.split("line 197")[1].split(" ")[-2])
    times[data][alg_name] = time_taken  
    
    comp_file = output_subdirectory + "reference_comparison_scores.txt"
    f3 = open(comp_file)
    spfn = 0
    spfp = 0
    tc_score = 0
    for line in f3.readlines(): 
        if len(line.split()) > 0:
            if line.split()[0] == "SP-Score":
                sp_ref = round(float(line.split()[1]), 3)
            if line.split()[0] == "SPFN":
                spfn = round(float(line.split()[1]), 3)
            if line.split()[0] == "SPFP":
                spfp = round(float(line.split()[1]), 3)
            if line.split()[0] == "TC":
                tc_score = round(float(line.split()[1]), 3)
    ref_comp[data][alg_name] = sp_ref
    sp_errors[data][alg_name] = (spfn + spfp)/2
    tc_scores[data][alg_name] = tc_score
    sp_scores[data][alg_name] = 1 - spfn
    modeler_scores[data][alg_name]  = 1 - spfp



    

    for K in Ks:  
    
        # compare upp_fp5
        
        
        output_subdirectory = final_output_directory+ data + "-fp5-B" + str(backbone) + "-K" + str(K) + "-frac-" + str(frac_B) + "-stc" + str(stc) + "-srnds" + str(rounds) + "-x" + str(processors) + "/"
    
        output_file = output_subdirectory + "output_alignment.fasta"
        
        output_subset_file = output_subdirectory + "output_subset_alignment.fasta"
        
        alg_name = "UPP-fp5" + "K" + str(K) 

        if first_run:
            f0 = open(output_file, 'r')
            f1 = open(output_subset_file, 'w')
            lines = f0.readlines()
            for i in range(len(lines)):
                line = lines[i]
                if len(line) > 0:
                    if line[0] == ">":
                        if line[1:].strip() in reference_record_ids:
                            f1.write(lines[i])
                            f1.write(lines[i+1])
            f1.close()
    
            os.system("script -c \"java -jar " + FastSP_folder + "FastSP.jar -r " + reference_file + " -e " + output_subset_file + " \" " + output_subdirectory + "reference_comparison_scores.txt")
            
                   
        terminal_file = output_subdirectory + "terminal_output.txt"
        f2 = open(terminal_file)
        for line in f2.readlines(): 
            if len(line.split("line 197")) > 1:
                time_taken = float(line.split("line 197")[1].split(" ")[-2])
        times[data][ alg_name ] = time_taken  
        
        comp_file = output_subdirectory + "reference_comparison_scores.txt"
        f3 = open(comp_file)
        spfn = 0
        spfp = 0
        tc_score = 0
        for line in f3.readlines(): 
            if len(line.split()) > 0:
                if line.split()[0] == "SP-Score":
                    sp_ref = round(float(line.split()[1]), 3)
                if line.split()[0] == "SPFN":
                    spfn = round(float(line.split()[1]), 3)
                if line.split()[0] == "SPFP":
                    spfp = round(float(line.split()[1]), 3)
                if line.split()[0] == "TC":
                    tc_score = round(float(line.split()[1]), 3)
        ref_comp[data][ alg_name ] = sp_ref
        sp_errors[data][ alg_name ] = (spfn + spfp)/2
        tc_scores[data][ alg_name ] = tc_score    
        sp_scores[data][alg_name] = 1 - spfn
        modeler_scores[data][alg_name]  = 1 - spfp




#print("")
#for data in datas:
	#print(data, num_seqs[data], ref_comp[data], times[data]) 
#    print(data, num_seqs[data], sp_errors[data], tc_scores[data], times[data]) 	
        
print("")
for data in datas:
	#print(data, num_seqs[data], ref_comp[data], times[data]) 
    #print(data, num_seqs[data], sp_errors[data], tc_scores[data], times[data]) 	
    print(data)
    print("num seqs:", num_seqs[data])
    print("max seq len:", seq_len_stats[data][1])
    print("avg seq len:", seq_len_stats[data][2])
    #print("avg seq len:", seq_len_stats[data][3])
    print("sp-errors:", sp_errors[data])
    print("tc score:", tc_scores[data])
    print("sp-score:", sp_scores[data])
    print("modeler score :", modeler_scores[data])
    print("time (s):", times[data])
    print("")
    
alg_to_sp = {}
alg_to_tc = {}
alg_to_time = {}
alg_to_sps = {}
alg_to_modeler = {}
for data in sp_errors:
    if data != "rhv":
        for alg in sp_errors[data]:
            if alg not in alg_to_sp:
                alg_to_sp[alg] = []
            alg_to_sp[alg].append(sp_errors[data][alg])
        
for data in tc_scores:
    if data != "rhv":
        for alg in tc_scores[data]:
            if alg not in alg_to_tc:
                alg_to_tc[alg] = []
            alg_to_tc[alg].append(tc_scores[data][alg])
        
for data in times:
    if data != "rhv":
        for alg in sp_scores[data]:
            if alg not in alg_to_time:
                alg_to_sps[alg] = []
            alg_to_sps[alg].append(sp_scores[data][alg])
        
for data in times:
    if data != "rhv":
        for alg in modeler_scores[data]:
            if alg not in alg_to_modeler:
                alg_to_modeler[alg] = []
            alg_to_modeler[alg].append(modeler_scores[data][alg])
            
for data in times:
    if data != "rhv":
        for alg in times[data]:
            if alg not in alg_to_time:
                alg_to_time[alg] = []
            alg_to_time[alg].append(times[data][alg])
        
for alg in alg_to_sp:
    print(alg)
    print("avg. sp error", np.mean(alg_to_sp[alg]))
    print("avg. tc score", np.mean(alg_to_tc[alg]))
    print("avg. sp score", np.mean(alg_to_sps[alg]))
    print("avg. modeler score", np.mean(alg_to_modeler[alg]))
    print("avg. time", np.mean(alg_to_time[alg]))
    print("")
    
    


avg_avg_seq_len = 0
avg_num_seqs = 0
for data in datas:
    if data != "rhv":
        avg_avg_seq_len += seq_len_stats[data][2] / (len(datas) - 1)
        avg_num_seqs += num_seqs[data] / (len(datas) - 1)

print("avg. avg. seq len", avg_avg_seq_len)
print("avg. num seqs", avg_num_seqs)
        
        
final_sp = {}
final_tc = {}
final_time = {}
final_sps = {}
final_modeler = {}
final_algs = ["UPP_bb", "UPP2_bb", "UPP-fp5K10"]
for alg in final_algs:
    final_sp[alg] = np.mean(alg_to_sp[alg])
    final_tc[alg] = np.mean(alg_to_tc[alg])   
    final_sps[alg] = np.mean(alg_to_sps[alg])
    final_modeler[alg] = np.mean(alg_to_modeler[alg])
    final_time[alg] = np.mean(alg_to_time[alg])   
    
print("avg sp-error : ", final_sp)
print("avg tc-score : ", final_tc)
print("avg time : ", final_time)
print("avg sp-score : ", final_sps)
print("avg modeler score : ", final_modeler)
