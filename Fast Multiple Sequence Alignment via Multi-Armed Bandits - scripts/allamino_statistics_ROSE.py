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
output_directory =  '/home/mazooji2/UPP-datasets-runs/ROSE/'              
input_directory =  '/home/mazooji2/UPP-datasets/ROSE/' #'16S/'

FastSP_folder = '/home/mazooji2/FastSP-code/FastSP/'

first_run = True



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
datas  = ["1000M1", "1000S1", "1000L1"] #, "1000S1", "1000L1"] #, "16S.T", "16S.3"]
#datas  = ["16S.T", "16S.3"]
#Ks = [i for i in range(5, 21)]
Ks = [15, 20]
frac_B = 0.2
rounds = 3
stc = 10    
backbone = 500
processors = 24



sample = "R0"

for data in datas:

    reference_file = input_directory + data + '/' + sample + '/rose.aln.true.fasta'
    
    #seqs_file = input_directory + data + '/model/initial.fas'
    
    

    f0 = open(reference_file)
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

    # compare upp
    
    
    times[data] = {}
    ref_comp[data] = {}
    sp_errors[data] = {}
    tc_scores[data] = {}
    sp_scores[data] = {}
    modeler_scores[data] = {}
    
    
    output_subdirectory = output_directory + data + "-" + sample + "-B" + str(backbone) + "-x" + str(processors) + "/"

    output_file = output_subdirectory + "output_alignment.fasta"
    
    
    #print("script -c \"java -jar " + FastSP_folder + "FastSP.jar -r " + reference_file + " -e " + output_file + " \" " + output_subdirectory + "reference_comparison_scores.txt")

    #os.system("script -c \"java -jar " + FastSP_folder + "FastSP.jar -r " + reference_file + " -e " + output_file + " \" " + output_subdirectory + "reference_comparison_scores.txt")

    
    if first_run:
        os.system("script -c \"java -jar " + FastSP_folder + "FastSP.jar -r " + reference_file + " -e " + output_file + " \" " + output_subdirectory + "reference_comparison_scores.txt")

    terminal_file = output_subdirectory + "terminal_output.txt"
	
    f2 = open(terminal_file)
    for line in f2.readlines(): 
        if len(line.split("line 197")) > 1:
            time_taken = float(line.split("line 197")[1].split(" ")[-2])
    times[data]["UPP"] = time_taken   

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
    ref_comp[data]["UPP"] = sp_ref
    sp_errors[data]["UPP"] = (spfn + spfp)/2
    tc_scores[data]["UPP"] = tc_score
    sp_scores[data]["UPP"] = 1 - spfn
    modeler_scores[data]["UPP"]  = 1 - spfp
    
    
    # compare upp2


    output_subdirectory = output_directory + data + "-" +  str(sample) +  "-UPP2" + "-B" + str(backbone) + "-x" + str(processors) + "/"

    output_file = output_subdirectory + "output_alignment.fasta"
    
    
    #print("script -c \"java -jar " + FastSP_folder + "FastSP.jar -r " + reference_file + " -e " + output_file + " \" " + output_subdirectory + "reference_comparison_scores.txt")

    #os.system("script -c \"java -jar " + FastSP_folder + "FastSP.jar -r " + reference_file + " -e " + output_file + " \" " + output_subdirectory + "reference_comparison_scores.txt")

    
    run_name = "UPP2"
    
    
    if first_run:
        os.system("script -c \"java -jar " + FastSP_folder + "FastSP.jar -r " + reference_file + " -e " + output_file + " \" " + output_subdirectory + "reference_comparison_scores.txt")

        
    terminal_file = output_subdirectory + "terminal_output.txt"
    f2 = open(terminal_file)
    for line in f2.readlines(): 
        if len(line.split("line 197")) > 1:
            time_taken = float(line.split("line 197")[1].split(" ")[-2])
    times[data][run_name] = time_taken  
    
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
    ref_comp[data][run_name] = sp_ref
    sp_errors[data][run_name] = (spfn + spfp)/2
    tc_scores[data][run_name] = tc_score    
    sp_scores[data][run_name] = 1 - spfn
    modeler_scores[data][run_name]  = 1 - spfp


    for K in Ks:  
    
        # compare upp_fp5
        
        output_subdirectory = final_output_directory+ data + "-fp5-B" + str(backbone) + "-K" + str(K) + "-frac-" + str(frac_B) + "-stc" + str(stc) + "-srnds" + str(rounds) + "-x" + str(processors) + "/"
    
        output_file = output_subdirectory + "output_alignment.fasta"
        
        
        run_name = "UPP-fp5-exact-" + "K" + str(K) 
        
        
        if first_run:
            os.system("script -c \"java -jar " + FastSP_folder + "FastSP.jar -r " + reference_file + " -e " + output_file + " \" " + output_subdirectory + "reference_comparison_scores.txt")
            #returned_output2 = subprocess.check_output("script -c \"java -jar " + FastSP_folder + "FastSP.jar -r " + reference_file + " -e " + output_file + " \" " + output_subdirectory + "reference_comparison_scores.txt",  shell=True, text=True)
            
        terminal_file = output_subdirectory + "terminal_output.txt"
        f2 = open(terminal_file)
        for line in f2.readlines(): 
            if len(line.split("line 197")) > 1:
                time_taken = float(line.split("line 197")[1].split(" ")[-2])
        times[data][run_name] = time_taken  
        
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
        ref_comp[data][run_name] = sp_ref
        sp_errors[data][run_name] = (spfn + spfp)/2
        tc_scores[data][run_name] = tc_score   
        sp_scores[data][run_name] = 1 - spfn
        modeler_scores[data][run_name]  = 1 - spfp



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
    print("")
    #print("avg seq len:", seq_len_stats[data][3])
    print("sp-errors:", sp_errors[data])
    print("")
    print("tc score:", tc_scores[data])
    print("")
    print("time (s):", times[data])
    print("")
    print("sp score:", sp_scores[data])
    print("")
    print("modeler score:", modeler_scores[data])
    print("")
    print("")
        