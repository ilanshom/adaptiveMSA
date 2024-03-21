#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 22:42:37 2023

@author: kmazooji
"""
#import Bio
#from Bio import Phylo
#from Bio import SeqIO
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord
import os


# For Homfam

output_directory = '/home/mazooji2/UPP-datasets-runs/compare_scores/'                     
input_directory =  '/home/mazooji2/UPP-datasets/homfam/'
tree_input_directory = '/home/mazooji2/UPP-datasets-runs/homfam/'

#datas = ["zf-CCHH", "tRNA-synt_2b", "aat_test-only", "p450_test-only"]
#datas = os.listdir(input_directory)
#datas = ["adh", "PDZ", "p450"]
#datas = os.listdir(input_directory)
datas = ["adh"] #, "PDZ", "p450"]
data_type = "amino"

#K = 25
#Ks = [i for i in range(5, 16)]
Ks = [5, 10, 15, 20]
backbone = 1000
#ktc = 50

processors = 24

exact = True

for data in datas:
    for K in Ks:
        backbone_tree_file = tree_input_directory + data +'-B1000-x24/output_pasta.fasttree'
        backbone_alignment_file = tree_input_directory + data +'-B1000-x24/output_pasta.fasta'
    
        data_file = input_directory + data + '/model/initial.fas'
    
        #FastSP_folder = '/home/mazooji2/FastSP-code/FastSP/'
    
        output_subdirectory = output_directory + data + "-cs-B" + str(backbone) + "-K" + str(K) + "/"
        #output_subdirectory = output_directory + data + "-fp3-B" + str(backbone) + "-K" + str(K) + "-ktc" + str(ktc) + "-stc" + str(stc) + "-srnds" + str(rounds) + "-iscr-x" + str(processors) + "/"
    
        os.system('mkdir ' + output_subdirectory)
    
        output_file = output_subdirectory + "output_alignment.fasta"
        
        #os.system("script -c  \"run_upp.py -A 10 -B " +  str(backbone) + " -M -1 -m " + data_type + " -s " + data_file + " --kmer " + str(K) + " --exact " + " -d " + output_subdirectory + " -x " +  str(processors) + " \" " + output_subdirectory + "terminal_output.txt")
        os.system("script -c  \"run_upp.py -A 10 -B " +  str(backbone) + " -M -1 -m " + data_type + " -s " + data_file + " --kmer " + str(K) +  " --exact " + " -d " + output_subdirectory + " -x " +  str(processors) + " -t " + backbone_tree_file  + " -a " + backbone_alignment_file + " --ignore-overlap \" " + output_subdirectory + "terminal_output.txt")

        #os.system("script -c  \"run_upp.py -A 10 -B " +  str(backbone) + " -M -1 -m " + data_type + " -s " + data_file + " --kmer " + str(K)  + " --exact -d " + output_subdirectory  + " \" " + output_subdirectory + "terminal_output.txt")

        #os.system("script -c  \"run_upp.py -A 10 -B " +  str(backbone) + " -M -1 -m rna -s " + data_file + " --kmer " + str(K) + " --kmerstc " + str(ktc) + " --setstc " + str(stc) + " --srnds " + str(rounds) + " -d " + output_subdirectory + " -x " +  str(processors) + " -t " + backbone_tree_file  + " -a " + backbone_alignment_file + " --ignore-overlap \" " + output_subdirectory + "terminal_output.txt")
        if not os.path.exists(output_subdirectory + "terminal_output_backup.txt"):
            os.system("cp " + output_subdirectory + "terminal_output.txt " + output_subdirectory + "terminal_output_backup.txt")
  