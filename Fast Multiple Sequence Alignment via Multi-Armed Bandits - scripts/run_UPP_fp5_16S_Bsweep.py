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

output_directory = '/home/mazooji2/UPP-datasets-runs/16S/Bsweep/'                     
input_directory =  '/home/mazooji2/UPP-datasets/'

#datas = ["zf-CCHH", "tRNA-synt_2b", "aat_test-only", "p450_test-only"]
#datas = os.listdir(input_directory)
datas = ["16S.B.ALL"]

#K = 25
#Ks = [5, 10, 15, 20, 25]
K = 20
backbone = 1000
stc = 10    
processors = 24

frac_Bs = [0.1, 0.2, 0.3]
Rs = [2, 3, 4]

for data in datas:
    for frac_B in frac_Bs:
        for rounds in Rs:  
            data_file = input_directory + data + '/R0/cleaned.alignment.TtoU.fasta'
        
            #FastSP_folder = '/home/mazooji2/FastSP-code/FastSP/'
        
            output_subdirectory = output_directory + data + "-fp5-B" + str(backbone) + "-K" + str(K) + "-frac-" + str(frac_B) + "-stc" + str(stc) + "-srnds" + str(rounds) + "-x" + str(processors) + "/"
            #output_subdirectory = output_directory + data + "-fp3-B" + str(backbone) + "-K" + str(K) + "-ktc" + str(ktc) + "-stc" + str(stc) + "-srnds" + str(rounds) + "-iscr-x" + str(processors) + "/"
        
            os.system('mkdir ' + output_subdirectory)
        
            output_file = output_subdirectory + "output_alignment.fasta"
            
            os.system("script -c  \"run_upp.py -A 10 -B " +  str(backbone) + " -M -1 -m rna -s " + data_file + " --kmer " + str(K) + " --usefracbatch --batchsizefrac " + str(frac_B) + " --minbatchsize 50 "  + " --setstc " + str(stc) + " --srnds " + str(rounds) + " -d " + output_subdirectory + " -x " +  str(processors) + " \" " + output_subdirectory + "terminal_output.txt")
        
            #os.system("script -c  \"run_upp.py -A 10 -B " +  str(backbone) + " -M -1 -m rna -s " + data_file + " --kmer " + str(K) + " --kmerstc " + str(ktc) + " --setstc " + str(stc) + " --srnds " + str(rounds) + " -d " + output_subdirectory + " -x " +  str(processors) + " -t " + backbone_tree_file  + " -a " + backbone_alignment_file + " --ignore-overlap \" " + output_subdirectory + "terminal_output.txt")
            if not os.path.exists(output_subdirectory + "terminal_output_backup.txt"):
                os.system("cp " + output_subdirectory + "terminal_output.txt " + output_subdirectory + "terminal_output_backup.txt")
            
        
            # output_subdirectory = output_directory + data + "-fp3-B" + str(backbone) + "-K" + str(K) + "-ktc" + str(ktc) + "-stc" + str(stc) + "-srnds" + str(srnds) + "-x" + str(processors)  + "/"
            # os.system('mkdir ' + output_subdirectory)
        
            # output_file = output_subdirectory + "output_alignment.fasta"
        
            # os.system("script -c  \"run_upp.py -A 10 -B " +  str(backbone) + " -M -1 -m amino -s " + data_file + " -kmer " + str(K) + " -kmerstc " + str(ktc) + " -setstc " + str(stc) + " -srnds " + str(srnds) + " -d " + output_subdirectory + " -x " +  str(processors) + " \" " + output_subdirectory + "terminal_output.txt")
        
            # os.system("muscle -spscore " + output_subdirectory + "output_alignment.fasta -log " + output_subdirectory + "output_spscore.txt")      



