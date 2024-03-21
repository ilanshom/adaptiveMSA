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

output_directory = '/home/mazooji2/UPP-datasets-runs/indelible/'                     
input_directory =  '/home/mazooji2/UPP-datasets/indelible/'
#tree_input_directory = '/home/mazooji2/UPP-datasets-runs/16S/'

#datas = ["zf-CCHH", "tRNA-synt_2b", "aat_test-only", "p450_test-only"]
#datas = os.listdir(input_directory)
#datas = ["10000M2", "10000M3", "10000M4"]
datas = ["10000M4"]
samples = ['0']


for data in datas:
    for sample in samples:
        backbone = 1000
        processors = 24
    
    
        #backbone_tree_file = tree_input_directory + data +'-B1000-x24/output_pasta.fasttree'
        #backbone_alignment_file = tree_input_directory + data +'-B1000-x24/output_pasta.fasta'
    
        data_file = input_directory + data + '/' + sample + '/model/initial.fas'

        if data == "10000M4":
            data_file = input_directory + data + '/' + sample + '/model/true.fasta'

        #FastSP_folder = '/home/mazooji2/FastSP-code/FastSP/'
    
        output_subdirectory = output_directory + data + "-UPP2-early_stop" + "-x" + str(processors) + "/"
        #output_subdirectory = output_directory + data + "-fp3-B" + str(backbone) + "-K" + str(K) + "-ktc" + str(ktc) + "-stc" + str(stc) + "-srnds" + str(rounds) + "-iscr-x" + str(processors) + "/"
    
        os.system('mkdir ' + output_subdirectory)
    
        output_file = output_subdirectory + "output_alignment.fasta"
    
        os.system("script -c  \"run_upp.py -A 10 -B " +  str(backbone) + " -M -1 --molecule dna -s " + data_file + " --decomp_only True --bitscore_adjust True --hier_upp True --early_stop  True " +  " -d " + output_subdirectory + " -x " +  str(processors) + " \" " + output_subdirectory + "terminal_output.txt")
   
        #os.system("script -c  \"run_upp.py -A 10 -B " +  str(backbone) + " -M -1 -m amino -s " + data_file + " -kmer " + str(K) + " -kmerstc " + str(ktc) + " -setstc " + str(stc) + " -srnds " + str(rounds) + " -d " + output_subdirectory + " -x " +  str(processors) + " -iscr \" " + output_subdirectory + "terminal_output.txt")
    
        #os.system("script -c  \"run_upp.py -A 10 -B " +  str(backbone) + " -M -1 -m amino -s " + data_file + " -kmer " + str(K) + " -kmerstc " + str(ktc) + " -setstc " + str(stc) + " -srnds " + str(rounds) + " -d " + output_subdirectory + " -x " +  str(processors) + " -t " + backbone_tree_file  + " -a " + backbone_alignment_file + " --ignore-overlap \" " + output_subdirectory + "terminal_output.txt")
    
    
    
    
        # output_subdirectory = output_directory + data + "-fp3-B" + str(backbone) + "-K" + str(K) + "-ktc" + str(ktc) + "-stc" + str(stc) + "-srnds" + str(srnds) + "-x" + str(processors)  + "/"
        # os.system('mkdir ' + output_subdirectory)
    
        # output_file = output_subdirectory + "output_alignment.fasta"
    
        # os.system("script -c  \"run_upp.py -A 10 -B " +  str(backbone) + " -M -1 -m amino -s " + data_file + " -kmer " + str(K) + " -kmerstc " + str(ktc) + " -setstc " + str(stc) + " -srnds " + str(srnds) + " -d " + output_subdirectory + " -x " +  str(processors) + " \" " + output_subdirectory + "terminal_output.txt")
    
        # os.system("muscle -spscore " + output_subdirectory + "output_alignment.fasta -log " + output_subdirectory + "output_spscore.txt")
    



