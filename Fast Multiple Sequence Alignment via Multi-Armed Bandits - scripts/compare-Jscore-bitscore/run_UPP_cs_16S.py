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


# output_directory = '/home/kmazooji/UPP-datasets-runs/compare_scores/'                     
# input_directory =  '/home/kmazooji/UPP-datasets/ROSE/'
#tree_input_directory = '/home/kmazooji/UPP-datasets-runs/homfam/'

output_directory = '/home/mazooji2/UPP-datasets-runs/compare_scores/'                     
input_directory =  '/home/mazooji2/UPP-datasets/'

#datas = ["zf-CCHH", "tRNA-synt_2b", "aat_test-only", "p450_test-only"]
datas = [ "16S.T"] #os.listdir(input_directory)
data_type = "rna"

Ks = [10] #[20, 15, 10]
backbone = 1000

for data in datas:
    for K in Ks:
    
        
        # backbone_tree_file = tree_input_directory + data +'-B1000-x24/output_pasta.fasttree'
        # backbone_alignment_file = tree_input_directory + data +'-B1000-x24/output_pasta.fasta'
        
        data_file = input_directory + data + '/R0/cleaned.alignment.TtoU.fasta'
    
    
        #output_subdirectory = output_directory + data + "-fp3-B" + str(backbone) + "-K" + str(K) + "-ktc" + str(ktc) + "-stc" + str(stc) + "-srnds" + str(rounds) + "-iscr-x" + str(processors) + "/"
    
        output_subdirectory = output_directory + data + "-cs-B" + str(backbone) + "-K" + str(K) + "/"

        #output_subdirectory = output_directory + data + "-cs-run2-B" + str(backbone) + "-K" + str(K) + "/"

    
    
        os.system('mkdir ' + output_subdirectory)
    
        output_file = output_subdirectory + "output_alignment.fasta"
    
        #os.system("script -c  \"run_upp.py -A 10 -B " +  str(backbone) + " -M -1 -m amino -s " + data_file + " -kmer " + str(K) + " -kmerstc " + str(ktc) + " -setstc " + str(stc) + " -srnds " + str(rounds) + " -d " + output_subdirectory + " -x " +  str(processors) + " -t " + backbone_tree_file  + " -a " + backbone_alignment_file + " -iscr --ignore-overlap \" " + output_subdirectory + "terminal_output.txt")
    
        os.system("script -c  \"run_upp.py -A 10 -B " +  str(backbone) + " -M -1 -m " + data_type + " -s " + data_file + " --kmer " + str(K)  + " --exact -d " + output_subdirectory  + " \" " + output_subdirectory + "terminal_output.txt")
    
        #os.system("script -c  \"run_upp.py -A 10 -B " +  str(backbone) + " -M -1 -m amino -s " + data_file + " -kmer " + str(K) + " -kmerstc " + str(ktc) + " -setstc " + str(stc) + " -srnds " + str(rounds) + " -d " + output_subdirectory + " -x " +  str(processors) + " \" " + output_subdirectory + "terminal_output.txt")
    
    
        #os.system("muscle -spscore " + output_subdirectory + "output_alignment.fasta -log " + output_subdirectory + "output_spscore.txt")


