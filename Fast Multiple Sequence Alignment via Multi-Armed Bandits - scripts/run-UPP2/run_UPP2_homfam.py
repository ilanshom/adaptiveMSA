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

output_directory = '/home/mazooji2/UPP-datasets-runs/homfam/'                     
input_directory =  '/home/mazooji2/UPP-datasets/homfam/'
tree_input_directory = '/home/mazooji2/UPP-datasets-runs/homfam/'

#datas = ["zf-CCHH", "tRNA-synt_2b", "aat_test-only", "p450_test-only"]
datas = os.listdir(input_directory)

for data in datas:
    backbone = 1000
    processors = 24

    data_file = input_directory + data + '/model/initial.fas'
    
    backbone_tree_file = tree_input_directory + data +'-B1000-x24/output_pasta.fasttree'
    backbone_alignment_file = tree_input_directory + data +'-B1000-x24/output_pasta.fasta'

    #FastSP_folder = '/home/mazooji2/FastSP-code/FastSP/'

    output_subdirectory = output_directory + data + "-UPP2-common_backbone-early_stop" + "-x" + str(processors) + "/"
    os.system('mkdir ' + output_subdirectory)

    output_file = output_subdirectory + "output_alignment.fasta"

    #os.system("script -c  \"run_upp.py -A 10 -B " +  str(backbone) + " -M -1 --molecule amino -s " + data_file + " --decomp_only True --bitscore_adjust True --hier_upp True --early_stop  True " +  " -d " + output_subdirectory + " -x " +  str(processors) + " \" " + output_subdirectory + "terminal_output.txt")
    os.system("script -c  \"run_upp.py " + " -M -1 --molecule amino -s " + data_file + " --decomp_only True --bitscore_adjust True --hier_upp True --early_stop  True " +  " -d " + output_subdirectory + " -x " +  str(processors) + " -t " + backbone_tree_file  + " -a " + backbone_alignment_file + " --ignore-overlap \" " +  output_subdirectory + "terminal_output.txt")


