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

output_directory = '/home/mazooji2/UPP-datasets-runs/16S/'                     
input_directory =  '/home/mazooji2/UPP-datasets/'

#datas = ["zf-CCHH", "tRNA-synt_2b", "aat_test-only", "p450_test-only"]
#datas = os.listdir(input_directory)
datas = ["16S.T", "16S.3", "16S.B.ALL"]

for data in datas:
    backbone = 1000
    processors = 24

    data_file = input_directory + data + '/R0/cleaned.alignment.TtoU.fasta'

    #FastSP_folder = '/home/mazooji2/FastSP-code/FastSP/'

    output_subdirectory = output_directory + data + "-rna-B" + str(backbone) + "-x" + str(processors) + "/"
    os.system('mkdir ' + output_subdirectory)

    output_file = output_subdirectory + "output_alignment.fasta"

    os.system("script -c  \"run_upp.py -A 10 -B " +  str(backbone) + " -M -1 -m rna -s " + data_file +  " -d " + output_subdirectory + " -x " +  str(processors) + " \" " + output_subdirectory + "terminal_output.txt")

    #os.system("muscle -spscore " + output_subdirectory + "output_alignment.fasta -log " + output_subdirectory + "output_spscore.txt")




