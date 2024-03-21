#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 01:58:43 2024

@author: kmazooji
"""
import os



# for 16S

output_directory = '/home/mazooji2/UPP-datasets-runs/final_experiments/'                     
input_directory =  '/home/mazooji2/UPP-datasets/'

#datas = ["zf-CCHH", "tRNA-synt_2b", "aat_test-only", "p450_test-only"]
#datas = os.listdir(input_directory)
datas = ["16S.T", "16S.3", "16S.B.ALL"]

#K = 25
#Ks = [5, 10, 15, 20, 25]
K = 20
backbone = 1000
stc = 10    
processors = 24
frac_B = 0.2
rounds = 3


for data in datas:
    data_file = input_directory + data + '/R0/cleaned.alignment.TtoU.fasta'

    output_subdirectory = output_directory + data + "-fp5-B" + str(backbone) + "-K" + str(K) + "-frac-" + str(frac_B) + "-stc" + str(stc) + "-srnds" + str(rounds) + "-x" + str(processors) + "/"

    os.system('mkdir ' + output_subdirectory)

    output_file = output_subdirectory + "output_alignment.fasta"
    
    os.system("script -c  \"run_upp.py -A 10 -B " +  str(backbone) + " -M -1 -m rna -s " + data_file + " --kmer " + str(K) + " --usefracbatch --batchsizefrac " + str(frac_B) + " --minbatchsize 50 "  + " --setstc " + str(stc) + " --srnds " + str(rounds) + " -d " + output_subdirectory + " -x " +  str(processors) + " \" " + output_subdirectory + "terminal_output.txt")

    if not os.path.exists(output_subdirectory + "terminal_output_backup.txt"):
        os.system("cp " + output_subdirectory + "terminal_output.txt " + output_subdirectory + "terminal_output_backup.txt")


input_directory =  '/home/mazooji2/UPP-datasets/indelible/'
datas = ["10000M2", "10000M3", "10000M4"]
sample = "0"

for data in datas:

    data_file = input_directory + data + '/' + sample + '/model/initial.fas'
    
    if data == "10000M4":
        data_file = input_directory + data + '/' + sample + '/model/true.fasta'


    output_subdirectory = output_directory + data + "-fp5-B" + str(backbone) + "-K" + str(K) + "-frac-" + str(frac_B) + "-stc" + str(stc) + "-srnds" + str(rounds) + "-x" + str(processors) + "/"
    
    os.system('mkdir ' + output_subdirectory)

    output_file = output_subdirectory + "output_alignment.fasta"
    
    os.system("script -c  \"run_upp.py -A 10 -B " +  str(backbone) + " -M -1 -m dna -s " + data_file + " --kmer " + str(K) + " --usefracbatch --batchsizefrac " + str(frac_B) + " --minbatchsize 50 "  + " --setstc " + str(stc) + " --srnds " + str(rounds) + " -d " + output_subdirectory + " -x " +  str(processors) + " \" " + output_subdirectory + "terminal_output.txt")

    if not os.path.exists(output_subdirectory + "terminal_output_backup.txt"):
        os.system("cp " + output_subdirectory + "terminal_output.txt " + output_subdirectory + "terminal_output_backup.txt")
    
    
input_directory =  '/home/mazooji2/UPP-datasets/rnasim/'
datas = ["10000", "50000", "100000"]

for data in datas:

    data_file = input_directory + data + '/1/model/true.fasta'

    output_subdirectory = output_directory + data + "-fp5-B" + str(backbone) + "-K" + str(K) + "-frac-" + str(frac_B) + "-stc" + str(stc) + "-srnds" + str(rounds) + "-x" + str(processors) + "/"
    
    os.system('mkdir ' + output_subdirectory)

    output_file = output_subdirectory + "output_alignment.fasta"
    
    os.system("script -c  \"run_upp.py -A 10 -B " +  str(backbone) + " -M -1 -m rna -s " + data_file + " --kmer " + str(K) + " --usefracbatch --batchsizefrac " + str(frac_B) + " --minbatchsize 50 "  + " --setstc " + str(stc) + " --srnds " + str(rounds) + " -d " + output_subdirectory + " -x " +  str(processors) + " \" " + output_subdirectory + "terminal_output.txt")

    if not os.path.exists(output_subdirectory + "terminal_output_backup.txt"):
        os.system("cp " + output_subdirectory + "terminal_output.txt " + output_subdirectory + "terminal_output_backup.txt")
