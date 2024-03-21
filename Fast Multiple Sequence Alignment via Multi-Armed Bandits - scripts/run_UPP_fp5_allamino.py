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


# Rose

output_directory = '/home/mazooji2/UPP-datasets-runs/final_experiments_amino/'                     
input_directory =  '/home/mazooji2/UPP-datasets/ROSE/'


#datas = ["zf-CCHH", "tRNA-synt_2b", "aat_test-only", "p450_test-only"]
#datas = os.listdir(input_directory)
datas = ["1000M1", "1000S1", "1000L1"]

#K = 25
#Ks = [i for i in range(16, 21)]
Ks = [15, 20]

backbones = [100, 500] #500
frac_B = 0.2
rounds = 3
stc = 10    
processors = 24

sample = "R0"

for backbone in backbones:
    for data in datas:
        for K in Ks:
            # backbone_tree_file = tree_input_directory + data +'-B1000-x24/output_pasta.fasttree'
            # backbone_alignment_file = tree_input_directory + data +'-B1000-x24/output_pasta.fasta'
        
            data_file = input_directory + data + '/' + sample + '/rose.aln.true.fasta'
        
            #FastSP_folder = '/home/mazooji2/FastSP-code/FastSP/'
        
            #output_subdirectory = output_directory + data + "-fp5-B" + str(backbone) + "-K" + str(K) + "-ktc" + str(ktc) + "-stc" + str(stc) + "-srnds" + str(rounds) + "-x" + str(processors) + "/"        
    
        
            output_subdirectory = output_directory + data + "-fp5-B" + str(backbone) + "-K" + str(K) + "-frac-" + str(frac_B) + "-stc" + str(stc) + "-srnds" + str(rounds) + "-x" + str(processors) + "/"
    
    
    
            os.system('mkdir ' + output_subdirectory)
        
            output_file = output_subdirectory + "output_alignment.fasta"
            
            #os.system("script -c  \"run_upp.py -A 10 -B " +  str(backbone) + " -M -1 -m amino -s " + data_file + " --kmer " + str(K) + " --kmerstc " + str(ktc) + " --setstc " + str(stc) + " --srnds " + str(rounds) + " -d " + output_subdirectory + " -x " +  str(processors) + " \" " + output_subdirectory + "terminal_output.txt")
    
    
        
            os.system("script -c  \"run_upp.py -A 10 -B " +  str(backbone) + " -M -1 -m amino -s " + data_file + " --kmer " + str(K) + " --usefracbatch --batchsizefrac " + str(frac_B) + " --minbatchsize 50 "  + " --setstc " + str(stc) + " --srnds " + str(rounds) + " -d " + output_subdirectory + " -x " +  str(processors) + " \" " + output_subdirectory + "terminal_output.txt")
        
    
    
            #os.system("script -c  \"run_upp.py -A 10 -B " +  str(backbone) + " -M -1 -m rna -s " + data_file + " --kmer " + str(K) + " --kmerstc " + str(ktc) + " --setstc " + str(stc) + " --srnds " + str(rounds) + " -d " + output_subdirectory + " -x " +  str(processors) + " -t " + backbone_tree_file  + " -a " + backbone_alignment_file + " --ignore-overlap \" " + output_subdirectory + "terminal_output.txt")
            if not os.path.exists(output_subdirectory + "terminal_output_backup.txt"):
                os.system("cp " + output_subdirectory + "terminal_output.txt " + output_subdirectory + "terminal_output_backup.txt")
            
        
            # output_subdirectory = output_directory + data + "-fp3-B" + str(backbone) + "-K" + str(K) + "-ktc" + str(ktc) + "-stc" + str(stc) + "-srnds" + str(srnds) + "-x" + str(processors)  + "/"
            # os.system('mkdir ' + output_subdirectory)
        
            # output_file = output_subdirectory + "output_alignment.fasta"
        
            # os.system("script -c  \"run_upp.py -A 10 -B " +  str(backbone) + " -M -1 -m amino -s " + data_file + " -kmer " + str(K) + " -kmerstc " + str(ktc) + " -setstc " + str(stc) + " -srnds " + str(srnds) + " -d " + output_subdirectory + " -x " +  str(processors) + " \" " + output_subdirectory + "terminal_output.txt")
        
            # os.system("muscle -spscore " + output_subdirectory + "output_alignment.fasta -log " + output_subdirectory + "output_spscore.txt")      





# homfam

input_directory =  '/home/mazooji2/UPP-datasets/homfam/'
tree_input_directory = '/home/mazooji2/UPP-datasets-runs/homfam/'

#datas = ["zf-CCHH", "tRNA-synt_2b", "aat_test-only", "p450_test-only"]
#datas = os.listdir(input_directory)
#datas = ["adh", "PDZ", "p450"]
datas = os.listdir(input_directory)

#K = 25
#Ks = [i for i in range(5, 16)]
Ks = [5, 10, 15, 20]
backbone = 1000
#ktc = 50

# frac_B = 0.2
# rounds = 3
# stc = 10    
# processors = 24

# exact = False

for data in datas:
    for K in Ks:
        backbone_tree_file = tree_input_directory + data +'-B1000-x24/output_pasta.fasttree'
        backbone_alignment_file = tree_input_directory + data +'-B1000-x24/output_pasta.fasta'
    
        data_file = input_directory + data + '/model/initial.fas'
    
        #FastSP_folder = '/home/mazooji2/FastSP-code/FastSP/'
    
        #output_subdirectory = output_directory + data + "-fp5-B" + str(backbone) + "-K" + str(K) + "-ktc" + str(ktc) + "-stc" + str(stc) + "-srnds" + str(rounds) + "-x" + str(processors) + "/"        

    
        output_subdirectory = output_directory + data + "-fp5-B" + str(backbone) + "-K" + str(K) + "-frac-" + str(frac_B) + "-stc" + str(stc) + "-srnds" + str(rounds) + "-x" + str(processors) + "/"



        os.system('mkdir ' + output_subdirectory)
    
        output_file = output_subdirectory + "output_alignment.fasta"
        
        #os.system("script -c  \"run_upp.py -A 10 -B " +  str(backbone) + " -M -1 -m amino -s " + data_file + " --kmer " + str(K) + " --kmerstc " + str(ktc) + " --setstc " + str(stc) + " --srnds " + str(rounds) + " -d " + output_subdirectory + " -x " +  str(processors) + " \" " + output_subdirectory + "terminal_output.txt")


    
        #os.system("script -c  \"run_upp.py -A 10 -B " +  str(backbone) + " -M -1 -m amino -s " + data_file + " --kmer " + str(K) + " --usefracbatch --batchsizefrac " + str(frac_B) + " --minbatchsize 50 "  + " --setstc " + str(stc) + " --srnds " + str(rounds) + " -d " + output_subdirectory + " -x " +  str(processors) + " \" " + output_subdirectory + "terminal_output.txt")
        os.system("script -c  \"run_upp.py -A 10 -B " +  str(backbone) + " -M -1 -m amino -s " + data_file + " --kmer " + str(K) + " --usefracbatch --batchsizefrac " + str(frac_B) + " --minbatchsize 50 "  + " --setstc " + str(stc) + " --srnds " + str(rounds) + " -d " + output_subdirectory + " -x " +  str(processors) + " -t " + backbone_tree_file  + " -a " + backbone_alignment_file + " --ignore-overlap \" " + output_subdirectory + "terminal_output.txt")



        #os.system("script -c  \"run_upp.py -A 10 -B " +  str(backbone) + " -M -1 -m rna -s " + data_file + " --kmer " + str(K) + " --kmerstc " + str(ktc) + " --setstc " + str(stc) + " --srnds " + str(rounds) + " -d " + output_subdirectory + " -x " +  str(processors) + " -t " + backbone_tree_file  + " -a " + backbone_alignment_file + " --ignore-overlap \" " + output_subdirectory + "terminal_output.txt")
        if not os.path.exists(output_subdirectory + "terminal_output_backup.txt"):
            os.system("cp " + output_subdirectory + "terminal_output.txt " + output_subdirectory + "terminal_output_backup.txt")
        
    
        # output_subdirectory = output_directory + data + "-fp3-B" + str(backbone) + "-K" + str(K) + "-ktc" + str(ktc) + "-stc" + str(stc) + "-srnds" + str(srnds) + "-x" + str(processors)  + "/"
        # os.system('mkdir ' + output_subdirectory)
    
        # output_file = output_subdirectory + "output_alignment.fasta"
    
        # os.system("script -c  \"run_upp.py -A 10 -B " +  str(backbone) + " -M -1 -m amino -s " + data_file + " -kmer " + str(K) + " -kmerstc " + str(ktc) + " -setstc " + str(stc) + " -srnds " + str(srnds) + " -d " + output_subdirectory + " -x " +  str(processors) + " \" " + output_subdirectory + "terminal_output.txt")
    
        # os.system("muscle -spscore " + output_subdirectory + "output_alignment.fasta -log " + output_subdirectory + "output_spscore.txt")      
