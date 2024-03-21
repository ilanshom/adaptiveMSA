# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 00:31:44 2024

@author: kmazo
"""
import os
import subprocess
import numpy as np

output_directory = '/home/mazooji2/UPP-datasets-runs/final_experiments_mem_usage/UPP/'
output_directory = '/home/mazooji2/UPP-datasets-runs/final_experiments_mem_usage/UPP2/'    
output_directory = '/home/mazooji2/UPP-datasets-runs/final_experiments_mem_usage/UPP-fp5/'                     

datas = os.listdir(output_directory)

for data in datas:
    if data[-4:] == ".dat":
        
        mem_file = output_directory + data
    
        reference_record_ids = {}
        f1 = open(mem_file, 'r')
        lines1 = f1.readlines()
        print(lines1[0])
        highest_seen = 0
        time_seen = None
        for line in lines1[1:]:
            if len(line) > 0:
                tokens = line.split()
                mem_used = float(tokens[1]) 
                if mem_used > highest_seen:
                    highest_seen = mem_used 
                    time_seen = float(tokens[2]) 
        print(highest_seen)
        print(time_seen)
        print("")