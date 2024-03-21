#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 13 16:04:10 2024

@author: kmazooji
"""
import os
import numpy as np
                   
# input_file =  '/home/mazooji2/UPP-datasets/16S.3/R0/cleaned.alignment.fasta'
# output_file = '/home/mazooji2/UPP-datasets/16S.3/R0/cleaned.alignment.TtoU.fasta'

#base = '/home/mazooji2/UPP-datasets/16S.3/R0/'
#base = '/home/mazooji2/UPP-datasets/16S.T/R0/'
base = '/home/mazooji2/UPP-datasets/16S.B.ALL/R0/'
input_file = base + 'cleaned.alignment.fasta'
output_file = base + 'cleaned.alignment.TtoU.fasta'

in_f = open(input_file, 'r')
o_f = open(output_file, 'w')

for line in in_f.readlines():
    if len(line) > 0:
        if line[0] == ">":
            o_f.write(line)
        if line[0] != ">":
            new_line = line.replace('T', 'U')
            o_f.write(new_line)
    else:
        o_f.write(line)

in_f.close()
o_f.close()
            
            





