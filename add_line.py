#! /home/ajan/anaconda/bin/python
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 12:55:40 2016

@author: ajan
"""

import glob				# used to list directory contents
import os		
import numpy

data_files=glob.glob("*.py")

target_folder="./pyshape_edited/"
os.makedirs(target_folder)

line_to_add="#! /home/ajan/anaconda/bin/python"

print data_files[0]

for files in data_files:
    with open(target_folder+files,'w') as outfile:
        outfile.write(line_to_add)
        outfile.write("\n")
        with open(files) as infile:
            code=infile.read().splitlines(True)
            if code[0][0]=="#":
                outfile.writelines(code[1:]) 
            else:
                print code[0][0]
                outfile.writelines(code)
             
#             

