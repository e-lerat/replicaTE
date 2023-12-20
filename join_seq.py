#! /usr/bin/python3
# -*- coding: Utf-8 -*-

# Copyright (C) 2021 Emmanuelle Lerat
 
# This file is part of replicaTE suite, developed by Van Anthony Le.

# replicaTE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# replicaTE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with replicaTE.  If not, see <http://www.gnu.org/licenses/>.

import click
import argparse 
import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np
import pandas as pd
import os
from re import sub
from variaTE import distrib_truncnorm, distrib_exponweib, simulate_Seq, generate_Seq
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC
import shlex, subprocess
import os

## Retrieve positions of genes and intergenic regions and concatenate all sequences
def concat_seq(df_inter_Seq, df_new_Gen):
    inter_Seq_list = df_inter_Seq.sequence.tolist()
    new_Gen_list = df_new_Gen.sequence.tolist()
    clean_Seq = ''
    Pos_Start_Gene = []
    Pos_End_Gene = []
    Pos_Start_intergenic = []
    Pos_End_intergenic = []
    print("         Final number of genes: ",len(new_Gen_list))
    
    for i in range(len(new_Gen_list)):
        
       # print(i)
        if i == 0:
            Pos_Start_intergenic.append(0)
            Pos_End_intergenic.append(len(inter_Seq_list[i])-1)

        else:
            Pos_Start_intergenic.append(Pos_End_Gene[-1]+1)
            Pos_End_intergenic.append(Pos_Start_intergenic[-1]+len(inter_Seq_list[i])-1)

        #print(new_Gen_list[i])#test
        Pos_Start_Gene.append(Pos_End_intergenic[-1]+1)
        #print (type(new_Gen_list[i]))#test
        if isinstance(new_Gen_list[i],str): #verify whether the cleaned gene has a sequence or not
            Pos_End_Gene.append(Pos_Start_Gene[-1]+len(new_Gen_list[i])-1)
        elif isinstance(new_Gen_list[i],float):
            Pos_End_Gene.append(Pos_Start_Gene[-1]-1)
        
        #print (str(Pos_Start_Gene)  + "\t" + str(Pos_End_Gene))#test
        #print(len(new_Gen_list[i]))#test
        
        clean_Seq = clean_Seq+inter_Seq_list[i]
        if isinstance(new_Gen_list[i],str): #verify whether the cleaned gene has a sequence or not
            clean_Seq = clean_Seq+new_Gen_list[i]
        elif isinstance(new_Gen_list[i],float):
            clean_Seq = clean_Seq

    clean_Seq = data_clean[0]

    ## Write clean_Seq (simulated genome with clean genes and simulated intergenic regions) in fasta format
    #if isinstance(seq_Id, list):
     #   seq_Id = ''.join(seq_Id)
        #print(seq_Id)
   # if len(seq_Desc) > 1:
    #    seq_Desc = organism
    #else:
     #   seq_Desc= ''.join(seq_Desc)
    clean_record = SeqRecord(Seq(clean_Seq, generic_dna),id = "toto")
    with open("test.fasta", "w") as output_handle:
        SeqIO.write(clean_record, output_handle, "fasta")
            
