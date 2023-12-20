#! /usr/bin/python3
# -*- coding: Utf-8 -*-

# Copyright (C) 2021 Emmanuelle Lerat
 
# This file is part of the replicaTE suite, developed by Van Anthony Le.

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
import time
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sps
import scipy.stats as stats
from random import choice, randint, random
from variaTE import pick_nt, distrib_exponweib, distrib_truncnorm, distrib_truncexp, generate_Seq
from Bio.Seq import Seq
#from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Create parameters and TE copy sequence in fasta.")
parser.add_argument("-m", "--min", type=int, help ="Minimal length of TE sequences (detecting tools can have min length requirement).", default = 80)
parser.add_argument("-n", "--nest", help ="Disable nested TEs ", action="store_false")
parser.add_argument("-d", "--div", type=int, help ="Maximal sequence divergence.", default = 20)
parser.add_argument("-c", "--cop", type=int, help ="Use simulated TE copy numbers (type -2).", default = -1)
parser.add_argument("-p", "--pst", type=float, help ="Probability to be on + strand.", default = 50)
parser.add_argument("-T", "--tsd",help ="Include TSD sequences", action="store_true") #by default, the parameter is set to "False"
parser.add_argument("-o", "--out",help ="Output name.", default="simulated_TEs")
args = parser.parse_args()

## Exponentiated Weibull distribution for number of copy
def simulate_copy(df_stat_TEs, mean_Copy):
    #data_cop = sorted(df_stat_TEs.nb_Copy.tolist())
    data_cop = sorted(df_stat_TEs.copy_nb.tolist())
    a, c, fit_loc, fit_scale = stats.exponweib.fit(data_cop, loc= min(data_cop), scale=np.mean(data_cop))
    distrib_Copy = distrib_exponweib(data_cop, a, c, 1, mean_Copy, 'y', "TE_copy_number")
    plt.hist(data_cop, bins = 200)
    plt.xlim(0, max(data_cop))
    plt.title("TE copy number distribution")
    plt.savefig("attachment/observed_copy_nb_distribution.png", bbox_inches='tight')
    #plt.show()
    plt.close()
    
    return distrib_Copy

## Generate parameters for the simulated TE copies and write them in a tsv file
def generate_param(min_Len_TEs, max_Div, mean_copy, prob_strand, presence_tsd):
#def generate_param(min_Len_TEs, max_Div, mean_copy, prob_strand):
    #print("min_Len_TEs " + str(min_Len_TEs))#test
    #print("max_Div " + str(max_Div))#test
    #print("mean_copy " + str(mean_copy))#test
    #print("prob_strand " + str(prob_strand))#test
    #print("presence_tsd " + str(presence_tsd)) #test
    
    ## Load data
    df_stat_TEs = pd.read_csv('data/stat_TEs_tab.tsv', sep ='\t')
    df_stat_TEs.rename(columns = {'Unnamed: 0':'id'}, inplace = True)
    df_intergenic = pd.read_csv('data/intergenic_sim_tab.tsv', sep ='\t')
    
    ## Initalize parameters
    #TEs_origin = df_stat_TEs.seq_Ancestral.tolist()
    TEs_origin = df_stat_TEs.ancestral_seq.tolist()
    length_TEs_origin = np.array([len(seq) for seq in TEs_origin])
    
    # Select TE families with a minimal length
    selected_TEs = np.where(length_TEs_origin >= min_Len_TEs)
    #seq_TEs_ori = df_stat_TEs.seq_Ancestral.values[selected_TEs]
    seq_TEs_ori = df_stat_TEs.ancestral_seq.values[selected_TEs]
    id_TEs_origin = df_stat_TEs.id.values[selected_TEs]
    mean_seq = df_stat_TEs.mean_length.values[selected_TEs]
    length_TEs_origin = length_TEs_origin[selected_TEs]
    nb_Fam = len(length_TEs_origin)
       
    #print("nb_Fam after selection of families of lenght > 80 bp" + "\t" + str(nb_Fam)) #test
    
    ## Retrieve copy number - either using the real copy number (-1 option) of by simulating the copy number (-2 option)
    if mean_copy == -1:
        #copy_per_Fam = np.array(df_stat_TEs.nb_Copy.tolist())
        copy_per_Fam = np.array(df_stat_TEs.copy_nb.tolist())
        #print("copy_per_Fam" + "\t" + str(copy_per_Fam)) #test

    if mean_copy == -2:
        #mean_copy = df_stat_TEs.nb_Copy.sum()/nb_Fam
        mean_copy = df_stat_TEs.copy_nb.sum()/nb_Fam
        if mean_copy > 40:
            mean_copy = mean_copy/4
        if mean_copy < 2:
            mean_copy = 2
        copy_per_Fam = simulate_copy(df_stat_TEs, mean_copy)        # Simulate copy number
        copy_per_Fam = np.random.choice(copy_per_Fam, nb_Fam, replace=False)        #take the selected TE family

    ## Truncated normal distribution for sequence copy divergence
    while max_Div >100 and max_Div <5:
        max_Div = int(input("Error: this divergence is incorrect, enter a new one: "))
        
    distrib_Div = list(distrib_truncnorm(copy_per_Fam.sum(),max_Div/2, 4, 0, max_Div, "TE_copy_divergence"))
    #print("distrib_Div" + "\t" + str(distrib_Div)) #test
    	
    len_copy_Fam = []
    div_per_Copy = []
    strand_copy_list = []
    TEs_TSD = []
    btm = 0
    iterator = 0
    strand_choice = ['+','-']
    #tsd_choice = [8, 7, 5, 4, 3, 2, 0]
    
    while prob_strand < 0 or prob_strand > 100:
         prob_strand = int(input("Error: this probability is incorrect, enter a new one: "))
            
    prob_strand = prob_strand/100 #50/100 => 0.5 with the default value
    
    if presence_tsd == False: 
        #tsd_size=0
        my_tsd_size = [0,0,0,0,0]
    else: #randomly select a TSD size among the proposed list
        my_tsd_size = [0,2,4,6,8]

    ## Attribute randomly length and divergence per copy
    for i in range(nb_Fam):

        # In case no copy
        if copy_per_Fam[i] == 0:
            len_copy_Fam.append("NA")
            div_per_Copy.append("NA")
            strand_copy_list.append("NA")
            TEs_TSD.append("NA")
            #TEs_TSD="NA"
        else:
            if copy_per_Fam[i] < 4:
                mean = length_TEs_origin[i]
                #mean = length_TEs_origin[i]/2 #commande originale
                #print(str(length_TEs_origin[i]) + "\t" + str(mean)) #test
                #print(str(copy_per_Fam[i]) + "\t" + str(mean)) #test
            else:
                mean = mean_seq[i]
                #print("mean nb copy > 4 \t" + str(mean)) #test
            btm = iterator
            iterator = iterator+(copy_per_Fam[i])
            
            pool_Length = distrib_truncexp(900, min_Len_TEs, length_TEs_origin[i]*1.025, mean).tolist()      # length can be > than length ancestral (fix = 2.5%)
            #print("avant "+ str(pool_Length)) #test
            add_len = [length_TEs_origin[i]]*100
            pool_Length = np.concatenate((pool_Length,add_len),axis=None)
            #print("apres " + str(pool_Length)) #test
            #length_Copy = np.random.choice(pool_Length, copy_per_Fam[i], replace=False)         # random select length
            length_Copy = np.random.choice(pool_Length, copy_per_Fam[i], replace=True)         # random select length #test avec replace =vrai
            #print(str(length_Copy)) #test
            strand_per_copy = list(np.random.choice(strand_choice, copy_per_Fam[i], p=[prob_strand, 1-prob_strand]))        # random select strand
            
            TEs_TSD_per_copy = list(np.random.choice(my_tsd_size, 1, p=[0.2, 0.2, 0.2, 0.2, 0.2]))*copy_per_Fam[i]      ## Generate TSD
            
            ## In case where the max copy length is too long, new random choice is made
            timeout = time.time() + 60*2

            while max(length_Copy) > 15000:
                #print(str(max(length_Copy)))
                #length_Copy = np.random.choice(pool_Length, copy_per_Fam[i], replace=False)
                length_Copy = np.random.choice(pool_Length, copy_per_Fam[i], replace=True)#test
                if time.time() > timeout:
                    print(" Error: Copy length is too long")
                    sys.exit()

            len_copy_Fam.append(list(length_Copy))
            div_per_Copy.append(distrib_Div[btm:iterator])
            strand_copy_list.append(strand_per_copy)
            TEs_TSD.append(TEs_TSD_per_copy)

    ## Test
    #print(copy_per_Fam.sum())# Test
    #print(len(distrib_Div))# Test
    #print(len(len_copy_Fam))# Test

###############
    ## Attribute randomly an intergenic sequence into which the TE copy will be inserted
    len_Inter = df_intergenic.index

    all_index = []
    #print("len_Inter " + str(len_Inter)) #test
    #print("len_copy_Fam " + str(len_copy_Fam)) #test
    for i in len_copy_Fam: 
        index_Inter = []
        for j in i:
            try:
                index_Inter.append(np.random.choice(len_Inter, replace=True))
                #print ("len_Inter " + str(len(len_Inter))) #test
            except ValueError:
            	try:
            	    index_Inter.append(np.random.choice(len_Inter, replace=False)) #test pour voir si c'est mieux
            	except:
            	     pass   
            	     
            #if randint(1,100) > 90 and args.nest:        ## Allows nested TEs => cree le probleme => comprendre pourquoi
             #   len_Inter = [x for x in len_Inter if x not in index_Inter]
            
        #print(str(len_Inter)) #test
        #print(str(len(index_Inter))) #test
        all_index.append(index_Inter)
        #print(str(all_index)) #test
    
    ## Output
    df_TEs_param = pd.DataFrame({'id': id_TEs_origin,'copy_number': copy_per_Fam, 'length_per_copy': len_copy_Fam, 'div_per_copy': div_per_Copy, 'loc_inter_seq': all_index, 'length_ancestral': length_TEs_origin, 'strand':strand_copy_list, 'TSD':TEs_TSD, 'Seq_ori': seq_TEs_ori})
    df_TEs_param.to_csv("data/param_TEs_tab.tsv",sep ='\t')
    return df_TEs_param

## Generate sequence for each copy
def generate_TEs_Seq(df_TEs_param):
    id_tot_TEs = []
    seq_tot_TEs = []
    test = []
    ## Loading bar
    print("  Generate simulated TE sequences")
    with click.progressbar (range(len(list(df_TEs_param.iterrows())))) as bar:
        for row in df_TEs_param.iterrows():
        
            bar.update(1)
            #if row[1].nb_Copy != 0:
            if row[1].copy_number != 0:
                Seq_ori = row[1].Seq_ori
                #len_origin = row[1].len_origin
                len_origin = row[1].length_ancestral
                gc_content = GC(Seq_ori)
                ## generate seq with param
                #for i in range(row[1].nb_Copy):
                for i in range(row[1].copy_number):
                    
                    #len_Copy = row[1].len_per_Copy[i]
                    len_Copy = row[1].length_per_copy[i]
                    id_Copy_fam = row[1].id+'_'+str(i+1) 
                    id_tot_TEs.append(id_Copy_fam)
                    if len_Copy  >= len_origin:
                        add_len_Seq = len_Copy - len_origin
                        add_seq = generate_Seq(add_len_Seq, gc_content) 
                        my_Seq = Seq_ori+add_seq

                    else:
                        pos_rand = randint(0, (len_origin - len_Copy))
                        my_Seq = Seq_ori[pos_rand:pos_rand+len_Copy]

                    ## Apply the divergence parser.add_argument("--out",help ="Name of the output.", default="output_clean_seq")
                    Divseq=[]
                    Pos=0
                    for nt in my_Seq:
                    
                        Pos+=1
                        #Div = row[1].div_per_Copy[i]/100
                        Div = row[1].div_per_copy[i]/100
                        Seed = random()
                        if Seed > Div :         # If the random number is above the diversity threshold, the nucleotide is picked up randomly
                            act=nt
                        else:
                            act=pick_nt(0.5)
                        Divseq.append(act)

                    new_seq = "".join(Divseq)       # Replace the not diverged sequence by the new one
                    
                    tsd_size= row[1].TSD[i]
                    if tsd_size == 2:
                        #TEs_TSD=str("AT")
                        new_seq="AT"+new_seq+"AT"
                        seq_tot_TEs.append(new_seq)
                    elif tsd_size == 4 :
                        #TEs_TSD=str("ATAT")
                        new_seq="ATAT"+new_seq+"ATAT"
                        seq_tot_TEs.append(new_seq)
                    elif tsd_size == 6 :
                        #TEs_TSD=str("ATATAT")
                        new_seq="ATATAT"+new_seq+"ATATATAT"
                        seq_tot_TEs.append(new_seq)
                    elif tsd_size == 8 :
                        #TEs_TSD=str("ATATATAT")
                        new_seq="ATATATAT"+new_seq+"ATATATAT"
                        seq_tot_TEs.append(new_seq)
                    elif tsd_size == 0:
                        seq_tot_TEs.append(new_seq)

    print("  Number of simulated TE copies: "+str(len(id_tot_TEs)))
    print("  Done."+'\n')
    return (id_tot_TEs, seq_tot_TEs)

def main():
    
    print('***************STEP2 generaTE: generate simulated TE sequences***************'+'\n')
    df_TEs_param = generate_param(args.min, args.div, args.cop, args.pst, args.tsd)
    #df_TEs_param = generate_param(args.min, args.div, args.cop, args.pst)
    data = generate_TEs_Seq(df_TEs_param)
    id_copy_TEs = data[0]
    seq_copy_TEs = data[1]
        
    ## Write clean_seq in fasta format
    records = []
    for seq in range(len(seq_copy_TEs)):
        records.append(SeqRecord(Seq(seq_copy_TEs[seq], generic_dna), id_copy_TEs[seq]+' '+id_copy_TEs[seq][:id_copy_TEs[seq].find('_')], description =''))

    with open("data/"+args.out+".fasta", "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")


if __name__ == "__main__":
    main()
