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

import os
import csv
import copy
import click
import argparse
import collections
import pandas as pd
import numpy as np
from Bio import SeqIO
from re import sub
from random import randint
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
from random import random, choice
from deleTE import concat_seq

parser = argparse.ArgumentParser(description="Insert TE sequences in the simulated genome")
parser.add_argument("-t", "--TEs", help ="fasta file of simulated TEs", default = "data/simulated_TEs.fasta")
parser.add_argument("-i", "--int", help ="Tsv file of simulated intergenic regions", default = "data/intergenic_sim_tab.tsv")
parser.add_argument("-g", "--gen", help ="Tsv file of genes", default = "data/gene_clean_tab.tsv")
parser.add_argument("-p", "--par", help ="Tsv file of parameters of simulated TEs", default = "data/param_TEs_tab.tsv")
parser.add_argument("-d", "--des", help ="get the description", default = "attachment/output_clean_seq.fasta")
parser.add_argument("-n", "--ind", help ="number of deleted genome (type int)", default = 1)
parser.add_argument("-m", "--met", help ="allows to have new insertions in the deleted genome (not present in the simulated genome)", action="store_true")
parser.add_argument("-o", "--out", help ="output file name", default="simulated_genome")

args = parser.parse_args()

## Insert TEs into the attributed intergenic region
def insert_TEs(generated_TES, df_Inter, df_gene, df_TEs):
    print("  Create the custom reference genome")
    ## load data
    records = list(SeqIO.parse(generated_TES, "fasta"))
    df_gene = df_gene.drop(columns=['Unnamed: 0'])
    df_Inter = df_Inter.drop(columns=['Unnamed: 0'])
    df_TEs = df_TEs.drop(columns=['Unnamed: 0'])
    df_TEs = df_TEs.dropna()


    insert = {}
    index = 0
    df_seq_inter = df_Inter['sequence']
    df_seq_length = df_Inter['length']
    df_inter = pd.concat([df_seq_length, df_seq_inter], axis=1)
    df_inter_copy = df_inter.copy(deep =True)
    df_Inter_inserted = pd.concat([df_seq_length, df_seq_inter], axis=1)
    df_Inter_inserted_masked = pd.concat([df_seq_length, df_seq_inter], axis=1)
    #ligne=0 #test
       
    for row in df_TEs.iterrows():
        loc_Inser = row[1].loc_inter_seq.split()
        div_copy = row[1].div_per_copy.split()
        tsd_copy = row[1].TSD.split()
        strand_copy = row[1].strand.split()
        iter = 0
        #ligne+=1 #test
        #print(ligne) #test
        #print(len(list(df_TEs.iterrows()))) #nb de familles

        #print("insert " + str(insert) + "loc_Inser "+ str(loc_Inser))#test
        #print("insert " + str(insert) + "strand "+ str(strand_copy))#test
        #print("insert " + str(insert) + "div "+ str(div_copy))#test
        for i in loc_Inser:
            ##print the nb of TE
            print(index) #test
            i = int(sub('[^0-9]', '', i)) #takes the number corresponding to a particular intergenic region
            strand_TE =sub('[^\+\-]', '', strand_copy[iter])
            div_TE = float(sub('[^0-9 \.]', '', div_copy[iter]))
            TSD_TE = int(sub('[^0-9 \.]', '', tsd_copy[iter]))
            seq_Inter = df_Inter_inserted.at[i,'sequence']
            #takes sequence of TE		
            seq_TE = str(records[index].seq)
            seq_TE_masked = 'N'*(len(seq_TE))
            ## Do the reverse complement
            if strand_TE =='-':
                seq_TE = Seq(seq_TE, generic_dna)
                seq_TE = str(seq_TE.reverse_complement())
            id_TE = records[index].id
            
            try:
                pos_insert = randint(TSD_TE, (len(seq_Inter)))
            except: 
                pos_insert = randint(TSD_TE, (len(seq_Inter)))
                
            pos_end_insert = pos_insert+(len(seq_TE)-1)
            #print(" pos_insert " + str(pos_insert) + " pos_end_insert " + str(pos_end_insert))

            ## Case  where (a) TE(s) are already inserted in the considered intergenic sequence ##
            if i in insert:
                inter = 0
                for j in range(len(insert[i])):
                    #print("j " + str(j))#test
                    start_TE_prev = insert[i][j][1]
                    end_TE_prev = insert[i][j][2]

                    ## check if the new TE is inserted before the previous TE ##
                    if pos_insert <= start_TE_prev:
                        #print("test1 "+str(insert[i])) #test
                        #print("test1 "+str(len(seq_TE))) #test
                        inter = j
                        ## Update all TE position upstream
                        while inter != len(insert[i]) :
                            insert[i][inter][1] = insert[i][inter][1]+(len(seq_TE))+ TSD_TE
                            insert[i][inter][2] = insert[i][inter][2]+(len(seq_TE))+ TSD_TE
                            inter+=1

                        df_Inter_inserted.at[i,'sequence'] = str(seq_Inter[0:pos_insert-1]+seq_TE+seq_Inter[(pos_insert-TSD_TE)-1:])
                        df_Inter_inserted_masked.at[i,'sequence'] = str(seq_Inter[0:pos_insert-1]+seq_TE_masked+seq_Inter[(pos_insert-TSD_TE)-1:])
                        df_Inter_inserted.at[i,'length'] = len(df_Inter_inserted.at[i,'sequence'])
                        insert[i].insert(j,[id_TE, pos_insert, pos_end_insert, div_TE, strand_TE, TSD_TE])
                        #print(insert[i])


                    ## check if the new TE is inserted after the previous TE ##
                    elif j == len(insert[i])-1 and pos_insert >= end_TE_prev :
                        df_Inter_inserted.at[i,'sequence'] = str(seq_Inter[0:pos_insert-1]+seq_TE+seq_Inter[(pos_insert-TSD_TE)-1:])
                        df_Inter_inserted_masked.at[i,'sequence'] = str(seq_Inter[0:pos_insert-1]+seq_TE_masked+seq_Inter[(pos_insert-TSD_TE)-1:])
                        df_Inter_inserted.at[i,'length'] = len(df_Inter_inserted.at[i,'sequence'])
                        insert[i].append([id_TE, pos_insert, pos_end_insert, div_TE, strand_TE, TSD_TE])
                        #print(insert[i])


                    ## Nested TEs ##
                    ## Check if the new TE will be inserted inside the previous TE ##
                    elif pos_insert < end_TE_prev and pos_insert > start_TE_prev:
                        #print("test nested "+str(insert[i])) #test
                        inter = j

                        insert[i][j][2] = end_TE_prev+(len(seq_TE))+ TSD_TE
                        if len(insert[i][j]) == 6:
                            insert[i][j].append('nested')
                        #print("test ori "+str(insert[i]))
                        #print(pos_insert)
                        inter = j
                        test1 = inter+1 
                        ## Check for all TEs which are already nested in the first TE and update positions ##
                        while inter != len(insert[i])-1 :
                            #print("test2 "+str(insert[i]))
                            inter+=1
                            start_TE_nest = insert[i][inter][1]
                            end_TE_nest = insert[i][inter][2]

                            ## If the new TE is inserted in the TE i
                            if pos_insert >= start_TE_nest and pos_insert <= end_TE_nest:
                                #print('test1')
                                insert[i][inter][2] = end_TE_nest+(len(seq_TE))+ TSD_TE
                                #print(insert[i][inter][2])
                                test1 = inter+1
                                if len(insert[i][inter]) <= 7:
                                    insert[i][inter].append('nested*')

                            ## If the new TE is inserted before the TE i (in or out of the nested TE) ##
                            elif pos_insert <= start_TE_nest:
                                #print('test2')
                                insert[i][inter][1] = insert[i][inter][1]+(len(seq_TE))+ TSD_TE
                                insert[i][inter][2] = insert[i][inter][2]+(len(seq_TE))+ TSD_TE
                                test1 = inter

                            ## If the new TE is inserted after the TE i ##
                            elif pos_insert > end_TE_nest:
                                test1 = inter

                        insert[i].insert(test1,[id_TE, pos_insert, pos_end_insert, div_TE, strand_TE, TSD_TE, '*'])
                        df_Inter_inserted.at[i,'sequence'] = str(seq_Inter[0:pos_insert-1]+seq_TE+seq_Inter[(pos_insert-TSD_TE)-1:])
                        df_Inter_inserted_masked.at[i,'sequence'] = str(seq_Inter[0:pos_insert-1]+seq_TE_masked+seq_Inter[(pos_insert-TSD_TE)-1:])
                        df_Inter_inserted.at[i,'length'] = len(df_Inter_inserted.at[i,'sequence'])
                        break


            ## Case where there is no TE in the considered intergenic sequence ## 
            else:
                insert[i] = [[id_TE, pos_insert, pos_end_insert, div_TE, strand_TE, TSD_TE]]
                df_Inter_inserted.at[i,'sequence'] = str(seq_Inter[0:pos_insert-1]+seq_TE+seq_Inter[(pos_insert-TSD_TE)-1:])
                df_Inter_inserted_masked.at[i,'sequence'] = str(seq_Inter[0:pos_insert-1]+seq_TE_masked+seq_Inter[(pos_insert-TSD_TE)-1:])
                df_Inter_inserted.at[i,'length'] = len(df_Inter_inserted.at[i,'sequence'])
            index+=1
            iter+=1
            #print(index, id_TE,df_Inter_inserted.iloc[i])
            
    ## Write TEs in each intergenic
    insert_ref = copy.deepcopy(insert)
    list_del = []
    for i in range(int(len(insert_ref)/10)):
        list_del.append(choice(list(insert_ref.keys())))
        insert_ref.pop(list_del[-1])
    od_insert_ref = collections.OrderedDict(sorted(insert_ref.items()))
    #with open('simulated_genome/insert_TEs.csv', 'w') as csv_file: 
    with open('simulated_genome/info_IR_inserted_TEs.csv', 'w') as csv_file: #filename changed
        writer = csv.writer(csv_file)
        for key, value in od_insert_ref.items():
            writer.writerow([key, value])

    df_Inter_inserted_ref = df_Inter_inserted.copy(deep =True)
    for i in list_del:
        df_Inter_inserted_ref.at[i,'sequence'] = df_inter_copy.at[i,'sequence']
        df_Inter_inserted_masked.at[i,'sequence'] = df_inter_copy.at[i,'sequence']
        df_Inter_inserted_ref.at[i,'length'] = len(df_Inter_inserted_ref.at[i,'sequence'])
        df_Inter_inserted_masked.at[i,'length'] = len(df_Inter_inserted_masked.at[i,'sequence'])
        
    ## Write data in file
    print("     Write gene and intergenic information information to file")
    data_insert = concat_seq(df_Inter_inserted_ref, df_gene)
    data_insert[1].to_csv("data/intergenic_inser_tab.tsv",sep ='\t')
    data_insert[2].to_csv("data/gene_inser_tab.tsv",sep ='\t')
    
    print("     Construct information for complete genome")
    data_insert_test = concat_seq(df_Inter_inserted, df_gene)
    clean_test = data_insert_test[0]
    clean_Insert = data_insert[0]
    df_Inter_inserted_ref = data_insert[1]

    print("     Construct information for masked genome")
    data_insert_masked = concat_seq(df_Inter_inserted_masked, df_gene)
    clean_Insert_masked = data_insert_masked[0]
    df_Inter_inserted_masked = data_insert_masked[1]
    
    return (df_Inter_inserted_ref, insert_ref, clean_Insert, clean_Insert_masked, df_Inter_inserted_masked, insert, df_Inter_inserted)
    
## Write TE annotation in file
def annotation_TEs(insert, df_Inter_inserted, clean_Insert, chr, title):
    annot_TEs = {}
    annot_TEs_correct = {}
    index = 0
    for key, value in insert.items():
        for i in value:
            add_len = df_Inter_inserted.at[key,'Pos_Start']
            if len(i) == 6:
                annot_TEs[i[0]] = (i[1]+add_len, i[2]+add_len, i[2]-i[1]+1 , i[3], 'none', i[4], i[5], key)
            if len(i) == 7:
                    if i[6] =='*':
                        annot_TEs[i[0]] = (i[1]+add_len, i[2]+add_len, i[2]-i[1]+1 , i[3],'insert*', i[4], i[5], key)
                    if i[6] =='nested':
                        annot_TEs[i[0]] = (i[1]+add_len, i[2]+add_len, i[2]-i[1]+1 , i[3],'nested', i[4], i[5], key)
            elif len(i) == 8:       
                annot_TEs[i[0]] = (i[1]+add_len, i[2]+add_len, i[2]-i[1]+1 , i[3],'insert&nested', i[4], i[5], key)
            index+=1

    df_annot_TEs = pd.DataFrame.from_dict(annot_TEs, orient='index', columns=['Pos_Start', 'Pos_End', 'length', 'Div', 'type', 'strand', 'TSD', 'index_intergenic'])
    df_annot_TEs = df_annot_TEs.sort_values(by=['Pos_Start'])
    
    ## Write distance between TEs
    Dist1 = []
    Dist2 = []
    GC_dist1 = []
    GC_dist2 = []
    pos_Start = df_annot_TEs.Pos_Start.tolist()
    pos_End = df_annot_TEs.Pos_End.tolist() 
    for i in range(len(pos_Start)):
        if i == 0:
            Dist1.append(0)
            Dist2.append((pos_Start[i+1]-pos_End[i])-1)
            GC_dist1.append(GC(clean_Insert[0:pos_Start[i]-1]))
            if Dist2[-1] <0:
                GC_dist2.append(0)
            else:
                GC_dist2.append(GC(clean_Insert[pos_End[i]+1:pos_Start[i+1]-1]))
                
        elif i == len(pos_Start)-1:
            Dist1.append((pos_Start[i]-pos_End[i-1])-1)
            Dist2.append(0)
            GC_dist2.append(GC(clean_Insert[pos_End[i]:-1]))
            if Dist1[-1] <0:
                GC_dist1.append(0)
            else:
                GC_dist1.append(GC(clean_Insert[pos_End[i-1]+1:pos_Start[i]-1]))
                
        else:
            Dist1.append((pos_Start[i]-pos_End[i-1])-1)
            Dist2.append((pos_Start[i+1]-pos_End[i])-1)
            if Dist2[-1] <0:
                GC_dist2.append(0)
            else:
                GC_dist2.append(GC(clean_Insert[pos_End[i]+1:pos_Start[i+1]-1]))
            if Dist1[-1] <0:
                GC_dist1.append(0)
            else:
                GC_dist1.append(GC(clean_Insert[pos_End[i-1]+1:pos_Start[i]-1]))
                                
    df_annot_TEs['Dist1'] = Dist1
    df_annot_TEs['Dist2'] = Dist2
    df_annot_TEs['GC1'] = GC_dist1
    df_annot_TEs['GC2'] = GC_dist2
    df_annot_TEs.to_csv(title+".tsv",sep ='\t')
    
    return df_annot_TEs

## Create simulated genome with half of the TE insertions compared to the complete simulated genome
def create_Del(df_Inter_inserted_raw, insert_raw, df_inter, df_gene, df_Inter_inserted_masked, chr, turn):
    print("  Generate half inserted genome(s)")
    with click.progressbar (range(int(turn))) as bar:
        for iter in bar:
            list_del = []
            insert_turn= {}
            insert_turn = copy.deepcopy(insert_raw)

            for i in range(int(len(insert_turn)/2)):
                test = choice(list(insert_turn.keys()))
                list_del.append(test)
                insert_turn.pop(list_del[-1])
            df_Inter_inserted_turn = df_Inter_inserted_raw.copy(deep =True)
            df_Inter_inserted_masked_turn =  df_Inter_inserted_masked.copy(deep =True)
            
            for j in list_del:
                # return the orginal intergenic region without TEs
                df_Inter_inserted_turn.at[j,'sequence'] = df_inter.at[j,'sequence']
                df_Inter_inserted_masked_turn.at[i,'sequence'] = df_inter.at[j,'sequence']
                
                df_Inter_inserted_turn.at[j,'length'] = len(df_Inter_inserted_turn.at[j,'sequence'])
                df_Inter_inserted_masked_turn.at[j,'length'] = len(df_Inter_inserted_masked_turn.at[j,'sequence'])
                
            #Deleted genome
            print("     Construct deleted genome")
            data_insert_del = concat_seq(df_Inter_inserted_turn, df_gene)
            clean_Insert_del = data_insert_del[0]
            df_Inter_inserted_del = data_insert_del[1]
            
            # Masked genome
            print("     Construct masked genome")
            data_insert_del_masked = concat_seq(df_Inter_inserted_masked_turn, df_gene)
            clean_Insert_del_masked = data_insert_del_masked[0]
            df_Inter_inserted_del_masked = data_insert_del_masked[1]
            
            # Write annotation
            df_annot_TEs_del = annotation_TEs(insert_turn, df_Inter_inserted_del, clean_Insert_del, chr, "deleted_genome/annot_TEs_del_"+str(iter+1))
            
            # Write sequence fasta
            clean_record = SeqRecord(Seq(clean_Insert_del, generic_dna),id = chr, description ='')
            insert_turn= {}
            with open("deleted_genome/"+args.out+"_Del_"+str(iter+1)+".fasta", "w") as output_handle:
                SeqIO.write(clean_record, output_handle, "fasta")
                
            # Write sequence fasta
            clean_record_masked = SeqRecord(Seq(clean_Insert_del_masked, generic_dna),id = chr, description ='')
            insert_turn= {}
            with open("deleted_genome/"+args.out+"_Del_"+str(iter+1)+"_masked.fasta", "w") as output_handle:
                SeqIO.write(clean_record_masked, output_handle, "fasta")
    print("  Done."+'\n')


def main():
    print('***************STEP3 inseraTE: generate simulated genomic sequences***************'+'\n')

    ## Load Data
    df_Inter = pd.read_csv(args.int, delimiter ='\t')
    df_gene = pd.read_csv(args.gen, delimiter ='\t')
    df_TEs = pd.read_csv(args.par, delimiter ='\t')
    data_insert = insert_TEs(args.TEs, df_Inter, df_gene, df_TEs)
    df_Inter_inserted = data_insert[0]
    insert = data_insert[1]
    clean_Insert = data_insert[2]
    clean_Insert_masked = data_insert[3]
    df_Inter_inserted_masked = data_insert[4]
    insert_ori = data_insert[5]
    df_Inter_inserted_ori = data_insert[6]
    ## find description 
    record = SeqIO.read(args.des, "fasta")
    desc_tot = record.description
    chr = desc_tot[desc_tot.find("chromosome")+11:]

    
    ## Write annotation
    df_annot_TEs = annotation_TEs(insert, df_Inter_inserted, clean_Insert, chr, "simulated_genome/annot_TEs")
    print("  Number of inserted TEs: "+str(len(df_annot_TEs)))
    print("  Done."+'\n')
    
    ## Create deleted genome
    if args.met: #if --met argument used, the deleted genome will also contain insertions not present in the simulated genome
        create_Del(df_Inter_inserted_ori, insert_ori, df_Inter, df_gene,df_Inter_inserted_masked, chr, args.ind)
    else:
        create_Del(df_Inter_inserted, insert, df_Inter, df_gene,df_Inter_inserted_masked, chr, args.ind)

    #Write masked sequence in fasta format
    clean_record1 = SeqRecord(Seq(clean_Insert_masked, generic_dna),id = chr, description = '')
    with open("simulated_genome/"+args.out+"_masked.fasta", "w") as output_handle:
        SeqIO.write(clean_record1, output_handle, "fasta")
    
    ## Write final sequence in fasta format
    clean_record2 = SeqRecord(Seq(clean_Insert, generic_dna),id = chr, description ='')
    with open("simulated_genome/"+args.out+".fasta", "w") as output_handle:
        SeqIO.write(clean_record2, output_handle, "fasta")
    
    ## Testing our density
    print("  Final estimated gene density: "+str((len(df_Inter)-1)/(len(clean_Insert)/10**6)))
    print("  Fasta file created.")
    
if __name__ == "__main__":
    main()
