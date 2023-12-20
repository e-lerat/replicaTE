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

## Extract data (postions, length, sequences) from "record" genbank file
def extract_data(record, Pos_chr):
    seq_Genome = str(record.seq)
    genes_Data = {}
    TEs_Data = {}
    n=0
    for feature in record.features:
        
        ## Extract gene data 
        if feature.type == "gene":
            
            if 'gene' in feature.qualifiers:
                gene_Name = feature.qualifiers['gene'][0]
            elif 'gene_id' in feature.qualifiers:
                gene_Name = feature.qualifiers['gene_id'][0]
            elif 'locus_tag' in feature.qualifiers:
                gene_Name = feature.qualifiers['locus_tag'][0]
                
            #print (gene_Name) #test
            location = str(feature.location)
            # Case of multiple locations
            if location.find("order") !=-1:
                print ("ignored gene due to multiple locations: ", gene_Name)
                #location = sub('[^0-9]', ' ', location[7:-5])
               # location = location.split('       ')
               # for i in range(len(location)):
                   # if location[i][0] ==' ':
                   #     location[i] = location[i][1:]
                   # pos_Start = int(location[i][:location[i].find(' ')])
                  #  pos_end = int(location[i][location[i].find(' ')+1:])
                  #  gene_Seq = seq_Genome[int(pos_Start):int(pos_End)]
                  #  pos_Start = pos_Start + Pos_chr
                 #   pos_end = pos_end + Pos_chr
                   # genes_Data[gene_Name+"_"+str(i)] = (pos_Start, pos_End, gene_length, gene_Seq)
            else:
                location = sub('[^0-9]', ' ', location)
                pos_Start = int(location[1:location[2:].find(' ')+2].replace(' ',''))
                pos_End = int(location[location[2:].find(' ')+2:].replace(' ',''))
                gene_Seq = seq_Genome[int(pos_Start):int(pos_End)]
                gene_length = len(gene_Seq)
                pos_Start = pos_Start + Pos_chr
                pos_End = pos_End-1 + Pos_chr
                genes_Data[gene_Name] = (pos_Start, pos_End, gene_length, gene_Seq)
            
            ## Extract TE data
            
            ##The name of the TE copies have to be different inside a given family in the genbank file (for example format like familly{}number)
        if feature.type == "mobile_element":
            ## retrieve a usable name
            if 'mobile_element_type' in feature.qualifiers:
                name_TE = feature.qualifiers['mobile_element_type'][0]
                name_TE = name_TE[name_TE.find(':')+1:]
            elif 'note' in feature.qualifiers:    
                name_TE = feature.qualifiers['note'][0]
                #name_TE = name_TE[name_TE.find('_'):-1]
            
            location = str(feature.location)
            strand = (sub('[^-+]', '', location))
            location = (sub('[^0-9]', ' ', location))
            pos_Start = int(location[1:location[2:].find(' ')+2].replace(' ',''))
            pos_End = int(location[location[2:].find(' ')+2:].replace(' ',''))
            #print(name_TE+"\t"+str(location)+"\t"+strand+"\t"+str(pos_Start)) #test
            if strand == '-':
                TE_Seq = str(Seq(seq_Genome[int(pos_Start):int(pos_End)]).reverse_complement())
            else:
                TE_Seq = seq_Genome[int(pos_Start):int(pos_End)]
            TE_length = len(TE_Seq)
            pos_Start = pos_Start + Pos_chr
            pos_End = pos_End + Pos_chr
            TEs_Data[name_TE] = (pos_Start, pos_End, TE_length, TE_Seq)
            n=n+1
            #print(TEs_Data[name_TE])#test
            
    #print(TEs_Data)#test
    #print("nb of extracted mobile element"+"\t"+str(n)) #test
    #print("taille TEs_Data = nombre de copies"+"\t"+str(len(TEs_Data))) #test
    df_Gene = pd.DataFrame.from_dict(genes_Data, orient='index', columns=['pos_Start', 'pos_End', 'length', 'sequence'])
    df_TE = pd.DataFrame.from_dict(TEs_Data, orient='index', columns=['pos_Start', 'pos_End', 'length', 'sequence'])
    #print("longueur df_TE fin extract_data" + "\t" + str(len(df_TE)))#test
    #print("longueur df_Gene fin extract_data" + "\t" + str(len(df_Gene)))#test

    return (df_Gene, df_TE)

## Read genbank file and extract data
def parser_Genbank(genbank_file):
    seq_Desc = list()
    Pos_chr = 0
    seq_genome = ''
    try:
        ## Mutiple record
        seq_Id = list()
        #print("dans multiple")#test
        df_Gene = pd.DataFrame(columns = ['pos_Start', 'pos_End', 'length', 'sequence'])
        df_TE = pd.DataFrame(columns = ['pos_Start', 'pos_End', 'length', 'sequence'])
        for record in SeqIO.parse(genbank_file, "genbank"):
            data_gen = extract_data(record, Pos_chr)
            df_iter_Gene = data_gen[0]
            df_iter_TE = data_gen[1]
            Pos_chr = Pos_chr + len(record.seq)
            df_Gene = df_Gene.append(df_iter_Gene)
            df_TE = df_TE.append(df_iter_TE)
            seq_genome = seq_genome + record.seq
            seq_Id.append(record.id)
            organism = record.annotations["source"]
            seq_Desc.append(record.description)

    except:
        ## One record
        print("dans one record")#test
        record = SeqIO.read(genbank_file, "genbank")
        data_gen = extract_data(record, Pos_chr)
        df_Gene = data_gen[0]
        df_TE = data_gen[1]
        seq_genome = record.seq
        seq_Id =record.id
        seq_Desc.append(record.description)
        organism = record.annotations["source"]
    
    #print(df_TE)#test
    df_Gene = df_Gene.rename_axis('id').reset_index()
    df_TE = df_TE.rename_axis('id').reset_index()
    #print(df_TE)#test
    #print("dans parser_genbank"+"\t"+str(len(df_TE))) #test
    print(seq_Id)
    return (df_Gene, df_TE, seq_genome, seq_Id, seq_Desc, organism)

## generate simulated genes based on the chromosome
def simulated_gene(df_Gene):
    np.seterr(divide='ignore', invalid='ignore')
    gc_content = [GC(x) for x in df_Gene.sequence]
    data_len = [len(x) for x in df_Gene.sequence if len(x) > 0 and len(x) !='NaN' and len(x) != 'Inf']
    data_len = sorted(data_len)
    size = len(df_Gene)

    ## Plot
    plt.hist(gc_content, bins= 50)
    plt.title("GC_content_distribution")
    plt.xlim(0, 100)
    plt.savefig("attachment/observed_gene_len_distribution.png", bbox_inches='tight')
    #plt.show()
    plt.close()
    
    plt.hist(data_len, bins= 100)
    plt.title("Length_distribution")
    plt.xlim(0, max(data_len))
    plt.savefig("attachment/observed_gene_gc_content_distribution.png", bbox_inches='tight')
    #plt.show()
    plt.close()
      
    mu, std = stats.norm.fit(gc_content)
    distrib_GC = distrib_truncnorm(size, mu, std, 15, max(gc_content), "gene_GC_content")
    
    a, c, fit_loc, fit_scale = stats.exponweib.fit(data_len, loc =min(data_len), scale=np.mean(data_len))
    distrib_len = distrib_exponweib(data_len, a, c, fit_loc, fit_scale, 'y', "gene_length")
#    distrib_len = [2500] * size        ## testing, make length of gene = 2500 bp
    print("  Simulate genes ")
    df_new_gene = simulate_Seq(distrib_len, distrib_GC)
    
    return df_new_gene



## Check if TEs are inside gene sequence and replace them with 'N'
def search_TEs_gene(df_TE, df_Gene):
    iter = 0
    count = 0
    nb_TE = len(df_TE)
    print("Number of extracted TEs: " + str(nb_TE)) #test
    ## loading bar
    print("  Check TEs in gene sequence and delete them")
    with click.progressbar (range(len(list(df_TE.iterrows())))) as bar:
        for row in df_TE.iterrows():
            
            bar.update(1)
            btm = row[1].pos_Start
            top = row[1].pos_End

            ## localize TEs with position with Start and End
            for row in df_Gene.iloc[iter:].iterrows():
                
                gene_Start = int(row[1].pos_Start)
                gene_End = int(row[1].pos_End)
                gene_Seq = str(row[1].sequence)

                if (btm > gene_Start) and (btm < gene_End):
                    if (top < gene_End):
                        sub_TE = 'N'*(top-btm)
                        df_Gene.at[iter,'sequence'] = gene_Seq.replace(gene_Seq[btm-gene_Start:top-gene_Start],sub_TE,1) 
                        
                        count+=1
                    else:
                        sub_TE = 'N'*(gene_End-btm)
                        df_Gene.at[iter,'sequence'] = gene_Seq.replace(gene_Seq[btm-gene_Start:],sub_TE,1)
                        count+=1
                    break

                elif (btm < gene_Start) and (top > gene_Start) and (top < gene_End):
                    sub_TE = 'N'*(top-gene_Start)
                    df_Gene.at[iter,'sequence'] = gene_Seq.replace(gene_Seq[0:top - gene_Start],sub_TE,1)
                    count+=1
                    break

                elif (top < gene_Start):
                    break
                
                else:
                    iter+=1
                    
    print("  Number of TEs annotated in genes / total number of TEs: "+str(count)+"/"+str(nb_TE))
    print("  Done."+'\n')
    return df_Gene

#Delete 'nested' genes
def search_nested(df_Gene):
    df_Gene.sort_values(by=['pos_Start'])
    gene_Start = df_Gene.pos_Start.tolist()
    gene_End = df_Gene.pos_End.tolist()
    gene_nested = []
    gene_name=df_Gene.id.tolist()
    i = 0
    while i < len(gene_End)-1:
        j=i+1
        inter = 0
        #print(str(len(gene_End))+"\t"+str(i)+"\t"+str(j))#test
       
        while gene_End[i] > gene_Start[j] and j < len(gene_End)-2:
            #print(gene_name[i]+"\t"+gene_name[j]+"\t"+str(i)+"\t"+str(j)+"\t"+str(gene_End[i])+"\t"+str(gene_Start[j]))#test
            gene_nested.append(j)
            inter+=1
            j=j+1
            
        i=i+inter+1
    df_Gene = df_Gene.drop(df_Gene.index[gene_nested]).reset_index(drop=True)
    return df_Gene

## Extract intergenic data (GC content, length, sequence)
def search_inter(df_Gene, seq_Genome):
    #print("len(df_Gene) " + str(len(df_Gene))) #test
    seq = str(seq_Genome)
    gene_Start = df_Gene.pos_Start.tolist()
    gene_End = df_Gene.pos_End.tolist()
    all_inter_Seq = []
    inter_gc = []
    ## Loading bar
    print("  Retrieve the intergenic sequences")
    with click.progressbar (range(len(gene_Start))) as bar:
        for i in bar: 
            
            if i == 0:
                inter_seq = seq[0:gene_Start[i]]
                all_inter_Seq.append(inter_seq)
                inter_gc.append(GC(inter_seq))
            else:
                inter_seq = seq[gene_End[i-1]:gene_Start[i]]
                all_inter_Seq.append(inter_seq)
                inter_gc.append(GC(inter_seq))
                 
    inter_seq = seq[gene_End[-1]:]
    all_inter_Seq.append(inter_seq)
    inter_gc.append(GC(inter_seq))
    inter_seq_len = [len(x) for x in all_inter_Seq]
    
    df_inter = pd.DataFrame({'GC_content':inter_gc, 'length': inter_seq_len, 'sequence': all_inter_Seq})
    #print("df_inter " + str(len(df_inter))) #test
    
    print("  Done."+'\n')
    return df_inter

## Delete masked TEs in genes
def delete_TEs(df_Gene):
    new_Gen_seq = []
    new_Gen_len = []
    gene_Seq = df_Gene.sequence.tolist()
    for i in range(len(gene_Seq)):
 
        new_Seq = gene_Seq[i].replace('N','')
        new_Gen_seq.append(new_Seq)
        new_Gen_len.append(len(new_Seq))
        
    ## Create a dataframe with the data
    df_new_Gen = pd.DataFrame({'id': df_Gene.id.tolist(), 'length': new_Gen_len, 'sequence': new_Gen_seq})
    return df_new_Gen

## Retrieve features from TEs (copy_nb, mean_length), do not consider families with a max length < 200 bp
def caract_Fam(df_TE):
    sorted_df_TEs = df_TE.sort_values(by=['id', 'length'], ascending=False)
    fam_id = [id[0:id.find('{')].replace(' ','') for id in sorted_df_TEs.id.tolist()]
    sorted_df_TEs['id'] = fam_id
    sorted_df_TEs = sorted_df_TEs.sort_values(by=['id', 'length'], ascending=False)
    nb_copy = 0
    Fam_Stat = {}
    length  = []
    seq  = []
    name = ''
    name_del = ''
    ligne=0#test
    #print(sorted_df_TEs)#test
    for row in sorted_df_TEs.iterrows():
        ligne+=1 #test
        if row[1].id != name:
            if nb_copy != 0:
#                Fam_Stat[name] = (nb_copy, sum(length)/len(length), generate_Seq(length[0], GC(seq)))
                Fam_Stat[name] = (nb_copy, sum(length)/len(length), seq)
                nb_copy = 0
                length = []
                name = ''
            else:
                name = ''

        if row[1].length > 200 and name =='':
            name = row[1].id
            length.append(row[1].length)
            seq = row[1].sequence

        if row[1].id == name:
            length.append(row[1].length)
            nb_copy+=1

        if row[1].length < 200 and name =='':
            name_del = row[1].id
            sorted_df_TEs = sorted_df_TEs[sorted_df_TEs.id != name_del]
            
        if ligne == len(sorted_df_TEs):#test to make sure it takes into account the last family
            Fam_Stat[name] = (nb_copy, sum(length)/len(length), seq)
            nb_copy = 0
            length = []
            name = ''
            
    stat_Chr_TEs = pd.DataFrame.from_dict(Fam_Stat, orient='index', columns=['copy_nb', 'mean_length', 'ancestral_seq'])
    #print(stat_Chr_TEs)#test
    return stat_Chr_TEs
    
## Give mean nb of intergenic sequences
def gene_Density(nb_Gene, len_Seq_tot):
    mean_Inter_Gene = nb_Gene/(len_Seq_tot/10**6)
    return  mean_Inter_Gene

## Generate simulated intergenic sequences
def simulate_inter(df_inter, gene_Dens, gc_mean, minlen, len_Gene_sum, len_Seq_tot):
    np.seterr(divide='ignore', invalid='ignore')
    data_len = df_inter.length.tolist()
    data_gc = df_inter.GC_content.tolist()
    size = len(data_len)
    ## Look if gc content is correct
    while gc_mean < 15 and gc_mean !=-1:
        gc_mean = float(input("Error: GC_content too low, enter a new GC_content: "))
        
    while gc_mean > 100: 
        gc_mean = float(input("Error: GC_content is too high (>100), enter a new GC_content: "))
        
    if gc_mean == -1:
        mu, std = stats.norm.fit(data_gc)  
    #simulation of the GC content distribution of intergenic regions    
    distrib_GC = distrib_truncnorm(size, mu, std, 15, max(data_gc), "intergenic_GC_content")
    
    ## Look if min length is correct
    while minlen < 1:
        minlen = int(input("Error: min length can't be negatif, enter a new min length: "))
        
    ## Look if density is correct
    mean_init = gene_Density(size-1, len_Seq_tot)
    min_Len = ((mean_init)*(size-1))*0.50
    length_Test = (((size-1)/gene_Dens)*10**6)-len_Gene_sum
    while length_Test <= min_Len and gene_Dens !=-1:
        gene_Dens = int(input("Error: this density creates negatif length, enter a new density: "))
        length_Test = (((size-1)/gene_Dens)*10**6)-len_Gene_sum
        
    while gene_Dens < 0 and gene_Dens !=-1:
        gene_Dens = int(input("Error: density can't be negatif, enter a new density: ")) 
        
    if gene_Dens == -1:
        gene_Dens = gene_Density(size-1, len_Seq_tot)
        print("  Initial gene density: "+str(gene_Dens))
        mean =((((size-1)/gene_Dens)*10**6)-len_Gene_sum)/(size-1)
        a, c, fit_loc, fit_scale = stats.exponweib.fit(data_len, scale=mean)
        distrib_len = distrib_exponweib(data_len, a, c, minlen, fit_scale, 'y', "intergenic_length")
#        distrib_len = [2] * size       ## testing, make length of IR = 2 bp

    else: 
        mean =((((size-1)/gene_Dens)*10**6)-len_Gene_sum)/(size-1)
        distrib_len = distrib_exponweib(data_len, 1, 0.89, minlen, mean, 'y', "intergenic_length")
        
    ## Plot observed intergenic GC% and length distribution found in the "attachment" directory
    plt.hist(data_gc, bins = 150)
    plt.xlim(0, 100)
    plt.title("Genomic intergenic GC content distribution")
    plt.savefig("attachment/observed_intergenic_GC_content_distribution.png", bbox_inches='tight')
    #plt.show()
    plt.close()
    plt.hist(data_len, bins = 150)
    #plt.xlim(0, max(distrib_len))
    plt.xlim(0, max(data_len))
    plt.title("Genomic intergenic length distribution")
    plt.savefig("attachment/observed_intergenic_length_distribution.png", bbox_inches='tight')
    #plt.show()df_Gene
    plt.close()
    
    #print("dans simulate_inter df_inter " + str(len(df_inter))) #test
    df_inter = simulate_Seq(distrib_len, distrib_GC) 
    
    return df_inter


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

        #clean_Seq = clean_Seq+new_Gen_list[i]

    clean_Seq = clean_Seq + inter_Seq_list[-1]
    Pos_Start_intergenic.append(Pos_End_Gene[-1]+1)
    Pos_End_intergenic.append(Pos_Start_intergenic[-1]+len(inter_Seq_list[-1])-1)
    df_inter_Seq["Pos_Start"] = Pos_Start_intergenic
    df_inter_Seq["Pos_End"] = Pos_End_intergenic
    ## Change order of the columns
    cols = df_inter_Seq.columns.tolist()
    cols = cols[-2:] + cols[:-2]
    df_inter_Seq = df_inter_Seq[cols]
    df_new_Gen["Pos_Start"] = Pos_Start_Gene
    df_new_Gen["Pos_End"] = Pos_End_Gene 
    cols = df_new_Gen.columns.tolist()
    cols = cols[-2:] + cols[:-2]
    df_new_Gen = df_new_Gen[cols]
    return(clean_Seq, df_inter_Seq, df_new_Gen)

def main():
    
    print('***************STEP1 deleTE: retrieve information from genes, intergenic regions and TEs***************'+'\n')
    #create directories containing result files
    os.makedirs('data', exist_ok=True)
    os.makedirs('attachment', exist_ok=True)
    os.makedirs('deleted_genome', exist_ok=True)
    os.makedirs('simulated_genome', exist_ok=True)

    ## Load data
    data = parser_Genbank(args.gb)
    df_Gene = data[0]
    df_Gene = search_nested(df_Gene)        # Delete 'insert' genes
    df_TE = data[1]
    seq_gen = data[2]
    len_Seq_tot = len(seq_gen)
    seq_Id = data[3]
    seq_Desc = data[4]
    organism = data[5]
    nb_Gene = len(df_Gene)
    
    df_TE.to_csv("attachment/ori_TEs_tab.tsv",sep ='\t')
    
    ## Check if we used simulated genes
    if args.sim == True:
        df_new_Gen = simulated_gene(df_Gene)
        sum_Gene_tot = df_new_Gen.length.sum()

    else:
        df_Gene = search_TEs_gene(df_TE, df_Gene)       # Dataframe with masked genes
        df_new_Gen = delete_TEs(df_Gene)        # Dataframe of 'clean' genes
        sum_Gene_tot = df_new_Gen.length.sum()      # Sum length of all sequences
        df_Gene.to_csv("attachment/gene_markedTEs_tab.tsv",sep ='\t')

    #print("df_Gene " + str(len(df_Gene))) #test
    #print("df_new_Gen " + str(len(df_new_Gen))) #test
    
    ## Generate the intergenic regions with fitted data
    print("  Simulate intergenic regions ")
    df_data_inter = search_inter(df_Gene, seq_gen)
    df_simulated_inter = simulate_inter(df_data_inter, args.dens, args.gc, args.min, sum_Gene_tot, len_Seq_tot)
    
    ## Delete too long TEs (> 11k bp)
    df_TE = df_TE.loc[df_TE['length'] < 11000]
    df_TE = df_TE.reset_index(drop=True)
    #print(df_TE)#test

    ## Check the statistic of TEs 
    df_stat_fam_TEs = caract_Fam(df_TE)
    #print(df_stat_fam_TEs)#test
    ## Check the statistic of all TE copy (if provided)
    if args.TEs != 'NA':
        length_TEs = []
        name_TEs = []
        Seq_TEs = []
        for record in SeqIO.parse(args.TEs, "fasta"):
            ident_correct = record.description[record.description.find('name=')+5:]
            name_TEs.append(ident_correct)
            length_TEs.append(len(record.seq))
            Seq_TEs.append(str(record.seq))
            
        df_TEs_genome = pd.DataFrame({'id':name_TEs, 'length': length_TEs, 'sequence': Seq_TEs})
        df_TEs_genome = df_TEs_genome.loc[df_TEs_genome['length'] < 11000]      # Delete TEs > 11 kb in case TE fasta file is provided
        df_stat_fam_TEs = caract_Fam(df_TEs_genome)
        df_stat_fam_TEs.to_csv("data/stat_TEs_tab.tsv",sep ='\t')
    else:
        df_stat_fam_TEs.to_csv("data/stat_TEs_tab.tsv",sep ='\t')

    ## Create clean sequences
    data_clean = concat_seq(df_simulated_inter, df_new_Gen)
    #print("data_clean " + str(data_clean[1])) #test
    
    ## Write dataframe containing the cleaned genes (without N) and simulated intergenic regions into files
    data_clean[1].to_csv("data/intergenic_sim_tab.tsv",sep ='\t')
    data_clean[2].to_csv("data/gene_clean_tab.tsv",sep ='\t')
    clean_Seq = data_clean[0]

    ## Write clean_Seq (simulated genome with clean genes and simulated intergenic regions) in fasta format
    if isinstance(seq_Id, list):
        seq_Id = ''.join(seq_Id)
        #print(seq_Id)
    if len(seq_Desc) > 1:
        seq_Desc = organism
    else:
        seq_Desc= ''.join(seq_Desc)
    clean_record = SeqRecord(Seq(clean_Seq, generic_dna),id = seq_Id, description = seq_Desc+" ")
    with open("attachment/"+args.out+".fasta", "w") as output_handle:
        SeqIO.write(clean_record, output_handle, "fasta")
        
    ## Testing our density
    print("  Estimated gene density: "+str(nb_Gene/(len(clean_Seq)/10**6)))
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulate a genome without TEs.")
    parser.add_argument("-g", "--gb", help ="Genbank file.")
    parser.add_argument("-s", "--sim", help ="Uses the characteristics of genes to create simulated genes.", action="store_true")
    parser.add_argument("-m", "--min", type=int, help ="Minimal length of intergenic sequences (can break the code).", default = 100)
    parser.add_argument("-d", "--dens", type=int, help ="Specify Gene density.", default = -1)
    parser.add_argument("-c", "--gc", type=int, help ="Specify GC content of intergenic regions.", default = -1)
    parser.add_argument("-t", "--TEs", help ="Fasta file of all TE copies from the genome.", default = 'NA')
    parser.add_argument("-o", "--out",help ="Output name.", default="output_clean_seq")
    args = parser.parse_args()
    main()
    
