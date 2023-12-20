#! /usr/bin/python3
# -*- coding: Utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sps
import scipy.stats as stats
import pandas as pd
import click
from Bio.SeqUtils import GC
from random import choice, uniform, gammavariate, random, randint

# Copyright (C) 2021 Emmanuelle Lerat
 
# This file is part of replicaTE suite, developed by Van Anthony Le.

# replicaTE is a free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# replicaTE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with replicaTE.  If not, see <http://www.gnu.org/licenses/>.

## Truncated exponential distribution 
def distrib_truncexp(N, lower, upper, scale):
    model_Exp = stats.truncexpon(b=(upper-lower)/scale, loc=lower, scale=scale)
    distrib_Length = model_Exp.rvs(N)
    ## Plot the length distribution
    #fig, ax = plt.subplots()
    #ax.hist(distrib_Length, density=True)
    #plt.show()
    
    return distrib_Length.astype(int)

## Truncated normal distribution
def distrib_truncnorm(nb, mu, sigma, lower, upper, title):
    distrib_Div = stats.truncnorm.rvs((lower-mu)/sigma, (upper-mu)/sigma, loc=mu, scale=sigma, size=nb)
    ## Plot the generated distribution
    plt.hist(distrib_Div, bins= 75, density=True, alpha=0.6, color='r')
    plt.xlim(0, max(distrib_Div)+5)
    plt.title(title+" using the truncated normal distribution")
    plt.savefig("attachment/simulated_"+title+"_truncated_normal_distribution.png", bbox_inches='tight')
    #plt.show() 
    plt.close()
    
    return np.around(distrib_Div, decimals = 2)

## Weibull distribution 
def distrib_exponweib(data_len, a, c, loc, scale, draw_plot, title):
    size = len(data_len)
    upper = max(data_len)
    distrib_Length = stats.exponweib.rvs(a, c, loc, scale, size)
    upper = max(distrib_Length)
    if draw_plot == 'y':
        ## Plot the length distribution
        plt.hist(distrib_Length, bins= 75, density=True, alpha=0.6, color='r')
        plt.xlim(0, upper)
        plt.title(title+" using the Weibull distribution")
        plt.savefig("attachment/simulated_"+title+"_Weibull_distribution.png", bbox_inches='tight')
        #plt.show() 
        plt.close()

    return distrib_Length.astype(int)

## Pick a random nucleotide depending of GC count
def pick_nt(GC): 
    RGC = random()
    GC = GC/100
    if GC < RGC :
        if randint(0,1) == 0 :
            return ("A")
        else :
            return ("T")
    else :
        if randint(0,1) == 0 :
            return ("G")
        else :
            return ("C")

## Generate DNA sequence with length and weights
def generate_Seq(length, GC_content):
    seq = []
    for i in range(length):
        nt = pick_nt(GC_content)
        seq.append(nt)
    seq = "".join(seq)
    
    return seq

## Create intergenic sequences with exponential law
def simulate_Seq(listLen, distrib_GC):
    print("  Create intergenic sequences using the Weibull distribution")
    seq = []
    
    with click.progressbar (range(len(listLen))) as bar:
        for i in bar:
            seq.append(generate_Seq(listLen[i], distrib_GC[i]))
            
    df_seq = pd.DataFrame({'length': listLen, 'GC_content':distrib_GC.tolist(),'sequence': seq})
    print("  Done."+'\n')
    
    return df_seq
