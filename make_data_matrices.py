#!/usr/bin/env python

import sys, os
import argparse
import pandas as pd
import numpy
from Bio import Phylo

parser = argparse.ArgumentParser(description='Script to create data matrices for input to PGLS')
parser.add_argument('-s', dest = 'SpFile', type = str, required=True,  help = 'Path to comma-delimited file with the names of all the species in the analysis separated by the state of interest (0 or 1/ C or S) and path to tab-separated expression matrix (of TPM values) for that species')
parser.add_argument('-f', dest= 'Dir', type = str, required=True, help ='Path to directory of gene family trees')
parser.add_argument('-e', dest= 'Ending', type = str, required=True, help ='Tree file ending')
args = parser.parse_args()

DIR = args.Dir
if DIR[-1] != "/": DIR += "/"

state_dict = {}
exp_dict = {}

with open(args.SpFile, 'rU') as f:
	for line in f:
		state_dict[line.split(',')[0]] = line.split(',')[1]
		exp_dict[line.split(',')[0]] = line.split(',')[2].strip('\n')
		
for file in os.listdir(DIR):
    if file.endswith(args.Ending):
        clust = file.split('.')[0]
        print clust
        output = open(DIR+clust+'.dat.csv', 'w')
        output.write(','.join(['Gene','State','ExpMean','ExpVar\n']))
        tree = Phylo.read(DIR+file, "newick")
        tips = []
        for leaf in tree.get_terminals():
            tips.append(leaf.name)
        for tip in tips:
            species = tip.split('@')[0]
            trans = tip.split('@')[1]
            dat = pd.read_csv(exp_dict[species], sep='\t', index_col=0)
            output.write(','.join([tip,state_dict[species],str(numpy.mean(dat.loc[trans])),str(numpy.var(dat.loc[trans])+0.000000000001)+'\n']))


'''
for line in args.SpFile
	make dict of species and state
	make another dict of species and path to expression matrix
	#just open the expression matrix file when needed and close each time using the dict

#select row with pandas
# df.loc['row_name']
# numpy.mean(df.loc['row_name'])
'''