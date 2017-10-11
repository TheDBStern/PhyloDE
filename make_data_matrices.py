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
parser.add_argument('--log2', dest= 'lg2', action ='store_true', default= False, help ='log2 transform expression values, default = False.')
args = parser.parse_args()

DIR = args.Dir
if DIR[-1] != "/": DIR += "/"

state_dict = {}
exp_dict = {}

with open(args.SpFile, 'rU') as f:
	for line in f:
		state_dict[line.split(',')[0]] = line.split(',')[1]
		exp_dict[line.split(',')[0]] = pd.read_csv(line.split(',')[2].strip('\n'), sep='\t', index_col=0)
			
for file in os.listdir(DIR):
	if file.endswith(args.Ending):
		clust = file.split('.')[0]
		print clust
		if args.lg2:
			output = open(DIR+clust+'.log2.dat.csv', 'w')
		else:
			output = open(DIR+clust+'.dat.csv', 'w')
		output.write(','.join(['Gene','State','ExpMean','ExpVar','ExpSE\n']))
		tree = Phylo.read(DIR+file, "newick")
		for leaf in tree.get_terminals():
			tip = leaf.name
			species = tip.split('@')[0]
			trans = tip.split('@')[1]
			dat = exp_dict[species]
			if args.lg2:
				output.write(','.join([tip,state_dict[species],str(numpy.mean(numpy.log2(dat.loc[trans]+1))),str(numpy.var(numpy.log2(dat.loc[trans]+1))+0.000000000001),str(numpy.std(numpy.log2(dat.loc[trans]+1))/numpy.sqrt(len(dat.loc[trans]))+0.000000000001)+'\n']))
			else:
				output.write(','.join([tip,state_dict[species],str(numpy.mean(dat.loc[trans])),str(numpy.var(dat.loc[trans])+0.000000000001),str(numpy.std(dat.loc[trans]/numpy.sqrt(len(dat.loc[trans])))+0.000000000001)+'\n']))
