#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 12:36:23 2021

@author: earl

ANI_genome_picking takes a matrix of average nucleotide identities and finds the X most divergent sequences
usage:
ANI_genome_picking.py PATH/TO/MATRIX number_of_output_genomes


Matrix should be in the form:
Genome  Seq1   Seq2   Seq3
Seq1    100    98.0   90
Seq2    98.0   100    93
Seq3    90     93     100

This file is created by orthophylo.1.6.slurm in the current git repo 
"""
from collections import Counter
import copy
#import sys
#print (int(sys.argv[2]))
#ANI_matrix=sys.argv[1]
ANI_matrix="../TEST/full_matrix"
#n_clusters=int(sys.argv[2])
n_clusters=10
full_dict = {}
f = open(ANI_matrix)
for line in f:    
    #print (line)
    l = line.replace('"',"")
    l = l.replace("'","")        
    l1 = l.replace("", "")
    l2 = l1.replace("\t\n", "")
    l2 = l2.replace("\n", "")
    l3 = l2.replace("\t\t", "\t")
    l4 = l3.split("\t") 
    #store the names of selected genomes in a dict with their avg ANI to other genomes
    if len(l4) > 1:
        if l4[0] == "Genomes":
            genomes = l4[1:]
        else:
            J=1
            full_dict[l4[0]] = {} #make a dict with genome as key (column 1)
            for I in l4[1:]:
                # for each data column, associate genome name with ANI calc
                full_dict[l4[0]][genomes[J-1]] = float(l4[J])
                J+=1
f.close()


#rm self entries
for key in full_dict:
    del full_dict[key][key]

full_dict_copy = copy.deepcopy(full_dict)

#dont remember what this is for...
def average_dict_values(in_dict):
    res = 0
    for val in in_dict.values():
        res += val
  
    # using len() to get total keys for mean computation
    mean = res / len(in_dict)
    return mean

def cluster(in_dict,first_iteration):
    ANI_max={}
    for key1,val in full_dict.items():
        #print (key1,val)
        max_key = max(val, key=val.get)        
        ANI_max[max_key,key1]=full_dict[key1][max_key]
    max_ANI_pair = max(ANI_max, key=ANI_max.get)
    sums = dict(Counter(full_dict[max_ANI_pair[0]]) + Counter(full_dict[max_ANI_pair[1]]))
    full_dict[max_ANI_pair] = {k: sums[k] / float((k in full_dict[max_ANI_pair[0]]) + (k in full_dict[max_ANI_pair[1]])) for k in sums}
    
    #add other values for other side of matrix for merged genomes
    for I in full_dict[max_ANI_pair]:
        full_dict[I][max_ANI_pair] = full_dict[max_ANI_pair][I]
    
    #print (full_dict.keys())
    for I in max_ANI_pair:
        del full_dict[I]
        #del full_dict[max_ANI_pair][I]
    for I in max_ANI_pair:
        for key in full_dict:
            if key != I:
                del full_dict[key][I]
    return full_dict, ANI_max, max_ANI_pair


first_iteration=0
while len(full_dict) > n_clusters:
    cluster(full_dict,first_iteration)
    first_iteration+=1

NN_distances={}
for I in full_dict:
    NN_distances[I]=str(I).replace(")","").replace("(","").replace(",,",",").replace("'","").replace(" ","").split(",")

dict_dist={}
for I in NN_distances:
    dict_dist[I]={}
    for J in NN_distances[I]:
        dict_dist[I][J]=average_dict_values(full_dict_copy[J])

species_list=[]
for I in dict_dist:
     species_list.append(max(dict_dist[I], key=dict_dist[I].get))

outfile=open("Species_shortlist", 'w')
for k in species_list:
    outline=str(k)
    outfile.write(outline + "\n")
outfile.close()
