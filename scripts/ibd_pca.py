# -*- coding: utf-8 -*-

import pandas as pd
import csv
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh

samples = pd.read_csv("IBD_GRM_t10.grm.cat.id", sep='\t', header=None)
N = len(set(samples[0]))
print(N,"samples are found.")
# crete a UKBB -> {0,..,N} mapping
samples = samples.set_index( samples[0])
samples[0] = range(N)

print("Loading as pandas dataframe...", end='')
df = pd.read_csv("IBD_GRM_t10.grm.gz", compression="gzip", header=None, sep='\t')
print("Creating the sparse matrix...", end='')
IBD_matrix = coo_matrix( (df[3], (samples.loc[ df[0], 0], samples.loc[ df[1], 0])), shape=(N,N) )
print("Erase the dataframe and free up memory...")
del df 
print("Calculate the eigenpairs..", end='')
eigvals, eigvec = eigsh( IBD_matrix + IBD_matrix.transpose(), 40)

fieldnames = ["sample"]
for i in range(1,101):
    fieldnames.append("PC" + str(i))

output_file = csv.DictWriter(open("IBD_PCs.txt", 'w'), fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
output_file.writeheader()

for i in range(0,N):

    r = {}
    r["sample"] = samples.index.values[i]

    for j in range(0,100):
        fieldname = "PC" + str(j + 1)
        r[fieldname] = eigvec[i][j]

    output_file.writerow(r)


