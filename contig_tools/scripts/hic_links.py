#!/usr/bin/env python
# cw21@sanger.ac.uk

import argparse

parser = argparse.ArgumentParser(
    description='Generate scaffold contact matrix')
parser.add_argument("--pairfile", help="input file with pairs", required=True)
parser.add_argument("--sizefile", help="input scaffold size file", required=True)
parser.add_argument("--outfile", help="Output .npy", required=True)
parser.add_argument("--outconn", help="Output .txt with bool", default = "None")
parser.add_argument("--outconnbp", help="Output .txt with bool", default = "None")

args = parser.parse_args()
print(args)

# required
infile = args.pairfile
outfile = args.outfile
sizefile = args.sizefile
# optional
connect_out = args.outconn
connect_bp = args.outconnbp

import numpy as np

# get scaffold number/lengths
scaff_lengths = [int(line.strip().split("\t")[1]) for line in open(sizefile)]
n_scaffs = len(scaff_lengths)

a = np.zeros((n_scaffs, n_scaffs))
links = open(infile, "r")

for line in links:
    line = line.split("\t")
    first = int(line[1].split("_")[1]) - 1
    second = int(line[5].split("_")[1]) - 1
    if second != first:
        a[first,second] += 1
        a[second, first] += 1
    else:
        a[first,first] += 1

# Save .npy with matrix
np.save(outfile, a)

# Write file encoding whether scaffolds have >0 connections if requested
if connect_out != "None":
    np.fill_diagonal(a, 0)
    is_conn = "\n".join([str((np.sum(a[i,:]) > 0).astype(int)) for i in range(len(a))])
    with open(connect_out, "w") as f:
        f.write(is_conn)
        f.close()
        
# Write file with (sum of connections/scaffold length) if requested
if connect_out != "None":
    np.fill_diagonal(a, 0)
    conn_base = "\n".join([str((np.sum(a[i,:]).astype(int)/scaff_lengths[i])) for i in range(len(a))])
    with open(connect_bp, "w") as f:
        f.write(conn_base)
        f.close()
