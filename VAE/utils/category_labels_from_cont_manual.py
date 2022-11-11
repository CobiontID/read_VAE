#!/usr/bin/env python
import numpy as np
import argparse

parser = argparse.ArgumentParser(
    description='Generate categorical integer labels by sequence ids from a continuous vector. Bin edges are specified manually')

parser.add_argument(
    "--idfiles", help="One or more files with sequence identifiers (one per line)", required=True)
parser.add_argument(
    "--feature", help="File with continuous vector to be binned (one per line)", required=True)
parser.add_argument("--labelled", help="Destination for output file with labels.", required=True)
parser.add_argument("--bins", help="List of bin edges.", required=True)


args = parser.parse_args()
print(args)

idfiles = args.idfiles.split()
feature = args.feature.split()
outfile = args.labelled
bins = [int(i) for i in args.bins.split()]

print(idfiles)

ids = [open(f).read().split() for f in idfiles]
ids = [y for x in ids for y in x]

cont = [open(f).read().split() for f in feature]
cont = np.array([float(y) for x in cont for y in x])


bins = [0.] + bins + [np.max(cont)]
digitized = np.digitize(cont, bins, right=True)

np.savetxt(outfile, digitized, fmt='%i', delimiter="\n")
