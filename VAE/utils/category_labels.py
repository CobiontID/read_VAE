#!/usr/bin/env python
import numpy as np
import argparse
from itertools import combinations

parser = argparse.ArgumentParser(
    description='Generate categorical integer labels by sequence ids. Unlabelled identifiers are assigned 0, the nth set of sequences is assigned the label n. Identifiers found in multiple sets are assigned to their own additional category if present.')

parser.add_argument(
    "--idfiles", help="One or more files with sequence identifiers (one per line)", required=True)
parser.add_argument(
    "--bins", help="One or more files listing identifiers of sequences belonging to a given bin (one per line)", required=True)
parser.add_argument("--labelled", help="Destination for file with labels.", required=True)

args = parser.parse_args()
print(args)

idfiles = args.idfiles.split()
binfiles = args.bins.split()
outfile = args.labelled

print(idfiles)
ids = [open(f).read().split() for f in idfiles]
ids = [y for x in ids for y in x]
mask = np.array(len(ids)*[0], dtype="int32")
print(len(mask))

all_sets = []
for i, cat in enumerate(binfiles):
    seq_set = set([j.strip() for j in open(cat)])
    all_sets.append(seq_set)
    np.put(mask, np.where([seqid in seq_set for seqid in ids]), [i+1])

# Add bin for identifiers that appear in multiple sets.
nt = lambda a, b: all_sets[a].intersection(all_sets[b])
in_multiple = set().union(*[nt(*j) for j in combinations(range(i+1), 2)])

if len(in_multiple) > 0:
    np.put(mask, np.where([seqid in in_multiple for seqid in ids]), [i+2])
    print("Adding extra bin containing intersect of sets: {}".format(i+2))

np.savetxt(outfile, mask, fmt='%i', delimiter="\n")
