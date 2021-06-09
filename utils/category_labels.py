import numpy as np
import argparse

parser = argparse.ArgumentParser(
    description='Generate categorical integer labels by sequence ids. Unlabelled identifiers are assigned 0, the nth set of sequences is assigned the label n')

parser.add_argument(
    "--idfiles", help="One or more files with sequence identifiers (one per line)")
parser.add_argument(
    "--bins", help="One or more files listing identifiers of sequences belonging to a given bin (one per line)")
parser.add_argument("--labelled", help="Destination for file with labels.")

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
if len(set.intersection(*all_sets)) > 0:
    np.put(mask, np.where([seqid in set.intersection(*all_sets) for seqid in ids]), [i+2])
    print("Adding extra bin containing intersect of sets: {}".format(i+2))

np.savetxt(outfile, mask, fmt='%i', delimiter="\n")
