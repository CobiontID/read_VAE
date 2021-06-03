import numpy as np
import argparse

parser = argparse.ArgumentParser(
    description='Generate categorical integer labels by sequence ids from a continuous vector. N specifies the number of quantiles to bin into.')

parser.add_argument(
    "--idfiles", help="One or more files with sequence identifiers (one per line)")
parser.add_argument(
    "--feature", help="File with continuous vector to be binned (one per line)")
parser.add_argument("--labelled", help="Destination for file with labels.")
parser.add_argument("--n", help="Number of bins.")

args = parser.parse_args()
print(args)

idfiles = args.idfiles.split()
feature = args.feature.split()
outfile = args.labelled
n = int(args.n)

print(idfiles)

ids = [open(f).read().split() for f in idfiles]
ids = [y for x in ids for y in x]

cont = [open(f).read().split() for f in feature]
cont = np.array([float(y) for x in cont for y in x])

bins = [np.quantile(cont, i) for i in np.linspace(0, 1, n+1)]
digitized = np.digitize(cont, bins[:-1]) - 1

np.savetxt(outfile, digitized, fmt='%i', delimiter="\n")
