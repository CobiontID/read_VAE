import numpy as np
import argparse

parser = argparse.ArgumentParser(
    description='Generate categorical integer labels by sequence from a continuous vector. N specifies the number of bins,')

parser.add_argument(
    "--idfiles", help="One or more files with sequence identifiers (one per line). Obsolete, retained for compatibility with old Snakefile.")
parser.add_argument(
    "--feature", help="File with continuous vector to be binned (one per line)", required=True)
parser.add_argument("--labelled", help="Destination for file with labels.", required=True)
parser.add_argument("--n", help="Number of bins.", default=10)
parser.add_argument("--linear", help="Use linear binning instead of quantiles", default="F")
parser.add_argument("--log", help="Use logarithmic binning instead of quantiles", default="F")

args = parser.parse_args()
print(args)

feature = args.feature.split()
outfile = args.labelled
n = int(args.n)
linear = args.linear
log = args.log

if (log == "T") & (linear == "T"):
    linear = "F"
    print("Can't set logarithmic and linear binning simultaneously, defaulting to log1p.")

cont = [open(f).read().split() for f in feature]
cont = np.array([float(y) for x in cont for y in x])

if linear == "T":
    print("Using linear binning")
    bins = np.linspace(np.min(cont), np.max(cont), n+1)
    digitized = np.digitize(cont, bins[:-1]) - 1
    edges = bins
elif log == "T":
    print("Using logarithmic binning")
    cont = np.log1p(cont)
    bins = np.linspace(np.min(cont), np.max(cont), n+1)
    digitized = np.digitize(cont, bins[:-1]) - 1
    edges = [np.exp(i)-1 for i in bins]
else:
    bins = [np.quantile(cont, i) for i in np.linspace(0, 1, n+1)]
    # Remove duplicate bin edges to ensure continuous numbering
    bins = sorted(list(set(bins)))
    edges = bins
    if len(bins) < n+1:
        print(f"Note: Collapsed {(n+1)-len(bins)} quantile bin(s). Kept {len(bins)-1}.")
    digitized = np.digitize(cont, bins[:-1]) - 1
    print(f'Counts per label: {" ".join(str(i) for i in np.bincount(digitized))}')

# Save labels
np.savetxt(outfile, digitized, fmt='%i', delimiter="\n")

# Save bin edges
np.savetxt(f"{outfile}.edges", edges, fmt='%2f', delimiter="\n")
