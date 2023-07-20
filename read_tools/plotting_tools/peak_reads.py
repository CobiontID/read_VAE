#!/usr/bin/env python
# cw21@sanger.ac.uk

import itertools
import argparse

parser = argparse.ArgumentParser(
    description='Get read sequences from local peaks in the latent space of the VAE.')

parser.add_argument(
    "--coords", help="2D array with 2D coordinates", required=True)
parser.add_argument(
    "--ids", help="Path to file with sequence ids (in same order as coordinates)", required=True)
parser.add_argument(
    "--fasta", help="Path to fasta file containing sequences to be extracted", required=True)
parser.add_argument("--sampleid", help="Name of sample",
                    default="sample", required=False)
parser.add_argument("--n_bins", help="Number of bins to use in 2D histogram",
                    type=int, default=100, required=False)
parser.add_argument("--n_reads", help="Number of reads per bin to extract",
                    type=int, default=2, required=False)
parser.add_argument("--outdir", help="Output directory",
                    default="", required=False)


args = parser.parse_args()
print(args)

coords = args.coords
ids = args.ids
fasta = args.fasta
sampleid = args.sampleid
n_bins = args.n_bins
outdir = args.outdir
n_reads = args.n_reads

#TODO: Check files exist

## Imports 
import subprocess
import numpy as np
import pandas as pd

from collections import defaultdict
from findpeaks import findpeaks
from holoviews import dim, opts
import holoviews as hv
from holoviews.operation.datashader import datashade
from datashader import transfer_functions as tf
hv.extension('bokeh')

# Load data

if coords.split(".")[-1] == "npy":
    a = np.load(coords)
else:
    a = np.loadtxt(coords, dtype=np.float32)

read_ids = np.array([i.strip() for i in open(ids)])


# Do binning
def make_histogram_bins(a, n_bins):
    h = np.histogram2d(a[:,0], a[:,1], bins=n_bins)

    bins_y = np.searchsorted(h[2], a[:,1])
    bins_x = np.searchsorted(h[1], a[:,0])

    d = defaultdict(list)

    for i in range(len(a)):
        d[(bins_x[i], bins_y[i])].append(i)

    return h, d

h, d = make_histogram_bins(a, n_bins)

eq = tf.eq_hist(h[0], nbins=n_bins)
print(type(eq))


# Get cells with peaks
def find_peak_cells(eq):
# Initialize peak finder
    fp = findpeaks(method='topology', window=20)
    results = fp.fit(eq)
    where = np.where(results['Xdetect'] > 0)
    peak_rows = [d[(i, j)][:n_reads] for i, j in zip(where[0], where[1])]
    return list(itertools.chain(*peak_rows))

flattened = find_peak_cells(eq)

def get_seq_file_batch(identifiers, fasta, sampleid):
    """Get batch of sequences"""
    grep = "zgrep" if fasta.split(".")[-1] == "gz" else "grep"
    id_string = "\n".join(identifiers)
    with open("tmp.ids", "w") as f:
        f.write(id_string)
    seq_cmd = f'{grep} -A1 -f tmp.ids {fasta} > {outdir}peaks.{sampleid}.fa'
    seq = subprocess.run(seq_cmd, capture_output=True, shell=True)
    return 1 if seq.returncode != 0 else 0


got_fasta = get_seq_file_batch(read_ids[flattened], fasta, sampleid)

if got_fasta == 1:
    print("Non-zero subprocess return code.")
    print(read_ids[flattened])

read_table = pd.DataFrame(data=np.hstack([read_ids[flattened][:, None], a[flattened]]), columns=["id", "x", "y"])
read_table.to_csv(f"{outdir}{sampleid}.peaks_reads_coords.tsv", sep="\t")

# Draw coordinates of points
samples = hv.Points(a[flattened]).opts(color="blue", size=5, show_grid=True)
points = datashade(hv.Points(a), cmap="grey").opts(width=800, height=800)

hv.save(points*samples, f"{outdir}{sampleid}.peaks.png", backend="bokeh", fmt='png')

