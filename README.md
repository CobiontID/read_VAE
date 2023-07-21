# kmer_decomposition

The documentation in this repository is currently still under construction.

## Workflows
### <a href="https://github.com/CobiontID/read_VAE/tree/main/read_tools">Read k-mer decomposition and visualisation</a>

Tallies k-mers in a read set, reduce to two dimensions and visualise read clusters (defaults to tetranucleotides). Annotates read plots with additional sequence features, such as estimated coding density, approximate k-mer coverage and sequence k-mer diversity.

#### Example data set: _Erannis defoliaria_
Decomposed read tetranucleotides from _Erannis defoliaria_ indicate the presence of bacteria in the sample (top). The reads are coloured by estimated coding density.

<img src="https://github.com/CobiontID/read_VAE/blob/main/ilEraDefo1_hexamer.2d_plot_labelled.png" width=500>

### <a href="https://github.com/CobiontID/read_VAE/tree/main/contig_tools">Contig/scaffold k-mer decomposition and visualisation</a>

As with the reads, tallies and reduces tetranucleotide composition to two dimensions and plots with annotations. In addition to estimated coding density and k-mer diversity, FastK provides a measure of repetitiveness, and coverage for primary Hifiasm assemblies can be extracted and used to annotate the plots. A selection tool allows sequences that are of interest to be selected and downloaded with their annotations. Where Hi-C data are available, a SALSA or YaHs pair file may also be provided to annotate plots with scaffold connectivity information. Take a look at an interactive version of the plot [here](https://cobiontid.github.io/examples.html#scaffold-tetranucleotide-visualisation).

#### Example data set: Hylocomiadelphus triquetrus
![image](https://user-images.githubusercontent.com/10507101/133108115-a3dbe6af-a602-47d9-a56c-27887d464084.png)


## Tools

### <a href="https://github.com/CobiontID/read_VAE/tree/main/contig_tools/VAE">Variational Autoencoder for k-mer decomposition</a>
#### vae.py
Read k-mer counts are reduced to two dimensions following the method of <a href="https://arxiv.org/abs/1312.6114">Kingma and Welling (2013)</a>. Outputs two-dimensional representation of the read set and a basic plot.

#### [Plotting tools for reads](https://github.com/CobiontID/read_VAE/tree/main/read_tools/plotting_tools)
Generate colour-coded plots of 2D representations learned by the VAE.

### [Visualisations for contigs](https://github.com/CobiontID/read_VAE/tree/main/contig_tools)
Workflow and utilities to generate interactive HTML file of decomposed tetranucleotide plots with binned annotations.

### Standalone tools used in workflows
#### [kmer-counter](https://github.com/CobiontID/kmer-counter)
Counts the number of occurences of each k-mer of size k for each record in a fasta file of nucleotide sequences (canonicalised or non-canonicalised). Implemented in Rust, runs approximately ten times faster than the equivalent code in Python.

#### [unique-kmer-counts](https://github.com/CobiontID/unique-kmer-counts)
Count the number of distinct k-mers of size k for each record in a fasta files of nucleotide sequences, and divide by sequence length. Implemented in Rust.

#### hexamer
Estimates the coding density using the sum of lengths of putative coding sequences divided by sequence length. The cobiont pipelines previously used a modified version of the old hexamer code. The relevant functionality is now available in an updated version of hexamer from
https://github.com/richarddurbin/hexamer (to extract the estimated density, pipe stdout to `awk '{ print $3/$2}'`)

#### [fastk-medians](https://github.com/CobiontID/fastk-medians)
Calculates the median number of times each k-mer of size k (in this case k = 31) occurs across the whole set of sequences. Provides an approximation of coverage for reads (provided they are not highly repetitive), or repetitiveness for assembled contigs or scaffolds.
