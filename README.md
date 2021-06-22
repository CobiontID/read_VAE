# kmer_decomposition

## Workflows
### <a href="https://github.com/CobiontID/kmer_decomposition/tree/main/readviz_pipeline">Read k-mer decomposition and visualisation</a>

### Contig/scaffold k-mer decomposition and visualisation
Tallies k-mers in a read set, reduce to two dimensions and visualise read clusters (defaults to tetranucleotides). Annotates read plots with additional sequence features, such as estimated coding density, approximate k-mer coverage and sequence k-mer diversity.

#### Example data set: _Erannis defoliaria_
Decomposed read tetranucleotides from _Erannis defoliaria_ indicate the presence of bacteria in the sample (top). The reads are coloured by estimated coding density.

<img src="https://github.com/CobiontID/kmer_decomposition/blob/main/ilEraDefo1_hexamer.2d_plot_labelled.png" width=500>

## Tools
### <a href="https://github.com/CobiontID/kmer_decomposition/tree/main/kmer-counter">kmer-counter</a>
Count the number of occurences of each k-mer of size k for each record in a fasta file of nucleotide sequences (canonicalised or non-canonicalised). Implemented in Rust, runs approximately ten times faster than the equivalent code in Python.
### <a href="https://github.com/CobiontID/kmer_decomposition/tree/main/unique-kmer-counts">unique-kmer-counter</a>
Count the number of distinct k-mers of size k for each record in a fasta files of nucleotide sequences, and divide by sequence length. Implemented in Rust.
