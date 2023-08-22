# VAE tools for read data

## <a href="https://github.com/CobiontID/read_VAE/blob/main/read_tools/Workflow.md">How to run the VAE workflow</a>
Example Snakemake workflow (configured for use on an LSF cluster), and instructions to run each step individually. Count tetranucleotides for read set, encode into two dimensions, and render plots annotated with estimated coding density, median number of 31-mers across the read set (approximates coverage), and the number of unique k-mers per base.

## [VAE](https://github.com/CobiontID/read_VAE/tree/main/read_tools/VAE)
Run the VAE model on a set of tetranucleotide counts, and render a plot of the latent space.

## [Tools for plotting 2D representations](https://github.com/CobiontID/read_VAE/tree/main/read_tools/plotting_tools)
Overview of scripts to plot 2D coordinates learned by the VAE with Datashader, with optional colour-coded annotations (these may be continuous variables or category labels). Also includes code to sample sequences from near local peaks in the latent space.

## [Interactive dashboard](https://github.com/CobiontID/read_VAE/tree/main/read_tools/dashboard)
Filter and explore 2D representations annotated with information about coding density and coverage, and class labels. Filter by class and/or coverage. Select sequences of interest to launch blast queries, display the distributions of coverage and coding density, and download statistics.
