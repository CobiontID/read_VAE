# Generate two-dimensional representations and visualisations for contig composition

## Running the pipeline

### Setup
This set-up assumes that the script will be run on an LSF cluster.

- Set up the configuration in <a href="https://github.com/CobiontID/kmer_decomposition/blob/main/contigviz_pipeline/config.yml">config.yml</a>
  - Set up the sample to be run: 
     - `sample_id`: The sample identifier (used to name output files)
     - `species_name`: The full species name 
     - `contig_file`: Path to the fasta (or fasta.gz) containing the contig sequences
     - `contig_noseq`: Path to hifiasm .noseq.gfa to calculate coverage (if not using hifiasm, set to None)
     - `assembler`: Assembler, e.g. hifiasm
     - `seq_type`: Type of sequence input for labelling purposes, e.g. "p_ctg" or "scaffolds"
     - `collapse_kmers`: Specify if k-mers should be canonicalised (collapsed) or not (uncollapsed)
     - `pair_file`: Pair file with scaffold contacts from SALSA, optional (if not using Hi-C, set None)
     - `size_file`: File with scaffold sizes from SALSA (if not using Hi-C, set None)

  - Base configuration:
     - `user_group`: User group to be used for LSF
     - `tetra_count_path`: Path to the k-mer counter executable
     - `fastk_path`: Base path to FastK and ProfMedianAll executables
     - `un_count_path`: Path to unique k-mer counter executable
     - `hexamer_path`: Path to hexamer executable (for estimating coding density)
     - `hextable_path`: Path to reference table for Hexamer
     - `reduced_plot_path`: Path to <a href="https://github.com/CobiontID/kmer_decomposition/blob/main/draw_contigs/Select_contigs_reduced_multi.py">script</a> to decompose k-mer counts and draw annotated contig selection graph.
     - `hic_link_path`: Path to <a href="https://github.com/CobiontID/kmer_decomposition/blob/main/Hi-C/utils/hic_links.py">script</a> to read SALSA pairs file and generate connectivity annotations.
     - `conda_tf`: The conda environment used for the decomposition and visualisation steps. Required packages are specified in in <a href="https://github.com/CobiontID/kmer_decomposition/blob/main/env_kmerviz.yaml">env_kmerviz.yaml</a>.

### Run
In a conda environment with <a href="https://snakemake.readthedocs.io/en/stable/">Snakemake</a> installed, run `Snakemake`
