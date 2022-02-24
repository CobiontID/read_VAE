# Generate two-dimensional representations and visualisations for read k-mer composition

## Running the pipeline

### Setup
This set-up assumes that the script will be run on an LSF cluster.

- Set up the configuration in <a href="https://github.com/CobiontID/kmer_decomposition/blob/main/readviz_pipeline/config.yml">config.yml</a>
  - Set up the sample to be run: 
     - `sample_id`: The sample identifier (used to name output files)
     - `species_name`: The full species name 
     - `read_file`: Path to the file containing the reads to be analysed
     - `fastk_tab`: Path to FastK .ktab file. Optional. If provided, FastK will profile the reads against the existing table. If not using, set to `None`
     
   - Base configuration:
     - `user_group`: User group to be used for LSF
     - `tetra_count_path`: Path to the k-mer counter executable
     - `fastk_path`: Base path to FastK and ProfMedianAll executables
     - `un_count_path`: Path to unique k-mer counter executable
     - `hexamer_path`: Path to hexsum executable (for estimating coding density)
     - `vae_path`: Path to Variational Autencoder-based k-mer decomposition script
     - `draw_vae_path`: Path to script for plotting decomposed read k-mers
     - `cat_labeller_path`: Path to script for generating bins for reads from a continuous vector
     - `hextable_path`: Path to reference table for Hexamer/Hexsum
     - `conda_tf`: The conda environment used for the decomposition and visualisation steps. Required packages are specified in in <a href="https://github.com/CobiontID/kmer_decomposition/blob/main/env_kmerviz.yaml">env_kmerviz.yaml</a>.

### Run
In a conda environment with <a href="https://snakemake.readthedocs.io/en/stable/">Snakemake</a> installed, run `Snakemake`

