# Workflow to generate two-dimensional representations and visualisations for contig composition

## Software requirements
- Install dependencies from source:
  - <a href="https://github.com/CobiontID/kmer-counter">kmer-counter</a> (required)
  - <a href="https://github.com/CobiontID/unique-kmer-counts">unique-kmer-counter</a>
  - <a href="https://github.com/CobiontID/fastk-medians">fastk-medians</a> (install `FastK` and `ProfMedianAll`)
  - <a href="https://github.com/richarddurbin/hexamer">Hexamer</a> (you will need `cds.worm.hex`).

- Set up a Conda environment according with <a href="https://github.com/CobiontID/read_VAE/blob/main/env_kmerviz.yaml">this</a> configuration (see <a href="https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file">here</a> for instructions).

- You will also need `Select_contigs_reduced_multi.py`from this repository and, if using scaffold contact information, `hic_links.py`.

## Running the Snakemake pipeline
The included example set-up assumes that the script will be run on an LSF cluster. For instructions how to run the steps individually, see below.

Note that the contig/scaffold workflow currently does not use a VAE. The number of sequences typically does not require it, and results can be inconsistent on small sets of sequences and require manual parameter tuning. For large numbers of very short contigs, consider using the [workflow designed for long reads](https://github.com/CobiontID/read_VAE/blob/main/read_tools/Workflow.md) instead. To avoid problems with global distances, there is a PCA option (described below). 

### Setup
This set-up assumes that the script will be run on an LSF cluster.

- Set up the configuration in <a href="https://github.com/CobiontID/read_VAE/blob/main/contig_tools/config.yml">config.yml</a>
  - Set up the sample to be run: 
     - `sample_id`: The sample identifier (used to name output files)
     - `species_name`: The full species name 
     - `contig_file`: Path to the fasta (or fasta.gz) containing the contig sequences
     - `contig_noseq`: Path to hifiasm .noseq.gfa to calculate coverage, optional (if not using hifiasm, set to None)
     - `assembler`: Assembler, e.g. hifiasm
     - `seq_type`: Type of sequence input for labelling purposes, e.g. "p_ctg" or "scaffolds"
     - `collapse_kmers`: Specify if k-mers should be canonicalised (collapsed) or not (uncollapsed)
     - `pair_file`: Pair file with scaffold contacts from SALSA or YaHs, optional (if not using Hi-C, set None). The relevant file will commonly be named `alignments_sorted.txt`.
     - `size_file`: File with scaffold sizes (if not using Hi-C, set None). The file will commonly be named `out_scaffolds_final.fa.chrom.sizes` or similar, and contains two columns (scaffold names and scaffold lengths).

  - Base configuration:
     - `user_group`: User group to be used for LSF
     - `tetra_count_path`: Path to the k-mer counter executable
     - `fastk_path`: Base path to FastK and ProfMedianAll executables
     - `un_count_path`: Path to unique k-mer counter executable
     - `hexamer_path`: Path to hexamer executable (for estimating coding density)
     - `hextable_path`: Path to reference table for Hexamer
     - `reduced_plot_path`: Path to <a href="https://github.com/CobiontID/read_VAE/blob/main/contig_tools/scripts/Select_contigs_reduced_multi.py">script</a> to decompose k-mer counts and draw annotated contig selection graph.
     - `hic_link_path`: Path to <a href="https://github.com/CobiontID/read_VAE/blob/main/contig_tools/scripts/hic_links.py">script</a> to read SALSA pairs file and generate connectivity annotations.
     - `conda_tf`: The conda environment used for the decomposition and visualisation steps.

### Run
In a conda environment with <a href="https://snakemake.readthedocs.io/en/stable/">Snakemake</a> installed, run `Snakemake`

### Run steps individually on the commandline
- Follow steps 1-4 outlined [here](https://github.com/CobiontID/read_VAE/blob/main/read_tools/Workflow.md) using the contig or scaffold sequences as input. Set k for `unique-kmers` to 15 to accommodate longer sequences.
- Extract coverage if using a suitable hifiasm assembly as input:
`grep "S.*ptg" contigs_sampleid.noseq.gfa | sed -n 's/.*rd:i:\([0-9]\)/\1/p' > sampleid.coverage.txt`
- If applicable, get connections between scaffolds
- `python hic_links.py --pairfile alignments_sorted.txt --sizefile out_scaffolds_final.fa.chrom.sizes --outfile sampleid.connections.npy --outconn  isconnected.sampleid.txt --outconnbp isconnected.sampleid.normbp.txt`
- Generate the plot:
  - Using a set of contigs as an example:`python Select_contigs_reduced_multi.py --seqtype p_ctg --infile sampleid.p_ctg.tetra.collapsed.npy --outfile sampleid_p_ctg_multi_select.html --seqidfile sampleid.p_ctg.ids.txt --annotfiles "sampleid.p_ctg.hexsum sampleid.median_31mer.txt sampleid.15_mers.txt sampleid.coverage.txt" --annotnames "Hexamer FastK Unique_15mers Coverage" --speciesname sampleid`
  - `--annotnames`contains a list of labels for each annotation vector that is supplied, and `--annotfiles` contains the corresponding list of files. If multiple files and labels are supplied, the lists must be encloded in quotation marks.
  - To include scaffold connections, add "isconnected.sampleid.txt isconnected.sampleid.normbp.txt" and "Is_Connected Connections_Base" to the lists. If no coverage information is available, omit "sampleid.coverage.txt" and "Coverage", and so on.
  - By default, the code will use UMAP to represent the tetranucleotide counts, but there is a commandline argument to use PCA instead (run the script with --help for details).
