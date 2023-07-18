# Generate two-dimensional representations and visualisations for read k-mer composition

## Setup
Install required tools from source:
- <a href="https://github.com/CobiontID/kmer-counter">kmer-counter</a>
- <a href="https://github.com/CobiontID/unique-kmer-counts">unique-kmer-counter</a>
- <a href="https://github.com/CobiontID/fastk-medians">fastk-medians</a>
- <a href="https://github.com/richarddurbin/hexamer">Hexamer</a>.

Set up a Conda environment according to https://github.com/CobiontID/read_VAE/blob/main/env_kmerviz.yaml

### Running the Snakemake pipeline

The included example set-up assumes that the script will be run on an LSF cluster.

- Set up the configuration in <a href="https://github.com/CobiontID/read_VAE/blob/main/read_tools/config.yml">config.yml</a>
  - Set up the sample to be run: 
     - `sample_id`: The sample identifier (used to name output files)
     - `species_name`: The full species name (this field is included for record-keeping purposes)
     - `read_file`: Path to the fasta file containing the reads to be analysed.
     - `fastk_tab`: Path to FastK .ktab file. Optional. If provided, FastK will profile the reads against the existing table. If not using, set to `None`. The default k-mer length is k = 31 as of Feb 2022.
     
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

#### Run
In a conda environment with <a href="https://snakemake.readthedocs.io/en/stable/">Snakemake</a> installed, run `Snakemake`. Depending on the input data, the memory requirements may need to be adjusted in the Snakefile.

### Steps to run without Snakemake

Count canonicalized tetranucleotides:
```
/software/kmer_counter/target/release/kmer-counter --file reads_sampleid.fa.gz --collapse 1 --ids sampleid.reads.ids.txt --klength 4 --out sampleid.reads.tetra.collapsed.npy
```

Count unique 8-mers per base:
```
/software/unique_kmers/target/release/unique-kmers --klength 8 --file reads_sampleid.fa.gz --out sampleid.8mer.txt'

```

Calculate estimated coding density with Hexamer (Piping is only necessary if the input fasta file is gzip-ed):
```
zcat reads_sampleid.fa.gz | /software/team301/hexamer/hexamer -T 20 -S /software/team301/user/cw21/hexamer/cds.worm.hex - | awk '{ print $3/$2 }' > sampleid.reads.hexsum
```

Get median number of times each k-mer of size 31 in each read occurs across the dataset:
```
cp -n reads_sampleid.fa.gz ./fastk/profile/reads_sampleid.fa.gz
/software/fastk/FastK -k31 -T8 -M14 -p ./fastk/profile/reads_sampleid.fa.gz
/software/fastk/ProfMedianAll ./fastk/profile/*prof > sampleid.median_31mer.txt
```
Run the VAE:
```
conda activate /software/conda/kmerviz
outdir=./vae/sampleid/
python ./VAE/vae.py --countfile sampleid.reads.tetra.collapsed.npy --fignames sampleid --kl 0.0025 --epochs 15 --outdir ${outdir}
```

Generate colour-coded plots. First, assign labels based on k-mer statistics:

```
python ./plotting_tools/category_labels_from_cont.py --feature sampleid.reads.hexsum --labelled hexsum.binned.sampleid --n 10
python ./plotting_tools/category_labels_from_cont.py --feature sampleid.median_31mer.txt --labelled 31mer.binned.sampleid --n 10
python ./plotting_tools/category_labels_from_cont.py --feature sampleid.8mer.txt --labelled 8mer.binned.sampleid --n 10
```

Draw the plots:
```
python ./VAE/vae_draw.py --zfile ${outdir}/sampleid.vae.out.2d.0 --outdir ${outdir}/ --fignames sampleid_31mer --labels 31mer.binned.sampleid --edges 31mer.binned.sampleid.edges --legend_y_label "Median 31-mer count"
python ./VAE/vae_draw.py --zfile ${outdir}/sampleid.vae.out.2d.0 --outdir ${outdir}/ --fignames sampleid_hexamer --labels hexsum.binned.sampleid --edges hexsum.binned.sampleid.edges --legend_y_label "Hexamer"
python ./VAE/vae_draw.py --zfile ${outdir}/sampleid.vae.out.2d.0 --outdir ${outdir}/ --fignames sampleid_8mer --labels 8mer.binned.sampleid --edges 8mer.binned.sampleid.edges --legend_y_label "Unique 8-mers/base"
```
