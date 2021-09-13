# Select_contigs_reduced_multi.py
Generates HTML file of decomposed tetranucleotide plots with binned annotations. See <a href="https://cobiontid.github.io/examples/cbHylTriq8_scaffolds_multi_select.html">example</a>.
## Usage
Arguments:
- `--infile`: Path to .npy file with k-mer counts
- `--outfile`: Path of .html file to write results to
- `--seqidfile`: Path to file with sequence identifiers (one identifier per line, identifiers in same order as in the .npy file)
- `--annotfiles`: Path(s) to file(s) with annotations to be binned, surrounded by quotes if more than one, e.g. `"sample.hexsum sample.coverage"`. 
- `--annotnames`: Label(s) for annotation(s), corresponding to annotation files provided, surrounded by quotes if more than one e.g. `"Density Coverage"`. Number of labels must match number of annotation files.
- `----bins`: Number of bins for annotations, default value = 5.
- `--discretize`: Discretize annotation by quantile (default if no value supplied), or linear (value other than "quantile" supplied).
- `--speciesname`: String with species or sample name (displayed in plot label).
- `--pca`: If set to `T`, reduce k-mers to two dimensions with PCA. Default to `F`(reduction with UMAP). 
- `--seqtype`: String with type of sequence, displayed in plot label (e.g. contigs or scaffolds).
