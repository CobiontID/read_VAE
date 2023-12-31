# Utilities for VAE plots

## Plotting: vae_draw.py
Generate a 2D scatter plot from a two-dimensional numpy array or text file containing x and y coordinates

### Arguments
#### Inputs and outputs
- `--zfile`: Path to a file containing an array of 2D coordinates. The input may be a tab delimited text file with no header or .npy.
  - `--npy`: Specify if input is in .npy format, otherwise text will be assumed (default = `F`)
- `--labels`: Path to a file containing integer labels for each row in the coordinate file (one integer per tow). The labels, which may be generated with the labelling utilities described below, are used to assign a colour to each category. If no labels are provided, the data will be shaded based on density.
- `--edges`: Specify labels for colourbar bin edges (Optional).
- `--fignames`: Prefix for the filename where the output figure is stored (default = "figure").
- `--outdir`: Output directory path (optional).
#### Plotting options
- `--legend_y_label`: Set the label for the colourbar legend
- `--canvas_size`: Set the size of the output (choose between 500, 1000 and 1500, default 1000)
- `--categorical`: Enable categorical colour maps by setting to `T` (default = `F`). This setting is suitable to highlight contrasts between different taxonomic classes. If categorical maps are not enabled, slices from a continuous colour map will be used.
- `--cmap`: Specify a custom categorical colour map as RGB, e.g. `"#d60000 #8c3bff #018700"` Useful for fine-tuning the colourscheme for a specific dataset.
- `--superimpose`: If set to `T`, all points with non-zero labels will be plotted *over* points labelled zero rather than shaded (default = `F`). This allows labelled points to be seen more easily when they are sparse.
- `--remove`: If set to `T` and `--labels` are provided, all points with non-zero labels will be filtered out (default = `F`). Allows reads that were not assigned to a class to be visualised by themselves, which can be useful for identifying false negatives (the default label generated by `category_labels.py` for reads that do not appear in any of the provided lists is 0).

### A note on accessibility
The default colour maps used in this tool to indicate different quantiles and categories are "Viridis" (matplotlib) and "Colorblind" ([Bokeh](https://docs.bokeh.org/en/latest/docs/reference/palettes.html)), which are designed to be colorblind safe. However, the latter provides a limited number of categories (eight), and larger accessible maps may not provide sufficient contrast.

Therefore limiting or consolidating classes may be advisable. In some cases, specifying a custom colour map (see above) may also be helpful. If more than eight classes are provided with the `--categorical` setting enabled, the code will print a warning and fall back on [glasbey_hv](https://colorcet.holoviz.org/user_guide/Categorical.html).

## Labelling
- **category_labels.py**: Generate a list of integer labels based on one of more lists of sequence identifiers. Can be used to colour read plots according to taxonomic classification from e.g. Kraken.
- **category_labels_from_cont.py**: Generate bins based on n quantiles from a continuous vector
- **category_labels_from_cont_manual.py**: Generate bins based on a predefined list of edges. Useful for e.g. slicing FastK median vectors according to the minima in a GenomeScope profile.

## Get sequences from local peaks: peak_reads.py

Detects local peaks in 2D space and samples a specified number of sequences near each peak. Returns a fasta file containing the sampled sequences, and a plot showing their coordinates. The fasta may then be used to query a reference database and identify components of the sample. Currently requires each read sequence to be in a single line in the fasta file, which may be gzip-ed.

### Arguments
`--coords`: Path to a file containing an array of 2D coordinates (required).
- `-- ids`: Path to file containing sequence identifiers (required).
- `--fasta`: Path to the fasta file containing the read sequences (required).
- `--sampleid`: Name of sample, used to name output files
- `--n_bins`: Number of bins to use in 2D histogram (default 100)
- `--n_reads`: Number of reads per peak  to extract (default 2)
- `--outdir`: Path to output directory

### Example output
<img src="https://github.com/CobiontID/read_VAE/assets/10507101/b5765e71-b213-4821-9be7-4c7b9dd740e1" width=500>
