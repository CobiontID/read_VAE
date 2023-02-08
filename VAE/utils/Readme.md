# Utilities for VAE plots

## Plotting
### vae_draw.py
Generate a 2D scatter plot from a two-dimensional numpy array or text file containing x and y coordinates

#### Arguments
##### Inputs and outputs
- `--zfile`: Path to a file containing an array of 2D coordinates. The input may be a tab delimited text file with no header or .npy.
  - `--npy`: Specify if input is in .npy format, otherwise text will be assumed (default = `F`)
- `--labels`: Path to a file containing integer labels for each row in the coordinate file (one integer per tow). The labels, which may be generated with the labelling utilities described below, are used to assign a colour to each category. If no labels are provided, the data will be shaded based on density.
- `--fignames`: Prefix for the filename where the output figure is stored (default = "figure").
- `--outdir`: Output directory path (optional).
##### Plotting options
- `--canvas_size`: Set the size of the output (choose between 500, 1000 and 1500, default 1500)
- `--categorical`: Enable categorical colour maps by setting to `T` (default = `F`). This setting is suitable to highlight contrasts between different taxonomic classes. If categorical maps are not enabled, slices from a continuous colour map will be used.
- `--cmap`: Specify a custom categorical colour map as RGB, e.g. `"#d60000 #8c3bff #018700"` Useful for fine-tuning the colourscheme for a specific dataset.
- `--superimpose`: If set to `T`, all points with non-zero labels will be plotted *over* points labelled zero rather than shaded (default = `F`). This allows labelled points to be seen more easily when they are sparse.

##### A note on accessibility
The default colour maps used in this tool to indicate different quantiles and categories in the data are "Viridis" (matplotlib) and "Colorblind" ([Bokeh](https://docs.bokeh.org/en/latest/docs/reference/palettes.html)), which are designed to be colorblind safe. However, the latter provides a limited number of categories (eight), and larger accessible maps may not provide sufficient contrast.

Therefore limiting or consolidating classes may be advisable. In some cases, specifying a custom colour map (see above) may also be helpful. If more than eight classes are provided with the `--categorical` setting enabled, the code will print a warning and fall back on [glasbey_hv](https://colorcet.holoviz.org/user_guide/Categorical.html).

## Labelling
- **category_labels.py**: Generate a list of integer labels based on one of more lists of sequence identifiers. Can be used to colour read plots according to taxonomic classification from e.g. Kraken.
- **category_labels_from_cont.py**: Generate bins based on n quantiles from a continuous vector
- **category_labels_from_cont_manual.py**: Generate bins based on a predefined list of edges. Useful for e.g. slicing FastK median vectors according to the minima in a GenomeScope profile.
