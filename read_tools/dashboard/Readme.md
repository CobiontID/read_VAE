# Read VAE dashboard

To explore the structure of a given readset, interactively filtering the data and cross-referencing different sequence features and taxonomic labels can be useful. The selection tool allows sequences of interest to be retrieved and blasted, and statistics for the selection may be displayed and downloaded. A demo for an example dataset is available [here](https://huggingface.co/spaces/cc7740/read_VAE) (screenshots below).

<img src="https://github.com/CobiontID/read_VAE/assets/10507101/42d5a687-4af4-4900-9843-e1b597856c17" width=600><br>
<img src="https://github.com/CobiontID/read_VAE/assets/10507101/1d3005f0-3654-438c-b9cf-5e00155e8b60f" width=600>


See https://cobiontid.github.io/reads.html for more details and a walkthrough of the example.

## Required inputs
To display a dataset with the dashboard, the following are usually required: 2D coordinates for each read, the corresponding sequence identifiers, coding density estimates, and k-mer coverage estimates, and the fasta file from which the data were generated.

All of these are provided by the workflow described [here](https://github.com/CobiontID/read_VAE/tree/main/read_tools). Optionally, lists of identifiers belonging to different sequence classes may be supplied (these could be based on taxonomic classifications or read mapping, for example). For simpler visualisations, consider using the tools designed for [static plots](https://github.com/CobiontID/read_VAE/tree/main/read_tools/plotting_tools) instead.

## Running the dashboard 
The server can be started with:

```panel serve --port=$port --address="0.0.0.0" --allow-websocket-origin=localhost:$port --show read_dash.py --args --config example.yml```

Here,`$port` is the number of an open port on the local machine. `example.yml` is the path to the sample configuration file. See below for more information on configuration. For especially large datasets, it may be necessary to set `--session-token-expiration` to a larger value. If the configuration file is updated after the app has been launched, refreshing the browser window will prompt all of the data to be reloaded using the updated file paths or class lists, if supplied.

If the dashboard is run on a remote machine, connecting with the SSH extension in Visual Studio Code and entering the above in the terminal may be convenient, as port forwarding will be handled automatically. When prompted "Do you want Code to open the external website?", clicking "open" will open the dashboard in a browser window.

### Setting up the configuration file
A basic example configuration file can be found [here](https://github.com/CobiontID/read_VAE/blob/main/read_tools/dashboard/example.yml).

#### Basic configuration
- `tolid`: Sample identifier (this should match the identifier used to run the pipeline, unless the paths to the required files are manually set)
- `fasta`: Path to the fasta file containing the sequences of the reads. Required to launch blast queries. It is assumed that each sequence record has a single line, unless indexing is enabled.
- `basedir`: Path to the base directory containing outputs from the read VAE workflow (but see details on how to override defaults).
- `class_lists`: One or more files containing the identifiers of reads belong to different classes. Each file should be on a separate line, preceded by a dash. If left blank, class labels will not be used.
- `remote_blast`: Specify `True` to use NCBI's servers for blast queries. Querying a local server may give faster turnaround.

#### Options to manually configure paths
- `k`: The k-mer size used with FastK in the pipeline (used to retrieve the median k-mer count file). If not specified, k defaults to 31.
- `vae_path`: A custom path pointing to the 2D coordinates to be used for plotting. Input files may either be in text format, e.g. ilApaMono1_16.vae.out.2d.0, or in .npy format (the latter is preferable for large datasets, since it allows faster loading).
- `read_ids_path`: Custom path to read identifiers
- `fastk_path`: Custom path pointing to the file containing median read k-mer counts
  
#### Use indexed fasta files
For large readsets, it may be useful to index the fasta file containing the reads with `samtools faidx`, in order to retrieve sequences more quickly. Indexing also makes retrieving records more convenient when sequences are split across multiple lines. To make use of an index, include these lines:
```
indexed: True
samtools_path: /software/samtools
```
The path of the .fai file is currently assumed to be the same as for the corresponding fasta file, plus the suffix.

#### Experimental settings
Adding ```no_annot: True``` to the yml file allows the interface to be used with datasets for which coding density and coverage estimates are missing (filtering by coverage and colouring by coding density will not work in this case).

## Caveats
This code was written for research purposes and comes without included batteries. Although it has been tested extensively on various outputs generated by the read workflow in this repository, exception handling is sparse.

### Large datasets
When the dataframe from which the plot is constructed is loaded, downsampling of the data is automatically triggered if there are more than 15 million rows to make interacting with the notebook more nimble. You may adjust this threshold to suit your computational resources.

The downsampling function preserves points in less dense areas of the plot by constructing a 2D histogram of the x and y coordinates, and only removing data from regions containing more points than the specified threshold value. Note that it does this without taking into account other features of the data, such as k-mer coverage. Reads that were dropped as a result of downsampling will, not be displayed in the downloaded table. In some instances, it may therefore be useful to pre-filter data based on other criteria.

The interface may take a while to load for large datasets, so keep an eye on the terminal if nothing appears to be happening in the browser.

### Blast query turnaround
When using the remote_blast option, it may take a very long time to get a result when the server is busy. For production use, sending queries to a local server may therefore be preferable. Alternatively, for preliminary exploration, consider sampling reads near local [histogram peaks](https://github.com/CobiontID/read_VAE/tree/main/read_tools/plotting_tools) and running a local query.

