# Read VAE dashboard

To explore the structure of a given readset, interactively filtering the data and cross-referencing different sequence features and taxonomic labels can be useful. A demo is available [here](https://huggingface.co/spaces/cc7740/read_VAE).

![image](https://github.com/CobiontID/read_VAE/assets/10507101/42d5a687-4af4-4900-9843-e1b597856c17)

See https://cobiontid.github.io/reads.html for more details and a walkthrough of the example.

## Running the dashboard 
The server can be started with:

```panel serve --port=$port --address="0.0.0.0" --allow-websocket-origin=localhost:$port --show read_dash.py --args --config example.yml```

Here,`$port` is the number of an open port on the local machine. `example.yml` is the path to the sample configuration file. See below for more information on configuration. For especially large datasets, it may be necessary to set `--session-token-expiration` to a larger value. If the configuration file is updated after the app has been launched, refreshing the browser window will prompt all of the data to be reloaded using the updated file paths or class lists, if supplied.

If the dashboard is run on a remote machine, connecting with the SSH extension in Visual Studio Code and entering the above in the terminal may be convenient, as port forwarding will be handled automatically. When prompted "Do you want Code to open the external website?", clicking "open" will open the dashboard in a browser window.

### Setting up the configuration file
#### Basic configuration
- `tolid`: Sample identifier (this should match the identifier used to run the pipeline, unless the paths to the required files are manually set)
- `fasta`: Path to the fasta file containing the sequences of the reads. Required to launch blast queries.
- `basedir`: Path to the base directory containing outputs from the read VAE workflow (but see details on how to override defaults).
- `class_lists`: One or more files containing the identifiers of reads belong to different classes. Each file should be on a separate line, preceded by a dash. If left blank, class labels will not be used. 

#### Options to manually configure paths
- `k`: The k-mer size used with FastK in the pipeline (used to retrieve the median k-mer count file). If not specified, k defaults to 31.
- `vae_path`: A custom path pointing to the 2D coordinates to be used for plotting. Input files may either be in text format, e.g. ilApaMono1_16.vae.out.2d.0, or in .npy format (the latter is preferable for large datasets, since it allows faster loading).
- `read_ids_path`: Custom path to read identifiers
- `fastk_path`: Custome path pointing to the file containing median read k-mer counts
  
#### To use indexed fasta files
For large readsets, it may be useful to index the fasta file containing the reads with `samtools faidx`. To make use of an index, include these lines:
```
indexed: True
samtools_path: /software/samtools
```
The path of the .fai file is currently assumed to be the same as for the corresponding fasta file, plus the suffix.
