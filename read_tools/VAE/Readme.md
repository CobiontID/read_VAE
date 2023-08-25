# VAE

To embed a dataset with the default configuration, run:

```python vae.py --countfile sampleid.tetranucs.npy --outdir ./VAE/sampleid/ --fignames sampleid```

Only the first argument is required. The necessary packages are included in the yml at the base of this repository.

## Configuration
### Basic input and output options
- `--countfile`: Path to 2D numpy array with k-mer counts to embed (required)
- `--outdir`: Output directory (default is the current working directory)
- `--latent_format`: Format to store latent representations in (defaults to txt for workflow compatibility, but npy is preferable for large datasets due to faster I/O)
- `--plot_latents`: Plot the latent space after training is complete (enabled by default, provided the number of latent dimensions equals 2)
- `--fignames`: Specify names for files with embeddings and figures (if no name is specified, "sample" will be used as a default)

### Training settings
- `--latent_dim` Number of latent dimensions (default = 2, other values will disable plotting functions)
- `--kl`, Weight for KL loss, corresponds to beta in a beta-VAE (beta == 1 corresponds to a standard VAE) (default = 0.0025)
- `--batch_size` Batch size for training (default = 256)
- `--epochs` Number of epochs to train for, provided early stopping is not triggered (default = 10)
- `--learning_rate` Set the learning rate (default = 0.001)
- `--reduce_beta` Reduce weight on KL loss (beta) if variance indicates posterior collapse at the end of each epoch (enabled by default). In some cases, this may not be sufficient to prevent latent collapse and manual tuning will be necessary.
- `--cb` Apply the CB correction (enabled by default). Omitting the correction may decrease runtime.

### Additional options
- `--model_out`: Directory path prefix to save model weights (optional)
- `--intermediate`: Draw embeddings for intermediate steps at the end of each epoch to track training (off by default, requires the number of latent dimensions to equal 2). Files will be labelled with "start" for the initial weights or the number of the epoch.
- `--extra_columns` Additional input feature data to load into matrix.
- `--remove` File specifying reads to remove before embedding (0: keep; other integer value: drop; one integer per line). May be useful for visualising low coverage reads, for example.
- `--verbose` Print extra model info, such as a summary of the encoder and decoder, which includes the shapes of the layers and the number of parameters (off by default).