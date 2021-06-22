# Utilities for VAE plots
## Labelling
- **category_labels.py**: Generate a list of integer labels based on one of more lists of sequence identifiers. Can be used to colour read plots according to taxonomic classification from e.g. Kraken.
- **category_labels_from_cont.py**: Generate bins based on n quantiles from a continuous vector
- **category_labels_from_cont_manual.py**: Generate bins based on a predefined list of edges. Useful for e.g. slicing FastK median vectors according to the minima in a GenomeScope profile.
