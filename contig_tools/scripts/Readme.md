## Draw contigs/scaffolds
**Select_contigs_reduced_multi.py**
Generate HTML file of decomposed tetranucleotide plots with binned annotations.


## Hi-C
**hic_links.py** Generate scaffold connectivity matrix from SALSA or YaHs pair file (see [here](https://github.com/c-zhou/yahs/blob/0f81dbf678abd86b601374f960bcb4e0b8d33426/scripts/run_yahs.sh#L48)) and store the resulting np array in a .npy file. Optionally returns a 1/0 vector for connected/uncnnected scaffolds, or (number of connections)/(scaffold length), for use with `Select_contigs_reduced_multi.py`.
