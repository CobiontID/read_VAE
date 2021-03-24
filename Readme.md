# Count the number of unique k-mers in a set of sequences

Counts distinct k-mers for each sequence record in a file and divides it by the total number of k-mers in that record, where 4^k >= record length.

## Usage

Input: Multi-entry fasta file

`unique-kmers --file <fasta file> --klength <k-mer length> --out <output file>`

Note that k-mers will not be canonicalised.
  
## Dependencies
Rust (see https://www.rust-lang.org/tools/install)

needletail = "0.4"
fnv = "1.0.7"
csv
clap
