extern crate needletail;
use needletail::{parse_fastx_file, Sequence};
use fnv::FnvHashMap;
use clap::{Arg, App};
extern crate csv;

use std::error::Error;
use csv::Writer;

fn main()  -> Result<(), Box<dyn Error>> {

    let matches = App::new("Umique k-mer counter")
        .version("0.1.0")
        .author("Claudia C. Weber <cw21@sanger.ac.uk>>")
        .about("Tally nucleotide counts in multi-entry fasta")
        .arg(Arg::with_name("file")
                 .short("f")
                 .long("file")
                 .takes_value(true)
                 .required(true)
                 .help("Fasta file to tally."))
        .arg(Arg::with_name("out")
                 .short("o")
                 .long("out")
                 .takes_value(true)
                 .required(true)
                 .default_value("counter_output.txt")
                 .help("Output file name."))
        .arg(Arg::with_name("klength")
                 .short("k")
                 .long("klength")
                 .takes_value(true)
                 .default_value("4")
                 .help("K-mer length"))
        .get_matches();
    
        let filename = matches.value_of("file").unwrap();
        let out = matches.value_of("out").unwrap();
        let k = matches.value_of("klength").unwrap();
        println!("K-mer length: {:#?}", k);
        println!("File: {:#?}", filename);
        println!("Output file: {:#?}", out);


    let mut reader = parse_fastx_file(&filename).expect("valid path/file");
    let mut wtr = Writer::from_path(out)?;

    while let Some(record) = reader.next() {
        let seqrec = record.expect("invalid record");
        let mut n_kmers = 0;
        let mut k_counts = FnvHashMap::default();
        let norm_seq = seqrec.normalize(false);
        //let rc = norm_seq.reverse_complement();
        for kmer in norm_seq.kmers(k.parse::<u8>().unwrap()) {
            n_kmers += 1;
            *k_counts.entry(kmer).or_insert(0) = 0;
        }

        let mut n_kmers_distinct = 0;
        for (_key, _value) in k_counts {
            n_kmers_distinct +=1;
        }

        wtr.write_record(&[(f64::from(n_kmers_distinct) / f64::from(n_kmers)).to_string()])?;

    }

    wtr.flush()?;
    Ok(())
}

