#[macro_use] extern crate log;

use clap::Parser;
use rust_htslib::bcf::{Header, HeaderRecord, Reader, Read, Writer, Format};
use rust_htslib::bcf::header;
use std::collections::{BTreeMap, BTreeSet};

/// Simple program to copy INFO fields to FORMAT tags
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
   /// Path to the input VCF or "-" to read from STDIN
   #[arg(short, long)]
   input: String,

   /// Path to the output VCF or "-" to write to STDOUT
   #[arg(short, long)]
   output: String,

   /// The INFO fields to copy over to FORMAT tag
   #[arg(short, long)]
   fields: Vec<String>,

   /// Transfer also the QUAL tag into FORMAT
   #[arg(short, long)]
   qual: bool,

   /// Show verbose output (sets log-level to debug or trace)
   #[arg(short, long, action = clap::ArgAction::Count)]
    verbose: u8,

   /// Indicates how often the tool should report the progress if verbose output is enabled
   #[arg(long = "verbose-report", default_value_t = 10000)]
    verbose_report: u32,
}

enum TagValue {
	Flag(i32),
	Integer(Vec<i32>),
	Float(Vec<f32>),
	String(Vec<Vec<u8>>)
}

fn main() {
  env_logger::init();

  trace!("Parsing arguments");
  let args = Args::parse();

  if args.verbose >= 3 {
    log::set_max_level(log::LevelFilter::Trace)
  } else if args.verbose >= 2 {
    log::set_max_level(log::LevelFilter::Debug)
  } else if args.verbose >= 1 {
    log::set_max_level(log::LevelFilter::Info)
  }

  if args.fields.len() == 0 && !args.qual {
		error!("No field for conversion identified. Use '-q' or '-f' options");
		return;
	}

  let mut fields : BTreeSet<String> = args.fields.iter()
		.map(|x| x.to_owned())
		.collect();

	// Open the BAM File and extract information from the header
	let mut input = if args.input == "-" {
    debug!("Reading from STDIN");
    Reader::from_stdin().expect("Can not open input stream")
  } else {
    debug!("Opening input VCF at {}", args.input);
    Reader::from_path(&args.input).expect("Can not open input VCF file")
  };

  trace!("Extracting header information");
	let header = input.header();
	if header.sample_count() != 1 {
		error!("input is not a single-sample VCF");
		return;
	}

  trace!("Building new header");
	let mut new_header = Header::from_template(header);
	// Store a map of field-names to data types
	let mut field_types = BTreeMap::new();
	for rec in header.header_records() {
	  if let HeaderRecord::Info { key: _, values: v } = rec {
			if let Some(id) = v.get("ID") {
				if fields.contains(id) {
					fields.remove(id);
					new_header.remove_info(id.as_bytes());

					let new_record = format!("##FORMAT=<ID={},Number={},Type={},Description={}>",
						v.get("ID").unwrap(), v.get("Number").unwrap(), v.get("Type").unwrap(), v.get("Description").unwrap()
					);
					trace!("Adding new FORMAT header record: {}", new_record);
					new_header.push_record(new_record.as_bytes());

					match v.get("Type").unwrap().as_str() {
						"Flag" => field_types.insert(id.to_owned(), header::TagType::Flag),
						"Integer" => field_types.insert(id.to_owned(), header::TagType::Integer),
						"Float" => field_types.insert(id.to_owned(), header::TagType::Float),
						"String" => field_types.insert(id.to_owned(), header::TagType::String),
						_ => {
							error!("Unknown tag type {}", v.get("Type").unwrap());
							return();
						}
					};

				}
			}
		}
	}
	trace!("Found {} INFO fields: {:?}", field_types.len(), field_types);
	for f in args.fields {
		if !field_types.contains_key(&f) {
			error!("Error: input VCF does not contain INFO tag '{}'", f);
			return;
		}
	}

	if args.qual {
	  let new_record = "##FORMAT=<ID=QUAL,Number=1,Type=Float,Description=Phred-scaled quality score for the assertion made in ALT>";
		trace!("Adding new FORMAT header record: {}", new_record);
		new_header.push_record(new_record.as_bytes());
	}


  let mut output = if args.output == "-" {
    debug!("Writing to STDOUT");
    Writer::from_stdout(&new_header, true, Format::Vcf).unwrap()
  } else {
    debug!("Opening output VCF at {}", args.output);
    Writer::from_path(&args.output, &new_header, true, Format::Vcf).unwrap()
  };

  // Cound the number of records that were processed
	let mut n_records = 0;
	for r in input.records() {
		// Extract the record
		let mut rec = match r {
		  Ok(i) => i,
		  Err(e) => {
		    error!("Malformed VCF record: {}", e);
		    return;
		  }
		};

    // Check verbose reporting
		n_records += 1;
	  if n_records % args.verbose_report == 0 {
      trace!(" - {n_records} processed")
	  }

    // Record the data for each field that is available
		let mut data = BTreeMap::new();
    for (tag, ttype) in &field_types {
			let tagb = tag.as_bytes();
			match ttype {
				header::TagType::Flag => {
					let v = rec.info(tagb).flag().unwrap();
					if v {
						data.insert(tag, TagValue::Flag(1));
						rec.clear_info_flag(tagb).expect("Can not remove INFO tag");
					} else {
						data.insert(tag, TagValue::Flag(0));
					}
				},
				header::TagType::Integer => {
					if let Some(v) = rec.info(tagb).integer().unwrap() {
						let x : Vec<i32> = v.iter().map(|x| *x).collect();
						data.insert(tag, TagValue::Integer(x));
						rec.clear_info_integer(tagb).expect("Can not remove INFO tag");
					}
				},
				header::TagType::Float => {
					if let Some(v) = rec.info(tagb).float().unwrap() {
						let x : Vec<f32> = v.iter().map(|x| *x).collect();
						data.insert(tag, TagValue::Float(x));
						rec.clear_info_float(tagb).expect("Can not remove INFO tag");
					}
				},
				header::TagType::String => {
					if let Some(v) = rec.info(tagb).string().unwrap() {
						let x : Vec<Vec<u8>> = v.iter().map(|x| x.to_vec()).collect();
						data.insert(tag, TagValue::String(x));
						rec.clear_info_string(tagb).expect("Can not remove INFO tag");
					}
				}
			}
		}

    // Replace the header for the record
		output.translate(&mut rec);

    // Re-Insert the data as Format tag
		for (tag, data) in &data {
			let tagb = tag.as_bytes();
			match data {
				TagValue::Flag(v) => {
					rec.push_format_integer(tagb, &[*v]).expect("Can not store FLAG-Integer in the format");
				},
				TagValue::Integer(v) => {
					rec.push_format_integer(tagb, &v).expect("Can not store Integer in the format");
				},
				TagValue::Float(v) => {
					rec.push_format_float(tagb, &v).expect("Can not store Float in the format");
				},
				TagValue::String(v) => {
					rec.push_format_string(tagb, &v).expect("Can not store String in the format");
				}
			}
		}
		if args.qual {
		  let q = [rec.qual()];
      rec.push_format_float("QUAL".as_bytes(), &q).expect("Can not store Float in the format");
		}


    // Write the record to the output stream
		output.write(&rec).expect("Can not write record to output stream");
	}
	info!("Finished transfering fields for {n_records} records")
}
