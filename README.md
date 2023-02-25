# vcf-info2format

A small utility to move INFO fields in a single-sample VCF into FORMAT fields

## Example Usage

If you have pre-existing single-sample VCF files or from a per-sample variant calling workflow, you
may want to copy some INFO fields to FORMAT before merging to retain the information for later use. 
In such cases, you can use this utility for move info-fields to format-fields and retain them like this:

```sh
bcftools call sampleA.bam \
	| vcf-info2format --fields FS --qual \
	> sampleA.vcf
bcftools call sampleB.bam \
	| vcf-info2format --fields FS --qual \
	> sampleB.vcf

bcftools merge sampleA.vcf sampleB.vcf > multi-sample.vcf
```

## Options 

```
Usage: vcf-info2format [OPTIONS] --input <INPUT> --output <OUTPUT>

Options:
  -i, --input <INPUT>
          Path to the input VCF or "-" to read from STDIN
  -o, --output <OUTPUT>
          Path to the output VCF or "-" to write to STDOUT
  -f, --fields <FIELDS>
          The INFO fields to copy over to FORMAT tag
  -q, --qual
          Transfer also the QUAL tag into FORMAT
  -v, --verbose...
          Show verbose output (sets log-level to debug or trace)
      --verbose-report <VERBOSE_REPORT>
          Indicates how often the tool should report the progress if verbose output is enabled [default: 10000]
  -h, --help
          Print help
  -V, --version
          Print version
```


## Installation

Install the Rust toolchain, e.g. using [rustup](https://rustup.rs), [conda](https://anaconda.org/conda-forge/rust), or [docker](https://hub.docker.com/_/rust), in version 1.67.0 or higher. 
Then run: 

```sh
cargo install --git https://github.com/Evotec-Bioinformatics/vcf-info2format.git
```

## License

[GNU General Public License v3.0](LICENSE)

