![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/wjwei-handsome/wgatools/ci.yml)
![GitHub repo size](https://img.shields.io/github/repo-size/wjwei-handsome/wgatools)

## *W*hole *G*enome *A*lignment **T**ools

## A Rust library and tools for whole genome alignment files

## TOOLS

### WHAT HAVE DONE
- [x] PAF file reader
- [x] MAF file reader
- [x] Chain file reader
- [x] CIGAR string parser
- [x] MAF2PAF
- [x] MAF2Chain
- [x] PAF2Chain
- [x] PAF2Blocks
- [x] PAF2MAF
- [x] Chain2MAF
- [x] Chain2PAF
- [x] Call Variants from MAF

### WHAT WILL DO IN FUTURE
- [ ] SAM file reader [really need?]
- [ ] SAM converter [really need?]
- [ ] Visualize genome alignment
- [ ] Call variants and statistics/visualize them
- [ ] Local improvement of alignment by re-alignment

### Install

```shell
git clone https://github.com/wjwei-handsome/wgatools.git
cd wgatools
cargo build --release
```

### Usages

```shell
> ./target/release/wgatools
wgatools -- a cross-platform and ultrafast toolkit for Whole Genome Alignment Files manipulation

Version: 0.1.0

Authors: Wenjie Wei <wjwei9908@gmail.com>

Usage: wgatools [OPTIONS] <COMMAND>

Commands:
  maf2paf    Convert MAF format to PAF format [aliases: m2p]
  maf2chain  Convert MAF format to Chain format [aliases: m2c]
  paf2maf    Convert PAF format to MAF format [aliases: p2m]
  paf2chain  Convert PAF format to Chain format [aliases: p2c]
  chain2maf  Convert Chain format to MAF format [aliases: c2m]
  chain2paf  Convert Chain format to PAF format [aliases: c2p]
  help       Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help (see more with '--help')
  -V, --version  Print version

GLOBAL:
  -o, --outfile <OUTFILE>  Output file ("-" for stdout) [default: -]
  -r, --rewrite            Bool, if rewrite output file [default: false]
  -t, --threads <THREADS>  [default: 1]
```

> NOTE: If you want to convert into MAF format, you should provide target and query genome sequence files in [.fa/.fa.gz].

## Library

Some simple reader and iterator for PAF, MAF and Chain files: 
```rust
use wgatools::parser::paf::PafReader;
use wgatools::parser::maf::MAFReader;
use wgatools::parser::chain::ChainReader;
fn main() {
    let mut mafreader = MAFReader::from_path("test.maf").unwrap();
    for record in mafreader.records() {
        let record = record.unwrap();
        println!("{:?}", record);
    }
    /// ...
}
```

### TODO for library
- [ ] Error detection and handling
- [ ] Test cases
- [ ] Documentations

[//]: # (> It should be extremely fast!![img]&#40;https://raw.githubusercontent.com/wjwei-handsome/wwjPic/main/img/20230706022535.png&#41;)

## Features

- use `nom` to parse CIGAR string
- use `rayon` to accelerate the speed of conversions
- ...

## Contributing

Feel free to dive in! [Open an issue](https://github.com/wjwei-handsome/GeneMap/issues/new) or submit PRs.

## License

GPL-3.0 © Weiwenjie