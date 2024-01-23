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
- [x] Visualize MAF file in terminal
- [x] Extract regions from MAF file
- [x] Build MAF index
- [x] Statistics of MAF/PAF file

### WHAT WILL DO IN FUTURE

- [ ] SAM converter [really need?]
- [ ] Local improvement of alignment by re-alignment
- [ ] MAF -> GAF -> HAL
- [ ] for BIG MAF, should optimize
- [ ] split & chop MAF file

### Install

```shell
git clone https://github.com/wjwei-handsome/wgatools.git
cd wgatools
cargo build --release
```

or just install from git:

```shell
cargo install --git https://github.com/wjwei-handsome/wgatools.git
```

#### Nix

A nix flake is also available. You can build from within the repo like this:

```shell
nix build .#wgatools
```

Or directly install from github:

```shell
nix profile install github:wjwei-handsome/wgatools
```

#### Docker and Singularity

Using nix, we can derive docker and singularity images:

```shell
nix build .#dockerImage
```

First, we load the docker image into the local daemon:

```shell
docker load < result
```

It's then possible to pack up a singularity image:

```shell
singularity build wgatools-$(git log -1 --format=%h --abbrev=8).sif docker-daemon://wgatools:latest
```

This can be useful when running on HPCs where it might be difficult to build wgatools.

### Usage

```shell
> wgatools
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
  maf-ext    Extract specific region from MAF file with index [aliases: me]
  call       Call Variants from MAF file [aliases: c]
  maf2sam    TEST: maf2sam [aliases: m2s]
  maf-index  Build index for MAF file [aliases: mi]
  tview      View MAF file in terminal [aliases: tv]
  stat       Statistics for Alignment file [aliases: st]
  help       Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help (see more with '--help')
  -V, --version  Print version

GLOBAL:
  -o, --outfile <OUTFILE>  Output file ("-" for stdout) [default: -]
  -r, --rewrite            Bool, if rewrite output file [default: false]
  -t, --threads <THREADS>  Threads, default 1 [default: 1]
  -v, --verbose...         Logging level [-v: Info, -vv: Debug, -vvv: Trace, defalut: Warn]
```

> NOTE: If you want to convert into MAF format, you should provide target and query genome sequence files in [.fa/.fa.gz].

### Examples

#### visualize MAF file in terminal

![example](./example.gif)

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

- [x] Error detection and handling
- [ ] Test cases
- [ ] Documentations

[//]: # "> It should be extremely fast!![img](https://raw.githubusercontent.com/wjwei-handsome/wwjPic/main/img/20230706022535.png)"

## Features

- use `nom` to parse CIGAR string
- use `rayon` to accelerate the speed of conversions
- use `ratatui` to visualize MAF file in terminal
- ...

## Contributing

Feel free to dive in! [Open an issue](https://github.com/wjwei-handsome/wgatools/issues/new) or submit PRs.

## License

[MIT License](./LICENSE) © WenjieWei
