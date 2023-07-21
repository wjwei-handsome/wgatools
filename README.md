## *W*hole *G*enome *A*lignment **T**ools

### !!! Early Stage !!!

### WHAT HAVE DONE
- [x] PAF format reader
- [x] MAF format reader
- [x] CIGAR string parser
- [x] MAF2PAF
- [x] MAF2Chain
- [x] MAF2BAM
- [x] PAF2Chain
- [x] PAF2Blocks
- [x] PAF2MAF
- [x] Call Variants from MAF

### WHAT WILL DO
- [ ] convert MAF to [blocks, delta]
- [ ] convert PAF to [bam, chain]
- [ ] convert BAM to [blocks, delta, MAF, PAF, chain]
- [ ] convert chains to [delta, MAF, PAF, blocks]
- [ ] convert delta to [chain, MAF, PAF, blocks]
- >NOTE: some conversions is unnecessary, but it's a good practice to implement them.
- [ ] visualize genome alignment
- [ ] call variants and statistics/visualize them

### A sample example of conversions

```rust
fn main() {
    let your_path = "/your/path/to/your.paf";
    let mut reader = PafReader::from_path(your_path).unwrap();
    reader.convert("out.chain", FileFormat::Chain);
    let mut reader = PafReader::from_path(your_path).unwrap();
    reader.convert("out.blocks", FileFormat::Blocks);
    let your_maf_path = "/your/path/to/your.maf";
    let mut reader = MafReader::from_path(your_path).unwrap();
    reader.convert("test.paf", FileFormat::Paf);
}
```
> It should be extremely fast!![img](https://raw.githubusercontent.com/wjwei-handsome/wwjPic/main/img/20230706022535.png)

### Features

- use `nom` to parse CIGAR string