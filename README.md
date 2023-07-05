## *W*hole *G*enome *A*lignment **T**ools

### !!! Early Stage !!!

### A sample example of conversion from PAF to Chain or Blocks

```rust
fn main() {
    let your_path = "/your/path/to/your.paf";
    let mut reader = PafReader::from_path(your_path).unwrap();
    reader.convert("out.chain", FileFormat::Chain);
    let mut reader = PafReader::from_path(your_path).unwrap();
    reader.convert("out.blocks", FileFormat::Blocks);
}
```
> It should be extremely fast!![img](https://raw.githubusercontent.com/wjwei-handsome/wwjPic/main/img/20230706022535.png)

### Features

- use `nom` to parse CIGAR string