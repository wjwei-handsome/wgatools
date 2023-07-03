use std::fs::File;
use wgatools::errors::ParseError;

fn main() -> Result<(), ParseError> {
    match File::open("ttt") {
        Ok(_file) => {
            println!("nothing");
            Ok(())
        }
        Err(e) => Err(e.into()),
    }
}
