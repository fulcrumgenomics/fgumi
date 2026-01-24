#### Summary

Using this crate is much faster than flate2 due to the use of ISA-L under the hood. 
However flate2 is _much_ more widely used. I suggest giving this a shot and provide feedback.

Criterion benchmark report available here: https://milesgranger.github.io/isal-rs/benches/criterion/report/

#### To run yourself:

- `cd benches`
- `tar -xvzf data.tar.gz`
- `cargo bench`
