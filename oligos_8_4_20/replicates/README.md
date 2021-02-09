Command for demultiplexing (based on [this webpage](https://community.nanoporetech.com/protocols/Guppy-protocol/v/gpb_2003_v1_revt_14dec2018/barcoding-demultiplexing)):
```
./ont-guppy/bin/guppy_barcoder --input_path guppy_out --save_path guppy_out_demultiplexed --config configuration.cfg --barcode_kits "EXP-NBD104 EXP-NBD114" -x cuda:2 -t 20
```
