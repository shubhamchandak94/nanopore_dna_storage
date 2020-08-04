18 experiments, first 9 with single CRC, next 9 with 2 CRC idea (index+CRC1+data+CRC2, each CRC 8 bit, CRC1 covering index and first half of data (in blocks of 16, first half includes lesser blocks when num RS blocks is odd), CRC2 covering index and second half of data). The first CRC is placed near index to protect it when second half is fully corrupted. The first half is kept smaller to make sure we recover it in more cases when index is preserved. We perform shuffling of oligo symbols to make sure the bias towards getting first half correct does not adversely affect RS performance. Experiments: m = 6, 8, 11; rate = 3/4, 5/6, 7/8. Number of RS blocks per oligo chosen so that oligo length <= 120+50 for each oligo. Barcodes were chosen using minimum edit distance, no repeating bases, and GC balance criteria (see generate_barcodes.py).

Other differences from previous experiment:
- focus on computationally practical codes: so try m=6 instead of m=14
- try very high rate code: rate = 7/8, remove very low rate code: rate = 1/2
- reduce RS slightly to get lower writing costs (0.25 as compared to 0.3)
- added some poems since we had more space due to increased pool size and lower RS
