'''
Usage:
    python3 compute_coverage_stats.py aligned_sequence_list_file total_num_oligos coverage
    aligned_sequence_list_file: file with each line having the sequence the read aligns to
'''

import sys
import random
from collections import Counter
import numpy as np
import scipy.stats

random.seed(999)

if len(sys.argv) != 4:
    print('Incorrect number of arguments')
    sys.exit(1)

with open(str(sys.argv[1])) as f:
    original_ref_list = [l.rstrip('\n') for l in f.readlines()]
num_oligos=int(sys.argv[2])
coverage=float(sys.argv[3])
num_reads = int(num_oligos*coverage)
print('num_oligos',num_oligos)
print('coverage',coverage)
print('num_reads',num_reads)

subsampled_ref_list = random.choices(original_ref_list,k=num_reads)
c = Counter(subsampled_ref_list)

count_list = []
for elem in c:
    count_list.append(c[elem])

num_seen = len(count_list)

count_list += (num_oligos-len(count_list))*[0]
print('mean',np.mean(count_list))
print('variance',np.var(count_list))
print('Poisson sampling variance',np.mean(count_list))
print('normalized variance',np.var(count_list)/np.mean(count_list))
print('number of oligos with 0 reads',num_oligos-num_seen)
print('fraction of oligos with 0 reads',1-num_seen/num_oligos)

bincount = np.bincount(count_list)
for b in bincount:
    print(b)
