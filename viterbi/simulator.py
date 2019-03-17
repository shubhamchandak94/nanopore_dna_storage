import sys
import numpy as np
from uuid import uuid4
from fast5_research import Fast5
import scrappy
import editdistance
import scipy.stats as st
import mappy as mp
import util

SYN_SUB_PROB = 0.0
SYN_DEL_PROB = 0.0
SYN_INS_PROB = 0.0

# from deepsimulator (probably more realistic dwell times)
def rep_rvs(size,a):
    # array_1 = np.ones(int(size*0.075)).astype(int)
    # samples = st.alpha.rvs(3.3928495261646932, 
    #     -7.6451557771999035, 50.873948369526737, 
    #     size=(size-int(size*0.075))).astype(int)
    # samples = np.concatenate((samples, array_1), 0)
    # samples[samples<1] = 2
    # print(a)
    a = a*5
    array_1 = np.ones(int(size*(0.075-0.015*a))).astype(int)
    samples = st.alpha.rvs(3.3928495261646932+a, 
        -7.6451557771999035+(2*a), 50.873948369526737, 
        size=(size-int(size*(0.075-0.015*a)))).astype(int)
    samples = np.concatenate((samples, array_1), 0)
    samples[samples<1] = 2
    np.random.shuffle(samples)
    return samples


# can handle fasta header, but also works without it
infile_seq = sys.argv[1]
f = open(infile_seq)
seq = f.readline().rstrip('\n')
if seq[0] == '>':
	seq = f.readline().rstrip('\n')
f.close()
len_seq = len(seq)
print('Length of seq: ', len_seq)

syn_seq = util.simulate_indelsubs(seq, sub_prob = SYN_SUB_PROB, del_prob = SYN_DEL_PROB, ins_prob = SYN_INS_PROB)

squiggle_array = scrappy.sequence_to_squiggle(syn_seq,rescale=True).data(as_numpy=True)
raw_data = np.array([])

# for dwell time use deepsimulator since the one provided by scrappie is way too good
# scrappie dwell gives around 3-4% edit distance while ds dwell gives around 15%
#ds_alpha = 0.1 # 0.1 is default parameter in deepsim
#squiggle_array[:,0] = rep_rvs(squiggle_array.shape[0], ds_alpha)

for squig in squiggle_array:
	mean = squig[1]
	stdv = squig[2]
	dwell = squig[0]
	raw_data = np.append(raw_data, np.random.laplace(mean, stdv/np.sqrt(2), int(round(dwell))))

print('Length of raw signal: ', len(raw_data))

file_raw = 'raw_data.txt'
f = open(file_raw,'w')
for val in raw_data:
	f.write(str(val)+'\n')
f.close()

# create fast5 (from https://nanoporetech.github.io/fast5_research/examples.html)
filename='my_new.fast5'

# example of how to digitize data
start, stop = int(min(raw_data - 1)), int(max(raw_data + 1))
rng = stop - start
digitisation = 8192.0
bins = np.arange(start, stop, rng / digitisation)
# np.int16 is required, the library will refuse to write anything other
raw_data_binned = np.digitize(raw_data, bins).astype(np.int16)

# The following are required meta data
channel_id = {
    'digitisation': digitisation,
    'offset': 0,
    'range': rng,
    'sampling_rate': 4000,
    'channel_number': 1,
    }
read_id = {
    'start_time': 0,
    'duration': len(raw_data),
    'read_number': 1,
    'start_mux': 1,
    'read_id': str(uuid4()),
    'scaling_used': 1,
    'median_before': 0,
}
tracking_id = {
    'exp_start_time': '1970-01-01T00:00:00Z',
    'run_id': str(uuid4()).replace('-',''),
    'flow_cell_id': 'FAH00000',
}
context_tags = {}

with Fast5.New(filename, 'w', tracking_id=tracking_id, context_tags=context_tags, channel_id=channel_id) as h:
    h.set_raw(raw_data_binned, meta=read_id, read_number=1)

#basecall = scrappy.basecall_raw(raw_data, model='rnnrf_r94')
#basecall_seq = basecall[0]

basecall_seq = scrappy.basecall_raw_python(raw_data)
#basecall_seq = scrappy.basecall_raw_python_no_homopolymer(raw_data)
#basecall_seq = scrappy.basecall_raw_python_sync(raw_data,period=5)
#basecall_seq = scrappy.basecall_raw_python_vocab(raw_data,['AAAAA','CCCCC','GGGGG','TTTTT'])
#basecall_seq = scrappy.basecall_raw_python_vocab(raw_data,['AGCT','ATCG'])
#basecall_seq = scrappy.basecall_raw_python_vocab(raw_data,['AGCT','ATCG','TGAC','CATG','ACGT','ATGC'])
#basecall_seq = scrappy.basecall_raw_python_vocab(raw_data,['AGGTTC','CAGACT','ATGACC','TGACGA','CTGCAT','TGTACG','ATTCGC','GATCTA'])
#basecall_seq = scrappy.basecall_raw_python_vocab(raw_data,['CGTAACCTGT','GGTCGCCGCC'])
#basecall_seq = scrappy.basecall_raw_python_vocab(raw_data,['CAGTTC','TATCGA'])
#basecall_seq = scrappy.basecall_raw_viterbi_conv(raw_data)
vocab_random_150 = [
'GGCCTACCGA',
'GGACCGCTTA',
'GTAGGACTGG',
'TGGCTGGGTT',
'AAGCCTTAAT',
'CAAGCGAAAA',
'TCCCCGGGAT',
'CTCCGAATAG',
'CTGCATTCAT',
'AACGATTTTC',
'ACGAAAGCAA',
'TGTTCTGGAC',
'ATAGGCCTAT',
'CATGAATTAG',
'ATGCTTGGAC',
]
#basecall_seq = scrappy.basecall_raw_python_vocab(raw_data,vocab_random_150)

vocab_random_150_30 = [
'GGCCTACCGAGGACCGCTTAGTAGGACTGG',
'TGGCTGGGTTAAGCCTTAATCAAGCGAAAA',
'TCCCCGGGATCTCCGAATAGCTGCATTCAT',
'AACGATTTTCACGAAAGCAATGTTCTGGAC',
'ATAGGCCTATCATGAATTAGATGCTTGGAC',
]
#basecall_seq = scrappy.basecall_raw_python_vocab(raw_data,vocab_random_150_30)

# find edit distance
edit_distance = editdistance.eval(basecall_seq,seq)
print('Edit distance: ', edit_distance)
print('Ratio: ', edit_distance/len_seq)

# generate fastq file
filename='my_new.fastq'
f = open(filename, 'w')
f.write('@my_new\n')
f.write(basecall_seq+'\n')
f.close()

# align to original and write to sam
filename='my_new.sam'
f = open(filename,'w')
a = mp.Aligner(seq = seq, preset='map-ont')
for hit in a.map(basecall_seq):
	f.write("{}\t{}\t{}\t{}\tNM:{}\n".format(hit.ctg, hit.r_st, hit.r_en, hit.cigar_str, hit.NM))
f.close()
