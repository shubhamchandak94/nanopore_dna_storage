import json
original_file = '/raid/nanopore/shubham/Nanopore-DNA-storeage/oligos_2_7_19/code/conv_m11_nosync_init_64payload_init10010110001_oligos_raptor_coded'
recon_file = 'tmpfile.m11_nosync.l_64.280743'
with open(recon_file) as f:
    m1 = json.loads(f.read())
with open(original_file) as f:
    m2 = json.loads(f.read())
print('Number of symbols in original:',len(m2['symbols']))
print('Number of symbols in reconstruction:',len(m1['symbols']))
d1 = {l[0]:l[1] for l in m1['symbols']}
d2 = {l[0]:l[1] for l in m2['symbols']}

to_be_removed = []
for k in d1:
    if k not in d2:
        # out of range
        to_be_removed.append(k)

for k in to_be_removed:
    del d1[k]
print('Number of symbols out of range:',len(to_be_removed))
print('Number of erroneous symbols in reconstruction:',sum([d1[k]!=d2[k] for k in d1]))
