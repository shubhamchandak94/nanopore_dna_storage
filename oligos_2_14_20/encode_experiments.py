import sys
sys.path.insert(1, '../')
import helper

barcode_start = [
        'GCTACATGTATACTGCGAGACAGAC',
        'TCTATCTACTCGTGCTCGCTAGCTG',
        'TGAGATCACAGCTACATAGTGAGAG',
        'AGCGTACACGACTGAGCACACTACG',
        'CATCAGCAGTAGAGAGTAGCGCGAT',
        'GAGTCTCTAGCGCTACGAGATATAT',
        'CACGAGATCTCAGTGTCGACACGTG',
        'CGCTGCAGTCTATCTCTCTGTACAT',
        'GTGACTCTGCATATGCTGTCTCGAT',
        'ATGAGAGAGAGCGCTCTCTCATGAG',
        'GTCGTACTGTGACTCTCACTCGTGC',
        'CGCACGCATATCACGATGAGTAGCT',
        'CGTATGCTCGTATGAGCATAGAGTG',
        'GTGTGAGCTATGCGAGCGACGATCT',
        'TATGTACTAGCTAGTCACGCACACA',
        'CTAGACGTGCGAGTATACTACTATG',
        'TGATCTATGATCGATGTACAGCGCG',
        'TCTCGCTCGAGCACAGAGATAGCGA',
        ]

barcode_end = [
        'CGATAGTCGCAGTCGCACATCACTC',
        'TGTCTGCACTGCACTAGTCGCATGT',
        'GATAGAGCACTGTCGTGTGCAGTCA',
        'ATGCGCTCACACGCATCTGCTAGCA',
        'TATCATCGACGCTAGCAGTGTCTGC',
        'AGCTCTGATGAGATCAGCAGACTGT',
        'AGTGATAGTGATCTCGACGCAGCTA',
        'GCTAGTAGAGATACATCTCGCGAGA',
        'CTGAGTGTGTGCACGTGAGATGCAC',
        'ACAGCGCGATAGATGTACATAGTAC',
        'TCTGCGTACGTGTCGACACATGTCT',
        'TGCGAGACATATAGTAGCTCTCTCA',
        'GTCGTCTAGATGTATCGATGTGATG',
        'CGTCGACGTCGTACACTGTACGTCA',
        'ACACTGCTGTGTCGTCAGTCTCACG',
        'ATGACACGCAGTCAGACGCGTGTGA',
        'TCTCACGCGCACGATCTCTGACACT',
        'ACATGTATCTGTAGATGCTGTGAGC',
         ]

RS_redundancy = [ 0.28 ]*18

conv_m = ([6, 8, 11]*3)*2

conv_r = ([3]*3+[5]*3+[7]*3)*2

pad = [
        True,
        False,
        False,
        False,
        False,
        False,
        True,
        False,
        True,
        ] + \
      [
        False,
        True,
        True,
        True,
        False,
        True,
        False,
        False,
        False,
        ]

bytes_per_oligo = ([18]*3 + [20]*3 + [22]*3) + ([16]*3+[20]*3+[20]*3)

DATA_PATH='./'
infile_name = [DATA_PATH+'data_files.tar.bz2.enc.1']*9 + [DATA_PATH+'data_files.tar.bz2.enc.2']*9
for i in range(18):
    print('-------------------------------------\n')
    print('i',i)
    print('bytes_per_oligo',bytes_per_oligo[i])
    print('RS_redundancy',RS_redundancy[i])
    print('conv_m',conv_m[i])
    print('conv_r',conv_r[i])
    print('pad',pad[i])

    if (i < 9):
        print('using encode with 1 CRC')
        helper.encode(data_file = infile_name[i], oligo_file = DATA_PATH+'reads.'+str(i), bytes_per_oligo = bytes_per_oligo[i], RS_redundancy = RS_redundancy[i], conv_m = conv_m[i], conv_r = conv_r[i], pad=pad[i])
    else:
        print('using encode with 2 CRCs')
        helper.encode_2crc(data_file = infile_name[i], oligo_file = DATA_PATH+'reads.'+str(i), bytes_per_oligo = bytes_per_oligo[i], RS_redundancy = RS_redundancy[i], conv_m = conv_m[i], conv_r = conv_r[i], pad=pad[i])

    with open(DATA_PATH+'reads.'+str(i)) as f_reads, open(DATA_PATH+'oligos_'+str(i)+'.fa','w') as f_oligos:
        for j, line in enumerate(f_reads):
            f_oligos.write('>oligos_'+str(i)+'_'+barcode_start[i]+'_'+barcode_end[i]+'_'+str(j)+'\n')
            f_oligos.write(barcode_start[i]+line.rstrip('\n')+barcode_end[i]+'\n')

