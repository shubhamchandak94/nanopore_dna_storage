import helper

barcode_start = [
        'CTGGCTCCTCTGTATGTTGGAGAAT',  
        'TGCGGATGCGGAAGTATGGTCCTCG', 
        'AGTAACGCCTATTGATAACGAAGCA',
        'CTGGCGGCCTTGGCCGACTATCTGC',  
        'TAGTCCGCGCTCGAATTCCGAGGCC',  
        'ATGTTCGGAACGTCAAGACCGAGGA', 
        'GCTAGTACGCGAACAGAGTGCAGTA',  
        'CACCTGTGCTGCGTCAGGCTGTGTC',  
        'CGTACAATCGTATTAGGCACCTTCC',  
        'GTATACATTCCTTGCCAACATAGTA',  
        'TATCGATTGCATGATACATCCGCAC',
        'GGCCTACCGAGGACCGCTTAGTAGG', 
        'GATACTATCGAGATTACTCCAAGTC',
        ]

barcode_end = [
         'CCTATATGTACCTCTATCGTAAGTC',
         'CACTAGAAGCATGTCGCTATCGAGT',
         'TAACCTTCGCTGCTAGGAACTGTCT',
         'ACCATGTCGTACAGTCGTTGTAACA',
         'TACAAGACTACGCAAGATCGCGCTA',
         'TGGCTCCATTATGCTACAATCACTA',
         'ACAGATGCAGTAATTCTCACGAACT',
         'GCTGTCCGTTCCGCATTGACACGGC',
         'GCGGACCTCCAGATCCACTTGTCTG',
         'TGAATCTGGATACGCGTTCCTCAAC',
         'GACCTGTGGAAGTTCCTCATTACTA',
         'CCTATCATGAATTAGATGCTTGGAC',
         'GCTAGTCGATCCTCTGCTGCAATCG',
         ]

RS_redundancy = [
        0.3,
        0.3,
        0.3,
        0.3,
        0.3,
        0.3,
        0.3,
        0.3,
        0.3,
        0.2,
        0.4,
        0.3,
        0.3,
        ]

conv_m = [
        8,
        11,
        14,
        8,
        11,
        14,
        8,
        11,
        14,
        11,
        11,
        11,
        11,
        ]

conv_r = [
        1,
        1,
        1,
        3,
        3,
        3,
        5,
        5,
        5,
        3,
        3,
        3,
        3,
        ]

pad = [
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        False,
        True,
        False,
        False,
        False,
        False,
        ]

bytes_per_oligo = [
        10,
        10,
        10,
        18,
        18,
        18,
        20,
        20,
        20,
        18,
        18,
        18,
        18,
        ]

DATA_PATH='../nanopore_dna_storage_data/'
infile_name = DATA_PATH+'encoded_file/data_files.tar.bz2.enc'
for i in range(13):
    print('i',i)
    print('bytes_per_oligo',bytes_per_oligo[i])
    print('RS_redundancy',RS_redundancy[i])
    print('conv_m',conv_m[i])
    print('conv_r',conv_r[i])
    print('pad',pad[i])
    helper.encode(data_file = infile_name, oligo_file = DATA_PATH+'oligos_files/reads.'+str(i), bytes_per_oligo = bytes_per_oligo[i], RS_redundancy = RS_redundancy[i], conv_m = conv_m[i], conv_r = conv_r[i], pad=pad[i])
    with open(DATA_PATH+'oligo_files/reads.'+str(i)) as f_reads, open(DATA_PATH+'oligo_files/oligos_'+str(i)+'.fa','w') as f_oligos:
        for j, line in enumerate(f_reads):
            f_oligos.write('>oligos_'+str(i)+'_'+barcode_start[i]+'_'+barcode_end[i]+'_'+str(j)+'\n')
            f_oligos.write(barcode_start[i]+line.rstrip('\n')+barcode_end[i]+'\n')

