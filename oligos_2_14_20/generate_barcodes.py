import numpy as np
import distance

next_base = {
            'A':['C','G','T'],
            'C':['A','G','T'],
            'G':['C','A','T'],
            'T':['C','G','A'],
            }

num_barcodes = 36 # 2 times because we need both start and end
barcode_len = 25
min_edit_distance = 12
min_AT_count = 11
max_AT_count = 14
np.random.seed(42)

barcodes = []
while True:
    if len(barcodes) == num_barcodes:
        break
    barcode = np.random.choice(['A','C','G','T'])
    for j in range(barcode_len-1):
        barcode += np.random.choice(next_base[barcode[-1]])
    if len(barcodes) == 0 or min([distance.levenshtein(barcodes[i],barcode) for i in range(len(barcodes))]) >= min_edit_distance:
        AT_count = barcode.count('A') + barcode.count('T')
        if min_AT_count <= AT_count <= max_AT_count:
            barcodes.append(barcode)

for barcode in barcodes:
    print(barcode)

edit_distances = np.zeros((num_barcodes,num_barcodes))

for i in range(num_barcodes):
    for j in range(num_barcodes):
        edit_distances[i][j] = distance.levenshtein(barcodes[i],barcodes[j])

print(edit_distances)
