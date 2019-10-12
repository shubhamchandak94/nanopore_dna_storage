import sys
import matplotlib.pyplot as plt

if len(sys.argv) != 3:
    print('Incorrect number of arguments')
    sys.exit(1)

prefix = str(sys.argv[1])
num_mapped = int(sys.argv[2])

f_csv = open(prefix+'.error_stats.csv','w')

subs_pos = []
subs_rate = []
ins_pos = []
ins_rate = []
del_pos = []
del_rate = []

f_csv.write('subs_pos,subs_rate\n')
with open(prefix+'.substitution_stats') as f:
    for line in f:
        l = line.split()
        subs_pos.append(int(l[0]))
        subs_rate.append(sum([int(a) for a in l[1:]]))
        subs_rate[-1] /= num_mapped
        f_csv.write(str(subs_pos[-1])+','+str(subs_rate[-1])+'\n')

f_csv.write('ins_pos,ins_rate\n')
with open(prefix+'.insertion_stats') as f:
    for line in f:
        l = line.split()
        ins_pos.append(int(l[0]))
        ins_rate.append(int(l[1]))
        ins_rate[-1] /= num_mapped
        f_csv.write(str(ins_pos[-1])+','+str(ins_rate[-1])+'\n')

f_csv.write('del_pos,del_rate\n')
with open(prefix+'.deletion_stats') as f:
    for line in f:
        l = line.split()
        del_pos.append(int(l[0]))
        del_rate.append(int(l[1]))
        del_rate[-1] /= num_mapped
        f_csv.write(str(del_pos[-1])+','+str(del_rate[-1])+'\n')

f_csv.close()
