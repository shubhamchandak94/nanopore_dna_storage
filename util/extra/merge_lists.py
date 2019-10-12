import shutil
import os

original_list_dir='/raid/nanopore/shubham/20190629_nanopore_data/decoded_lists/11_3_RS_0.2/'
merged_list_dir='/raid/nanopore/shubham/20190629_nanopore_data/decoded_lists/11_3_RS_0.2/merged/'
num_thr = 15

num_done = 0
f_info = open(merged_list_dir+'info.txt','w')

for i in range(num_thr):
    with open(original_list_dir+'info_'+str(i)+'.txt') as f:
        lines = f.readlines()
        num_lines = len(lines)
        for j in range(num_lines):
            f_info.write(lines[j])
            src = original_list_dir+'list_'+str(i)+'_'+str(j)
            dst = merged_list_dir + 'list_'+str(num_done+j)
            if os.path.isfile(src):
                shutil.copy(src,dst)
        num_done += num_lines
