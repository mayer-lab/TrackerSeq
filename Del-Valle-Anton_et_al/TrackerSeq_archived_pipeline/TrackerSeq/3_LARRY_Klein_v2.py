
import sys
import numpy as np, matplotlib.pyplot as plt, networkx as nx, pickle, json, gzip
import pandas as pd
import itertools


N_READS = 10
N_UMIS = 9
N_HAMMING = 6
CELL_BCS_PATH = "<file path of whitelist>"
LIB_NAMES_PATH = "<file path of lib.txt>"



# cell_bcs = open(CELL_BCS_PATH).read().strip('\n').split('\n') #if the barcode is extracted from Seurat object
cell_bcs = pd.read_csv(CELL_BCS_PATH, sep='\t', header = None)
cell_bcs = cell_bcs.loc[:,0].tolist()
lib_names = open(LIB_NAMES_PATH).read().strip('\n').split('\n')



counts = {}
f = open(str(sys.argv[3]))
l = f.readline().strip('\n')
current_tag = []
i = 0
print('Reading in all barcodes')
while not (l == '' and len(current_tag)==0):
    i += 1
    if i % (3*10**6)==0: print('Processed '+repr(int(i/3))+' reads')
    if l == '':
        current_tag = []
    elif l[0] == '>':
        current_tag = l[1:].split(',')
    elif l != '' and len(current_tag)==3:
        current_tag.append(l)
        current_tag = tuple(current_tag)
        if not current_tag in counts: counts[current_tag] = 0
        counts[current_tag] += 1
        
    l = f.readline().strip('\n')



num_reads = [v for k,v in counts.items()]
plt.hist(np.log(num_reads)/np.log(10), bins=50)
plt.plot([np.log(N_READS)/np.log(10),np.log(N_READS)/np.log(N_READS)],[1,10**5],'-k',linewidth=2)
plt.xticks(range(6),np.logspace(0,5,6))
plt.text(np.log(N_READS)/np.log(10)*1.1,10**5*.8,'N_READS cutoff', fontsize=12)
plt.yscale('log')

counts_filtered = {k:v for k,v in counts.items() if v >= N_READS}
print('Retaining '+repr(len(counts_filtered))+ ' out of '+repr(len(counts))+' (Sample,Cell-BC,UMI,GFP-BC) combinations')



import time
from Levenshtein import distance

all_gfp_bcs = sorted(set([k[3] for k in counts_filtered]))

def collapse_barcodes(all_bcs, N_HAMMING=1):
    # Convert the barcodes to a numpy array of strings
    all_bcs_arr = np.array(all_bcs, dtype='U')
    
    # Initialize the list of good barcodes and the barcode map
    good_bcs = []
    bc_map = {}
    
    # Loop through all barcodes
    for i, bc1 in enumerate(all_bcs_arr):
        # Print status
        if i % 500 == 0:
            print('Mapping {} out of {} barcodes'.format(i, len(all_bcs_arr)))
            elapsed_time = time.time() - start_time
            remaining_time = elapsed_time * (len(all_bcs_arr) - i) / (i + 1)
            print('Estimated remaining time: {:.2f} seconds'.format(remaining_time))

        # Try to map the barcode to an existing good barcode
        mapped = False
        for bc2 in good_bcs:
            if distance(bc1, bc2) <= N_HAMMING:
                mapped = True
                bc_map[bc1] = bc2
                break
        if not mapped:
            good_bcs.append(bc1)
    
    return good_bcs, bc_map

# Example usage
print('Collapsing barcodes with Levenshtein distance {}...'.format(N_HAMMING))
start_time = time.time()
good_gfp_bcs, bc_map = collapse_barcodes(all_gfp_bcs, N_HAMMING=N_HAMMING)
end_time = time.time()
print('Done. Elapsed time: {:.2f} seconds'.format(end_time - start_time))

print('\nCollapsed '+repr(len(bc_map))+' barcodes')
for bc in good_gfp_bcs: bc_map[bc] = bc

cell_data = {}
for lib,cell in zip(lib_names,cell_bcs):
    cell_data[(lib,cell)] = {}

for lib,cell,umi,BC in counts_filtered.keys():
    if (lib,cell) in cell_data:
        if not BC in cell_data[(lib,cell)]:
            cell_data[(lib,cell)][BC] = 0
        cell_data[(lib,cell)][BC] += 1

BC_lists = []
for i in range(1,10):
    BC_list = []
    for lib,cell in zip(lib_names,cell_bcs):
        bc_counts = cell_data[(lib,cell)]
        valid_bcs = [bc_map[k] for k,v in bc_counts.items() if v >= i]
        BC_list.append(','.join(sorted(valid_bcs))) # added a comma to separate multiple-barcode 
    BC_lists.append(BC_list)

efficiency = np.array([len([ll for ll in l if len(ll)>0]) for l in BC_lists]) / len(cell_bcs)
plt.plot(range(1,10),efficiency)
plt.plot([N_UMIS,N_UMIS],[np.min(efficiency),np.max(efficiency)],'-k',linewidth=2)
plt.text(N_UMIS*1.1,np.max(efficiency)*.95,'UMI cutoff',fontsize=14)

final_BCs = BC_lists[N_UMIS-1]
print('\nFinal annotation has '+repr(len(set(final_BCs)))+' clones in '+ repr(len([l for l in final_BCs if len(l)>0]))+' cells')



BC_set = sorted(set([bc for bc in final_BCs if bc != '']))
BC_set = list(map(lambda x: x.split(','), BC_set))
BC_set = list(itertools.chain.from_iterable(BC_set))
clone_mat = np.zeros((len(final_BCs),len(BC_set)))
for i,bc in enumerate(final_BCs):
    if bc != '':
        multi_bc = bc.split(',')
        for n in range(len(multi_bc)):
            j = BC_set.index(multi_bc[n])
            clone_mat[i,j] = 1
clone_mat = np.array(clone_mat,dtype=int)
np.savetxt(str(sys.argv[4]),clone_mat,delimiter=',',fmt='%i');
np.save(str(sys.argv[5]),clone_mat);
open(str(sys.argv[6]),'w').write('\n'.join(final_BCs));
