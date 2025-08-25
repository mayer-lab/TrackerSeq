import sys
import numpy as np, matplotlib.pyplot as plt, networkx as nx, pickle, json, gzip
import pandas as pd
import itertools
import statistics as stats
#import seaborn as sns
import os



# INPUTS
CELL_BCS_PATH = str(sys.argv[1]) #input cell barcode whitelist from preprocessing 
BAR_REF = str(sys.argv[2]) # input reformatted lineage barcode file
BARTENDER_BARCODES_PCR = str(sys.argv[3]) # input raw barcodes and cluster ids
BARTENDER_CLUSTERS_PCR = str(sys.argv[4]) # input cluster centers and quality scores

# OUTPUTS
CLONE_MAT_CSV = str(sys.argv[5])+'.csv' # output for cell - lineage barcode matrix in csv
CLONE_MAT_NPY = str(sys.argv[5])+'.npy' # output for cell - lineage barcode matrix in native npy
LBC_TXT = str(sys.argv[6]) + '_LBC_list.txt' # output for lineage barcode sequences as txt
CBC_CSV = str(sys.argv[6]) + '_cellbc_list.csv'# output for lineage barcode sequences in csv
BAR_MAP = str(sys.argv[6]) + '_barcode_map.npy' # ouptut for barcode mapping file

N_READS = int(sys.argv[7])
N_UMIS = int(sys.argv[8])
N_HAMMING = 5 
SMALL_DATASET = False
CLUSTER_QUALITY_CUTOFF = 0.8


# cell_bcs = open(CELL_BCS_PATH).read().strip('\n').split('\n') #if the barcode is extracted from Seurat object
cell_bcs = pd.read_csv(CELL_BCS_PATH, sep='\t', header = None)
cell_bcs = cell_bcs.loc[:,0].tolist()

# if inputs are available, read in the bartender clustering file outputs (bartender_single_com) for collapsing low hamming distance barcodes
if os.path.exists(BARTENDER_BARCODES_PCR):
    bartender_barcodes = pd.read_csv(BARTENDER_BARCODES_PCR, sep=',', header = 'infer')
if os.path.exists(BARTENDER_CLUSTERS_PCR):
    bartender_clusters = pd.read_csv(BARTENDER_CLUSTERS_PCR, sep=',', header = 'infer')

# For each unique barcode in the bartender_barcodes file, match it with the barcode center sequence based on cluster id
# This generates a whitelist (default hamming = 3) which can be used to map barcodes later in this pipeline step
# This approach is superior than the current implementation in this pipeline for mapping barcodes since it uses umis and does iterative clustering

# Filter for quality threshold, -P*log2(P) - (1 - P)*log2(1-P) where P is the % majority sequence at each NT position
# 0.5 is a decent default as this means ~9/10 instances of a barcode have a NT value in a given position
# 0.8 would be 1/4
# 0.9 1/3... 
# likely doesnt matter below this level since umi filtering will remove anything below (1/N_UMI) and N_READS filtering will help as well

bartender_whitelist = {}
for i, unique_bc in enumerate(bartender_barcodes['Unique.reads']): 
    cluster = bartender_barcodes['Cluster.ID'][i]
    cluster_quality = bartender_clusters['Cluster.Score'][bartender_clusters['Cluster.ID']==cluster].values[0]
    if cluster_quality < CLUSTER_QUALITY_CUTOFF:
        center = bartender_clusters['Center'][bartender_clusters['Cluster.ID']==cluster].values[0]
        bartender_whitelist[unique_bc] = center

print('\nWhitelist maps '+repr(len(set(bartender_whitelist.keys()))) + ' barcodes to ' +repr(len(set(bartender_whitelist.values()))) + ' barcodes')

# Save the whitelist to the BAR_MAP output file
np.save(BAR_MAP, bartender_whitelist)


# Process reads into memory while checking for valid data
counts = {}
f = open(str(BAR_REF))
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

# format of counts as dict {(experiment, cellbc, umibc, lineagebc): umi_count}
all_lbcs = [k[3] for k,v in counts.items()]
all_cells = [k[1] for k in counts.keys()]
low_counts = {k:v for k,v in counts.items() if v <= 2}
low_count_lbcs = [k[3] for k,v in low_counts.items()]

num_reads = [v for k,v in counts.items()]
# #Plot readcount histogram - this is the # of reads for each unique umi, cellbc, lineage bc combination
# plt.hist(num_reads, bins = 3000)
# plt.plot([N_READS, N_READS],[1,10**6],'-k',linewidth=2)
# plt.xticks()
# plt.yticks()
# plt.text(10*1.5,10**5*.8,'N_READS cutoff', fontsize=12)
# plt.yscale('log')
# plt.xscale('log')

# Create function to calculate hamming distance between pair of barcodes
def hamming(bc1,bc2): return np.sum([x1 != x2 for x1,x2 in zip(bc1,bc2)])

# Import Levenshtein distance function -- probably a better function for calculating distances because it can handle all reasonable types of edits
import time
#from Levenshtein import distance

# Barcode reads filtered for minimum number of reads, sweeping across reasonable values
READ_sweep = range(1, 31)
Retained = list()
counts_filtered_sweep = list()
for r in range(len(READ_sweep)):
    counts_filtered = {k:v for k,v in counts.items() if v >= READ_sweep[r]}
    counts_filtered_sweep.append((counts_filtered))
    Retained.append(len(counts_filtered))
    print("With "+ repr(READ_sweep[r]) +', retaining '+repr(len(counts_filtered))+ ' out of '+repr(len(counts))+' (Sample,cellBC,UMI,GFP-LBC) combinations')


# # Plot the retention of barcodes based on N_reads cutoff
# plt.plot(READ_sweep,Retained)
# plt.plot([N_READS,N_READS],[np.min(Retained),np.max(Retained)],'-k',linewidth=2)
# plt.text(10*1.1,np.max(Retained)*.8,'N_READS cutoff',fontsize=12)
# plt.title("N_READs sweep - # Barcodes retained")
# plt.yscale('log')

# Depending on computational complexity, set the counts_filtered list to the barcodes retained using imputed N_reads cutoff (addl reporting will be available in log file)
if SMALL_DATASET == True: 
    # If dataset size is smaller, can compute totals and distance metrics for all reads prior to filtering (at BC merge, removal steps)
    counts_filtered = counts_filtered_sweep[1-min(READ_sweep)]
else: 
    # Cutoff all reads below N_READS threshold to reduce computational load, which can be significant! Default N_READS = 10
    counts_filtered = counts_filtered_sweep[3-min(READ_sweep)]



#Convert counts filtered to data frame with named columns
temp_counts = []
for a,b,c,d in counts_filtered:
    temp_counts.append([a,b,c,d,counts_filtered[(a,b,c,d)]])
counts_filtered_df = pd.DataFrame(temp_counts, columns = ['experiment', 'cellbc', 'umi', "LBC", "umi_count"])

# Generate sums for unique lineage barcodes (LBCs)
all_LBC = np.unique(counts_filtered_df["LBC"])
lbc_totals_list = []
lbc_total_reads = {}
lbc_total_umis = {}
reads_x_umis = {}
cells_per_lbc = {}
for lbc in all_LBC:
    lbc_array = counts_filtered_df[counts_filtered_df["LBC"]==lbc]
    umis = len(lbc_array)
    lbc_total_umis[lbc] = umis
    totes = sum(lbc_array["umi_count"])
    lbc_total_reads[lbc] = totes
    lbc_totals_list.append([lbc, umis, totes])
    reads_x_umis[lbc] = umis*totes
    cells_per_lbc[lbc] = len(np.unique(lbc_array['cellbc']))


lbc_read_counts = [v for k, v in lbc_total_reads.items()] 
# #Plot histogram of reads per lbc
# plt.hist(lbc_read_counts, bins =10000)
# plt.yscale('log')
# #plt.xscale('log')
# #plt.xlim(0,200)
# plt.xlim(0,100000)


lbc_cell_counts = [v for k,v in cells_per_lbc.items() if v<1000]
# #Plot histogram of cells per lbc
# plt.hist(lbc_cell_counts, bins=1000)
# plt.yscale('log')
# #plt.xscale('log')
# #plt.xlim(0,200)
# plt.xlim(0,100)

umi_totals = [v for k,v in lbc_total_umis.items() if v<1000]
# #Plot histogram of umis per lbc
# plt.hist(umi_totals, bins=1000)
# plt.yscale('log')
# #plt.xscale('log')
# #plt.xlim(0,200)
# plt.xlim(0,100)

# #This section performs Levenstein distance calculation and maps barcode pairs with low lev distance to each other based on the barcode with more umis
# #Note this is not the current approach for mapping barcodes, but does provide the initial pairwise distance calculations for graphing barcode similarity
# #If implemented, this would require an iterative mapping approach in order to adequately cluster barcodes 
# #(e.g. the same block of code seen twice below would be repeated n times until all barcodes met the minimum distance threshold)
unique_barcodes = np.unique(counts_filtered_df['LBC'])
good_gfp_bcs = []
lev_list = []
lev_values_list = []
initial_bc_map = {}
# Examine all pairwise barcode sequences to determine similarity -- the stats are still used from this block, 
# however the barcode map is deprecated, instead the bartender generated clustering files are utilized
for i,bc1 in enumerate(unique_barcodes):
    if i > 0 and i % 500 == 0: print('Mapped '+repr(i)+' out of '+repr(len(unique_barcodes))+' barcodes')
    mapped = False
    for bc2 in good_gfp_bcs:
        lev_dist = hamming(bc1,bc2)
        lev_list.append([bc1, bc2, lev_dist])
        lev_values_list.append(lev_dist)
        if lev_dist <= N_HAMMING:
            #mapped = True  #removed this since it can lead to 2 pairwise barcodes not being compared
            if lbc_total_umis[bc1]>lbc_total_umis[bc2]:
                initial_bc_map[bc2] = bc1
            elif lbc_total_umis[bc1]<lbc_total_umis[bc2]:
                initial_bc_map[bc1] = bc2
            else:
                initial_bc_map[bc1] = bc2
            break
    if not mapped:
        good_gfp_bcs.append(bc1)

mapped_barcodes = list(set([v for v in initial_bc_map.values()]))
second_barcodes = []
lev_list_mapped = []
lev_values_list_mapped = []
second_bc_map = {}
for i,bc1 in enumerate(mapped_barcodes):
    if i > 0 and i % 500 == 0: print('Mapped '+repr(i)+' out of '+repr(len(mapped_barcodes))+' barcodes')
    mapped = False
    for bc2 in second_barcodes:
        lev_dist = hamming(bc1,bc2)
        lev_list.append([bc1, bc2, lev_dist])
        lev_values_list.append(lev_dist)
        if lev_dist <= N_HAMMING:
            #mapped = True  #removed this since it can lead to 2 pairwise barcodes not being compared
            if lbc_total_umis[bc1]>lbc_total_umis[bc2]:
                second_bc_map[bc2] = bc1
            elif lbc_total_umis[bc1]<lbc_total_umis[bc2]:
                second_bc_map[bc1] = bc2
            else:
                second_bc_map[bc1] = bc2
            break
    if not mapped:
        second_barcodes.append(bc1)


# # Plot histogram of the Levenshtein distance metrics for each bc pair
# plt.hist(lev_values_list, bins=25)
# plt.xticks()
# plt.yticks()
# plt.title('levenshtein distances for ' + str(len(lev_values_list)) + ' pairwise comparisons')
# plt.yscale('log')


# Filter for N_READS and N_UMIS as a cross product
# Filter for N_UMIS alone
# Combine the retention list to both conditions are met

retained_barcodes_reads = [k for k,v in lbc_total_reads.items() if v>N_READS] # e.g. N_reads = 5 * N_umi = 3 = 15
retained_barcodes_umis = [k for k,v in lbc_total_umis.items() if v>N_UMIS] # e.g. N_reads = 5 * N_umi = 3 = 15
retained_barcodes = [x for x in retained_barcodes_reads if x in retained_barcodes_umis]


# Implementation to merge LBCs at or below the hamming threshold into a single barcode based on clustering performed in bartender
# Filters the counts dictionary to include only barcodes meeting read + umi thresholds BEFORE MAPPING BARCODES
# Uses bartender_whitelist produced from bartender_single_com
counts_filtered_retained = {k:v for k,v in counts_filtered.items() if k[3] in retained_barcodes}
print('\nFiltered '+repr(len(set(counts_filtered.keys()))) + ' reads to ' +repr(len(set(counts_filtered_retained.keys()))) + ' reads meeting N_READ and N_UMI standards')
collapsed_counts = {}
for k,v in counts_filtered_retained.items():
    if k[3] in bartender_whitelist.keys(): 
        new_bc = bartender_whitelist[k[3]]
        new_key = [k[0], k[1], k[2], new_bc]
        new_key = tuple(new_key)
        collapsed_counts[new_key] = v
    else: 
        collapsed_counts[k] = v

collapsed_bcs = list(set([k[3] for k,v in collapsed_counts.items()]))
collapsed_cells = list(set([k[1] for k,v in collapsed_counts.items()]))
print('\nCollapsed '+repr(len(set(retained_barcodes))) + ' barcodes to ' +repr(len(set(collapsed_bcs))) + ' barcodes')
print('\nThere are now ' + repr(len(collapsed_bcs)) + ' barcodes in ' + repr(len(collapsed_cells)) + ' cells') 



# From the now mapped and filtered counts dictionary, extract the names of the sample libraries and the unique cells
lib_names_filtered = [k[0] for k in collapsed_counts.keys()]
cell_bcs_filtered = [k[1] for k in collapsed_counts.keys()]

# Initialize dicts for cell-specific data and filtered cell-specific data for outputs
cell_data = {}
cell_data_umi_filtered = {}
for lib,cell in zip(lib_names_filtered, cell_bcs_filtered):
    if cell !='cbc_not_found':
        cell_data[(lib,cell)] = {}
        cell_data_umi_filtered[(lib,cell)] = {}


# Aggregate UMI counts per cellbc - LBC combination
for lib,cell,umi,BC in collapsed_counts.keys():
    if (lib,cell) in cell_data:
        if not BC in cell_data[(lib,cell)]:
            cell_data[(lib,cell)][BC] = 0
        cell_data[(lib,cell)][BC] += 1

# Filter each cell based on min N_UMIS
for lib,cell in cell_data.keys():
        for BC in cell_data[(lib,cell)].keys():
            if cell_data[(lib,cell)][BC] > N_UMIS:
                cell_data_umi_filtered[(lib,cell)][BC] = cell_data[(lib,cell)][BC]

# Sort the LBC list based on the order of the overall data frames for export - the LBC list order will correspond to columns of the clone matrix
i = 0
clone_mat = np.zeros((len(cell_data_umi_filtered.keys()),len(set(collapsed_bcs))))
final_cell_bcs = [c[1] for c in cell_data_umi_filtered.keys()]
for i, k in enumerate(cell_data_umi_filtered.keys()):
    lbcs_in_cell = [h for h in cell_data_umi_filtered[k].keys()]
    for j in lbcs_in_cell:
        lbc_idx = collapsed_bcs.index(j)
        clone_mat[i, lbc_idx] = 1
clone_mat_df = pd.DataFrame(clone_mat, index=final_cell_bcs, columns=collapsed_bcs, dtype=int)
#clone_mat = np.array(clone_mat,dtype=int)

# Write out
clone_mat_df.to_csv(CLONE_MAT_CSV, index=True, header=True, sep=',')
np.save(CLONE_MAT_NPY, clone_mat_df)

open(LBC_TXT,'w').write('\n'.join(retained_barcodes));
np.savetxt(CBC_CSV, final_cell_bcs , delimiter=',' , fmt = "%s")




