
import sys
import numpy as np, pandas as pd
from matplotlib import pyplot as plt
from matplotlib import colors
from matplotlib import figure


CLONE_MAT = str(sys.argv[1]) + '.csv'
CLONE_NPY = str(sys.argv[1]) + '.npy'
CELL_BCS = str(sys.argv[2])
DATASET = str(sys.argv[3]) + '.csv'


# load csv file as pandas and numpy file to do boolean filtering
clone_mat = pd.read_csv(CLONE_MAT, on_bad_lines='skip', header = 0, index_col = 0)
clone_mat_np = np.load(CLONE_NPY)
# read in cell bc file
# cell_bcs = open(CELL_BCS).read().strip('\n').split('\n') #use this to read in barcodes that are already in .txt file format
cell_bcs = pd.read_csv(CELL_BCS, sep='\t', header = None)


# examine dimensions of the matrix
clone_mat.shape


# get rid of columns that label no clones. OR take out any columns where sum < 1 (i.e. does not label cell) 
multi_cc = np.sum(clone_mat_np, axis = 0)
multi_cc_bool = multi_cc >= 1
clone_mat = clone_mat.iloc[:,multi_cc_bool]
clone_mat_np = clone_mat_np[:,multi_cc_bool]
clone_mat.head()
clone_mat.shape


# only keep rows (cells) that sum to up to 1 or more (i.e. must contain a lineage barcode)
sum_bc = np.sum(clone_mat_np, axis = 1)
sum_bc_bool = sum_bc >= 1
clone_mat = clone_mat.iloc[sum_bc_bool, :]
clone_mat_np = clone_mat_np[sum_bc_bool,:]
clone_mat.shape
clone_mat.head()


clone_mat.shape


clone_mat.to_csv(DATASET, index_label = 'cellbc')




