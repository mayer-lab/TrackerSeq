import sys
import re
import numpy as np
import pandas as pd

# This script extracts read information from the fastq out of the bbduk run (which combines R2 lineage barcode reads with R1 cellbc and umi information)
# and saves into 3 useful formats for subsequent pipeline analysis. 


# INPUTS: 
EXTRACTED_FASTQ = str(sys.argv[1]) # path to extracted fastq files containing lineage barcodes

# OUTPUTS
# Fasta handle without suffix, now added within the code due to 3 different ouputs!
REFORMATTED_FASTA = str(sys.argv[2]) # path to reformatted fasta fileRE
LIB = str(sys.argv[3]) # string indicating which library the file came from - does not require a comma, as in prev versions!

def find_cbc(s):
    try:
        pattern_cbc = r'_[ACTG]{16}_'
        pattern_extract_cbc = r'[ACTG]{16}'
        cbc_w_handles = re.findall(pattern_cbc, s)
        if len(cbc_w_handles) == 1:
            cbc_extract = re.findall(pattern_extract_cbc, cbc_w_handles[0])
        else: 
            cbc_extract = ['cbc_not_found']
        return cbc_extract[0]
    except ValueError:
        return ""


def find_umi(s):
    try:
        pattern_umi = r'_[ACTG]{10} |[ACTG]{12} '
        pattern_extract_umi = r'[ACTG]{10}|[ACTG]{12}'
        umi_w_handles = re.findall(pattern_umi, s)
        if len(umi_w_handles) == 1:
            umi_extract = re.findall(pattern_extract_umi, umi_w_handles[0])
        elif len(umi_w_handles) > 1:
            umi_extract = ['multiple_umi_matches_check_format']
        else: 
            umi_extract = ['umi_not_found']
        return umi_extract[0]
    except ValueError:
        return ""



def find_lbc(s):
    try:
        pattern_extract_lbc = r'[ACTG]{37}'
        lbc_extract = re.findall(pattern_extract_lbc, s)
        if lbc_extract==[]: 
            lbc_extract = ['lbc_not_found']
        return lbc_extract[0]
    except ValueError:
        return ""


file = open(EXTRACTED_FASTQ, 'r')
barcodeList = file.readlines()

temp_list = []
barcodeList_reformatted = barcodeList
for i in range(0, len(barcodeList), 4):
    temp_list.append([(">" + LIB), find_cbc(barcodeList[i]), find_umi(barcodeList[i]), find_lbc(barcodeList[i+1])])
    barcodeList_reformatted[i] = ">" + LIB + "," + find_cbc(barcodeList[i]) + "," + find_umi(barcodeList[i])


barcode_df_to_csv = pd.DataFrame(temp_list, columns = ['experiment', 'cellbc', 'umi', "LBC"])
# remove read rows without a lbc or cell-bc
clean_barcode_df_to_csv = barcode_df_to_csv[barcode_df_to_csv['cellbc']!="cbc_not_found"]
clean_barcode_df_to_csv = clean_barcode_df_to_csv[clean_barcode_df_to_csv['umi']!="umi_not_found"]
clean_barcode_df_to_csv = clean_barcode_df_to_csv[clean_barcode_df_to_csv['LBC']!="lbc_not_found"]


# write the legacy text output in .fa format for pipeline step 3
with open(REFORMATTED_FASTA + '.fa',"w") as f:
    for i in range(0,len(barcodeList_reformatted), 4):
        f.write(str(barcodeList_reformatted[i]) + "\n" + str(barcodeList_reformatted[i+1]))


# write the numpy df as a .csv for alternative analysis
np.savetxt(REFORMATTED_FASTA +'.csv', clean_barcode_df_to_csv, delimiter=",", fmt='%s')


# reorder columns and write only LBCs, umis to a .csv for bartender clustering (whitelisting)
barcode_df_for_bartender = clean_barcode_df_to_csv[['LBC', 'umi']]

# write the numpy df as a .csv for alternative analysis (e.g. bartender v1.1 umi clustering)
np.savetxt(REFORMATTED_FASTA + '_for_bartender.csv', barcode_df_for_bartender, delimiter=",", fmt='%s')



















