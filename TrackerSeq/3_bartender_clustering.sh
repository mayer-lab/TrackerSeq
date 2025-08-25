# Note that bartender does not run on all systems -- particularly redhat linux

reformat_csv=$1
output=$2

### 1. clustering lineage barcodes based on umi and NT position information, using Bartender v1.1

bartender_single_com -f reformat_csv -o output -d 3 -z 5

