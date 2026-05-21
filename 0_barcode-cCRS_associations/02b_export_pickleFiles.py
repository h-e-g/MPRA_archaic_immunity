#source /pasteur/helix/projects/evo_immuno_pop/MPRA/MPRAflow/conda/envs/mpraflow_py36/bin/deactivate

import sys
import pickle
import pandas as pd
from collections import Counter
from collections import defaultdict
import pyarrow
from matplotlib import pyplot as plt
import seaborn
import re

runName=sys.argv[1]
prefix='/pasteur/helix/projects/evo_immuno_pop/MPRA/MPRAflow/outs/'+runName+'/'+runName
####### export non filtered associations to a txt file
Associations=prefix+'_coords_to_barcodes.pickle'

dic=pickle.load(open(Associations,'rb'))
AssociationsTxt=prefix+'_coords_to_barcodes.tsv'
f=open(AssociationsTxt,'w')
for (oligo_name,sublist) in dic.items():
  for barcode in sublist:
    f.write(oligo_name+'\t'+barcode+'\n')

f.close()

####### export filtered associations to a txt file
## Deactivated: MPRAflow crashed before doing the filtering. filtering will be done directly inside R

# Associations=prefix+'_filtered_coords_to_barcodes.pickle'
#
# dic=pickle.load(open(Associations,'rb'))
# AssociationsTxt=prefix+'_filtered_coords_to_barcodes.tsv'
# f=open(AssociationsTxt,'w')
# for (oligo_name,sublist) in dic.items():
#   for barcode in sublist:
#     f.write(oligo_name+'\t'+barcode+'\n')
#
# f.close()


####### export non filtered associations to a txt file

# count associations
# raw_counts=(pd.Series(dic, name = 'n_barcodes').str.len().rename_axis('coord').reset_index())
# print('filter')
# print(raw_counts)

# all_bcs=(dic.values())

# allpairs = [(name,item) for (name,sublist) in dic.items() for item in sublist]
# flat_BC = [item for sublist in all_bcs for item in sublist]
