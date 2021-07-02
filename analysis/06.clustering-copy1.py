#!/usr/bin/env python
# coding: utf-8
# UMAP-LEARN NEEDS AN OLDER ENV - scf_terra (for example)
# ENV - routine not compatible with this

import datetime
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import logging

import sys
sys.path.insert(0, '/cndd2/fangming/projects/SingleCellRoutines')
from __init__scr import * 
import utils
import clst_utils

log = utils.create_logger()

np.random.seed(0)
today = datetime.date.today()

# samples = [
#     'Slice1_Replicate1',
#     'Slice1_Replicate2',
#     'Slice1_Replicate3',
    
#     'Slice2_Replicate1',
#     'Slice2_Replicate2',
#     'Slice2_Replicate3',
    
#     'Slice3_Replicate1',
#     'Slice3_Replicate2',
#     'Slice3_Replicate3',
# ]

file = '../data/processed_vizgen_merfish_may1_2021.h5ad'
sample = 'Slice2_Replicate1'
output = '../data/clustering_embedding_{}_{}.tsv.gz'.format(sample, today)
# clustering embedding
npc = 50
leiden_knn = 30
leiden_resolution = 1
umap_knn = 30
umap_min_dist = 0.1

logging.info(
    '''
    sample: {}
    Number of PCs: {}
    Leiden kNN: {}
    Leiden resolution: {}
    UMAP kNN: {}
    UMAP min dist: {}
    '''.format(
        sample,
        npc, 
        leiden_knn, 
        leiden_resolution,
        umap_knn,
        umap_min_dist,
        )
    )

thedata = pd.read_hdf(file, 'data_'+sample)
# informations
metacols = thedata.columns[:11]
genes = thedata.columns[11:]
gmat = thedata.iloc[:,11:]
print(metacols, len(genes), gmat.shape)

logging.info("dataset readed in")

# # gmats are normalized
# print(gmat.sum(axis=1))
X = (PCA(n_components=npc, 
        svd_solver='randomized', 
        random_state=0)
    .fit_transform(gmat.values)
    )
logging.info("PCA done")

df_umap = clst_utils.run_umap_lite(
    X, gmat.index.values,
    n_neighbors=umap_knn,
    min_dist=umap_min_dist,
    random_state=0,
    )
logging.info("UMAP done")

df_clst = clst_utils.clustering_routine(
    X, gmat.index.values, 
    leiden_knn, 
    resolution=leiden_resolution,
   )
logging.info("Clustering done")
df_clst = 'C'+df_clst.astype(str)
df_res = (df_clst.join(df_umap))
df_res.index.name = 'cell'

# save results
df_res.to_csv(output, sep='\t')
logging.info("Output file saved {}".format(output))