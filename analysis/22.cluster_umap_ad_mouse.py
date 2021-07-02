#!/usr/bin/env python
# coding: utf-8
# UMAP-LEARN NEEDS AN OLDER ENV - scf_terra (for example)
# ENV - routine not compatible with this

"""FUTURE: 
- explore zscore genes 
"""

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

all_samples = [
    'Xulab_2_5_region_0',
    'Xulab_2_5_region_1',
    'Xulab_2_6_region_0',
    'Xulab_2_6_region_1',
]

samples_list = [
    all_samples,
]

for samples in samples_list:
    samples_shortname = "_".join([
        sample.replace('Slice', 'S')
              .replace('_Replicate', 'R')
              .replace('region_', 'R')
              .replace('Xulab_', 'Xu')
        for sample in samples
    ])

    input = '../data/processed_merfish_ad_mouse_june2_2021.hdf5'
    output = '../data/clustering_embedding_{}_{}.tsv.gz'.format(samples_shortname, today)
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
            samples,
            npc, 
            leiden_knn, 
            leiden_resolution,
            umap_knn,
            umap_min_dist,
            )
        )

    logging.info("reading datasets")
    gmat = []
    meta = []
    for sample in samples:
        _gmat = pd.read_hdf(input, 'mat_'+sample)
        gmat.append(_gmat)
        _meta = pd.read_hdf(input, 'meta_'+sample)
        _meta['sample'] = sample 
        meta.append(_meta)
    gmat = pd.concat(gmat)
    meta = pd.concat(meta)

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
    df_res['sample'] = meta.loc[df_res.index, 'sample']

    # save results
    df_res.to_csv(output, sep='\t')
    logging.info("Output file saved {}".format(output))