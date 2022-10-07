import os
import scanpy as sc
import tangram as tg
import pandas as pd
import numpy as np
from pathlib import Path

Path("data/gene_exprs/SlideSeq").mkdir(parents=True, exist_ok=True)
Path("data/gene_exprs/STARmap").mkdir(parents=True, exist_ok=True)

for dataset in ['SlideSeq', 'STARmap']:
    if dataset == 'SlideSeq':
        ##### Load the spatial data:  Slide-seq data, 9852 spatial voxels, 24518 genes
        print('Reading SlideSeq spatial data ...')
        path_to_sp = 'data/slideseq_MOp_1217.h5ad'
    else:
        ##### Load the spatial data: STARmap 1549 voxels with 1020 genes
        print('Reading STARmap spatial data ...')
        path_to_sp = 'data/STARmap.h5ad'
    ad_sp = sc.read_h5ad(path_to_sp)
    print(ad_sp)

    if dataset == 'SlideSeq':
        ##### Load single cell data, MOp 10Xv3 dataset, 26431 profiled cells with 27742 genes
        print('Reading MOp 10Xv3 single cell data ...')
        path_to_sc = 'data/mop_sn_tutorial.h5ad'
        print('Reading marker genes ...')
        df_genes = pd.read_csv('data/MOp_markers.csv', index_col=0)
        markers = np.reshape(df_genes.values, (-1,))
        markers = list(markers)
        print(len(markers))
    else:
        ##### Load SMARTSeq2 sn data, 45768 profiled cells with 23178 genes
        print('Reading MOp 10Xv3 single cell data ...')
        path_to_sc = 'data/SMARTSeq2.h5ad'
        markers = None
    ad_sc = sc.read_h5ad(path_to_sc)
    print(ad_sc)
    # normalize the number of counts within each cell to a fixed  number
    print('Normalizing counts ...')
    sc.pp.normalize_total(ad_sc)
    ##### Prepare to map
    print('Preparing for mapping ...')
    tg.pp_adatas(ad_sc, ad_sp, genes=markers)

    #Load mapping
    path_maps = f"data/ad_maps/{dataset}/"
    all_objects = os.listdir(path_maps)
    path_ge = f"data/gene_exprs/{dataset}/"
    print(f"Getting mappings from {path_maps}, saving gene expression predictions to {path_ge}")
    for i in range(len(all_objects)):
        print(i)
        map_obj = sc.read_h5ad(path_maps + all_objects[i])
        ad_ge = tg.project_genes(adata_map=map_obj, adata_sc=ad_sc)
        ad_ge.write(filename=path_ge + f"ge{i}.h5ad")