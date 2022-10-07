import numpy as np
import pandas as pd
import scanpy as sc
import tangram as tg
from pathlib import Path

Path("data/ad_maps/SlideSeq").mkdir(parents=True, exist_ok=True)
Path("data/ad_maps/STARmap").mkdir(parents=True, exist_ok=True)

# we follow the Tangram Tutorial without squidpy

##### Load the spatial data:  Slide-seq data, 9852 spatial voxels, 24518 genes
print('Reading SlideSeq spatial data ...')
ad_sp = sc.read_h5ad('data/slideseq_MOp_1217.h5ad')
print(ad_sp)

##### Load single cell data, MOp 10Xv3 dataset, 26431 profiled cells with 27742 genes
print('Reading MOp 10Xv3 single cell data ...')
ad_sc = sc.read_h5ad('data/mop_sn_tutorial.h5ad')
print(ad_sc)
# normalize the number of counts within each cell to a fixed  number
print('Normalizing counts ...')
sc.pp.normalize_total(ad_sc)

##### Prepare to map
print('Reading marker genes ...')
df_genes = pd.read_csv('data/MOp_markers.csv', index_col=0)
markers = np.reshape(df_genes.values, (-1, ))
markers = list(markers)
print(len(markers))
print('Preparing for mapping ...')
tg.pp_adatas(ad_sc, ad_sp, genes=markers)

##### Map
print('Mapping single cells onto space ...')
for i in range(9):
    # either using cpu or using GPU, I strongly recommend GPU, either from the server or from Google Collab
    ad_map = tg.map_cells_to_space(
        adata_sc=ad_sc,
        adata_sp=ad_sp,
        #device='cpu',
        device='cuda:0',
    )
    print('Saving map')
    ad_map.write(filename=f'data/ad_maps/SlideSeq/ad_map{i}.h5ad')

##### Load the spatial data: STARmap 1549 voxels with 1020 genes
print('Reading STARmap spatial data ...')
ad_sp = sc.read_h5ad('data/STARmap.h5ad')
print(ad_sp)

##### Load SMARTSeq2 sn data, 45768 profiled cells with 23178 genes
print('Reading MOp 10Xv3 single cell data ...')
ad_sc = sc.read_h5ad('data/SMARTSeq2.h5ad')
print(ad_sc)
# normalize the number of counts within each cell to a fixed  number
print('Normalizing counts ...')
sc.pp.normalize_total(ad_sc)
print('Preparing for mapping ...')
tg.pp_adatas(ad_sc, ad_sp, genes=None)

##### Map
print('Mapping single cells onto space ...')
for i in range(9):
    # either using cpu or using GPU, I strongly recommend GPU, either from the server or from Google Collab
    ad_map = tg.map_cells_to_space(
        adata_sc=ad_sc,
        adata_sp=ad_sp,
        device='cuda:0',
        # device='cpu',
        mode='constrained',
        density_prior='uniform',
        num_epochs=500,
        target_count=ad_sp.shape[1],
        lambda_f_reg=1,
        lambda_count=1
    )
    print('Saving map')
    ad_map.write(filename=f'data/ad_maps/STARmap/ad_map{i}.h5ad')