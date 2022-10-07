import numpy as np
import pandas as pd
import scanpy as sc
import tangram as tg

# we follow the Tangram Tutorial without squidpy

##### Load and visualize the spatial data:  Slide-seq data, 9852 spatial voxels, 24518 genes
ad_sp = sc.read_h5ad('data/slideseq_MOp_1217.h5ad')
print(ad_sp)

##### Load single cell data, MOp 10Xv3 dataset, approximately 26k profiled cells with 28k genes
ad_sc = sc.read_h5ad('data/mop_sn_tutorial.h5ad')
print(ad_sc)
# normalize the number of counts within each cell to a fixed  number
sc.pp.normalize_total(ad_sc)

##### Prepare to map
df_genes = pd.read_csv('data/MOp_markers.csv', index_col=0)
markers = np.reshape(df_genes.values, (-1, ))
markers = list(markers)
print(len(markers))

tg.pp_adatas(ad_sc, ad_sp, genes=markers)

##### Map
# either using cpu or using GPU, I strongly recommend GPU, either from the server or from Google Collab
ad_map = tg.map_cells_to_space(
    adata_sc=ad_sc,
    adata_sp=ad_sp,
    #device='cpu',
    device='cuda:0',
)
ad_map.write(filename='data/ad_maps/ad_map1.h5ad')