import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
sp_path = "data/20180505_BY3_1kgenes/"

##### Spatial Data (STARmap)
print('Reading in spatial data ... ')
cellxgene_df = pd.read_csv(sp_path + "cell_barcode_count.csv", header=None)
genes = pd.read_csv(sp_path + "genes.csv", header=None)
cellxgene_df.columns = genes[0].values
centroids = pd.read_csv(sp_path + "centroids.tsv", header=None, sep="\t")
n_obs = len(centroids)
print('Making metadata ...')
obs_meta = pd.DataFrame({
    'x':centroids[0].values,
    'y':centroids[1].values
    },
    index=np.arange(n_obs, dtype=int).astype(str),
)
print('Making AnnData object ...')
ad_sp = ad.AnnData(cellxgene_df, obs=obs_meta)
print('Save STARmap h5ad object')
ad_sp.write("data/STARmap.h5ad")

xs = ad_sp.obs.x.values
ys = ad_sp.obs.y.values
plt.axis('off')
plt.scatter(xs, ys, s=3)
plt.title('STARmap Spatial Data')
plt.show()

##### SN Data
print('Reading metadata ...')
metadata = pd.read_csv("data/metadata.csv")
metadata.index = metadata['sample_name']
metadata = metadata.drop(columns=['sample_name'])

print('Reading sn data ...')
cellxexons_df = pd.read_csv("data/cells_exons_VISp.csv")
cellxexons_df = cellxexons_df.set_index('Unnamed: 0').T

metadata = metadata.loc[cellxexons_df.index, :]
print('Making AnnData object ..')
ad_sc = ad.AnnData(cellxexons_df, obs=metadata)
print('Save SMARTSeq2 h5ad object')
ad_sc.write("data/SMARTSeq2.h5ad")
