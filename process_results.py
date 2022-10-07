import os
import scanpy as sc
import pandas as pd
import random
import numpy as np
from scipy.stats import pearsonr


def generate_ge_train_test_dfs(ad_ge, dataset):
    df = ad_ge.to_df()
    if dataset == 'STARmap':
        markers = ad_ge.uns['training_genes']
    else:
        df_genes = pd.read_csv('data/MOp_markers.csv', index_col=0)
        markers = np.reshape(df_genes.values, (-1,))
        markers = list(markers)
    markers_lc = [marker.lower() for marker in markers]
    df_training = df[df.columns.intersection(markers_lc)]
    df_test = df[df.columns.difference(markers_lc)]
    return df_training, df_test


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
    #Load mapping
    path_maps = f"data/ad_maps/{dataset}/"
    all_objects = os.listdir(path_maps)

    all_pearsons = []
    mapping_df = pd.DataFrame({'x': ad_sp.obs.x.values, 'y': ad_sp.obs.y.values})
    for i in range(len(all_objects)):
        path = all_objects[i]
        print('Reading object')
        ad_map = sc.read_h5ad(path_maps + path)
        random.seed(1234)
        for j in range(50):
            rnd_number = random.randint(0, ad_map.X.shape[0])
            mapping = ad_map.X.T[:, rnd_number]
            mapping_df[f'{i}_{rnd_number}'] = mapping
        for j in range(i + 1, len(all_objects)):
            path2 = all_objects[j]
            ad_map2 = sc.read_h5ad(path_maps + path2)
            x1 = np.array(ad_map.to_df()).flatten()
            x2 = np.array(ad_map2.to_df()).flatten()
            print(f"Calculating Pearson ...")
            r, pval_p = pearsonr(x1, x2)
            print(f"#### Mapping {i} vs. {j} ####")
            print(f'Pearson normal: {r}')
            all_pearsons.append(r)
    mapping_df = pd.wide_to_long(mapping_df, ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9'], i=['x', 'y'], j='cell',
                                 sep='_')
    mapping_df = mapping_df.reset_index()
    mapping_df = pd.melt(mapping_df, id_vars=['x', 'y', 'cell'], var_name='run', value_name='probability')
    mapping_df.to_csv(f"data/mapping_df_{dataset}.csv")

    df = pd.DataFrame({'all_pearsons': all_pearsons})
    df.to_csv(f"data/correlations_mappings_{dataset}.csv")

    # Gene expression data
    path_ge = f"data/gene_exprs/{dataset}/"
    all_objects = os.listdir(path_ge)
    all_pearsons_train = []
    all_pearsons_test = []

    for i in range(len(all_objects)):
        for j in range(i + 1, len(all_objects)):
            print(f"#### GE {i} vs. {j} ####")
            path1 = all_objects[i]
            path2 = all_objects[j]
            print("Reading data ... ")
            ad_ge = sc.read_h5ad(path_ge + path1)
            ad_ge2 = sc.read_h5ad(path_ge + path2)

            print("Splitting predictions ...")
            df_training1, df_test1 = generate_ge_train_test_dfs(ad_ge, dataset)
            df_training2, df_test2 = generate_ge_train_test_dfs(ad_ge2, dataset)

            print("Calculating coefficients for training data ...")
            x1 = np.array(df_training1).flatten()
            x2 = np.array(df_training2).flatten()
            r, pval_p = pearsonr(x1, x2)
            all_pearsons_train.append(r)

            print("Calculating coefficients for test data ...")
            x1 = np.array(df_test1).flatten()
            x2 = np.array(df_test2).flatten()
            r, pval_p = pearsonr(x1, x2)
            all_pearsons_test.append(r)
    df = pd.DataFrame({
        'all_pearsons_train': all_pearsons_train,
        'all_pearsons_test': all_pearsons_test
    })

    df.to_csv(f"data/correlations_expr_{dataset}.csv")