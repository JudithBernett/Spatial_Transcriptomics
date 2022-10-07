# Robustness Tests for Tangram and CellTrek

## Tangram
1. Setup your environment using one of the .yml files
2. Make a data folder and download the Tangram sample data by running:
```
# MOp 10Xv3 dataset
wget https://storage.googleapis.com/tommaso-brain-data/tangram_demo/mop_sn_tutorial.h5ad.gz -O data/mop_sn_tutorial.h5ad.gz
# SlideSeq data
wget https://storage.googleapis.com/tommaso-brain-data/tangram_demo/slideseq_MOp_1217.h5ad.gz -O data/slideseq_MOp_1217.h5ad.gz
wget https://storage.googleapis.com/tommaso-brain-data/tangram_demo/MOp_markers.csv -O data/MOp_markers.csv
gunzip -f data/mop_sn_tutorial.h5ad.gz
gunzip -f data/slideseq_MOp_1217.h5ad.gz
# Smart-Seq2_VISp_snRNAseq
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115746/suppl/GSE115746_cells_exon_counts.csv.gz -O data/cells_exons_VISp.csv.gz
gunzip data/cells_exons_VISp.csv.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115746/suppl/GSE115746_complete_metadata_28706-cells.csv.gz -O data/metadata.csv.gz
gunzip data/metadata.csv.gz
```
3. For STARMAP download https://zenodo.org/record/3967291#.Yz7aAexBxJU and transfer the Spatial/Starmap/visual_1020/20180505_BY3_1kgenes/ folder to data
4. Convert the Smart-Seq2 and STARMAP datasets to h5ad files with ```convert_to_h5ad.py```
5. Run ```obtain_tangram_ad_maps``` to get all mappings
6. Run ```obtain_tangram_expressions``` to get all gene expression predictions