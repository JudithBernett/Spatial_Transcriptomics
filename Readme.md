# Robustness Tests for Tangram and CellTrek

## Tangram
1. Setup your environment using one of the .yml files
2. Make a data folder and download the Tangram sample data by running:
```
wget https://storage.googleapis.com/tommaso-brain-data/tangram_demo/mop_sn_tutorial.h5ad.gz -O data/mop_sn_tutorial.h5ad.gz
wget https://storage.googleapis.com/tommaso-brain-data/tangram_demo/slideseq_MOp_1217.h5ad.gz -O data/slideseq_MOp_1217.h5ad.gz
wget https://storage.googleapis.com/tommaso-brain-data/tangram_demo/MOp_markers.csv -O data/MOp_markers.csv
gunzip -f data/mop_sn_tutorial.h5ad.gz
gunzip -f data/slideseq_MOp_1217.h5ad.gz
```
3. Make a directory data/ad_maps/ and run obtain_tangram_ad_maps to get all mappings