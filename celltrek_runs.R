#library(devtools)
#install_github("navinlabcode/CellTrek")
library(data.table)
library(CellTrek)
library(Seurat)

# Read in Data
print('Reading ST data ...')
brain_st_cortex <- readRDS("data/brain_st_cortex.rds")
print(dim(brain_st_cortex))
print('Reading SC data ...')
brain_sc <- readRDS("data/brain_sc.rds")
print(dim(brain_sc))
# Subsample
print('Subsampling SC data ...')
brain_sc <-
  brain_sc[, sample(colnames(brain_sc), size = 500, replace = F),]
print(dim(brain_sc))
# Rename the cells/spots with syntactically valid names
print('Renaming cells/spots')
brain_st_cortex <-
  RenameCells(brain_st_cortex, new.names = make.names(Cells(brain_st_cortex)))
brain_sc <-
  RenameCells(brain_sc, new.names = make.names(Cells(brain_sc)))

# Coembed ST and scRNA-seq datasets using traint
brain_traint <- CellTrek::traint(
  st_data = brain_st_cortex,
  sc_data = brain_sc,
  sc_assay = 'RNA',
  cell_names = 'cell_type'
)

# After coembedding, we can chart single cells to their spatial locations.
# Here, we use the non-linear interpolation (intp = T, intp_lin=F) approach to augment the ST spots.
# Run this 9 times

coord_lists <- list()
for (i in c(1:9)) {
  print(i)
  brain_celltrek <- CellTrek::celltrek(
    st_sc_int = brain_traint,
    int_assay = 'traint',
    sc_data = brain_sc,
    sc_assay = 'RNA',
    reduction = 'pca',
    intp = T,
    intp_pnt = 5000,
    intp_lin = F,
    nPCs = 30,
    ntree = 1000,
    dist_thresh = 0.55,
    top_spot = 1,
    spot_n = 5,
    repel_r = 20,
    repel_iter = 20,
    keep_model = T
  )$celltrek
  
  tmp <- data.table(
    id = rownames(brain_celltrek@meta.data),
    coord_x = brain_celltrek@meta.data$coord_x,
    coord_y = brain_celltrek@meta.data$coord_y
  )
  brain_celltrek$cell_type <-
    factor(brain_celltrek$cell_type, levels = sort(unique(brain_celltrek$cell_type)))
  
  tmp <- tmp[, cell_type := brain_celltrek$cell_type[id]]
  coord_lists <- append(coord_lists, list(i = tmp))
}

names(coord_lists) <- c(1:9)
saveRDS(coord_lists, "coord_lists_topspot1.rds")
