CellTrek Notebook
================

``` r
library(data.table)
library(ggplot2)
library(CellTrek)
```

    ## Warning: vorhergehender Import 'data.table::last' durch 'dplyr::last' während
    ## des Ladens von 'CellTrek' ersetzt

    ## Warning: vorhergehender Import 'data.table::first' durch 'dplyr::first' während
    ## des Ladens von 'CellTrek' ersetzt

    ## Warning: vorhergehender Import 'MASS::select' durch 'dplyr::select' während des
    ## Ladens von 'CellTrek' ersetzt

    ## Warning: vorhergehender Import 'data.table::between' durch 'dplyr::between'
    ## während des Ladens von 'CellTrek' ersetzt

    ## Warning: vorhergehender Import 'dplyr::union' durch 'igraph::union' während des
    ## Ladens von 'CellTrek' ersetzt

    ## Warning: vorhergehender Import 'dplyr::as_data_frame' durch
    ## 'igraph::as_data_frame' während des Ladens von 'CellTrek' ersetzt

    ## Warning: vorhergehender Import 'dplyr::groups' durch 'igraph::groups' während
    ## des Ladens von 'CellTrek' ersetzt

    ## Warning: vorhergehender Import 'igraph::groups' durch 'plotly::groups' während
    ## des Ladens von 'CellTrek' ersetzt

    ## Warning: vorhergehender Import 'magrittr::set_names' durch 'purrr::set_names'
    ## während des Ladens von 'CellTrek' ersetzt

    ## Warning: vorhergehender Import 'data.table::transpose' durch 'purrr::transpose'
    ## während des Ladens von 'CellTrek' ersetzt

    ## Warning: vorhergehender Import 'igraph::compose' durch 'purrr::compose' während
    ## des Ladens von 'CellTrek' ersetzt

    ## Warning: vorhergehender Import 'igraph::simplify' durch 'purrr::simplify'
    ## während des Ladens von 'CellTrek' ersetzt

    ## Warning: vorhergehender Import 'purrr::partial' durch 'randomForestSRC::partial'
    ## während des Ladens von 'CellTrek' ersetzt

    ## Warning: vorhergehender Import 'data.table::melt' durch 'reshape2::melt' während
    ## des Ladens von 'CellTrek' ersetzt

    ## Warning: vorhergehender Import 'data.table::dcast' durch 'reshape2::dcast'
    ## während des Ladens von 'CellTrek' ersetzt

    ## Warning: vorhergehender Import 'purrr::discard' durch 'scales::discard' während
    ## des Ladens von 'CellTrek' ersetzt

    ## Warning: vorhergehender Import 'igraph::as_data_frame' durch
    ## 'tibble::as_data_frame' während des Ladens von 'CellTrek' ersetzt

    ## Warning: vorhergehender Import 'magrittr::extract' durch 'tidyr::extract'
    ## während des Ladens von 'CellTrek' ersetzt

    ## Warning: vorhergehender Import 'igraph::crossing' durch 'tidyr::crossing'
    ## während des Ladens von 'CellTrek' ersetzt

``` r
library(dplyr)
```

    ## 
    ## Attache Paket: 'dplyr'

    ## Die folgenden Objekte sind maskiert von 'package:data.table':
    ## 
    ##     between, first, last

    ## Die folgenden Objekte sind maskiert von 'package:stats':
    ## 
    ##     filter, lag

    ## Die folgenden Objekte sind maskiert von 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(Seurat)
```

    ## Attaching SeuratObject

    ## Attaching sp

``` r
library(RColorBrewer)
```

``` r
brain_st_cortex <- readRDS("data/brain_st_cortex.rds")
print(dim(brain_st_cortex))
```

    ## [1] 31053  1075

``` r
brain_sc <- readRDS("data/brain_sc.rds")
print(dim(brain_sc))
```

    ## [1] 34617  4785

``` r
## Rename the cells/spots with syntactically valid names
brain_st_cortex <- RenameCells(brain_st_cortex, new.names=make.names(Cells(brain_st_cortex)))
brain_sc <- RenameCells(brain_sc, new.names=make.names(Cells(brain_sc)))
```

``` r
SpatialDimPlot(brain_st_cortex)
```

![](celltrek_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
DimPlot(brain_sc, label = T, label.size = 4.5)
```

![](celltrek_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

``` r
coord_lists <- readRDS("data/coord_lists_topspot1.rds")
all_coords <- rbindlist(coord_lists, idcol="run")
all_coords_wide <- dcast(all_coords, id ~ run, drop=FALSE, value.var=c("coord_x", "coord_y"))

rowDistsMean <- function(tmprow) {
  tmpx <- tmprow[2:11]
  tmpy <- tmprow[12:21]
  return(mean(dist(t(rbind(tmpx, tmpy))), na.rm = TRUE))
}
rowDistsMedian <- function(tmprow) {
  tmpx <- tmprow[2:11]
  tmpy <- tmprow[12:21]
  return(median(dist(t(rbind(tmpx, tmpy))), na.rm = TRUE))
}
dists <- data.table(id=all_coords_wide$id, mean=apply(all_coords_wide, 1, rowDistsMean))
dists$median <- apply(all_coords_wide, 1, rowDistsMedian)
dists <- melt(dists, measure.var = c("mean", "median"), variable.name = "euclidean distance", value.name="value")
ggplot(dists, aes(x = `euclidean distance`, y = value))+
  geom_boxplot()+
  theme_minimal()
```

    ## Warning: Removed 862 rows containing non-finite values (stat_boxplot).

![](celltrek_files/figure-gfm/unnamed-chunk-4-1.png)<!-- --> After cell
charting, we can interactively visualize the CellTrek result using the
code from celltrek_vis

``` r
plot_celltrek <- function(coords){
  coords <- coords[!is.na(cell_type)]
  img_temp <- brain_st_cortex@images$anterior1@image
  scale_factor <- brain_st_cortex@images$anterior1@scale.factors$lowres
  pnt_colors <- colorRampPalette(brewer.pal(9, 
                "Set1"))(length(levels(coords$cell_type)))
  plotly::plot_ly(d = coords, 
                  x = ~coord_y * scale_factor, 
                  y = ~ dim(img_temp)[1] - coord_x * scale_factor, customdata = ~id, 
                  color = ~cell_type,
                  colors = pnt_colors,
                  type = "scatter", 
                  mode = "markers", 
                  marker = list(line = list(color = "rgb(1,1,1)", 
                  width = 0.5), size = 8, opacity = 0.8)) %>% 
    plotly::layout(xaxis = list(range = c(0, dim(img_temp)[2]), showgrid = FALSE, showline = FALSE), 
                   yaxis = list(range = c(0, dim(img_temp)[1]), showgrid = FALSE, showline = FALSE), 
                   images = list(source = plotly::raster2uri(as.raster(img_temp)),
                                 x = 0, y = 0, sizex = dim(img_temp)[2], 
                                 sizey = dim(img_temp)[1], xref = "x", yref = "y", 
                                 xanchor = "left", yanchor = "bottom", 
                                 layer = "below", sizing = "stretch")
                   )
}
```

``` r
p1 <- plot_celltrek(as.data.table(coord_lists[[1]]))
p2 <- plot_celltrek(as.data.table(coord_lists[[2]]))
p3 <- plot_celltrek(as.data.table(coord_lists[[3]]))
p4 <- plot_celltrek(as.data.table(coord_lists[[4]]))
p5 <- plot_celltrek(as.data.table(coord_lists[[5]]))
p6 <- plot_celltrek(as.data.table(coord_lists[[6]]))
p7 <- plot_celltrek(as.data.table(coord_lists[[7]]))
p8 <- plot_celltrek(as.data.table(coord_lists[[8]]))
p9 <- plot_celltrek(as.data.table(coord_lists[[9]]))
```

``` r
plotly::subplot(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrows=3)
```

![](celltrek_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
top_20_var_cells <- data.table(dists[`euclidean distance` == 'median'][order(value, decreasing = T)][1:20, c("id", "value")])
top_20_var_cells$top <- factor(paste("top", c(1:20), sep="_"), levels = paste("top", c(1:20), sep="_"))

plot_celltrek_var <- function(coord) {
  tmp <- merge(as.data.table(coord)[id %in% top_20_var_cells$id], top_20_var_cells)
  img_temp <- brain_st_cortex@images$anterior1@image
  scale_factor <- brain_st_cortex@images$anterior1@scale.factors$lowres
  pnt_colors <- colorRampPalette(brewer.pal(12, "Set1"))(length(levels(tmp$top)))
  plotly::plot_ly(d = tmp, 
                  x = ~coord_y * scale_factor, 
                  y = ~ dim(img_temp)[1] - coord_x * scale_factor, 
                  customdata = ~id, 
                  color = ~top,
                  text = ~top,
                  textposition = 'bottom center',
                  textfont = list(color = '#fff000', size = 8),
                  colors = pnt_colors,
                  type = "scatter", 
                  mode = "markers+text", 
                  marker = list(line = list(color = "rgb(1,1,1)", 
                  width = 0.5
                  ), size = 8, opacity = 0.8)) %>% 
    plotly::layout(xaxis = list(range = c(0, dim(img_temp)[2]), showgrid = FALSE, showline = FALSE), 
                   yaxis = list(range = c(0, dim(img_temp)[1]), showgrid = FALSE, showline = FALSE), 
                   images = list(source = plotly::raster2uri(as.raster(img_temp)),
                                 x = 0, y = 0, sizex = dim(img_temp)[2], 
                                 sizey = dim(img_temp)[1], xref = "x", yref = "y", 
                                 xanchor = "left", yanchor = "bottom", 
                                 layer = "below", sizing = "stretch")
                   )
}
```

``` r
p1 <- plot_celltrek_var(as.data.table(coord_lists[[1]]))
```

    ## Warning in brewer.pal(12, "Set1"): n too large, allowed maximum for palette Set1 is 9
    ## Returning the palette you asked for with that many colors

``` r
p2 <- plot_celltrek_var(as.data.table(coord_lists[[2]]))
```

    ## Warning in brewer.pal(12, "Set1"): n too large, allowed maximum for palette Set1 is 9
    ## Returning the palette you asked for with that many colors

``` r
p3 <- plot_celltrek_var(as.data.table(coord_lists[[3]]))
```

    ## Warning in brewer.pal(12, "Set1"): n too large, allowed maximum for palette Set1 is 9
    ## Returning the palette you asked for with that many colors

``` r
p4 <- plot_celltrek_var(as.data.table(coord_lists[[4]]))
```

    ## Warning in brewer.pal(12, "Set1"): n too large, allowed maximum for palette Set1 is 9
    ## Returning the palette you asked for with that many colors

``` r
p5 <- plot_celltrek_var(as.data.table(coord_lists[[5]]))
```

    ## Warning in brewer.pal(12, "Set1"): n too large, allowed maximum for palette Set1 is 9
    ## Returning the palette you asked for with that many colors

``` r
p6 <- plot_celltrek_var(as.data.table(coord_lists[[6]]))
```

    ## Warning in brewer.pal(12, "Set1"): n too large, allowed maximum for palette Set1 is 9
    ## Returning the palette you asked for with that many colors

``` r
p7 <- plot_celltrek_var(as.data.table(coord_lists[[7]]))
```

    ## Warning in brewer.pal(12, "Set1"): n too large, allowed maximum for palette Set1 is 9
    ## Returning the palette you asked for with that many colors

``` r
p8 <- plot_celltrek_var(as.data.table(coord_lists[[8]]))
```

    ## Warning in brewer.pal(12, "Set1"): n too large, allowed maximum for palette Set1 is 9
    ## Returning the palette you asked for with that many colors

``` r
p9 <- plot_celltrek_var(as.data.table(coord_lists[[9]]))
```

    ## Warning in brewer.pal(12, "Set1"): n too large, allowed maximum for palette Set1 is 9
    ## Returning the palette you asked for with that many colors

``` r
plotly::subplot(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrows=3)
```

![](celltrek_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->
