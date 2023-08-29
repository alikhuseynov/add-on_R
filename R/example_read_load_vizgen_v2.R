##--------------------------------------------------
### make sure all libs are installed a priori!
# load libs ----
suppressPackageStartupMessages({
  library(ggplot2)
  library(Seurat)
  library(dplyr)
  library(magrittr)
  library(BiocParallel)
  library(progressr)
  library(sf)
})

# helper function to return Seurat object metadata
callmeta <- function (object = NULL) { 
  return(object@meta.data) 
  }

# set directory to read from
dir_use <- "./merfish_seurat_debug/" # or something like "./region_0/"

##
start.time <- Sys.time()
obj <-
  LoadVizgen(data.dir = dir_use,  
             fov = "merfish.test", 
             assay = "Vizgen",
             metadata = c("volume", "fov"), # add cell volume info
             type = c("segmentations", "centroids"), # type of cell spatial coord matrices
             z = 3L,
             add.zIndex = TRUE, # add z slice section to a cell
             update.object = TRUE,
             use.BiocParallel = TRUE,
             workers.MulticoreParam = 14, # for `BiocParallel` processing
             min.area = 5, # minimal polygon area to use as a threshold for filtering segmentaion geometries
             add.molecules = TRUE, # if to add "molecules" coordinates to FOV of the object
             verbose = TRUE
             )
end.time <- Sys.time()
message("Time taken to load object = ", 
        round(end.time - start.time, digits = 2), " ", 
        attr(c(end.time - start.time), "units"))
obj
obj %>% callmeta %>% str
