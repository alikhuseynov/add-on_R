##--------------------------------------------------
### make sure all libs are installed a priori!
## load libs ----
    suppressPackageStartupMessages({
    library(ggplot2)
    library(ggdark)
    library(ggridges)
    library(ggpubr)
    library(SingleCellExperiment)
    library(scater)
    library(Seurat)
    library(dplyr)
    library(magrittr)
    library(png)
    library(cowplot)
    library(parallel)
    library(harmony)
    library(cetcolor)
    library(gridExtra)
    library(plotly)
    library(ggridges)
    library(jpeg)
    library(fields)
    library(spatstat)
    library(patchwork)
    library(Matrix)
    #library(BayesSpace)
    library(progressr)    
    #library(SpotClean)
    #library(SpatialExperiment)  
    #library(STutility)  
    #library(terra)
    #library(Giotto)
    library(BiocParallel)  
})

# set object params & names ----
assay.name <- "Vizgen"
project.name <- "Medulloblastoma - Vizgen"
sample.name <- "my_sample_name"

# function to return metadata of Seurat obj
callmeta <- function(object = NULL) {
    return(object@meta.data) 
}
# eg `callmeta(obj) %>% str`

# eg, load 2 samples ----
dir_use <- c("vizgen_data/xxx/region_0",
             "vizgen_data/xxx/region_1")
# each of path "region_" folders should have:
	# detected_transcripts.csv # molecule coords
	# "Cellpose_" folder which should have these files:
		# cellpose_mosaic_space.parquet
		# cellpose_micron_space.parquet
		# cellpose_cell_by_gene.csv
		# cellpose_cell_metadata.csv

# load optimized `LoadVizgen()` & `ReadVizgen()`
source("read_load_vizgen_v2.R")
# ie https://github.com/alikhuseynov/add-on_R/blob/develop/R/example_read_load_vizgen_v2.R

# load objects using lapply to a list - using `BiocParallel`
# set params
mol.type <- "microns"
coord.space <- "micron"
z.stack <- 3L # z stack to use
fov <- c("region_0", "region_1")

# create Seurat obj
message("Using optimized `LoadVizgen_opt()` function")
obj_vz_list <- 
  lapply(dir_use %>% seq, function(i) {
    message("Creating/Loading object(s): ", "\n", fov[i])
    gc() %>% invisible() # collect garbage
    # renders progress of function
    with_progress(
        LoadVizgen_opt(data.dir = dir_use[i],
                       fov = fov[i], 
                       assay = "Vizgen",
                       use.cellpose.out = TRUE, # if TRUE and ./Cellpose dir exists
                       #filter = "^Blank-", # filter gene names starting with specific pattern
                       metadata = c("volume", "fov"), # add cell volume info
                       mol.type = mol.type, # molecule coords in µm space
                       coord.space = coord.space, 
                       type = c("segmentations", "centroids", "boxes"), # type of cell spatial coord matrices
                       z = z.stack,
                       add.zIndex = TRUE, # add z slice section to a cell
                       update.object = TRUE,
                       use.BiocParallel = TRUE,
                       workers.total = 15, # for `BiocParallel` processing 
                       DTthreads.pct = NULL # percentage of total threads to use for `data.table::fread`
        )
    )
    }
    )

end.time <- Sys.time()
message("Time taken to create/load object = ", 
        round(end.time - start.time, digits = 2), " ", 
        attr(c(end.time - start.time), "units"))

# list of objects
obj_vz_list
