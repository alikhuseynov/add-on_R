
##=======================
## using specific Seurat branch
message("Using Seurat 'feat/imaging' branch")
# remotes::install_github(repo = 'satijalab/seurat', ref = 'feat/imaging') # install 'feat/imaging' branch

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
})

##=======================
# set object params & names
assay.name <- "Vizgen"
project.name <- "Medulloblastoma - Vizgen"
sample.name <- "my_sample_name"

# set root dir ----
dir_working <- "./vizgen_data/analyzed_data/" 
setwd(dir_working)
## this folder should be the output from MERSCOPE Instrument
## eg, it should contain something like that:
#── HuBrainTumor_S1-xxx_VSxx_QC_Cellbound_VAxx_Vxx_XX_x-x-xx
#│   ├── region_0
#│   ├── region_1
#│   ├── region_2

#── HuBrainTumor_S2-xxx_VSxx_QC_Cellbound_VAxx_Vxx_XX_x-x-xx
#│   ├── region_0
#│   ├── region_1
#│   ├── region_2

## then one can select to read/load each region at a time
# set dir to use - specific sample ----
dir_use <- list.files(dir_working, pattern = "_S1|region_1", full.names = TRUE) %>% 
    list.files(pattern = "region_1", full.names = TRUE)
dir_use
setwd(dir_use)
# system(paste0("ls -lth ", dir_use), intern = TRUE) # check dir content.

# function to return metadata of Seurat obj
callmeta <- function(object = NULL) {
    return(object@meta.data) 
}
# eg `callmeta(obj) %>% str`

##=======================
# read dataset - generates data.frame to create an obj
# optimized function - mclapply based..
# load optimized `LoadVizgen()` & `ReadVizgen()`
source("read_load_vizgen_v1.R")
# ie https://github.com/alikhuseynov/add-on_R/blob/develop/R/read_load_vizgen_v1.R

# set params
mol.type = "microns"
coord.space = "micron"
z.stack = 3L # z stack to use

# optimized function - mclapply based..
message("Using optimized `ReadVizgen_opt()` function")
gc() %>% invisible()
start.time <- Sys.time()
message("Reading Vizgen data: ", sample.name)
setwd(dir_use)

# renders progress of function
# [vignette("progressr-intro")](https://progressr.futureverse.org/articles/progressr-intro.html)
handlers("progress")
with_progress(
  data <- ReadVizgen_opt(data.dir = dir_use,
                         use.cellpose.out = TRUE, # if TRUE and ./Cellpose dir exists
                         #transcripts = "./Cellpose/cellpose_cell_by_gene.csv", # count matrix
                         #molecules = "./detected_transcripts.csv", # molecule spatial coord matrices
                         #spatial = "./Cellpose/cellpose_cell_metadata.csv", # cell spatial coord matrices
                         #filter = "^Blank-", # filter gene names starting with specific pattern
                         metadata = c("volume", "fov"), # add cell volume info
                         mol.type = mol.type, # molecule coords in µm space
                         coord.space = coord.space, 
                         type = c("segmentations", "centroids", "boxes"), # type of cell spatial coord matrices
                         z = z.stack, 
                         use.parallel = TRUE, mc.cores.total = 20, # for `parallel` processing 
                         mc.cores.portion = 2, # denominator for `mc.cores.total`, ie `mc.cores.total / mc.cores.portion` for faster extraction of cell boundaries
                         DTthreads.pct = NULL # percentage of total threads to use for `data.table::fread`
  )
)
data %>% str

end.time <- Sys.time()
message("Time taken to read data = ", 
        round(end.time - start.time, digits = 2), " ",
        attr(c(end.time - start.time), "units"))

##=======================
# create Seurat obj
# set params
mol.type = "microns"
coord.space = "micron"
z.stack = 3L # z stack to use

# create Seurat obj
gc() %>% invisible()
message("Using optimized `LoadVizgen_opt()` function")

start.time <- Sys.time()
## name to store FOV
# eg sample name or fusud sample name with `region.0`
use.fused.name <- TRUE # set to FALSE if to use sample name only!
region.name <- 
  stringr::str_split(dir_use, pattern = "/", simplify = TRUE) %>% 
  grep("region", ., value = TRUE) %>% 
  gsub(pattern = "_", replacement = ".", .)

if (use.fused.name) {
  fov <- paste0(tolower(sample.name), ".", region.name)
} else {
  fov <- tolower(sample.name)
}
message("Creating/Loading object(s): ", "\n", 
        fov)

# renders progress of function
with_progress(
  obj <- 
    LoadVizgen_opt(data.dir = dir_use,
                   fov = fov, 
                   assay = assay.name,
                   use.cellpose.out = TRUE, # if TRUE and ./Cellpose dir exists
                   #transcripts = "./Cellpose/cellpose_cell_by_gene.csv", # count matrix
                   #molecules = "./detected_transcripts.csv", # molecule spatial coord matrices
                   #spatial = "./Cellpose/cellpose_cell_metadata.csv", # cell spatial coord matrices
                   #filter = "^Blank-", # filter gene names starting with specific pattern
                   metadata = c("volume", "fov"), # add cell volume info
                   mol.type = mol.type, # molecule coords in µm space
                   coord.space = coord.space, 
                   type = c("segmentations", "centroids", "boxes"), # type of cell spatial coord matrices
                   z = z.stack,
                   Update.object = TRUE,
                   use.parallel = TRUE, mc.cores.total = 24, # for `parallel` processing 
                   mc.cores.portion = 2,
                   DTthreads.pct = NULL # percentage of total threads to use for `data.table::fread`
    )
)

end.time <- Sys.time()
message("Time taken to create/load object = ", 
        round(end.time - start.time, digits = 2), " ", 
        attr(c(end.time - start.time), "units"))
obj

