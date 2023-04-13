#' @import dplyr
#' @importFrom magrittr %>% %<>%
#' @importFrom progressr progressor
NULL

#' Optimized functions solution to \code{subset()}:
#' New arguments:
#' @param use.BiocParallel If to use \code{BiocParallel::bplapply()}, default is \code{TRUE}, if \code{FALSE}, uses \code{future} library
#' @param workers.total Number of cores to use for \code{BiocParallel::bplapply()}
#' @param DTthreads.pct Set percentage eg \code{50} of total threads to use for \code{data.table::fread}, if set to \code{NULL} will use default setting as in \code{data.table::getDTthreads(verbose = T)}
#' @param mol.type.use In which space to use molecules coords: "microns" or "pixels". Default to "microns". This arg is for creating an object \code{LoadVizgen_opt()}
#' @param coord.space = "micron" Which coordinate space to use: "microns" or "mosaic" (pixel space). Default to "microns"
#' @param use.cellpose.out If TRUE, and ./Cellpose folder exists, will load results from current MERSCOPE Instrument output. Default to TRUE. Set to FALSE if to use previous outputs (ie. non-Cellpose).
#' @param add.zIndex If to add \code{z} slice index to a cell
#' @param update.object If to update final object, default to TRUE.
#' @param ... Arguments passed to \code{ReadVizgen_opt()}

## TODO:
    # - add modifications to cloned gihub repo https://github.com/satijalab/seurat/tree/develop new branch `feat/vizgen`
    # - add an option to read image data, store image as a raster similar to Visium `Read10X_Image`


## ==============================================================
## Reading data - output is a list of dataframes
## ==============================================================

ReadVizgen_opt <-
  function (data.dir, transcripts = NULL, spatial = NULL, molecules = NULL, type = "segmentations", 
            use.cellpose.out = TRUE, mol.type = "microns", coord.space = "micron", 
            metadata = NULL, filter = NA_character_ , z = 3L, use.BiocParallel = TRUE, 
            workers.total = 12, DTthreads.pct = data.table::getDTthreads(verbose = F)
  )
    
  {   
    if (!requireNamespace("data.table", quietly = TRUE)) {
      stop("Please install 'data.table' for this function", call. = FALSE)
    }
      
    if (!requireNamespace("BiocParallel", quietly = TRUE)) {
      stop("Please install 'BiocParallel' for parallelization", call. = FALSE)
    }
    
    if (!requireNamespace("sfarrow", quietly = TRUE)) {
      stop("Please install 'sfarrow' for reading cell boundaries from `.parquet` files ",
		   call. = FALSE)
    }

    # setting workers to use for parallel computing - `parallel`
    if (use.BiocParallel) {
      message("Using parallelization with: `BiocParallel`")
      if (is.null(workers.total)) { 
        workers.total <- quantile(BiocParallel::multicoreWorkers() %>% seq)[4] %>% ceiling
        message(workers.total, " of total workers available will be used")
      } else { message("Setting total workers to: ", workers.total) }   
    } else { message("Using parallelization with: `future`", "\n",
                     "NOTE: set workers for parallelization, eg: `future::plan('multisession', workers = 10)`")
           }
      
    # setting cores to use for parallel computing - `data.table`
    if (!is.null(DTthreads.pct)) {
      message("Using parallelization with: `data.table`", "\n",
              "..for `data.table::fread`")
      data.table::setDTthreads(threads = 0) # all cores
      DTthreads <- data.table::getDTthreads() # max cores
      DTthreads <- c((DTthreads * DTthreads.pct) / 100) %>% ceiling # percentage from total threads
      message("Setting DTthreads to: ", DTthreads, " (", paste0(DTthreads.pct, "%"), ")")
      data.table::setDTthreads(threads = DTthreads) # set
    }
    
    gc() %>% invisible() # collect garbage
    
    hdf5 <- requireNamespace("hdf5r", quietly = TRUE)
    type <- match.arg(arg = type, choices = c("segmentations",
                                              "centroids", "boxes"), several.ok = TRUE)
    mol.type <- match.arg(arg = mol.type, choices = c("pixels",
                                                      "microns"), several.ok = TRUE)
    if (!is.null(x = metadata)) {
      metadata <- match.arg(arg = metadata, 
                            choices = c("volume", "fov"), several.ok = TRUE)
    }
    if (!z %in% seq.int(from = 0L, to = 6L)) {
      stop("The z-index must be in the range [0, 6]")
    }
    use.dir <- all(vapply(X = c(transcripts, spatial, molecules),
                          FUN = function(x) {
                            return(is.null(x = x) || is.na(x = x))
                          }, FUN.VALUE = logical(length = 1L)))
    if (use.dir && !dir.exists(paths = data.dir)) {
      stop("Cannot find Vizgen directory ", data.dir)
    }

    # use Cellpose output
    if (use.cellpose.out && grep("Cellpose", list.files(data.dir)) %>% any) {
      message("Cellpose output will be used..")
      files <- c(transcripts = transcripts %||% paste0(data.dir, 
           paste0("/", list.files(data.dir, pattern = "Cellpose"), 
                  "/cellpose_cell_by_gene.csv")),
                 
                 spatial = spatial %||% paste0(data.dir, 
           paste0("/", list.files(data.dir, pattern = "Cellpose"), 
                  "/cellpose_cell_metadata.csv")),
                 
                 molecules = molecules %||% "detected_transcripts[_a-zA-Z0-9]*.csv")
    } else { 
      files <- c(transcripts = transcripts %||% "cell_by_gene[_a-zA-Z0-9]*.csv",
                 spatial = spatial %||% "cell_metadata[_a-zA-Z0-9]*.csv",
                 molecules = molecules %||% "detected_transcripts[_a-zA-Z0-9]*.csv")
    }
    
    files[is.na(x = files)] <- NA_character_
    h5dir <- file.path(ifelse(test = dirname(path = files["spatial"]) ==
                                ".", yes = data.dir, no = dirname(path = files["spatial"])),
                       "cell_boundaries")
    zidx <- paste0("zIndex_", z)
    files <- vapply(X = files, FUN = function(x) {
      x <- as.character(x = x)
      if (isTRUE(x = dirname(path = x) == ".")) {
          fnames <- list.files(path = data.dir, pattern = x,
                               recursive = FALSE, full.names = TRUE)
          return(sort(x = fnames, decreasing = TRUE)[1L])
      }
      else {
        return(x)
      }
    }, FUN.VALUE = character(length = 1L), USE.NAMES = TRUE)
    files[!file.exists(files)] <- NA_character_
    if (all(is.na(x = files))) {
      stop("Cannot find Vizgen input files in ", data.dir)
    }
    
    if (!is.na(x = files[["spatial"]])) {
      pprecoord <- progressor()
      pprecoord(message = "Preloading cell spatial coordinates",
                class = "sticky", amount = 0)
      sp <- data.table::fread(file = files[["spatial"]], sep = ",",
                              data.table = FALSE, verbose = FALSE)
      pprecoord(type = "finish")
      rownames(x = sp) <- as.character(x = sp[, 1])
      sp <- sp[, -1, drop = FALSE]
      if ("segmentations" %in% type) {
        poly <- if (isFALSE(x = hdf5) && !use.cellpose.out) {
          warning("Cannot find hdf5r; unable to load segmentation vertices",
                  immediate. = TRUE)
          FALSE
        }
        else if (!dir.exists(paths = h5dir) && !use.cellpose.out) {
          warning("Cannot find cell boundary H5 files",
                  immediate. = TRUE)
          FALSE
        } 
        else if (use.cellpose.out) { # added for Cellpose output
          files2scan <-
            list.files(data.dir, pattern = ".parquet$", 
                       all.files = TRUE, 
                       full.names = TRUE, 
                       recursive = TRUE)
          if (length(files2scan)) { 
            message("Cell segmentations are found in `.parquet` file(s)", "\n", 
                    "..using ", coord.space, " space coordinates")
          }
          TRUE
        }
        else {
          TRUE
        }
        if (isFALSE(x = poly)) {
          type <- setdiff(x = type, y = "segmentations")
        }
      }
      spatials <- rep_len(x = files[["spatial"]], length.out = length(x = type))
      names(x = spatials) <- type
      files <- c(files, spatials)
      files <- files[setdiff(x = names(x = files), y = "spatial")]
    }
    else if (!is.null(x = metadata)) {
      warning("metadata can only be loaded when spatial coordinates are loaded",
              immediate. = TRUE)
      metadata <- NULL
    }
    
    if (!is.na(x = files[["molecules"]])) {
      ppremol <- progressor()
      ppremol(message = "Preloading molecule coordinates",
              class = "sticky", amount = 0)
      mx <- data.table::fread(file = files[["molecules"]],
                              sep = ",", verbose = FALSE)
      mx <- mx[mx$global_z == z, , drop = FALSE]
      if (!is.na(x = filter)) {
        ppremol(message = paste("Filtering molecules with pattern",
                                filter), class = "sticky", amount = 0)
        mx <- mx[!grepl(pattern = filter, x = mx$gene), ,
                 drop = FALSE]
      }
      ppremol(type = "finish")
      mols <- rep_len(x = files[["molecules"]], length.out = length(x = mol.type))
      names(x = mols) <- mol.type
      files <- c(files, mols)
      files <- files[setdiff(x = names(x = files), y = "molecules")]
    }
    files <- files[!is.na(x = files)]
    outs <- vector(mode = "list", length = length(x = files))
    names(x = outs) <- names(x = files)
    if (!is.null(metadata)) {
      outs <- c(outs, list(metadata = NULL))
    }
    
    gc() %>% invisible() # collect garbage
    
    for (otype in names(x = outs)) {
      outs[[otype]] <- 
        switch(EXPR = otype, 
               transcripts = {
                 ptx <- progressor()
                 ptx(message = "Reading counts matrix", class = "sticky",
                     amount = 0)
                 tx <- data.table::fread(file = files[[otype]], sep = ",",
                                         data.table = FALSE, verbose = FALSE)
                 rownames(x = tx) <- as.character(x = tx[, 1])
                 tx <- t(x = as.matrix(x = tx[, -1, drop = FALSE]))
                 if (!is.na(x = filter)) {
                   ptx(message = paste("Filtering genes with pattern",
                                       filter), class = "sticky", amount = 0)
                   tx <- tx[!grepl(pattern = filter, x = rownames(x = tx)),
                            , drop = FALSE]
                 }
                 ratio <- getOption(x = "Seurat.input.sparse_ratio",
                                    default = 0.4)
                 if ((sum(tx == 0)/length(x = tx)) > ratio) {
                   ptx(message = "Converting counts to sparse matrix",
                       class = "sticky", amount = 0)
                   tx <- as.sparse(x = tx)
                 }
                 ptx(type = "finish")
                 tx
                 
                 }, centroids = {
                   pcents <- progressor()
                   pcents(message = "Creating centroid coordinates",
                          class = "sticky", amount = 0)
                   pcents(type = "finish")
                   data.frame(x = sp$center_x, y = sp$center_y, cell = rownames(x = sp),
                              stringsAsFactors = FALSE)
                   
                   # use Cellpose segmentations
                   }, segmentations = {
                     if (use.cellpose.out) { 
                       files2scan <-
                         list.files(data.dir, pattern = ".parquet$", 
                                    all.files = TRUE, 
                                    full.names = TRUE, 
                                    recursive = TRUE)
                       if (length(files2scan)) {
                         file2read <- files2scan %>%
                           grep(coord.space, ., value = TRUE) %>% 
                           grep(".parquet$", ., value = TRUE)
                       }
                       
                       # read .parquet file..
                       parq <- sfarrow::st_read_parquet(file2read)
                       gc() %>% invisible() # free up memory
                         
                       # get all cell segmentations
                       segs <- filter(parq, ZIndex == z) %>% pull(Geometry)
                       # check if any cells have > 1 segmentation boundary
                       test.segs <- 
                         lapply(segs %>% seq, function(i) segs[[i]][[1]] %>% length)
                       if (which(unlist(test.segs) > 1) %>% any) {
                         segs.art.index <- which(unlist(test.segs) > 1)
                         message(segs.art.index %>% length, 
                                 " Cells have > 1 segmentaion boundary artifacts", "\n",
                                 "..removing artifacts", "\n",
                                 "..keeping cell boundaries with maximum coords")
                         # usually artifacts have small boundaries/coords
                         # find cell boundaries with maximum coords
                         for (i in segs.art.index %>% seq) { 
                           dims <- lapply(segs[[segs.art.index[i]]][[1]] %>% seq(), 
                                          function(d) { dim(segs[[segs.art.index[i]]][[1]][[d]])[1] } )
                           # get & keep boundaries with maximum coords
                           maxs.segs <- which(unlist(dims) == max(unlist(dims)))  
                           segs[[segs.art.index[i]]][[1]] <- segs[[segs.art.index[i]]][[1]][maxs.segs]
                         }
                       } else { message("All cells have 1 segmentaion boundary (no artifacts)") }
                       # add cell names 
                       names(segs) <- filter(parq, ZIndex == z) %>% pull(EntityID) %>% as.character       
                       
                       # TODO: (optionally) resample & make cell boundaries equidistant!  
                       if (use.BiocParallel) {
                           # extract cell boundaries per cells (faster)
                           gc() %>% invisible() # free up memory
                           message("Using fast extraction of cell boundaries..")  
                           segs_list <-
                           BiocParallel::bplapply(segs %>% seq,
                                    function(i) {
                                        segs[[i]][[1]] %>% 
                                            data.table::as.data.table(x = .) %>%
                                            as.data.frame %>%
                                            mutate(cell = names(segs)[i]) },
                                    BPPARAM = BiocParallel::MulticoreParam(workers.total, tasks = 100L, 
                                                                           force.GC = FALSE, 
                                                                           progressbar = TRUE)
                                   )
                       } else {
                       # extract cell boundaries per cells
                       message("Extracting of cell boundaries..")    
                       segs_list <-
                           lapply(segs %>% seq,
                                  function(i) {
                                      segs[[i]][[1]] %>%
                                          data.table::as.data.table(x = .) %>%
                                          as.data.frame %>% 
                                          mutate(cell = names(segs)[i])
                                  }
                                 )
                       }
                                
                       #segs_list %>% length
                       # df of all cell segmentations
                       gc() %>% invisible() # free up memory
                       #segs <- do.call("rbind", segs_list)
                       segs <- data.table::rbindlist(segs_list)
                       names(segs)[1:2] <- c("x", "y") 
                       message("All cell boundaries (with no artifacts) are loaded..")        
                       segs
                     } else { # use non-Cellpose segmentations
                       ppoly <- progressor(steps = length(x = unique(x = sp$fov)))
                       ppoly(message = "Creating polygon coordinates", class = "sticky",
                             amount = 0)
                       # use parallel or future
                       if (use.BiocParallel) { 
                         pg <- BiocParallel::bplapply(X = unique(x = sp$fov), FUN = function(f, ...) {
                           fname <- file.path(h5dir, paste0("feature_data_", f, ".hdf5"))
                           if (!file.exists(fname)) {
                             warning("Cannot find HDF5 file for field of view ",
                                     f, immediate. = TRUE)
                             return(NULL)
                           }
                           # reading hdf5 files
                           hfile <- hdf5r::H5File$new(filename = fname, 
                                                      mode = "r")
                           on.exit(expr = hfile$close_all())
                           cells <- rownames(x = subset(x = sp, subset = fov == f))     
                           # creating df for cell boundaries     
                           df <- BiocParallel::bplapply(X = cells, FUN = function(x) {
                             return(tryCatch(expr = {
                               cc <- hfile[["featuredata"]][[x]][[zidx]][["p_0"]][["coordinates"]]$read()
                               cc <- as.data.frame(x = t(x = cc))
                               colnames(x = cc) <- c("x", "y")
                               cc$cell <- x
                               cc
                             }, error = function(...) {
                               return(NULL)
                             }))
                           }, BPPARAM = BiocParallel::MulticoreParam(workers.total, tasks = 100L, 
                                                                     force.GC = FALSE, 
                                                                     progressbar = TRUE))   
                             gc() %>% invisible()
                           ppoly()
                           return(do.call(what = "rbind", args = df))
                         }, BPPARAM = BiocParallel::MulticoreParam(workers.total, tasks = 100L, 
                                                                   force.GC = FALSE, 
                                                                   progressbar = TRUE))
                           gc() %>% invisible()
                       } else { 
                         pg <- 
                           future.apply::future_lapply(X = unique(x = sp$fov), 
                                                       FUN = function(f, ...) {
                                                         fname <- file.path(h5dir, paste0("feature_data_",
                                                                                          f, ".hdf5"))
                                                         if (!file.exists(fname)) {
                                                           warning("Cannot find HDF5 file for field of view ",
                                                                   f, immediate. = TRUE)
                                                           return(NULL)
                                                         }
                                                         hfile <- hdf5r::H5File$new(filename = fname,
                                                                                    mode = "r")
                                                         on.exit(expr = hfile$close_all())
                                                         cells <- rownames(x = subset(x = sp, subset = fov == f))
                                                         df <- lapply(X = cells, FUN = function(x) {
                                                           return(tryCatch(expr = {
                                                             cc <- hfile[["featuredata"]][[x]][[zidx]][["p_0"]][["coordinates"]]$read()
                                                             cc <- as.data.frame(x = t(x = cc))
                                                             colnames(x = cc) <- c("x", "y")
                                                             cc$cell <- x
                                                             cc
                                                           }, error = function(...) {
                                                             return(NULL)
                                                           }))
                                                         })
                                                         ppoly()
                                                         return(do.call(what = "rbind", args = df))
                                                         }
                                                       )
                         }
                       gc() %>% invisible() # collect garbage
                       ppoly(type = "finish")
                       # cell polygons
                       pg <- do.call(what = "rbind", args = pg)
                       npg <- length(x = unique(x = pg$cell))
                       if (npg < nrow(x = sp)) {
                         warning(nrow(x = sp) - npg, " cells missing polygon information",
                                 immediate. = TRUE)
                         }
                       pg # final segmentaions out..
                       }
                     }, boxes = {
                       pbox <- progressor(steps = nrow(x = sp))
                       pbox(message = "Creating box coordinates", class = "sticky",
                            amount = 0)
                       # use parallel or future
                       if (use.BiocParallel) {
                         message("Creating box coordinates..")
                         bx <- BiocParallel::bplapply(X = rownames(x = sp), FUN = function(cell) {
                           row <- sp[cell, ]
                           df <- expand.grid(x = c(row$min_x, row$max_x),
                                             y = c(row$min_y, row$max_y), cell = cell, KEEP.OUT.ATTRS = FALSE,
                                             stringsAsFactors = FALSE)
                           df <- df[c(1, 3, 4, 2), , drop = FALSE]
                           pbox()
                           return(df)
                           }, BPPARAM = BiocParallel::MulticoreParam(workers.total, tasks = 100L, 
                                                                     force.GC = FALSE, progressbar = TRUE)
                           )
                           gc() %>% invisible()
                         } else {
                           bx <- future.apply::future_lapply(X = rownames(x = sp), FUN = function(cell) {
                             row <- sp[cell, ]
                             df <- expand.grid(x = c(row$min_x, row$max_x),
                                               y = c(row$min_y, row$max_y), cell = cell, KEEP.OUT.ATTRS = FALSE,
                                               stringsAsFactors = FALSE)
                             df <- df[c(1, 3, 4, 2), , drop = FALSE]
                             pbox()
                             return(df)
                           })
                           }
                       gc() %>% invisible() # collect garbage
                       pbox(type = "finish")
                       do.call(what = "rbind", args = bx)
                       }, metadata = {
                         pmeta <- progressor()
                         pmeta(message = "Loading metadata", class = "sticky",
                               amount = 0)
                         pmeta(type = "finish")
                         sp[, metadata, drop = FALSE]
                         }, pixels = {
                           ppixels <- progressor()
                           ppixels(message = "Creating pixel-level molecule coordinates",
                                   class = "sticky", amount = 0)
                           df <- data.frame(x = mx$x, y = mx$y, gene = mx$gene,
                                            stringsAsFactors = FALSE)
                           ppixels(type = "finish")
                           df
                           }, microns = {
                             pmicrons <- progressor()
                             pmicrons(message = "Creating micron-level molecule coordinates",
                                      class = "sticky", amount = 0)
                             df <- data.frame(x = mx$global_x, y = mx$global_y,
                                              gene = mx$gene, stringsAsFactors = FALSE)
                             pmicrons(type = "finish")
                             df
                           }, stop("Unknown MERFISH input type: ", type))
      }
    
    # add z-slice index for cells
    outs$zIndex <- data.frame(z = rep_len(z, length.out = outs$centroids$cell %>% length), cell = outs$centroids$cell)
    
    return(outs)
    gc() %>% invisible()
    
  }

## ==============================================================
## Loading data - output is a Seurat object
## ==============================================================

LoadVizgen_opt <- 
  function (data.dir, fov = "vz", assay = "Vizgen", 
            add.zIndex = TRUE, update.object = TRUE,
            ...)
  {
    data <- ReadVizgen_opt(data.dir = data.dir, ...)
    
    gc() %>% invisible()  
    # if "segmentations" are not present, use cell bounding boxes instead
    if (!"segmentations" %in% names(data)) { 
      bound.boxes <- CreateSegmentation(data[["boxes"]])
      cents <- CreateCentroids(data[["centroids"]])
      bound.boxes.data <- list(centroids = cents, boxes = bound.boxes)
      message("Creating FOVs", "\n",
              "..using box coordinates instead of segmentations")  
      coords <- CreateFOV(coords = bound.boxes.data, type = c("boxes",
                                                              "centroids"), molecules = data[[mol.type]], assay = assay)
      message("Creating Seurat object")  
      obj <- CreateSeuratObject(counts = data[["transcripts"]], assay = assay)
      coords <- subset(x = coords, 
                       cells = intersect(x = Cells(x = coords[["boxes"]]),
                                         y = Cells(x = obj)))
    } else {
      segs <- CreateSegmentation(data[["segmentations"]])
      cents <- CreateCentroids(data[["centroids"]])
      segmentations.data <- list(centroids = cents, segmentation = segs)
      message("Creating FOVs", "\n", 
              "..using segmentations")    
      coords <- CreateFOV(coords = segmentations.data, type = c("segmentation",
                                                                "centroids"), molecules = data[[mol.type]], assay = assay)
      message("Creating Seurat object..")  
      obj <- CreateSeuratObject(counts = data[["transcripts"]], assay = assay)
      coords <- subset(x = coords,
                       cells = intersect(x = Cells(x = coords[["segmentation"]]),
                                         y = Cells(x = obj)))    
    }
    
    # add z-stack index for cells
    if (add.zIndex) { obj$z <- data$zIndex$z }
   
    # add metadata vars
    message("..adding metadata infos")  
    if (c("metadata" %in% names(data))) {
      metadata <- match.arg(arg = "metadata", choices = names(data), several.ok = TRUE)
      meta.vars <- names(data[[metadata]])
      for (i in meta.vars %>% seq) {
        obj %<>% AddMetaData(metadata = data[[metadata]][[meta.vars[i]]], 
                             col.name = meta.vars[i])
      }
    }
    
    # sanity on fov name
    fov %<>% gsub("_|-", ".", .) 
     
    message("..adding FOV")
    obj[[fov]] <- coords
      
    if (update.object) { 
        message("Updating object..")
        obj %<>% UpdateSeuratObject() } 
    
    message("Object is ready!")    
    return(obj)
      
    gc() %>% invisible()
  }

