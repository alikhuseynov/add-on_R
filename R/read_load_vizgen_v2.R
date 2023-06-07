#' @import dplyr
#' @import tibble
#' @import purrr
#' @import furrr
#' @import magrittr
#' @importFrom spatstat.geom is.empty
#' @importFrom progressr progressor
NULL

#' @param data.dir Path to the directory with Vizgen MERFISH files; requires at
#' least one of the following files present:
#' \itemize{
#'  \item \dQuote{\code{cell_by_gene.csv}}: used for reading count matrix
#'  \item \dQuote{\code{cell_metadata.csv}}: used for reading cell spatial
#'  coordinate matrices
#'  \item \dQuote{\code{detected_transcripts.csv}}: used for reading molecule
#'  spatial coordinate matrices
#' }
#' @param transcripts Optional file path for counts matrix; pass \code{NA} to
#' suppress reading counts matrix
#' @param spatial Optional file path for spatial metadata; pass \code{NA} to
#' suppress reading spatial coordinates. If \code{spatial} is provided and
#' \code{type} is \dQuote{segmentations}, uses \code{dirname(spatial)} instead of
#' \code{data.dir} to find HDF5 files
#' @param molecules Optional file path for molecule coordinates file; pass
#' \code{NA} to suppress reading spatial molecule information
#' @param type Type of cell spatial coordinate matrices to read; choose one
#' or more of:
#' \itemize{
#'  \item \dQuote{segmentations}: cell segmentation vertices; requires
#'  \href{https://cran.r-project.org/package=hdf5r}{\pkg{hdf5r}} to be
#'   installed and requires a directory \dQuote{\code{cell_boundaries}} within
#'   \code{data.dir}. Within \dQuote{\code{cell_boundaries}}, there must be
#'   one or more HDF5 file named \dQuote{\code{feature_data_##.hdf5}}
#'  \item \dQuote{centroids}: cell centroids in micron coordinate space
#'  \item \dQuote{boxes}: cell box outlines in micron coordinate space
#' }
#' @param mol.type Type of molecule spatial coordinate matrices to read;
#' choose one or more of:
#' \itemize{
#'  \item \dQuote{pixels}: molecule coordinates in pixel space
#'  \item \dQuote{microns}: molecule coordinates in micron space
#' }
#' @param metadata Type of available metadata to read;
#' choose zero or more of:
#' \itemize{
#'  \item \dQuote{volume}: estimated cell volume
#'  \item \dQuote{fov}: cell's fov
#' }
#' @param filter A character to filter molecules by, pass \code{NA} to skip
#'  molecule filtering
#' @param z Z-index to load; must be between 0 and 6, inclusive
#' @param use.BiocParallel If to use \code{BiocParallel::bplapply}, 
#' default is to \code{TRUE}, if \code{FALSE}, uses \code{future} library
#' @param workers.MulticoreParam Number of cores to use for \code{BiocParallel::bplapply}
#' @param DTthreads.pct Optional, set percentage eg \code{50} of total threads to use for \code{data.table::fread}, 
#' if set to \code{NULL} will use default setting as in \code{data.table::getDTthreads(verbose = T)}
#' @param use.furrr When working with lists, use \code{furrr} with \code{future} parallelization;
#'  if \conde{FALSE} standard \code{purrr} will be used
#' @param verbose If to print processing messages using \code{message()}, default is to \code{FALSE}
#'
#' @return \code{ReadVizgen}: A list with some combination of the
#' following values:
#' \itemize{
#'  \item \dQuote{\code{transcripts}}: a
#'  \link[Matrix:dgCMatrix-class]{sparse matrix} with expression data; cells
#'   are columns and features are rows
#'  \item \dQuote{\code{segmentations}}: a data frame with cell polygon outlines in
#'   three columns: \dQuote{x}, \dQuote{y}, and \dQuote{cell}
#'  \item \dQuote{\code{centroids}}: a data frame with cell centroid
#'   coordinates in three columns: \dQuote{x}, \dQuote{y}, and \dQuote{cell}
#'  \item \dQuote{\code{boxes}}: a data frame with cell box outlines in three
#'   columns: \dQuote{x}, \dQuote{y}, and \dQuote{cell}
#'  \item \dQuote{\code{microns}}: a data frame with molecule micron
#'   coordinates in three columns: \dQuote{x}, \dQuote{y}, and \dQuote{gene}
#'  \item \dQuote{\code{pixels}}: a data frame with molecule pixel coordinates
#'   in three columns: \dQuote{x}, \dQuote{y}, and \dQuote{gene}
#'  \item \dQuote{\code{metadata}}: a data frame with the cell-level metadata
#'   requested by \code{metadata}
#' }

## TODO:
# - add modifications to cloned gihub repo https://github.com/satijalab/seurat/tree/develop new branch `feat/vizgen`
# - add an option to read image data, store image as a raster similar to Visium `Read10X_Image`

## ==============================================================
## Reading data - output is a list of dataframes
## ==============================================================
ReadVizgen_opt <- function(data.dir,
                           transcripts = NULL,
                           spatial = NULL,
                           molecules = NULL,
                           type = c('centroids', 'segmentations'),
                           mol.type = 'microns',
                           metadata = NULL,
                           filter = NA_character_,	
                           z = 3L,
                           use.BiocParallel = TRUE,
                           workers.MulticoreParam = 12,
                           DTthreads.pct = NULL,
                           use.furrr = TRUE,
                           verbose = FALSE
)
{
  # TODO: handle multiple segmentations per z-plane
  # NOTE: this is only needed when segmentations differ between z-planes
  
  # packages that needs to be installed a priori
  pkgs <- c("data.table", "arrow", "sfarrow",
            "tidyverse", "furrr", "BiocParallel")
  lapply(pkgs %>% length %>% seq, function(i) 
  { !requireNamespace(pkgs[i], quietly = TRUE) } ) %>% 
    unlist %>% 
    { if (c(which(.) > 0) %>% any()) 
    { c("Please install ->", "\n",
        paste0("'", pkgs[which(.)], "'", collapse = ", "), " for this function") %>% 
        stop(., call. = FALSE) } }
  
  # setting workers to use for parallelization - `BiocParallel` ----
  if (use.BiocParallel) {
    if (verbose) { message("Using parallelization with: `BiocParallel`") }
    if (is.null(workers.MulticoreParam)) { 
      workers.MulticoreParam <- quantile(BiocParallel::multicoreWorkers() %>% seq)[4] %>% ceiling
      if (verbose) { message(workers.MulticoreParam, " of total workers available will be used") }
    } else { if (verbose) { message("Setting total workers to: ", workers.MulticoreParam) } }   
  } else { if (verbose) { message("Using parallelization with: `future`", "\n",
                                  "NOTE: set workers for parallelization, eg: `future::plan('multisession', workers = 10)`") }
  }
  # support parallelization on unix and windows
  # credit to https://github.com/Bioconductor/BiocParallel/issues/98
  BPParam <-
    if (.Platform$OS.type == "windows") {
      if (verbose) { message("Using parallelization for Windows with: `BiocParallel::SerialParam`") }
      BiocParallel::SerialParam(force.GC = FALSE, progressbar = TRUE)
    } else {
      BiocParallel::MulticoreParam(workers.MulticoreParam, tasks = 50L, 
                                   force.GC = FALSE, progressbar = TRUE)
    }
  
  # (optional) 
  # setting additional cores to use for parallelization in `data.table`
  # might allow to read large files faster
  if (!is.null(DTthreads.pct)) {
    if (verbose) { message("Using parallelization with: `data.table`", "\n",
                           "..for `data.table::fread`") }
    data.table::setDTthreads(threads = 0) # all cores
    DTthreads <- data.table::getDTthreads() # max cores
    DTthreads <- c((DTthreads * DTthreads.pct) / 100) %>% ceiling # percentage from total threads
    if (verbose) { message("Setting DTthreads to: ", DTthreads, " (", paste0(DTthreads.pct, "%"), ")") }
    data.table::setDTthreads(threads = DTthreads) # set
  }
  
  # hdf5r is only used for loading polygon boundaries
  # Not needed for all Vizgen input
  hdf5 <- requireNamespace("hdf5r", quietly = TRUE)
  # Argument checking
  type <- match.arg(
    arg = type,
    choices = c('segmentations', 'centroids', 'boxes'),
    several.ok = TRUE
  )
  mol.type <- match.arg(
    arg = mol.type,
    choices = c('pixels', 'microns'),
    several.ok = TRUE
  )
  if (!is.null(x = metadata)) {
    metadata <- match.arg(
      arg = metadata,
      choices = c("volume", "fov"),
      several.ok = TRUE
    )
  }
  if (!z %in% seq.int(from = 0L, to = 6L)) {
    stop("The z-index must be in the range [0, 6]")
  }
  use.dir <- all(vapply(
    X = c(transcripts, spatial, molecules),
    FUN = function(x) {
      return(is.null(x = x) || is.na(x = x))
    },
    FUN.VALUE = logical(length = 1L)
  ))
  if (use.dir && !dir.exists(paths = data.dir)) {
    stop("Cannot find Vizgen directory ", data.dir)
  } else { 
    if (verbose) { message("Reading data from:" , "\n", data.dir) }
  }
  
  # Use segmentation output from ".parquet" file ----
  # check if file exists
  parq <-
    # look in the current directory
    list.files(data.dir, 
               pattern = ".parquet$",
               full.names = TRUE) %>%
    { if (is.empty(.)) { 
      # look in the sub directory (if nothing is found)
      list.files(data.dir,
                 pattern = ".parquet$",
                 full.names = TRUE,
                 recursive = TRUE)
    } else { (.) }}
  # set to use .parquet" file if present
  use.parquet <- 
    parq %>% length() %>% any
  
  if (use.parquet) {
    if (verbose) { message("Cell segmentations are found in `.parquet` file", "\n", 
                           ">>> using ", gsub("s", "", mol.type), " space coordinates") }}
  # Identify input files..
  # if no files are found in the current directory..
  #..look for them in the sub directory
  files <- c(transcripts = transcripts %||% 
               list.files(data.dir, 
                          pattern = "cell_by_gene",
                          full.names = TRUE) %>%
               { if (is.empty(.)) { 
                 list.files(data.dir,
                            pattern = "cell_by_gene",
                            full.names = TRUE,
                            recursive = TRUE)
               } else { (.) }},
             
             spatial = spatial %||% 
               list.files(data.dir, 
                          pattern = "cell_metadata",
                          full.names = TRUE) %>%
               { if (is.empty(.)) { 
                 list.files(data.dir,
                            pattern = "cell_metadata",
                            full.names = TRUE,
                            recursive = TRUE)
               } else { (.) }},
             
             molecules = molecules %||% 
               list.files(data.dir, 
                          pattern = "detected_transcripts",
                          full.names = TRUE) %>%
               { if (is.empty(.)) { 
                 list.files(data.dir,
                            pattern = "detected_transcripts",
                            full.names = TRUE,
                            recursive = TRUE)
               } else { (.) }}
  )
  
  files[is.na(x = files)] <- NA_character_
  h5dir <- file.path(
    ifelse(
      test = dirname(path = files['spatial']) == '.',
      yes = data.dir,
      no = dirname(path = files['spatial'])
    ),
    'cell_boundaries'
  )
  zidx <- paste0('zIndex_', z)
  files <- vapply(
    X = files,
    FUN = function(x) {
      x <- as.character(x = x)
      if (isTRUE(x = dirname(path = x) == '.')) {
        fnames <- list.files(
          path = data.dir,
          pattern = x,
          recursive = FALSE,
          full.names = TRUE
        )
        return(sort(x = fnames, decreasing = TRUE)[1L])
      } else {
        return(x)
      }
    },
    FUN.VALUE = character(length = 1L),
    USE.NAMES = TRUE
  )
  files[!file.exists(files)] <- NA_character_
  if (all(is.na(x = files))) {
    stop("Cannot find Vizgen input files in ", data.dir)
  }
  # Checking for loading spatial coordinates
  if (!is.na(x = files[['spatial']])) {
    pprecoord <- progressor()
    pprecoord(
      message = "Preloading cell spatial coordinates",
      class = 'sticky',
      amount = 0
    )
    sp <- data.table::fread(
      file = files[['spatial']],
      sep = ',',
      data.table = FALSE,
      verbose = FALSE
      # showProgress = progressr:::progressr_in_globalenv(action = 'query')
      # showProgress = verbose
    )
    pprecoord(type = 'finish')
    rownames(x = sp) <- as.character(x = sp[, 1])
    #sp <- sp[, -1, drop = FALSE]
    if ((names(sp) == "transcript_count") %>% any) {
      if (verbose) { message(">>> filtering `cell_metadata` - keep cells with `transcript_count` > 0") }
      sp %<>% select(-1) %>% filter(transcript_count > 0)
    } else { sp %<>% select(-1) }
    
    # Check to see if we should load segmentations
    if ('segmentations' %in% type) {
      poly <- if (isFALSE(x = hdf5) && !use.parquet) {
        warning(
          "Cannot find hdf5r; unable to load segmentation vertices",
          immediate. = TRUE
        )
        FALSE
      } else if (!dir.exists(paths = h5dir) && !use.parquet) {
        warning("Cannot find cell boundary H5 files", immediate. = TRUE)
        FALSE
      } else if (use.parquet) { # for non .hdf5 files
        if (length(parq)) {
        }
        TRUE
      }
      else {
        TRUE
      }
      if (isFALSE(x = poly)) {
        type <- setdiff(x = type, y = 'segmentations')
      }
    }
    spatials <- rep_len(x = files[['spatial']], length.out = length(x = type))
    names(x = spatials) <- type
    files <- c(files, spatials)
    files <- files[setdiff(x = names(x = files), y = 'spatial')]
  } else if (!is.null(x = metadata)) {
    warning(
      "metadata can only be loaded when spatial coordinates are loaded",
      immediate. = TRUE
    )
    metadata <- NULL
  }
  # Check for loading of molecule coordinates
  if (!is.na(x = files[['molecules']])) {
    ppremol <- progressor()
    ppremol(
      message = "Preloading molecule coordinates",
      class = 'sticky',
      amount = 0
    )
    mx <- data.table::fread(
      file = files[['molecules']],
      sep = ',',
      verbose = FALSE
      # showProgress = verbose
    )
    mx <- mx[mx$global_z == z, , drop = FALSE]
    if (!is.na(x = filter)) {
      ppremol(
        message = paste("Filtering molecules with pattern", filter),
        class = 'sticky',
        amount = 0
      )
      mx <- mx[!grepl(pattern = filter, x = mx$gene), , drop = FALSE]
    }
    ppremol(type = 'finish')
    mols <- rep_len(x = files[['molecules']], length.out = length(x = mol.type))
    names(x = mols) <- mol.type
    files <- c(files, mols)
    files <- files[setdiff(x = names(x = files), y = 'molecules')]
  }
  files <- files[!is.na(x = files)]
  
  # Read input data ----
  outs <- vector(mode = 'list', length = length(x = files))
  names(x = outs) <- names(x = files)
  if (!is.null(metadata)) {
    outs <- c(outs, list(metadata = NULL))
  }
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
               # avoid converting to dense matrix?
               #tx <- t(x = as.matrix(x = tx[, -1, drop = FALSE]))
               
               # keep cells with `transcript_count` > 0
               if (verbose) { message(">>> filtering `cell_by_gene` - keep cells with counts > 0") }
               tx %<>% select(-1) %>%
                 filter_all(any_vars(. > 0)) %>% {
                   if ((names(sp) == "transcript_count") %>% any) {
                     # match cells to filtered data from spatial (cell_metadata)
                     filter(., rownames(.) %in% rownames(sp)) 
                   } else { (.) } } %>% # return filtered count data
                 t() %>% as.sparse()
               
               # filter cell metadata df
               # match filtered cell IDs from count matrix to cell metadata df
               if (!(names(sp) == "transcript_count") %>% any) {
                 sp %<>% filter(rownames(.) %in% colnames(tx))
               }
               
               ratio <- getOption(x = "Seurat.input.sparse_ratio",
                                  default = 0.4)
               if ((sum(tx == 0)/length(x = tx)) > ratio) {
                 ptx(message = "Counts are converted to sparse matrix",
                     class = "sticky", amount = 0)
                 #tx <- as.sparse(x = tx)
               }
               if (!is.na(x = filter)) {
                 ptx(message = paste("Filtering genes with pattern",
                                     filter), class = "sticky", amount = 0)
                 tx <- tx[!grepl(pattern = filter, x = rownames(x = tx)),
                          , drop = FALSE]
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
               
               # use segmentations from ".parquet"
             }, segmentations = {
               if (use.parquet) {
                 if (length(parq) > 1) {
                   # eg, if two files are present:
                   # `cellpose_micron_space.parquet`
                   # `cellpose_mosaic_space.parquet`
                   parq %<>% { 
                     if (mol.type == "pixels") {
                       # for `cellpose_mosaic_space.parquet`
                       grep("mosaic", ., value = TRUE)
                     } else if (mol.type == "microns") {
                       # for `cellpose_micron_space.parquet`
                       grep(gsub("s", "", mol.type), ., value = TRUE)
                     }
                     } %>%
                     { if (is.empty(.)) {
                       # only if single ".parquet" file present
                       # eg, `cell_boundaries.parquet`
                       parq %>%
                         grep(".parquet$", ., value = TRUE)
                     } else { (.) } }
                 }
                 
                 # Read .parquet file
                 parq %<>% sfarrow::st_read_parquet(.)
                 
                 # get all cell segmentations
                 segs <- filter(parq, ZIndex == z) %>% 
                   pull(Geometry) %>%
                   # add cell ID
                   set_names(filter(parq, ZIndex == z) %>% 
                               pull(EntityID) %>% as.character)
                 
                 ## Sanity checks on segmentation polygons - part 1 ----
                 test.segs <- lapply(segs %>% seq, 
                                     function(i) segs[[i]] %>% length) %>% unlist()
                 if (any(test.segs > 1)) {
                   segs.art.index <- which(test.segs > 1)
                   segs.empty.index <- which(test.segs < 1)
                   if (verbose) { message("Sanity checks on cell segmentation polygons:", "\n", 
                                          ">>> found ", segs.art.index %>% length,
                                          " cells with > 1 (nested) polygon lists:", "\n",
                                          ">>> flattening polygon lists", "\n",
                                          if (c(segs.empty.index %>% length) > 0) {
                                            paste0(">>> removing ", 
                                                   segs.empty.index %>% length,
                                                   " empty polygon lists") }
                   ) }
                   
                   # step 1 - find polygon nested lists with > 1 length
                   # get original planned session if exists
                   if (is(future::plan(), "multisession")) {
                     orig.plan <- future::plan()
                   }
                   segs.1 <- segs %>% 
                     purrr::keep(., c(purrr::map(., length) %>% unlist > 1)) %>%
                     # flatten each list
                     { 
                       if (use.furrr) {
                         # set temporary workers
                         f.plan <- future::plan("multisession", workers = 4L, gc = TRUE)
                         #on.exit(f.plan %>% future::plan()) # exiting doesn't to work with furrr
                         # using furrr (faster purrr with future) 
                         furrr::future_map(., purrr::flatten)
                       } else {
                         purrr::map(., purrr::flatten)
                       }
                     } %>% suppressPackageStartupMessages()
                   # set originally plannned session back
                   if (exists("orig.plan")) { 
                     future::plan(orig.plan)
                   } else { 
                     # or plan sequential
                     future::plan("sequential")
                   }
                   
                   # step 2 - apply multiple filtering  
                   segs.2 <- 
                     segs %>%
                     # remove empty elements               
                     purrr::keep(., !c(purrr::map(., length) %>% unlist < 1)) %>%
                     # keep lists with length == 1 polygon information      
                     purrr::keep(., c(purrr::map(., length) %>% unlist == 1)) %>%
                     # collapse to a list
                     purrr::flatten(.)
                   
                   # step 3 - combine segmentaion lists     
                   segs <- c(segs.1, segs.2) 
                 }
                 
                 ## Sanity checks on segmentation polygons - part 2 ----
                 # check if any cells have > 1 segmentation boundary
                 # TODO: (optionally) optimize sanity code using purrr?
                 test.segs <- lapply(segs %>% seq, function(i) segs[[i]] %>% length) %>% unlist()
                 if (any(test.segs %>% unlist > 1)) {
                   segs.art.index <- which(test.segs %>% unlist > 1)
                   if (verbose) { 
                     message("Found ", segs.art.index %>% length,
                             c(" cells with > 1 polygon artifacts:",
                               ">>> removing artifacts",
                               ">>> keeping cell boundary with maximum coords") %>% 
                               paste0(., "\n"))
                   }
                   # usually artifacts have small boundaries/coords
                   # find cell boundaries with maximum coords
                   
                   # TODO: (optionally)
                   # - calculate geometric properties, eg circularity:
                   # - keep single polygon with high circularity values (likely a cell?)
                   
                   for (i in segs.art.index %>% seq) { 
                     dims <- lapply(segs[[segs.art.index[i]]] %>% seq(), 
                                    function(d) { dim(segs[[segs.art.index[i]]][[d]])[1] } )
                     # get & keep boundaries with maximum coords
                     maxs.segs <- which(unlist(dims) == max(unlist(dims)))
                     segs[[segs.art.index[i]]] <- segs[[segs.art.index[i]]][maxs.segs]
                   }
                 } else { if (verbose) { message("All cells have single polygon boundary (no artifacts)") } }                 
                 
                 # some cells might have > 1 polygon boundary with identical lenght
                 #..in this case, keep only the 1st polygon boundary
                 test.segs <- lapply(segs %>% seq, function(i) segs[[i]] %>% length) %>% unlist()                
                 if (any(test.segs > 1)) {
                   segs.art.index <- which(test.segs > 1)
                   if (verbose) { message("Additionally found ", segs.art.index %>% length,
                                          c(" cells with > 1 polygons (identical length):",
                                            ">>> only the 1st polygon boundary will be kept") %>% 
                                            paste0(., "\n")) }
                   for (i in segs.art.index %>% seq) { 
                     # TODO: (optionally)
                     # - calculate geometric properties, eg circularity:
                     # - keep single polygon with high circularity values (likely a cell?)
                     
                     segs[[segs.art.index[i]]] <- segs[[segs.art.index[i]]][1]
                   }
                 }              
                 
                 ## Extract cell boundaries per cells
                 # TODO: (optionally) resample & make cell boundaries equidistant? 
                 # TODO: (optionally) optimize with purrr::map to get df per list instance?
                 if (use.BiocParallel) {
                   gc() %>% invisible() # free up memory
                   if (verbose) { message("Extracting cell segmentations - using `BiocParallel`") } 
                   segs_list <-
                     BiocParallel::bplapply(segs %>% seq,
                                            function(i) {
                                              segs[[i]][[1]] %>% 
                                                data.table::as.data.table(.) %>%
                                                mutate(cell = names(segs)[i]) },
                                            BPPARAM = BPParam)
                 } else {
                   if (verbose) { message("Extracting cell segmentations - using `future`") }    
                   segs_list <-
                     future.apply::future_lapply(segs %>% seq,
                                                 function(i) {
                                                   segs[[i]] %>%
                                                     data.table::as.data.table(.) %>%
                                                     mutate(cell = names(segs)[i])
                                                 }
                     )
                 }
                 
                 # df of all cell segmentations
                 segs <- 
                   data.table::rbindlist(segs_list) %>% 
                   data.table::setnames(c("x", "y", "cell"))
                 #names(segs)[1:2] <- c("x", "y")
                 if (verbose) { message("All cell segmentations are loaded..") }      
                 segs
               } else { 
                 # else use ".hdf5" files from ./cell_boundaries (older version)
                 ppoly <- progressor(steps = length(x = unique(x = sp$fov)))
                 ppoly(message = "Creating polygon coordinates", class = "sticky",
                       amount = 0)
                 # use `BiocParallel` or `future`
                 if (use.BiocParallel) {
                   if (verbose) { message("Reading '.hdf5' files..") }
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
                     #return(do.call(what = "rbind", args = df))
                     return(data.table::rbindlist(df))
                   }, BPPARAM = BPParam)
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
                 if (verbose) { message("Creating box coordinates..") }
                 bx <- BiocParallel::bplapply(X = rownames(x = sp), FUN = function(cell) {
                   row <- sp[cell, ]
                   # faster version for grid construction
                   df <- data.table::CJ(x = c(row$min_x, row$max_x),
                                        y = c(row$min_y, row$max_y), 
                                        cell = cell) %>% 
                     slice(c(1, 3, 4, 2))
                   #df <- expand.grid(x = c(row$min_x, row$max_x),
                   #                 y = c(row$min_y, row$max_y), cell = cell, KEEP.OUT.ATTRS = FALSE,
                   #                stringsAsFactors = FALSE)
                   #df <- df[c(1, 3, 4, 2), , drop = FALSE]
                   pbox()
                   return(df)
                 }, BPPARAM = BPParam)
               } else {
                 bx <- future.apply::future_lapply(X = rownames(x = sp), FUN = function(cell) {
                   row <- sp[cell, ]
                   df <- data.table::CJ(x = c(row$min_x, row$max_x),
                                        y = c(row$min_y, row$max_y), cell = cell) %>%
                     slice(c(1, 3, 4, 2))
                   #df <- expand.grid(x = c(row$min_x, row$max_x),
                   #y = c(row$min_y, row$max_y), cell = cell, KEEP.OUT.ATTRS = FALSE,
                   #stringsAsFactors = FALSE)
                   #df <- df[c(1, 3, 4, 2), , drop = FALSE]
                   pbox()
                   return(df)
                 })
               }
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
  
  # add z-slice index for cells ----
  outs$zIndex <- 
    data.frame(z = rep_len(z, length.out = outs$centroids %>% pull(cell) %>% length), 
               cell = outs$centroids %>% pull(cell))
  
  return(outs)
}

## ==============================================================
## Loading data - output is a Seurat object
## ==============================================================

#' @return \code{LoadVizgen}: A \code{\link[SeuratObject]{Seurat}} object
#'
#' @param fov Name to store FOV as
#' @param assay Name to store expression matrix as
#' @param add.zIndex If to add \code{z} slice index to a cell
#' @param update.object If to update final object, default to TRUE
#' @param ... Arguments passed to \code{ReadVizgen}
#'
#' @importFrom SeuratObject Cells CreateCentroids CreateFOV
#' CreateSegmentation CreateSeuratObject
#' @import dplyr
#'
#' @export
#'
#' @rdname ReadVizgen

LoadVizgen_opt <- function(data.dir, 
                           fov = 'vz', 
                           assay = 'Vizgen',
                           mol.type = 'microns',
                           filter = '^Blank-',
                           z = 3L,
                           add.zIndex = TRUE, 
                           update.object = TRUE,
                           verbose,
                           ...)
{
  # reading data..
  data <- ReadVizgen_opt(data.dir = data.dir,
                         mol.type = mol.type,
                         filter = filter,
                         z = z,
                         verbose = verbose,
                         ...)
  
  if (verbose) { message("Creating Seurat object..") }  
  obj <- CreateSeuratObject(counts = data[["transcripts"]], assay = assay)
  
  # in case no segmentation is present, use boxes
  if (!"segmentations" %in% names(data)) {
    if ("boxes" %in% names(data)) {
      bound.boxes <- CreateSegmentation(data[["boxes"]])
      cents <- CreateCentroids(data[["centroids"]])
      bound.boxes.data <- list(centroids = cents, 
                               boxes = bound.boxes)
      if (verbose) { 
        message("Creating FOVs..", "\n",
                ">>> using box coordinates instead of segmentations") 
      }
      coords <- CreateFOV(coords = bound.boxes.data, 
                          type = c("boxes", "centroids"),
                          molecules = data[[mol.type]], 
                          assay = assay)
    } else { 
      # in case no segmentation & no boxes are present, use centroids only
      cents <- CreateCentroids(data[["centroids"]])
      if (verbose) { 
        message("Creating FOVs..", "\n",
                ">>> using only centroids") 
      }
      coords <- CreateFOV(coords = list(centroids = cents), 
                          type = c("centroids"),
                          molecules = data[[mol.type]], 
                          assay = assay)
      coords <- subset(x = coords, 
                       cells = intersect(x = Cells(x = coords[["centroids"]]),
                                         y = Cells(x = obj))) 
    }
  } else if ("segmentations" %in% names(data)) {
    segs <- CreateSegmentation(data[["segmentations"]])
    cents <- CreateCentroids(data[["centroids"]])
    segmentations.data <- list(centroids = cents, segmentation = segs)
    if (verbose) { 
      message("Creating FOVs..", "\n", 
              ">>> using segmentations") 
    }
    coords <- CreateFOV(coords = segmentations.data, 
                        type = c("segmentation", "centroids"), 
                        molecules = data[[mol.type]], 
                        assay = assay)
    # only consider the cells we have counts and a segmentation.
    # Cells which don't have a segmentation are probably found in other z slices.
    coords <- subset(x = coords,
                     cells = intersect(x = Cells(x = coords[["segmentation"]]),
                                       y = Cells(x = obj)))
  }
  
  # add z-stack index for cells
  if (add.zIndex) { obj$z <- data$zIndex %>% pull(z) }
  
  # add metadata vars
  if (verbose) { message(">>> adding metadata infos") }
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
  
  if (verbose) { message(">>> adding FOV") }
  obj[[fov]] <- coords
  
  ## filter - keep cells with counts > 0
  # helper function to return metadata
  callmeta <- function (object = NULL) { return(object@meta.data) }
  nCount <- grep("nCount", callmeta(obj) %>% names, value = TRUE)
  if (any(obj[[nCount]] == 0)) {
    if (verbose) { message(">>> filtering object - keeping cells with counts > 0") }
    obj %<>% subset(subset = !!base::as.symbol(nCount) > 0)
  } else { if (verbose) { message(">>> all counts are > 0") } }
  
  if (update.object) { 
    if (verbose) { message("Updating object:") 
      obj %<>% UpdateSeuratObject()
    } else { 
      obj %<>% 
        UpdateSeuratObject() %>% 
        suppressMessages() } }
  
  if (verbose) { message("Object is ready!") } 
  return(obj)
  
}
