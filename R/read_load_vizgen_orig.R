

ReadVizgen <- 
function (data.dir, transcripts = NULL, spatial = NULL, molecules = NULL, 
    type = "segmentations", mol.type = "microns", metadata = NULL, 
    filter = NA_character_, z = 3L) 
{
    if (!requireNamespace("data.table", quietly = TRUE)) {
        stop("Please install 'data.table' for this function")
    }
    hdf5 <- requireNamespace("hdf5r", quietly = TRUE)
    type <- match.arg(arg = type, choices = c("segmentations", 
        "centroids", "boxes"), several.ok = TRUE)
    mol.type <- match.arg(arg = mol.type, choices = c("pixels", 
        "microns"), several.ok = TRUE)
    if (!is.null(x = metadata)) {
        metadata <- match.arg(arg = metadata, choices = c("volume", 
            "fov"), several.ok = TRUE)
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
    files <- c(transcripts = transcripts %||% "cell_by_gene[_a-zA-Z0-9]*.csv", 
        spatial = spatial %||% "cell_metadata[_a-zA-Z0-9]*.csv", 
        molecules = molecules %||% "detected_transcripts[_a-zA-Z0-9]*.csv")
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
            poly <- if (isFALSE(x = hdf5)) {
                warning("Cannot find hdf5r; unable to load segmentation vertices", 
                  immediate. = TRUE)
                FALSE
            }
            else if (!dir.exists(paths = h5dir)) {
                warning("Cannot find cell boundary H5 files", 
                  immediate. = TRUE)
                FALSE
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
    for (otype in names(x = outs)) {
        outs[[otype]] <- switch(EXPR = otype, transcripts = {
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
        }, segmentations = {
            ppoly <- progressor(steps = length(x = unique(x = sp$fov)))
            ppoly(message = "Creating polygon coordinates", class = "sticky", 
                amount = 0)
            pg <- future_lapply(X = unique(x = sp$fov), FUN = function(f, 
                ...) {
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
                cells <- rownames(x = subset(x = sp, subset = fov == 
                  f))
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
            })
            ppoly(type = "finish")
            pg <- do.call(what = "rbind", args = pg)
            npg <- length(x = unique(x = pg$cell))
            if (npg < nrow(x = sp)) {
                warning(nrow(x = sp) - npg, " cells missing polygon information", 
                  immediate. = TRUE)
            }
            pg
        }, boxes = {
            pbox <- progressor(steps = nrow(x = sp))
            pbox(message = "Creating box coordinates", class = "sticky", 
                amount = 0)
            bx <- future_lapply(X = rownames(x = sp), FUN = function(cell) {
                row <- sp[cell, ]
                df <- expand.grid(x = c(row$min_x, row$max_x), 
                  y = c(row$min_y, row$max_y), cell = cell, KEEP.OUT.ATTRS = FALSE, 
                  stringsAsFactors = FALSE)
                df <- df[c(1, 3, 4, 2), , drop = FALSE]
                pbox()
                return(df)
            })
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
    return(outs)
}


LoadVizgen <- 
function (data.dir, fov, assay = "Vizgen", z = 3L) 
{
    data <- ReadVizgen(data.dir = data.dir, filter = "^Blank-", 
        type = c("centroids", "segmentations"), z = z)
    segs <- CreateSegmentation(data$segmentations)
    cents <- CreateCentroids(data$centroids)
    segmentations.data <- list(centroids = cents, segmentation = segs)
    coords <- CreateFOV(coords = segmentations.data, type = c("segmentation", 
        "centroids"), molecules = data$microns, assay = assay)
    obj <- CreateSeuratObject(counts = data$transcripts, assay = assay)
    coords <- subset(x = coords, cells = intersect(x = Cells(x = coords[["segmentation"]]), 
        y = Cells(x = obj)))
    obj[[fov]] <- coords
    return(obj)
}
