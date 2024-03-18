#' Custom Crop function to remove any tissue area
#' @param x A data.frame with \code{x}, \code{y} and eg \code{cell} (cell barcode or ID) variable
#' @param object A \code{Seurat} object. NOTE, currently works on object with single FOV only. 
#' @param col_id A \code{character} vector specifing which variable has the cell ids
#' @param xy_pts A data.frame of xy point coordinates, must have 3 or more xy points, 
#'  it can also be a list of data frames. This will generate a convex hull for cropping
#' @param c_hull_include Everything under (\code{TRUE}) convex hull polygon is included or cropped out (\code{FALSE})
#' @param crop_molecules When \code{Seurat} is present, if to crop molecule cooridnates
#'  NOTE, this can take time especially when there are many molecules.
#' @param BPPARAM description
#' @return Returns cropped data.frame with same variables as input.
#'  Or cropped \code{Seurat} object.
#' @import sf
#' @import BiocParallel
#' @import SeuratObject
#' @importFrom magrittr %>% %<>%
#' @importFrom dplyr filter pull transmute
#' @examples
#' # NOTE, `sf` package must be installed!
#' # provide xy coords to make convex hull polygon
#' xy_pts <- 
#' data.frame("x" = c(5400, 11000, 5400, 6100, 7500, 11000),
#' "y" = c(0, 0, 1500, 2000, 4100, 4100))
#' # run crop function on Seurat object (image-based spatial object, like Xenium, Vizgen, etc..)
#' obj_crop <- 
#' Crop_custom(object = seurat_object, col_id = c("cell"), xy_pts = xy_pts,
#' c_hull_include = TRUE,
#' crop_molecules = TRUE,
#' BPPARAM = BiocParallel::MulticoreParam(5, tasks = 10L,force.GC = FALSE, progressbar = TRUE))
#' # run crop function on any dataframe that contains xy coords
#' df_crop <- 
#' Crop_custom(x = df_mols, col_id = c("molecule"), xy_pts = xy_pts,
#' c_hull_include = TRUE,
#' BPPARAM = BiocParallel::MulticoreParam(5, tasks = 10L,force.GC = FALSE, progressbar = TRUE))
#'
Crop_custom <- 
  function(x = NULL, 
           object = NULL,
           col_id = c("cell"),
           xy_pts = NULL,
           c_hull_include = TRUE,
           crop_molecules = TRUE,
           BPPARAM = BiocParallel::SerialParam()) {
    
    # check packages
    pkgs <- c("data.table", "sf", "BiocParallel", 
              "tidyverse", "magrittr")
    lapply(pkgs %>% seq(), function(i)
    { !requireNamespace(pkgs[i], quietly = TRUE) } ) %>% 
      unlist() %>% 
      { if (c(which(.) > 0) %>% any()) 
      { c("Please install ->", "\n",
          paste0("'", pkgs[which(.)], "'", collapse = ", "), " for this function") %>% 
          stop(., call. = FALSE) } }
    
    # check inputs
    if (is.null(xy_pts)) {
      stop(">>> Please provide xy point coordinates in `xy_pts`")
    }
    
    if (!is.null(object)) {
      if (is(object, "Seurat")) {
        message(">>> Using `Seurat` object")
        x <- NULL
        df_xy <- object[[Images(object)[1]]] %>% GetTissueCoordinates()
      }
    } else if (!is.null(x)) {
      # check x
      is_x <-
      inherits(x = x, 
               what = c("data.frame", "data.table", "tibble", "matrix"))
      #grep("data.frame|data.table|tibble|matrix", 
      #       class(x)) %>% any()
      if (is_x) { 
        df_xy <- x
        object <- NULL
      }
    } else if (is.null(object) && is.null(x)) {
      stop(">>> Please provide either `object` or `x` data.frame")
    }
    
    # make sf data.frame
    sf_df <- st_as_sf(df_xy, coords = c("x", "y"))
    # make convex hull
    if (is(xy_pts, "list")) {
      xy_pts <- data.table::rbindlist(xy_pts)
    }
    c_hull <- 
      st_as_sf(xy_pts, coords = c("x", "y")) %>% 
      st_combine() %>% st_convex_hull()
    crop_df <-
      st_intersection(sf_df, c_hull)
    
    # subset data.frame given the cell ids
    # use col_id instead of $cell
    df_xy %<>% {
      if (c_hull_include) { 
        filter(., !!as.symbol(col_id) %in% pull(crop_df, col_id))
      } else {
        filter(., !(!!as.symbol(col_id)) %in% pull(crop_df, col_id))
      }
    }
    
    # TODO make more cleaner code for cropping df ----
    # output data.frame ----
    if (is.null(object) && !is.null(x)) {
      message(">>> Cropping `data.frame`")
      if (c_hull_include) {
        # faster with `st_join`
        mols <-
          st_join(x = sf_df, 
                  join = st_within, 
                  left = FALSE,
                  y = st_sf(geometry = c_hull))
      } else {
        #mols <- st_difference(sf_df, c_hull)
        mols <-
          st_join(x = sf_df, 
                  join = st_disjoint,
                  left = FALSE,
                  y = st_sf(geometry = c_hull))
      }
      
      genes <- mols %>% pull(col_id) %>% unique()
      mols <-
        bplapply(genes %>% seq(), function(i) {
          mols %>%
            
            # TODO: use col_id instead of molecule, eg !!as.symbol(col_id)
            filter(molecule == genes[i]) %>%
            st_geometry() %>%
            st_coordinates() %>%
            as.data.frame() %>%
            transmute(x = X, y = Y, 
                      molecule = genes[i])
        }, BPPARAM = BPPARAM) %>%
        data.table::rbindlist()
      message(">>> Return output: `data.frame`")
      if (crop_molecules) {
        return(mols)
      } else { return(df_xy) }
    }
    
  
    # cropping Seurat object ----
    # output Seurat
    if (!is.null(object)) {
      message(">>> Return output - `Seurat` object")
      # subset object, centroids and segmentations as well
      object %<>%
        subset(x = ., 
               cells = intersect(x = colnames(x = .),
                                 y =  pull(df_xy, col_id)))
      
      # crop molecules ----
      fov <- object[[Images(object)[1]]]
      if (crop_molecules && 
          grep("molecules", names(fov)) != 0) {
        message(">>> Cropping molecule coordinates - might take time!")
        # using previously made convex hull
        sf_df_mols <- st_as_sf(fov[["molecule"]] %>% GetTissueCoordinates(), 
                               coords = c("x", "y"))
        if (c_hull_include) {
          #mols <- st_intersection(sf_df_mols, c_hull)
          # faster with `st_join`
          mols <-
            st_join(x = sf_df_mols, 
                    join = st_within, 
                    left = FALSE,
                    y = st_sf(geometry = c_hull))
        } else {
          #mols <- st_difference(sf_df_mols, c_hull)
          mols <-
            st_join(x = sf_df_mols, 
                    join = st_disjoint,
                    left = FALSE,
                    y = st_sf(geometry = c_hull))
          }
        genes <- mols$molecule %>% unique()
        mols <-
          bplapply(genes %>% seq(), function(i) {
            mols %>%
              filter(molecule == genes[i]) %>%
              st_geometry() %>%
              st_coordinates() %>%
              as.data.frame() %>%
              transmute(x = X, y = Y, 
                        gene = genes[i])
          }, BPPARAM = BPPARAM)
        # Create Molecule FOV only
        mols %<>% 
          data.table::rbindlist() %>%
          CreateMolecules()
        # replace and add to FOV of the object
        object[[Images(object)[1]]][["molecule"]] <- mols
        
        # TODO: make sure that cropped mols are added to object? ----
        #..and no mols are present in cropped out regions, since few mols were still present.
        # test with GetTissueCoordinates() and plot them!
        
        
        
      }
      validObject(object = object)
      return(object)
    }
  }
