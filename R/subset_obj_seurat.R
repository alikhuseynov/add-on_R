#' @importFrom magrittr %>% %<>%
NULL

#' Intermediate solution to \code{subset()}:
#' subset FOVs/centroids if selected cells are NOT found in each FOV
#' NOTE: some code parts and args are takes from SeuratObject

#' Function params/args:
#' @param object An S4 object
#' @param cells A vector of cells to keep; if \code{NULL}, defaults to all cells
#' @param update.slots If to update slots of an object
#' @param idents A vector of identity classes to keep
#' @param ... Arguments passed to \code{subset()} and other methods

subset_opt <- 
function (object = NULL, subset, 
          cells = NULL, idents = NULL, 
          update.slots = FALSE, ...)

{

if (update.slots) { 
    message("Updating object slots..")
    object <- UpdateSlots(object = object) 
}

message("Cloing object..")
obj_subset <- object

if (!missing(x = subset)) {
    subset <- enquo(arg = subset)
    # cells to keep in the object
    message("Extracting cells matched to `subset`, `idents`")
    cells <-
    WhichCells(object = obj_subset, 
               cells = cells,
               idents = idents, 
               expression = subset, 
               return.null = TRUE, ...)
}
    
# check if cells are present in all FOV
message("Matching cells in FOVs..")
cells_check <-
lapply(Images(obj_subset) %>% seq, 
       function(i) { 
           any(obj_subset[[Images(obj_subset)[i]]][["centroids"]] %>% Cells %in% cells) 
       }) %>% unlist
if (all(cells_check)) { 
    message("Cell subsets are found in all FOVs!", "\n",
            "Subsetting object..") 
    obj_subset %<>% base::subset(cells = cells, ...)
}

# if cells are present only in one or several FOVs:
# subset FOVs
fovs <- 
lapply(Images(obj_subset) %>% seq, function(i) {
    if (any(obj_subset[[Images(obj_subset)[i]]][["centroids"]] %>% Cells %in% cells)) {
        message("Cell subsets are found only in FOV: ", "\n", Images(obj_subset)[i])
        message("Subsetting Centroids..")
        base::subset(x = obj_subset[[Images(obj_subset)[i]]], cells = cells)
}
})
    
# replace subsetted FOVs, and remove FOVs with no matching cells
message("removing FOVs in which no cells are found: ", "\n", 
        paste0(Images(object)[which(!cells_check == TRUE)], "\n"), "\n",
        "Subsetting cells..")
for (i in fovs %>% seq) { obj_subset[[Images(object)[i]]] <- fovs[[i]] }  

# subset final object
obj_subset %<>% base::subset(cells = cells, ...)
    
UpdateSeuratObject(obj_subset)
return(obj_subset)
}
