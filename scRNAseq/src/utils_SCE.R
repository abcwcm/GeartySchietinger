#' Add additional entries to the colData of an SCE
#' @param SCE.in SCE object
#' @param to_add data.frame with additional colData entries
#' @param merge_on Indicate the column of \code{to_add}, which should contain the
#' entries corresponding to the colnames(SCE.in). Default: NULL, which will use
#'  rownames of the \code{to_add} df.
#' @author Friederike Dündar
#' @return SCE object with amended colData entries
#' @import SingleCellExperiment
.add_colData <- function(SCE.in, to_add, merge_on = NULL){

    to_add <- as.data.frame(to_add)

    if(is.null(merge_on)){
        cnames <- rownames(to_add)
    }else{
        cnames <- to_add[, merge_on, drop=TRUE]
    }
    if(!all(cnames %in% colnames(SCE.in))){stop("The rownames/merge_on entries
        of the data.frame to be added aren't all part of the colnames of the SCE")}
    if(!all(colnames(SCE.in) %in% cnames)){message("There are more cells in the
        SCE object than entries in the df; there will be NAs in the resulting colData")}

    cd <- colData(SCE.in)
    to_add$MERGE <- cnames
    if(!is.null(merge_on) && merge_on %in% names(cd)){
        to_add[, merge_on] <- NULL
    }
    cdnames <- names(to_add)
    dblnms <- cdnames[cdnames %in% names(cd)]
    if(length(dblnms > 0)){
        cdnames[cdnames %in% dblnms] <- paste0(dblnms, ".y")
        names(to_add) <- cdnames
    }
    cd$MERGE <- colnames(SCE.in)
    cd.out <- merge(cd, to_add, by = "MERGE", all.x = TRUE)
    rownames(cd.out) <- cd.out$MERGE
    cd.out$MERGE <- NULL

    sce.out <- SCE.in
    colData(sce.out) <- DataFrame(cd.out[colnames(sce.out),])
    return(sce.out)

}


.extract_colData <- function(SCE_in, clonotype_col){
    met_df <- as.data.frame(colData(SCE_in))
    met_df <- tibble::as_tibble(met_df, rownames = ".cell_id")
    met_df <- dplyr::filter(met_df, !is.na(!!sym(clonotype_col)))

    return(met_df)
}

.wrestle_tcr_df <- function(tcr_df, clonotype_col){

    met_df <- tibble::as_tibble(tcr_df, rownames = ".cell_id")
    met_df <- dplyr::filter(met_df, !is.na(!!sym(clonotype_col)))

    return(met_df)
}
