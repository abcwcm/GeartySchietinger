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

#' Prepare a data.table in skinny format suitable for ggplotting
#'
#' @param object \code{SCE} object
#' @param exprs_values what type of expression values should be used from the
#' \code{SCE} object (default: "logcounts"); can be more than one name
#' (e.g., \code{names(assay(object))}).
#' @param features vector of feature names (e.g. gene names, HTOs, ADTs) 
#' for which the \code{data.table} should be made
#' (optional; default: NULL)
#' @param include_metaData indicate whether the resulting \code{data.table}
#' should also include specified entries of colData(object). Cacn include specific
#' gene or other feature names. Default: NULL. 
#'
#'
#' @return \code{data.table} in skinny format with additional columns
#' "cell" and whatever was specified via \code{include_metaData}
#' 
#' @author Friederike Dündar
#'
#' @examples \dontrun{
#' # get a long data.table with values from all assays
#' exp.dt <- make_long_dt(sce, exprs_values = names(assays(sce)),
#'                       features = c("SHOX2", "MYH6"),
#'                       include_metaData= c("SampleName", "Condition","Replicate"))
#' }
#'
#' @seealso \code{\link{fx.extract_exprs}}
#' 
#' @import data.table
#' @import scater
#' 
#' @export
#' 
make_long_dt <- function (object, exprs_values = "logcounts", features = NULL, 
    include_metaData = NULL){
    
    ## define genes
    if (!is.null(features)) {
        features <- unique(features)
        genes <- features[features %in% rownames(object)]
        obj <- object[genes, ]
        gn <- genes
    }
    else {
        obj <- object
        gn <- rownames(object)
    }
    
    ## extract the values for all indicated exprs_value types -------------------
    dt_list <- lapply(exprs_values, function(x){
        fx.extract_exprs(obj, x, gn) } )
    out.dt <- rbindlist(dt_list)
    
    ## check altExp for remaining features
    if(length(features[!features %in% rownames(object)]) > 0){
        
        for(EXP in altExpNames(object)){
            
            obj.alt <- altExp(object, EXP)
            gn.alt <- features[features %in% rownames(obj.alt)]
            
            if(length(gn.alt) > 0){
                ev <- exprs_values[exprs_values %in% assayNames(obj.alt)]
                out.alt <- rbindlist(lapply(ev, function(x){
                    fx.extract_exprs(obj.alt, x, gn.alt)}))
                out.dt <- rbindlist(list(out.dt, out.alt))
            }
        }
    }
    
    out.dt <- data.table::dcast(data = out.dt, feature_name + cell ~ value_type, value.var = "value")
    
    ## add colData -----------------------------------------------------
    if(!is.null(include_metaData)){
        found_mt <- FALSE
        varLabs <- names(colData(obj))
        if (any(include_metaData %in% varLabs )) {
            inc_mt <- include_metaData[include_metaData %in%  varLabs]
            meta_info <- as.data.table(colData(obj)[inc_mt])
            meta_info$cell <- colnames(obj)
            out.dt <- meta_info[out.dt, on = "cell"]
            found_mt <- TRUE
        }
        
        ## add gene values as meta data
        gnNames <- rownames(object)
        if (any(include_metaData %in% gnNames)) {
            inc_mt <- include_metaData[include_metaData %in% gnNames]
            if(!exists("out.dt")){out.dt <- NULL}
            out.dt <- fx.features_as_metaData(object, exprs_values = exprs_values,
                which_feat = inc_mt, meta.dt = out.dt)
            found_mt <- TRUE
        }
        
        ## add additional features' values as meta data
        for(EXP in altExpNames(object)){
            obj.alt <- altExp(object, EXP)
            gn.alt <- rownames(obj.alt)
            if(any(include_metaData %in% gn.alt)){
                inc_mt <- include_metaData[include_metaData %in% gn.alt]
                if(!exists("out.dt")){out.dt <- NULL}
                ev <- exprs_values[exprs_values %in% assayNames(obj.alt)]
                out.dt <- fx.features_as_metaData(obj.alt, exprs_values = ev,
                    which_feat = inc_mt, meta.dt = out.dt)
                found_mt <- TRUE
            }
        }
        
        if (!found_mt) {
            if (!(all(include_metaData == FALSE))) {
                warning("Not all of the labels indicated by you via include_metadata are part of the sce.object. There will be no metadata returned.")
            }
        }
    }
    return(out.dt)
}

#' Features as meta data
#' @description Extracts expression values for individual features and adds
#' them as a column 
#' @param object SCE object
#' @param exprs_values 
#' @param which_feat feature name(s)
#' @param meta.dt data.table with already extracted metadata
#' 
#' @author Friederike Dündar
#' 
fx.features_as_metaData <- function(object,exprs_values, which_feat, meta.dt = NULL){
    
    mt_gn <- lapply(exprs_values, function(x) {
        tmout <- as.data.frame(assay(object[unique(which_feat), ], x))
        tmout$feature_name <- row.names(tmout)
        tmout <- as.data.table(tmout)
        tmout$value_type <- x
        return(tmout)
    })
    
    mt_gn <- melt.data.table(
        rbindlist(mt_gn), 
        variable.name = "cell", 
        id.vars = c("value_type", "feature_name"))
    mt_gn <- dcast.data.table(mt_gn, 
        cell ~ feature_name + value_type, sep = ".")
    
    ## merge with other metadata
    if(!is.null(meta.dt)){
        outt <- mt_gn[meta.dt, on = "cell"]
    }else{
        outt <- mt_gn
    }
    
    return(outt)
}

#' Extract expression values from SCE object
#'
#' @description Wrapper function for extracting expression values into a data.table
#' from an SCE object including checking the output.
#'
#' @author Friederike Dündar
#' 
#' @param object SCE object
#' @param exprs_values which type of \code{assay} values to extract
#' @param genes vector of gene names (subset of rownames)
#'
#' @return molten data.table with one row per gene and cell and the corresponding
#' \code{exprs_values} and their type in a third and fourth column.
#' @seealso \code{\link{make_long_dt}}
#'
#' @import magrittr
#' @import data.table
#' 
fx.extract_exprs <- function (object, exprs_values, genes){
    
    long.dt <- assay(object[unique(genes), ], exprs_values) %>% as.matrix %>%
        data.table
    if (is.null(long.dt)) {
        stop("There were no values to extract. Check the name for exprs_values.")
    }
    
    long.dt$feature_name <- unique(genes)
    long.dt <- melt(long.dt, id.vars = "feature_name", variable.name = "cell")
    long.dt$value_type <- exprs_values
    return(long.dt)
}
