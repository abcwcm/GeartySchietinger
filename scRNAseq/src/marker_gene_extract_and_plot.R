#' Extract marker gene results from findMarkers into a data.table
#' 
#' @author Friederike Dündar
#' @description Wrestling the output of \code{scran::findMarkers()} into a more palatable 
#' data.table. CAUTION! It is expected that \code{findMarkers()} was run with \code{direction = "up"}
#' 
#' @param sce.object SCE object
#' @param marker_search_result result of \code{scran::findMarkers}, e.g. l.mks$XvsY;
#' typically a list where the names indicate the comparison. This function here
#' assumes that findMarkers was run with \code{direction = "up"}.
#' @param FDR_thresh Default: 0.01
#' @param rank_thresh maximum rank that the extracted genes must have in addition
#' to \code{FDR_thresh}
#' @param logFC_thresh minimum logFC that the extracted genes must have (default: 0)
#' @param comp_name add a name for the comparison that underlies the marker search
#' result, e.g. "DB_clst"; if NULL (default), just the name from the \code{marker_search_result}
#' will be used.
#' 
#' 
#' @import data.table
#' @import magrittr
#' 
#' @return data.table with three columns:
#' \itemize{
#' \item \code{gene_symbol}
#' \item \code{up_in}: in which condition(s) the gene was upregulated
#' \item \code{classify}: whether the gene was found in just one comparison or 
#' multiple ones
#' }
#' 
#' @examples \dontrun{
#' gns_XvsY <- extract_markers(sf, marker_search_result = l.mks$XvsY,
#'                                FDR_thresh = 0.01, rank_thresh = 35)
#' plot_marker_heatmap(sf, gns = gns_XvsY$gene_symbol,
#'                     cells = colnames(sf[,sf$condition != "Z"]),
#'                     exprs_values = "logcounts",show_HKgenes = FALSE,
#'                     gns_as_cellAnnotation = c("Malat1"), 
#'                     n_quant_breaks = 150,
#'                     fontsize_row = 6, 
#'                     main = "Markers when comparing X to Y cells")
#' }
#' 
#' @seealso \code{\link{plot_marker_heatmap}}
#' 
#' 
#' @export
#' 
extract_markers <- function(sce.object, marker_search_result, FDR_thresh = 0.01,
                            rank_thresh = 35, logFC_thresh = 0, 
                            comp_name = NULL){
  
  if(!all(is.numeric(FDR_thresh) & is.numeric(logFC_thresh) & is.numeric(rank_thresh))){
    stop("One of your filters (FDR, rank, logFC) is not numeric.")
  }
  
  ## extract genes from the marker list-----------------------------------------
  gns <- lapply(seq_along(marker_search_result), function(i){
    ## use those genes that pass the defined thresholds
    tmp <- as.matrix(subset(marker_search_result[[i]],
                            FDR <= FDR_thresh & Top <= rank_thresh))
    
    ## at least one logFC value (there may be multiple per gene) must pass the threshold
    tmp2 <- tmp[,grep("logFC",colnames(tmp), value = TRUE), drop=FALSE]
    keep <- apply(tmp2, 1, function(x) any(x > logFC_thresh))
    goi <- rownames(tmp[keep,])
    
    ## turn the results into a data.table
    if(!is.null(goi)){
      out <- data.table(gene_symbol = goi)
      out$cell_type <- paste0(comp_name, names(marker_search_result[i]))
      return(out)
    }else{
      message(paste("No genes passed the thresholds for", names(marker_search_result[i]) ))
    }
  }) %>% rbindlist
  
  if(dim(gns)[1] > 0){
    gns <- gns[,paste(cell_type, collapse = ","), gene_symbol]
    setnames(gns,"V1", "up_in")
    gns[, classify := ifelse(grepl(",", up_in), "shared", "unique")][] # [] is needed for the data.table being printed immediately when the function is used
  
      # > head(gns)
      #gene_symbol      up_in classify
      #1:       Rps27 DB_clstX   unique
      #2:       Rpl41 DB_clstX   unique
      #3:       Rps29 DB_clstX   unique
      #4:      Malat1 DB_clstX   unique
      #5:       Rps28 DB_clstX   unique
      #6:       Pcsk2 DB_clstX   unique
    return(gns)
  }
  
}


##plot pheatmap=====================
#' Plot heatmap of markers
#' 
#' @author Friederike Dündar
#' @description A wrapper function around \code{pheatmap} that helps with tuning
#' the color scheme to what is usually used for scRNA-seq data, selecting specific
#' cell groups, keeping cells sorted according to specific colData entries and so on.
#' CAUTION! Default settings will exclude housekeeping genes as identified by "^mt-" and "^Rp[sl]"!
#' 
#' @param sce.object SCE object
#' @param gns specify which genes to use for the heatmap, default: NULL will lead
#' to all genes of the \code{sce.object} being used
#' @param cells specify which cells to use for the heatmap, default: NULL will lead
#' to all cells of the \code{sce.object} being used
#' @param exprs_values what type of assay data to use. Default: "logcounts"
#' @param show_HKgenes whether housekeeping genes should be shown or not. If 
#' FALSE (default), HK genes will be removed. If TRUE, \emph{only} HK genes will
#' be shown. If "both", the heatmap will contain all of the genes defined via 
#' \code{gns}
#' @param mito_pattern default: "^mt-" (case will be ignored)
#' @param ribo_pattern default: "^Rp[sl]" (case will be ignored)
#' @param gns_as_cellAnnotation specify genes (e.g. "Malat1") that are better displayed
#' as an annotation to the cells rather than as part of the heatmap. Per default
#' setting of \code{exclude_genes}, these genes will be
#' removed from the matrix and won't influence the clustering any more. If you 
#' want to keep them in the matrix, set \code{exclude_genes} to \code{NULL} or just
#' specify those genes that should be excluded.
#' @param exclude_genes Whether or not genes that are used as
#' cell annotation should be excluded from the actual matrix (i.e. they won't
#' influence the clustering etc. and won't be shown in the heatmap). The default
#' setting will exclude all \code{gns_as_cellAnnotation} entries from the matrix.
#' If you want to exclude just specific genes, specify their names here, e.g. "Tox".
#' @param n_quant_breaks if the colors should be assigned to the quantiles, specify
#' the numnber here that will be used to determine the percentiles. \code{quant_breaks = 100}
#' will return 100 percentiles. To turn this off completely, set to NULL (Default).
#' @param add_cell_annotation Add a \code{data.frame} with annotation for the cells, e.g.
#' \code{data.frame(row.names = colnames(sf.tmp),cluster = sf.tmp$beta.lgct_k0.3pct)}.
#' Default: NULL. 
#' @param sort_cells_by If you want to use a pre-specified order of the cells, you can 
#' either indicate an entry of \code{colData} that should be used to sort the cells (e.g. "beta.lgct_k0.3pct"),
#' or you can supply a vector of cell names in the order that you would like them to appear.
#' Note that the clustering by column will be turned off!
#' @param cluster_cells_within_sortGroup set this to TRUE if you want to cluster the columns within each group
#' defined by \code{sort_cells_by}. Default: FALSE. (Will only work if \code{sort_cells_by} is part of \code{colData(sce.object)}).
#' @param annotation_colors pheatmap color list. Default:
#' \code{list(condition = c(`X` = "gray25",`Y` = "limegreen"))}
#' @param scale whether the values should be shown as is ("none"; default) or 
#' as z-scores based on "row" or "column" scaling
#' @param col_palette default: \code{rev(brewer.pal(n = 7, name = "RdYlBu"))},
#' alternative, for example: \code{c("gray88","cornsilk","firebrick3")}
#' @param ... additional parameters for \code{pheatmap}
#' 
#' @examples \dontrun{
#' gns_XvsY <- extract_markers(sf, marker_search_result = l.mks$XvsY,
#'                                FDR_thresh = 0.01, rank_thresh = 35)
#' plot_marker_heatmap(sf, gns = gns_XvsY$gene_symbol,
#'                     cells = colnames(sf[,sf$condition != "DB"]),
#'                     exprs_values = "logcounts", show_HKgenes = FALSE,
#'                     gns_as_cellAnnotation = c("Malat1"), 
#'                     n_quant_breaks = 150,
#'                     fontsize_row = 6, 
#'                     main = "Markers when comparing X to Y cells")
#' 
#' ## a more involved example
#' plot_marker_heatmap(sf, 
#'                     gns = gns_XvsY[classify == "unique"]$gene_symbol, 
#'                     cells = colnames(sf[,as.character(sf$k0.3pct) %in% c("2","6","7","8")]),
#'                     exprs_values = "logcounts", 
#'                     show_HKgenes = FALSE,
#'                     gns_as_cellAnnotation = c("Cst3"), # will be shown on top of the heatmap as a color bar
#'                     n_quant_breaks = 150,
#'                     fontsize_row = 6,
#'                     sort_cells_by = "k0.3pct", # this turns off the clustering of the cols
#'                     add_cell_annotation = data.frame(row.names = colnames(sf),
#'                                                      cluster = sf$k0.3pct),
#'                     annotation_row = data.frame(row.names = gns_XvsY[classify=="unique"]$gene_symbol,
#'                                                 up_in = gns_XvsY[classify=="unique"]$up_in),
#'                     define_anno_cols = list(condition = c(`X` = "gray25",
#'                                                           `Y` = "limegreen"),
#'                                            cluster = c(`2` = "gold",
#'                                                        `8` = "dodgerblue", 
#'                                                        `7` = "violet", 
#'                                                        `6` = "lightcyan1"),
#'                                            Cst3 = scales::brewer_pal("seq", "Purples")(5)[1:4]),
#'                      main = "Markers for different cell clusters\n(unique to one cluster)")
#' }
#' 
#' @seealso \code{\link{extract_markers}}
#' 
#' @import scater
#' @import pheatmap
#' @import RColorBrewer
#' 
#' @export
#' 
plot_marker_heatmap <- function(sce.object, gns = NULL, cells=NULL,
                                exprs_values = "logcounts",
                                show_HKgenes = FALSE,  mito_pattern = "^mt-", ribo_pattern = "^Rp[sl]",
                                gns_as_cellAnnotation = NULL, exclude_genes = gns_as_cellAnnotation,
                                n_quant_breaks = NULL,
                                add_cell_annotation = NULL,
                                define_anno_cols = NULL,
                                sort_cells_by = NULL, cluster_cells_within_sortGroup = FALSE,
                                cluster_cols = TRUE, scale= "none",
                                col_palette = rev(brewer.pal(n = 7, name = "RdYlBu")),
                                ...){
  
  stmp <- sce.object
  
  ## extract genes -----------  
  if(!is.null(gns)){
    if( !any(gns %in% rownames(stmp)) ){
      warning("None of the genes seems to be found in the rownames of the SCE object.")
    }
    stmp <- stmp[unique(gns),]
  }
  
  ## extract cells-------------
  if(!is.null(cells)){
    if( !any(cells %in% colnames(stmp)) ){
      warning("None of the cells seems to be found in the colnames of the SCE object.")
    }
    stmp <- stmp[,unique(cells)]
  }
  
  ## get matrix of logcounts-----
  mat <-  as.matrix(assay(stmp, exprs_values))
  
  ## mitochondrial genes----------
  hk_genes <- unique(c(grep(mito_pattern, rownames(mat), value=TRUE, ignore.case = TRUE),  # mitochondrial genes
                       grep(ribo_pattern, rownames(mat), value=TRUE,ignore.case = TRUE))) # ribosomal genes
  
  ## define genes that should be excluded from the heatmap-------
  exclude <- exclude_genes
 
  ## figure out how to handle HK genes---------------------------
  if(show_HKgenes == FALSE){
    exclude <- unique(c(exclude, hk_genes))
  }
  if(show_HKgenes == TRUE){
    exclude <- unique(c(exclude, 
                        rownames(mat)[which( !(rownames(mat) %in% hk_genes) )]
    ))
  }
  
  ## define annotation for cells-------------------------------------------------
  
  if(!is.null(gns_as_cellAnnotation)){
    cell_anno <- data.frame(t(mat[make.names(gns_as_cellAnnotation), , drop=FALSE]))#,
                            #condition = stmp$condition)
  }else{
    cell_anno <- NULL
    #cell_anno <- data.frame(row.names = colnames(stmp),
    #                        condition = stmp$condition)
  }
  if(!is.null(add_cell_annotation)){
    
    if(!is.null(cell_anno)){
      cell_anno <- merge(cell_anno, add_cell_annotation, by = "row.names")
      row.names(cell_anno) <- cell_anno$Row.names
      cell_anno$Row.names <- NULL
      }else{
        cell_anno <- add_cell_annotation
      }
  }
  
  ## specify the matrix that will actually be used for the heatmap -------------
  if(!is.null(exclude)){
    mat_plt <- mat[ !rownames(mat) %in% exclude,] 
  }else{
    mat_plt <- mat
  }
  
  ## assign heatmap colors based on quantiles-----------------------------------
  if(!is.null(n_quant_breaks)){
    if(scale != "none"){
      mat4breaks <- pheatmap:::scale_mat(mat_plt, scale)
    }else{
      mat4breaks <- mat_plt
    }
    ## define breaks for color assignment
    mat_breaks <- quantile(mat4breaks,
                           probs = seq(0, 1, length.out = n_quant_breaks))
    mat_breaks <- mat_breaks[!duplicated(mat_breaks)]
    cols <- colorRampPalette(col_palette)(length(mat_breaks) - 1)
  }else{
    cols = colorRampPalette(col_palette)(100)
    #if(is.null(breaks)){
    mat_breaks <- NULL
    #}else{mat_breaks <- breaks}
  }
  
  ## sort cells by a specified vector, then turn off the clustering------------
  if( !is.null(sort_cells_by) ){
    
    ## pre-cluster if that's indicated
    if(cluster_cells_within_sortGroup == TRUE){
      if( sort_cells_by %in% names(colData(stmp))){
        
        clstl <- list() ## for each value of sort_cells_by, do a separate clustering of the columns
        for (CLUST in sort(unique(as.character(colData(stmp)[, sort_cells_by ]))) ){
          cells_tmp <- as.character(colData(stmp)[, sort_cells_by]) == CLUST
          
          x <- pheatmap(mat_plt[, cells_tmp], cluster_cols = TRUE, silent = TRUE)
          clstl[[CLUST]] <- x$tree_col$labels[x$tree_col$order] ## saving cell order
        }
        preclustered_cells <- unlist(clstl)
        sort_cells_by <- preclustered_cells
      }else{warning("Couldn't pre-cluster the cells because the variable given to 'sort_cells_by' is not part of colData(sce.object).")}
    }
    
    ## sort the matrix accordingly
    if( all( colnames(stmp) %in% sort_cells_by) ){
      sort_cells_by <- sort_cells_by[which(sort_cells_by %in% colnames(stmp))]
      mat_plt <- mat_plt[, sort_cells_by]
      cluster_cols <- FALSE
      }else if(sort_cells_by %in% names(colData(stmp)) ){
        sort_cells_by <- colData(stmp[,colnames(mat_plt)])[,sort_cells_by]
        mat_plt <- mat_plt[, order(sort_cells_by)]
        cluster_cols <- FALSE
      }
    else{ 
      warning(paste(sort_cells_by, "is not part of colData of the SCE object. Cells will not be pre-sorted."))
      if( is.null(cluster_cols) ){cluster_cols <- TRUE}
    }
  }
  
  ## heatmap! ------------------------------------------------------------------
  pheatmap(mat_plt,
           show_colnames = FALSE,
           breaks = mat_breaks,
           color = cols,
           annotation_col = cell_anno,
           #    annotation_row = data.frame(row.names = goi[classify !="unique" & !(gene_symbol %in% exclude)]$gene_symbol, 
           #                               up_in = goi[classify !="unique" & !(gene_symbol %in% exclude)]$up_in),
           annotation_colors = define_anno_cols,
           cluster_cols = cluster_cols,
           scale = scale,
           ...)
  
}

