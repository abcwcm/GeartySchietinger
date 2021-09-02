#' Download and locally cache an RDA or RDS file that's stored on a website
#' 
#' @author Friederike Dündar
#' @description Using the library \code{BiocFileCache}, this function takes as 
#' input the direct link to a RDA (or RDS) file, downloads it and stashes it in a 
#' local directory from where it can be retrieved (using this exact function) in
#' a later session. I.e. if the data has already been downloaded, it won't be down-
#' loaded again.
#' 
#' @details Based on instructions from the BiocFileCache vignette: \url{https://bioconductor.org/packages/release/bioc/vignettes/BiocFileCache/inst/doc/BiocFileCache.html#local-cache-of-an-internet-resource}. This is a slightly more versatile form of
#' the function \code{load_RDSdata_from_Box()} as it also allows to load RDA files
#' directly into the global environment by setting the option \code{load_rds} to FALSE.
#' 
#' @param shared_link The URL for the file as a string. NOTE: This should be the
#' \emph{direct} link, i.e. it should end with ".rda" or ".rdata", "rds" or the like
#' @param data_name A customized name that may be used to more easily refer to 
#' the cached data (\code{rname} parameter of \code{BiocFileCache} functions),
#' e.g. "SCE2019" -- any kind of string will work, but it's important
#' to make a note of this name for reloading the data next time, otherwise the 
#' data set will be downloaded again if there's a mismatch with the name.
#' If set to NULL (default), the \code{rname} in the cache will correspond to the
#' url indicated via \code{shared_link}.
#' @param cache_path If wanted, this allows to specify a path where the data
#' should be cached. if NULL (default), the default setting of \code{BiocFileCache()}
#' will be used.
#' @param check_for_update Indicate whether the data should be downloaded again.
#' 
#' @return The final command here is \code{readRDS()}, so you should be assigning
#' the result of this function to an object name in your environment.
#'
#' @examples 
#' 
#' ## load RDS requires the user to define an object name (here: testingsce)
#' direct_link <- "https://wcm.box.com/shared/static/4887keeu77eskmdotve1sm7n6u7udckg.rds"
#' testingsce <- load_data_from_Box(shared_link = direct_link, data_name = "testRDS")
#' 
#' ## loading RDA
#' sl <- "https://wcm.box.com/shared/static/mmf98o464n94fdh78x14sbwo7ts56ukq.rda"
#' load_data_from_Box(sl, load_rds = FALSE) # should add tiny_sce to global env.
#' ls()
#' 
#' @export
#'
load_data_from_Box <- function(shared_link = "https://wcm.box.com/shared/static/mmf98o464n94fdh78x14sbwo7ts56ukq.rda",
  data_name = NULL, cache_path = NULL, check_for_update = FALSE, load_rds = TRUE){
  
  if(grepl("\\.rds$", shared_link, ignore.case = TRUE)){
    load_rds <- TRUE
  }
  
  suppressPackageStartupMessages({ library(BiocFileCache) })
  load_new <- FALSE
  
  if(is.null(data_name)){
    data_name <- shared_link
  }
  
  ## load the cache
  if(is.null(cache_path)){
    bfc <- BiocFileCache(ask=F)
  }else{
    bfc <- BiocFileCache(cache_path, ask=FALSE)
  }
  
  ## check if url is being tracked
  bfc_res <- bfcquery(bfc, shared_link)
  
  if(bfccount(bfc_res) == 0L){
    load_new <- TRUE
  }else{
    ## if it is in cache, get path to load
    rid <- bfc_res[bfc_res$fpath == shared_link & bfc_res$rname == data_name,]$rid 
    
    ## double-check we're looking at what we want
    if(length(rid) > 0){
      if(length(rid) > 1){
        message(paste0("Multiple entries found for", data_name, shared_link,". We'll use only the first one."))
        rid <- rid[1]
      }
      
      ans <- bfcrpath(bfc, rids = rid)
      
      if(check_for_update){
        ## check to see if the resource needs to be updated
        check <- bfcneedsupdate(bfc, rids = rid)
        if( is.na(check)) upd <- TRUE ## 'check' can be NA if it cannot be determined, choose how to handle
        if(upd) ans <- bfcdownload(bfc, rid = rid)
      }
    }else{
      load_new <- TRUE
    }
  }
  
  ## if it wasn't in the cache, add 
  if(load_new){
    ans <- bfcadd(bfc, rname = data_name, fpath = shared_link)
  }
  
  message(paste("Cached data here:", ans))
  
  if(load_rds){
    readRDS(ans)
  }else{
    load(ans, .GlobalEnv)
  }
  
}


#' Download and locally cache an RDS file that's stored on a website
#' 
#' @author Friederike Dündar
#' @description Using the library \code{BiocFileCache}, this function takes as 
#' input the direct link to a RDS (!) file, downloads it and stashes it in a 
#' local directory from where it can be retrieved (using this exact function) in
#' a later session. I.e. if the data has already been downloaded, it won't be down-
#' loaded again.
#' 
#' @details Based on instructions from the BiocFileCache vignette: \url{https://bioconductor.org/packages/release/bioc/vignettes/BiocFileCache/inst/doc/BiocFileCache.html#local-cache-of-an-internet-resource}
#' 
#' @param shared_link The URL for the file as a string. NOTE: This should be the
#' \emph{direct} link, i.e. it should end with ".rds" (!)
#' @param data_name A customized name that may be used to more easily refer to 
#' the cached data (\code{rname} parameter of \code{BiocFileCache} functions),
#' e.g. "SCE2019" -- any kind of string will work, but it's important
#' to make a note of this name for reloading the data next time, otherwise the 
#' data set will be downloaded again if there's a mismatch with the name.
#' If set to NULL (default), the \code{rname} in the cache will correspond to the
#' url indicated via \code{shared_link}.
#' @param cache_path If wanted, this allows to specify a path where the data
#' should be cached. if NULL (default), the default setting of \code{BiocFileCache()}
#' will be used.
#' @param check_for_update Indicate whether the data should be downloaded again.
#' 
#' @return The final command here is \code{readRDS()}, so you should be assigning
#' the result of this function to an object name in your environment.
#'
#' @examples 
#' direct_link <- "https://wcm.box.com/shared/static/4887keeu77eskmdotve1sm7n6u7udckg.rds"
#' testingsce <- load_RDSdata_from_Box(shared_link = direct_link, data_name = "testRDS")
#'
#' @export
#'
load_RDSdata_from_Box <- function(shared_link, data_name = NULL, 
    cache_path = NULL, check_for_update = FALSE){
  
  suppressPackageStartupMessages({ library(BiocFileCache) })
  load_new <- FALSE
  
  if(is.null(data_name)){
    data_name <- shared_link
  }
  
  ## load the cache
  if(is.null(cache_path)){
    bfc <- BiocFileCache(ask=F)
  }else{
    bfc <- BiocFileCache(cache_path, ask=FALSE)
  }
  
  ## check if url is being tracked
  bfc_res <- bfcquery(bfc, shared_link)
  
  if(bfccount(bfc_res) == 0L){
    load_new <- TRUE
  }else{
    ## if it is in cache, get path to load
    rid <- bfc_res[bfc_res$fpath == shared_link & bfc_res$rname == data_name,]$rid 
    
    ## double-check we're looking at what we want
    if(length(rid) > 0){
       if(length(rid) > 1){
         message(paste0("Multiple entries found for", data_name, shared_link,". We'll use only the first one."))
         rid <- rid[1]
       }
      
      ans <- bfcrpath(bfc, rids = rid)
      
      if(check_for_update){
        ## check to see if the resource needs to be updated
        check <- bfcneedsupdate(bfc, rids = rid)
        if( is.na(check)) upd <- TRUE ## 'check' can be NA if it cannot be determined, choose how to handle
        if(upd) ans <- bfcdownload(bfc, rid = rid)
      }
    }else{
      load_new <- TRUE
    }
  }
  
  ## if it wasn't in the cache, add 
  if(load_new){
    ans <- bfcadd(bfc, rname = data_name, fpath = shared_link)
  }
  
  message(paste("Cached data here:", ans))
  readRDS(ans)
}
