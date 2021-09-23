# Processing single-cell TCR sequencing data
> Friederike Dündar, 9/8/2021

In order to determine whether clones from the pancreatic lymph node (pLN) ended up feeding the supply of beta-cell-destroying cells in the pancreas, we performed single-cell TCR-sequencing for CD8 T cells from the pLN as well as from the pancreas of the same mouse. 
(For detail of the cell sorting, see the publication)

The code here is meant to illustrate the steps we took to generate the R object
containing the VDJ data that can be downloaded from the Box.
The raw data can be obtained from GEO.

In brief, the processing entails:

- reading in CellRanger output,
- adding the mouse IDs, which are inferred from the scRNA-seq data results (via HTOs),
- focussing on clonotypes that are found in pLN and pancreas samples from the same mouse ("paired"),
- adding a clonotype ID that is (a) unique, (b) captures the same clonotype across the different samples from the same mouse, and (c) allows for the possibility that instances of cells where we captured just a single TRAx or a single TRBy (instead of the proper (TRAx;TRBy pair) are highly likely to belong to the same clonotype 

The final object can be downloaded from the [Box](https://wcm.box.com/shared/static/5wj8h9qssai2u0ya2sre2tsudahvvoni.rds).

```
library(SingleCellExperiment); library(magrittr); library(data.table)
data_dir <- "Sofia_CITEseq/data/"

## read in SCE to have the info about paired samples (stored in the form of
## the Mouse IDs in the colData
scel <- readRDS(file = paste0(data_dir, "sce-list_allCells_filtered_noHTODoubs.rds"))
paired_samples <- lapply(scel, function(x) data.table(Mouse = unique(x$Mouse))) %>%
    rbindlist(., idcol = "Sample") %>% 
    .[, .N, Mouse] %>% .[N>1] %>% .$Mouse

## read in the TCR info from CellRanger ======================================
smps <- names(scel)
tcrfile <- "outs/vdj_t/filtered_contig_annotations.csv"
tcr <- lapply(smps, function(x) fread(paste0(data_dir, x,"/",tcrfile)))
names(tcr) <- smps

## combine the TCR info from the individual samples
tcr.all <-  rbindlist(tcr, idcol="Sample")
setnames(tcr.all, "barcode", "Barcode")

## extract those cells that came from paired samples
cddt <-  lapply(scel, function(x){
    as.data.table(colData(x[, x$Mouse %in% paired_samples]))
    }) %>% rbindlist
tcr <- tcr.all[cddt, on =c( "Barcode","Sample")] %>% .[Mouse %in% paired_samples]
##!save(tcr, file = paste0(data_dir, "tcr_contigs_pairedSamples.rda"))

## consensus clonotype IDs from CellRanger ==================================
 consens_annot <- fread("../data/panc_1/outs/vdj_t/consensus_annotations.csv")

tcr[, Tissue := gsub("_.*","",Sample)]
tcr[ , Sample := gsub(".*_", "s", Sample)]

## Since `raw_clonotype_id == "clonotype1"` in Sample X is not the same as
## `raw_clonotype_id == "clonotype1"` in Sample Y, I need to find a **unified 
## clonotype ID per paired sample**
## That's what I now use the `djvdjs` import function for
## see https://github.com/friedue/SCEdjvdj
source("SCEdjvdj/R/import_vdj.R")
source("SCEdjvdj/R/utils_plots.R")
source("SCEdjvdj/R/plots_SCE.R")
source("SCEdjvdj/R/calc_abundance.R")
source("SCEdjvdj/R/utils_SCE.R") # .add_colData
vdj.dirs <- paste0(data_dir, c("panc_1","panc_3","panc_4", "pLN_1", "pLN_3","pLN_4"), "/")
names(vdj.dirs) <- c("panc_1","panc_3","panc_4", "pLN_1", "pLN_3","pLN_4")

## get mouse information
mice <- fread("mouse_per_cell.csv")
# remove cells without TCR info
#mice <- mice[Barcode %in% tcr.df$Barcode]
mice[, Tissue:= gsub("[0-9]_.*", "", cell)]
mice[, Sample:= gsub("(panc|pLN)([0-9])_.*", "\\1_\\2", cell)]

## determine which mice were profiled in both tissues
paired_samples <- mice[, c("Mouse","Tissue")] %>% unique %>% .[, .N, by=c("Mouse")] %>% .[N>1] %>% .$Mouse
mice[ , pairedExp := ifelse(Mouse %in% paired_samples, TRUE, FALSE)]

tcr.df <- import_vdj(vdj_dir = vdj.dirs)
tcr.df$Barcode <- gsub(".*_[134]_","", rownames(tcr.df))
tcr.df$Sample <- gsub("(.*_[134])_.*","\\1", rownames(tcr.df))

# remove cells without mouse info
tcr.df <- subset(tcr.df, Barcode %in% mice$Barcode)

tcr.cells <- as.data.table(tcr.df) %>% mice[., on =c("Barcode","Sample")]
tcr.cells <- tcr.cells[!is.na(Tissue)]

cl.ids <- tcr.cells[, c("cdr3_nt","Mouse"), with=FALSE] %>% unique
cl.ids[, clono.ID := paste0("cid",1:nrow(cl.ids))]
tcr.cells <- cl.ids[tcr.cells, on = c("cdr3_nt","Mouse")]

## Assign new clonotype IDs that collapse single-TRAs and single-TRBs
## with cells that show the same TRA or TRB in a proper pair
## because it is highly likely that we've just missed one of the chains
## for some cells [also, I'll focus on single-TRAs/TRBs that aren't found in
## multiple pairings ("monogamous") and are also not present in "wrong" pairs 
## such as TRA;TRA;TRB 
## This strategy isn't 1000% clean and possibly too strict, but it should give
## us a less restricted set than completely ignoring the fact that the single 
## TRA/TRB instances are probably represented in proper TRA;TRB pairs, too
traStatus <- lapply(unique(tcr.cells$Mouse), function(x){
    tmp2 <- unique(tcr.cells[Mouse==x & chains %in% c("TRA","TRB","TRA;TRB") & TRA != "noProperPair",
        c("Mouse","chains","cdr3","cdr3_nt","TRA","TRB")])
    tmp <- tcr.cells[Mouse==x, c("Mouse","chains","cdr3_nt")] %>% unique
    data.table(
    TRA = tmp2$TRA,
    Mouse = tmp2$Mouse,
    inFringeCases.TRA = unlist(lapply(tmp2$TRA, function(y) any(grepl(y, unique(tmp[!chains %in% c("TRA","TRB","TRA;TRB")]$cdr3_nt)))))) %>% unique
}) %>% rbindlist

trbStatus <- lapply(unique(tcr.cells$Mouse), function(x){
  #  trbtmp <- trbPairCount[Mouse==x]
#    tmp <- tcr.cells[Mouse==x]
     trbtmp <- unique(tcr.cells[Mouse==x & chains %in% c("TRA","TRB","TRA;TRB") & TRB != "noProperPair",
        c("Mouse","chains","cdr3","cdr3_nt","TRA","TRB")])
    tmp <- tcr.cells[Mouse==x, c("Mouse","chains","cdr3_nt")] %>% unique
    data.table(
    TRB = trbtmp$TRB,
    Mouse = trbtmp$Mouse,
    inFringeCases.TRB = unlist(lapply(trbtmp$TRB, function(x) any(grepl(x, unique(tmp[!chains %in% c("TRA","TRB","TRA;TRB")]$cdr3_nt)))))) %>% unique
}) %>% rbindlist

traPairCount <- traPairCount[traStatus, on = c("Mouse","TRA")]
trbPairCount <- trbPairCount[trbStatus, on = c("Mouse","TRB")]

traPairCount[, monogamous.tra := ifelse(!is.na(n_pairs) & n_pairs > 1, FALSE, !inFringeCases.TRA)]
trbPairCount[, monogamous.trb := ifelse(!is.na(n_pairs) & n_pairs > 1, FALSE, !inFringeCases.TRB)]

## add info to original dt; "none" will be used for non-TRA;TRB pairs
tcr.cells <- traPairCount[, -c("n_pairs","inFringeCases.TRA"), with=FALSE] %>% .[tcr.cells, on = c("Mouse","TRA")]
tcr.cells <- trbPairCount[, -c("n_pairs","inFringeCases.TRB"), with=FALSE] %>% .[tcr.cells, on = c("Mouse","TRB")]

cid <- tcr.cells[chains %in% c("TRA","TRB","TRA;TRB"), 
    c("TRA","TRB","monogamous.tra","monogamous.trb","chains","Mouse","cdr3_nt","clono.ID"), with=FALSE] %>% unique
cidProper <- cid[ chains=="TRA;TRB"]
cidProper[, new_cid := paste0("clone",1:nrow(cidProper))]

## monogamous, single TRAs and TRBs get the new and matching cids
cidTRA <- cid[
    chains == "TRA" & !is.na(monogamous.tra) & monogamous.tra == TRUE,
    c("chains","TRA","Mouse","cdr3_nt","clono.ID")] %>% 
    cidProper[., on = c("TRA","Mouse")] %>% .[, c("i.chains","i.cdr3_nt","i.clono.ID","Mouse","new_cid"),with=FALSE] %>% .[!is.na(new_cid)]
setnames(cidTRA, names(cidTRA), gsub("^i\\.", "", names(cidTRA)))
cidTRB <- cid[
    chains == "TRB" & !is.na(monogamous.trb) & monogamous.trb == TRUE,
    c("chains","TRB","Mouse","cdr3_nt","clono.ID")] %>% 
    cidProper[., on = c("TRB","Mouse")] %>% .[, c("i.chains","i.cdr3_nt","i.clono.ID","Mouse","new_cid"),with=FALSE] %>% .[!is.na(new_cid)]
setnames(cidTRB, names(cidTRB), gsub("^i\\.", "", names(cidTRB)))

## the remaining single-chain instances get their own new cid
cidsSingle <- cid[!(clono.ID %in% unique(c(cidTRB$clono.ID, cidTRA$clono.ID))) & chains %in% c("TRA","TRB"),
    c("chains","cdr3_nt","clono.ID","Mouse")]
cidsSingle[, new_cid := paste0("cloneS.", 1:nrow(cidsSingle))]

## combining all single-chain instances with their new cids
cidsSingle <- rbind(cidsSingle, cidTRA, cidTRB)

newCids <- rbind(cidProper[, c("chains","Mouse","cdr3_nt","new_cid","clono.ID")], cidsSingle)

## add back into tcr.cells
tcr.cells <- newCids[tcr.cells, on = c("chains","Mouse","cdr3_nt","clono.ID")]
##!saveRDS(tcr.cells[pairedExp==TRUE], file = paste0(data_dir,"tcr_perCell_pairedSamples_newCID_2021-05.rds"))
```
