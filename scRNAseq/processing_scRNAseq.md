# Processing scRNA-seq
> Friederike Dündar; 9/8/2021

For the initial processing, we mostly wanted to identify barcodes that were 
clearly representing suboptimal droplets based on mostly technical parameters such
as gene content, mitochondrial content, and whether we captured reasonable values 
from the antibody-derived tags (ADT) and the hash-tagged oligos (HTOs).

>The final object `sce_integrated_pLNSamples_filtered.rds` can be downloaded from [here](https://wcm.box.com/shared/static/2wvbdfs3vja2cnlckcqpk20eu1693o8o.rds).

The following code chunk shows how we:

- removed barcodes with very high mitochondrial content and very low gene content,
- removed cells without ADT or HTO information:
    * `ADT`: 7 "functional" **antibodies** characterizing the cell types (*CD44* through *CD39*)
    * `HTO`: varying numbers of **hash tags** (`HT0[1-9]`) marking the mice from which a cell originated
- removed HTO-based doublets (detected with `CiteFuse`).

```r
library(SingleCellExperiment)
library(magrittr);
data_dir <- "Sofia_CITEseq/data/"
smpls <- c("pLN_1","pLN_3","pLN_4") 

## load the CellRanger output into R =======================================
scel <- lapply(smpls, function(x){
    print(x)
    out.sce <- DropletUtils::read10xCounts(
        samples = paste0(data_dir, x,
            "/outs/count/filtered_feature_bc_matrix"),
        sample.names = x,
       version = "auto")
    
    print("Adding alternative experiments")
    # add ADT and HTO info as alternative Experiments to the SCE
    # this will also allow for different numbers of HTO per sample
    print("HTO")
    is.hto <- grepl("^HTO", rownames(out.sce))
    out.sce <- splitAltExps(out.sce, ifelse(is.hto, "HTO", "gene"))
    counts(altExp(out.sce, "HTO")) <- as.matrix(counts(altExp(out.sce, "HTO")))
    
    print("ADT")
    is.ab <- !grepl("^ENSMUS|^HTO", rownames(out.sce))
    out.sce <- splitAltExps(out.sce, ifelse(is.ab, "ADT", "gene"))
    counts(altExp(out.sce, "ADT")) <- as.matrix(counts(altExp(out.sce, "ADT")))
    
    # add colnames
    ncells <- ncol(out.sce)
    colnames(out.sce) <- paste(gsub("_", "", x), 1:ncells, sep = "_")

    return(out.sce)
})

names(scel) <- gsub("_","",smpls)

## QC metrics =================================================================
## To identify low-quality cells, it is easier to look at all samples together.
full_data <- do.call(cbind, lapply(scel, counts))
cell_info <- do.call(rbind, lapply(scel, colData))
gene_info <- rowData(scel[[1]])

sce.all <- SingleCellExperiment(
    list(counts = full_data), 
    rowData = gene_info,
    colData = cell_info, 
    metadata = list(Samples = names(scel))
)

## filtering based on mito content and gene content ---------------------------
is.mito <- grepl("mt-", ignore.case = TRUE, rowData(sce.all)$Symbol)
sce.all <- scuttle::addPerCellQC(sce.all, subsets=list(mitochondrial=is.mito))
sce.all$mito.discard <- isOutlier(sce.all$subsets_mitochondrial_percent, type="higher")
sce.all <- sce.all[, !sce.all$mito.discard]
sce.all$gene.discard <- log10(sce.all$detected) <= 2.5
sce.all <- sce.all[, !sce.all$gene.discard]

## remove cells without ADT detection -----------------------------------------
sce.all$ADT.discard <- unlist(lapply(scel, function(x){
    kp <- colnames(x) %in% colnames(sce.all)
    rm_low_ADTcells(x[, kp])
}))

sce.all <- sce.all[, !sce.all$ADT.discard]

## remove genes without any coverage -------------------------------------
gnszero <- Matrix::rowSums(counts(sce.all)) == 0
sce.all <- sce.all[!gnszero, ]
## Seurat scTransform removes genes that are expressed in less than 5 cells;
## for consistency's sake, let's do that, too
gnslow <- Matrix::rowSums( counts(sce.all) > 0) < 5 # more than 4 cells with min. 1 UMI
sce.all <- sce.all[!gnslow, ]
sce.all$Sample <- gsub("_", "", sce.all$Sample)

## HTOs (mouse labels) =====================================================
#Now that I've removed the low-quality cells, HTO demuxing should be done on the 
#individual samples because the hash tags aren't shared across samples.
#I'll return to separate SCE instances for each sample to 
#ease the separate processing.
## remove filtered out cells and genes (see above)
scel <- lapply(scel, function(x){
    kpc <- colnames(x) %in% colnames(sce.all)
    kpg <- rowData(x)$ID %in% rowData(sce.all)$ID
    out <- x[kpg, kpc]
    rownames(out) <- uniquifyFeatureNames(rowData(out)$ID, rowData(out)$Symbol)
    return(out)
    }
)

library(data.table)
inf <- fread(paste0(data_dir,"samples.txt"))
setnames(inf, names(inf), make.names(names(inf)))
setnames(inf, "Sample.name", "Sample")
inf[, HTO := paste0("HTO",Hashtags)]

## CiteFuse approach for doublet identification
## CiteFuse returns two colData entries:
#    - `doubletClassify_between_label`: the HTO indeces
#    - `doubletClassify_between_class`: doublet/multiplet, negative, Singlet 
scel <- lapply(scel, function(x){
    out.sce <- x
   out.sce$HTO.label <- out.sce$doubletClassify_between_label
    out.sce$HTO.class <- out.sce$doubletClassify_between_class
    out.sce$out.sce$doubletClassify_between_label <- NULL
    out.sce$out.sce$doubletClassify_between_class <- NULL
    out.sce <- scater::runTSNE(out.sce, altexp = "HTO", name = "TSNE_HTO", pca = TRUE)
    #out.sce <- scater::runTSNE(out.sce, altexp = "ADT", name = "TSNE_ADT", pca = TRUE)
    return(out.sce)
}
)
## remove doublets etc.
scel.filt <- lapply(scel, function(x){
    out.sce <- x[, x$HTO.class == "Singlet"] ## only keeping singlets
    out.sce$HTO.class <- NULL
    
    ## get the HTO names
    htos <- data.frame(HTO = rownames(altExp(out.sce, "HTO")))
    
    ## the labels are the indeces of the HTOs
    out.sce$HTO <- htos[(out.sce$HTO.label), , drop = TRUE]
    #out.sce$doubletClassify_between_label <- NULL
    
    cd <- colData(out.sce)
    cdd <- merge(
        data.frame(cd, cell=rownames(cd)),
        data.frame(inf[Sample == as.character(unique(out.sce$Sample)), -c("Sample","Hashtags")]),
        by.x = "HTO", by.y="HTO", all = TRUE ) %>% DataFrame
    rownames(cdd) <- cdd$cell
    colData(out.sce) <- cdd[colnames(out.sce), ]
    out.sce$Mouse <- paste0("m.", out.sce$Mouse.numbers)
    
    return(out.sce)
})

## Normalizing ADT values ======================================================
scel.filt <- lapply(names(scel.filt), function(x){
    
    print(x)
    out.sce <- scel.filt[[x]]
    
    # add sizefactor for ADT ------------------------------------------------
    Exp <- "ADT"
    
    baseline <- DropletUtils::inferAmbience(counts(altExp(out.sce, Exp)))
    sizeFactors(altExp(out.sce, Exp)) <- scuttle::medianSizeFactors(altExp(out.sce, Exp), reference=baseline)
    
    # lognorm --------------------------------------------------------------
    print("Log-normalization")
    ## need to calculate the size factors for the RNA-seq
    quiclu <- scran::quickCluster(out.sce)
    out.sce <- scran::computeSumFactors(out.sce, cluster=quiclu, min.mean=0.1)
    
    ## logNorm will transform both expression and altExps
    out.sce <- scuttle::logNormCounts(out.sce, use.altexps="ADT") 
    
    return(out.sce)
})
names(scel.filt) <- lapply(scel.filt, function(x) unique(x$Sample)) %>% unlist

##!saveRDS(scel.filt, file = paste0(data_dir,"sce-list_allCells_filtered_noHTODoubs.rds"))
```

Once the filtering based on technical properties was done, we integrated all
samples, performed clustering, dimensionality reduction, marker gene detection,
cell cycle score calculations and so on.
During that first round, we found that specific subpopulations and small clusters of cells
represented droplets that we weren't interested in for subsequent analyses, including CD4 T cells, B cells, and likely doublets.
We're not showing the iterative process of determining those cells here, just the
final code that was used (a) to remove the barcodes representing CD4 T cells etc.
and (b) to integrate the remaining, cleaned cells into one composite SingleCellExperiment object.

As a reminder, these were the types of cells that were removed in the following chunk of code:

- removed CD4,
- removed B cells,
- removed two small clusters of cells with high G1 scores and elevated numbers of genes compared to all other cells

```r
library(magrittr)
library(batchelor)
library(scran)
library(scater)
library(patchwork)
fnmerged <- "sceMultiout_integrated_pLNSamples.rds"

# 1. List of SCE
all.sce <- readRDS("sce-list_allCells_filtered_noHTODoubs.rds")
rmcd4 <- read.table("rm_cd4.txt") ## read in the barcodes of cells to be removed
rmb <- read.table("rm_bcells.txt")## read in the barcodes of cells to be removed
all.sce <- lapply(all.sce, function(x){
    rmcells <- unique(c(rmcd4$V1, rmb$V1))
    rmcells <- rmcells[rmcells %in% colnames(x)]
    kpcells <- colnames(x)[!colnames(x) %in% rmcells]
    return(x[, kpcells])
})

## filter genes
rd_qc <- lapply(all.sce,perFeatureQCMetrics)
for(x in seq_along(all.sce)){rowData(all.sce[[x]])$qc.mean <- rd_qc[[x]]$mean}
all.sce <- lapply(all.sce, function(x){ x[rowData(x)$qc.mean > 0.001,]})
## combine
universe <- Reduce(intersect, lapply(all.sce, rownames))
all.sce2 <- lapply(all.sce, "[", i=universe,)
normed.sce <- do.call(multiBatchNorm, all.sce2) # returns a list
# Identifying a set of HVGs using stats from all batches
all.dec <- lapply(all.sce2, modelGeneVar)
combined.dec <- do.call(combineVar, all.dec)
combined.hvg <- getTopHVGs(combined.dec, n=2500)
# Merge with MNN ----------------------------
## prep
combined <- noCorrect(normed.sce)
assayNames(combined) <- "logcounts"
combined$Sample <- combined$batch
combined$Tissue <- gsub("_.*","",combined$batch) 
combined$Sample.ID <- gsub(".*_","s",combined$batch)
combined$batch <- NULL
set.seed(1010100)
## progressively merge cells from each sample in each batch until all cells 
## are mapped onto a common coordinate space
multiout <- fastMNN(combined, batch=combined$Sample, subset.row=combined.hvg)
# Renaming metadata fields for easier communication later.
multiout$Sample <- multiout$batch
multiout$Tissue <- gsub("_.*","",multiout$batch)
multiout$Sample.ID <- gsub(".*_","s",multiout$batch)
# clustering ----------------------------
g <- buildSNNGraph(multiout, use.dimred="corrected", k = 20)
clusters <- igraph::cluster_louvain(g)
colLabels(multiout) <- factor(clusters$membership)
table(colLabels(multiout), multiout$batch)
multiout$cluster.20 <- colLabels(multiout)
saveRDS(multiout,file = fnmerged)
## UMAP------ ----------------------------
set.seed(10101010)
multiout <- runUMAP(multiout, dimred="corrected")
saveRDS(multiout,file = fnmerged)

# generate composite file for the combined file of all shared genes ============
comb.mat <- lapply(all.sce, function(x) counts(x)) %>% do.call(cbind, .)
colnames(comb.mat) <- unlist(lapply(all.sce, function(x) colnames(x)))

### rowData
rd <- rowData(all.sce[[1]])[, c("ID","Symbol")]
rd <- rd[rownames(comb.mat),]

### colData
cd <- lapply(all.sce, function(x) colData(x)[, c("Sample","Barcode","cell","Mouse","sizeFactor","Diabetic.status","Age..wk.")]) %>% 
    do.call(rbind, .)
#### get cluster information from multiout
cd2 <- as.data.frame(colData(multiout))
cd2$cell <- rownames(cd2)
cd2 <- cd2[, c("cell","cluster.20")]
cd <- merge(cd, cd2, by = "cell")

#### ADTs in future colData
adts <- lapply(all.sce, function(x){
    assay(altExp(x, "ADT"), "logcounts")
}) %>% do.call(cbind, .) %>% t %>% as.data.frame
adts$cell <- rownames(adts)
cd <- merge(cd, adts, by = "cell")

rownames(cd) <- cd$cell
cd <- cd[colnames(comb.mat),]

sce.pLN <- SingleCellExperiment(assays = list(counts = comb.mat), colData = cd, rowData = rd)

# add redDims from the merged data set
rdu <- reducedDim(multiout, "UMAP") 
reducedDim(sce.pLN, "UMAP") <- rdu[colnames(sce.pLN),]
reducedDim(sce.pLN, "PCA_corr") <- reducedDim(multiout, "corrected")

# add log-counts
qckclst <- quickCluster(sce.pLN, method = "igraph", min.mean = 0.1)
sce.pLN <- computeSumFactors(sce.pLN, min.mean=0.1, cluster = qckclst)
sce.pLN <- scater::logNormCounts(sce.pLN)

# more colData details
sce.pLN$Tissue <- gsub("_.*","",sce.pLN$Sample) 
sce.pLN$Sample.ID <- gsub(".*_","s",sce.pLN$Sample)

## remove high G1-Phase and high-gene-content clusters (1, 3) -----------------
sln <- sce.pLN[, as.character(sce.pLN$cluster.20) %in% c("2","4","5")]
rd <- reducedDim(sceln, "UMAP") %>% as.data.frame
rm <- subset(rd, (V1 < 0 & V1 > -2) | V2 >8) %>% rownames
sln <- sln[, !(colnames(sln) %in% rm)]
colData(sln) <-  colData(sln)[, 1:26]
sln$population <- ifelse(as.character(sln$cluster.20) == "2", "AIC", 
    ifelse(as.character(sln$cluster.20) == "4", "AMC", "transition"))

saveRDS(sln, file = "sce_integrated_pLNSamples_filtered.rds") 
```

The final object `sce_integrated_pLNSamples_filtered.rds` can be downloaded from [here](https://wcm.box.com/shared/static/2wvbdfs3vja2cnlckcqpk20eu1693o8o.rds).
