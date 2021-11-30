[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5525867.svg)](https://doi.org/10.5281/zenodo.5525867)


# Bioinformatic methods for Gearty et al. (2021)

>Gearty, S.V., Dündar, F., Zumbo, P. et al. An autoimmune stem-like CD8 T cell population drives type 1 diabetes. Nature (2021). https://doi.org/10.1038/s41586-021-04248-x

In Gearty et al. (2021) we describe CD8 T cell populations that contribute to the continuous destruction of pancreatic beta cells leading to type 1 diabetes.
The data described here entails bulk RNA-seq samples of beta-cell-specific CD8 T cells from murine pancreatic lymph nodes (pLN) and pancreas as well as single-cell RNA-seq data coupled with antibody-derived tags (CITE-seq) and single-cell TCR sequencing.

* [Bulk RNA-seq](#bulk-rna-seq)
* [Single-cell VDJ and CITE-seq data processing](#scvdj-and-scrna-seq-data-processing)
* [Integration of public scRNA-seq data sets](#integration-with-public-scRNA-seq-data-sets)
* [References](#references)
* [Software versions](#package-versions)
* [Data for download](#data-for-download)

All scripts were written by Paul Zumbo and Friederike Dündar.
Samples were prepared by Sofia Gearty and the sequencing facilities at MSKCC and Weill Cornell Medicine.
Don't hesitate to [get in touch](https://abc.med.cornell.edu/) with questions related to the code.

![](WCM_MB_LOGO_HZSS1L_CLR_RGB.png)

## Bulk RNA-seq

### Sample Preparation 

For pLN and pancreas samples, NOD mice aged 14-20 weeks were used.
Samples were isolated as follows: samples were stained with NRP-V7 tetramer, Live/dead Zombie dye, and antibodies against CD8alpha, CD45.1, CD44, and Ly108. For the pancreas, samples were enriched for islets using methods adapted from published protocols.
1,000-2,000 cells were sorted into FCS, washed with PBS, and resuspended in TRIzol LS reagent and stored at -80C.
For pLN samples, 500-3,000 cells were sorted directly into TRIzol LS reagent.
For naive pLN controls, NOD mice aged 6 weeks were used and cells were stained with CD8alpha, NRP-V7 tetramer, CD45.1, Live/dead Zombie dye, and CD44. 5,000-10,000 cells were sorted directly intro TRIzol LS reagent.
RNA extraction was extracted with TRIzol LS.
After RiboGreen quantification and quality control by Agilent BioAnalyzer, RNA underwent amplification using the SMART-Seq v4 Ultra Low Input RNA Kit, with 12 cycles of amplification.
Subsequently, 0.7-3ng of amplified cDNA was used to prepare libraries with the KAPA Hyper Prep Kit using 8 cycles of PCR. 
Samples were barcoded and run on a HiSeq 4000 in a 50bp/50bp paired end run, using the HiSeq 3000/4000 SBS Kit (Illumina). 

### Processing 

Sequencing reads were mapped with `STAR v2.6.0c` with default parameters to the mouse reference genome89 (GRCm38.p6).
Fragments per gene were counted with `featureCounts v1.6.2` with respect to Gencode vM17 comprehensive gene annotations. 
Differentially expressed genes between pairwise comparisons were identified by Wald tests using `DESeq2 v1.26.092`, and only Benjamini–Hochberg corrected P-values <0.10 were considered statistically significant. 
Likelihood ratio tests (LRT) were used to determine genes which vary between the four conditions in any way (adj. P < 0.01). Genes identified by the LRT method were clustered based on their changes in expression across the four conditions with the degPatterns function from the `DEGreport v1.22.0` R package. 
Base-2 log-transformed counts per million (CPM) values were used for heatmap plots of bulk RNA-seq data, which were centered and scaled by row. 
For gene set enrichment analysis, the `fgseaMultilevel` function from the R-package `fgsea v1.12.0` was run in pre-ranked mode (genes were ranked based on `DESeq2`'s Wald statistic).
Gene ontology over representation analysis was performed for biological processes ontology domain with `clusterProfiler v3.14.3`.


## scVDJ and scRNA-seq data processing

In principle, we followed the ideas and workflows discussed by Amezquita et al. in the excellent online book ["Orchestrating single-cell analysis"](https://osca.bioconductor.org/) describing single-cell analyses carried out with packages of the Bioconductor environment.

### Sample prep

Female NOD mice aged 13-20 weeks were used. Samples were prepared as described above for RNA-seq, except before pooling samples were incubated with hashtags (Biolegend Total-Seq C0301-10) to enable multiplexing.
Samples were stained with NRP-V7 tetramer, Live/dead Zombie dye, and antibodies against CD8alpha and CD45.1; in addition, Biolegend feature barcoding antibodies against CD44, CD62L, CD127, CD73, PD1, CD38, and CD39 were added for CITE-seq analysis. 
2,000-16,000 cells were sorted into PBS+0.04% BSA, immediately analyzed for viability, and processed. 
For details of the library preparation, see [`wetLabPrep_details.md`](https://github.com/abcwcm/GeartySchietinger/blob/master/scRNAseq/wetLabPrep_details.md).

### Processing and analysis

VDJ and single-cell RNA-seq including hash-tagged oligos (HTO; used to label cells coming from the same mouse), and antibody-derived tags (ADT; used to label cells based on the expression of selected surface markers) data were processed following the recommendations by [Amezquita et al](https://bioconductor.org/books/release/OSCA) (version: 2020-11-13). In brief, raw read files were processed and aligned using the `CellRanger` pipeline (`cellranger-5.0.0` with `refdata-gex-mm10-2020-A` and `vdj_GRCm38_alts_ensembl-5.0.0`; for annotation supplied by 10X Genomics: <https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest>).
Read counts representing gene expression were used for removing barcodes representing empty droplets or droplets with high amounts of mitochondrially encoded gene products: using the `isOutlier()` function of the scater package, all cells whose fractions of mitochondrial reads were higher than 3x median absolute deviation (MAD) were flagged and removed as well as cells with fewer than 10^2.5 detected genes.
Following the first round of filtering, HTO demultiplexing and doublet detection were done with [`CiteFuse`](http://www.bioconductor.org/packages/release/bioc/html/CiteFuse.html), removing more droplets and identifying the donor mouse for each remaining droplet. After removal of low-quality droplets and barcodes associated with HTO-inferred doublets, ADT values were normalized using
`DropletUtils::inferAmbience()` and `scuttle::medianSizeFactors()`.
Finally, droplets for which less than 50% of the ADT were detected compared to other droplets were removed and only genes that were expressed in at least 5 cells per sample were kept for downstream analyses.

VDJ sequencing data was imported into R using the `import_vdj()` function of the [`djvdj` package](https://github.com/rnabioco/djvdj), which was adjusted to work with `SingleCellExperiment` objects (<https://github.com/friedue/SCEdjvdj>).

Size-factor normalized logcounts were obtained via `scran::computeSumFactors()` and `scater::logNormCounts()` (Lun2016, McCarthy2017).
Cells from the three different technical replicates were integrated with `batchelor::fastMNN()` using the top 2500 most variable genes with min. mean normalized expression of 0.001 (Haghverdi2018).
Clustering was performed with `igraph::cluster_louvain()`, dimensionality reductions were done with `scater` functions (`runUMAP()`) using the batch-corrected values and `destiny::diffusionMap()` (Lun2016clustering, Angerer2015).
Marker genes were determined with `scran::findMarkers()`.
The `TSCAN package (v.1.28.0)` was used to calculate pseudotime values and trajectories as well as genes associated with the pseudotime gradients.
GO term enrichments were calculated with the `clusterProfiler` package's functions `compareCluster()` and `enrichGO()` after excluding ribosomal genes from the gene lists of interest.
All plots were generated using `ggplot2` packages and the `pheatmap` package for heatmaps.

For details, see [processing_scRNAseq.md](scRNAseq/processing_scRNAseq.md), [processingVDJseq.md](scRNAseq/processingVDJseq.md), [figures_scRNAseq.Rmd](scRNAseq/figures_scRNAseq.Rmd) and [figures_VDJseq.Rmd](scRNAseeq/figures_VDJseq.Rmd).

## Integration with public scRNA-seq data sets

Gene count matrices for day 7 CD8 T cells from acute and chronic infection were obtained from GEO (GSE119940); using cell labels provided by Chen Yao we extracted the data for barcodes corresponding to memory precursor and memory-like cells as described in Yao et al. (2019).
scRNA-seq data for Tcm were downloaded from GEO (day 129 following acute infection: GSM3732587) where they had been uploaded to by Schauder et al.
scRNA-seq of Tcms, MPECS, and Tpex were subsequently integrated with the pLN data set using `batchelor::multiBatchNorm()` and `batchelor::fastMNN()` as described above.
For global comparisons of the different populations of T cells, we created pseudo-bulk samples by aggregating the read counts per gene across cells of the same population.
These were then cpm-normalized via `edgeR::calcNormFactors` (Robinson 2010) and subsequently analyzed and visualized via PCA and hierarchical clustering using base R functions as well as `pcaExplorer::hi_loadings()` and the `dendextend` package (Marini 2019, Galili 2015).
All other analyses were done with the same principles and packages as described above.

For details see [Schauder2021.Rmd](scRNAseq/Schauder2021.Rmd), [Yao2019.Rmd](scRNAseq/Yao2019.Rmd), [figures_public_scRNAseq.Rmd](scRNAseq/figures_public_scRNAseq.Rmd) and the corresponding PDF/HTML files.

## References

* Amezquita, Robert, Lun, Aaron, Hicks, Stephanie, and Gottardo, Raphael (version 1.0.6). "Orchestrating Single-Cell Analysis with Bioconductor". <https://bioconductor.org/books/release/OSCA/>. Nat Methods. 2020 Feb;17(2):137-145. doi: 10.1038/s41592-019-0654-x. PMID: 31792435.
* **scater**: McCarthy DJ, Campbell KR, Lun ATL, Willis QF (2017). "Scater: pre-processing, quality control, normalisation and
visualisation of single-cell RNA-seq data in R." _Bioinformatics_, *33*, 1179-1186. <https://doi.org/10.1093/bioinformatics/btw777>.
* **scran**: L. Lun, A. T., Bach, K., & Marioni, J. C. (2016). Pooling across cells to normalize single-cell RNA sequencing data with many zero counts. Genome Biology, 17(1), 75. <https://doi.org/10.1186/s13059-016-0947-7>
* **MNN correction**:
    - Haghverdi L, Lun ATL, Morgan MD, Marioni JC (2018). "Batch effects in single-cell RNA-sequencing data are corrected
by matching mutual nearest neighbors." _Nat. Biotechnol._, *36*(5), 421-427. <https://doi.org/10.1038/nbt.4091>
    - Lun, Aaron. Further MNN algorithm development. 2018. https://MarioniLab.github.io/FurtherMNN2018/theory/description.html
* **clustering etc**: Lun ATL, McCarthy DJ, Marioni JC (2016). "A step-by-step workflow for low-level analysis of single-cell RNA-seq data
with Bioconductor." _F1000Res._, *5*, 2122. <https://doi.org/10.12688/f1000research.9501>
* **destiny**: Angerer, P., Haghverdi, L., Büttner, M., Theis, F. J., Marr, C., & Buettner, F. (2016). Destiny: Diffusion maps for large-scale single-cell data in R. Bioinformatics, 3(8), 1241–1243. <https://doi.org/10.1093/bioinformatics/btv715>.
* **TSCAN**: Zhicheng Ji, Hongkai Ji, TSCAN: Pseudo-time reconstruction and evaluation in single-cell RNA-seq analysis, _Nucleic Acids Research_, Volume 44, Issue 13, 27 July 2016, Page e117, <https://doi.org/10.1093/nar/gkw430>
* Galili, Tal (2015). dendextend: an R package for visualizing, adjusting, and comparing trees of hierarchical clustering. Bioinformatics. DOI: 10.1093/bioinformatics/btv428
* Marini, Federico and Binder, Harald (2019). pcaExplorer: an R/Bioconductor package for interacting with RNA-seq principal components. doi: 10.1186/s12859-019-2879-1
* Robinson MD, Oshlack A (2010). A scaling normalization method for differential expression analysis of RNA-seq data. Genome Biology 11, R25.
* Schauder DM, Shen J, Chen Y, Kasmani MY et al. E2A-regulated epigenetic landscape promotes memory CD8 T cell differentiation. Proc Natl Acad Sci U S A 2021 Apr 20;118(16). PMID: 33859041
* Yao C, Sun HW, Lacey NE, Ji Y et al. Single-cell RNA-seq reveals TOX as a key regulator of CD8<sup>+</sup> T cell persistence in chronic infection. Nat Immunol 2019 Jul;20(7):890-901. PMID: 31209400

## Package versions

- STAR v2.6.0c
- CellRanger v.5.0.0 with mm10-2020-A

>R version 4.0.3 (2020-10-10)

- AnnotationDbi_1.50.3 
- batchelor_1.6.0
- CiteFuse_1.2.0
- clusterProfiler_3.16.1 (v3.14.3 for bulk RNA-seq)
- DEGreport v1.22.0
- DESeq2 v1.26.092
- destiny_3.4.0
- DropletUtils_1.10.0
- edgeR_3.30.3
- featureCounts v1.6.2 
- fgsea v1.12.0 
- GenomeInfoDb_1.26.0 
- GenomeInfoDbData_1.2.4
- ggalluvial_0.12.3
- ggplot2_3.3.3
- GOstats_2.54.0
- igraph_1.2.6
- limma_3.46.0
- Matrix_1.3-2
- pcaExplorer v2.14.2
- TSCAN_1.28.0
- scater_1.18.0
- scuttle_1.0.0
- SingleCellExperiment_1.12.0
- SummarizedExperiment_1.20.0
- scran_1.18.0

## Data for download

The raw data (fastq files, read counts from CellRanger) can be downloaded from GEO (GSE151652).
Single-cell data are stored under accession no. GSM5641677 through GSM5641682.
The remaining samples are bulk RNA-seq.

For the single-cell data, there are multiple data set types per sample:

* **draining pancreatic lymph node** (pLN) - three technical replicates (pLN_1, pLN_3, pLN_4)
    * single-cell expression values (barcodes, features, matrix)
    * hash-tagged oligos to mark which individual mouse a given barcode/cell belongs to (see [mouse_per_cell.csv](https://github.com/abcwcm/GeartySchietinger/blob/master/scRNAseq/data/mouse_per_cell.csv) for the processed file; [samples_Gearty_singleCell_hashtagBarcodes.csv](https://github.com/abcwcm/GeartySchietinger/blob/master/scRNAseq/data/samples_Gearty_singleCell_hashtagBarcodes.csv) contains the actual HTO barcodes)
    * tagged antibodies (ADT) against CD44, CD62L, CD127, CD73, PD-1, CD38, CD39 (see [ADTBarcodes.csv](https://github.com/abcwcm/GeartySchietinger/blob/master/scRNAseq/data/ADTBarcodes.csv))
    * V(D)J sequencing 
* **pancreas** - three technical replicates (panc_1, panc_3, panc_4)
    * hash-tagged oligos 
    * V(D)J sequencing

>The easiest way to get started is to use the processed data provided here.

For the single-cell data, some of the data can be downloaded from Box in the form of RDS (load into R via `in_data <- readRDS()`) or RDA objects (load into R via `load()`).
The `data/` directory in the scRNA-seq directory contains some text files that contain just the cell labels and the mouse labels for individual cells.

| File_name |	Robject_type	| Details |
|--------------|-------------|------------|
| [sce_integrated_pLNSamples_filtered.rds](https://wcm.box.com/shared/static/2wvbdfs3vja2cnlckcqpk20eu1693o8o.rds) | RDS |	`SingleCellExperiment` object with logcounts, reduced dimensionality coordinates, labels and so on for the integrated pLN data set that is the basis for most of the scRNA-seq figures in Gearty et al.; code details can be found in `scRNA_seq/processing_scRNA-seq.md` and `scRNA_seq/figures_scRNAseq.Rmd` |
| [GOenrich_BP_pLNmarkers.rda](https://wcm.box.com/shared/static/8k2cudvwqxw05iyywz8d36oof2qavjb8.rda) |	RDA	| Result of `clusterProfiler`'s `goEnrich()` on marker genes of pLN populations with GO terms biological processes |
| [GOenrich_MF_pLNmarkers.rda](https://wcm.box.com/shared/static/3zstp6yohqvrc31t26tzs8ox5pr5ukm6.rda) |	RDA	| Result of `clusterProfiler`'s `goEnrich()` on marker genes of pLN populations with GO terms molecular functions |
| [tcr_perCell_pairedSamples_newCID_2021-05.rds](https://wcm.box.com/shared/static/5wj8h9qssai2u0ya2sre2tsudahvvoni.rds) | RDS	| clonotype IDs and labels for cells that were found in pLN as well as pancreas samples from the same mice ("paired") |
| [tscan_lineData_DiffMap.rda](https://wcm.box.com/shared/static/rakg00duv7j5v18e9jw2jq93px187l3h.rda) |	RDA | `TSCAN` results |
| [sce_integrated_with_YaoProgsMpecs-and-Schauder.rds](https://wcm.box.com/shared/static/ixtkzzow5b4r2ck3u6zit2jtep7dunzv.rds) | RDS	| SingleCellExperiment object with logcounts and labels following the integration of Gearty's pLN data as well as MPECs and ProgLikes from [Yao et al. (2019)](https://dx.doi.org/10.1038/s41590-019-0403-4) and d129-Tcm from [Schauder et al (2021)](https://dx.doi.org/10.1073/pnas.2013452118)!

