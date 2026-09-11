# PIPseq analysis
# April 7th, Tuesday 2026
# ---------------------
# Lecture, basic scRNA-seq analysis: https://satijalab.org/seurat/reference/read10x
# Lecture, Hashtag analysis: https://satijalab.org/seurat/articles/hashing_vignette


# libraries
library(Seurat)
library(patchwork)
library(ggplot2)
library(dittoSeq)
library(future)
library(scales)

##################################################################################################
# settings
# setwd("~/RStudio_files/PHF19/scRNAseq_Eileen")
dir <- dirname(rstudioapi::getSourceEditorContext()$path); print(dir) # [1] "/Users/lees130/Library/CloudStorage/OneDrive-NYULangoneHealth/N18_PipSeqMouse_KernHast/06_SeuratObject_Pipseq"
setwd(dir);
create_UMAP <- TRUE

PipseqFileDir <- "/Users/lees130/Library/CloudStorage/OneDrive-NYULangoneHealth/N18_PipSeqMouse_KernHast/05_ProcessedCellRangerOutput"
##################################################################################################

# load data
# if(create_UMAP){

      ########## ================ ########## ================ ########## ================ ########## ================ 
      ## Step1. Read RNA and HTO matrices.   Make a seurat object (RNA)
      ########## ================ ########## ================ ########## ================ ########## ================ 
      ######### 1) cDNA expression
      ## paste0(PipseqFileDir, "/01_Pipseeker_cDNA/") has the "barcodes.tsv.gz", "features.tsv.gz", and "matrix.mtx.gz" files of cDNA expression
      cDNA_CountMatrix <- Read10X(data.dir=paste0(PipseqFileDir, "/01_Pipseeker_cDNA/"),gene.column=2, unique.features=TRUE,strip.suffix=FALSE); class(cDNA_CountMatrix); # list. It takes about 2 min
      dim(cDNA_CountMatrix) # 32227gene 34119cell 
      saveRDS(cDNA_CountMatrix, "cDNA_CountMatrix.rds")

      # create seurat obj
      cDNAExpSeuratObj <- CreateSeuratObject(counts=cDNA_CountMatrix,  project="MousePipseq")
      saveRDS(cDNAExpSeuratObj, "cDNAExpSeuratObj.rds")
      
      ########## 2) HTO
      ## paste0(PipseqFileDir, "/02_Pipseeker_HTO/") has the "barcodes.tsv.gz", "features.tsv.gz", and "matrix.mtx.gz" files of library2. 
      HTO.data <- Read10X(data.dir=paste0(PipseqFileDir, "/02_Pipseeker_HTO/")); class(HTO.data) # list # This takes 5 min.
      names(HTO.data) # [1] "Gene Expression" "Cell Hashing"   
      saveRDS(HTO.data, "HTO.data.rds")
      
      ########## ================ ########## ================ ########## ================ ########## ================ 
      ## Step2. Add HTO as a second assay
      ########## ================ ########## ================ ########## ================ ########## ================   
      ## Extract the HTO matrix
      HTO.data <- HTO.data[["Cell Hashing"]]
      dim(HTO.data) # 4 34119
      head(colnames(HTO.data)) # Barcodes: [1] "AAAAAAAAAACGGTCA" "AAAAAAAAAAGCCGGA" "AAAAAAAAACACGCGA" "AAAAAAAAACATGACG" "AAAAAAAAACGGTCAA" "AAAAAAAAACGGTTGA"
      
      cDNAExpSeuratObj[["HTO"]] <- CreateAssayObject(counts = HTO.data)
      
      # Make sure barcodes match:
      all(colnames(cDNAExpSeuratObj) == colnames(HTO.data)) # TRUE
      
      ########## ================ ########## ================ ########## ================ ########## ================ 
      ## Step3. Normalize HTO data (CRITICAL STEP) # HTO uses CLR normalization, not log-normalization:
      ########## ================ ########## ================ ########## ================ ########## ================       
      cDNAExpSeuratObj <- NormalizeData(cDNAExpSeuratObj, assay = "HTO", normalization.method = "CLR")      
      
      ########## ================ ########## ================ ########## ================ ########## ================ 
      ## Step4. Demultiplex cells
      ########## ================ ########## ================ ########## ================ ########## ================       
      cDNAExpSeuratObj <- HTODemux(cDNAExpSeuratObj, assay = "HTO", positive.quantile = 0.99)      
      # Cutoff for TotalSeq-A0301 : 143 reads
      # Cutoff for TotalSeq-A0302 : 133 reads
      # Cutoff for TotalSeq-A0303 : 163 reads
      # Cutoff for TotalSeq-A0304 : 2025 reads
      
      #### Inspect the results
      table(cDNAExpSeuratObj$HTO_classification.global)
      # Doublet Negative  Singlet    ## singlet: good cells. ## doublet: multiple tags ## Negative: no clear tag
      # 2106    16261    15752 
       
      table(cDNAExpSeuratObj$HTO_classification) # counts per tag
      # Negative                TotalSeq-A0301 TotalSeq-A0301_TotalSeq-A0302 TotalSeq-A0301_TotalSeq-A0303 TotalSeq-A0301_TotalSeq-A0304                TotalSeq-A0302 
      # 16261                          4420                           540                           575                           108                          4437 
      # TotalSeq-A0302_TotalSeq-A0303 TotalSeq-A0302_TotalSeq-A0304                TotalSeq-A0303 TotalSeq-A0303_TotalSeq-A0304                TotalSeq-A0304 
      # 613                           116                          5736                           154                          1159 
      
      ########## ================ ########## ================ ########## ================ ########## ================ 
      ## Step6. Visualize HTO separation (VERY important) or Feature scatter plot 
      ########## ================ ########## ================ ########## ================ ########## ================       
      HTOseparationPlot <-  RidgePlot(cDNAExpSeuratObj, assay = "HTO", features = rownames(cDNAExpSeuratObj[["HTO"]]))
      # Picking joint bandwidth of 0.0826
      # Picking joint bandwidth of 0.0793
      # Picking joint bandwidth of 0.0739
      # Picking joint bandwidth of 0.0744
      OutHTOSepFileName <- "Out_HTOseparationPlot.pdf"  # Line 80 determine decreasing or increasing.
      ggsave(filename=OutHTOSepFileName, HTOseparationPlot, width=12, height=10)
      
     
      FeatureScatterPlot<- FeatureScatter(cDNAExpSeuratObj, feature1 = "TotalSeq-A0301", feature2 = "TotalSeq-A0302")
      OutFeatScatterFileName <- "Out_FeatureScatterPlot.pdf"  # Line 80 determine decreasing or increasing.
      ggsave(filename=OutFeatScatterFileName, FeatureScatterPlot, width=6, height=5)
# } else {
#       # load previous UMAP
#       phf19.seurat <- readRDS(paste0(getwd(),"/results/phf19.seurat.RDS"))
# }

      ########## ================ ########## ================ ########## ================ ########## ================ 
      ## Step7. Proceed with RNA analysis:  Normalize RNA, Find variable genes, SCTransform & PCA, Clust
      ###       From now, you can refer to # Lecture, basic scRNA-seq analysis: https://satijalab.org/seurat/reference/read10x
      ########## ================ ########## ================ ########## ================ ########## ================   

