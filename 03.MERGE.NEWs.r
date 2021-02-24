rm(list=ls())
library(Seurat)
library(Seurat)
library(ggplot2)
library(readr)
library(readxl)
library(Matrix)
library(tidyr)
library(dplyr)
library(readxl)
library(cowplot) 
library(knitr)
library(harmony)
library(markdown) 
setwd("/rsrch3/scratch/lym_myl_rsch/qcai1/SC145")
dataDir  =  "/rsrch3/scratch/lym_myl_rsch/qcai1/SC145/Data.RDS/"
samplePDX <- read_xlsx("/rsrch3/scratch/lym_myl_rsch/qcai1/SC145/Data/PDX_scRNA_WES_info.20201023cqs.xlsx", sheet = 1)
unique(as.character(samplePDX$scSeqName))
data.frame(samplePDX)

allSampleName = list.files(dataDir)
#allSampleName <- allSampleName[!allSampleName %in% "scripts"]
allSampleName  <- gsub(".RDS","",allSampleName)

initialSampleName = "S10"
allSampleName = allSampleName[!allSampleName %in% initialSampleName]


i <- initialSampleName
Object1 <- readRDS(paste0("./Data.RDS/", i , ".RDS"))
Object1@meta.data$Pt <-    samplePDX[samplePDX$scSeqName %in% i,]$Pt   ######################################
Object1@meta.data$sample <-    samplePDX[samplePDX$scSeqName %in% i,]$scSeqName
Object1@meta.data$Newname <-    samplePDX[samplePDX$scSeqName %in% i,]$Newname
Object1@meta.data$HumanMouse <-    samplePDX[samplePDX$scSeqName %in% i,]$HumanMouse
Object1@meta.data$PDXmodel <-    samplePDX[samplePDX$scSeqName %in% i,]$PDXmodel
Object1 <- NormalizeData(Object1)
Object1 <- FindVariableFeatures(Object1, selection.method = "vst", nfeatures= 2000)
Object1 <- ScaleData(Object1)
Object1 <- subset(Object1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)               #REQUESTED BY LINGHUA 20


Object.combined = c()
for ( i in allSampleName ) {
    print(i)
    j = i
    Object32 <- readRDS(paste0("./Data.RDS/", i , ".RDS"))
    
Object32@meta.data$Pt <-    samplePDX[samplePDX$scSeqName %in% i,]$Pt   ######################################
Object32@meta.data$sample <-    samplePDX[samplePDX$scSeqName %in% i,]$scSeqName
Object32@meta.data$Newname <-    samplePDX[samplePDX$scSeqName %in% i,]$Newname 
Object32@meta.data$HumanMouse <-    samplePDX[samplePDX$scSeqName %in% i,]$HumanMouse
Object32@meta.data$PDXmodel <-    samplePDX[samplePDX$scSeqName %in% i,]$PDXmodel
Object32[["percent.mt"]] <- PercentageFeatureSet(Object32, pattern = "^MT-")

Object32 <- NormalizeData(Object32)
Object32 <- FindVariableFeatures(Object32, selection.method = "vst", nfeatures= 2000)
Object32 <- ScaleData(Object32)
Object32 <- subset(Object32, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)   #REQUESTED BY LINGHUA 20
Object32    
   Object.combined = c(Object.combined, Object32)
  }
 

Objectall <- merge(Object1, y =   Object.combined, add.cell.ids =c(initialSampleName, allSampleName)  , project = "SC145")
 
 
 
ObjectPX@meta.data$Symbol <- "NNN"
ObjectPX@meta.data$Generation <- "NNN"
ObjectPX@meta.data$Finalname <- "NNN"
ObjectPX@meta.data$Pt <- "NNN"
for (i in unique(Objectall@meta.data$orig.ident)) {
Objectall@meta.data[Objectall@meta.data[,"orig.ident"] %in%  i,]$Symbol     <- samplePDX[samplePDX$scSeqName %in% i,]$Symbol
Objectall@meta.data[Objectall@meta.data[,"orig.ident"] %in%  i,]$Pt         <- samplePDX[samplePDX$scSeqName %in% i,]$Pt
Objectall@meta.data[Objectall@meta.data[,"orig.ident"] %in%  i,]$Finalname  <- samplePDX[samplePDX$scSeqName %in% i,]$Finalname
Objectall@meta.data[Objectall@meta.data[,"orig.ident"] %in%  i,]$Generation <- samplePDX[samplePDX$scSeqName %in% i,]$Generation
} 


 saveRDS(Objectall, "./Data/Object_ABCDEF.H&M.0223.2021.rds")
 
 
 
 pt1 <-"ABCDEF"
 part <- ".all_cell"  # ".Bcell_only"  "all"
 #ObjectPX <- readRDS("./Data/Object_ABCDEF.H&M.0223.2021.rds") # include all
 #ObjectPX <- subset(ObjectPX,       subset = Pt %in% c(pt1, "P") )  
ObjectPX <- Objectall


 ObjectPX <- NormalizeData(object = ObjectPX, normalization.method = "LogNormalize", scale.factor = 10000)  #, verbose = FALSE)
 ObjectPX <- FindVariableFeatures(object = ObjectPX, selection.method = "vst", nfeatures = 2000)
 ObjectPX <- ScaleData(ObjectPX,verbose = FALSE)
 ObjectPX <- RunPCA(ObjectPX, features = VariableFeatures(object = ObjectPX))
 ObjectPX <- FindNeighbors(ObjectPX, dims = 1:10)
 ObjectPX <- FindClusters(ObjectPX, resolution = 0.5)
 #ObjectPX <- RunUMAP(ObjectPX, dims = 1:30)
 ObjectPX <- RunTSNE(ObjectPX, dims = 1:30)
 


 pt1 <-"ABCDEF"
 roww = 4 
 part <- ".all_cell"                             # ".Bcell_only" 
 knit("/rsrch3/scratch/lym_myl_rsch/qcai1/SC145/report/TSNE.RMD", paste0(pt1, "_PBMC", part,".RMD"))
 markdownToHTML(paste0(pt1, "_PBMC", part,".RMD"), paste0("./report/",pt1, "&PBMC", part,".html"))

 
