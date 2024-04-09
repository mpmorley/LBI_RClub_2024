library(Seurat)
library(tidyverse)
library(UCell)


#Install Tthe SCP library 
#https://zhanghao-njmu.github.io/SCP/
library(SCP)
#Install teh UCell library
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("UCell")
library(UCell)



# Functions ---------------------------------------------------------------


named_group_split <- function(.tbl, ...) {
  grouped <- group_by(.tbl, ...)
  names <- rlang::inject(paste(!!!group_keys(grouped), sep = " / "))
  
  grouped %>% 
    group_split() %>% 
    rlang::set_names(names)
}



# Create a Seurat Object --------------------------------------------------


counts <- Read10X_h5('Filtered.h5')
scrna <- CreateSeuratObject(counts = counts)



# Access functions -------------------------------------------------------

scrna

#Get Cells and features(genes)
colnames(scrna)
rownames(scrna)

#Functions do the same as above
Cells(scrna)
Features(scrna)

#Using grep to find genes

grep('Rps|Rpl',rownames(scrna))

grep('Rps|Rpl',rownames(scrna),value = T)



# Meta Data ---------------------------------------------------------------


#Get meta table 
scrna[[]]


#Direct referencing to "slots" should be avoided since the authors could change name and or remove. I do it a lot anyway
scrna@meta.data


#Add some new fields

scrna$samplename <- 'CD1_control_P5_male'
scrna$treatment <- 'control'
scrna$strain <- 'CD1'
scrna$age <- 'P5'
scrna$sex <- 'M'

scrna[[]]


#Compute other QC values

scrna[["percent.mito"]] <- PercentageFeatureSet(scrna, pattern = "^Mt-")
scrna[["percent.ribo"]] <- PercentageFeatureSet(scrna, pattern = "^Rps|^Rpl")

VlnPlot(scrna, features = c("nFeature_RNA", "nCount_RNA", "percent.mito","percent.ribo"), ncol = 4)

scrna[[]]


# Standard Analysis  ------------------------------------------------------


scrna <- NormalizeData(object = scrna)
scrna <- FindVariableFeatures(object = scrna)
scrna <- ScaleData(object = scrna)
scrna <- RunPCA(object = scrna)
scrna <- FindNeighbors(object = scrna, dims = 1:30)
scrna <- FindClusters(object = scrna)
scrna <- RunUMAP(object = scrna, dims = 1:30)


# SCtransform Version -----------------------------------------------------


scrna <- SCTransform(object = scrna)
scrna <- RunPCA(object = scrna)
scrna <- FindNeighbors(object = scrna, dims = 1:30)
scrna <- FindClusters(object = scrna)
scrna <- RunUMAP(object = scrna, dims = 1:30)


scrna

#We can use pipes!!! 
scrna <- SCTransform(scrna,vars.to.regress = c('nFeature_RNA')) %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters() %>%
  RunUMAP(dims = 1:30)


scrna

scrna[[]]



# Basic QC stats ----------------------------------------------------------

#Using Seurat to plot

#Seurats Identy class

#Revist the QC vlnplot, it's different! Seurat updated the Identy after we clustered 
VlnPlot(scrna, features = c("nFeature_RNA", "nCount_RNA","percent.ribo"), ncol = 4)

Idents(scrna)

#We can set the Idents to anything in meta.data

#Set sample_name. 

Idents(scrna) <- 'samplename'

VlnPlot(scrna, features = c("nFeature_RNA", "nCount_RNA","percent.ribo"), ncol = 4)

# Seurat plotting functions have group.by argument to define 

VlnPlot(scrna, features = c("nFeature_RNA", "nCount_RNA","percent.ribo"), group.by = 'seurat_clusters',ncol = 4)

#Change Ident back to cluster info
Idents(scrna) <- 'seurat_clusters'




# Seurat Visualtions ------------------------------------------------------

# Plot is for discrete categorical data
#Using the ident
DimPlot(scrna)

#Same but defineing it,which is safer. 

DimPlot(scrna,group.by = 'seurat_clusters')

#We can do some cleaning of this plot

DimPlot(scrna,group.by = 'seurat_clusters') + coord_equal() + theme_void()

#Featureplot is for continuous data, like gene expression

FeaturePlot(object=scrna,features='Ager')

FeaturePlot(object=scrna,features=c('Ager','Lamp3'))

FeaturePlot(object=scrna,features='Ager',order = T) + coord_equal() + theme_void()


#Let's build a quick plot using patchwork. 


d1 <- DimPlot(scrna,group.by = 'seurat_clusters') + coord_equal() + theme_void()
f1 <- FeaturePlot(object=scrna,features='Ager',order = T) + coord_equal() + theme_void()


d1|f1




# Nicer Plots using the SCP pacakge ---------------------------------------

#Let's revist the vlnplot for QC
#Let's make a variable called wcfeatires so we don;'t need to type this in for everplot.
qcfeatures <- c("nCount_RNA","nFeature_RNA","percent.ribo")

FeatureStatPlot(scrna, stat.by =qcfeatures , group.by = "seurat_clusters",ncol = 1)

#Cleanup with stacking and no legend
v1 <- FeatureStatPlot(scrna, stat.by =qcfeatures , group.by = "seurat_clusters",stack = T,legend.position ='none')
v1



#Basic default plot, SCP does not use IdentsQ  

CellDimPlot(scrna, group.by = "seurat_clusters")

#We can do a bit better, let's label clusters and make pretty axes. 


d1 <- CellDimPlot(scrna, group.by = "seurat_clusters", reduction = "UMAP",theme_use = 'theme_blank', label = T,label_insitu = T)
d1


f1 <- FeatureDimPlot(scrna,qcfeatures ,theme_use = 'theme_blank',ncol=1)

f1

(d1|f1)/v1



# Run Findmarkers ---------------------------------------------------------
fm <- FindAllMarkers(scrna,only.pos = T) %>% mutate(pctdiff = pct.1-pct.2)


#Use dplyr to take the top 7 genes with the largest % diff. 
diffgenes <-
  fm %>% group_by(cluster) %>% top_n(n = 7, wt = pctdiff)

#Use SCP to make a heatmap. 
ht1 <- GroupHeatmap(
  scrna,
  features = diffgenes$gene,
  group.by = 'seurat_clusters',
  show_row_names = T,
  row_names_side = 'left',
  cluster_rows = F
)

ht1

#Create a list of data frames by spltting the large data frame by cluster

dge.list <- fm %>% named_group_split(cluster)
dge.list

#The write_xlsx can take a list of data frames and then save them to seperate sheets in an excel workbook. 
writexl::write_xlsx(dge.list,path = paste0(plotdir,'FindMarkers.xlsx'))


# AT2 Score ---------------------------------------------------------------


gene.sets <- list(AT2=c('Sftpc', 'Lamp3', 'Slc34a2', 'Sftpb', 'Sftpa1', 'Cxcl15', 'Npc2', 'Dram1', 'Sfta2', 'Bex4', 'Ppp1r14c', 'Abca3', 'Sftpd'))
scrna <- AddModuleScore_UCell(scrna, features = gene.sets,name='_score')
f3 <- FeatureDimPlot(scrna,'AT2_score',theme_use = 'theme_blank')
v3 <- FeatureStatPlot(scrna, stat.by = gene.sets$AT2, group.by = "seurat_clusters",stack = T,xlab = 'Cluster',legend.position = 'none')

(d1|f3)/v2








