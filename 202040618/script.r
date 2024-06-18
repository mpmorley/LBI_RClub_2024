library(Seurat)
library(tidyverse)
library(patchwork)
library(SCP)


## Download the Seurat object from Dropbox, see my email

scrna.all <- readRDS("SeuratTest.RDS")


crna.all@meta.data <- scrna.all[[]] |> select(-var_sample, -seurat_clusters,-var_cluster,-starts_with('SCT_snn'))

head(scrna.all[[]])

#We don't have a meta.data field for which sample, all we have is our library ID. We can use mapvalues function to gor from lib ID to sample name
table(scrna.all$orig.ident)

scrna.all$sample <- plyr::mapvalues(scrna.all$orig.ident,from = c("EEM-scRNA-305","EEM-scRNA-309","EEM-scRNA-313"),to=c('Cdc42KO_d14','Cdc42KO_d28','Ctrl'))
                                    
table(scrna.all$sample)

#Let's make some UMAP DimPlots. 

#Seurat version
DimPlot(scrna.all,group.by = 'sample')

#Nicer Seurat version
DimPlot(scrna.all,group.by = 'sample') + theme_void() + coord_equal()


#SCP version
CellDimPlot(scrna.all, group.by = "sample", reduction = "UMAP",theme_use = 'theme_blank', label = T,label_insitu = T,label_repulsion = T)

#Let's look at celltypes as well, We see some celtypes we are not intested in

CellDimPlot(scrna.all, group.by = "predicted.celltype", reduction = "UMAP",theme_use = 'theme_blank', label = T,label_insitu = T,label_repulsion = T)

#We see some celtypes we are not interested in, so let's remove them

scrna <- subset(scrna.all,subset = predicted.celltype %in% c('Cap1','Cap2','Cap_postflu','Cap1_2'))

#We need to reprocess some steps with those missing cells now. 

head(scrna[[]])


scrna <- FindVariableFeatures(scrna) |>
  RunPCA() %>%
  FindNeighbors(dims = 1:20) |>
  FindClusters() |>
  RunUMAP(dims = 1:20)



head(scrna[[]])


#Let's plot again.. 

CellDimPlot(scrna, group.by = "sample", reduction = "UMAP",theme_use = 'theme_blank', label = T,label_insitu = T,label_repulsion = T)
#We still have big difference between Ctrl and KO samples

#What genes are drving this, we can compare the samples using findallmarkers, we first to change sample to ident to do so

Idents(scrna) <- 'sample'

fm <- FindAllMarkers(scrna,only.pos = T) |> mutate(diffpct=pct.1-pct.2)


plotgenes <- fm |> slice_max(order_by=avg_log2FC,n = 10) |> pull(gene)


FeatureStatPlot(scrna,stat.by = plotgenes,
                legend.position = "top", legend.direction = "horizontal",
                group.by = "sample", bg.by = "sample", stack = TRUE
)

#What if sort by %diff?


plotgenes <- fm |> slice_max(order_by=diffpct,n = 10) |> pull(gene)
  FeatureStatPlot(scrna,stat.by = plotgenes,
                  legend.position = "top", legend.direction = "horizontal",
                  group.by = "sample", bg.by = "sample", stack = TRUE
  )

#Seems like some ribosome genes too.. Let's plot those.   
  plotgenes <- fm |> filter(grepl('Rps|Rpl',gene)) |> slice_max(order_by=diffpct,n = 10) |> pull(gene)
  
  FeatureStatPlot(scrna,stat.by = plotgenes,
                  legend.position = "top", legend.direction = "horizontal",
                  group.by = "sample", bg.by = "sample", stack = TRUE
  )
  
  
  
#We don't use genes directly but rather PCs for UMAP and cluster, we can see which genes most influence each PC
  
VizDimLoadings(scrna)
  

#Let's instead use map to plot the first 10 PCs

map(1:10,function(pcnum){
  print(str_glue("Plotting PC{pcnum}\n"))
    VizDimLoadings(scrna,dims = pcnum)
})


#Let's lok at some of these genes for PC1

FeatureDimPlot(scrna,features = c('Mgp','Cldn5'),legend.position = 'none',theme_use = 'theme_blank')



#What if now removed the highly expressed X genes, ribo and mito genes.  

scrna <- scrna[c(grep('Rps|Rpl|Xist|Tsix|mt',rownames(scrna),value = T,invert = T)),]

scrna <- FindVariableFeatures(scrna) |>
  RunPCA() %>%
  FindNeighbors(dims = 1:8) |>
  FindClusters() |>
  RunUMAP(dims = 1:8)
scrna$cluster <- Idents(scrna)



d1 <- CellDimPlot(scrna, group.by = "orig.ident", reduction = "UMAP",theme_use = 'theme_blank', label = T,label_insitu = T,label_repulsion = T)
d2 <- CellDimPlot(scrna, group.by = "predicted.celltype", reduction = "UMAP",theme_use = 'theme_blank', label = T,label_insitu = T,label_repulsion = T,palette = 'Categorical.12')
d3 <- CellDimPlot(scrna, group.by = "seurat_clusters", reduction = "UMAP",theme_use = 'theme_blank', label = T,label_insitu = T,label_repulsion = T,palette = 'Categorical.12')
d1 + d2

b1 <- CellStatPlot(scrna,stat.by = 'orig.ident',group.by='cluster')

d3+b1

head(scrna[[]])


# Eval cluster resolutions -------------------------------------------------

#Findclusters can take a vector of resiltions to perform
scrna <- scrna |>
FindNeighbors(dims = 1:8) |>
  FindClusters(resolution = seq(from=0.2,to=1.2,by=.2)) 

#View all the new resolutions
head(scrna[[]])


#Pick one to plot. 
CellDimPlot(scrna, group.by = "SCT_snn_res.0.2", reduction = "UMAP",theme_use = 'theme_blank', label = T,label_insitu = T,label_repulsion = T,palette = 'Categorical.12')


#Let's construct a plot of all the different res. We will first construct a vector of all the resolutions names. 

clsres <- paste0('SCT_snn_res.',seq(from=0.2,to=1.2,by=.2))

map(clsres,\(name){
  CellDimPlot(scrna, group.by = name, reduction = "UMAP",theme_use = 'theme_blank', label = T,label_insitu = T,label_repulsion = T,palette = 'Categorical.12',title = name)
})

#We will now "map" over all teh resultion names and make a plot, the map function will return a list of ggplot objects. 
plotlist <- map(clsres,\(name){
  CellDimPlot(scrna, group.by = name, reduction = "UMAP",theme_use = 'theme_blank', label = T,label_insitu = T,label_repulsion = T,legend.position = 'none', palette = 'Categorical.12',title = name)
})

#We can now use the wrap_plots function in patchwork to make a single plot
wrap_plots(plotlist)


