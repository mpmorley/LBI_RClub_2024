library('Seurat')
library(tidyverse)
source('bin/ColorPal.R')

outdir <- 'Results/Step2_subset/'
dir.create(outdir,showWarnings = F)

all <- readRDS('Results/Step1_merge/Seurat.RDS')



# Subset for Endos --------------------------------------------------------


#Subset both on Cluster and predicted celltype. 
scrna <- subset(all, subset=seurat_clusters %in% c(0,3,4,23) & predicted.celltype %in% c('Cap_postflu', 'Cap1','Cap1_2','Cap2'))
scrna <- FindVariableFeatures(scrna)
scrna <- RunPCA(scrna)
scrna <- RunUMAP(scrna,dims = 1:10)
scrna <- FindNeighbors(scrna, dims = 1:10)
scrna <- FindClusters(scrna, resolution = 0.4)



d1 <- CellDimPlot(scrna,group.by = 'predicted.celltype',label = TRUE, label_insitu = TRUE, label_repel = TRUE,theme_use = "theme_blank", label.fg = "black",label.bg.r = 0,legend.position = 'none')
d2 <- CellDimPlot(scrna,group.by = 'seurat_clusters',label = TRUE, label_insitu = TRUE, label_repel = TRUE,theme_use = "theme_blank", label.fg = "black",label.bg.r = 0,legend.position = 'bottom')
d1|d2


# QC plots for review -----------------------------------------------------

qc.plots <- 
  list(
    FeatureStatPlot(scrna, stat.by = c("nFeature_RNA", "nCount_RNA", "percent.mito", 'predicted.celltype.score','hybrid_score', 'scrublet_score'), 
                    group.by = "seurat_clusters", add_box = TRUE, stack = TRUE, legend.position = 'none'),
    CellStatPlot(scrna, group.by = "seurat_clusters", stat.by = "predicted.celltype", plot_type = "bar"),
    CellStatPlot(scrna, group.by = "seurat_clusters", stat.by = c("sample"), plot_type = "bar"),
    CellStatPlot(scrna, group.by = "seurat_clusters", stat.by = "Phase", plot_type = "bar")
  )

qc.plots <- setNames(qc.plots, paste0('QC_subset_', 1:length(qc.plots)))
imap(qc.plots, ~ggsave(filename = str_glue('{outdir}/QCPlots/{.y}.png'),.x,height=8,width=10))



# Plot Cap cell markers ---------------------------------------------------


genelist <- c('Car4','Ednrb', 'Gpihbp1')

f1 <- FeatureDimPlot(scrna,genelist,theme_use = "theme_blank",ncol = 1)

v1 <- FeatureStatPlot(scrna, stat.by =genelist, 
                group.by = "seurat_clusters",add_box = TRUE, stack = TRUE, legend.position = 'none')

d2|f1


# Cell Annotation ---------------------------------------------------------
#Read in the cell anno, I added a column called "order" this is the order I want the cells to appear on plots etc. 
cellanno <- readxl::read_xlsx(str_glue('{outdir}/data/CellAnno.xlsx')) |> arrange(order)

scrna$celltype <- plyr::mapvalues(scrna$seurat_clusters,from=cellanno$cluster,to=cellanno$celltype)
scrna$celltype <- factor(scrna$celltype,levels = c('Cap1','Cap1_Cap2','Cap2'))
scrna$celltype.level2 <- plyr::mapvalues(scrna$seurat_clusters,from=cellanno$cluster,to=cellanno$celltype.level2)
scrna$celltype.level2 <- factor(scrna$celltype.level2 ,levels = cellanno$celltype.level2)

# #Explore some of the Cap1 and Cap1_Cap2 subclusters.  -------------------

#make a list of what to compare. 

comparison.list <- list(
  'Cap1_a' = list(a='Cap1_a',b=c('Cap1_b','Cap1_c')),
  'Cap1_b' = list(a='Cap1_b',b=c('Cap1_a','Cap1_c')),
  'Cap1_c' = list(a='Cap1_c',b=c('Cap1_a','Cap1_b')),
  'Cap1_Cap2_a' = list(a='Cap1_Cap2_a',b=c('Cap1_Cap2_b','Cap1_Cap2_c')),
  'Cap1_Cap2_b' = list(a='Cap1_Cap2_b',b=c('Cap1_Cap2_a','Cap1_Cap2_c')),
  'Cap1_Cap2_c' = list(a='Cap1_Cap2_c',b=c('Cap1_Cap2_a','Cap1_Cap2_b'))
)


#We have to set the seurat celltype.level2 to the Seurat Ident to be used with Findmarkers 
Idents(scrna) <- 'celltype.level2'

comparison.res <- imap(comparison.list,\(c,name){
  FindMarkers(scrna,ident.1 = c[['a']],ident.2 = c[['b']],only.pos=T) |>
    rownames_to_column('gene') |>
    mutate(comparison=name,
         pct.diff=pct.1-pct.2
         )
})



# Plot the top 5 genes per cell type --------------------------------------

comparison.res |> bind_rows() |> group_by(comparison) |>
  slice_max(pct.diff,n=5) |> select(gene) |> distinct() |> pull(gene) |>
  FeatureStatPlot(scrna, stat.by =_, 
                  group.by = "celltype.level2",add_box = TRUE, stack = TRUE, 
                  legend.position = 'none',
                  palcolor = pal.list$endo.lvl2,
                  bg.by='celltype',
                  bg_palcolor = pal.list$endo
  )


#Appears we might have some ambient RNAs  
genelist = c('Scgb1a1','Sftpc','Lamp3')


FeatureStatPlot(scrna, stat.by =genelist, 
                group.by = "celltype.level2",add_box = TRUE, stack = TRUE, legend.position = 'none',split.by='sample',
                palcolor = pal.list$sample,
                bg.by='celltype',
                bg_palcolor = pal.list$endo
                )


#We also seem to have Hb genes defining Cap2_Cap2_c as well. 
genelist = c('Hba-a1','Hbb-bs','Hbb-bt')
FeatureStatPlot(scrna, stat.by =genelist, 
                group.by = "celltype.level2",add_box = TRUE, stack = TRUE, legend.position = 'none',split.by='sample',
                palcolor = pal.list$sample,
                bg.by='celltype',
                bg_palcolor = pal.list$endo
)



# Refine Dims for UMAP and clustering -------------------------------------


#if we look back, 1A1 is affecting PC3 and the Hbs PC7,PC8. 


VizDimLoadings(scrna, dims = 1:4, reduction = "pca")
VizDimLoadings(scrna, dims = 5:8, reduction = "pca")
VizDimLoadings(scrna, dims = 9:12, reduction = "pca")
VizDimLoadings(scrna, dims = 13:16, reduction = "pca")



dims <-  c(1:2,4:6,9:12)
scrna <- RunUMAP(scrna,dims = dims)
scrna <- FindNeighbors(scrna, dims = dims)
scrna <- FindClusters(scrna, resolution = 0.4)






p1 <- CellDimPlot(scrna,group.by = 'seurat_clusters',label = TRUE, 
            label_insitu = TRUE, label_repel = TRUE,
            theme_use = "theme_blank",
            legend.position = 'none')

p2 <- FeatureStatPlot(scrna, stat.by =c('Gpihbp1','Car4'),
                group.by = "seurat_clusters",
                add_box = TRUE, 
                stack = TRUE, 
                legend.position = 'none')

p3 <- CellStatPlot(scrna, group.by = "seurat_clusters", 
                   stat.by = c("sample"), 
                   plot_type = "bar",
                   palcolor = pal.list$sample)


(p1|p3)/p2


genelist <- c('Gpihbp1','Car4','Ednrb','Hba-a1','Hbb-bs','Scgb1a1','Sftpc','Lamp3')
FeatureStatPlot(scrna, stat.by =genelist,
                group.by = "seurat_clusters",
                add_box = TRUE, 
                stack = TRUE, 
                legend.position = 'none',
                #split.by='sample',
                #palcolor = pal.list$sample
)

#Set celltype without a file. 
scrna$celltype <- plyr::mapvalues(scrna$seurat_clusters,from=0:5,to=c('Cap1','Cap1','Cap2','Cap1','Cap1_Cap2','Cap1_Cap2'))
scrna$celltype <- factor(scrna$celltype,levels = c('Cap1','Cap1_Cap2','Cap2'))

#Just remember that the celltype.level2 is set to the older clusters, for now let's set it NULL so we don't use it. 
scrna$celltype.level2 <- NULL
scrna$celltype.level2 

p4 <- CellDimPlot(scrna,group.by = 'celltype',label = TRUE, 
                  label_insitu = TRUE, label_repel = TRUE,
                  theme_use = "theme_blank",
                  legend.position = 'none',
                  palcolor = pal.list$endo)

p1|p4


saveRDS(scrna,str_glue({'{outdir}/Seurat.RDS'}))

