library(Seurat)
library(tidyverse)
library(patchwork)

#We'll use this later 

BiocManager::install("MAST")

#Read in the data from 
scrna <- readRDS("SeuratTest.RDS")


#In the meta table the column sample contains the grouping variable of our samples, we have 3 groups. 
table(scrna$sample)

#Let's make the default ident to sample 
Idents(scrna) <- 'sample'


#We can test all the  Cdc42KO_d28 cells vs the the ctrl and Cdc42KO_d14 cells. In this case the other cells are the ref, so +logFC is an increase in Cdc42KO_d28
dge.1 <- FindMarkers(scrna, ident.1 = "Cdc42KO_d28") 
  

#Notice that genenames are rownames, let's move that to a column and create a new column with Pct diff          
dge.2 <- FindMarkers(scrna, ident.1 = "Cdc42KO_d28") |> 
  rownames_to_column('gene') |>
  mutate(pctdiff= pct.1 - pct.2)


#We can only compare 2 groups, but we can specify the 2nd group. Here we'll test d28 vs ctrl. 

dge.3 <- FindMarkers(scrna, ident.1 = "Cdc42KO_d28", ident.2 = 'Ctrl') |> 
  rownames_to_column('gene') |>
  mutate(pctdiff= pct.1 - pct.2)


#The default method for Diff Express is non-parametic wilcox test. We can try MAST.  
dge.4 <- FindMarkers(scrna, ident.1 = "Cdc42KO_d28", ident.2 = 'Ctrl',test.use = 'MAST') |> 
  rownames_to_column('gene') |>
  mutate(pctdiff= pct.1 - pct.2)

#With MAST we can use latent variables in the model, such as % mito, nFeatures Cell Cycle score, to adjust any effect they may have. 
dge.5 <- FindMarkers(scrna, ident.1 = "Cdc42KO_d28", ident.2 = 'Ctrl',test.use = 'MAST',latent.vars = c("percent.mito", "Phase")) |> 
  rownames_to_column('gene') |>
  mutate(pctdiff= pct.1 - pct.2)


#The entire we've been testing all the celltypes togather, but what if we want tested a single cell type. 

#IMPORTNAT We need to make celltype the ident 
Idents(scrna) <- 'var_celltype'
table(Idents(scrna))        

dge.cap2 <- FindMarkers(scrna, ident.1 = "Cdc42KO_d28", ident.2 = 'Ctrl', group.by = 'sample', subset.ident = 'Cap2') |> 
  rownames_to_column('gene') |>
  mutate(pctdiff= pct.1 - pct.2)

#A couple of ways to plot this using Seurat's VlnPlot. 
VlnPlot(scrna,c('Sftpc','Xist'),split.by = 'sample')
VlnPlot(scrna,c('Sftpc','Xist'),idents = 'Cap2',group.by = 'sample')
VlnPlot(scrna,c('Sftpc','Xist'),idents = 'Cap2',split.by = 'sample')



#What if we want to plot all the celltypes? We can use map to loop over. 

#We can use base R. 
celltypes.list <- unique(scrna$var_celltype)

#Or we can filter out celltypes that don't have many cells, This is just on all the cells regardless of sample, try figuring out to do 
#improve this filter. 


celltypes.list <- scrna[[]] |> group_by(var_celltype,sample) |> summarise(n=n()) |> filter(n>100) |> distinct(var_celltype) |> pull(var_celltype)

# DiffExp in single cell space --------------------------------------------


diffexp.d28 <- map(celltypes.list,\(celltype){
  FindMarkers(scrna, ident.1 = "Cdc42KO_d28", ident.2 = 'Ctrl', group.by = 'sample', subset.ident = celltype) |> 
    rownames_to_column('gene') |>
    mutate(pctdiff= pct.1 - pct.2)
}) |> set_names(celltypes.list)


#Let's merge all the results into a single sheet. 

all.diffexp <- bind_rows(diffexp.d28)


#However there's no celltype!! 

#We can loop the list of results and add it, we need to use imap. The list element names are celltype, so we can get the celltypename from there.

diffexp.d28 <- imap(diffexp.d28, \(res,celltype){
  res |> mutate(celltype=celltype)
})

all.diffexp <- bind_rows(diffexp.d28)

#We could also modify the Diff Exp code and add it there.  
diffexp.d28 <- map(celltypes.list,\(celltype){
  FindMarkers(scrna, ident.1 = "Cdc42KO_d28", ident.2 = 'Ctrl', group.by = 'sample', subset.ident = celltype) |> 
    rownames_to_column('gene') |>
    mutate(pctdiff= pct.1 - pct.2,
           celltype=celltype
           )
}) |> set_names(celltypes.list)

all.diffexp <- bind_rows(diffexp.d28)



#We can do some quick summaries. 


all.diffexp |> 
  filter(p_val_adj<0.05) |>
  group_by(celltype) |> 
  summarise(sig_genes = n())

#We can use if_else to classify up/down genes 
all.diffexp |> 
  mutate(direction = if_else(avg_log2FC > 0,'Up','Down')) |>
  filter(p_val_adj<0.05) |>
  group_by(celltype,direction) |> 
  summarise(sig_genes = n())


#We use tidyr to then make a little table. 
all.diffexp |> 
  mutate(direction = if_else(avg_log2FC > 0,'Up','Down')) |>
  filter(p_val_adj<0.05) |>
  group_by(celltype,direction) |> 
  summarise(sig_genes = n()) |>
  pivot_wider(names_from = direction,values_from = sig_genes)



# How find genes alterred in all celltypes, we have 5 celltypes

genelist <- all.diffexp  |> filter(p_val_adj < 0.05 & abs(pctdiff) > .25 ) |>
  group_by(gene) |> summarise(count = n()) |> filter(count==5) |> pull(gene)


all.diffexp  |> filter(p_val_adj<0.05 & abs(pctdiff) > .25 & gene %in% genelist) |>
  view()


# How find genes different in all celltypes

uniquegenes <- all.diffexp  |> filter(p_val_adj < 0.05 & abs(pctdiff) > .25 ) |>
  group_by(gene) |> summarise(count = n()) |> filter(count==1) |> pull(gene)


all.diffexp  |> filter(p_val_adj<0.05 & abs(pctdiff) > .25 & gene %in% uniquegenes) |>
  view()
  








