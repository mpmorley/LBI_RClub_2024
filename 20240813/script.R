library(Seurat)
library(tidyverse)
library(patchwork)
library(SCP)
library(clusterProfiler)


#We'll create a function that takes a data.frame splits into a named list of data frames one per item of the 
# group variable.

named_group_split <- function(.tbl, ...) {
  grouped <- group_by(.tbl, ...)
  names <- rlang::inject(paste(!!!group_keys(grouped), sep = " / "))
  
  grouped %>% 
    group_split() %>% 
    rlang::set_names(names)
}

#Read in the data from 
scrna <- readRDS("SeuratTest.RDS")

#Here we're using var_celltype, your celltype annotation name could be different. THis is a simple violinPlot using seurat. 
VlnPlot(scrna,'Car4',group.by = 'var_celltype')

#Here's the same using the SCP package
FeatureStatPlot(srt = scrna, group.by = 'var_celltype',stat.by = c("Car4"))


#I don't like the order of celltytpes in this plot, to change them we need to make a factor and define the levels. 
# Let's reorder the celltypes 
scrna$var_celltype <- factor(scrna$var_celltype, c("Cap1","Cap1_2",  "Cap2", 'Cap_postflu',"Arterial_endothelium","Venous_endothelium","Lymphatic_endothelium"))

#Let's replot. 
FeatureStatPlot(srt = scrna, group.by = 'var_celltype',stat.by = c("Car4"))


#Seurat using the special variable called ident to define a cell, this can be cluster, sample name etc. Certain Seurat functions will use this Ident
# We are going to set it now to var_celltype. 
Idents(scrna) <- "var_celltype"
DefaultAssay(scrna) <- 'RNA'

DiffExp <- FindAllMarkers(scrna) |> mutate(pctdiff = pct.1 - pct.2)




d.tmp <- DiffExp |> filter(cluster=='Cap1' & p_val_adj < 0.05)



enrichGO(
  gene=d.tmp$gene,
  OrgDb         = org.Mm.eg.db::org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05,
  readable = T
)

#What happened? Be defaul enrichGO expects entrezID as input not gene symbol. 
#We can either convert symbols to entrezID or we can tell enrichGO to use symbol


enrichGO(
  gene=d.tmp$gene,
  OrgDb         = org.Mm.eg.db::org.Mm.eg.db,
  keyType = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05,
  readable = T
)


#So how do we go about doing this for every celltype?

#Let's first start by re-visiting the group_by function. Remeber cluster is the name in FindallMarkers output
DiffExp |> group_by(cluster)

#Let's find the max pct diff gene for each Celltype
DiffExp |> group_by(cluster) |> slice_max(pctdiff,n=1)

#In fact we can pipe this function into a plot command as well, our plot command however takes a vector of genes, not a table. 
# We can use the dplyr pull to get vector of genes. 


DiffExp |> group_by(cluster) |> slice_max(pctdiff,n=1) |> pull(gene)

#Now we can pipe, however the stat.by is not the first argument, so we have to direct our pipe to the proper argument with a placeholder, 
# We use the "_" as the palceholder. 
DiffExp |> group_by(cluster) |> slice_max(pctdiff,n=1) |> pull(gene) |> 
  FeatureStatPlot(srt = scrna, group.by = 'var_celltype',stat.by = _)

#If we're using magrittr pipes we would use the "." 

DiffExp %>% group_by(cluster) %>% slice_max(pctdiff,n=1) %>%pull(gene) %>%
  FeatureStatPlot(srt = scrna, group.by = 'var_celltype',stat.by = .)


#This is ugly but we can use the stack parameter. 
DiffExp |> group_by(cluster) |> slice_max(pctdiff,n=1) |> pull(gene) |> 
  FeatureStatPlot(srt = scrna, group.by = 'var_celltype',stat.by = _, stack=T)


#Back to performing enrichment across all cells. 


#The group_split function will return back a list of data.frames for every unique level in the group variable, in our case celltype. However it's not a named list. 
DiffExp |> group_by(cluster) |> group_split()

#using the above function we'll now get a named list. 
DiffExp |> named_group_split(cluster)

#We can now use map to "loop" over the data for each celltype. Let's do some trivial like taking the first row of each subset

DiffExp |> named_group_split(cluster) |> map(\(d){
 head(d,n=1) 
})



#Let's now put this all togather. 


go.results <- 
DiffExp |> filter(p_val_adj < 0.05) |> 
  named_group_split(cluster) |>
  map(\(d){
    enrichGO(
      gene=d$gene,
      OrgDb = org.Mm.eg.db::org.Mm.eg.db,
      keyType = "SYMBOL",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.01,
      qvalueCutoff  = 0.05,
      readable = T
    )
  })
  
#We can save the results to an excel sheet. just recall that the results are 
# stored in the results slot of the enrichGO object 


go.results$Cap1@result  
  
#What we need to do then is extract out the results from each object, we can use map again

map(go.results,\(x){
    x@result
    }) |> 
writexl::write_xlsx('GOresults.xlsx')
  

#We can also make plots in bulk as well. 

barplots.list <- map(go.results,\(x){
  barplot(x, showCategory=6,label_format = 1) 
})

barplots.list <- imap(go.results,\(x,name){
  barplot(x, showCategory=6,label_format = 1) + ggtitle(name)
})

#We can plot them as one plot.. 
wrap_plots(barplots.list)


#Or we can do one plot at a time. We'll use iwalk this time

#first let's create an output dir for this plots,take note of the 

plotdir<-'plots'
dir.create(plots)

#Take note of the "/" in the name to define dir path

iwalk(barplots.list,\(p,name){
  ggsave(str_glue("{plotdir}/{name}_barplot.png"),p)
})

#Let's make a function that makes 2 plots 


MakePlots <- function(ego,name){
  p1 <- barplot(ego, showCategory=6,label_format = 1) + ggtitle(name)
  p2 <- cnetplot(ego,showCategory = 5)
  p1|p2
}

#We can plot and save it on one step. 
iwalk(go.results,\(x,name){
  p <- MakePlots(x,name)
  ggsave(str_glue("{plotdir}/{name}_multiplot.png"),p,bg = 'white')
})


