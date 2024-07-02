library(Seurat)
library(tidyverse)
library(patchwork)
library(SCP)


#More rehash on the purrr library, which is the map family of functions 

alist <- list('name1' = 'value1',
              'name2' = 'value2',
              'name3' = 'value3'
)

alist
# We'll use map to loop over
map(alist,function(element_value){
  print(str_glue("The current element value is {element_value}"))
})
  
  
#Map always returns a list that is the same size as in the input list, in the above case we just printing and we don't need to return a new list, We can use
# the walk function for this. 

walk(alist,function(element_value){
  print(str_glue("The current element value is {element_value}"))
})
  
#We can use the shorthand for function which is "\" 

walk(alist,\(element_value){
  print(str_glue("The current element value is {element_value}"))
})


#What if we need is use the name of the list? We can use the imap or iwalk functions. 
#The first argument will always be the value and the 2nd will be the name in the list. 
iwalk(alist,\(element_value, element_name){
  print(str_glue("The current element is {element_name} and the value is {element_value}"))
})

## Download the Seurat object from Dropbox, see my email

scrna <- readRDS("SeuratTest.RDS")

#Seurat version
DimPlot(scrna,group.by = 'predicted.celltype')

#Let's use SCP package to make nicer plots
d1 <- CellDimPlot(scrna, group.by = "predicted.celltype", reduction = "UMAP",theme_use = 'theme_blank', 
                  label = T,legend.position = 'none', label_insitu = T,label_repulsion = T,palette = 'Categorical.12')
f1 <- FeatureDimPlot(scrna,features='Car4',reduction = "UMAP",theme_use = 'theme_blank')
v1 <- FeatureStatPlot(srt = scrna, group.by = 'var_celltype',stat.by = c("Car4"), add_box = TRUE)
(d1+f1)/v1


# Let's reorder the celltypes 
scrna$predicted.celltype <- factor(scrna$predicted.celltype, c("Cap1","Cap1_2",  "Cap2", 'Cap_postflu',"Arterial_endothelium","Venous_endothelium","Lymphatic_endothelium"))

#Let's remake the violinplot now

v1 <- FeatureStatPlot(srt = scrna, group.by = 'predicted.celltype',stat.by = c("Car4"), add_box = TRUE)
v1

#and now a 3 panel figure. 
(d1+f1)/v1

#What if we want plot many genes. First let's turn this into a fucntion. We need a seaurt object and a genename as args. 

myGenePlot <- function(object,gene){
  d1 <- CellDimPlot(object, group.by = "predicted.celltype", reduction = "UMAP",theme_use = 'theme_blank', label = T,legend.position = 'none', label_insitu = T,label_repulsion = T,palette = 'Categorical.12')
  f1 <- FeatureDimPlot(object,features=gene,reduction = "UMAP",theme_use = 'theme_blank')
  v1 <- FeatureStatPlot(srt = object, group.by = 'predicted.celltype',stat.by = gene, add_box = TRUE)
  (d1+f1)/v1
}

#let's test it. 
myGenePlot(object=scrna,gene='Atf3')

#Let's do a few genes 

genestoplot <- c('Car4','Atf3','Prox1','Gpihbp1')


#This output the plots to rstudio plots window
map(genestoplot,\(gene){
  myGenePlot(object=scrna,gene=gene)
})
    
#We can now save the output to a list. 
plots <- map(genestoplot,\(gene){
  myGenePlot(object=scrna,gene=gene)
})

plots[[1]]
plots[[2]]

#It would be very convenient if the set the names in the list
# to the gene. set_names will aloow us to do that. 

plots <- map(genestoplot,\(gene){
  myGenePlot(object=scrna,gene=gene)
}) |> set_names(genestoplot)


names(plots)

plots$Prox1

#Let's save these figures to a file, one file per gene. 
#Since we are only saving files, we don't need a return object so we will use walk. 
#Are list cleverly uses the gene as the name, so we then can use iwalk to get gene name for the 
# file name as well. First we create a new dir to store are plots, we will use a variable for the plot dir
#so we can use it later. 

plotdir <- 'plots'
dir.create(plotdir)

iwalk(plots,\(plot,genename){
  ggsave(filename=str_glue('{plotdir}/{genename}.png'),plot = plot,width=12,height = 10)  
})




iwalk(genestoplot,\(gene){
  plot <- myGenePlot(object=scrna,gene=gene)
  ggsave(filename=str_glue('{plotdir}/{gene}.png'),plot = plot,width=12,height = 10) 
})





