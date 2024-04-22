library(Seurat)
library(tidyverse)
library(UCell)
library(SeuratData)
library(SCP)
library(patchwork)

# Control structures if else ------------------------------------------------------

a <- 10
b <- 20

if (b > a) { #These brackets define a code block
  print("b is greater than a")
} 




a <- 20
b <- 10

if (b > a) {
  print("b is greater than a")
} else {
  print("a is greater than b")
}



#What if there equal>?
a <- 20
b <- 20

if (b > a) {
  print("b is greater than a")
} else {
  print("a is greater than b")
}

#We can use if else

if (b > a) {
    print("b is greater than a")
  } else if (a== b ) {
    print("a is equal to b")
  } else {
  print("a is greater than b")
}

#######################


#paste function

c1 <- 'AT1'
c2 <- 'AT2'

paste(c1,c2)

paste(c1,c2,sep=',')

cells <- c('AT1','AT2','CAP1')

#Sep will return a vector same size as largest input vector. 
paste('Celltype',cells,sep=':')

#Collapse will reurn a string 
paste('Celltype',cells,collapse =':')

#paste0 is a shirtcut function, mostly used for concatenation operations. 
paste0()

paste0('The first cell is ',c1,',the second cell is ',c2,' in the data')

#We can use glue instrad for better concateantion, either using hte library glue or str_glue from tidyverse
#We create a statment in which block our variables using {}
str_glue('The first cell is {c1},the second cell is {c2} in the data')

#What if need to use "{}" in your string? Just use double {{ or }}

str_glue('This is not a {{varible}} but this is: {c1}')




# Functions ---------------------------------------------------------------


my_function <- function() { # create a function with the name my_function
  print("Hello World!")
} 

my_function()


#Function with arguments. 
addme <- function(x,y){
  x+y
}

addme(2,4)




# Putting it all togather -------------------------------------------------
data("pbmc_small")


scrna <- pbmc_small
scrna
scrna$seurat_clusters <- scrna$RNA_snn_res.0.8 



d1 <- CellDimPlot(scrna, group.by = "seurat_clusters", reduction = "UMAP",theme_use = 'theme_blank', label = T,label_insitu = T,legend.position = 'none')
f1 <- FeatureDimPlot(scrna,'SPON2' ,theme_use = 'theme_blank',ncol=1)
v1 <- FeatureStatPlot(scrna, stat.by ="SPON2" , group.by = "seurat_clusters",legend.position ='none')


(d1|f1)/v1

(d1|f1)/v1 + plot_layout(heights = c(5,1))

pw <- (d1|f1)/v1 + plot_layout(heights = c(5,1))

ggsave('SPON2_plot.png',pw,height=10,width=10)


#This is great but what if I wanted to plot another gene! or 100. 


myGenePlot <- function(gene,savefile=FALSE){
  d1 <- CellDimPlot(scrna, group.by = "seurat_clusters", reduction = "UMAP",theme_use = 'theme_blank', label = T,label_insitu = T,legend.position = 'none')
  f1 <- FeatureDimPlot(scrna,gene ,theme_use = 'theme_blank',ncol=1)
  v1 <- FeatureStatPlot(scrna, stat.by =gene , group.by = "seurat_clusters",legend.position ='none')
  pw <- (d1|f1)/v1 + plot_layout(heights = c(5,1))
  if(savefile){
    ggsave(str_glue('{gene}_plot.png'),pw,height=10,width=10)
  }else{
    pw
  }
  
}

myGenePlot('PTPN22')
myGenePlot('CCR7')
myGenePlot('SIT1',savefile = TRUE)
  
  
  








