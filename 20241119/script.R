library(Seurat)
library(tidyverse)
library(patchwork)
library(biomaRt)
library(UCell)


# Functions ---------------------------------------------------------------

named_group_split <- function(.tbl, ...) {
  grouped <- group_by(.tbl, ...)
  names <- rlang::inject(paste(!!!group_keys(grouped), sep = " / "))
  
  grouped %>% 
    group_split() %>% 
    rlang::set_names(names)
}



#Read in seurat object
scrna <- readRDS("SeuratTest.RDS")
#set Ident to celltype so it's the default group when we plot
Idents(scrna) <- 'predicted.celltype'


# Manual genelist ---------------------------------------------------------


CAP1genes <- list('CAP1'=c("Gpihbp1",'Sema3c','Kit','Cd93','Glp1r','Hmcn1','Cxcl12','Cd74','Vegfa'))

#Use Seurats function AddModuleScore function
  
scrna <- AddModuleScore(object=scrna,features =CAP1genes,name='CAP1_')

#We have a new column in the meta table CAP1_1
names(scrna[[]])

#We can plot the score 

VlnPlot(object = scrna,features = 'CAP1_1')


scrna.sub <- subset(scrna, downsample=200)
scrna.sub <- AddModuleScore(object=scrna.sub,features =CAP1genes,name='CAP1_sub')

FeatureScatter(scrna.sub,feature1 = 'CAP1_1',feature2 = 'CAP1_sub1')


# Use Ucell https://github.com/carmonalab/UCell ---------------------------------------------------------------


#Remeber features is a list and not a vector!!!
scrna <- AddModuleScore_UCell(scrna, features =CAP1genes)

# A new column appears in Meta data with "_UCell suffix
names(scrna[[]])


VlnPlot(object = scrna,features = 'CAP1_UCell')

#Similar but different.. UCELL ranges from 0-1. 
FeatureScatter(scrna,feature1 = 'CAP1_1',feature2 = 'CAP1_UCell')

scrna.sub <- subset(scrna, downsample=200)
scrna.sub <- AddModuleScore_UCell(scrna.sub,features =list('CAP1_sub' = CAP1genes[[1]]))

#UCell does not change when subsetted. 

FeatureScatter(scrna.sub,feature1 = 'CAP1_UCell',feature2 = 'CAP1_sub_UCell')


# Another way to generate a gene list -------------------------------------


geneset2 <- list('ribo_genes' = grep("^Rps|^Rpl",rownames(scrna),value = T))
scrna <- AddModuleScore_UCell(scrna, features =geneset2)

VlnPlot(object = scrna,features = 'ribo_genes_UCell')



# Retrieve a genelist from the web -----------------------------------------


#this list is human genes so we need to convert to mouse. 

hs.stressgenes <- read_table('https://raw.githubusercontent.com/alexmascension/revisit_reynolds_fb/master/papers_genes_bad_quality/stress_genes.txt',col_names = 'gene') |> pull(gene)



#From JAX we download a table of human and mouse genes, table is organized as rows of mouse/genes. So we have to split the table ny organism, then join the 2 tables. 
mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t") |>
  mutate(Common.Organism.Name = gsub(', laboratory', '',Common.Organism.Name)) |> #The mouse field has laboratory, let's remove it
  dplyr::select(DB.Class.Key,Common.Organism.Name,Symbol) |> #We only need a few cols of data
  group_split(Common.Organism.Name) |> #This will split the data.frame in into a list of dataframes one per organism
  reduce(inner_join,by="DB.Class.Key") |> #This will join all lsu elements.  
  rename(human_symbol = Symbol.x, #CVlean up the names
         mouse_symbol = Symbol.y)

#We can use this table to convert the stressgenes

mm.stressgenes <-  mouse_human_genes |> filter(human_symbol %in% hs.stressgenes) |> pull(mouse_symbol)

scrna <- AddModuleScore_UCell(scrna, features = list('stress' = mm.stressgenes))
VlnPlot(object = scrna,features = 'stress_UCell')



# Create a gene list From FindallMarkers -------------------------------------------


fam.res <- FindAllMarkers(scrna,only.pos = T) |>
  mutate(pct.diff = pct.1-pct.2) |> 
  filter(cluster %in%  c('Cap1','Cap2','Arterial_endothelium','Venous_endothelium','Lymphatic_endothelium'))

cellmarkers <- 
  fam.res |>
  group_by(cluster) |> 
  top_n(12,pct.diff) |> 
  named_group_split(cluster) |>
  map(\(d)pull(d,gene))



scrna <- AddModuleScore_UCell(scrna, features = cellmarkers)
signature.names <- paste0(names(cellmarkers), "_UCell")

VlnPlot(object = scrna,features = signature.names,stack = T,flip = T) + NoLegend()



# Get Genes From GO Term --------------------------------------------------


go_term <- "GO:0001935"  # endothelial cell proliferation

# Connect to the Ensembl BioMart database
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")  # 

# Retrieve the list of genes associated with the GO term
mart_results <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "go_id", "name_1006","mgi_symbol"), 
  filters = "go", 
  values = go_term, 
  mart = mart
) 

genelist4 <- mart_results |> dplyr::select(mgi_symbol) |> distinct() |> pull(mgi_symbol)

scrna <- AddModuleScore_UCell(scrna, features = list('EndoProlif' = genelist4))

VlnPlot(object = scrna,features = "EndoProlif_UCell",pt.size = 0) + NoLegend()

