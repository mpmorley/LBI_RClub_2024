#############################################################################
# MPM 20240604 R club Script for Gene Set Enrichment 
# 
# https://yulab-smu.top/biomedical-knowledge-mining-book/index.html
#

# Install clusterprofiler from BioConductor. 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
BiocManager::install("enrichplot")
BiocManager::install("pathview")


# We also need to install species specific Gene Annotation database
# We'll install the human one,
# see http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
BiocManager::install("org.Hs.eg.db")


#Load some libraries 

library(tidyverse)
library(clusterProfiler)

#Let's read in some Diff Expression results, we're going to sort by Foldchange as well,
# we need to do this for later. 


diff <- read.csv('DiffExRresults.csv') |>
  mutate(ENTREZID=as.character(ENTREZID)) |>
  arrange(-fc)



diff |> select(SYMBOL)

#PROBLEM!!! 

#Look at the output from loading this library, it say it masks select from dplyr
library(org.Hs.eg.db)

#Now we try using dplyr select and it fails, it's using the AnnotationDbi version
diff |>select(SYMBOL)

#We can fix it a few ways, we can define the full name 

diff |> dplyr::select(SYMBOL)

#we can remask 

select <- dplyr::select

#Or we can not load org.Hs.eg.db and then we can define the full name were we need it, this what I do 
#Since this package imports so many other pacakges it's best to restart R.

org.Hs.eg.db::org.Hs.eg.db




# Let's do some analysis --------------------------------------------------

#First thing is look at your gene annotation, Do you have offical symbols
# refseq iD, Ensembl IDs, EntrezID etc. 


gene.df <- bitr(diff$ENTREZID, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db::org.Hs.eg.db)

gene.df




#We start basic Gene Ontology over representation analysis 
#For this we need create a vector of genes we'd like to test.. 
# Say Upregulated genes with a FC > 2

genes <- diff  |>
  filter(adj.P.Val < 0.05 & fc > 2) |>
  pull(ENTREZID)

ego<-  enrichGO(
  gene=genes,
  OrgDb         = org.Hs.eg.db::org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05,
  readable = T
)


#Though we really should provide a background set of genes, the first test 
#tested all of the genes in the genome, however we really should test the 
# gene expressed in our experiment. 


ego.up <- enrichGO(gene       = genes,
                   universe      = diff$ENTREZID  ,
                   OrgDb         = org.Hs.eg.db::org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE
)





#We can use symbol as well, let's also use a pipe here too. 

ego.up <- diff  |>
  filter(adj.P.Val < 0.05 & fc > 2) |>
  pull(SYMBOL) |>
  enrichGO(
    universe      = diff$SYMBOL,
    keyType       = 'SYMBOL',
    OrgDb         = org.Hs.eg.db::org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.01,
    qvalueCutoff  = 0.05,
    readable=F
  )



#Let's do the down reg genes 

ego.down <- diff  |>
  filter(adj.P.Val < 0.05 & fc < -2) |>
  pull(SYMBOL) |>
  enrichGO(
    universe      = diff$SYMBOL,
    keyType       = 'SYMBOL',
    OrgDb         = org.Hs.eg.db::org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.01,
    qvalueCutoff  = 0.05,
    readable=F
  )



ego.both <- diff  |> filter(adj.P.Val < 0.05 & abs(fc) > 2) |>
  pull(SYMBOL) |>
  enrichGO(
    universe      = diff$SYMBOL,
    keyType       = 'SYMBOL',
    OrgDb         = org.Hs.eg.db::org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.01,
    qvalueCutoff  = 0.05,
    readable=F
  )





# Wikipathways ------------------------------------------------------------


WP <- enrichWP(genes,
               organism = "Homo sapiens"
               ) 
WP@result

# Gene Set Enrichment -----------------------------------------------------


#We need to create a "named" vector of Foldchanges. Remeber GSEA takes as input all the genes, we don't filter
geneList <- setNames(diff$fc,diff$ENTREZID)


gse <- gseGO(geneList     = geneList,
             OrgDb        = org.Hs.eg.db::org.Hs.eg.db,
             ont          = "BP",
             nPerm        = 1000,
             minGSSize    = 100,
             maxGSSize    = 500,
             pvalueCutoff = 0.05,
             verbose      = FALSE)




# Results tables ----------------------------------------------------------

#All these functions return a "enrichRresult" object

view(WP@result)

WP.symbols <- setReadable(WP, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType="ENTREZID")

view(WP.symbols@result)

#Write out results for a file. 
write_csv(WP.symbols@result,file='WP_results.csv')



# Plot the results  ------------------------------------------------------
library(enrichplot)
library(patchwork)
#Barplot and dotplots are very common ways to show GO data
barplot(ego.up, showCategory=20)

p1 <- barplot(ego.up, showCategory=10,label_format = 7) + ggtitle('Up Regulated Genes')
p2 <- barplot(ego.down, showCategory=10,label_format = 7) + ggtitle('Down Regulated Genes')
p1+p2



dotplot(ego.up, showCategory=30) 

#The Cnet plot might be useful for certain 


geneFC <- setNames(diff$fc,diff$SYMBOL)

cnetplot(ego.both,showCategory = 5)

#can use 
cnetplot(ego.both, foldChange=geneFC)

#Heatmap plot
heatplot(ego.both,showCategory = 20)
heatplot(ego.both,showCategory = 20,foldChange = geneFC)




ego.pair <- pairwise_termsim(ego.both)

treeplot(ego.pair)




upsetplot(ego.both) + theme(plot.margin = unit(c(.5, .5, .5, 3), 
                                            "inches"))


#You can use dplyr verbs and ggplot as well, such as filter 

ego.up.filter <- ego.up |> filter(ID %in% c('GO:0001501','GO:0022407'))

cnetplot(ego.up.filter)

ego.up.sort <- ego.up |> arrange(-Count)


#You can also ggplot
ggplot(ego.up.sort,showCategory = 20, 
       aes(Description,Count,fill=-log10(pvalue))) + 
  geom_bar(stat='Identity') + 
  coord_flip() + 
  scale_fill_viridis_b(option='C') + 
  theme_bw()



# Plot GSEA ---------------------------------------------------------------
#Let's look at results

view(gse@result)

gseaplot(gse, geneSetID = 1, by = "runningScore", title = gse$Description[1])


gseaplot2(gse, geneSetID = 15, title = gse$Description[15])


gseaplot2(gse, geneSetID = c(1,4,8))


gseaplot2(gse, geneSetID = c(1,4,8),pvalue_table = T,ES_geom = 'dot')




# KEGG analysis -----------------------------------------------------------


kk <- diff  |>
  filter(adj.P.Val < 0.05 & fc > 2) |>
  pull(ENTREZID) |>
  enrichKEGG(
    organism     = 'hsa',
    universe = diff$ENTREZID,
    pvalueCutoff = 0.05,
    
  )

view(kk@result)

#Plot the kegg results

browseKEGG(kk, kk.result$ID[1])

#Plot using pathview 
library(pathview)


pathview(gene.data  = geneList,
         pathway.id =  "hsa04340",
         species    = "hsa",
         limit      = list(gene=max(abs(geneList)), cpd=1))


