library(Seurat)
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(msigdbr)

scrna <- readRDS("SeuratTest.RDS")




#Make a mouse/human ortholog lookup table. We'll use this later. 
mh_data <- read_tsv("https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt") |> 
  group_split(`Common Organism Name`) |> 
  map(~select(.x,`DB Class Key`,Symbol,`EntrezGene ID`)) |>
  list(l=_) |>
  with(inner_join(x=l[[1]],y=l[[2]],by = "DB Class Key",suffix = c('_human','_mouse')))


#Many of the enrichments require ENTREZ id and not gene symbol, we'll create a table of gene symbol and ENTREZ ID. We use
# the bitr helper function from ClusterProfiler. 

eg = bitr(rownames(scrna), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Mm.eg.db::org.Mm.eg.db)

# Run Diffexpressed --------------------------------------


#We need to set an extra parameter, min.pct=0. For GSEA we need the results for all genes in the data set. 

Idents(scrna) <- 'celltype'
dge <- FindMarkers(scrna, ident.1 = "Cdc42KO_d14", ident.2 = 'Ctrl', group.by = 'sample', subset.ident = 'Cap1',min.pct=0) |>
  rownames_to_column('gene') |>
  mutate(pctdiff= pct.1 - pct.2) |>
  inner_join(x=_,y=eg,by=c('gene'='SYMBOL'))


# over-representation analysis --------------------------------------------

#For this we need to select a set of genes of interest. There's numerous ways to do this

# Adj P and FC
dge.genes <- dge |> filter(p_val_adj < 0.05 & abs(avg_log2FC) > log2(1.5) ) |> pull(ENTREZID)
length(dge.genes)


er.BP <- enrichGO(
  gene=dge.genes,
  OrgDb         = org.Mm.eg.db::org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05,
  readable = T
)

er.BP_universe <- enrichGO(
  gene=dge.genes,
  universe = eg$ENTREZID,
  OrgDb         = org.Mm.eg.db::org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05,
  readable = T
)

#Compare 
p1 <- barplot(er.BP, showCategory=10,label_format = 1) 
p2 <- barplot(er.BP_universe, showCategory=10,label_format = 1) 

p1|p2



# We can use pctdiff as a cutoff as well. 
dge.genes.pctdiff <- dge |> filter(p_val_adj < 0.05 & abs(pctdiff) > .15) |> pull(ENTREZID)
length(dge.genes.pctdiff)

er.BP_universe.pctdiff <- enrichGO(
  gene=dge.genes.pctdiff,
  universe = eg$ENTREZID,
  OrgDb         = org.Mm.eg.db::org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05,
  readable = T
)

p3 <- barplot(er.BP_universe.pctdiff, showCategory=10,label_format = 1) 
p2|p3


#For GSEA we don't pre select the genes but rather we rank them by a metric, we'll use FC.  
# We need to create a named vector, the names being ENTREZID and the value being the FC sorted in descending order. 

ranks <- deframe(dge |> arrange(-avg_log2FC) |> select(ENTREZID,avg_log2FC))


gseGO <- gseGO(geneList     = ranks, 
        OrgDb        = org.Mm.eg.db::org.Mm.eg.db,
        ont          = 'BP',
        minGSSize    = 50,
        maxGSSize    = 500,
        pvalueCutoff = 1,
        verbose      = FALSE)



#There isn't a barplot for the gseaResult, we'll use the dotplot. 
d1 <- dotplot(er.BP_universe, showCategory=10,label_format = 1) 
d2 <- dotplot(gseGO, showCategory=10,label_format = 1) 

d1|d2



# use MsigDB -------------------------------------------------------------
#https://www.gsea-msigdb.org/gsea/msigdb

#We can see what is avaible. 
msigdbr_collections(db_species='MM')

#Let's do hallmark gene set. 

gs.MH <- msigdbr(db_species='MM', species = "Mus musculus", collection = "MH") |>
  dplyr::select(gs_name, ncbi_gene)

gsea>MH <- GSEA(ranks, TERM2GENE = gs.MH)



msiglist <- c('M1','M2','M3','M5','MH')
msig_genesets <- map(msiglist,~msigdbr(db_species='MM',species = "Mus musculus", collection = .x) |>
                       dplyr::select(gs_name, ncbi_gene)
) |> setNames(msiglist)

gse.list <- map(msig_genesets,\(gs){
  GSEA(ranks, TERM2GENE = gs)
})


# GSE on KEGG and KEGG module  ---------------------------------------------------------------------

gse.list[['KEGG']] <- gseKEGG(
  geneList     = ranks,
  organism     = 'mmu',
  minGSSize    = 120,
  pvalueCutoff = 1,
  verbose      = FALSE)

gse.list[['MKEGG']] <- gseMKEGG(geneList = ranks,
                                organism = 'mmu',
                                pvalueCutoff = 1)

# GSE on Wikipathways -----------------------------------------------------

gse.list[['WP']]  <- gseWP(ranks, organism = "Mus musculus")



# Download a file from Enrichr and perform GSEA --------------------------------


panther_gs <- read.gmt("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=Panther_2016") |>
  inner_join(mh_data,by=c('gene'='Symbol_human')) |> select(term,`EntrezGene ID_mouse`)
gse.list[['PANTHER']]  <- GSEA(ranks, TERM2GENE = panther_gs)

gse.list[['PANTHER']]@result


#Convert the ENTREZ ID back to gene symbols. 
gse.list <- map(gse.list, ~setReadable(.x,OrgDb = org.Mm.eg.db::org.Mm.eg.db, keyType="ENTREZID"))

map(gse.list,~.x@result) |>
  writexl::write_xlsx(path='GSE_results.xlsx')




