## Go term enrichment

##Authors:
# Fran Romero Campero
# Ana Belén Romero Losada
# Pedro de los Reyes Rodríguez

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("clusterProfiler")
# BiocManager::install("org.At.tair.db")

library(clusterProfiler)
library(org.At.tair.db)

## Auxiliary function to compute enrichments for GO table
compute.enrichments <- function(gene.ratios, bg.ratios)
{
  gene.ratios.eval <- sapply(parse(text=gene.ratios),FUN = eval)
  bg.ratios.eval <- sapply(parse(text=bg.ratios),FUN = eval)
  enrichments <- round(x=gene.ratios.eval/bg.ratios.eval,digits = 2)
  enrichments.text <- paste(enrichments, " (", gene.ratios, "; ", bg.ratios, ")",sep="")
  
  return(enrichments.text)  
}

atha.universe <- unique(select(org.At.tair.db,columns = c("GO"),keys=keys(org.At.tair.db,keytype = "TAIR"))[["TAIR"]])
length(atha.universe)

target.genes <- read.table(file="../tablas/targets_and_upregulated_genes.txt")
target.genes <- as.vector(target.genes[[1]])
length(target.genes)


enrich.go <- enrichGO(gene          = target.genes ,
                      universe      = atha.universe,
                      OrgDb         = org.At.tair.db,
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = FALSE,
                      keyType = "TAIR")

## Generate ouput table
enrich.go.result <- as.data.frame(enrich.go)

## GO term Description P-value Q-value Enrichment (SetRatio, BgRatio) Genes
go.term.enrichments <- compute.enrichments(gene.ratios = enrich.go.result$GeneRatio,
                                           bg.ratios = enrich.go.result$BgRatio)

go.result.table <- data.frame(enrich.go.result$ID, enrich.go.result$Description,
                              enrich.go.result$pvalue, enrich.go.result$qvalue,
                              go.term.enrichments, 
                              gsub(pattern = "/",replacement = " ",x = enrich.go.result$geneID),
                              stringsAsFactors = FALSE)

colnames(go.result.table) <- c("GO ID", "Description", "p-value", "q-value",
                               "Enrichment (Target Ratio; BG Ration)","Genes")

# write.table(go.result.table, file = "go_terms_table.txt", sep = "\t", row.names = F,
            # quote = F)

goplot(enrich.go,showCategory = 10)
barplot(enrich.go,drop=TRUE,showCategory = 40)
emapplot(enrich.go)
cnetplot(enrich.go)
dotplot(enrich.go, showCategory =20)

go.result.table[73,]
go.result.table[74,]

## Genero na imagen tiff grande con alta resolución para que se vean bien
## todas las categorías

tiff(file="barplo_goterms.tiff",
     width     = 8,
     height    = 10,
     units     = "in",
     res       = 300
)

barplot(enrich.go,drop=TRUE,showCategory = 30)

dev.off()

### Plotting specific GO terms #####

goplot(enrich.go, term=c("GO:0007623", "GO:0048511"))

plotGOgraph(enrich.go)

key.terms <- c("GO:0007623", "GO:0048511", "GO:0010218", "GO:0010114",
               "GO:0009637", "GO:0010224", "GO:0009639", "GO:0009411",
               "GO:0048574", "GO:0042752", "GO:0048571")

key.terms <- c("GO:0007623","GO:0048511","GO:0071482")

keep <- go.result.table$`GO ID`[!go.result.table$`GO ID` %in% key.terms]
circadian.go <- dropGO(enrich.go, term = keep )

goplot(circadian.go)
tiff('cnetplot__nolabels.tiff', units="in", width=4, height=4, res=300, compression = 'lzw')
cnetplot(circadian.go, cex_label_gene=0.4,
         cex_label_category= 0.5, node_label="category")
dev.off()
emapplot(circadian.go)

## Reducing the redundancy of GO using rrvgo package

# get the similarity matrix between terms
simMatrix <- calculateSimMatrix(go.result.table$`GO ID`,
                                orgdb="org.At.tair.db",
                                ont="BP",
                                method="Rel")

# group terms based on similarity
scores <- setNames(-log10(go.result.table$`q-value`), go.result.table$`GO ID`)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.At.tair.db")

# Plot similarity matrix as a heatmap
heatmapPlot(simMatrix,
            reducedTerms,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=6)

# Plot GO terms as scattered points
scatterPlot(simMatrix, reducedTerms)

# Treemaps are space-filling visualization of hierarchical structures
treemapPlot(reducedTerms)

# wordcloud
wordcloudPlot(reducedTerms, min.freq=1, colors="black")
