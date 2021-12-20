##Author: Pedro de los Reyes Rodríguez / Francisco J. Romero-Campero
##Email: pedderod@alum.us.es / fran@us.es
##Date: November 2016

arguments <- commandArgs(trailingOnly = TRUE)

if(length(arguments) == 0)
{
  print("Usage: Rscript go_enrichment.R text_file_with_gene_set ouput_file_for_image output_file_for_table species")
}

## First argument consists of a text file with the gene set
tf.targets <- arguments[1] 
## Biological process is the only analyzed ontology
ontology.type <- "BP"
## Second argument consists of the file name to store the graphical 
## representation of the enrichment as a DAG
tf.image.GO <- arguments[2] 
## The third argument consists of the file name to store the table summing up
## the significantive GO terms detected
tf.table.GO <- arguments[3] 
## The fourth argument consists of the species id to analyze 
## (ota - Ostreoccocus tauri)
## (cre - Chlamydomonas reinhardtii)
## (atha - Arabidopsis thaliana)
species <- arguments[4] #"cre"

## Loading library topGO
library("topGO")
## The p-value threshold is set to 0.05
p.value <- 0.05

##Reading in gene set file
target.genes <- read.table(tf.targets)
target.genes <- as.vector(target.genes[[1]])


##Loading GO annotation
if(species == "atha")
{
  map.file <- "../data/GO_terms/athaliana_go/athaliana_go.map"
} else if(species == "ota")
{
  map.file <- "../data/GO_terms/otauri_go/otauri_GO_v2.map"
  
} else if(species == "cre")
{
  map.file <- "../data/GO_terms/creinhardtii_go/creinhardtii_go.map"
}

geneID2GO <- readMappings(file=map.file)

##Setting the background
gene.names <- attributes(geneID2GO)[[1]]
gene.background <- rep(1, length(gene.names))
names(gene.background) <- gene.names

gene.background[target.genes] <- 0

##Setting target genes (those assigned a 0)
ath.gene.selec <- function(gene.list)
{
  return(gene.list == 0)
}

##Generating topGOdata object
sampleGOdata <- new("topGOdata",
                    description = "Arabidopsis session", ontology = ontology.type,
                    allGenes = gene.background, geneSel = ath.gene.selec,
                    nodeSize = 10,
                    annot = annFUN.gene2GO, gene2GO = geneID2GO)

resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")

##Extracting the top 100 most significative GO terms
allRes <- GenTable(sampleGOdata, classicFisher = resultFisher,
                   ranksOf = "classicFisher", topNodes = 100,numChar=100)


##Extracting p-value column
all.pvalues <- allRes[["classicFisher"]]

##Converting p-values to numeric
corrected.all.pvalues <- vector(length=length(all.pvalues), mode="numeric")

for (i in 1:length(all.pvalues))
{
  if (is.na(as.numeric(all.pvalues[i])))
  {
    corrected.all.pvalues[i] <- 0
  }else 
  {
    corrected.all.pvalues[i] <- as.numeric(all.pvalues[i])
  }
}

##Keeping only GO.terms with a p.value < 0.05
filter.Res <- allRes[corrected.all.pvalues < p.value, ]


res.with.genes <- data.frame(filter.Res,genes=vector(length=nrow(filter.Res),mode="character"))


##Extracting GO terms
target.genes.annotation <- geneID2GO[target.genes]
significative.GO.terms <- as.vector(filter.Res[["GO.ID"]])

##Getting genes associated to each GO
offspring <- as.list(GOBPOFFSPRING)

genes.for.each.GO <- vector(length=length(significative.GO.terms),mode="character")

for (i in 1:length(significative.GO.terms))
{
  GO.genes <- vector(mode="character")
  k <- 1
  go.offspring <- offspring[[significative.GO.terms[i]]]

  for (j in 1:length(target.genes.annotation))
    {
      if(length(intersect(go.offspring,target.genes.annotation[[j]]))!=0)
      {
        GO.genes[k] <- names(target.genes.annotation[j])
        k <- k+1
      } 
    }

  GO.genes.string <- paste(unique(GO.genes), collapse=" ")
  
  genes.for.each.GO[i] <- GO.genes.string
  print(i)
}


res.with.genes[["genes"]] <- genes.for.each.GO

##Generating image with DAG
png(file=tf.image.GO,
    width     = 10,
    height    = 10,
    units     = "in",
    res       = 600
)

showSigOfNodes(sampleGOdata, score(resultFisher), firstSigNodes = 50, useInfo = 'all')

dev.off()

##Writing table with significative GO terms
write.table(res.with.genes,file=tf.table.GO, row.names=FALSE)

