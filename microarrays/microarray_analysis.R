
###### --- Microarray analysis --- ####

# Authors
# Francisco J. Romero-Campero
# Pedro de los Reyes Rodríguez

# Load libraries
library(affy)
library(simpleaffy)
library(affyPLM)
library(limma)
library(annaffy)
library(moe430a.db) #it depends on the analyzed organism 
library(gplots)

## Read the .CEL files 

microarray.data <- ReadAffy(verbose=TRUE) 
microarray.data
cdfName(microarray.data) #getting info          

######### Visualization of microarray data #####
image(microarray.data[,1],col=rainbow(100))
image(microarray.data[,2],col=rainbow(100))
image(microarray.data[,3],col=rainbow(100))


######### Quality analysis #########################
quality.analysis <- qc(microarray.data)
plot(quality.analysis)

# Debido a que la beta-actina y GAPDH se expresan es la mayoría de los tipos de células
# y son genes relativamente largos, Affymetrix los usa como controles de la 
# calidad del ARN. Tres conjuntos de sondas están diseñados en 3 regiones de estos genes
# (5', mid (llamado M) y 3'). Intensidades similares para sus 3 regiones indican que las 
# transcripciones no fueron truncadas y están etiquetadas igualmente a lo largo de la 
# secuencia. Además, puesto que la degradación del ARN comienza desde el extremo 5' de 
# la molécula, es común que la intensidad del probeset (conjunto de sondas) en ese extremo 
# sea ligeramente inferior.

# Para un array de buena calidad, Affymetrix recomienda que la proporción 3 '/ 5'
# no exceda de:
#   - 3 para la beta-actina
# - 1,25 para GAPDH
# Estos valores se establecieron para las muestras humanas, por lo que la proporción 
# puede ser ligeramente superior para otras especies.


# Boxplot previsualization to see if the samples are comparables

boxplot(microarray.data,col=rainbow(14),las=2,ylab="Luminescencia")
hist(microarray.data,col=rainbow(14))

# Normalization

microarray.processed.data <- rma(microarray.data)

# Boxplot visualization post normalization

boxplot(microarray.processed.data,col=rainbow(8),las=2,ylab="Luminiscencia")
hist(microarray.processed.data,col=rainbow(8))


########### Average expression levels estimation ############

expression.level <- exprs(microarray.processed.data)
dim(expression.level)
head(expression.level)
sampleID <- c("Rb_nocis_1","Rb_nocis_2","Rb_nocis_3",
              "WT_nocis_1", "WT_nocis_2", "WT_nocis_3",
              "Rb_cis_1","Rb_cis_2","Rb_cis_3",
              "WT_cis_1","WT_cis_2","WT_cis_3")
colnames(expression.level) <- sampleID
head(expression.level)

Rb_nocis <- (expression.level[,"Rb_nocis_1"] + expression.level[,"Rb_nocis_2"]
             + expression.level[,"Rb_nocis_3"])/3

Rb_cis <- (expression.level[,"Rb_cis_1"] + expression.level[,"Rb_cis_2"]
             + expression.level[,"Rb_cis_3"])/3

WT_nocis <- (expression.level[,"WT_nocis_1"] + expression.level[,"WT_nocis_2"]
             + expression.level[,"WT_nocis_3"])/3

WT_cis <- (expression.level[,"WT_cis_1"] + expression.level[,"WT_cis_2"]
           + expression.level[,"WT_cis_3"])/3

probe.names <- names(Rb_nocis)
mean.expression <- matrix(c(Rb_nocis,WT_nocis,Rb_cis,WT_cis),ncol=4)
conditions.id <- c("Rb_nocis","WT_nocis","Rb_cis","WT_cis")
rownames(mean.expression) <- probe.names
colnames(mean.expression) <- conditions.id
head(mean.expression)


#######################################################################################
##############       Differential gene expression analysis            #################
#######################################################################################


# Comparations

plot(WT_nocis,Rb_nocis, col="grey",xlab="WT no tratadas",ylab="Rb-/- no tratadas",
     pch=19,cex=0.5)

plot(WT_nocis, Rb_nocis)

plot(WT_cis,Rb_cis,col="grey",xlab="WT tratadas",ylab="Rb-/- tratadas",pch=19,cex=0.5)

plot(WT_nocis, WT_cis, col="grey",xlab="WT no tratadas",ylab="WT tratadas",pch=19,cex=0.5)

plot(Rb_nocis, Rb_cis, col="grey",xlab="Rb-/- no tratadas",ylab="Rb-/- tratadas",pch=19,cex=0.5)



# Differential gene expression analysis

experimental.design <- model.matrix(~ -1+factor(c(1,1,1,2,2,2,3,3,3,4,4,4)))
colnames(experimental.design) <- c("Rb_nocis","WT_nocis","Rb_cis","WT_cis")
linear.fit <- lmFit(expression.level, experimental.design)
##lmFit fits to a linear model

contrast.matrix <- makeContrasts(Rb_nocis-WT_nocis, Rb_cis-WT_cis,
                                 WT_nocis-WT_cis, Rb_nocis-Rb_cis,
                                 levels=c("Rb_nocis","WT_nocis","Rb_cis","WT_cis"))
head(contrast.matrix)
contrast.linear.fit <- contrasts.fit(linear.fit, contrast.matrix)
## contrasts.fit: Given a linear model fit to the microarray data,
## it computes estimated coefficients and standard errors for a given set of comparisons. 
contrast.results <- eBayes(contrast.linear.fit)


# Differential gene expression tables

Rb_nocis.WT_nocis <- topTable(contrast.results, number=22690,coef=1,sort.by="logFC")
head(Rb_nocis.WT_nocis)

Rb_cis.WT_cis <- topTable(contrast.results, number=22690,coef=2,sort.by="logFC")
head(Rb_cis.WT_cis)

WT_nocis.WT_cis <- topTable(contrast.results, number=22690,coef=3,sort.by="logFC")
head(WT_nocis.WT_cis)

Rb_nocis.Rb_cis <- topTable(contrast.results, number=22690,coef=4,sort.by="logFC")
head(Rb_nocis.Rb_cis)

## Comparison 1

fold.change.Rb_nocis.WT_nocis <- Rb_nocis.WT_nocis[["logFC"]]
genes.ids.Rb_nocis.WT_nocis <- rownames(Rb_nocis.WT_nocis)

activated.genes.Rb_nocis.WT_nocis <- genes.ids.Rb_nocis.WT_nocis[fold.change.Rb_nocis.WT_nocis > log2(2)]
repressed.genes.Rb_nocis.WT_nocis <- genes.ids.Rb_nocis.WT_nocis[fold.change.Rb_nocis.WT_nocis < - log2(2)]

length(activated.genes.Rb_nocis.WT_nocis)
length(repressed.genes.Rb_nocis.WT_nocis)

plot(WT_nocis,Rb_nocis,pch=19,cex=0.5,col="grey",xlab="WT no tratadas",ylab="Rb-/- no tratadas",cex.lab=1.5)
points(WT_nocis[activated.genes.Rb_nocis.WT_nocis],Rb_nocis[activated.genes.Rb_nocis.WT_nocis],pch=19,cex=0.5,col="green")
points(WT_nocis[repressed.genes.Rb_nocis.WT_nocis],Rb_nocis[repressed.genes.Rb_nocis.WT_nocis],pch=19,cex=0.5,col="red")
text(WT_nocis["1417019_a_at"]+1,Rb_nocis["1417019_a_at"]+1,"Cdc6", col="black", cex=1)
text(WT_nocis["22535_at"]-0.5,Rb_nocis["1422535_at"]-0.5,"Ciclina E2", col="black", cex=1)
text(WT_nocis["1450223_at"]+0.5,Rb_nocis["1450223_at"]+0.5,"Apaf1", col="black", cex=1)

## Comparison 2

fold.change.Rb_cis.WT_cis <- Rb_cis.WT_cis[["logFC"]]
genes.ids.Rb_cis.WT_cis <- rownames(Rb_cis.WT_cis)

activated.genes.Rb_cis.WT_cis <- genes.ids.Rb_cis.WT_cis[fold.change.Rb_cis.WT_cis > log2(2)]
repressed.genes.Rb_cis.WT_cis <- genes.ids.Rb_cis.WT_cis[fold.change.Rb_cis.WT_cis < - log2(2)]

length(activated.genes.Rb_cis.WT_cis)
length(repressed.genes.Rb_cis.WT_cis)

plot(WT_cis,Rb_cis,pch=19,cex=0.5,col="grey",xlab="WT tratadas",ylab="Rb-/- tratadas",cex.lab=1.5)
points(WT_cis[activated.genes.Rb_cis.WT_cis],Rb_cis[activated.genes.Rb_cis.WT_cis],pch=19,cex=0.5,col="green")
points(WT_cis[repressed.genes.Rb_cis.WT_cis],Rb_cis[repressed.genes.Rb_cis.WT_cis],pch=19,cex=0.5,col="red")
text(WT_cis["1417019_a_at"]+1,Rb_cis["1417019_a_at"]+1,"Cdc6", col="black", cex=1)
text(WT_cis["1422535_at"]-0.5,Rb_cis["1422535_at"]-0.5,"Ciclina E2", col="black", cex=1)
text(WT_cis["1450223_at"]+0.5,Rb_cis["1450223_at"]+0.5,"Apaf1", col="black", cex=1)


## Annotate and Export data

activated.genes.Rb_nocis.WT_nocis.table <- aafTableAnn(activated.genes.Rb_nocis.WT_nocis, "moe430a.db", aaf.handler())
saveHTML(activated.genes.Rb_nocis.WT_nocis.table, file="activated_genes_Rb_nocis-WT_nocis.html")
saveText(activated.genes.Rb_nocis.WT_nocis.table, file="activated_genes_Rb_nocis-WT_nocis.txt")

repressed.genes.Rb_nocis.WT_nocis.table <- aafTableAnn(repressed.genes.Rb_nocis.WT_nocis, "moe430a.db", aaf.handler())
saveHTML(repressed.genes.Rb_nocis.WT_nocis.table, file="repressed_genes_Rb_nocis-WT_nocis.html")
saveText(repressed.genes.Rb_nocis.WT_nocis.table, file="repressed_genes_Rb_nocis-WT_nocis.txt")

activated.genes.Rb_cis.WT_cis.table <- aafTableAnn(activated.genes.Rb_cis.WT_cis, "moe430a.db", aaf.handler())
saveHTML(activated.genes.Rb_cis.WT_cis.table, file="activated_genes_Rb_cis-WT_cis.html")
saveText(activated.genes.Rb_cis.WT_cis.table, file="activated_genes_Rb_cis-WT_cis.txt")

repressed.genes.Rb_cis.WT_cis.table <- aafTableAnn(repressed.genes.Rb_cis.WT_cis, "moe430a.db", aaf.handler())
saveHTML(repressed.genes.Rb_cis.WT_cis.table, file="repressed_genes_Rb_cis-WT_cis.html")
saveText(repressed.genes.Rb_cis.WT_cis.table, file="repressed_genes_Rb_cis-WT_cis.txt")


## Heatmap

complete.DEGs <- c(activated.genes.Rb_cis.WT_cis,repressed.genes.Rb_cis.WT_cis,
                   activated.genes.Rb_nocis.WT_nocis,repressed.genes.Rb_nocis.WT_nocis)

complete.DEGs <- unique(complete.DEGs)
length(complete.DEGs)
DEG.expression <- mean.expression[complete.DEGs,]

normalized.DEG.expression <- t(scale(t(DEG.expression)))

heatmap.2(normalized.DEG.expression,Colv = FALSE,dendrogram="row",labRow=c(""),density.info="none",trace="none",col=redgreen(100),margins = c(8,8),cexCol=1.2)

