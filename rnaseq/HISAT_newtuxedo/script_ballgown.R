####################################################################
## Análisis de los resultados de HISAT2 y STRINGTIE               ##
## usando los paquetes de R ballgown y limma.                     ##  
##                                                                ##
## Autores: Francisco J. Romero-Campero fran@us.es                ##
##          Pedro de los Reyes Rodriguez pedro.reyes@ibvf.csic.es ##
####################################################################

## El paquete de bioconductor ballgown proporciona las funciones necesarias para 
## realizar un análisis de expresión génica diferencial y visualizar los resultados
## a partir de procesamiento de los datos brutos de secuenciación realizados con 
## hisat2 y stringtie. 

## Para ejecutar con éxito este script es necesario descargar la carpeta samples
## completa a tu ordenador, mover este script a la carpeta samples y fijar el 
## Working Directory To Source File Location. 

## Instalación y carga de los paquetes necesarios. Sólo es necesario instalar los
## paquetes la primera vez que se ejecuta este script en un ordenador el resto de las
## veces bastará cargar los paquetes simplemente. 

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("ballgown")


library(ballgown)
library(genefilter)

## Para cargar los datos es necesario crear previamente un fichero tabular
## que contenga como primera columna los nombres de las carpetas donde se guarda
## cada muesra típicamente sample1, sample2, etc ... El resto de columnas
## haran referencia al genotipo, tratamiento y demás caracteriśticas de cada muestra. 
pheno.data <- read.csv("pheno_data.csv")
pheno.data

## La función ballgown se usa para cargar o leer los datos. Es necesario especificar
## el directorio donde se encuentra las muestras. En nuestro caso especificamos .
## para indicar que se encuentran en el actual directorio. 
bg.data <- ballgown(dataDir = ".", samplePattern = "sample", pData=pheno.data)
bg.data
sampleNames(bg.data)

## where did you merge from
bg.data@dirs
## when did you merge?
bg.data@mergedDate
## what did you import?
bg.data@meas

##Look at the exon, intron and transcript  data
structure(bg.data)$exon
structure(bg.data)$intron
structure(bg.data)$trans

## La función gexpr extrae los niveles de expresión génicos en cada muestra
## medidos como FPKM (Fragments Per Kilobase of exon and Million of mapped reads)
gene.expression <- gexpr(bg.data)
head(gene.expression)
dim(gene.expression)

gene.expression[1:10,]
gene.expression[50:70,]

gene.expression["AT5G15840",] #co
gene.expression["AT5G24470",] #prr5
gene.expression["AT1G65480",] #ft
gene.expression["AT1G01060",] #lhy




## La función texpr extrae los niveles de expresión génicos en cada muestra
## medidos como FPKM para transcritos
transcript.expression <- texpr(bg.data)
head(transcript.expression)
dim(transcript.expression)

rownames(transcript.expression) <- transcriptNames(bg.data)
# transcript.expression["AT2G28670.1",]
# transcript.expression["AT2G28670.2",]



## Nombramos las columnas con los nombres de nuestras muestras. 
colnames(transcript.expression) <- c("35SCO_1", "35SCO_2",
                                     "Col0_1","Col0_2")
colnames(gene.expression) <- c("35SCO_1", "35SCO_2",
                               "Col0_1","Col0_2")

##Note: to avoid log 0, add 1 to log values

## Previsualizamos la similitud entre las réplicas
# plot(log2(gene.expression[,1]+1),log2(gene.expression[,2]+1),pch=19,cex=0.7,xlab="ablated_1",ylab=substitute(italic("ablated_2")),cex.lab=1.25)
# plot(log2(gene.expression[,1]+1),log2(gene.expression[,2]+1),pch=19,cex=0.7,xlab="log2(FPKM + 1) wta",ylab="log2(FPKM + 1) wtb",cex.lab=1.25)
# plot(log2(gene.expression[,3]+1),log2(gene.expression[,4]+1),pch=19,cex=0.7,xlab="log2(FPKM + 1) m2a",ylab="log2(FPKM + 1) m2b",cex.lab=1.25)
# plot(log2(gene.expression[,5]+1),log2(gene.expression[,6]+1),pch=19,cex=0.7,xlab="log2(FPKM + 1) b91a",ylab="log2(FPKM + 1) b91b",cex.lab=1.25)
# plot(log2(gene.expression[,7]+1),log2(gene.expression[,8]+1),pch=19,cex=0.7,xlab="log2(FPKM + 1) 902a",ylab="log2(FPKM + 1) 902b",cex.lab=1.25)


## Construimos un boxplot para comprobar que las distribuciones globales de las
## muestras son similares y comparables.
boxplot(log2(gene.expression+1),col=rainbow(ncol(gene.expression)),ylab="log2(FPKM + 1)",cex.lab=1.5)


# Normalización Cuantílica Superior

upper.quantiles <- vector(mode="numeric",length=ncol(gene.expression))
normalized.expression <- matrix(ncol = ncol(gene.expression), nrow = nrow(gene.expression))

for(i in 1:ncol(gene.expression))
{
  upper.quantiles[i] <- quantile(gene.expression[,i],probs=0.75)
}

## División de cada columna por su cuartil superior. 
for(i in 1:ncol(gene.expression))
{
  normalized.expression[,i] <- gene.expression[,i] / upper.quantiles[i]
}
head(normalized.expression)
dim(normalized.expression)
colnames(normalized.expression) <- colnames(gene.expression)
rownames(normalized.expression) <- rownames(gene.expression)
head(normalized.expression)

normalized.expression["AT5G15840",] #co
normalized.expression["AT5G24470",] #prr5
normalized.expression["AT1G65480",] #ft
normalized.expression["AT1G01060",] #lhy


## Transformación logarítmica y visualización de los datos normalizados.
log.normalized.expression <- log2(normalized.expression+1)
head(log.normalized.expression)
dim(log.normalized.expression)

log.normalized.expression["AT5G15840",] #co
log.normalized.expression["AT5G24470",] #prr5
log.normalized.expression["AT1G65480",] #ft
log.normalized.expression["AT1G01060",] #lhy

# Previsualización de los Datos con Boxplot
boxplot(log.normalized.expression,col=rainbow(ncol(normalized.expression)),ylab="log2(FPKM + 1)",
        cex.lab=1, outline=F)

# write.table(file="tablas/logaritmic_normalized_expression.txt", x = log.normalized.expression,
#             col.names = TRUE, row.names = TRUE, sep = "\t")


# Mean expression matrix
x35SCO.fpkm <- (log.normalized.expression[,1] + log.normalized.expression[,2])/2
col0.fpkm <- (log.normalized.expression[,3] + log.normalized.expression[,4])/2

mean.expression <- matrix(c(x35SCO.fpkm, col0.fpkm), ncol=2)
colnames(mean.expression) <- c("35SCO", "Col0")
rownames(mean.expression) <- names(x35SCO.fpkm)
head(mean.expression)
mean.expression["AT5G15840",]
mean.expression["AT5G24470",]

# write.table(mean.expression, file="tablas/log_mean_expression.txt",
# col.names = TRUE, row.names = TRUE, sep = "\t")


        

#--------------------------------------------------------------------------------------#
##### ------Limma (expresión diferencial) ----- #####
##El paquete **limma** (Linear Models for Microarray Analysis) proporciona las 
##funciones necesarias para determinar los genes expresados de forma 
##diferencial (DEGs). 

library(limma)

## Especificamos el diseño experimental

experimental.design <- model.matrix(~ -1+factor(c(1,1,2,2)))
colnames(experimental.design) <- c("x35SCO", "Col0")

##A continuación, ajustamos la estimación de los niveles de expresión de cada
##gen a un modelo lineal teniendo en cuenta el diseño experimental. Ene este paso
##fundamentalmente se calcula la media de las réplicas en cada condición.

linear.fit <- lmFit(log.normalized.expression, experimental.design)
linear.fit
head(linear.fit$coefficients)
linear.fit$coefficients["AT5G15840",]


##Para especificar los constrastes a realizar utilizamos la función
##*makeContrasts* que recibe como entrada los contrastes a realizar separados 
##por comas y especificados con los nombres de las dos condiciones 
##correspondientes separadas por un guión -. También recibe el argumento 
##levels, un vector con el nombre de las condiciones:

contrast.matrix <- makeContrasts(x35SCO-Col0, levels=c("x35SCO", "Col0"))
head(contrast.matrix)

##Calculamos el fold-change y los p-valores correspondientes para cada gen en
##cada uno de los constrastes especificados utilizando las funciones *constrasts.fit* 
##y *eBayes*.

contrast.linear.fit <- contrasts.fit(linear.fit, contrast.matrix)
contrast.results <- eBayes(contrast.linear.fit)

nrow(log.normalized.expression)

##NOTA: lmFit assumes your values are already log-transformed, 
## so the log-fold change is just: (rep1 + rep2)/2 - (rep1+rep2)/2

#Contraste 35SCO vs wt
x35SCO.wt <- topTable(contrast.results, number=33602,coef=1,sort.by="logFC")
# x35SCO.wt <- topTable(contrast.results, number=33602,coef=1,sort.by="logFC", p.value=0.05)#, lfc = 1)
head(x35SCO.wt)
x35SCO.wt["AT1G65480",]
x35SCO.wt["AT5G24470",]
nrow(x35SCO.wt)

fold.change.35SCO.wt <- x35SCO.wt$logFC
genes.ids.35SCO.wt <- rownames(x35SCO.wt)
p.value <- x35SCO.wt$adj.P.Val

activated.genes.35SCO <- genes.ids.35SCO.wt[fold.change.35SCO.wt > log2(2)]# & p.value < 0.05] 
repressed.genes.35SCO <- genes.ids.35SCO.wt[fold.change.35SCO.wt < -log2(1.5)]# & p.value < 0.05]

length(activated.genes.35SCO)
length(repressed.genes.35SCO)

# write.table(activated.genes.35SCO, file="tablas/fc_1.5/activated_genes_35SCO_wt_EH.txt", sep="\n", quote = FALSE, col.names = FALSE, row.names = FALSE)
# write.table(repressed.genes.35SCO, file="tablas/fc_1.5/repressed_genes_35SCO_wt.txt", sep="\n", quote = FALSE, col.names = FALSE, row.names = FALSE)
# write.table(activated.genes.35SCO, file="tablas/fc_2/activated_genes_35SCO_wt.txt", sep="\n", quote = FALSE, col.names = FALSE, row.names = FALSE)
# write.table(repressed.genes.35SCO, file="tablas/fc_2/repressed_genes_35SCO_wt.txt", sep="\n", quote = FALSE, col.names = FALSE, row.names = FALSE)

## Previsualizamos el efecto de la mutación en un scatterplot.

plot(x35SCO.fpkm, col0.fpkm, pch=19,cex=0.7,xlab="35SCO",ylab=substitute(italic("wt")),cex.lab=1.25)
points(x35SCO.fpkm[activated.genes.35SCO], col0.fpkm[activated.genes.35SCO], pch=19,cex=0.5,col="red")
points(x35SCO.fpkm[repressed.genes.35SCO], col0.fpkm[repressed.genes.35SCO], pch=19,cex=0.5,col="blue")


#Contraste SUC2 vs wt
suc2.wt <- topTable(contrast.results, number=33602,coef=2,sort.by="logFC")
# suc2.wt <- topTable(contrast.results, number=33602,coef=2,sort.by="logFC", p.value = 0.05)

head(suc2.wt)
suc2.wt["AT5G24470",]
suc2.wt["AT5G60100",]

fold.change.suc2.wt <- suc2.wt$logFC
genes.ids.suc2.wt <- rownames(suc2.wt)
p.value <- suc2.wt$adj.P.Val

activated.genes.suc2 <- genes.ids.suc2.wt[fold.change.suc2.wt > log2(1.5)] # & p.value < 0.05]
repressed.genes.suc2 <- genes.ids.suc2.wt[fold.change.suc2.wt < -log2(1.5)] #& p.value < 0.05]

length(activated.genes.suc2)
length(repressed.genes.suc2)

# write.table(activated.genes.suc2, file="tablas/activated_genes_suc2_wt_fc2.txt", sep="\n", quote = FALSE, col.names = FALSE, row.names = FALSE)
# write.table(repressed.genes.suc2, file="tablas/repressed_genes_suc2_wt_fc2.txt", sep="\n", quote = FALSE, col.names = FALSE, row.names = FALSE)

write.table(activated.genes.suc2, file="tablas/fc_1.5/activated_genes_suc2_wt_fc1.5.txt", sep="\n", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(repressed.genes.suc2, file="tablas/fc_1.5/repressed_genes_suc2_wt_fc1.5.txt", sep="\n", quote = FALSE, col.names = FALSE, row.names = FALSE)

plot(suc2.fpkm, col0.fpkm ,pch=19,cex=0.7,xlab="SUC2",ylab=substitute(italic("wt")),cex.lab=1.25)
points(suc2.fpkm[activated.genes.suc2], col0.fpkm[activated.genes.suc2], pch=19,cex=0.5,col="red")
points(suc2.fpkm[repressed.genes.suc2], col0.fpkm[repressed.genes.suc2], pch=19,cex=0.5,col="blue")


#Contraste co10 vs wt
co10.wt <- topTable(contrast.results, number=33602,coef=3,sort.by="logFC")
# co10.wt <- topTable(contrast.results, number=33602,coef=3,sort.by="logFC", p.value = 0.05)
head(co10.wt)
co10.wt["AT5G15840",]

fold.change.co10.wt <- co10.wt$logFC
genes.ids.co10.wt <- rownames(co10.wt)

# activated.genes.co10 <- genes.ids.co10.wt[fold.change.co10.wt > log2(2)]
# repressed.genes.co10 <- genes.ids.co10.wt[fold.change.co10.wt < -log2(2)]
activated.genes.co10 <- genes.ids.co10.wt[fold.change.co10.wt > log2(1.5)]
repressed.genes.co10 <- genes.ids.co10.wt[fold.change.co10.wt < -log2(1.5)]


length(activated.genes.co10)
length(repressed.genes.co10)

# write.table(activated.genes, file="tablas/activated_genes_co10_wt.txt", sep="\n", quote = FALSE, col.names = FALSE, row.names = FALSE)
# write.table(repressed.genes, file="tablas/repressed_genes_co10_wt.txt", sep="\n", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(activated.genes.co10, file="tablas/fc_1.5/activated_genes_co10_wt.txt", sep="\n", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(repressed.genes.co10, file="tablas/fc_1.5/repressed_genes_co10_wt.txt", sep="\n", quote = FALSE, col.names = FALSE, row.names = FALSE)

plot(co10.fpkm, col0.fpkm ,pch=19,cex=0.7,xlab="co10",ylab=substitute(italic("wt")),cex.lab=1.25)
points(co10.fpkm[activated.genes.co10], col0.fpkm[activated.genes.co10], pch=19,cex=0.5,col="red")
points(co10.fpkm[repressed.genes.co10], col0.fpkm[repressed.genes.co10], pch=19,cex=0.5,col="blue")


#--------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------#

#######-------------HEATMAP#############################################
#draw heatmap for DEGS
head(log.normalized.expression)

degs <- c(activated.genes.35SCO , repressed.genes.35SCO, 
          activated.genes.suc2, repressed.genes.suc2,
          activated.genes.co10, repressed.genes.co10)

length(degs)
degs <- unique(degs)
length(degs)

degs.table <- log.normalized.expression[degs,]
colnames(degs.table) <- colnames(gene.expression)
rownames(pheno.data) <- colnames(gene.expression)
pheno.data$sample <- NULL
head(pheno.data)
head(degs.table)

library(pheatmap)
library(dplyr)
library(readr)
library(viridis) #Default Color Maps
library(gplots)

tiff("heatmap_degs_FC_pvalue_CLUSTER.tiff", height = 4, width = 6, 
     units = 'in', res=300, compression="lzw")
pheatmap(as.matrix(degs.table), cluster_rows = T, cluster_cols = T,
         scale = "row", clustering_distance_rows = "correlation", 
         clustering_method = "complete", #annotation_col = pheno.data, 
         main="DEGs",fontsize_col=14, fontsize_row = 6, color = greenred(28))
dev.off()


#Heatmap for selected genes

# selected.genes <- read.table(file="selected_genes/list.txt", sep = "\t", header = T)
# head(selected.genes)
# 
# up.selected <- selected.genes$up
# length(up.selected)
# up.selected <- up.selected[1:22]
# up.table <- filtered.data[as.vector(up.selected),]
# tiff("selected_genes/genes_up_cluster.tiff", height = 4, width = 6, 
#      units = 'in', res=300, compression="lzw")
# pheatmap(as.matrix(up.table), cluster_rows = F, cluster_cols = T,
#          scale = "row", clustering_distance_rows = "correlation", 
#          clustering_method = "complete", #annotation_col = pheno.data, 
#          main="DEGs",fontsize_col=14, fontsize_row = 6, color = greenred(28))
# dev.off()
# 
# down.selected <- selected.genes$down
# length(down.selected)
# down.table <- filtered.data[as.vector(down.selected),]
# tiff("selected_genes/genes_down_cluster.tiff", height = 4, width = 6, 
#      units = 'in', res=300, compression="lzw")
# pheatmap(as.matrix(down.table), cluster_rows = F, cluster_cols = T,
#          scale = "row", clustering_distance_rows = "correlation", 
#          clustering_method = "complete", #annotation_col = pheno.data, 
#          main="DEGs",fontsize_col=14, fontsize_row = 6, color = greenred(28))
# dev.off()





#--------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------#
##########-----------PCA#############################################

## Draw PCA plot
# transpose the data and compute principal components
# pca.data <- prcomp(t(degs.table)) #Only for DEGs

# Log transform the data before performing PCA
#the log-transform is used to reduce the influence 
#of extreme values or outliers. 


pca.data <- prcomp(t(log.normalized.expression))
nrow(pca.data$rotation)


# Calculate PCA component percentages
pca.data.perc <- round(100*pca.data$sdev^2/sum(pca.data$sdev^2),1)

# Extract 1 and 2 principle components and create a data frame with sample names, first and second principal components and group information
df.pca.data <- data.frame(PC1 = pca.data$x[,1], PC2 = pca.data$x[,2], sample = colnames(gene.expression), 
                          condition = c(rep("35SCO",3), rep("Col0",3), rep("SUC2",3), rep("co10",3)))
head(df.pca.data)

## We will use ggplot2 and ggrepel packages to draw the PCA plots.
# color by sample
library(ggplot2)
library(ggrepel)
ggplot(df.pca.data, aes(PC1,PC2, color = sample))+
  geom_point(size=8)+
  labs(x=paste0("PC1 (",pca.data.perc[1],")"), y=paste0("PC2 (",pca.data.perc[2],")"))


# color by condition/group
# tiff("images/PCA.tiff", height = 8, width = 14, 
     # units = 'in', res=600, compression="lzw")
ggplot(df.pca.data, aes(PC1,PC2, color = condition))+
  geom_point(size=8)+
  labs(x=paste0("PC1 (",pca.data.perc[1],")"), y=paste0("PC2 (",pca.data.perc[2],")"))+
  geom_text_repel(aes(label=sample),point.padding = 0.75)
# dev.off()




#--------------------------------------------------------------------------------------#


###################TRANSCRITOS##################

## En este punto ballgown no realiza ninguna normalización de los datos. 
## Utilizamos el paquete de R Normalyzer para esta tarea. Para ello es neceario generar
## un fichero con un formato específico.

normalyzer.t.data <- data.frame(rownames(transcript.expression), transcript.expression)
dim(normalyzer.t.data)
head(normalyzer.t.data)

colnames(normalyzer.t.data) <- NULL
rownames(normalyzer.t.data) <- NULL

normalyzer.table <- rbind(c(0,rep(1:4,each=2)),                                  
                          rbind(c("Gene",colnames(transcript.expression)),normalyzer.t.data))

head(normalyzer.table)

write.table(normalyzer.table,file="normalyzer_transcripts_table.tab",col.names = F, quote = F, row.names = F, sep="\t")

# biocLite(c("PerformanceAnalytics","ape","raster"))
# install.packages("Normalyzer_1.1.1.tar.gz",repos=NULL,type="source")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("NormalyzerDE")

# I have to learn how to use the updated normalyzerDE package, but, in the meantime
# I am going to use the older Normalyzer

library(Normalyzer)
library(grid)
normalyzer(datafile = "normalyzer_transcripts_table.tab", getjob = "data_normalization_transcripts")

normalized.data <- read.table(file="data_normalization_transcripts//Quantile-normalized.txt", header=T)
head(normalized.data)
normalized.data[is.na(normalized.data)] <- 0

normalized.data[normalized.data$Gene=="AT1G01115.1",]
normalized.data[normalized.data$Gene=="AT1G01600.1",]
normalized.data[normalized.data$Gene=="AT2G28670.1",]
normalized.data[normalized.data$Gene=="AT2G28671.1",]


##La normalización hace log2 de los FPKM, por eso, si quiero obtener los FPKMs normalizados
## elevo a 2 los datos. Pero hay un problema, después de la normalización, los NAs los
## convertimos a ceros, y 2 elevado a 0 es 1. POr tanto, en la tabla vamos a tener genes
## con valor de FPKM=1 que en realidad es 0. Hay que reconvertirlos a 0
normalized.fpkm.data <- 2^normalized.data 
normalized.fpkm.data[normalized.fpkm.data == 1] <- 0
normalized.fpkm.data$Gene <- normalized.data$Gene
head(normalized.fpkm.data)
normalized.fpkm.data[normalized.fpkm.data$Gene=="AT1G01115.1",]
normalized.fpkm.data[normalized.fpkm.data$Gene=="AT5G52300.1",]
write.table(normalized.fpkm.data, file="expression_no_filtered_by_replicates_TRANSCRIPTS.txt", col.names = TRUE, row.names = TRUE, sep = "\t")


## Testeamos la normalización
boxplot(normalized.data[,2:9],col=rep(c("red","blue", "green", "yellow"),each=2))

## Filtramos los datos normalizados (quitamos los genes que tengan menos de
## 1 FPKM en todas las muestras)

## NOTA:   Normalyzer ha hecho previamente el log2 de la expresión


head(normalized.data)
nrow(normalized.data)
genes.ids <- normalized.data$Gene
length(genes.ids)



# help(subset)
genes.to.remove <- subset(normalized.fpkm.data, wta < 1 & wtb < 1 & m2a < 1 & m2b < 1
                          & b91a < 1 & b91b < 1 & x902a < 1 & x902b < 1)
head(genes.to.remove)

genes.ids.to.remove <- genes.to.remove$Gene
length(genes.ids.to.remove) #17646

genes.to.keep <- genes.ids %in% genes.ids.to.remove


filtered.data <- normalized.fpkm.data[!genes.to.keep,]
head(filtered.data)
nrow(filtered.data) #24025
# write.table(filtered.data, file="expression_filtered_by_replicates_TRANSCRIPTS.txt", col.names = TRUE, row.names = TRUE, sep = "\t")


## Calculamos la matrix de expresión media con los datos normalizados. (y en log2) 
wt <- (normalized.data[,2] + normalized.data[,3])/2
m2 <- (normalized.data[,4] + normalized.data[,5])/2
b91 <- (normalized.data[,6] + normalized.data[,7])/2
x902 <- (normalized.data[,8] + normalized.data[,9])/2

mean.expression <- matrix(c(wt, m2, b91, x902),ncol=4)
colnames(mean.expression) <- c("wt","m2", "b91", "902")
rownames(mean.expression) <- normalized.data$Gene
head(mean.expression)
mean.expression["AT1G01600.1",]

write.table(mean.expression, file="mean_expression.txt", col.names = TRUE, row.names = TRUE, sep = "\t")

## Calculamos la matrix de expresión media con los datos filtrados y en FPKM 
wt <- (filtered.data[,2] + filtered.data[,3])/2
m2 <- (filtered.data[,4] + filtered.data[,5])/2
b91 <- (filtered.data[,6] + filtered.data[,7])/2
x902 <- (filtered.data[,8] + filtered.data[,9])/2

mean.expression <- matrix(c(wt, m2, b91, x902),ncol=4)
colnames(mean.expression) <- c("wt","m2", "b91", "902")
rownames(mean.expression) <- filtered.data$Gene
head(mean.expression)
mean.expression["AT1G01600.1",]

write.table(mean.expression, file="mean_expression_filtered_FPKM_TRANSCRIPTS.txt", col.names = TRUE, row.names = TRUE, sep = "\t")



###Barplot individual genes#####
library(ggplot2)
library(wesanderson)
# CO
co.expression <- data.frame(expression=c(log.normalized.expression["AT5G15840",]),
                            genotype= c("35SCO", "35SCO", "Col-0", "Col-0"))
sd.35sco <- sd(subset(co.expression, genotype == "35SCO")$expression)
sd.col0 <- sd(subset(co.expression, genotype == "Col-0")$expression)
mean.35sco <- mean(subset(co.expression, genotype == "35SCO")$expression)
mean.col0 <- mean(subset(co.expression, genotype == "Col-0")$expression)

co.mean <- data.frame(mean=c(mean.35sco, mean.col0),
                      sd=c(sd.35sco, sd.col0),
                      genotype=c("35SCO", "Col-0"))
head(co.mean)

tiff("barplots/co.tiff", height = 6, width = 4, 
     units = 'in', res=300, compression="lzw")
ggplot(data = co.mean, aes(x=genotype, y=mean, fill=genotype)) + 
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.3) + 
  # scale_fill_manual(values=c("#c72a00","#39c700"))
  # scale_fill_brewer(palette="Dark2")
  scale_fill_manual(values=wes_palette(n=3, name="Moonrise2")) +
  ggtitle("CONSTANS") +
  xlab("") + 
  ylab("Log2(FPKM+1)")
dev.off()

# PRR5
prr5.expression <- data.frame(expression=c(log.normalized.expression["AT5G24470",]),
                            genotype= c("35SCO", "35SCO", "Col-0", "Col-0"))
sd.35sco <- sd(subset(prr5.expression, genotype == "35SCO")$expression)
sd.col0 <- sd(subset(prr5.expression, genotype == "Col-0")$expression)
mean.35sco <- mean(subset(prr5.expression, genotype == "35SCO")$expression)
mean.col0 <- mean(subset(prr5.expression, genotype == "Col-0")$expression)

prr5.mean <- data.frame(mean=c(mean.35sco, mean.col0),
                      sd=c(sd.35sco, sd.col0),
                      genotype=c("35SCO", "Col-0"))
head(prr5.mean)

tiff("barplots/prr5.tiff", height = 6, width = 4, 
     units = 'in', res=300, compression="lzw")
ggplot(data = prr5.mean, aes(x=genotype, y=mean, fill=genotype)) + 
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.3) + 
  scale_fill_manual(values=wes_palette(n=3, name="Moonrise2")) +
  ggtitle("PRR5") +
  xlab("") + 
  ylab("Log2(FPKM+1)")
dev.off()

#FT
