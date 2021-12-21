## This script generates a Venn diagram representation of the overlap 
## between sets of  genes. Additonally, this scripts determines the 
## significance of this overlap using Fisher's exact test.

## Loading the library VennDiagram
library(VennDiagram)

#Load data
tf1.file <- "tf1.txt" ## File with the targets of a TF
tf2.file <- "tf2.txt" ## File with the targets of a TF
tf1.name <- "CO" ## Name of the TF1
tf2.name <- "PRR5"

## Reading in activated genes in 35SCO at ZT16
tf1.targets <- as.vector(read.table(file=tf1.file)[[1]])
length(tf1.targets)
tf2.targets <- as.vector(read.table(file=tf2.file)[[1]])
length(tf2.targets)


## Intersection 
tf1.tf2.overlap <- intersect(tf1.targets, tf2.targets)
length(tf1.tf2.overlap)


## Clearing the graphical panel and generating the venn diagram
tiff(file=paste0(tf1.name,"_", tf2.name, ".tiff"),
     width     = 4,
     height    = 4,
     units     = "in",
     res       = 600
)
grid.newpage()
draw.pairwise.venn(area1 = length(tf1.targets),
                   area2 = length(tf2.targets),
                   cross.area = length(tf1.tf2.overlap),
                   # category = c(paste0(tf1.name, " targets"),
                                # paste("\n   ", tf.name, "targets")),
                   fill=c("mediumseagreen","darkred"), 
                   cat.pos = c(0,190), cat.dist = rep(0.05, 1.3), cat.cex=1.6,
                   lwd=3,cex=1.4, fontfamily="arial",
                   cat.fontfamily="arial",fontface="bold",
                   cat.fontface="bold")
dev.off()


## Computing the significance of the overlap using fisher exact test
fisher.test(matrix(c(33600,length(tf1.targets),
                     length(tf2.targets),length(tf1.tf2.overlap)),ncol=2),alternative="greater")


######Three sets###############################################

tf3.targets <- as.vector(read.table(file="HY5_targets_2000bp.txt")[[1]])
tf3.name <- "HY5"

triple.intersection <- intersect(tf1.tf2.overlap, tf3.targets)
length(triple.intersection)

n12 <- length(tf1.tf2.overlap)
n13 <- length(intersect(tf1.targets, tf3.targets))
n23 <- length(intersect(tf2.targets, tf3.targets))
n123 <- length(triple.intersection)

## Clearing the graphical panel and generating the venn diagram
# help(draw.triple.venn)


grid.newpage()

overrideTriple=T

draw.triple.venn(area1 = length(tf1.targets),
                 area2 = length(tf2.targets),
                 area3 = length(tf3.targets),
                 n12 = n12, n13 = n13, n23 = n23, n123 = n123,
                 euler.d= TRUE, scaled = TRUE,
                    
                   # cross.area = length(activated.35sco.tf.targets),
                   # category = c("    35SCO Activated \n Genes at ZT16",
                   #              paste("\n   ", tf.name, "targets"), "HY5 targets"),
                   fill=c("mediumseagreen","darkred", "grey"), 
                   cat.pos = c(0,45,90), cat.dist = c(0.02,0.02,0.02), cat.cex=1.6,
                   lwd=3,cex=1.4, fontfamily="arial",
                   cat.fontfamily="arial",fontface="bold",
                   cat.fontface="bold")
                 # ,label.col = 0) ##transparencia al 0 para quitar los números



################ Tiff image ##################33


tiff(file=paste0(tf1.name,"_", tf2.name, "_",
                 tf3.name,".tiff"),
     width     = 4,
     height    = 4,
     units     = "in",
     res       = 600
)

overrideTriple=T

draw.triple.venn(area1 = length(tf1.targets),
                 area2 = length(tf2.targets),
                 area3 = length(tf3.targets),
                 n12 = n12, n13 = n13, n23 = n23, n123 = n123,
                 euler.d= TRUE, scaled = TRUE,
                 
                 # cross.area = length(activated.35sco.tf.targets),
                 # category = c("    35SCO Activated \n Genes at ZT16",
                 #              paste("\n   ", tf.name, "targets"), "HY5 targets"),
                 fill=c("mediumseagreen","darkred", "grey"), 
                 cat.pos = c(0,45,90), cat.dist = c(0.02,0.02,0.02), cat.cex=1.6,
                 lwd=3,cex=1.4, fontfamily="arial",
                 cat.fontfamily="arial",fontface="bold",
                 cat.fontface="bold",
                 label.col = 0) ##transparencia al 0 para quitar los números

dev.off()
