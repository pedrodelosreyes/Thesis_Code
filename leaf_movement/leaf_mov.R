## This script takes movement data from TRIP software,
## plot it and analyze it with circacompare package

##Author: Pedro de los Reyes
##Date: 2020

# Read the files with the leaf movement data. Output from
# TRIP software

co10_1 <- read.table(file= "crop_co10_1.csv", sep = ",")
co10_2 <- read.table(file= "crop_co10_2.csv", sep = ",")
co10_3 <- read.table(file= "crop_co10_3.csv", sep = ",")
co10_4 <- read.table(file= "crop_co10_4.csv", sep = ",")


head(co10_1)
nrow(co10_1)
tail(co10_1)

time <- c()

time <- seq(from=6, to =nrow(co10_1)*6, by=6)
horas <- time/60

# Adjust the initial recording time
final.horas <- horas + 2.25  

co10_1.data <- data.frame(final.horas, co10_1$V1)
head(co10_1.data)

co10_2.data <- data.frame(final.horas, co10_2$V1)
head(co10_2.data)

co10_3.data <- data.frame(final.horas, co10_3$V1)
head(co10_3.data)

co10_4.data <- data.frame(final.horas, co10_4$V1)
head(co10_4.data)


plot(co10_1.data, type ="l")


##Plot the leaf movement

plot(co10_1.data, type="l", lwd=2, ylim=c(-0.2, 0.2), 
     ylab="Vertical motion", xlab="ZT (hours in continuous light)")
lines(co10_2.data, col="blue", lwd=2)
lines(co10_3.data, col="green", lwd=2)
lines(co10_4.data, col="red", lwd=2)


## Add photoperiod rectangle
polygon(x=c(0:16,16:0),y=c(rep(-0.14,17),rep(-0.2,17)))
polygon(x=c(16:24,24:16),y=c(rep(-0.14,9),rep(-0.2,9)),col = "grey90")

polygon(x=c(24:40,40:24),y=c(rep(-0.14,17),rep(-0.2,17)))
polygon(x=c(40:48,48:40),y=c(rep(-0.14,9),rep(-0.2,9)),col = "grey90")

polygon(x=c(48:64,64:48),y=c(rep(-0.14,17),rep(-0.2,17)))
polygon(x=c(64:72,72:64),y=c(rep(-0.14,9),rep(-0.2,9)),col = "grey90")

polygon(x=c(72:88,88:72),y=c(rep(-0.14,17),rep(-0.2,17)))
polygon(x=c(88:96,96:88),y=c(rep(-0.14,9),rep(-0.2,9)),col = "grey90")

polygon(x=c(96:112,112:96),y=c(rep(-0.14,17),rep(-0.2,17)))
polygon(x=c(112:120,120:112),y=c(rep(-0.14,9),rep(-0.2,9)),col = "grey90")

polygon(x=c(120:136,136:120),y=c(rep(-0.14,17),rep(-0.2,17)))
polygon(x=c(136:144,144:136),y=c(rep(-0.14,9),rep(-0.2,9)),col = "grey90")


# legend("topright",
#        legend=c("Col-0", expression(italic("co10")),expression(italic("prr5")),expression(paste("35S:",italic("CO")))),#, expression(italic("co-10"))),
#        col=c("black","blue","green","red"),
#        pch=c(0,1),
#        lwd=3)



##Smoothing

# Aplying a smoothing spline
# df=12 I sete 12 degrees of freedom
# Ref: https://www.nature.com/articles/s41598-020-65372-8
df <- 12
smooth.co10_1.data <- smooth.spline(co10_1.data, df=df) #df=12 pongo 12 grados de libertad
smooth.co10_2.data <- smooth.spline(co10_2.data, df=df) #df=12 pongo 12 grados de libertad
smooth.co10_3.data <- smooth.spline(co10_3.data, df=df) #df=12 pongo 12 grados de libertad
smooth.co10_4.data <- smooth.spline(co10_4.data, df=df) #df=12 pongo 12 grados de libertad

# smooth24.co10_a.data <- smooth.spline(co10_a.data, df=24) #df=24 pongo 24 grados de libertad


plot(smooth.co10_1.data, type="l", lwd=2, ylim=c(-0.20, 0.20), 
     ylab="Vertical motion", xlab="ZT (hours in continuous light)")
# lines(smooth12.co10_a.data, col="grey", lwd=2)
lines(smooth.co10_2.data, col="blue", lwd=2)
lines(smooth.co10_3.data, col="red", lwd=2)
lines(smooth.co10_4.data, col="green", lwd=2)


## Add photoperiod rectangle
polygon(x=c(0:16,16:0),y=c(rep(-0.095,17),rep(-0.2,17)))
polygon(x=c(16:24,24:16),y=c(rep(-0.095,9),rep(-0.2,9)),col = "grey90")

polygon(x=c(24:40,40:24),y=c(rep(-0.095,17),rep(-0.2,17)))
polygon(x=c(40:48,48:40),y=c(rep(-0.095,9),rep(-0.2,9)),col = "grey90")

polygon(x=c(48:64,64:48),y=c(rep(-0.095,17),rep(-0.2,17)))
polygon(x=c(64:72,72:64),y=c(rep(-0.095,9),rep(-0.2,9)),col = "grey90")

polygon(x=c(72:88,88:72),y=c(rep(-0.095,17),rep(-0.2,17)))
polygon(x=c(88:96,96:88),y=c(rep(-0.095,9),rep(-0.2,9)),col = "grey90")

polygon(x=c(96:112,112:96),y=c(rep(-0.095,17),rep(-0.2,17)))
polygon(x=c(112:120,120:112),y=c(rep(-0.095,9),rep(-0.2,9)),col = "grey90")

polygon(x=c(120:136,136:120),y=c(rep(-0.095,17),rep(-0.2,17)))
polygon(x=c(136:144,144:136),y=c(rep(-0.095,9),rep(-0.2,9)),col = "grey90")

legend("topleft",
       legend=c("Col-0", expression(italic("co10")),expression(italic("prr5")),expression(paste("35S:",italic("CO")))),#, expression(italic("co-10"))),
       col=c("black","blue","green","red"),
       pch=c(0,1),
       lwd=3)





##########-----Rhythmic analysis---#####


# install.packages("devtools")
# devtools::install_github("RWParsons/circacompare")
# install.packages("nlstools")

library(circacompare)
# help(circacompare)


# circa_single() is used to analyse a single rhythm and provide 
# estimates of its mesor, amplitude and phase.

# circacompare() is used to analyse a dataset with two groups of 
# rhythmic data. It fits a model to estimate and statistically 
# support differences in mesor, amplitude and phase between the two groups.

#Examples

# create an example data frame with a phase shift between groups of 6 hours.
df <- make_data(hours_diff = 6)
# call circacompare.
out <- circacompare(x = df, col_time = "time", col_group = "group", col_outcome = "measure")
# view the graph.
out[[1]]
# view the results.
out[[2]]
# view the nls object summary and confidence intervals for parameters.  To understand the parameters in this object, it would be useful to read the circacompare publication (https://doi.org/10.1093/bioinformatics/btz730)
summary(out[[3]])
nlstools::confint2(out[[3]])
head(df)

co10_1.data <- cbind(co10_1.data, rep(x = "co10_1",times=1405))
co10_2.data <- cbind(co10_2.data, rep(x = "co10_2",times=1405))
co10_3.data <- cbind(co10_3.data, rep(x = "co10_3",times=1405))
co10_4.data <- cbind(co10_4.data, rep(x = "co10_4",times=1405))



colnames(co10_1.data) <- c("time", "measure", "plant")
colnames(co10_2.data) <- c("time", "measure", "plant")
colnames(co10_3.data) <- c("time", "measure", "plant")
colnames(co10_4.data) <- c("time", "measure", "plant")


head(co10_1.data)
head(co10_2.data)

nrow(co10_1.data)

only.co10_1 <- circa_single(x=co10_1.data, col_time="time", col_outcome="measure", period=24)
only.co10_1[[3]]
only.co10_1[[1]]
only.co10_1[[2]]
co10_1.phase <- only.co10_1[[2]]$phase_radians*3.8197
# Multiplying by 3.8917 ti convert radians to hours
#Para pasar de radianes a ángulos de hora lo multiplico por 3.8197
#phase refers to the time at which the response variable peaks


only.co10_2 <- circa_single(x=co10_2.data, col_time="time", col_outcome="measure", period=24)
only.co10_2[[3]]
only.co10_2[[1]]
only.co10_2[[2]]
co10_2.phase <- only.co10_2[[2]]$phase_radians*3.819


only.co10_3 <- circa_single(x=co10_3.data, col_time="time", col_outcome="measure", period=24)
only.co10_3[[3]]
only.co10_3[[1]]
only.co10_3[[2]]
co10_3.phase <- only.co10_3[[2]]$phase_radians*3.8197


only.co10_4 <- circa_single(x=co10_4.data, col_time="time", col_outcome="measure", period=24)
only.co10_4[[3]]
only.co10_4[[1]]
only.co10_4[[2]]
co10_4.phase <- only.co10_4[[2]]$phase_radians*3.8197






# co10_.prr5.data <- rbind(co10_.data, prr5.data)
# out <- circacompare(x = co10_.prr5.data, col_time = "time", col_group = "plant", col_outcome = "measure", period = 24)
# out[[1]]
# out[[2]]
# summary(out[[3]])

# Phase difference estimate  -1.086989e-01
# P-value for difference in phase   5.057888e-01




#####Notes

# When to use what
# 
# If you are looking to estimate the rhythmic parameters of a single group, 
# use circa_single(). If you are looking to estimate the differences between
# two rhythmic datasets, use circacompare()
# 
# If your data has a hierarchical structure, a mixed model may be more 
# appropriate. This may be the case if you have repeated measurements 
# from the same subjects/tissues over time, for example. In this case, 
# consider the equivalents of the above: circa_single_mixed() and 
# circacompare_mixed(). In addition to what has been described, 
# these mixed models require the user to specify which parameters 
# ought to have a random effect and the identifying column (col_id) 
# for this hierarchical structure.


# Function to get the change of period across days (i) in continuous light

get.periods <- function(mov.data)
{
  periods <- c()
  for (i in 1:5)
  {
    circa.data <- circa_single(x=mov.data[1:(i*240)+1,], col_time="time", col_outcome="measure", period=24)
    periods[i] <- (circa.data[[2]]$phase_radians*3.8197)*2
  }
  return(periods)
}

co10_.periods <- get.periods(co10_.data)
prr5.periods <- get.periods(prr5.data)
co10.periods <- get.periods(co10.data)
x35sco.periods <- get.periods(x35sco.data)

##Plot

plot(co10_.periods, type="l", lwd=3, ylim=c(25, 30), 
     ylab="Period (hours)", xlab="Days in continuous light")
lines(co10.periods, col="blue", lwd=3)
lines(prr5.periods, col="green", lwd=3)
lines(x35sco.periods, col="red", lwd=3)

legend("bottomleft",
       legend=c("Col-0", expression(italic("co10")),expression(italic("prr5")),expression(paste("35S:",italic("CO")))),#, expression(italic("co-10"))),
       col=c("black","blue","green","red"),
       pch=c(0,1),
       lwd=3)




#Referencia

# CircaCompare v0.1.0 (29) was used to detect a phase shift in 
# the expression pattern of genes cycling in two or more tissues. 
# For genes cycling in three or four tissues, all possible pairs 
# of tissues were considered. A phase shift was considered 
# statistically significant if the difference in the peak 
# times estimated by CircaCompare exceeded 2 hours with an 
# FDR < 0.05. GO enrichment for the phase-shifted genes was 
# performed with the R package topGO v.2.34.0. 





##############################################################
#Plotting

##ggplot for shaded area with raw data
# https://www.r-graph-gallery.com/104-plot-lines-with-error-envelopes-ggplot2.html
all.data <- data.frame(co10_1.data$measure, co10_2.data$measure, co10_3.data$measure, co10_4.data$measure)
head(all.data)

all.data$time <- co10_1.data$time
all.data$mean <- rowMeans(all.data[,1:4])
all.data$max <- apply(all.data[,1:4], 1, max)
all.data$min <- apply(all.data[,1:4], 1, min)
all.data$sd <- apply(all.data[,1:4], 1, sd)

all.data$type <- "co10"
head(all.data)
# write.table(all.data, file = "../all_plants/co10_data.txt", sep = "\t")


night1 <- data.frame (xmin=16, xmax=24, ymin=-Inf, ymax=Inf)
night2 <- data.frame (xmin=40, xmax=48, ymin=-Inf, ymax=Inf)
night3 <- data.frame (xmin=64, xmax=72, ymin=-Inf, ymax=Inf)
night4 <- data.frame (xmin=88, xmax=96, ymin=-Inf, ymax=Inf)
night5 <- data.frame (xmin=112, xmax=120, ymin=-Inf, ymax=Inf)
night6 <- data.frame (xmin=136, xmax=144, ymin=-Inf, ymax=Inf)

# generate break positions
breaks = c(0,16,24,40,48,64,72,88,96,112,120,136,144)
# and labels
labels = as.character(breaks)


# ggplot(data=all.data, aes(x=time, y=mean, ymin=min, ymax=max, fill=type, linetype=type)) + 
ggplot(data=all.data, aes(x=time, y=mean, ymin=mean-sd, ymax=mean+sd, fill=type, linetype=type)) + 
  geom_line() + 
  geom_ribbon(alpha=0.5) + 
  # scale_x_log10() + 
  # scale_y_log10() + 
  xlab("Hours in continuous light") + 
  ylab("Relative vertical motion") +
  geom_rect(data=night1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="black", alpha=0.1, inherit.aes = FALSE) +
  geom_rect(data=night2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="black", alpha=0.1, inherit.aes = FALSE) +
  geom_rect(data=night3, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="black", alpha=0.1, inherit.aes = FALSE) +
  geom_rect(data=night4, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="black", alpha=0.1, inherit.aes = FALSE) +
  geom_rect(data=night5, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="black", alpha=0.1, inherit.aes = FALSE) +
  geom_rect(data=night6, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="black", alpha=0.1, inherit.aes = FALSE) +
  
  scale_x_continuous(limits = c(1, 144), breaks = breaks, labels = labels,
                     name = "Time in constant light (h)") +
  guides(fill=guide_legend(title=NULL)) +
  guides(linetype=guide_legend(title=NULL))




##ggplot for shaded area with smoothed data
all.data <- data.frame(smooth.co10_1.data$y, smooth.co10_2.data$y, smooth.co10_3.data$y,
                       smooth.co10_4.data$y)

all.data$time <- co10_1.data$time
all.data$mean <- rowMeans(all.data[,1:4])
all.data$max <- apply(all.data[,1:4], 1, max)
all.data$min <- apply(all.data[,1:4], 1, min)
all.data$sd <- apply(all.data[,1:4], 1, sd)
all.data$type <- "co10"
head(all.data)
# write.table(all.data, file = "../all_plants/smoothed.co10_data.txt", sep = "\t")


library(ggplot2)
# ggplot(data=all.data, aes(x=time, y=mean, ymin=min, ymax=max, fill=type, linetype=type)) + 
ggplot(data=all.data, aes(x=time, y=mean, ymin=mean-sd, ymax=mean+sd, fill=type, linetype=type)) + 
  geom_line() + 
  geom_ribbon(alpha=0.5) + 
  # scale_x_log10() + 
  # scale_y_log10() + 
  xlab("Hours in continuous light") + 
  ylab("Vertical motion") +
  geom_rect(data=night1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="black", alpha=0.1, inherit.aes = FALSE) +
  geom_rect(data=night2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="black", alpha=0.1, inherit.aes = FALSE) +
  geom_rect(data=night3, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="black", alpha=0.1, inherit.aes = FALSE) +
  geom_rect(data=night4, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="black", alpha=0.1, inherit.aes = FALSE) +
  geom_rect(data=night5, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="black", alpha=0.1, inherit.aes = FALSE) +
  geom_rect(data=night6, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="black", alpha=0.1, inherit.aes = FALSE) 







#Calculating the period as the distance between peaks
#for each day. Then calculate the mean of time distances. 
#IN SMOOTHED DATA!!! using the mean measures
peak1 <- subset(all.data, mean == (max(all.data[1:241,"mean"])))$time
peak2 <- subset(all.data, mean == (max(all.data[241:481,"mean"])))$time
peak3 <- subset(all.data, mean == (max(all.data[481:721,"mean"])))$time
peak4 <-subset(all.data, mean == (max(all.data[721:961,"mean"])))$time
peak5 <-subset(all.data, mean == (max(all.data[961:1201,"mean"])))$time
peak6 <- subset(all.data, mean == (max(all.data[1201:1405,"mean"])))$time

co10.period <- mean(c(peak2-peak1, peak3-peak2, peak4-peak3, peak5-peak4, peak6-peak5))
#25.94 horas

head(all.data)
