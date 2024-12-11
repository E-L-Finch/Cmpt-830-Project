library(dplyr)
library(data.table)
library(ggplot2)
library(lattice)
library(caret)
library (tidyverse)
library(splines)



datafile <- read.csv("SAMPLE LCMS DATAFile.csv")
#this is the trial sample set also available on GitHub

#pulling out the metabolite values from each column from .csv file

Gly <- (select(datafile,GLY))
MH <- (select(datafile,X1MH))
ARG <- (select(datafile,ARG))
SN <- (select(datafile,Sample.Name))
class <- select (datafile, class)
Tm <- select (datafile, TM)

#creating the hydration status normalized values for each metabolite
GlyTm <- Gly / Tm
MhistTm <- MH / Tm
ArgTm <- ARG / Tm

#create a file of the hydration status normalized values to keep for reference. 
TmNormData <- data.frame (SN, GlyTm, MhistTm, ArgTm, Tm, class, 
                          stringsAsFactors = FALSE)
sink('Tm Normalized Data.csv')
write.csv (TmNormData)
sink()
#this step will create a file in your working directory as a .csv with the TM
#normalized data. This data will be used for further processing. 
# this intermediate data is used for lab use


# Sample Processing for MS Run
Raw <- TmNormData
raw <- as.data.table(Raw)
#both a table and data frame format are required for different queries 


#Test of the coding for each metabolite first to see if the code is working
#staring with Glycine metabolites (Gly)
QCname <- raw[class ==2]$Sample.Name
#this is selecting out the sample name (this name must be a numerical value)
#this will select all the sample names which have a corresponding class value=2
GlyQC <- raw[class ==2]$GLY
#this is selecting out the gly values for the QC samples 
#this is done by selecting the rows with a class = 2
#this can be used to create the LOESS values using only the QCs. 



# create an LOESS curve with Optimized Span 

dfmetab <- data.frame (x=c(QCname),
                       y=c(GlyQC))

#we can now apply the data frame to LOESS equation 
#note the number of your data points will impact the choice of span parameter 
#that can be used. This example only looks at 0.75 and 0.9 span options
#this can and should be optimized based on each data set and can be done with 
#k-fold optimization not shown here

metabloess75 <- loess(y ~ x, data=dfmetab, span=.75)
smooth75 <- predict(metabloess75) 

metabloess90 <- loess(y ~ x, data=dfmetab, span=.9)
smooth90 <- predict(metabloess90) 



# Plot the LOESS Curve(s) 

plot(dfmetab$x, dfmetab$y, xlab = "Sample Name", ylab = "Concentration (ng/ml)", pch=19, main='Loess Regression Models for Glycine QCs')
points(smooth75, x=dfmetab$x, col='purple')
lines(smooth75, x=dfmetab$x, col='purple')
points(smooth90, x=dfmetab$x, col='darkgreen')
lines(smooth90, x=dfmetab$x, col='darkgreen')
legend('bottomright', legend=c('.75', '.9'),
       col=c('purple', 'darkgreen'), pch=19, title='Smoothing Span')
#here the data indicates a span value of .9 will give adequate smoothing

# Pull out LOESS Curve Values for Spline Interpolation 
# predict(metabloess90)
# these values ca be used as intermediate data as a check stop



# Create a cubic spline interpolation between LOESS values for every datapt.

qcSp <- spline(QCname, predict(metabloess90), xmin = 1, xmax = 49, n = 49)
plot(qcSp, main = "Cubic Spline Interpolation of Glycine QCs", 
     xlab = "Sample Name", ylab = "Concentration (ng/ml)")
points(smooth90, x=dfmetab$x, col='blue', pch=19)
lines(smooth90, x=dfmetab$x, col='blue')

resid <- c(qcSp$y[1], qcSp$y) - c(qcSp$y, qcSp$y[49])

#gives list of residual to be subtracted/added to sample data in the correct
#sample order
# Apply the Cubic Spline to the relevant samples from the orginal data 
# View data points Raw$GLY [1:49]
GlyCorr <- Raw$GLY + resid
plot(GlyCorr, main = "Corrected Concentrations of Urine Samples")

#Now for 1-methylhistamine MH
MHQC <- raw[class == 2]$X1MH
QCname <- raw[class ==2]$Sample.Name


# create an LOESS curve with Optimized Span 

dfmetab2 <- data.frame (x=c(QCname),
                        y=c(MHQC))


metabloess75 <- loess(y ~ x, data=dfmetab2, span=.75)
smooth75 <- predict(metabloess75) 

metabloess90 <- loess(y ~ x, data=dfmetab2, span=.9)
smooth90 <- predict(metabloess90) 

# Plot the LOESS Curve(s) 

plot(dfmetab2$x, dfmetab2$y, pch=19,xlab = "Sample Name", ylab = "Concentration", main='Loess Regression Models for 1-MH QCs')
lines(smooth75, x=dfmetab2$x, col='purple')
points(smooth75, x=dfmetab2$x, col='purple')
lines(smooth90, x=dfmetab2$x, col='darkgreen')
points(smooth90, x=dfmetab2$x, col='darkgreen')
legend('bottomright', legend=c('.75', '.9'),
       col=c('purple', 'darkgreen'), pch=19, title='Smoothing Span')


# Pull out LOESS Curve Values for Spline 
#predict(metabloess90)

# Create a cubic spline interpolation between LOESS values for every datapt.

qcSp <- spline(QCname, predict(metabloess90), xmin = 1, xmax = 49, n = 49)
plot(qcSp, main= "Cubic Spline Interpolation of 1-Methyl Histamine",  
     xlab = "Sample Name", ylab = "Concentration")
points(smooth90, x=dfmetab$x, col='blue', pch=19)
lines(smooth90, x=dfmetab$x, col='blue')

resid <- c(qcSp$y[1], qcSp$y) - c(qcSp$y, qcSp$y[49])
#gives list of residual to be subtracted/added to sample data in that order

#determine the final corrected values for MH metabolite and plot them
# view data Raw$X1MH[1:49]
MHCorr <- Raw$X1MH + resid
plot(MHCorr, main= 'Plot of Corrected Samples for 1-Methylhistamine Metabolite')

#exact same process but with Arginine now

ARGQC <- raw[class == 2]$ARG

# create an LOESS curve with Optimized Span 

dfmetab3 <- data.frame (x=c(QCname),
                        y=c(ARGQC))

metabloess75 <- loess(y ~ x, data=dfmetab3, span=.75)
smooth75 <- predict(metabloess75) 

metabloess90 <- loess(y ~ x, data=dfmetab3, span=.9)
smooth90 <- predict(metabloess90) 

# Plot the LOESS Curve(s) 

plot(dfmetab3$x, dfmetab3$y, pch=19,xlab = "Sample Name", ylab = "Concentration", main='Loess Regression Models for Arginine QCs')
lines(smooth75, x=dfmetab3$x, col='purple')
points(smooth75, x=dfmetab2$x, col='purple')
lines(smooth90, x=dfmetab3$x, col='darkgreen')
points(smooth90, x=dfmetab2$x, col='darkgreen')
legend('bottomright', legend=c('.75', '.9'),
       col=c('purple', 'darkgreen'), pch=19, title='Smoothing Span')

# Pull out LOESS Curve Values for Spline 

#predict(metabloess90)

# Create a cubic spline interpolation between LOESS values for every datapt.

qcSp <- spline(QCname, predict(metabloess90), xmin = 1, xmax = 49, n = 49)
plot(qcSp, main= "Cubic Spline Interpolation of Arginine",
     xlab = "Sample Name", ylab = "Concentration")
points(smooth90, x=dfmetab$x, col='blue', pch=19)
lines(smooth90, x=dfmetab$x, col='blue')

resid <- c(qcSp$y[1], qcSp$y) - c(qcSp$y, qcSp$y[49])
#gives list of residual to be subtracted/added to sample data in that order

#make the final corrected values for Arginine

ArgCorr <- c(Raw$ARG + resid)
plot(ArgCorr, main= "Plot of Corrected Samples for Arginine Metabolite")





# Export fully corrected and normalized values 
#(from BOTH the TM and the Cubic Spline normalization)

ProcessedData <- data.frame (SN, GlyCorr, MHCorr, ArgCorr, Tm, class, 
                             stringsAsFactors = FALSE)
ProcessedData
sink('Processed LC-MS Data.csv')
write.csv (ProcessedData)
sink()
#this step will create a file in your working directory as a .csv as the 
#final corrected data



#now if you have several metabolites a loop can be used to automate correction


#here the values 2:4 can be adjusted in future depending on the number of 
#metabolites of interest and can be confirmed you are looking at the correct 
#columns by looking at the data table using the head() function

#loop to normalize hydration status 
datafile <- read.csv("SAMPLE LCMS DATAFile.csv")
#this is the trial sample set also available on GitHub

#pulling out the metabolite values from each column from .csv file

metabs <- colnames(datafile)[2:4]
#this will allow you to select each individual metabolite (generalized for expansion)
#metab <- metabs[#]

TmNormData <- datafile
for (metab in metabs) {
  SN <- (select(datafile,Sample.Name))
  class <- select (datafile, class)
  Tm <- select (datafile, TM)
  Metab <- select(datafile, all_of(metab))
  #creating the hydration status normalized values for each metabolite
  MetabTm <- Metab/Tm
  TmNormData[metab] <- MetabTm
}

sink( "Tm Norm Data")
write.csv(TmNormData)
sink

metab <- metabs[1]

#this is the loop for the Drift normalization
ProcessedData <- TmNormData
for (metab in metabs) {
  #this part is selecting the QC values by pulling samples with a class=2 value
  metabQC <- dplyr::filter(Raw, class == 2)[, metab]
  #this part is making the data for each metabolite value a DF
  dfMetab <- data.frame (x=c(QCname),
                       y=c(metabQC))
  #using this dataframe we will create an LOESS curve 
  #here we are using a generalized span=0.75 
  #further span parameters and optmization can be added as required
  Metabloess75 <- loess(y ~ x, data=dfMetab, span=.75)
  Metabsmooth75 <- predict(Metabloess75) 
  #this plots the LOESS Line on a graph
  plot(dfMetab$x, dfMetab$y, pch=19, xlab = "Sample Name", ylab = "Concentration",
       main= paste("LOESS regression of", metab, "Sample Data"))
  lines(Metabsmooth75, x=dfMetab$x, col='purple')
  points(smooth75, x=dfmetab2$x, col='purple')
  legend('bottomright', legend=c('.75'), pch=19, title='Smoothing Span')
  #this pulls out the new values from the optimized span LOESS for QCs
  predict(Metabloess75)
  #this part builds a cubic spline for the equivalent length of SAMPLE dataset
  qcSp <- spline(QCname, predict(Metabloess75), xmin = 1, xmax = 49, n = 49)
  #note the xmin and xmax values will correspond to the number of samples
  plot(qcSp, main =paste("Cubic Spline Interpolation for", metab),
       xlab = "Sample Name", ylab = "Concentration" )
  points(Metabsmooth75, x=dfmetab$x, col='pink', pch=19)
  lines(Metabsmooth75, x=dfmetab$x, col='darkgreen')
  #this part will determine the residuals for each point which can be applied to
  #the corresponding sample and used to correct the final data. 
  resid <- c(qcSp$y[1], qcSp$y) - c(qcSp$y, qcSp$y[49])
  resid <- resid[1:49]
  #then we will put it together to make the corrected sample amounts
  Metab <- select(TmNormData, all_of(metab))
  metabcorr <- Metab+resid
  ProcessedData[metab] <- metabcorr
}

ProcessedData
TmNormData
sink('Processed LC-MS Data.csv')
write.csv (ProcessedData)
sink()
#this step will create a file in your working directory as a .csv as the 
#final corrected data


#note the above does the same as the individual, but is applicable for the 
#addition of more metabolites to the loop by adjusting the column parameters
#also note that the spline section requires the exact number of samples in the 
#file to be able to work correctly
# if there is a specific metabolite of interest for the LOESS or Spline graphs
#these can be obtained by using the individual coding aspect from above the loop
# or running the lines of the loop individually for specific metabolite. 
# Note: the individual algorithms use a span =0.9 and the loop uses 0.75
#this was to demonstrate the effect of a smoothing parameter and will be
# optimized based on number and amount of variation in true data. 
#Pink and Green were chosen to differentiate between graphs that used a span 
# of 0.75 compared to those that used span=0.9 which is blue




#Statistical Processing: PCA plot generation 

# Create PCA plot 
library(ggplot2)
library(ggfortify)


data.matrix <- ProcessedData[c(1:6)]
results <- na.omit(data.matrix)

pca <- prcomp(results)
plot(pca$x[,1], pca$x[,2], xlab = "PC 1", ylab = "PC 2", main = "PCA of Processed LC-MS Data")

pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")
#lovely little bar plot which shows how much variation explained by each PC


pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])
pca.data

ggplot(data=pca.data, aes(x=X, y=Y, label=Sample)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("Corrected Sample PCA with Sample Name")
#this will give a plot where each point is labelled with the associated sample
#number and can be useful to ensure all points are on the graph

ggplot(data=pca.data, aes(x=X, y=Y, colour = as.factor(results$class))) +
  geom_point() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_classic() +
  ggtitle("Corrected Sample PCA")
#this graph will create a PCA where the points are coloured depending on the 
#value in the class column
#this is ideal for easy visualization of how data is being separated from the 
#model based on how they relate to class
#note you will want the QCs to be clustered together

