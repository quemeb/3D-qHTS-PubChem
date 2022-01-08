# THE PURPOSE OF THIS PROGRAM IS TO TAKE A .CSV FILE AND CREATE A 3D WATERFALL PLOT

#Important settings
keyword_1 <- "Active"
keyword_2 <- #"Inactive"
alpha_1 <- 1

#Color of points in graph ------------------------------------------------------
col_1 <- c("royalblue3") #active points
col_2 <- c("gray")   # noise  or non-actives
col_3 <- c("darkgree")   #keyword_1 curve color
col_4 <- c("darkgreen") #keyword_2 curve color
# all R colors can be found at http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf


#Checking for existence of RGL Package -----------------------------------------
avail <- installed.packages()
pack <- avail[,1]
if(is.element("rgl",pack) == FALSE){
  install.packages("rgl", type = "binary")
}
if(is.element("stringr", pack) == FALSE){
  install.packages("stringr")
}
if(is.element("dplyr", pack) == FALSE){
  install.packages("stringr")
}
if(is.element("tidyr", pack) == FALSE){
  install.packages("stringr")
}

#Loading libraries -------------------------------------------------------------
library(rgl)
library(stringr)
library(dplyr)
library(tidyr)

#Asking to select file ---------------------------------------------------------
ifile <- file.choose()

#Reading first line of CSV file to detect data type ----------------------------
cdata <- read.csv(ifile, header=TRUE, na.string="null",nrows = 1)
heads <- colnames(cdata)

#Checking if input file is from pubchem ----------------------------------------
if("PUBCHEM_ACTIVITY_OUTCOME" %in% heads){
  pubchem_data <- 1
}else{
  pubchem_data <- 0
}

#Identifying data headers to load ----------------------------------------------
if(pubchem_data == 1){
  #matching variable names
  fit <- "Fit_Output"
  readout <- "PUBCHEM_ACTIVITY_OUTCOME"
  lac50 <- "Fit_LogAC50"
  hill <- "Fit_HillSlope"
  inf <- "Fit_InfiniteActivity"
  zero <- "Fit_ZeroActivity"
  
  #finding the index to start reading
  finder <- read.csv(ifile, header = TRUE, na.strings = NULL, nrows = 10)
  skip <- which(finder[1]==1)
  cdata <- read.csv(ifile, col.names = heads, na.string="null", skip=skip-1)

  heads_len = length(heads)
  conc <- NULL
  conc_loc <- NULL
  a <- 1 #counting
  
  #Looking for concentration values
  for(i in 1:heads_len) {
    if(str_starts(heads[i], "Activity.at") == TRUE & 
       is.na(cdata[1,i]) == FALSE){
      conc[a] <- heads[i]
      conc_loc[a] <- i
      a <- a+1
    } 
  }
  titrations <- length(conc)
  x <- NULL
  y <- NULL
  
  #Extracting concentration values
  for(i in 1:titrations) {
    s <- unlist(str_split(conc[i], "\\."))
    x[i] <- str_c(s[3], ".", s[4])
    y[i] <- s[5]
  }
  
  #converting to numeric
  x <- as.numeric(x)
  
}else{
  cdata <- read.csv(ifile, header=TRUE, na.string="null")
  #matching variable names
  fit <- "Fit_Output"
  readout <- "Sample.Data.Type"
  lac50 <- "Log.AC50..M."
  hill <- "Hill.Coef"
  inf <- "Inf.Activity"
  zero <- "Zero.Activity"
}


#SETTING LOG MOLAR CONCENTRATIONS ----------------------------------------------
#if it is pubchem data, the program automatically finds and converts it
#if it is data from somewhere else, please input the data after the "else" 
#statement
#REMINDER: it is the log([concentration(M)])

if(pubchem_data == 1){
  conc <- NULL
  #log(molar concentration)
  l <- length(y)
  for(i in 1:l){
    conc[i] <- switch(y[i],
      "uM" = log(x[i]*1e-06),
      "nM" = log(x[i]*1e-09),
      "mM" = log(x[i]*1e-03),
      "pM" = log(x[i]*1e-12),
      "cM" = log(x[i]*1e-02),
      "dM" = log(x[i]*1e-01),
      "fM" = log(x[i]*1e-15))
  }
}else{ #please INPUT your concentration here if the data is NOT from PubChem
    conc <- c(
      -9.231,
      -9.11,
      -8,
      -7,
      -6,
      -5,
      -4)
}

conc <- format(round(conc, 2), nsmall = 2)
conc <- as.numeric(conc)

titrations <- length(conc)
lowerBound <- min(conc) 
upperBound <- max(conc)

#creating a Fit_Output if not given one ----------------------------------------
if(fit %in% heads == FALSE){
  l <- nrow(cdata)
  for (i in 1:l){
    if(cdata[i,readout]==keyword_1 || cdata[i,readout]==keyword_2){
      cdata[i,fit] = 1
    }else{
      cdata[i,fit] = 0
    }
  }
}

#creating smaller data sets with only needed data-------------------------------
if(pubchem_data == 1){
  cdata_points <- cdata %>% select(readout, conc_loc)
  cdata_curves <- cdata %>% select(fit, readout, lac50, hill, inf, zero)
}else{
  dataCols <- NULL
  for(i in 1:titrations){
    dataCols[i] <- c(paste("Data",(i-1), sep=""))
    }
  cdata_points <- cdata[,c(readout,dataCols)]
  cdata_curves <- cdata[,c(fit, readout, lac50, hill, inf, zero)]
}

names(cdata_curves) <- c("Fit_Output", "readout", "LAC50", "HILL", "INF", "ZERO")



#titration points to be put on the graph ---------------------------------------
l <- nrow(cdata_points)
for(i in 1:l){
  cdata_points$z[i] <- i
}
colnames(cdata_points) <- c("readout", conc, "z")

l <- length(colnames(cdata_points))
mainMatrix <- pivot_longer(cdata_points, cols = 2:(l-1), names_to = "x", values_to = "y")


#correcting data type-----------------------------------------------------------

mainMatrix$x <- as.double(mainMatrix$x)
mainMatrix$y <- as.double(mainMatrix$y)
mainMatrix$z <- as.double(mainMatrix$z)

waterfall_POINTS_data <- mainMatrix

#separating the data------------------------------------------------------------
waterfall_POINTS_data_1 <- data.frame(x=double(), y=double(), z=double())
waterfall_POINTS_data_2 <- data.frame(x=double(), y=double(), z=double())

l <- nrow(waterfall_POINTS_data)
for(i in 1:l){
  if(waterfall_POINTS_data$readout[i] == keyword_1 || waterfall_POINTS_data$readout[i] == keyword_2) {
    waterfall_POINTS_data_1 <- bind_rows(waterfall_POINTS_data_1, data.frame(x=waterfall_POINTS_data$x[i],
                                                                         y=waterfall_POINTS_data$y[i],
                                                                         z=waterfall_POINTS_data$z[i]))
  }else{
    waterfall_POINTS_data_2 <- bind_rows(waterfall_POINTS_data_2, data.frame(x=waterfall_POINTS_data$x[i],
                                                                         y=waterfall_POINTS_data$y[i],
                                                                         z=waterfall_POINTS_data$z[i]))
  }
  
}

#taking care of missing data ---------------------------------------------------
cdata_curves$LAC50[is.na(cdata_curves$LAC50)] <- log10(10)
cdata_curves$HILL[ is.na(cdata_curves$HILL)] <- 1 


#functions to recreate titration curves ----------------------------------------
interleave <- function(x) {
  unlist(lapply(1:(length(x)-1), function(i) c(x[i], x[i+1])))
}

f <- function(params, concs, interleave=TRUE) {
  xx <- seq(min(concs), max(concs), length=25)
  yy <- with(params, ZERO + (INF-ZERO)/(1 + 10^((LAC50-xx)*HILL)))
  return(data.frame(x=xx, y=yy))
}


#recreating titration curves keyword_1------------------------------------------
mainMatrix <- data.frame(x=double(),y=double(),z=double())
rowIndex = 0;
l <- nrow(cdata_curves)
for (i in 1:l) {
  if(cdata_curves[i,"Fit_Output"]==1 && cdata_curves[i,"readout"]==keyword_1) {
    rowIndex = rowIndex+1
    d1 <- data.frame(f(cdata_curves[i,], c(lowerBound, upperBound)), z=i)
    
    #add multiple rows
    mainMatrix <- bind_rows(mainMatrix, data.frame(x=d1[,1], z=i, y=d1[,2]))
    #needed for break mechanic when graphing
    mainMatrix <- bind_rows(mainMatrix, data.frame(x=NA, z=NA, y=1))
    
  }
}

#correcting data type if needed
mainMatrix$x <- as.double(mainMatrix$x)
mainMatrix$y <- as.double(mainMatrix$y)
mainMatrix$z <- as.double(mainMatrix$z)

waterfall_LINES_data_1 <- mainMatrix

#recreating titration curves keyword_2------------------------------------------
mainMatrix <- data.frame(x=double(),y=double(),z=double())
rowIndex = 0;
l <- nrow(cdata_curves)
for (i in 1:l) {
  if(cdata_curves[i,"Fit_Output"]==1 && cdata_curves[i,"readout"]==keyword_2) {
    rowIndex = rowIndex+1
    d1 <- data.frame(f(cdata_curves[i,], c(lowerBound, upperBound)), z=i)
    
    #add multiple rows
    mainMatrix <- bind_rows(mainMatrix, data.frame(x=d1[,1], z=i, y=d1[,2]))
    #needed for break mechanic when graphing
    mainMatrix <- bind_rows(mainMatrix, data.frame(x=NA, z=NA, y=1))
    
  }
}

#correcting data type if needed
mainMatrix$x <- as.double(mainMatrix$x)
mainMatrix$y <- as.double(mainMatrix$y)
mainMatrix$z <- as.double(mainMatrix$z)

waterfall_LINES_data_2 <- mainMatrix

#3D Graphing -------------------------------------------------------------------

#CHANGING POP-UP WINDOW PARAMETERS
# Changing parameters defaults from open3d() using par3d() for what I consider
# an initial good view of graph and window size
# The open3() window can be manually enlarged and the graph can be rotated 
# using user's mouse

# graph view coordinates
newMatrix = t(matrix(c(-0.6416085, 0.01792431, -0.7668231,  0, -0.1346225, 0.98157728,  0.1355843,  0,
                       0.7551258, 0.19022357, -0.6273754,  0, 0.0000000, 0.00000000,  0.0000000,  1), 
                     ncol = 4, nrow =4))

# Matrix was transposed since f(matrix) always reads y1,y2,y3,y4 
# col first then rows 

# Pop-up window size
newWindowRect = c(138, 161, 886, 760)

#SCALING THE 3D WINDOW ---------------------------------------------------------
#scale=c(concentration, %, #samples),  # 
open3d(userMatrix = newMatrix, windowRect = newWindowRect)

#PLOTTING POINTS IN 3D GRAPH ---------------------------------------------------
for(i in waterfall_POINTS_data_1$x)
{
  points3d(x=waterfall_POINTS_data_1$x[i], y = waterfall_POINTS_data_1$y[i],
           z=waterfall_POINTS_data_1$z[i], col = col_1,
           size = 2)
  break
}

for(i in waterfall_POINTS_data_2$x)
{
  points3d(x=waterfall_POINTS_data_2$x[i],y=waterfall_POINTS_data_2$y[i],
           z=waterfall_POINTS_data_2$z[i], col = col_2,
           size = 0.5)
  break
}

#PLOTTING LINES IN 3D GRAPH ----------------------------------------------------
for(i in waterfall_LINES_data_1$x)
{ 
  lines3d(waterfall_LINES_data_1[i],col = col_3, alpha = alpha_1)
  break
}

for(i in waterfall_LINES_data_2$x)
{ 
  lines3d(waterfall_LINES_data_2[i],col = col_3, alpha = alpha_1)
  break
}

#CREATE BOX AROUND EDGES FOR THE GRAPH -----------------------------------------
axes3d(expand = 1.03, box = FALSE, xunit = 'pretty', yunit = "pretty", zunit = 'pretty')

# Adding Grid lines
grid3d(c("x","y","z+"))

aspect3d(1,1,3)

#EXPORTING THE IMAGE FILES -----------------------------------------------------

## the following code is put as a comment to allow the user to position the interactive 
## graph however they please
## NAME THE FILE HOWEVER YOU LIKE WITHIN THE "" MARKS BELOW

# rgl.snapshot(filename = ".png")
# rgl.postscript(".svg", fmt="svg")

