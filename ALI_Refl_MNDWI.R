# File:   ALI_Refl_MNDWI.R                                  ("`-''-/").___..--''"`-._
# Author: George Panteras (gxp37@psu.edu)                   (`6_ 6  )   `-.  (     ).`-.__.`)
#         Geoinformatics and Earth Observation Laboratory   (_Y_.)'  ._   )  `._ `. ``-..-' 
#         (http://geoinf.psu.edu)                             _ ..`--'_..-_/  /--'_.' ,'
#         Department of Geography                           (il),-''  (li),'  ((!.-'
#         The Pennsylvania State University                 WE ARE ... PENN STATE!               
#    
# Functionality: - Radiance to reflectance conversion
#                - Modified normalized-difference water index (MNDWI) for EO-1 ALI.
# 
# The following NDWI is based on Li et al. 2013.
# 
# Hint: Li, W.; Du, Z.; Ling, F.; Zhou, D.; Wang, H.; Gui, Y.; Sun, B.; Zhang, X.	
#       A Comparison of Land Surface Water Mapping Using the Normalized Difference 
#       Water Index from TM, ETM+ and ALI. Remote Sens. 2013, 5, 5530-5549.
#
# Also: Chander, G., Markham, B.L., and Helder, D.L. (2009). Summary of Current 
#       Radiometric Calibration Coefficients for Landsat MSS, TM, ETM+, and EO-1 
#       ALI Sensors. Remote Sensing of Environment 113 (2009) 893–903.
################################################################################instal

########### libraries / packages

# install.packages("raster")
library(raster)

########### Load Data 

# Assign the pathname
path = "~/Downloads/ALI/EO1A0160372015283110KF"

## Create a RasterStack of the image 
# exclude the PAN since we cannot stack bands of different spatial res
img = list.files (path, pattern = '*Refl.tif$') [1:9]
img # check that you have only the (9) appropriate bands

# Stack the 9 bands in one
ali.img = stack(img)

# Check the layers
ali.img@layers
# Assign band names to each layer
names(ali.img) = c('b2','b3','b4','b5','b6','b7','b8','b9', 'b10')

########### Read Data and create the necessary properties for ALI

# Solar Spectral Irradiances. Bands 1-10 respectively (W/m^2*um). Modified from Chander et al., 2009
esun = c(1724, 1857, 1996, 1807, 1536, 1145, 955.8, 452.3, 235.1, 82.38)

# Scaling factor. Bands 1-10 respectively. Modified from Chander et al., 2009; USGS EO-1 (2011)
Grescale = c(0.045, 0.043, 0.028, 0.018, 0.011, 0.0091, 0.0083, 0.0028, 0.00091)

# Offset. Bands 1-10 respectively. Modified from Chander et al., 2009; USGS EO-1 (2011)
Brescale = c(-3.4, -4.4, -1.9, -1.3, -0.85, -0.65, -1.3, -0.6, -0.21)

metExtract = function( met, val ) {
  return( as.numeric(met[grep(val,met)[1]+2]))
}

# Basic conversions
deg2rad <- function(deg) return(deg*pi/180)
rad2deg <- function(rad) return(rad*180/pi)

## Read the data
filename = dir(path=path,pattern="MTL_L1GST.TXT",full.names=T)
met      = scan(filename,"char")

date     = met[grep("ACQUISITION_DATE",met)[1]+2]
sunelev  = metExtract(met, "SUN_ELEVATION")

ESdist = function(adate) 
{
  edist <- julian(as.Date(adate), origin = 
                    as.Date(paste(substring(adate, 1, 4), "12", "31", sep = "-")))[[1]]
  edist <- 1 - 0.016729 * cos((2 * pi) * (0.9856 * (edist - 4)/360))
  return( edist )
}

########### Radiance to reflectance

## Calculate the Reflectance
ali.reflectance = function( band, Grescale, Brescale, date, esun, sunelev ) {
  val = values(band)
  val[val==0]=NA  # to remove the frame around the layer
  val = ( val * Grescale + Brescale ) 
  if ( !is.na(esun) ) {
    val = pi * val * ( (ESdist(date)^2) / (esun * sin(deg2rad(sunelev))) ) 
  }
  val[val<0]=0
  
  values(band) = val
  
  return(band)
}

## Convert all the bands from radiance to reflectance

# Initialize an list to store all the reflectance images
img.list = list()

for (i in 1:9){
  img.list[[i]] = ali.reflectance(ali.img[[i]], Grescale[i], Brescale[i], date, esun[i], sunelev)
  
  # Optional: Comment the following line if you don't want to write each returned reflectance image 
  writeRaster(refl.b,file=paste0(strtrim(img[2], 24),i+1,"_Refl.TIF"),format="GTiff", overwrite=TRUE)
}

# Constract a raster brick with the converted bands
img.brick = brick(img.list)


########### Water Index

## MNDWI calculation
# The input should be the reflectance bands 4 and 6 corresponding to Green and NIR
ali.mndwi = (img.brick[[3]] - img.brick[[5]]) / (img.brick[[3]] + img.brick[[5]])

## Alternative valid band combinations for MNDWI
# mndwi = (b4 - b7) / (b4 + b7)   # mndwi = (b4 - b8) / (b4 + b8)
# mndwi = (b4 - b9) / (b4 + b9)   # mndwi = (b4 - b10) / (b4 + b10)

# Write the MNDWI image
writeRaster(ali.mndwi,filename=paste0(strtrim(img[2], 22),"_MNDWI.TIF"),format="GTiff", overwrite=TRUE)

