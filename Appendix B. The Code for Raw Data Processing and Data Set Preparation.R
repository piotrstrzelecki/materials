# 1. Required packages and functions --------------------------------------

# Install libraries; remove '#' to execute the line
# install.packages("DescTools")

# Load libraries
library(DescTools)

# Load functions
RotationMatrix <- function(x, y, z) {
  # Creates rotation matrix for 3D vectors in Euclidean space
  # x - rotation around x axis in degrees
  # y - rotation around y axis in degrees
  # z - rotation around z axis in degrees
  x <- -x / 180 * pi
  y <- -y / 180 * pi
  z <- -z / 180 * pi
  R <- cbind(c(1, 0, 0), c(0, cos(x), sin(x)), c(0, -sin(x), cos(x))) %*%
    cbind(c(cos(y), 0, -sin(y)), c(0, 1, 0), c(sin(y), 0, cos(y))) %*%
    cbind(c(cos(z), sin(z), 0), c(-sin(z), cos(z), 0), c(0, 0, 1))
  return(R)
}

ImageJtoGeo <- function(x, y, z) {
  # Function for coordinate system conversion:
  # ImageJ (Cartesian) to Geographical (spherical) coordinates
  # x - the x component of a 3D vector
  # y - the y component of a 3D vector
  # z - the z component of a 3D vector
  temp <- as.data.frame(CartToSph(x, y, z, up = T))
  temp[, 2:3] <- temp[, 2:3] * 180 / pi
  for (i in 1:nrow(temp)) {
    if (temp$theta[i] < 0) {
      if (temp$phi[i] < 90) {
        temp$phi[i] <- temp$phi[i] + 270
      }
      else {
        temp$phi[i] <- temp$phi[i] - 90
      }
    }
    else {
      temp$phi[i] <- temp$phi[i] + 90
    }
  }
  temp$theta <- abs(temp$theta)
  return(temp[, 2:3])
}

# 2. Load data ------------------------------------------------------------

# Set the directory
setwd("C:/Folder/...") # COMPLETE !

# Load raw data
Intensity <- read.csv("intensity.csv", header = T, sep = ",", dec = ".")
Region <- read.csv("regions.csv", header = T, sep = ",", dec = ".")
Geodesic <- read.csv("geodesic.csv", header = T, sep = ",", dec = ".")
Vectors <- read.csv("eigenvectors.csv", header = T, sep = ",", dec = ".")
Measure <- read.csv("measure.csv", header = T, sep = ",", dec = ".")

# Set a resolution of the measurements. To preserve pixel units set to 1
resolution <- 1 # e.g. 1 pix = 0.012 mm

# 3. Prepare a final dataset ----------------------------------------------

# (OPTIONAL) Rotation of components.
# Relevant if different position of the sample orientation is required
RZ <- RotationMatrix(0, 0, 0) # Rotation around z axis in ImageJ
RY <- RotationMatrix(0, 0, 0) # Rotation around y axis in ImageJ
RX <- RotationMatrix(0, 0, 0) # Rotation around y axis in ImageJ
# Rotation
Vectors[, 2:4] <- as.matrix(Vectors[, 2:4]) %*% RZ %*% RY %*% RX # rotate semi-axis a
Vectors[, 5:7] <- as.matrix(Vectors[, 5:7]) %*% RZ %*% RY %*% RX # rotate semi-axis b
Vectors[, 8:10] <- as.matrix(Vectors[, 8:10]) %*% RZ %*% RY %*% RX# rotate semi-axis c

# Prepare the final dataset (parameters in pixel units)
data <- data.frame(row.names = 1:nrow(Intensity))[1:nrow(Intensity), ]
{
  data$`mean intensity` <- Intensity$Mean
  data$`intensity sd` <- Intensity$StdDev
  data$`max intensity` <- Intensity$Max
  data$`min intensity` <- Intensity$Min
  data$`median intensity` <- Intensity$Median
  data$`mode intensity` <- Intensity$Mode
  data$`intensity skewness` <- Intensity$Skewness
  data$`intensity kurtosis` <- Intensity$Kurtosis
  data$`volume object` <- Measure$Vol..pix.
  data$`volume ellipsoid` <- Measure$Vol..pix. / Measure$RatioVolEllipsoid
  data$`volume box` <- Measure$VolBounding..pix.
  data$`volume ball` <- Geodesic$Radius^3 * 4 / 3 * pi
  data$`surface area` <- Measure$SurfCorr..pix.
  data$`Feret diameter` <- Measure$Feret..unit. / resolution
  data$`mean breadth` <- Region$MeanBreadth / resolution
  data$`geodesic diameter` <- Geodesic$Geod..Diam.
  data$`ball radius` <- Geodesic$Radius
  data$`semi-axis a` <- Region$Elli.R1 / resolution
  data$`semi-axis b` <- Region$Elli.R2 / resolution
  data$`semi-axis c` <- Region$Elli.R3 / resolution
  data$`mean radius` <- replace(Measure$DCMean..unit., is.na(Measure$DCMean..unit.), 0) / resolution
  data$`radius sd` <- replace(Measure$DCSD..unit., is.na(Measure$DCSD..unit.), 0) / resolution
  data$`max radius` <- replace(Measure$DCMax..unit., is.na(Measure$DCMax..unit.), 0) / resolution
  data$`min radius` <- replace( Measure$DCMin..unit., is.na(Measure$DCMin..unit.), 0) / resolution
  data$`compactness` <- Measure$CompCorr..pix.
  data$`discrete compactness` <- Measure$CompDiscrete
  data$`sphericity` <- Measure$SpherCorr..pix.
  data$`Euler number` <- Region$EulerNumber
  data$`geodesic elongation` <- Geodesic$Geod..Elong.
  data$`semi-axes ratio a/b` <- Region$Elli.R1.R2
  data$`semi-axes ratio a/c` <- Region$Elli.R1.R3
  data$`semi-axes ratio b/c` <- Region$Elli.R2.R3
  data$`object vol./ellipsoid vol. ratio` <- Measure$RatioVolEllipsoid
  data$`object vol./box vol. ratio` <- Measure$RatioVolbox
  data$`ball vol./object vol. ratio` <- Geodesic$Radius^3 * 4 / 3 * pi / Measure$Vol..pix.
  data$`x coordinate` <- Measure$CX..pix.
  data$`y coordinate` <- Measure$CY..pix.
  data$`z coordinate` <- Measure$CZ..pix.
  data$`semi-axis a trend` <- ImageJtoGeo(Vectors[, 2], Vectors[, 3], Vectors[, 4])[, 2]
  data$`semi-axis a plunge` <- ImageJtoGeo(Vectors[, 2], Vectors[, 3], Vectors[, 4])[, 1]
  data$`semi-axis b trend` <- ImageJtoGeo(Vectors[, 5], Vectors[, 6], Vectors[, 7])[, 2]
  data$`semi-axis b plunge` <- ImageJtoGeo(Vectors[, 5], Vectors[, 6], Vectors[, 7])[, 1]
  data$`semi-axis c trend` <- ImageJtoGeo(Vectors[, 8], Vectors[, 9], Vectors[, 10])[, 2]
  data$`semi-axis c plunge` <- ImageJtoGeo(Vectors[, 8], Vectors[, 9], Vectors[, 10])[, 1]
}

# Save the dataset in pixel units
write.csv(data, "dataset_pixel.csv", row.names = F)

# Save the dataset in desired unit
resolution <- 0.012
unit <- "mm"

# Prepare the final dataset (parameters in desired units)
data[, 9:12] <- data[, 9:12] * resolution^3
data[, 13] <- data[, 13] * resolution^2
data[, c(14:24, 36:38)] <- data[, c(14:24, 36:38)] * resolution

# Save the dataset in mm units
write.csv(data, paste0("dataset_", unit, ".csv"), row.names = F)

# Clear the workspace
rm(list = ls())
gc()