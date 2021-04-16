# Let’s get started as normal. 
objects()
rm(list=objects())
if (!require("pacman")) install.packages("pacman")
pacman::p_load(httr,EcoHydRology,curl,data.table,multisensi)
setwd("~/Week11/")
list.files(all.files = T)
objects()   # Should be empty.
vignette("multisensi-vignette")

# Class example
verhulst <- function(K, Y0, a, t) {
  output <- K/(1 + (K/Y0 - 1) * exp(-a * t))
  return(output)
}
T <- seq(from = 5, to = 100, by = 5)
verhulst2 <- function(X, t = T) {
  out <- matrix(nrow = nrow(X), ncol = length(t), NA)
  for (i in 1:nrow(X)) {
    out[i, ] <- verhulst(X$K[i], X$Y0[i], X$a[i], t)
  }
  out <- as.data.frame(out)
  names(out) <- paste("t", t, sep = "")
  return(out)
}
n <- 10
set.seed(1234)
X <- data.frame(K = runif(n, min = 100, max = 1000), Y0 = runif(n, min = 1,
                                                                max = 40), a = runif(n, min = 0.05, max = 0.2))
Y <- verhulst2(X)
par(cex.axis = 0.7, cex.lab = 0.8)
plot(T, Y[1, ], type = "l", xlab = "Time", ylab = "Population size",
     ylim = c(0, 1000))
for (i in 2:n) {
  lines(T, Y[i, ], type = "l", col = i)
}
library(multisensi)
verhulst.seq <- multisensi(model=verhulst2, reduction=NULL, center=FALSE,
                           design.args = list( K=c(100,400,1000), Y0=c(1,20,40), a=c(0.05,0.1,0.2)))
## [*] Design
print(verhulst.seq, digits = 2)
plot(verhulst.seq, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)
title(xlab = "Time in half-decades.")
plot(verhulst.seq, normalized = FALSE, color = terrain.colors, gsi.plot = FALSE)
title(xlab = "Time in half-decades.")

################################################################# Lab 11
# Read the man page and look at the example for the Solar() function
# and become familiar with the variables and parameters passed in.
?EcoHydRology::Solar
#
# We start by defining your objective function, same as the “Use Case” on 
# page 2 of the vignette. Your function is defined, though you need to extend
# it for Julian days, and of course, Tmin cannot be greater than Tmax, so you 
# will have a Tmin and dT (with Tmax=Tmin+dT)
#
J <- seq(from = 1, to = 365, by = 5)
# Solar(lat, Jday, Tx, Tn, albedo=0.2, forest=0, slope=0, aspect = 0,
#      units="kJm2d")
# Note that the EcoHydRology::Solar() function is for specific days, 
# as such, we will want to create a function to loop through our period
# of interest:
Solar_Looped <- function(X, Jday = J) {
  out <- matrix(nrow = nrow(X), ncol = length(Jday), NA)
  for (i in 1:nrow(X)) {
    out[i, ] <- Solar(lat=X$lat[i],
                      Jday=Jday, Tx=X$Tx[i], 
                      Tn=(X$Tx[i]-X$Trange[i]), 
                      X$slope[i],X$aspect[i],units="Wm2")
  }
  out <- as.data.frame(out)
  names(out) <- paste("Jday", Jday, sep = "")
  return(out)
}
# A sample of with graph, from the vignette, we continue to build a 
# dataframe for our specific case with random uniform numbers for the 
# Tx, Tn (Tx - Trange), slope, and aspect.
# 
n <- 10
set.seed(1234)
X <- data.frame(Tx = runif(n, min = 5, max = 30), Trange = runif(n, min = 2,
                                                                 max = 16), slope = runif(n, min = 0.0, max = 0.2),
                aspect = runif(n, min = 0.0, max = 0.2),
                lat=runif(n, min = 0.0, max = 1.1))  # 1.1 radians lat is where?
# 
# Look at what you are passing into your new Solar_Looped() function
View(X)
#
Y <- Solar_Looped(X,Jday = J)
#
# You can ignore all the warnings, remember Errors=bad, warnings=not so much 
# So lets move on and build our summary graph
par(cex.axis = 0.7, cex.lab = 0.8)
plot(J, Y[1, ], type = "l", xlab = "Day of Year", 
     ylab = "Surface Short Wave Rad(W/m^2)")
for (i in 2:n) {
  lines(J, Y[i, ], type = "l", col = i)
}
# 
# Well, that is kewl, yet expected
#
# Multisensitivities
# 3 Sequential univariate sensitivity analyses
# 3.1 Calculation of sensitivity indices
Solar_Looped.seq <- multisensi(model=Solar_Looped, reduction=NULL, center=FALSE,
                               design.args = list( Tx = c(5,15,25), 
                                                   Trange = c(2,9,16), 
                                                   slope = c(0.1,0.2,0.3),
                                                   aspect = c(0.1,.5,1.0),
                                                   lat=c(0.1,.77,1.1)))

print(Solar_Looped.seq, digits = 2)
#
# 3.2 Graphical representation of sensitivity indices
#
dev.off() # Clean up previous par()
plot(Solar_Looped.seq, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)#normalized the upper subplot shows the extreme (tirets), #inter-quartile (grey) and median (bold line) output values
title(xlab = "Days of the Year.")
plot(Solar_Looped.seq, normalized = FALSE, color = terrain.colors, gsi.plot = FALSE)
title(xlab = "Days of the Year.")


#
# Take note of section 3.3 Calculating simulations apart, as there are 
# several ways to invoke the multisensi() function.
# First by building our X (design input) and Y (model run on X). Here 
# here is us doing the same sensitivities by passing previously run 
# models (Y)
X <- expand.grid(Tx = c(5,15,25), 
                 Trange = c(2,9,16), 
                 slope = c(0.1,0.2,0.3),
                 aspect = c(0.1,.5,1.0),
                 lat=c(0.1,.77,1.1))
# Look at our input
head(X,10)
Y <- Solar_Looped(X,Jday=J) # can be performed outside R if necessary
# Look at our model output
head(Y,10)
# Notice on the next line that “model=Solar_Looped” is replaced with our 
# model output “Y” and we add in the input into “design=X”. This 
# is exactly the same as above, though with the model run 
# external to the multisensi() function:
Solar_Looped.seq <- multisensi(design=X, model=Y, reduction=NULL, center=FALSE) 



# 4 Multivariate sensitivity analysis based on PCA
# read the vignette, though note we are using the multisensi() function
# to run our model (i.e. no “design” variable, and model=Solar_Looped)
Solar_Looped.pca <- multisensi(model=Solar_Looped, reduction=basis.ACP, scale=FALSE,
                               design.args = list( Tx = c(5,15,25), 
                                                   Trange = c(2,9,16), 
                                                   slope = c(0.1,0.2,0.3),
                                                   aspect = c(0.1,.5,1.0),
                                                   lat=c(0.1,.77,1.1)))

summary(Solar_Looped.pca, digits = 2)
# 4.2 Graphical representation for PCA based analysis with 
# explanation in vignette. These graphs require the plot window to be larger
# and might give "Error in plot.new() : figure margins too large". 
# If so expand the plot window.
dev.off()
plot(Solar_Looped.pca, graph = 1)
plot(Solar_Looped.pca, graph = 2)
plot(Solar_Looped.pca, graph = 3)
#
# 5.1 Polynomial reduction of the multivariate output
# Skip 5.1 Polynomial reduction for now and move on to
# 6 Alternative methods of sensitivity analysis
# 6.1 With Sobol2007 implemented in the package sensitivity
# 
library(sensitivity)
m <- 10000
Xb <- data.frame(Tx = runif(m, min = 5, max = 30), 
                 Trange = runif(m, min = 2,max = 16), 
                 slope = runif(m, min = 0.0, max = 0.2),
                 aspect = runif(m, min = 0.0, max = 0.2),
                 lat=runif(m, min = 0.0, max = 1.1))

Solar_Looped.seq.sobol <- multisensi(design = sobol2007, model = Solar_Looped,
                                     reduction = NULL, analysis = analysis.sensitivity, center = TRUE,
                                     design.args = list(X1 = Xb[1:(m/2), ], X2 = Xb[(1 + m/2):m, ], nboot = 100),
                                     analysis.args = list(keep.outputs = FALSE))
#
# Note, this is a good time time get a drink of water and/or pee as 
# it is running the function m=10,000 times (a few minutes).
#
print(Solar_Looped.seq.sobol, digits = 2)
dev.off()
plot(Solar_Looped.seq.sobol, normalized = TRUE, color = terrain.colors)
dev.off()  # this also cleans the graphics device. 
#
# 6.2 With fast99 implemented in the package sensitivity
#
Solar_Looped.seq.fast <- multisensi(design = fast99, model = Solar_Looped,
                                    center = FALSE, reduction = NULL, analysis = analysis.sensitivity,
                                    design.args=list( factors=c("Tx","Trange","slope","aspect","lat"), 
                                                      n=1000, q = "qunif",
                                                      q.arg = list(list(min=5, max=30), 
                                                                   list(min=2, max=16),
                                                                   list(min=0, max=.2),
                                                                   list(min=0, max=.2),
                                                                   list(min = 0.0, max = 1.1))),
                                    analysis.args=list(keep.outputs=FALSE))

print(Solar_Looped.seq.fast,digits=2)
plot(Solar_Looped.seq.fast, normalized = TRUE, color = terrain.colors)

########################## HW1
objects()
rm(list=objects())
J <- seq(from = 1, to = 365, by = 5)
# Solution for PET_fromTemp
# Trick is you have to notice that "lat_radians" has replaced "lat" and
# there is no "units" variable... and... notice that the function has to
# be fixed to allow Jday to be a vector of different size than Tmax and Tmin
PET_fromTemp <- function (Jday, Tmax_C, Tmin_C, lat_radians, AvgT = (Tmax_C + Tmin_C)/2, albedo = 0.18, TerrestEmiss = 0.97, aspect = 0, slope = 0, forest = 0, PTconstant=1.26, AEparams=list(vp=NULL, opt="linear"))
{
  cloudiness <- EstCloudiness(Tmax_C, Tmin_C)
  DailyRad <- NetRad(lat_radians, Jday, Tmax_C, Tmin_C, albedo, forest, slope, aspect, AvgT, cloudiness, TerrestEmiss, AvgT, AEparams=AEparams)
  potentialET <- PTpet(DailyRad, AvgT, PTconstant)
  potentialET[which(potentialET < 0)] <- 0
  potentialET[which(Tmax_C == -999 | Tmin_C == -999)] <- (-999)
  return(potentialET)
}
PET_fromTemp_Looped <- function(X, Jday = J) {
  out <- matrix(nrow = nrow(X), ncol = length(Jday), NA)
  for (i in 1:nrow(X)) {
    out[i, ] <- PET_fromTemp(lat_radians=X$lat_radians[i],
                             Jday=Jday, 
                             Tmax_C=X$Tmax_C[i],
                             Tmin_C=(X$Tmax[i]-X$Trange[i]),
                             X$slope[i],
                             X$aspect[i])
  }
  out <- as.data.frame(out)
  names(out) <- paste("Jday", Jday, sep = "")
  return(out)
}
n <- 10
set.seed(1234)
X <- data.frame(Tmax_C = runif(n, min = 5, max = 30),
                Trange = runif(n, min=2, max = 16),
                slope = runif(n, min = 0.0, max = 0.2),
                aspect = runif(n, min = 0.0, max = 0.2),
                lat_radians=runif(n, min = 0.0, max = 1.1))
Y <- PET_fromTemp_Looped(X,Jday = J)
#
# You can ignore all the warnings, remember Errors=bad, warnings=not so much 
# So lets move on and build our summary graph
par(cex.axis = 0.7, cex.lab = 0.8)
plot(J, Y[10, ], type = "l", xlab = "Day of Year", 
     ylab = "PET")
for (i in 2:n) {
  lines(J, Y[i, ], type = "l", col = i)
}
# Multisensitivities
# 3 Sequential univariate sensitivity analyses
# 3.1 Calculation of sensitivity indices
PET_fromTemp_Looped.seq <- multisensi(model=PET_fromTemp_Looped, reduction=NULL, center=FALSE,
                                      design.args = list( Tmax_C = c(5,15,25),
                                                          Trange = c(2,9,16),
                                                          slope = c(0.1,0.2,0.3),
                                                          aspect = c(0.1,.5,1.0),
                                                          lat_radians=c(0.1,.77,1.1)))

print(PET_fromTemp_Looped.seq, digits = 2)
#
# 3.2 Graphical representation of sensitivity indices
#
dev.off() # Clean up previous par()
plot(PET_fromTemp_Looped.seq, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)#normalized the upper subplot shows the extreme (tirets), #inter-quartile (grey) and median (bold line) output values
title(xlab = "Days of the Year.")
plot(PET_fromTemp_Looped.seq, normalized = FALSE, color = terrain.colors, gsi.plot = FALSE)
title(xlab = "Days of the Year.")

X <- expand.grid(Tmax_C = c(5,15,25),
                 Trange = c(2,9,16),
                 slope = c(0.1,0.2,0.3),
                 aspect = c(0.1,.5,1.0),
                 lat_radians=c(0.1,.77,1.1))
# Look at our input
head(X,10)
Y <- PET_fromTemp_Looped(X,Jday=J) # can be performed outside R if necessary
# Look at our model output
head(Y,10)
# Notice on the next line that “model=Solar_Looped” is replaced with our 
# model output “Y” and we add in the input into “design=X”. This 
# is exactly the same as above, though with the model run 
# external to the multisensi() function:
PET_fromTemp_Looped.seq <- multisensi(design=X, model=Y, reduction=NULL, center=FALSE) 

# 4 Multivariate sensitivity analysis based on PCA
# read the vignette, though note we are using the multisensi() function
# to run our model (i.e. no “design” variable, and model=Solar_Looped)
PET_fromTemp_Looped.pca <- multisensi(model=PET_fromTemp_Looped, reduction=basis.ACP, scale=FALSE,
                                      design.args = list( Tmax_C = c(5,15,25),
                                                          Trange = c(2,9,16),
                                                          slope = c(0.1,0.2,0.3),
                                                          aspect = c(0.1,.5,1.0),
                                                          lat_radians=c(0.1,.77,1.1)))

summary(PET_fromTemp_Looped.pca, digits = 2)
# 4.2 Graphical representation for PCA based analysis with 
# explanation in vignette. These graphs require the plot window to be larger
# and might give "Error in plot.new() : figure margins too large". 
# If so expand the plot window.
dev.off()
plot(PET_fromTemp_Looped.pca, graph = 1)
plot(PET_fromTemp_Looped.pca, graph = 2)
plot(PET_fromTemp_Looped.pca, graph = 3)

# 5.1 Polynomial reduction of the multivariate output
# Skip 5.1 Polynomial reduction for now and move on to
# 6 Alternative methods of sensitivity analysis
# 6.1 With Sobol2007 implemented in the package sensitivity
# 
library(sensitivity)
m <- 10000
Xb <- data.frame(Tmax_C = runif(m, min = 5, max = 30),
                 Trange = runif(m, min = 2, max = 16),
                 slope = runif(m, min = 0.0, max = 0.2),
                 aspect = runif(m, min = 0.0, max = 0.2),
                 lat_radians=runif(m, min = 0.0, max = 1.1))

PET_fromTemp_Looped.seq.sobol <- multisensi(design = sobol2007, model = PET_fromTemp_Looped,
                                            reduction = NULL, analysis = analysis.sensitivity, center = TRUE,
                                            design.args = list(X1 = Xb[1:(m/2), ], X2 = Xb[(1 + m/2):m, ], nboot = 100),
                                            analysis.args = list(keep.outputs = FALSE))
#
# Note, this is a good time time get a drink of water and/or pee as 
# it is running the function m=10,000 times (a few minutes).
#
print(PET_fromTemp_Looped.seq.sobol, digits = 2)
dev.off()
plot(PET_fromTemp_Looped.seq.sobol, normalized = TRUE, color = terrain.colors)
dev.off()  # this also cleans the graphics device. 

# 6.2 With fast99 implemented in the package sensitivity
#
PET_fromTemp_Looped.seq.fast <- multisensi(design = fast99, model = PET_fromTemp_Looped,
                                           center = FALSE, reduction = NULL, analysis = analysis.sensitivity,
                                           design.args=list( factors=c("Tmax_C","Trange","slope","aspect","lat_radians"), 
                                                             n=1000, q = "qunif",
                                                             q.arg = list(list(min=5, max=30),
                                                                          list(min=2, max=16),
                                                                          list(min=0, max=.2),
                                                                          list(min=0, max=.2),
                                                                          list(min = 0.0, max = 1.1))),
                                           analysis.args=list(keep.outputs=FALSE))

print(PET_fromTemp_Looped.seq.fast,digits=2)
plot(PET_fromTemp_Looped.seq.fast, normalized = TRUE, color = terrain.colors)

############################################################# NetRad()
# NetRad(lat, Jday, Tx, Tn, albedo = 0.18, forest = 0, slope = 0, aspect = 0, airtemp = (Tn+Tx)/2, cloudiness = "Estimate", surfemissivity = 0.97, surftemp = (Tn+Tx)/2, units = "kJm2d", AEparams=list(vp=NULL, opt="linear"))
NetRad_Looped <- function(X, Jday = J) {
  out <- matrix(nrow = nrow(X), ncol = length(Jday), NA)
  for (i in 1:nrow(X)) {
    out[i, ] <- NetRad(lat=X$lat[i],
                       Jday=Jday,
                       Tx=X$Tx[i], 
                       Tn=(X$Tx[i]-X$Trange[i]), 
                       X$slope[i],
                       X$aspect[i],
                       units="kJm2d")
  }
  out <- as.data.frame(out)
  names(out) <- paste("Jday", Jday, sep = "")
  return(out)
}
n <- 10
set.seed(1234)
X <- data.frame(Tx = runif(n, min = 5, max = 30),
                Trange = runif(n, min = 2, max = 16), 
                slope = runif(n, min = 0.0, max = 0.2),
                aspect = runif(n, min = 0.0, max = 0.2),
                lat=runif(n, min = 0.0, max = 1.1))
Y <- NetRad_Looped(X,Jday = J)
#
# You can ignore all the warnings, remember Errors=bad, warnings=not so much 
# So lets move on and build our summary graph
par(cex.axis = 0.7, cex.lab = 0.8)
plot(J, Y[1, ], type = "l", xlab = "Day of Year", 
     ylab = "Net Rad (kJm2d)")
for (i in 2:n) {
  lines(J, Y[i, ], type = "l", col = i)
}
# Multisensitivities
NetRad_Looped.seq <- multisensi(model=NetRad_Looped, reduction=NULL, center=FALSE,
                                design.args = list( Tx = c(5,15,25), 
                                                    Trange = c(2,9,16), 
                                                    slope = c(0.1,0.2,0.3),
                                                    aspect = c(0.1,.5,1.0),
                                                    lat=c(0.1,.77,1.1)))

print(NetRad_Looped.seq, digits = 2)
#
# 3.2 Graphical representation of sensitivity indices
#
dev.off() # Clean up previous par()
plot(NetRad_Looped.seq, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)#normalized the upper subplot shows the extreme (tirets), #inter-quartile (grey) and median (bold line) output values
title(xlab = "Days of the Year.")
plot(NetRad_Looped.seq, normalized = FALSE, color = terrain.colors, gsi.plot = FALSE)
title(xlab = "Days of the Year.")

X <- expand.grid(Tx = c(5,15,25), 
                 Trange = c(2,9,16), 
                 slope = c(0.1,0.2,0.3),
                 aspect = c(0.1,.5,1.0),
                 lat=c(0.1,.77,1.1))
Y <- NetRad_Looped(X,Jday=J) # can be performed outside R if necessary
NetRad_Looped.seq <- multisensi(design=X, model=Y, reduction=NULL, center=FALSE) 

# 4 Multivariate sensitivity analysis based on PCA
NetRad_Looped.pca <- multisensi(model=NetRad_Looped, reduction=basis.ACP, scale=FALSE,
                                design.args = list( Tx = c(5,15,25), 
                                                    Trange = c(2,9,16), 
                                                    slope = c(0.1,0.2,0.3),
                                                    aspect = c(0.1,.5,1.0),
                                                    lat=c(0.1,.77,1.1)))

summary(NetRad_Looped.pca, digits = 2)
# 4.2 Graphical representation for PCA based analysis
dev.off()
plot(NetRad_Looped.pca, graph = 1)
plot(NetRad_Looped.pca, graph = 2)
plot(NetRad_Looped.pca, graph = 3)

# 6.1 With Sobol2007 implemented in the package sensitivity
# 
library(sensitivity)
m <- 10000
Xb <- data.frame(Tx = runif(m, min = 5, max = 30), 
                 Trange = runif(m, min = 2,max = 16), 
                 slope = runif(m, min = 0.0, max = 0.2),
                 aspect = runif(m, min = 0.0, max = 0.2),
                 lat=runif(m, min = 0.0, max = 1.1))

NetRad_Looped.seq.sobol <- multisensi(design = sobol2007, model = NetRad_Looped,
                                      reduction = NULL, analysis = analysis.sensitivity, center = TRUE,
                                      design.args = list(X1 = Xb[1:(m/2), ], X2 = Xb[(1 + m/2):m, ], nboot = 100),
                                      analysis.args = list(keep.outputs = FALSE))
print(NetRad_Looped.seq.sobol, digits = 2)
dev.off()
plot(NetRad_Looped.seq.sobol, normalized = TRUE, color = terrain.colors)
dev.off()  # this also cleans the graphics device. 
#
# 6.2 With fast99 implemented in the package sensitivity
#
NetRad_Looped.seq.fast <- multisensi(design = fast99, model = NetRad_Looped,
                                     center = FALSE, reduction = NULL, analysis = analysis.sensitivity,
                                     design.args=list( factors=c("Tx","Trange","slope","aspect","lat"), 
                                                       n=1000, q = "qunif",
                                                       q.arg = list(list(min=5, max=30), 
                                                                    list(min=2, max=16),
                                                                    list(min=0, max=.2),
                                                                    list(min=0, max=.2),
                                                                    list(min = 0.0, max = 1.1))),
                                     analysis.args=list(keep.outputs=FALSE))

print(NetRad_Looped.seq.fast,digits=2)
plot(NetRad_Looped.seq.fast, normalized = TRUE, color = terrain.colors)
