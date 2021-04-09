list.files(all.files = T)
objects()   # Should be empty.
#

if (!require("pacman")) install.packages("pacman")
pacman::p_load(httr,EcoHydRology,curl,elevatr,raster,rgdal,
               data.table,foreign,maptools,dataRetrieval,gdistance)
##############################################
# 0205551460 LICK RUN ABOVE PATTON AVENUE AT ROANOKE, VA
##############################################
# Make a function:
make_usgs_gage_list=function(siteNo = "0205551460", parameterCd = c("00060","00065"), start.date = "2017-05-01", end.date = "2017-11-01")
{
  #siteNo = "0205551460"
  #parameterCd = c("00060","00065")
  #start.date = "2017-05-01"  # Not frozen to not frozen
  #end.date = "2017-11-01"    # to still not frozen
  
  gagelist=list()   # Organize the data in a nice list as in previous labs
  gagelist[["flowdata"]]<- readNWISuv(siteNumbers = siteNo,parameterCd = parameterCd,startDate = start.date,endDate = end.date)
  #head(gagelist$flowdata)
  # Convert units.
  gagelist$flowdata$depth_m=gagelist$flowdata$X_00065_00000*0.3048
  # m/ft depth
  gagelist$flowdata$cms=gagelist$flowdata$X_00060_00000*.02832
  # m3/ft3 flow
  # Let's add in the USGS gage site information to the list and inspect
  gagelist[["site"]]=readNWISsite(siteNo)
  head(gagelist$site)
  class(gagelist$site$dec_lat_va)
  #
  # Set the Manning Coefficient in the USGS Gage's Site Table
  #
  gagelist$site$man_n=.035/1.49
  #
  # Create a SpatialPointsDataFrame out of the site dataframe in the USGS list
  coordinates(gagelist$site)=~dec_long_va+dec_lat_va
  
  return (gagelist)
}

# Run the function for the other sites.
USGS0205551460=make_usgs_gage_list(siteNo= "0205551460")
USGS02056000=make_usgs_gage_list(siteNo= "02056000")
USGS02055100=make_usgs_gage_list(siteNo= "02055100")      
USGS02054530=make_usgs_gage_list(siteNo= "02054530")
USGS02055000=make_usgs_gage_list(siteNo= "02055000")

# DEM data
ab_ll=rbind(USGS02056000$site,
            USGS0205551460$site,
            USGS02055100$site,
            USGS02055000$site,
            USGS02054530$site)
class(ab_ll)
ab_ll@proj4string
proj4_utm = paste0("+proj=utm +zone=",
                   trunc((180+coordinates(USGS02055000$site)[1])/6+1), 
                   " +datum=WGS84 +units=m +no_defs")
print(proj4_utm)
# Lat/Lon (_ll) is much easier!
proj4_ll = "+proj=longlat"
crs_ll=CRS(proj4_ll)
crs_utm=CRS(proj4_utm)
proj4string(ab_ll)=proj4_ll
ab_utm=spTransform(ab_ll,crs_utm)
ab_utm@coords
mydem=get_aws_terrain(locations=ab_utm@coords, 
                      z = 12, prj = proj4_utm,expand=1)
#
# Lets plot the DEM and the gage locations so we can guess 
# what gages connect with what gages
#
plot(mydem)
plot(ab_utm,add=T)
text(ab_utm, labels=ab_utm@data$site_no, cex=0.6, font=2,pos=1)
# From Lab02, I know I can get an overview of streams with the 
# USGS H
url="https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/NHD/HU8/HighResolution/Shape/NHD_H_03010101_HU8_Shape.zip"
curl_download(url,"NHD_H_03010101_HU8_Shape.zip")
unzip("NHD_H_03010101_HU8_Shape.zip",exdir="03010101")
streams=readOGR("03010101/Shape/NHDFlowline.dbf")
streams_utm=spTransform(streams,crs_utm)
plot(streams_utm,col="blue",add=T)


################################################### USGS0205551460
# USGS gage height for USGS02056000
USGS02056000$flowdata=USGS02056000$flowdata[,c(1,2,3,4,5,8,10)]
USGS02056000[["rating"]]=readNWISrating(USGS02056000$site$site_no)
plot(USGS02056000$rating$DEP,USGS02056000$rating$INDEP,xlab="DEP",ylab="INDEP")

USGS02056000$flowdata$X_00065_00000=approx(USGS02056000$rating$DEP,
                                           USGS02056000$rating$INDEP, xout = USGS02056000$flowdata$X_00060_00000, ties = min)$y
points(USGS02056000$flowdata$X_00060_00000,USGS02056000$flowdata$X_00065_00000,
       col="red")
#
USGS02056000$flowdata$depth_m=(USGS02056000$flowdata$X_00065_00000*0.3048)-0.8
# m/ft depth

# Distance between points.
A=SpatialPoints(USGS0205551460$site)# Up gradient site Lick Run
B=SpatialPoints(USGS02056000$site) # Down gradient site ROA River atNiagara
proj4string(A)=proj4_ll
proj4string(B)=proj4_ll
A_utm=spTransform(A,crs_utm)
B_utm=spTransform(B,crs_utm)
# Cut the DEM down to a more manageable size
cropmydem=crop(mydem,extend(extent(ab_utm),600))
cropmydem=trim(cropmydem)
cropmydem=cropmydem*1000.0
plot(cropmydem)
plot(ab_utm,add=T)
# Set up the weighting functions
altDiff <- function(x){x[2] - x[1]}
hd <- transition(cropmydem, altDiff, 8, symm=FALSE)
slope <- geoCorrection(hd)
adj <- adjacent(cropmydem, cells=1:ncell(cropmydem), pairs=TRUE, directions=8)
speed <- slope
speed[adj] <- 6 * exp(-3.5 * abs(slope[adj] + 0.05))
Conductance <- geoCorrection(speed)
# Find and plot the flow path
AtoB <- shortestPath(Conductance, A_utm, B_utm, output="SpatialLines")
plot(AtoB,add=T)
plot(streams_utm,col="blue",add=T)
plot(AtoB,add=T)
SpatialLinesLengths(AtoB)
USGS0205551460$site$L=SpatialLinesLengths(AtoB) # km to m
USGS0205551460$site$L # reach length in m

# Slope of stream
USGS0205551460$site$slope=(extract(mydem,A_utm)-
                             extract(mydem,B_utm))/USGS0205551460$site$L
USGS0205551460$site$slope

USGS0205551460$flowdata$B=(USGS0205551460$site$man_n*
                             USGS0205551460$flowdata$cms)/(USGS0205551460$flowdata$depth_m^(5/3)*
                                                             sqrt(USGS0205551460$site$slope))
head(USGS0205551460$flowdata)
# Lets look at how B changes with flow.    
plot(USGS0205551460$flowdata$dateTime,USGS0205551460$flowdata$B, main="LICK RUN TO ROANOKE RIVER AT NIAGARA, VA")
# Does this seem reasonable (...like order of magnitude reasonable)? You can 
# perform a quick and dirty check using google earth and measuring the channel 
# width in a few places.
#
plot(USGS0205551460$flowdata$cms,USGS0205551460$flowdata$depth_m, main="LICK RUN TO ROANOKE RIVER AT NIAGARA, VA")

# ck
# USGS0205551460$flowdata$ck = ???
  # ANS
  USGS0205551460$flowdata$ck =
  5/3*sqrt(USGS0205551460$site$slope)/USGS0205551460$site$man_n*
  (USGS0205551460$flowdata$depth_m^(2/3))
#
# USGS0205551460$flowdata$dt = ???
  USGS0205551460$flowdata$dt =
  USGS0205551460$site$L/USGS0205551460$flowdata$ck

plot(USGS0205551460$flowdata$dateTime,USGS0205551460$flowdata$dt)
USGS0205551460$flowdata$outTime=USGS0205551460$flowdata$dateTime+
  USGS0205551460$flowdata$dt

# Find beginning of  Waves
USGS0205551460$flowdata$newwave=
  USGS0205551460$flowdata$cms *1.1 <
  data.table::shift(USGS0205551460$flowdata$cms)
summary(USGS0205551460$flowdata$newwave)
# Add plot of the point found
len=length(USGS0205551460$flowdata$newwave)
USGS0205551460$flowdata$newwave[is.na(USGS0205551460$flowdata$newwave)]=F
# Removes repeated finds by going through loop backwords
for (i in seq(len,2)){
  print(i)
  if(USGS0205551460$flowdata$newwave[i]==T &
     USGS0205551460$flowdata$newwave[i-1]==T){
    USGS0205551460$flowdata$newwave[i]=F
  }
}
plot(USGS0205551460$flowdata$dateTime,USGS0205551460$flowdata$cms,type="l")
points(USGS0205551460$flowdata$dateTime[USGS0205551460$flowdata$newwave],
       USGS0205551460$flowdata$cms[USGS0205551460$flowdata$newwave],col=2)

# Find the time locations where waves begin
which(USGS0205551460$flowdata$newwave == TRUE)
plot(USGS0205551460$flowdata$dateTime,USGS0205551460$flowdata$cms,
     type="l",xlim=c(USGS0205551460$flowdata$dateTime[1109],
                     USGS0205551460$flowdata$dateTime[1109+200]), xlab="time", ylab="flow (m3/s)")
lines(USGS0205551460$flowdata$outTime,USGS0205551460$flowdata$cms,col=2)


################################################### USGS02055100
# USGS gage height for USGS02056000
# Distance between points.
A=SpatialPoints(USGS02055100$site)# Up gradient site
B=SpatialPoints(USGS02056000$site) # Down gradient site ROA River atNiagara
proj4string(A)=proj4_ll
proj4string(B)=proj4_ll
A_utm=spTransform(A,crs_utm)
B_utm=spTransform(B,crs_utm)
# Cut the DEM down to a more manageable size
cropmydem=crop(mydem,extend(extent(ab_utm),600))
cropmydem=trim(cropmydem)
cropmydem=cropmydem*1000.0
# Set up the weighting functions
altDiff <- function(x){x[2] - x[1]}
hd <- transition(cropmydem, altDiff, 8, symm=FALSE)
slope <- geoCorrection(hd)
adj <- adjacent(cropmydem, cells=1:ncell(cropmydem), pairs=TRUE, directions=8)
speed <- slope
speed[adj] <- 6 * exp(-3.5 * abs(slope[adj] + 0.05))
Conductance <- geoCorrection(speed)
# Find and plot the flow path
AtoB <- shortestPath(Conductance, A_utm, B_utm, output="SpatialLines")
plot(AtoB,add=T)
plot(streams_utm,col="blue",add=T)
plot(AtoB,add=T)
SpatialLinesLengths(AtoB)
USGS02055100$site$L=SpatialLinesLengths(AtoB) # km to m
USGS02055100$site$L # reach length in m

# Slope of stream
USGS02055100$site$slope=(extract(mydem,A_utm)-
                           extract(mydem,B_utm))/USGS02055100$site$L
USGS02055100$site$slope

USGS02055100$flowdata$B=(USGS02055100$site$man_n*
                           USGS02055100$flowdata$cms)/(USGS02055100$flowdata$depth_m^(5/3)*
                                                         sqrt(USGS02055100$site$slope))
head(USGS02055100$flowdata)
# Lets look at how B changes with flow.    
plot(USGS02055100$flowdata$dateTime,USGS02055100$flowdata$B, main="LICK RUN TO ROANOKE RIVER AT NIAGARA, VA")
# Does this seem reasonable (...like order of magnitude reasonable)? You can 
# perform a quick and dirty check using google earth and measuring the channel 
# width in a few places.
#
plot(USGS02055100$flowdata$cms,USGS02055100$flowdata$depth_m, main="LICK RUN TO ROANOKE RIVER AT NIAGARA, VA")


# ck
# USGS02055100$flowdata$ck = ???
# ANS
USGS02055100$flowdata$ck =
  5/3*sqrt(USGS02055100$site$slope)/USGS02055100$site$man_n*
  (USGS02055100$flowdata$depth_m^(2/3))
#
# USGS02055100$flowdata$dt = ???
USGS02055100$flowdata$dt =
  USGS02055100$site$L/USGS02055100$flowdata$ck

plot(USGS02055100$flowdata$dateTime,USGS02055100$flowdata$dt)
USGS02055100$flowdata$outTime=USGS02055100$flowdata$dateTime+
  USGS02055100$flowdata$dt

# Find beginning of  Waves
USGS02055100$flowdata$newwave=
  USGS02055100$flowdata$cms *1.1 <
  data.table::shift(USGS02055100$flowdata$cms)
summary(USGS02055100$flowdata$newwave)
# Add plot of the point found
len=length(USGS02055100$flowdata$newwave)
USGS02055100$flowdata$newwave[is.na(USGS02055100$flowdata$newwave)]=F
# Removes repeated finds by going through loop backwords
for (i in seq(len,2)){
  print(i)
  if(USGS02055100$flowdata$newwave[i]==T &
     USGS02055100$flowdata$newwave[i-1]==T){
    USGS02055100$flowdata$newwave[i]=F
  }
}
plot(USGS02055100$flowdata$dateTime,USGS02055100$flowdata$cms,type="l")
points(USGS02055100$flowdata$dateTime[USGS02055100$flowdata$newwave],
       USGS02055100$flowdata$cms[USGS02055100$flowdata$newwave],col=2)

# Find the time locations where waves begin
which(USGS02055100$flowdata$newwave == TRUE)
plot(USGS02055100$flowdata$dateTime,USGS02055100$flowdata$cms,
     type="l",xlim=c(USGS02055100$flowdata$dateTime[1109],
                     USGS02055100$flowdata$dateTime[1109+200]), xlab="time", ylab="flow (m3/s)")
lines(USGS02055100$flowdata$outTime,USGS02055100$flowdata$cms,col=2)

################################################### USGS02054530
# USGS gage height for USGS02056000
# Distance between points.
A=SpatialPoints(USGS02054530$site)# Up gradient site
B=SpatialPoints(USGS02056000$site) # Down gradient site ROA River atNiagara
proj4string(A)=proj4_ll
proj4string(B)=proj4_ll
A_utm=spTransform(A,crs_utm)
B_utm=spTransform(B,crs_utm)
# Cut the DEM down to a more manageable size
cropmydem=crop(mydem,extend(extent(ab_utm),600))
cropmydem=trim(cropmydem)
cropmydem=cropmydem*1000.0
# Set up the weighting functions
altDiff <- function(x){x[2] - x[1]}
hd <- transition(cropmydem, altDiff, 8, symm=FALSE)
slope <- geoCorrection(hd)
adj <- adjacent(cropmydem, cells=1:ncell(cropmydem), pairs=TRUE, directions=8)
speed <- slope
speed[adj] <- 6 * exp(-3.5 * abs(slope[adj] + 0.05))
Conductance <- geoCorrection(speed)
# Find and plot the flow path
AtoB <- shortestPath(Conductance, A_utm, B_utm, output="SpatialLines")
plot(AtoB,add=T)
plot(streams_utm,col="blue",add=T)
plot(AtoB,add=T)
SpatialLinesLengths(AtoB)
USGS02054530$site$L=SpatialLinesLengths(AtoB) # km to m
USGS02054530$site$L # reach length in m

# Slope of stream
USGS02054530$site$slope=(extract(mydem,A_utm)-
                           extract(mydem,B_utm))/USGS02054530$site$L
USGS02054530$site$slope

USGS02054530$flowdata$B=(USGS02054530$site$man_n*
                           USGS02054530$flowdata$cms)/(USGS02054530$flowdata$depth_m^(5/3)*
                                                         sqrt(USGS02054530$site$slope))
head(USGS02054530$flowdata)
# Lets look at how B changes with flow.    
plot(USGS02054530$flowdata$dateTime,USGS02054530$flowdata$B, main="LICK RUN TO ROANOKE RIVER AT NIAGARA, VA")
# Does this seem reasonable (...like order of magnitude reasonable)? You can 
# perform a quick and dirty check using google earth and measuring the channel 
# width in a few places.
#
plot(USGS02054530$flowdata$cms,USGS02054530$flowdata$depth_m, main="LICK RUN TO ROANOKE RIVER AT NIAGARA, VA")


# ck
# USGS02054530$flowdata$ck = ???
# ANS
USGS02054530$flowdata$ck =
  5/3*sqrt(USGS02054530$site$slope)/USGS02054530$site$man_n*
  (USGS02054530$flowdata$depth_m^(2/3))
#
# USGS02054530$flowdata$dt = ???
USGS02054530$flowdata$dt =
  USGS02054530$site$L/USGS02054530$flowdata$ck

plot(USGS02054530$flowdata$dateTime,USGS02054530$flowdata$dt)
USGS02054530$flowdata$outTime=USGS02054530$flowdata$dateTime+
  USGS02054530$flowdata$dt

# Find beginning of  Waves
USGS02054530$flowdata$newwave=
  USGS02054530$flowdata$cms *1.1 <
  data.table::shift(USGS02054530$flowdata$cms)
summary(USGS02054530$flowdata$newwave)
# Add plot of the point found
len=length(USGS02054530$flowdata$newwave)
USGS02054530$flowdata$newwave[is.na(USGS02054530$flowdata$newwave)]=F
# Removes repeated finds by going through loop backwords
for (i in seq(len,2)){
  print(i)
  if(USGS02054530$flowdata$newwave[i]==T &
     USGS02054530$flowdata$newwave[i-1]==T){
    USGS02054530$flowdata$newwave[i]=F
  }
}
plot(USGS02054530$flowdata$dateTime,USGS02054530$flowdata$cms,type="l")
points(USGS02054530$flowdata$dateTime[USGS02054530$flowdata$newwave],
       USGS02054530$flowdata$cms[USGS02054530$flowdata$newwave],col=2)

# Find the time locations where waves begin
which(USGS02054530$flowdata$newwave == TRUE)
plot(USGS02054530$flowdata$dateTime,USGS02054530$flowdata$cms,
     type="l",xlim=c(USGS02054530$flowdata$dateTime[1109],
                     USGS02054530$flowdata$dateTime[1109+200]), xlab="time", ylab="flow (m3/s)")
lines(USGS02054530$flowdata$outTime,USGS02054530$flowdata$cms,col=2)


################################################### USGS02055000
# USGS gage height for USGS02056000
# Distance between points.
A=SpatialPoints(USGS02055000$site)# Up gradient site
B=SpatialPoints(USGS02056000$site) # Down gradient site ROA River atNiagara
proj4string(A)=proj4_ll
proj4string(B)=proj4_ll
A_utm=spTransform(A,crs_utm)
B_utm=spTransform(B,crs_utm)
# Cut the DEM down to a more manageable size
cropmydem=crop(mydem,extend(extent(ab_utm),600))
cropmydem=trim(cropmydem)
cropmydem=cropmydem*1000.0
# Set up the weighting functions
altDiff <- function(x){x[2] - x[1]}
hd <- transition(cropmydem, altDiff, 8, symm=FALSE)
slope <- geoCorrection(hd)
adj <- adjacent(cropmydem, cells=1:ncell(cropmydem), pairs=TRUE, directions=8)
speed <- slope
speed[adj] <- 6 * exp(-3.5 * abs(slope[adj] + 0.05))
Conductance <- geoCorrection(speed)
# Find and plot the flow path
AtoB <- shortestPath(Conductance, A_utm, B_utm, output="SpatialLines")
plot(AtoB,add=T)
plot(streams_utm,col="blue",add=T)
plot(AtoB,add=T)
SpatialLinesLengths(AtoB)
USGS02055000$site$L=SpatialLinesLengths(AtoB) # km to m
USGS02055000$site$L # reach length in m

# Slope of stream
USGS02055000$site$slope=(extract(mydem,A_utm)-
                           extract(mydem,B_utm))/USGS02055000$site$L
USGS02055000$site$slope

USGS02055000$flowdata$B=(USGS02055000$site$man_n*
                           USGS02055000$flowdata$cms)/(USGS02055000$flowdata$depth_m^(5/3)*
                                                         sqrt(USGS02055000$site$slope))
head(USGS02055000$flowdata)
# Lets look at how B changes with flow.    
plot(USGS02055000$flowdata$dateTime,USGS02055000$flowdata$B, main="LICK RUN TO ROANOKE RIVER AT NIAGARA, VA")
# Does this seem reasonable (...like order of magnitude reasonable)? You can 
# perform a quick and dirty check using google earth and measuring the channel 
# width in a few places.
#
plot(USGS02055000$flowdata$cms,USGS02055000$flowdata$depth_m, main="LICK RUN TO ROANOKE RIVER AT NIAGARA, VA")


# ck
# USGS02055000$flowdata$ck = ???
# ANS
USGS02055000$flowdata$ck =
  5/3*sqrt(USGS02055000$site$slope)/USGS02055000$site$man_n*
  (USGS02055000$flowdata$depth_m^(2/3))
#
# USGS02055000$flowdata$dt = ???
USGS02055000$flowdata$dt =
  USGS02055000$site$L/USGS02055000$flowdata$ck

plot(USGS02055000$flowdata$dateTime,USGS02055000$flowdata$dt, type="l")
USGS02055000$flowdata$outTime=USGS02055000$flowdata$dateTime+
  USGS02055000$flowdata$dt

# Find beginning of  Waves
USGS02055000$flowdata$newwave=
  USGS02055000$flowdata$cms *1.1 <
  data.table::shift(USGS02055000$flowdata$cms)
summary(USGS02055000$flowdata$newwave)
# Add plot of the point found
len=length(USGS02055000$flowdata$newwave)
USGS02055000$flowdata$newwave[is.na(USGS02055000$flowdata$newwave)]=F
# Removes repeated finds by going through loop backwords
for (i in seq(len,2)){
  print(i)
  if(USGS02055000$flowdata$newwave[i]==T &
     USGS02055000$flowdata$newwave[i-1]==T){
    USGS02055000$flowdata$newwave[i]=F
  }
}
plot(USGS02055000$flowdata$dateTime,USGS02055000$flowdata$cms,type="l")
points(USGS02055000$flowdata$dateTime[USGS02055000$flowdata$newwave],
       USGS02055000$flowdata$cms[USGS02055000$flowdata$newwave],col=2)

# Find the time locations where waves begin
which(USGS02055000$flowdata$newwave == TRUE)
plot(USGS02055000$flowdata$dateTime,USGS02055000$flowdata$cms,
     type="l",xlim=c(USGS02055000$flowdata$dateTime[1109],
                     USGS02055000$flowdata$dateTime[1109+200]), xlab="time", ylab="flow (m3/s)")
lines(USGS02055000$flowdata$outTime,USGS02055000$flowdata$cms,col=2)















