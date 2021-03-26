# Cleaning up
objects()
rm(list=objects())
#
# Build a working directory for this weeks lab and change working dir
setwd("~/Week08")
url="https://raw.githubusercontent.com/vtdrfuka/BSE5304/main/Lab07Sol.R"
download.file(url,"Lab07Sol.R")
file.edit("Lab07Sol.R")

# Start of Lab 8
# Make a poly with raster library (slow)
# or from thee command line gdal (fast)
# gdal_polygonize.py -8 mydemw.tif mydemw_poly_gdal.shp   <<< run this in terminal.
# mydemw_poly=rasterToPolygons(mydemw,dissolve = T,na.rm = T)  <<< takes too long, use ^^^.
mydemw_poly=readOGR("mydemw_poly_gdal.shp")
plot(mybasindem)
plot(mydemw_poly,add=T,border="blue")
mydemw_poly
writeOGR(mydemw_poly,dsn=".",layer="mydemw",driver="ESRI Shapefile", overwrite_layer=TRUE)
# We will use this ESRI shape file, a zipped group of it, to download 
# our soil extent from the WebSoilSurvey Website
zip("mydemw.zip",list.files(pattern="mydemw[:.:]"))
# Download to your local machine mydemw.zip from the "Files" tab
# Open the WebSoilSurvey site to: 
browseURL("https://websoilsurvey.sc.egov.usda.gov/App/WebSoilSurvey.aspx")
# "Creat AOI from a zipped shapefile"
# Open "Download Soils Data" Tab
# "Create Download Link" in lower right hand corner
# Right-Click on download link and "Copy Link Address" and 
# paste into a url object
url="https://websoilsurvey.sc.egov.usda.gov/DSD/Download/AOI/zyt2p4ubwccybyhhko301equ/wss_aoi_2021-03-25_18-25-12.zip"
download.file(url,"wss_aoi.zip")
unzip("wss_aoi.zip")
mysoil_ll=readOGR("wss_aoi_2021-03-25_18-25-12/spatial/soilmu_a_aoi.shp")
View(mysoil_ll@data)

# Notice that MUKEY is all capital but, in lab 4 everything was lower case
mysoil_ll@data$mukey=mysoil_ll@data$MUKEY # put mukey in lower case in new column heading
View(mysoil_ll@data)

mysoil_utm = spTransform(mysoil_ll,crs_utm)
plot(mysoil_utm,add=T)
#
# Note that in the rasterize command, the TIC raster layer is being used as 
# a reference raster. Any of the matching UTM raster layers could be used in 
# its place. 
rmysoil_utm=rasterize(mysoil_utm,TIC,field=as.numeric(mysoil_utm$mukey))
#
ratify(TIC)
ratify(round(mybasinslp*10))
ratify(TIC*10^9)

unique(ratify(round(mybasinslp*10+1)))
unique(ratify((rmysoil_utm*10^3)))
#
# Now build an HRU table with the combination of the 1) raster Soils, 2) TIC,
# and 3) slope layers. 
#
hru=ratify(TIC*10^9 + (rmysoil_utm*10^3) + round(mybasinslp*10+1))
unique(values(hru))
length(unique(values(hru)))
pacman::p_load(circlize)
plot(hru,col=rand_color(length(unique(values(hru)))))

# Build an HRU attribute table
hru_table = levels(hru)[[1]]
origID = hru_table$ID # preserve data order for later reordering
# metadata parameters from a string... this will make more sense
# after the next "head()" command
hru_table$TIclass = as.numeric(substr(sprintf("%10.0f", hru_table$ID), 1,1))
hru_table$mukey = as.numeric(substr(sprintf("%10.0f", hru_table$ID), 2,7))
hru_table$slp = (as.numeric(substr(sprintf("%10.0f", 
                                             hru_table$ID), 8,10))-1)/10
#
# Calculate the area for each unique soil (mukey) X TIClass combination
# using res(raster) for x and y resolution in units of m
hru_table$areaSQKM = as.vector(round(res(hru)[1]*res(hru)[2]*
                                         table(values(hru))/10^6, 3))
#
# To better understand what happened, look at the new hru_table
head(hru_table)
summary(hru_table)

rm("mu2co")
for (mymukey in unique(hru_table$mukey)){
  print(mymukey)
  mukey_statement = format_SQL_in_statement(mymukey)
  q_mu2co = paste("SELECT mukey,cokey FROM component 
           WHERE mukey IN ", mukey_statement, sep="")
  if(!exists("mu2co")){
    mu2co=SDA_query(q_mu2co)} 
  else{
    mu2co=rbind(mu2co,SDA_query(q_mu2co))
  } 
}
View(mu2co)
# Second associate cokey with ksat_r,awc_r,hzdepb_r from chorizon
# cokey_statement = format_SQL_in_statement(unique(mu2co$cokey))
# q_co2ch = paste("SELECT cokey,ksat_r,awc_r,hzdepb_r,frag3to10_r  
# FROM chorizon WHERE cokey IN ", cokey_statement, sep="")
# co2ch = SDA_query(q_co2ch)
rm("co2ch")
for (mycokey in unique(mu2co$cokey)){
  print(mycokey)
  cokey_statement = format_SQL_in_statement(mycokey)
  q_co2ch = paste("SELECT cokey,ksat_r,awc_r,hzdepb_r,frag3to10_r FROM chorizon WHERE cokey IN ", cokey_statement, sep="")
  print(q_co2ch)
    if(!exists("co2ch")){
    co2ch=SDA_query(q_co2ch)
  } else{
    try((co2ch=rbind(co2ch,SDA_query(q_co2ch))))
  } 
}
View(co2ch)
rm("co2co")
for (mycokey in unique(mu2co$cokey)){
  print(mycokey)
  cokey_statement = format_SQL_in_statement(mycokey)
  q_co2co = paste("SELECT cokey,slopelenusle_r FROM component WHERE cokey IN ", cokey_statement, sep="")
  print(q_co2co)
  if(!exists("co2co")){
    co2co=SDA_query(q_co2co)} 
  else{
    try((co2co=rbind(co2co,SDA_query(q_co2co))))} 
}
View(co2co)
# Last, bring them back together, and aggregate based on max values
# of ksat_r,awc_r, and hzdepb_r
mu2ch=merge(mu2co,co2ch)
mu2ch=merge(mu2ch,co2co)
View(mu2ch)

# Merge and then aggregate spatially derived hru_table with
# SSURGO summary table
 MUSLE_mrg=merge(hru_table,mu2ch)   # Grad Homework 4 base
 MUSLE_mrg$ksat_r=as.numeric(MUSLE_mrg$ksat_r)
 MUSLE_mrg$awc_r=as.numeric(MUSLE_mrg$awc_r)
 MUSLE_mrg$hzdepb_r=as.numeric(MUSLE_mrg$hzdepb_r)
 MUSLE_mrg$slopelenusle_r=as.numeric(MUSLE_mrg$slopelenusle_r)
 MUSLE_mrg$frag3to10_r=as.numeric(MUSLE_mrg$frag3to10_r)
 MUSLE=aggregate(MUSLE_mrg,list(MUSLE_mrg$TIclass),mean,na.rm=T)
#
# Easiest first! Eq. 4:1.1.15 Course Fragment Factor
 MUSLE$CFRG=exp(-0.053*MUSLE$frag3to10_r)
 MUSLE
#
# LSusle is calculated using eq. 4.1.12
 MUSLE$alpha=atan(MUSLE$slp/100)
 MUSLE$LSm=.6*(1-exp(-35.835*MUSLE$slp/100))
 MUSLE$LS=(MUSLE$slopelenusle_r/22.1)^MUSLE$LSm * (65.41*sin(MUSLE$alpha)^2+4.56*sin(MUSLE$alpha)+0.065)
#
# Pusle
 MUSLE$Pusle=.50
#
# Cusle
 MUSLE$Cusle=.20
#
# Kusle
 MUSLE$Kusle=0.28
#
# Build a constant for those we are not changing day to day
 attach(MUSLE)
 MUSLE$KCPLSCFRG118=11.8*Kusle*Cusle*Pusle*LS*CFRG
 detach(MUSLE)
 MUSLE # Make sure values look correct, Pusle, Cusle, Kusle
#
# Now we need to use each of the TIClass Q solutions from Lab06 to calculate
# peak flows (qpeak) and complete the MUSLE Sediment Loss for each class.
# Run Model
#
# Now we need to use each of the TIClass Q solutions from Lab07 to calculate
# peak flows (qpeak) and complete the MUSLE Sediment Loss for each class.
# Run Model
download.file("https://github.com/vtdrfuka/BSE5304/raw/main/grabdata.R",
                "grabdata.R")
download.file("https://github.com/vtdrfuka/BSE5304/raw/main/functions.R",
              "functions.R")
 pacman::p_load(data.table)
# Use the estimated S for our watershed (Lab06)
 Sest = 157
# We will split into 5 VSA areas represented by 5 TI Classes
 nTIclass=5
 VSAsol=data.table(TIClass=seq(from=nTIclass,to=1),
                    As=seq(1:nTIclass)*(1/nTIclass),Wetfrac=(1/nTIclass))
 VSAsol[,sSratio:=2*(sqrt(1-shift(As))-sqrt(1-As))/Wetfrac-1]
#
 VSAsol$sSratio[1]=2*(sqrt(1-0)-sqrt(1-VSAsol$As[1]))/VSAsol$Wetfrac[1]-1
# Calculate TI Class localized sigma and Curve Number
 VSAsol[,sigma:=Sest*sSratio]
 VSAsol[,CN:=25400/(sigma+254)]
 VSAsol

 VSAParams=merge(VSAsol,MUSLE,by.x="TIClass",by.y="TIclass")
 View(VSAParams)
 modeldata$HillslopeAboveExcess=0
 TIC01=modeldata
 TIC02=modeldata
 TIC03=modeldata
 TIC04=modeldata
 TIC05=modeldata
 
# For TIC01 CNavg=VSAParams$CN[1] but confirm
 TIC01 = CN_Model(fnc_CNModel = TIC01, CNavg=VSAParams$CN[1])
 TIC01$qpeak=TIC01$Qpred/3600/24/1000*myflowgage$area/nTIclass*10^6 #m^3/sec
 TIC01$sed=(TIC01$Qpred*TIC01$qpeak*myflowgage$area/nTIclass*100)^.56*MUSLE$KCPLSCFRG118[1]    # Eq. 4:1.1.1 SWAT Theory
# Route the water and continue
 TIC02$HillslopeAboveExcess=TIC01$Qpred
 TIC02 = CN_Model(fnc_CNModel = TIC02, CNavg=VSAParams$CN[2])
 TIC02$qpeak=TIC02$Qpred/3600/24/1000*myflowgage$area/nTIclass*10^6
 TIC02$sed=(TIC02$Qpred*TIC02$qpeak*myflowgage$area/nTIclass*100)^.56*MUSLE$KCPLSCFRG118[2]
# TIC03
 TIC03$HillslopeAboveExcess=TIC02$Qpred
 TIC03 = CN_Model(fnc_CNModel = TIC03, CNavg=VSAParams$CN[3])
 TIC03$qpeak=TIC03$Qpred/3600/24/1000*myflowgage$area/nTIclass*10^6
 TIC03$sed=(TIC03$Qpred*TIC03$qpeak*myflowgage$area/nTIclass*100)^.56*MUSLE$KCPLSCFRG118[3]
# TIC04
 TIC04$HillslopeAboveExcess=TIC03$Qpred
 TIC04 = CN_Model(fnc_CNModel = TIC04, CNavg=VSAParams$CN[4])
 TIC04$qpeak=TIC04$Qpred/3600/24/1000*myflowgage$area/nTIclass*10^6
 TIC04$sed=(TIC04$Qpred*TIC04$qpeak*myflowgage$area/nTIclass*100)^.56*MUSLE$KCPLSCFRG118[4]
# TIC05
 TIC05$HillslopeAboveExcess=TIC04$Qpred
 TIC05 = CN_Model(fnc_CNModel = TIC05, CNavg=VSAParams$CN[5])
 TIC05$qpeak=TIC05$Qpred/3600/24/1000*myflowgage$area/nTIclass*10^6
 TIC05$sed=(TIC05$Qpred*TIC05$qpeak*myflowgage$area/nTIclass*100)^.56*MUSLE$KCPLSCFRG118[5]

 # Plot sediment over time for each TI
 plot(TIC01$date, TIC01$sed)
 plot(TIC02$date, TIC02$sed)
 plot(TIC03$date, TIC03$sed)
 plot(TIC04$date, TIC04$sed)
 plot(TIC05$date, TIC05$sed)
 