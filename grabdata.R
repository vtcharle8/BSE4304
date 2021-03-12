# Save your functions and build up your working "modeldata" dataframe
# so you can save a package.
objects()
rm(list=objects())
dir.create("~/Week06")
setwd("~/Week06")
#
# We will add to the mix a handy date manipulation package
# which you will use the month() function to isolate non-snow months
# below
#

if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table,httr,EcoHydRology,curl,elevatr,raster,soilDB,
               rgdal,lubridate, BSEHydroModels)



myflowgage_id="0205551460"
myflowgage=get_usgs_gage(myflowgage_id,begin_date = "2015-01-01",end_date = "2021-03-01")

# Note that flow returned is in m3/day, but we want mm/day for the basin
myflowgage$flowdata$Qmm = myflowgage$flowdata$flow/myflowgage$area/10^3
stns=meteo_distance(
  station_data=ghcnd_stations(),
  lat=myflowgage$declat,
  long=myflowgage$declon,
  units = "deg",
  radius = 30,
  limit = NULL
)
# We are looking for stations with elements that have PRCP, TMAX and TMIN 
# and current data (i.e. Year 2021). 
WXStn=stns[stns$element=="TMAX"&stns$last_year>=2020,]$id[2]
WXData=meteo_pull_monitors(
  monitors=WXStn,
  keep_flags = FALSE,
  date_min = "2016-01-01",
  date_max = NULL,
  var = c("TMAX","TMIN","PRCP") 
)
summary(WXData)
# Create an aligned modeldata data frame to build our model in
modeldata=merge(WXData,myflowgage$flowdata,by.x="date",by.y="mdate")
modeldata$MaxTemp=modeldata$tmax/10 # Converting to C
modeldata$MinTemp=modeldata$tmin/10 # Converting to C
modeldata$P=modeldata$prcp/10 # Converting to mm
# Compare your precipitation to the flow out of your basin
mean(modeldata$Qmm)
mean(modeldata$P)
modeldata$P[is.na(modeldata$P)]=0
modeldata$MinTemp[is.na(modeldata$MinTemp)]=0
modeldata$MaxTemp[is.na(modeldata$MaxTemp)]=
  modeldata$MinTemp[is.na(modeldata$MaxTemp)] +1
modeldata$MaxTemp[modeldata$MaxTemp<=modeldata$MinTemp]=
  modeldata$MinTemp[modeldata$MaxTemp<=modeldata$MinTemp]+1
modeldata$AvgTemp=(modeldata$MaxTemp+modeldata$MaxTemp)/2.0



#
# We will explore
#
url="https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/NHD/HU8/HighResolution/Shape/NHD_H_03010101_HU8_Shape.zip"
curl_download(url,"NHD_H_03010101_HU8_Shape.zip")
unzip("NHD_H_03010101_HU8_Shape.zip",exdir="03010101")
streams=readOGR("03010101/Shape/NHDFlowline.dbf")
mystream=subset(streams,GNIS_ID=="01478950")
plot(mystream,col="red")
#
# Use the spatial extents from our stream to download elevation raster.
#
proj4_ll = "+proj=longlat"
proj4string(mystream) = proj4_ll
mydem=get_aws_terrain(locations=mystream, 
                      z = 11, prj =proj4string(mystream) ,src ="aws",clip="bbox")
#
# Pretty pictures of our area help ease the frustration
#
plot(mydem)
lines(mystream,col="blue",lwd=4)
points(myflowgage$declon,myflowgage$declat,pch = 24, cex=2, col="blue", bg="red", lwd=2)


#
# For initializing slopes, we store the summary stats for terrain slope
#
plot(terrain(mydem, opt='slope',unit = "radians"))
lines(mystream,col="blue",lwd=4)
points(myflowgage$declon,myflowgage$declat,pch = 24, cex=2, col="blue", bg="red", lwd=2)
slope_sum=summary(terrain(mydem, opt='slope',unit = "radians"))
#
# And after our initialization is done, we are ready to get our estimated 
# dP from our TMWB model. Remember that dP = P - ET - SnowFall + SnowMelt
# 
# We are building on our prior lab solutions, need to build out our previous 
# TMWB model. We will grab functions from the solutions from Week 4’s Lab  
# Lab05

modeldata$HillslopeAboveExcess=0
BasinTMWB=modeldata
BasinTMWB = TMWB_Model(fnc_TMWB = BasinTMWB,fnc_slope=0, 
                       fnc_aspect=0,func_DAWC=.3,
                       func_z=500,fnc_fcres=.3)


attach(BasinTMWB)
plot(date,AW)
plot(dP,Qmm)
detach(BasinTMWB)
#
# But, we know that our systems behave differently during snowy winter
# months, so we will isolate our June ( month>5) - October ( < 11 ) data (_JO)
#
BasinTMWB_JO=BasinTMWB[(month(BasinTMWB$date) > 5 
                        & month(BasinTMWB$date) < 11),]
attach(BasinTMWB_JO)
plot(dP,Qmm)
detach(BasinTMWB_JO)

attach(BasinTMWB_JO)
NSE=function(Yobs,Ysim){
  return(1-sum((Yobs-Ysim)^2, na.rm=TRUE)/sum((Yobs-mean(Yobs, na.rm=TRUE))^2, na.rm=TRUE))
}

# Loop through different S values to optimize NSE.
for (myS in 50:350){
  print(paste0(myS, NSE(Qmm,dP^2/(dP+myS))))
}

# Estimate S=157 from the for-loop.
Sest=157
plot(dP,Qmm)
points(dP,dP^2/(dP+Sest),col="red")
detach(BasinTMWB_JO)


# We will split into 5 VSA areas represented by 5 TI Classes
nTIclass=5
VSAsol=data.table(WetClass=seq(from=nTIclass,to=1),
                  As=seq(1:nTIclass)*(1/nTIclass),Wetfrac=(1/nTIclass))
VSAsol[,sSratio:=2*(sqrt(1-shift(As))-sqrt(1-As))/Wetfrac-1]

#Fill in first value
VSAsol$sSratio[1]=2*(sqrt(1-0)-sqrt(1-VSAsol$As[1]))/VSAsol$Wetfrac[1]-1

# Calculate TI Class localized sigma and Curve Number
VSAsol[,sigma:=Sest*sSratio]
VSAsol[,CN:=25400/(sigma+254)]
VSAsol
plot(VSAsol$As,VSAsol$sigma)
lines(VSAsol$As,VSAsol$sigma)
plot(VSAsol$As,VSAsol$CN)
lines(VSAsol$As,VSAsol$CN)

# Use the CN_Model Function from last week, but set up 5 TI Classes instead
# of 3 hillslope locations
#
#
# Initialize the TI Class objects from top to bottom of slope
TIC05=modeldata
# For TIC05 CNavg=VSAsol$CN[1]
TIC05 = CNModel(fnc_CNModel = TIC05, CNavg=VSAsol$CN[1],
                 func_DAWC=.3,IaFrac=0.05,
                 func_z=1000,fnc_fcres=.3)
TIC04=modeldata
TIC04$HillslopeAboveExcess=TIC05$Qpred
TIC04 = CNModel(fnc_CNModel = TIC04, CNavg=VSAsol$CN[2],
                func_DAWC=.3,IaFrac=0.05,
                func_z=1000,fnc_fcres=.3)
TIC03=modeldata
TIC03$HillslopeAboveExcess=TIC04$Qpred
TIC03 = CNModel(fnc_CNModel = TIC03, CNavg=VSAsol$CN[3],
                func_DAWC=.3,IaFrac=0.05,
                func_z=1000,fnc_fcres=.3)
TIC02=modeldata
TIC02$HillslopeAboveExcess=TIC03$Qpred
TIC02 = CNModel(fnc_CNModel = TIC02, CNavg=VSAsol$CN[4],
                func_DAWC=.3,IaFrac=0.05,
                func_z=1000,fnc_fcres=.3)
TIC01=modeldata
TIC01$HillslopeAboveExcess=TIC02$Qpred
TIC01 = CNModel(fnc_CNModel = TIC01, CNavg=VSAsol$CN[5],
                func_DAWC=.3,IaFrac=0.05,
                func_z=1000,fnc_fcres=.3)

# HW1 plots of runoff.
plot(TIC05$date,TIC05$Qpred,type="l")
plot(TIC04$date,TIC04$Qpred,type="l")
plot(TIC03$date,TIC03$Qpred,type="l")
plot(TIC02$date,TIC02$Qpred,type="l")
plot(TIC01$date,TIC01$Qpred,type="l")

# HW2 plots of AW.
plot(TIC05$date, TIC05$AW, type="l")
plot(TIC04$date, TIC04$AW, type="l")
plot(TIC03$date, TIC03$AW, type="l")
plot(TIC02$date, TIC02$AW, type="l")
plot(TIC01$date, TIC01$AW, type="l")

# HW2 plots of ET.
plot(TIC05$date, TIC05$ET, type="l")
plot(TIC04$date, TIC04$ET, type="l")
plot(TIC03$date, TIC03$ET, type="l")
plot(TIC02$date, TIC02$ET, type="l")
plot(TIC01$date, TIC01$ET, type="l")
