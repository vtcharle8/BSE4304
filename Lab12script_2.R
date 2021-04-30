objects()
rm(list=objects())

if (!require("pacman")) install.packages("pacman")
pacman::p_load(EcoHydRology,curl,httr,rnoaa,DEoptim,rstudioapi)
setwd("~/Week12/")
vignette("DEoptim")
?DEoptim

url="https://r-forge.r-project.org/scm/viewvc.php/*checkout*/pkg/EcoHydRology/R/FillMissWX.R?root=ecohydrology"
source(url)
url="https://r-forge.r-project.org/scm/viewvc.php/*checkout*/pkg/EcoHydRology/man/FillMissWX.Rd?root=ecohydrology"
download.file(url,"FillMissWX.Rd")
rstudioapi::previewRd("FillMissWX.Rd")

# Weather data
# Using the latitude and longitude of USGS 04282650 
# LITTLE OTTER CREEK AT FERRISBURG, VT.
## Not run: 
flowgage_id="04282650" 
flowgage=get_usgs_gage(flowgage_id,begin_date = "2010-01-01",
                       end_date = "2022-01-01")
LOCFer=FillMissWX(declat=flowgage$declat, declon=flowgage$declon,
                  StnRadius=20,date_min="2010-01-01",date_max="2021-01-01")
plot(LOCFer$date,LOCFer$prcpDis, xlab="Date", 
     ylab = "Distance (km)",ylim=c(0,30))
points(LOCFer$date,LOCFer$tmaxDis-0.5,col="red",pch=3,cex=.3)
points(LOCFer$date,LOCFer$tminDis+.5,col="blue",pch=4,cex=.3)
legend("topleft",legend = c("prcpDis","tmaxDis-0.5","tminDis+0.5"),
       col = c("black","red","blue"),lty = 1:2, cex = 0.7)

save(LOCFer, file="LOCFer.Rd")
url="https://github.com/vtdrfuka/BSE5304/raw/main/LOCFer.Rd"
download.file(url, "LOCFer2.Rd")
load("LOCFer2.Rd")
download.file("https://github.com/vtdrfuka/BSE5304/raw/main/functions.R", "functions.R")
source("functions.R")


flowgage$flowdata$Qmm = flowgage$flowdata$flow/flowgage$area/10^3
modeldata=merge(LOCFer,flowgage$flowdata,by.x="date",by.y="mdate")
#
# Remember our Qmm vs dP plots (Flow vs Delta Precipitation)?
# Run either CN or TMWB Model to get our dP. 
#
modeldata$HillslopeAboveExcess=0
Qmm_dP=CN_Model(modeldata)
#
# Simple Calibration based on our earliest CN estimation model
# from Lab06(ish)
#
View(Qmm_dP)
Qmm=Qmm_dP$Qmm
dP=Qmm_dP$dP
plot(dP,Qmm)
points(dP,dP^2/(dP+45),col="red")  # S guestimates in bold
points(dP,dP^2/(dP+260),col="blue")# S guestimates in bold
#
# Build our function to optimize NSE based on Qmm and dP
#
Sestfun <- function(x){
  x <- x[1]
  NSE1=NSeff(Qmm,dP^2/(dP+x))
  return(abs(NSE1-1))
}
# Test and explore our estimation of S (Sest)
Sestfun(45)
Sestfun(260)



lower <- c(30)
upper <- c(700)
# run DEoptim and set a seed first for replicability
set.seed(1234)
DEoptim(Sestfun, lower, upper)
outDEoptim=DEoptim(Sestfun, lower, upper)
outDEoptim$optim$bestmem[1]
outDEoptim$optim$bestval
# NSE = 1- bestval
1-outDEoptim$optim$bestval
BestS=outDEoptim$optim$bestmem[1]
CN=(25.4*1000/BestS)+10
CN
rm(list=c("dP","Qmm"))


CN_ModelFun <- function(x){
  #x=outDEoptim$optim$bestmem   # This will be useful for building Flow vs. Date
  CNavg <- x[1]
  IaFrac = x[2]
  fnc_slope=x[3] 
  fnc_aspect=x[4]
  func_DAWC=x[5]
  func_z=x[6]
  fnc_fcres=x[7]
  TempBias=x[8]
  Qmm_dP=CN_Model(modeldata, 
                  CNavg = CNavg,IaFrac = IaFrac,fnc_slope=fnc_slope, 
                  fnc_aspect=fnc_aspect,func_DAWC=func_DAWC,func_z=func_z,
                  fnc_fcres=fnc_fcres,TempBias=TempBias,declat=flowgage$declat)
  NSE1=NSeff(Qmm_dP$Qmm,Qmm_dP$Qpred)
  return(abs(NSE1-1))
}


lower <- c(30,0.01,0.0,0.0,0.2, 500,0.1,-3.0)
upper <- c(99,0.20,0.1,0.1,0.4,3000,0.5, 3.0)
# run DEoptim and set a seed first for replicability
set.seed(1234)

if (!require("pacman")) install.packages("pacman")
pacman::p_load(EcoHydRology,curl,httr,rnoaa,DEoptim,rstudioapi)
setwd("~/Week12/")
load("~/Week12/lab12env.RData")

system.time((outDEoptim=DEoptim(CN_ModelFun, lower, upper,
                                DEoptim.control(strategy = 6,NP = 16,itermax=100,parallelType = 1,
                                                packages = c("EcoHydRology"),parVar=c("CN_Model","modeldata","flowgage")))))

outDEoptim$optim$bestmem
# NSE = 1- bestval
1-outDEoptim$optim$bestval
CN_ModelFun(outDEoptim$optim$bestmem)
names(outDEoptim$member$lower)=c("CNavg","IaFrac","fnc_slope", 
                                 "fnc_aspect","func_DAWC","func_z",
                                 "fnc_fcres","TempBias")
plot(outDEoptim)
plot(outDEoptim,plot.type="bestvalit")
plot(outDEoptim,plot.type="bestmemit")
dev.off()
plot(1-outDEoptim$member$bestvalit,ylab="NSE")
dev.off()

# HW 2 plot for CN range [30,99]
dev.off()
plot(Qmm_dP$date, Qmm_dP$Qpred, type="l", col="blue", xlab="Date", ylab="Flow")
lines(Qmm_dP$date, Qmm_dP$Qmm, type="l", col="red")
