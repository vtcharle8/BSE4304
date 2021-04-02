# Note we have the solutions from last weeks lab in objects TIC05,
# TIC04, etc. 
ls(pattern = "TIC*")

# 
# Build a dataframe for your DP Load model in TI Class 05,(you will need to do 
#this 5 times for TIC 1-4 where you will need the date, Qpred, and Average 
# Temperature
DPTI05=data.frame(date=TIC05$date,
                  Rt=TIC05$Qpred/1000, 
                  Tavg=(TIC05$MaxTemp+TIC05$MinTemp)/2)
# Take it to volumetric water content rather # than volumetric soil content, 
# should be in m of water
tau=9.3  # days
dt=1     # days time step
kF=.015  # Table 2
# Initialize MF and DF
DPTI05$MF=0
DPTI05$DF=0
# Spread your P Fertilizer on ~May 1, ~August 1, and ~October 1
DPTI05$MF[(format(DPTI05$date,"%j") %in% c(121,213,274))]=5.4*10^-4
# Remember what we say about attaching! 
attach(DPTI05)
#
# Loop to solve MF and DF
for (i in 2:length(date)){
  if(MF[i]<=MF[i-1]){
    MF[i]=MF[i-1]*exp(-dt/tau)-DF[i-1]
  }
  DF[i]=MF[i]*(kF*MF[i]*Rt[i]/(1+kF*MF[i]*Rt[i]))
}
DPTI05$MF=MF
DPTI05$DF=DF
detach(DPTI05)
rm(list=c("MF","DF")) # Clean up the environment 
dev.off() #reset graphics device
plot(DPTI05$date,DPTI05$MF)
plot(DPTI05$date,DPTI05$DF)

# Calculate your Export Coef from Easton et al 2007 Figure 2 using TI Class 5
# and note that the bold gives the TI Class. This figure gives a range of 
# 0 - 520 micrograms/litre 
# For TIC=5
muTS_TI05=(((520-0)/5)*VSAsol$TIClass[1]+0) # use VSAsol$TIClass table to 
# assign TIC to calc (remember TIC05 is in location 1 in the VSAsol$TIClass 
# table
# Setting range of Soil P values (MS) using the range given on page 7 of 
# Easton et al. 2007 (3.7-18.5mg/kg), then 
# Moore 1993… assume soil density is 2000kg/m^3 of soil
MS_TI05=(((18.5-3.7)/5)*VSAsol$TIClass[1]+3.7)*2000  # range from Easton et 
# al. 2007, pg 7. Moore suggests linear relationship
# We will take care of all of TIClass 05 now as will so 
# it makes sense when you repeat for TI Classes 1-4
# You will use muTS_TI01 thru muTS_TI04 to finish this lab
#
QS= 3.0 # A guess using the middle of the range 1-5
TR=20   # reference Temperature from Table 2.
DPTI05$muS= muTS_TI05*QS^((DPTI05$Tavg-TR)/10)  # Eq. 5
DPTI05$DS=(DPTI05$muS*MS_TI05*DPTI05$Rt)/10^6          # Eq. 4
plot(DPTI05$date,DPTI05$DS)

# DP TI04
DPTI04=data.frame(date=TIC04$date,
                  Rt=TIC04$Qpred/1000, 
                  Tavg=(TIC04$MaxTemp+TIC04$MinTemp)/2)
# Initialize MF and DF
DPTI04$MF=0
DPTI04$DF=0
# Spread your P Fertilizer on ~May 1, ~August 1, and ~October 1
DPTI04$MF[(format(DPTI04$date,"%j") %in% c(121,213,274))]=5.4*10^-4
attach(DPTI04)
#
# Loop to solve MF and DF
for (i in 2:length(date)){
  if(MF[i]<=MF[i-1]){
    MF[i]=MF[i-1]*exp(-dt/tau)-DF[i-1]
  }
  DF[i]=MF[i]*(kF*MF[i]*Rt[i]/(1+kF*MF[i]*Rt[i]))
}
DPTI04$MF=MF
DPTI04$DF=DF
detach(DPTI04)
rm(list=c("MF","DF")) # Clean up the environment 
dev.off() #reset graphics device
plot(DPTI04$date,DPTI04$MF)
plot(DPTI04$date,DPTI04$DF)
muTS_TI04=(((520-0)/5)*VSAsol$TIClass[2]+0) # use VSAsol$TIClass table to 
MS_TI04=(((18.5-3.7)/5)*VSAsol$TIClass[2]+3.7)*2000  # range from Easton et 
QS= 3.0 # A guess using the middle of the range 1-5
TR=20   # reference Temperature from Table 2.
DPTI04$muS= muTS_TI04*QS^((DPTI04$Tavg-TR)/10)  # Eq. 5
DPTI04$DS=(DPTI04$muS*MS_TI04*DPTI04$Rt)/10^6  

# DP TI03
DPTI03=data.frame(date=TIC03$date,
                  Rt=TIC03$Qpred/1000, 
                  Tavg=(TIC03$MaxTemp+TIC03$MinTemp)/2)
# Initialize MF and DF
DPTI03$MF=0
DPTI03$DF=0
# Spread your P Fertilizer on ~May 1, ~August 1, and ~October 1
DPTI03$MF[(format(DPTI03$date,"%j") %in% c(121,213,274))]=5.4*10^-4
attach(DPTI03)
#
# Loop to solve MF and DF
for (i in 2:length(date)){
  if(MF[i]<=MF[i-1]){
    MF[i]=MF[i-1]*exp(-dt/tau)-DF[i-1]
  }
  DF[i]=MF[i]*(kF*MF[i]*Rt[i]/(1+kF*MF[i]*Rt[i]))
}
DPTI03$MF=MF
DPTI03$DF=DF
detach(DPTI03)
rm(list=c("MF","DF")) # Clean up the environment 
dev.off() #reset graphics device
plot(DPTI03$date,DPTI03$MF)
plot(DPTI03$date,DPTI03$DF)
muTS_TI03=(((520-0)/5)*VSAsol$TIClass[3]+0) # use VSAsol$TIClass table to 
MS_TI03=(((18.5-3.7)/5)*VSAsol$TIClass[3]+3.7)*2000  # range from Easton et 
DPTI03$muS= muTS_TI03*QS^((DPTI03$Tavg-TR)/10)  # Eq. 5
DPTI03$DS=(DPTI03$muS*MS_TI03*DPTI03$Rt)/10^6  

# DP TI02
DPTI02=data.frame(date=TIC02$date,
                  Rt=TIC02$Qpred/1000, 
                  Tavg=(TIC02$MaxTemp+TIC02$MinTemp)/2)
# Initialize MF and DF
DPTI02$MF=0
DPTI02$DF=0
# Spread your P Fertilizer on ~May 1, ~August 1, and ~October 1
DPTI02$MF[(format(DPTI02$date,"%j") %in% c(121,213,274))]=5.4*10^-4
attach(DPTI02)
#
# Loop to solve MF and DF
for (i in 2:length(date)){
  if(MF[i]<=MF[i-1]){
    MF[i]=MF[i-1]*exp(-dt/tau)-DF[i-1]
  }
  DF[i]=MF[i]*(kF*MF[i]*Rt[i]/(1+kF*MF[i]*Rt[i]))
}
DPTI02$MF=MF
DPTI02$DF=DF
detach(DPTI02)
rm(list=c("MF","DF")) # Clean up the environment 
dev.off() #reset graphics device
plot(DPTI02$date,DPTI02$MF)
plot(DPTI02$date,DPTI02$DF)
muTS_TI02=(((520-0)/5)*VSAsol$TIClass[4]+0) # use VSAsol$TIClass table to 
MS_TI02=(((18.5-3.7)/5)*VSAsol$TIClass[4]+3.7)*2000  # range from Easton et 
DPTI02$muS= muTS_TI02*QS^((DPTI02$Tavg-TR)/10)  # Eq. 5
DPTI02$DS=(DPTI02$muS*MS_TI02*DPTI02$Rt)/10^6  

# DP TI01
DPTI01=data.frame(date=TIC01$date,
                  Rt=TIC01$Qpred/1000, 
                  Tavg=(TIC01$MaxTemp+TIC01$MinTemp)/2)
# Initialize MF and DF
DPTI01$MF=0
DPTI01$DF=0
# Spread your P Fertilizer on ~May 1, ~August 1, and ~October 1
DPTI01$MF[(format(DPTI01$date,"%j") %in% c(121,213,274))]=5.4*10^-4
attach(DPTI01)
#
# Loop to solve MF and DF
for (i in 2:length(date)){
  if(MF[i]<=MF[i-1]){
    MF[i]=MF[i-1]*exp(-dt/tau)-DF[i-1]
  }
  DF[i]=MF[i]*(kF*MF[i]*Rt[i]/(1+kF*MF[i]*Rt[i]))
}
DPTI01$MF=MF
DPTI01$DF=DF
detach(DPTI01)
rm(list=c("MF","DF")) # Clean up the environment 
dev.off() #reset graphics device
plot(DPTI01$date,DPTI01$MF)
plot(DPTI01$date,DPTI01$DF)
muTS_TI01=(((520-0)/5)*VSAsol$TIClass[5]+0) # use VSAsol$TIClass table to 
MS_TI01=(((18.5-3.7)/5)*VSAsol$TIClass[5]+3.7)*2000  # range from Easton et 
DPTI01$muS= muTS_TI01*QS^((DPTI01$Tavg-TR)/10)  # Eq. 5
DPTI01$DS=(DPTI01$muS*MS_TI01*DPTI01$Rt)/10^6  


# We will assume a simple base flow model, where the stream baseflow 
# B equals the minimum measured flow of the basin. We will build the 
# Base Flow solution directly into the Stream Load dataframe (DPLT).
# Note that we only consider TIC05 here because TIC05 is the class that 
# directly contributes to the stream.
DPLT=data.frame(date=TIC05$date,
                Rt=TIC05$Qpred/1000.0,
                Tavg=(TIC05$MaxTemp+TIC05$MinTemp)/2)
DPLT$B=min(TIC05$Qmm)*myflowgage$area*1000*1000/1000 # m^3/day
muTB=2.1*10^(-5) # Easton Table 2
QB=2.2           # Easton Table 2
TB=17            # Easton Table 2
DPLT$muB=muTB*QB^((DPLT$Tavg-TB)/10)  # Easton eq. 10
DPLT$LB=DPLT$muB*DPLT$B     # Easton eq. 9
plot(DPLT$date,DPLT$LB)

TIC05$area=myflowgage$area/5
TIC04$area=myflowgage$area/5
TIC03$area=myflowgage$area/5
TIC02$area=myflowgage$area/5
TIC01$area=myflowgage$area/5

DPLT$LT=DPLT$LB +
  TIC05$area*(DPTI05$DF + DPTI05$DS)+
  TIC04$area*(DPTI04$DF + DPTI04$DS)+
  TIC03$area*(DPTI03$DF + DPTI03$DS)+
  TIC02$area*(DPTI02$DF + DPTI02$DS)+
  TIC01$area*(DPTI01$DF + DPTI01$DS)
plot(DPLT$date,DPLT$LT, type="l", ylim=c(1,10))

# HW 1
TotalP = (TIC05$area*DPTI05$MF)+(TIC04$area*DPTI04$MF)+(TIC03$area*DPTI03$MF)+(TIC02$area*DPTI02$MF)+(TIC01$area*DPTI01$MF)
TotalP
plot(DPLT$date, cumsum(TotalP), type="l")