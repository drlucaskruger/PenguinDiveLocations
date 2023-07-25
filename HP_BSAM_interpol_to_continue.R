library(track2KBA)
library(lubridate)
library(dplyr)
library(ggplot2)
library(sp)
library(adehabitatLT)
library(bsam)
library(ggplot2)
library(patchwork)
library(data.table)
library(tidyverse)
library(move)
library(lfstat)

crs=CRS("+proj=longlat +datum=WGS84 +no_defs")
spo=CRS('+init=epsg:3031')


####---------Load data and run BSAM DCRW over GPS data-------------



#### ---------2019/20-----------

### load GPS and dive data

setwd("D:/UruguayPenguins/Harmony/GPS/2019/")
# get a vector of all file names
myfiles <- list.files("D:/UruguayPenguins/Harmony/GPS/2019/")

# loop over files names, reading in and saving each data.frame as an element in a list
n <- length(myfiles )
datalist <- vector(mode="list", length=n)
for(i in 1:n) {
  cat("importing file", i, ":", myfiles[i], "\n")
  datalist[[i]] <- read.csv(myfiles[i])
}

# assign list elements the file names
names(datalist) <- myfiles 

# combine all data.frames in datalist, use idcol argument to assign original file name
GPS19 <- data.table::rbindlist(datalist, idcol=TRUE,fill=T)


setwd("D:/UruguayPenguins/Harmony/TDR/2019/")
# get a vector of all file names
myfiles <- list.files("D:/UruguayPenguins/Harmony/TDR/2019/")

# loop over files names, reading in and saving each data.frame as an element in a list
n <- length(myfiles )
datalist <- vector(mode="list", length=n)
for(i in 1:n) {
  cat("importing file", i, ":", myfiles[i], "\n")
  datalist[[i]] <- read.csv(myfiles[i])
}

# assign list elements the file names
names(datalist) <- myfiles 

# combine all data.frames in datalist, use idcol argument to assign original file name
TDR19 <- data.table::rbindlist(datalist, idcol=TRUE,fill=T)

TDR19$id<-TDR19$'.id'

head(TDR19)

TDR19$begdesc<-as.POSIXct(strptime(TDR19$begdesc, format="%Y-%m-%d %H:%M:%S", tz="GMT"))

#save GPS and TDR data for later

saveRDS(GPS19,"D:/UruguayPenguins/Harmony/GPS/2019/GPS2019.Rds")
saveRDS(TDR19,"D:/UruguayPenguins/Harmony/TDR/2019/TDR2019.Rds")

###-------------Process GPS data------------------
head(GPS19)

GPS19$date<-as.POSIXct(strptime(GPS19$Timestamp, format="%d/%m/%Y %H:%M:%S", tz="GMT"))

df19 <- GPS19[!duplicated(GPS19[, c("Lat", "Long","date")]), ]
df19<-subset(df19,Lat<(-62) & Lat>(-63))

#check if there are potentially wrong GPS fixes
ggplot(df19,aes(Long,Lat))+geom_point()+facet_wrap(.id~.)  # PS: ID INC mean the last 5 days of december; ID CR means first half of ja

#speed filter PS: when I applied the speed/distance filter, the SSM did not work
#so I am moving forward without filtering for thre reasons: 
# - the SSM might eliminate most aberrant locations
# - inspection of prediction of each individual can indicate issues with aberrant locations
# - if necessary, it is still possible to apply a posteriori filters


#spdf19<-SpatialPointsDataFrame(coords  = coordinates(df19[,c("Long","Lat")]),df19,proj4string = crs)
#spdf19<-spTransform(spdf19,CRSobj=spo,center=FALSE)
#plot(spdf19)
#traj19<-as.ltraj(xy=coordinates(spdf19),date=spdf19$date,id=spdf19$id,
#                typeII = TRUE,infolocs = spdf19,
#                slsp = c("missing"),
                
#                proj4string = spo)
#plot(traj19)
#tdf<-ld(traj19)
#summary(tdf$dist)
#ggplot(tdf,aes(dist))+geom_histogram()+xlim(0,1000)
#tdf<-subset(tdf,dist<501)
#ggplot(subset(df19,id=="HP18_INC19"),aes(Long,Lat))+geom_point()+facet_wrap(.id~.)+theme_bw()+
#ggplot(subset(tdf,id=="HP18_INC19"),aes(x,y))+geom_point()+facet_wrap(id~.)+theme_bw()

dfMP19<-data.frame(id=df19$id,date=df19$date,lc=c(3),lon=df19$Long,lat=df19$Lat) # using the maximum accuracy in argos data

# run the SSM. PS. When I did a pre-filter of positions within the colony site, the prediction ability of the model
# decreased substantially, so I eliminated the points estimated in the colony area after the SSM
# With this parameters the SSM is time consuming. With these parameterization, each animal takes between 
# 4 and 6 minutes to be processed. I tested separately a few individuals, and the better models were those with 
# samples above 1000, large thin value and very low span value. tsteps of 0.0065 corresponds to a 10min interval between locations

bs19<-fit_ssm(dfMP19,model="DCRW",tstep = 0.0065, adapt = 50, samples = 2000, 
              thin = 5, span = 0.05)      

diag_ssm(bs19)  # diagnostic statistics for each animal 
                #trace and density of different posteriors distribution should be similar, indicating model convergence
                #GRB shrink factor should not be too far from 1.   

#map_ssm(bs19) # all animals together
map_ssm(bs19,onemap = F) #animals separated. HP09_INC19, HP18_INC19 and HP16_INC19 needs further processing

ssm19<-get_summary(bs19) # extract model output

saveRDS(ssm19,"D:/UruguayPenguins/Harmony/SSM19.Rds")

rm(list=ls(all.names=TRUE))
gc()
.rs.restartR()


###------------2021/22------------------


crs=CRS("+proj=longlat +datum=WGS84 +no_defs")
spo=CRS('+init=epsg:3031')

### load GPS and dive data

setwd("D:/UruguayPenguins/Harmony/GPS/2021/")
# get a vector of all file names
myfiles <- list.files("D:/UruguayPenguins/Harmony/GPS/2021/")

# loop over files names, reading in and saving each data.frame as an element in a list
n <- length(myfiles )
datalist <- vector(mode="list", length=n)
for(i in 1:n) {
  cat("importing file", i, ":", myfiles[i], "\n")
  datalist[[i]] <- read.csv(myfiles[i])
}

# assign list elements the file names
names(datalist) <- myfiles 

# combine all data.frames in datalist, use idcol argument to assign original file name
GPS21 <- data.table::rbindlist(datalist, idcol=TRUE,fill=T)


setwd("D:/UruguayPenguins/Harmony/TDR/2021/")
# get a vector of all file names
myfiles <- list.files("D:/UruguayPenguins/Harmony/TDR/2021/")

# loop over files names, reading in and saving each data.frame as an element in a list
n <- length(myfiles )
datalist <- vector(mode="list", length=n)
for(i in 1:n) {
  cat("importing file", i, ":", myfiles[i], "\n")
  datalist[[i]] <- read.csv(myfiles[i])
}

# assign list elements the file names
names(datalist) <- myfiles 

# combine all data.frames in datalist, use idcol argument to assign original file name
TDR21 <- data.table::rbindlist(datalist, idcol=TRUE,fill=T)

TDR21$id<-TDR21$'.id'

TDR21$begdesc<-as.POSIXct(strptime(TDR21$begdesc, format="%Y-%m-%d %H:%M:%S", tz="GMT"))

head(TDR21)


#save GPS and TDR data for later

saveRDS(GPS21,"D:/UruguayPenguins/Harmony/GPS/2021/GPS2021.Rds")
saveRDS(TDR21,"D:/UruguayPenguins/Harmony/TDR/2021/TDR2021.Rds")
###-------------Process GPS data------------------
head(GPS21)

GPS21$date<-as.POSIXct(strptime(GPS21$Timestamp, format="%d/%m/%Y %H:%M:%S", tz="GMT"))

df21 <- GPS21[!duplicated(GPS21[, c("location.lat", "location.lon","date")]), ]
df21<-subset(df21,location.lat<(-62) & location.lat>(-63))

#check if there are potentially wrong GPS fixes
ggplot(df21,aes(location.lon,location.lat))+geom_point()+facet_wrap(.id~.)  # PS: ID INC mean the last 5 days of december; ID CR means first half of ja

dfMP21<-data.frame(id=df21$'.id',date=df21$date,lc=c(3),lon=df21$location.lon,lat=df21$location.lat) # using the maximum accuracy in argos data
head(dfMP21)

bs21<-fit_ssm(dfMP21,model="DCRW",tstep = 0.0065, adapt = 50, samples = 2000, 
              thin = 5, span = 0.05)      

diag_ssm(bs21)  # diagnostic statistics for each animal 

#map_ssm(bs21) # all animals together
map_ssm(bs21,onemap = F) #animals separated.

ssm21<-get_summary(bs21) # extract model output

saveRDS(ssm21,"D:/UruguayPenguins/Harmony/SSM21.Rds")

rm(list=ls(all.names=TRUE))
gc()
.rs.restartR()


#-------------2022/23------------------


crs=CRS("+proj=longlat +datum=WGS84 +no_defs")
spo=CRS('+init=epsg:3031')

### -------GPS data

setwd("D:/UruguayPenguins/Harmony/GPS/2022/")
# get a vector of all file names
myfiles <- list.files("D:/UruguayPenguins/Harmony/GPS/2022/")

# loop over files names, reading in and saving each data.frame as an element in a list
n <- length(myfiles )
datalist <- vector(mode="list", length=n)
for(i in 1:n) {
  cat("importing file", i, ":", myfiles[i], "\n")
  datalist[[i]] <- read.csv(myfiles[i])
}

# assign list elements the file names
names(datalist) <- myfiles 

# combine all data.frames in datalist, use idcol argument to assign original file name
GPS22 <- data.table::rbindlist(datalist, idcol=TRUE,fill=T)

###--------TDR data
setwd("D:/UruguayPenguins/Harmony/TDR/2022/")
# get a vector of all file names
myfiles <- list.files("D:/UruguayPenguins/Harmony/TDR/2022/")

# loop over files names, reading in and saving each data.frame as an element in a list
n <- length(myfiles )
datalist <- vector(mode="list", length=n)
for(i in 1:n) {
  cat("importing file", i, ":", myfiles[i], "\n")
  datalist[[i]] <- read.csv(myfiles[i])
}

# assign list elements the file names
names(datalist) <- myfiles 

# combine all data.frames in datalist, use idcol argument to assign original file name
TDR22 <- data.table::rbindlist(datalist, idcol=TRUE,fill=T)

TDR22$id<-TDR22$'.id'

TDR22$begdesc<-as.POSIXct(strptime(TDR22$begdesc, format="%Y-%m-%d %H:%M:%S", tz="GMT"))

head(TDR22)


#save GPS and TDR data for later

saveRDS(GPS22,"D:/UruguayPenguins/Harmony/GPS/2022/GPS2022.Rds")
saveRDS(TDR22,"D:/UruguayPenguins/Harmony/TDR/2022/TDR2022.Rds")

###-------------Process GPS data------------------
head(GPS22)

GPS22$date<-as.POSIXct(strptime(paste(GPS22$Date,GPS22$Time), format="%d/%m/%Y %H:%M:%S", tz="GMT")) # in this set date and time were separated
head(GPS22)

df22 <- GPS22[!duplicated(GPS22[, c("location.lat", "location.lon","date")]), ]
df22<-subset(df22,location.lat<(-62) & location.lat>(-63))

#check if there are potentially wrong GPS fixes
ggplot(df22,aes(location.lon,location.lat))+geom_point()+facet_wrap(.id~.)  

dfMP22<-data.frame(id=df22$'.id',date=df22$date,lc=c(3),lon=df22$location.lon,lat=df22$location.lat) # using the maximum accuracy in argos data
head(dfMP22)

bs22<-fit_ssm(dfMP22,model="DCRW",tstep = 0.0065, adapt = 50, samples = 2000, 
              thin = 5, span = 0.05)      

diag_ssm(bs22)  # diagnostic statistics for each animal. Some fits were not too good...

#map_ssm(bs22) # all animals together
map_ssm(bs22,onemap = F) #animals separated.

#bad fits: 2022_KI-5_R4_2bV; 2022_KI-20_R2_6bV; 2022_KI-15_R2_6bA; 2022_DI-13_R2_5bA
#          2022_CP-07_R3_3A;2022_CP-07_R2_4V; 2022_CP-06_R2_3A; 2022_CP-05_R2_3V;
#          2022_CP-02_R2_1A; 2022_CP-01_R2_1V  - so out of 32 tracking events, 22 had useful data. It is consistent with other years

ssm22<-get_summary(bs22) # extract model output

saveRDS(ssm22,"D:/UruguayPenguins/Harmony/SSM22.Rds")


rm(list=ls(all.names=TRUE))
gc()
.rs.restartR()


####----------Merge BSAM DCRW and GPS------------------

#### merge DCRW predicted locations and observed locations, eliminate NA values and filter using moving average


#previously saved GPS and DCRW data

gps19<-readRDS("D:/UruguayPenguins/Harmony/GPS/2019/GPS2019.Rds")
gps21<-readRDS("D:/UruguayPenguins/Harmony/GPS/2021/GPS2021.Rds")
gps22<-readRDS("D:/UruguayPenguins/Harmony/GPS/2022/GPS2022.Rds")

ssm19<-readRDS("D:/UruguayPenguins/Harmony/SSM19.Rds")
ssm21<-readRDS("D:/UruguayPenguins/Harmony/SSM21.Rds")
ssm22<-readRDS("D:/UruguayPenguins/Harmony/SSM22.Rds")


#------2019----------
gps19$date<-as.POSIXct(strptime(gps19$Timestamp, format="%d/%m/%Y %H:%M:%S", tz="GMT"))
gps19$id<-substring(gps19$'.id',first=1,last=9)
ssm19$date<-as.POSIXct(strptime(ssm19$date, format="%Y-%m-%d %H:%M:%S", tz="GMT"))

merge19<-merge(gps19,ssm19,by=c("id","date"),all=T)
head(merge19)

merge19$LongM<-ifelse(is.na(merge19$Long),merge19$lon,merge19$Long)
merge19$LatM<-ifelse(is.na(merge19$Lat),merge19$lat,merge19$Lat)

# moving average
merge19$lonmean<-ma(merge19$LongM,n=5,sides=2)
merge19$latmean<-ma(merge19$LatM,n=5,sides =2)


#-------2021-------------
head(gps21)
head(ssm21)

gps21$date<-as.POSIXct(strptime(gps21$Timestamp, format="%d/%m/%Y %H:%M:%S", tz="GMT"))
gps21$id<-(gps21$'.id')
ssm21$date<-as.POSIXct(strptime(ssm21$date, format="%Y-%m-%d %H:%M:%S", tz="GMT"))

merge21<-merge(gps21,ssm21,by=c("id","date"),all=T)
head(merge21)

merge21$LongM<-ifelse(is.na(merge21$location.lon),merge21$lon,merge21$location.lon)
merge21$LatM<-ifelse(is.na(merge21$location.lat),merge21$lat,merge21$location.lat)

# moving average
merge21$lonmean<-ma(merge21$LongM,n=5,sides=2)
merge21$latmean<-ma(merge21$LatM,n=5,sides =2)

####----------2022/23-------------

head(gps22)
head(ssm22)

gps22$date<-as.POSIXct(strptime(paste(gps22$Date,gps22$Time), format="%d/%m/%Y %H:%M:%S", tz="GMT"))
gps22$id<-(gps22$'.id')
ssm21$date<-as.POSIXct(strptime(ssm21$date, format="%Y-%m-%d %H:%M:%S", tz="GMT"))

merge22<-merge(gps22,ssm22,by=c("id","date"),all=T)
head(merge22)

merge22$LongM<-ifelse(is.na(merge22$location.lon),merge22$lon,merge22$location.lon)
merge22$LatM<-ifelse(is.na(merge22$location.lat),merge22$lat,merge22$location.lat)

# moving average
merge22$lonmean<-ma(merge22$LongM,n=5,sides=2)
merge22$latmean<-ma(merge22$LatM,n=5,sides =2)

###-------application of speed filters -----------

library(SDLfilter)

sdl19<-data.frame(id=merge19$id,DateTime=merge19$date,lat=merge19$LatM,lon=merge19$LongM,qi=c(5))
vmax19<-vmax(sdl19, qi = 5, method = "ML", prob = 0.5)
dd19<-ddfilter_speed(sdl19,vmax=vmax19)

ggplot(subset(merge19,id=="HP03_CR19"),aes(Long,Lat))+geom_point()+theme_bw()+ggtitle(label="a.Raw locations")+
ggplot(subset(merge19,id=="HP03_CR19"),aes(lon,lat))+geom_point()+theme_bw()+ggtitle(label="b.DCRW locations")+
ggplot(subset(merge19,id=="HP03_CR19"),aes(LongM,LatM))+geom_point()+theme_bw()+ggtitle(label="c.Reconstructed locations")+
ggplot(subset(dd19,id=="HP03_CR19"),aes(lon,lat))+geom_point()+theme_bw()+ggtitle(label="d.Filtered locations")

head(merge21)

sdl21<-data.frame(id=merge21$id,DateTime=merge21$date,lat=merge21$LatM,lon=merge21$LongM,qi=c(5))
vmax21<-vmax(sdl21, qi = 5, method = "ML", prob = 0.5)
dd21<-ddfilter_speed(sdl21,vmax=vmax21)

ggplot(subset(merge21,id=="HP03_CR21.csv"),aes(location.lon,location.lat))+geom_point()+theme_bw()+ggtitle(label="a.Raw locations")+
  ggplot(subset(merge21,id=="HP03_CR21.csv"),aes(lon,lat))+geom_point()+theme_bw()+ggtitle(label="b.DCRW locations")+
  ggplot(subset(merge21,id=="HP03_CR21.csv"),aes(LongM,LatM))+geom_point()+theme_bw()+ggtitle(label="c.Reconstructed locations")+
  ggplot(subset(dd21,id=="HP03_CR21.csv"),aes(lon,lat))+geom_point()+theme_bw()+ggtitle(label="d.Filtered locations")



sdl22<-data.frame(id=merge22$id,DateTime=merge22$date,lat=merge22$LatM,lon=merge22$LongM,qi=c(5))
vmax22<-vmax(sdl22, qi = 5, method = "ML", prob = 0.5)
dd22<-ddfilter_speed(sdl22,vmax=vmax22)

ggplot(subset(merge22,id=="2022_KI-3_R4_1V.csv"),aes(location.lon,location.lat))+geom_point()+theme_bw()+ggtitle(label="a.Raw locations")+
  ggplot(subset(merge22,id=="2022_KI-3_R4_1V.csv"),aes(lon,lat))+geom_point()+theme_bw()+ggtitle(label="b.DCRW locations")+
  ggplot(subset(merge22,id=="2022_KI-3_R4_1V.csv"),aes(LongM,LatM))+geom_point()+theme_bw()+ggtitle(label="c.Reconstructed locations")+
  ggplot(subset(dd22,id=="2022_KI-3_R4_1V.csv"),aes(lon,lat))+geom_point()+theme_bw()+ggtitle(label="d.Filtered locations")




#-------------- Dive Locations Interpolations

crs=CRS("+proj=longlat +datum=WGS84 +no_defs")



head(dd19)
#dd19$id<-paste(dd19$id,c("csv"),sep=".") # ids in 2019 had not the ".csv"

ddf<-rbind(dd19,dd21,dd22)
head(ddf)
ddf<-subset(ddf,lat>(-63))


#### join dive data together

tdr19<-readRDS("D:/UruguayPenguins/Harmony/TDR/2019/TDR2019.Rds")
tdr21<-readRDS("D:/UruguayPenguins/Harmony/TDR/2021/TDR2021.Rds")
tdr22<-readRDS("D:/UruguayPenguins/Harmony/TDR/2022/TDR2022.Rds")

head(tdr19)

tdr19b<-data.frame(id=tdr19$id,begdesc=tdr19$begdesc,tdr19[,6:38])

head(tdr21)
tdr21b<-data.frame(id=tdr21$id,begdesc=tdr21$begdesc,tdr21[,6:38])

head(tdr22)
tdr22b<-data.frame(id=tdr22$id,begdesc=tdr22$begdesc,tdr22[,6:38])

tdr<-rbind(tdr19b,tdr21b,tdr22b)

head(tdr)

summary(as.factor(ddf$id))  ### some IDs had errors, no idea why. so lets correct it

ddf$id[ddf$id=="HP03_INC1.csv"]<-"HP03_INC19.csv"
ddf$id[ddf$id=="HP06_INC1.csv"]<-"HP06_INC19.csv"
ddf$id[ddf$id=="HP07_INC1.csv"]<-"HP07_INC19.csv"
ddf$id[ddf$id=="HP09_INC1.csv"]<-"HP09_INC19.csv"


ddf <- ddf[!duplicated(ddf[, c("DateTime","id")]), ]


spdf<-SpatialPointsDataFrame(coords  = coordinates(ddf[,c("lon","lat")]),ddf,proj4string = crs)

#plot(spdf)

peng_track <- as.ltraj (xy = coordinates(spdf), date = spdf$DateTime,#typeII=T,
                          id = spdf$id,
                          #burst = spdf19$id,
                          proj4string = crs) 

mtrack<-move(peng_track)
plot(mtrack)
summary(as.factor(tdr$id))

combined_TDR <- data.frame()
for (i in unique(tdr$id)) {
  # Check if the ID exists in mtrack19
  if (i %in% mtrack@idData$id) {
    # Subset TDR_filtered data for the current burst based on the minimum and maximum dates in mtrack19
    TDR_filtered <- subset(tdr, id == i & begdesc >= min(mtrack@timestamps) & begdesc <= max(mtrack@timestamps))
    
    # Combine the data for the current animal with the existing combined_TDR data
    combined_TDR <- rbind(combined_TDR, TDR_filtered)
    
  } else {
    # Handle the case when the ID doesn't exist in mtrack19 
    warning(paste("No corresponding GPS data for animal ID", i))
  }
}


summary(as.factor(combined_TDR$id))

tdr_traj<- as.ltraj (xy = coordinates(combined_TDR[,c("desctim","asctim")]), 
                     date = combined_TDR$begdesc,#typeII=T,
                     id = combined_TDR$id)

tdr_move<-move::move(tdr_traj)

animal_ids <- unique(tdr_move@idData$id)  # Use animal_ids from tdr_move


# Interpolate GPS positions for each animal in mtrack19
track_interpolated <- data.frame()

animal_ids <- unique(tdr_move@idData$id)  # Use animal_ids from tdr_move

# Filter mtrack19 based on animal_ids from tdr_move
mtrack_filtered <- mtrack[mtrack@idData$id %in% animal_ids]

# Initialize an empty list to store the interpolated tracks
interpolated_tracks <- list()

for (animal_id in animal_ids) {
  # Get the GPS track for the current animal
  track_sel <- mtrack_filtered[[which(mtrack_filtered@idData$id == animal_id)]]
  
  # Get the corresponding diving data for the current animal
  tdr_data <- tdr_move[[which(tdr_move@idData$id == animal_id)]]
  
  # Check the time range of GPS track and diving data
  gps_start_time <- min(track_sel@timestamps)
  gps_end_time <- max(track_sel@timestamps)
  
  diving_start_time <- min(tdr_data@timestamps)
  diving_end_time <- max(tdr_data@timestamps)
  
  cat("Animal ID:", animal_id, "\n")
  cat("GPS Start Time:", gps_start_time, "GPS End Time:", gps_end_time, "\n")
  cat("Diving Start Time:", diving_start_time, "Diving End Time:", diving_end_time, "\n")
  cat("\n")
  
  # Use tryCatch to handle the error and continue with the next animal
  tryCatch({
    # Interpolate the positions based on the diving data
    t_interp <- move::interpolateTime(track_sel, time = tdr_data@timestamps, spaceMethod = 'euclidean')
    
    # Store the interpolated track in the list
    interpolated_tracks[[animal_id]] <- t_interp
  }, error = function(e) {
    # Handle the error (optional)
    # For example, you can print a warning message:
    warning(paste("No corresponding diving data for animal ID", animal_id))
  })
}

# Combine all the interpolated tracks into a single move stack
track_interpolated <- move::moveStack(interpolated_tracks)

summary(track_interpolated)

trackdf<-as.data.frame(track_interpolated)
head(trackdf)
summary(as.factor(trackdf$id))

head(tdr)
head(trackdf)

tdr$timestamps<-tdr$begdesc


divelocs<-merge(trackdf,tdr,by=c("id","timestamps"))
summary(as.factor(divelocs$id))

write.csv(divelocs,"D:/UruguayPenguins/Harmony/DiveLocs.csv")

saveRDS(divelocs,"D:/UruguayPenguins/Harmony/DiveLocs.Rds")


####--------classify trips-------------

dfSM<-data.frame(id=divelocs$id,date=divelocs$timestamps,lon=divelocs$x,lat=divelocs$y)

dfSM$date_gmt<-paste(year(dfSM$date),month(dfSM$date),day(dfSM$date),sep="-")
dfSM$time_gmt<-paste(hour(dfSM$date),minute(dfSM$date),second(dfSM$date),sep=":")


dataGroup <- formatFields(
  dataGroup = dfSM, 
  fieldID   = "id", 
  fieldDate = "date_gmt", 
  fieldTime = "time_gmt",
  fieldLon  = "lon", 
  fieldLat  = "lat"
)
str(dataGroup)

colony <- dataGroup %>% 
  summarise(
    Longitude = first(Longitude), 
    Latitude  = first(Latitude)
  )


trips <- tripSplit(
  dataGroup  = dataGroup,
  colony     = colony,
  innerBuff  = 3,      # kilometers
  returnBuff = 10,
  duration   = 1,      # hours
  rmNonTrip  = TRUE
)

write.csv(trips,"D:/UruguayPenguins/Harmony/DiveLocs_trips.csv")
saveRDS(trips,"D:/UruguayPenguins/Harmony/DiveLocs_trips.Rds")

#calculate trips caracteristics

tripsy <- subset(trips, trips$Returns == "Yes")

sumTrips <- tripSummary(trips = tripsy, colony = colony)

sumTrips


write.csv(sumTrips,"D:/UruguayPenguins/Harmony/summary_of_trips.csv")
saveRDS(sumTrips,"D:/UruguayPenguins/Harmony/summary_of_trips.Rds")


tracks <- projectTracks( dataGroup = trips, projType = 'azim', custom=TRUE )
class(tracks)

hVals <- findScale(
  tracks   = tracks,
  scaleARS = TRUE,
  sumTrips = sumTrips,
  res=1)

hVals # smoothing parameter

tracks <- tracks[tracks$ColDist > 3, ] # remove trip start and end points near colony

KDE <- estSpaceUse(
  tracks = tracks, 
  scale = 1, 
  levelUD = 75, 
  polyOut = TRUE,
  res=1
)

mapKDE(KDE = KDE$UDPolygons, colony = colony)

repr <- repAssess(
  tracks    = tracks, 
  KDE       = KDE$KDE.Surface, 
  levelUD   = 50,
  iteration = 30, 
  bootTable = FALSE)

Site <- findSite(
  KDE = KDE$KDE.Surface,
  represent = repr$out,
  levelUD = 70,
  popSize = 40000,     
  polyOut = TRUE
)

class(Site)

potSite <- Site %>% dplyr::filter(.data$potentialSite==TRUE) %>% 
  summarise(
    max_animals = max(na.omit(N_animals)), # maximum number of animals aggregating in the site
    min_animals = min(na.omit(N_animals))  # minimum number using the site
  )


mapSite(Site)


### Exploratory Analysis of Trips Characteristics

#to continue