######################################################################
# NOTES, PACKAGES, AND DATA
######################################################################
## This code is used to aggregate surveys and calculate prevalence at the AD level
## Read coviariate code FIRST to access relevant shapefiles

## PACKAGES
library(dplyr); library(sp)
library(spatialEco); library(rgdal)

## DATA
## WWARN = Worldwide Antimalarial Resistance Network 
## Dataset adapted from the PDMS database in WWARN
setwd("~/desktop/MalDrugRes_SSA")
wwarn <- read.csv('Survey_MolecMarker_Data.csv')
wwarn$MidYear <- round((wwarn$StartYr+wwarn$EndYr)/2,0)

######################################################################
# Creating spatial points for surveys
######################################################################

## Choose marker of interest(uncomment)
wwarn <- wwarn[wwarn$MrkrType=="pfcrt 76T",]
#wwarn <- wwarn[wwarn$MrkrType=="pfmdr1 86Y",]
#wwarn <- wwarn[wwarn$MrkrType=="pfmdr1 184F",]
#wwarn <- wwarn[wwarn$MrkrType=="pfmdr1 1246Y",]

## Define time period of interest (uncomment)
wwarn <- wwarn[wwarn$MidYear<=2009,]
#wwarn <- wwarn[wwarn$MidYear>2009,]

wwarn.points <- as.data.frame(cbind(x=as.numeric(wwarn$Lon),y=as.numeric(wwarn$Lat),
                ID=as.character(wwarn$ID),sampyear=as.numeric(wwarn$MidYear),pubyear=as.numeric(wwarn$PubYear),
                country=as.character(wwarn$Country),tested=as.numeric(wwarn$Tested),pres=as.numeric(wwarn$MixedPres)))
wwarn.points <- wwarn.points[!is.na(wwarn.points$x),]
wwarn.points <- wwarn.points[wwarn.points$country!="Comoros",] ## excluded in CAR analysis
wwarn.points <-distinct(wwarn.points, x, y, ID, sampyear, .keep_all= TRUE)
wwarn.points$pt.ids <- 1:length(wwarn.points$x)
wwarn.points.xy <- wwarn.points[,c("x","y","sampyear","tested","pres","pt.ids")]
wwarn.points.xy$x <- as.numeric(as.character(wwarn.points.xy$x))
wwarn.points.xy$y <- as.numeric(as.character(wwarn.points.xy$y))

######################################################################
# Intersecting points and polygons (to calculate mean observed prevelance)
######################################################################

## Reminder to read in covariate/shapefile code first 
wwarn.points.sp <- wwarn.points.xy
coordinates(wwarn.points.sp) <- ~x+y
wwarn.points.sp = SpatialPoints(wwarn.points.sp,proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
plot(wwarn.points.sp,pch=16, cex=0.5)
pts <- point.in.poly(wwarn.points.sp, ctry_shps)
pts.df <- as.data.frame(pts)
wwarn.merge <- merge(wwarn.points.xy,pts.df,by='pt.ids')
wwarn.merge$tested = as.numeric(as.character(wwarn.merge$tested))
wwarn.merge$pres = as.numeric(as.character(wwarn.merge$pres))
wwarn.merge$NAME_1 = as.character(wwarn.merge$NAME_1)
wwarn.merge$group_id <- wwarn.merge %>% group_indices(NAME_1) 

n = length(unique(wwarn.merge$group_id))
df.poly <- data.frame("group_id"=unique(sort(wwarn.merge$group_id)),"N"=rep(NA,n),"Mut"=rep(NA,n))
for (i in 1:n) {
  df.poly$N[i]=sum(wwarn.merge$tested[wwarn.merge$group_id==i])
  df.poly$Mut[i]=sum(wwarn.merge$pres[wwarn.merge$group_id==i])}
df.poly$prev <- df.poly$Mut/df.poly$N
df.poly$NAME_1 <- wwarn.merge$NAME_1[match(df.poly$group_id, wwarn.merge$group_id)]
df.poly$NAME_0 <- wwarn.merge$NAME_0[match(df.poly$group_id, wwarn.merge$group_id)]

## Merging binomial prev. data with spatial polygon data
ctry_shps@data <- left_join(ctry_shps@data, df.poly, 
                  by=c("NAME_1" = "NAME_1", "NAME_0" = "NAME_0"), all.x = TRUE)

##########################################################
# Save observed data aggregated by admin. division
##########################################################

setwd("~/desktop/MalDrugRes_SSA/Observed_Prevalence")

## If rerunning code or using different markers/time periods, save agg. R data using: 
writeOGR(obj=ctry_shps, dsn="XXX", layer="ctry_shps1", driver="ESRI Shapefile")
## Replace XXX with name of marker and time period, e.g. crt_time1
## Delete layers if re-running code for same marker/time; will not overwrite!


