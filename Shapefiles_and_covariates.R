######################################################################
# NOTES, PACKAGES, AND DATA
######################################################################

## This code is used to download relevant shapefiles and covariate data
## Run this code before merging molecular marker survey data 

library(GADMTools); library(dplyr)
library(raster); library(sp); library(sf)
library(spdep); library(spatstat); library(rgdal); library(rgeos)
library(magrittr); library(spatialEco); library(largeList)
library(maptools); library('malariaAtlas')
library(reshape2); library(PrevMap); library(tmap)
library(ggplot2); library(RColorBrewer)

######################################################################
# Reading in the Shapefiles
######################################################################

setwd("~/desktop/MalDrugRes_SSA/gadm")
ssa = list()
ssa$countries = c("AGO","BDI","BEN","BFA","BWA","CAF","CIV","CMR","COD","COG","ERI",
                  "ETH","GAB","GHA","GIN","GMB","GNB","GNQ","KEN","LBR","MDG","MLI",
                  "MOZ","MRT","MWI","NAM","NER","NGA","RWA","SDN","SEN","SLE","SOM",
                  "SSD","SWZ","TCD","TGO","TZA","UGA","ZAF","ZMB","ZWE","LSO")
ctry_shps =     do.call("bind", lapply(ssa$countries, 
                function(x) raster::getData('GADM', country=x, level=1)))
ctry_outlines = do.call("bind", lapply(ssa$countries, 
                function(x) raster::getData('GADM', country=x, level=0)))

## Countries and regions excluded from analysis 
other = list()
other$countries = c("MAR","DZA","CPV","DJI","EGY","LBY","MUS","REU","STP","SYC","TUN","COM")
other_shps =      do.call("bind", lapply(other$countries, function(x) 
                  raster::getData('GADM', country=x, level=0)))
# COM (Comoros) excluded because 0 links (thereby not usable in CAR model)
# All other countries are generally not considered part of sub-Saharan Africa 

ctry_shps = ctry_shps[ctry_shps$NAME_1!="Bolama",]
ctry_shps = ctry_shps[ctry_shps$NAME_1!="Annobón",]
# Islands with 0 links

######################################################################
# Adding covariate data
######################################################################

setwd("~/desktop/MalDrugRes_SSA/covariates")

## Final variables included in model after removing correlated variables are below
## E.g. excluded housing quality (source=MAP), population count/density (source=GPW Columbia)

## Example of downloading data from MAP 
## Uncomment to run
#shp <- getShp(ISO = ssa$countries, admin_level = "admin0")
#PfPR2_10_2006 <- getRaster(surface = "Plasmodium falciparum PR2-10",shp = shp, year = 2006)
#v_06 <- extract(PfPR2_10_2006, ctry_shps) 
#saveList(object = v_06, file = "pfPR2_10_2006.RData")

## COVAR= Pf-PR_2-10 in 2006 and 2014; SOURCE= Malaria Atlas Project (MAP)
v_06 <- getList('pfPR2_10_2006.RData', compress = TRUE, verbose = FALSE, truncate = FALSE)
  means <- rep(NA, length(v_06))
  for (i in 1:length(v_06)){ means[i] <- mean(v_06[[i]],na.rm=T)}
  ctry_shps$meanPfPRrate_06 <- means
v_14 <- getList('pfPR2_10_2014.RData', compress = TRUE, verbose = FALSE, truncate = FALSE)
  means <- rep(NA, length(v_14))
  for (i in 1:length(v_14)){ means[i] <- mean(v_14[[i]],na.rm=T)}
  ctry_shps$meanPfPRrate_14 <- means

## COVAR= Insecticide-Treated Nets (ITNs) coverage in 2006 and 2014; SOURCE= MAP
itns_06 <- getList('itns_2006.RData', compress = TRUE, verbose = FALSE, truncate = FALSE)
  means <- rep(NA, length(itns_06))
  for (i in 1:length(itns_06)){ means[i] <- mean(itns_06[[i]],na.rm=T)}
  ctry_shps$itns_06 <- means
itns_14 <- getList('itns_2014.RData', compress = TRUE, verbose = FALSE, truncate = FALSE)
  means <- rep(NA, length(itns_14))
  for (i in 1:length(itns_14)){ means[i] <- mean(itns_14[[i]],na.rm=T)}
  ctry_shps$itns_14 <- means

## COVAR= City accessibility; SOURCE= MAP
access <- getList('city_travel.RData', compress = TRUE, verbose = FALSE, truncate = FALSE)
  means <- rep(NA, length(access))
  for (i in 1:length(access)){ means[i] <- mean(access[[i]],na.rm=T)}
  ctry_shps$access <- means
  ctry_shps$log.access <- log(ctry_shps$access)

## COVAR= ACT coverage in 2006 and 2014; SOURCE= MAP
a <- readRDS('ACTcov_2006.RData')
  means <- rep(NA, length(a))
  for (i in 1:length(a)){ means[i] <- mean(a[[i]],na.rm=T)}
  ctry_shps$ACTcov06 <- means
b <- readRDS('ACTcov_2014.RData')
  means <- rep(NA, length(b))
  for (i in 1:length(b)){ means[i] <- mean(b[[i]],na.rm=T)}
  ctry_shps$ACTcov14 <- means

## COVAR= P. vivax in 2006 and 2014; SOURCE= MAP
pv_06 <- getList('pv_06.RData', compress = TRUE, verbose = FALSE, truncate = FALSE)
  pv_06_2 <- vector("list")
  for (i in 1:length(pv_06)) {
          ex <- pv_06[[i]]
          ex[ex==-1] <- 0
          pv_06_2[[i]] <- ex }
  means <- rep(NA, length(pv_06_2))
  for (i in 1:length(pv_06_2)){ means[i] <- mean(pv_06_2[[i]],na.rm=T)}
  ctry_shps$meanPv_06 <- means
pv_14 <- getList('pv_14.RData', compress = TRUE, verbose = FALSE, truncate = FALSE)
  pv_14_2 <- vector("list")
  for (i in 1:length(pv_14)) {
          ex <- pv_14[[i]]
          ex[ex==-1] <- 0
          pv_14_2[[i]] <- ex }
  means <- rep(NA, length(pv_14_2))
  for (i in 1:length(pv_14_2)){ means[i] <- mean(pv_14_2[[i]],na.rm=T)}
  ctry_shps$meanPv_14 <- means

## COVAR= First-line drug policy; SOURCE= WHO
ctry_shps$drug_pol.f <- NA
x <- c("Angola","Benin","Botswana","Central African Republic","Chad","Ethiopia",
       "Gambia","Guinea-Bissau","Kenya","Malawi","Mozambique","Namibia","Niger",
       "Rwanda","South Africa","Swaziland","Uganda",
       "Zambia","Zimbabwe","Tanzania")
y <- c("Burundi","Cameroon","Republic of Congo","Côte d'Ivoire",
       "Democratic Republic of the Congo",
       "Equatorial Guinea","Eritrea","Gabon","Guinea","Liberia","Madagascar",
       "Sierra Leone","South Sudan")
z <- c("Burkina Faso","Mali","Togo","Nigeria","Sudan","Senegal","Ghana",
       "Mauritania","Somalia","Lesotho")
for (i in x) {  ctry_shps$drug_pol.f[which(ctry_shps$NAME_0==i)] <- 1} #AL policy
for (i in y) {  ctry_shps$drug_pol.f[which(ctry_shps$NAME_0==i)] <- 2} #ASAQ policy
for (i in z) {  ctry_shps$drug_pol.f[which(ctry_shps$NAME_0==i)] <- 3} #AL & ASAQ (or other) policy
ctry_shps$drug_pol.f <- as.factor(ctry_shps$drug_pol.f)
ctry_shps$drug_pol.al <- ifelse(ctry_shps$drug_pol.f=="1",1,0)
ctry_shps$drug_pol.asaq <- ifelse(ctry_shps$drug_pol.f=="2",1,0)
ctry_shps$drug_pol.multi <- ifelse(ctry_shps$drug_pol.f=="3",1,0)

## COVAR= Rainfall seasonality; SOURCE= WorldClim
#invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE))
#library(raster)
#climate <- getData('worldclim', var='bio', res=2.5)
#seasonalprecip <- climate$bio15
#seasonalprecip.inter <- extract(seasonalprecip, ctry_shps) 
#saveRDS(object = seasonalprecip.inter, file = 'seasonalprecip.RData')
seasonality <- readRDS('seasonalprecip.RData')
means <- rep(NA, length(seasonality))
for (i in 1:length(seasonality)){ means[i] <- mean(seasonality[[i]],na.rm=T)}
ctry_shps$seasonality <- means

