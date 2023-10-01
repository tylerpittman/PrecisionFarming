## Takes farm combine data and fits spatial quintiles
## Tyler Pittman, 15 September 2018
# Rscript --no-save /Users/tylerpittman/Farm/STEP_11_COMBINE_contour_polygons_2018_15September2018.R

#Sys.setenv("pypath"="/Library/Frameworks/GDAL.framework/Programs");
#Sys.which("gdal_polygonize.py");

setwd("/Users/tylerpittman/Farm/input");

library(TeachingDemos);
library(sgeostat);
#library(car);
library(sp);
library(shapefiles);
library(PBSmapping);
library(mapdata);
library(MASS);
library(maptools);
library(spBayes); 	## important gives Bayesian spatial modeling functions
library(RColorBrewer);
library(classInt);     	# finds class intervals for continuous variables
library(MBA);		# make sure to have `boost headers for C++' installed on linux system	
#library(gpclib);
library(rgeos);
library(raster);
library(fossil);
library(akima);
library(mapproj);
library(rgdal);
library(fields);
gpclibPermit() 		##need to use this to permit open license for unionSpatialPolygons feature
#library(RSurvey);	# requires x11 to be installed in ubuntu
#library(rgdal);
library(spatstat); 	## important gives point process modeling functions for kernel density
library(gstat); 	### Must load this after spatstat or will receive error message!!
library(RMySQL);
#library(spTimer);
#library(imputation);
library(xtable);
library(spam);
library(parallel);	# Does multicore processing when using mclapply WORKS on Linux!
source("shape2poly.R"); # Reads in shape2poly function and others
source("polygonizer.R"); # Reads in polygonizer function

# Set seeding crops below here in field for each year with given input boundary and combine yield data from previous year;
key <- matrix(c( 
900001, "Canary","Canola","Wheat","Lentils","Wheat", 
#900002,"Canary","Lentils","Canary","Lentils","Wheat",
900003,"Wheat","Lentils","Wheat","Lentils","Canola",
900004,"Lentils","Wheat","Lentils","Wheat","Lentils",
900005,"Wheat","Lentils","Canary","Canola","Wheat",
900006,"Lentils","Wheat","Lentils","Wheat","Lentils",
900007,"Wheat","Lentils","Wheat","Lentils","Canola",
900008,"Canola","Wheat","Lentils","Canola","Wheat",
900009,"Canary","Canola","Wheat","Lentils","Canola",
#900010,"Wheat","Canola","Wheat","Lentils","Wheat",
900011,"Lentils","Wheat","Lentils","Canola","Wheat",
900012,"Wheat","Lentils","Wheat","Lentils","Canola",
900013,"Lentils","Wheat","Lentils","Wheat","Lentils",
900014,"Wheat","Lentils","Wheat","Lentils","Canola",
#"900014b","Wheat","Lentils","Wheat","Barley",
900015,"Lentils","Wheat","Lentils","Wheat","Lentils",
900016,"Canola","Wheat","Lentils","Wheat","Lentils",
900017,"Canary","Lentils","Wheat","Canola","Wheat",
900018,"Lentils","Wheat","Lentils","Canola","Wheat",
900020,"Lentils","Wheat","Lentils","Barley","Lentils",
900022,"Lentils","Wheat","Lentils","Canola","Wheat",
900024,"Lentils","Wheat","Canary","Lentils","Wheat",
900025,"Lentils","Wheat","Canary","Lentils","Wheat",
900026,"Wheat","Lentils","Wheat","Canary","Lentils",
900027,"Wheat","Lentils","Wheat","Canary","Lentils",
900028,"Lentils","Wheat","Lentils","Wheat","Lentils",
#800028,"Wheat","Wheat","Lentils","Barley","Barley",
900029,"Lentils","Wheat","Lentils","Wheat","Canola",
900030,"Wheat","Wheat","Lentils","Wheat","Lentils",
900031,"Wheat","Lentils","Wheat","Lentils","Wheat"
), ncol=6, byrow=T);
key <- as.data.frame(key);
colnames(key) <- c("field", "Crop2014", "Crop2015", "Crop2016", "Crop2017", "Crop2018");
key$Crop2014 <- as.character(levels(key$Crop2014))[key$Crop2014];
key$Crop2015 <- as.character(levels(key$Crop2015))[key$Crop2015];
key$Crop2016 <- as.character(levels(key$Crop2016))[key$Crop2016];
key$Crop2017 <- as.character(levels(key$Crop2017))[key$Crop2017];
key$Crop2018 <- as.character(levels(key$Crop2018))[key$Crop2018];

dataColName <- c("DryYield", "HarvestM", "FuelRate", "Elevatio");
combineOperation <- c("yield", "harvestm", "fuelrate", "elevation");
seasonOperation <- c("Harvest Yield", "Harvest Moisture", "Harvest Fuel Rate", "Elevation");
plotOperation <- c("Dry Yield Quintiles (in Bushels/Acre) \n From", "Harvest Moisture Quintiles (in % wet basis) \n From", "Fuel Rate Quintiles (in l/h) \n From", "Elevation Quintiles (in Feet) \n From");
colorPlot <- c("RdYlGn", "BrBG", "RdGy", "PuOr");
operationName <- c("Yield", "Moisture", "FuelRate", "Elevation");

####################### All Looping starts below here ###############################;
counterCombine <- length(combineOperation);
	combineLoop <- function(z) {
	k <- z;
####################### Looping starts below here ###############################;
counterField <- length(key$field);

	fieldLoop <- function(y) {
	j <- y;

#k <- 1;
#j <- 24;
year = 2018;
dc = dataColName[k];
op = combineOperation[k];
so = seasonOperation[k];
po = plotOperation[k];
on = operationName[k];
cp = colorPlot[k];
fb = key$field[j];
fy = key$field[j];
key$Crop <- eval(parse(text=paste("key$Crop", year, sep="")));
type = key$Crop[j];
SP_ID = key$field[j];
#fb=900010;
#fy=900010;
#type = "Lentils";
setwd(paste("/Users/tylerpittman/Farm/input/", year, "/boundaries", sep=""));
fieldBoundary <- readShapeSpatial(paste(fb, "_", year, "_boundary", ".shp", sep=""));
fieldBoundary.dbf <- read.dbf(paste(fb, "_", year, "_boundary", ".dbf", sep=""), header=TRUE);
#fieldYield <- readShapeSpatial(paste(fy, "_2014_yield", ".shp", sep=""));
setwd(paste("/Users/tylerpittman/Farm/input/", year, "/yield", sep=""));
fieldYield <- readOGR(dsn=".", layer=paste(fy, "_", year, "_yield", sep=""));
fieldYield.dbf <- read.dbf(paste(fy, "_", year, "_yield", ".dbf", sep=""), header=TRUE);
#10/5/2014 7:36:51 PM
timeFull <- strptime(fieldYield@data$Time, format="%m/%d/%Y %I:%M:%S %p");
CombDate <- median(as.Date(timeFull));
CombDate <- format(CombDate, format="%B %d");

proj4string(fieldBoundary) <- CRS("+proj=utm +zone=13 +ellps=WGS84")
fieldBoundary <- spTransform(fieldBoundary, CRS("+proj=longlat"))

### This fixes invalid polygon boundaries;
#library(sf);
#st_is_valid(fieldBoundary, reason = TRUE);
#st_make_valid(fieldBoundary);
###
fieldYield.simp <- gSimplify(fieldYield, tol = 0.000001, topologyPreserve=FALSE);
fieldYield.simp <- gSimplify(fieldYield.simp, tol = 0.000001, topologyPreserve=FALSE);
fieldBoundary.simp <- gSimplify(fieldBoundary, tol = 0.000001, topologyPreserve=FALSE);
fieldBoundary.simp <- gSimplify(fieldBoundary.simp, tol = 0.000001, topologyPreserve=FALSE);
#gIsValid(fieldYield, reason = TRUE);
#gIsValid(fieldBoundary, reason = TRUE);
#gIsValid(fieldYield.simp, reason = TRUE);
#gIsValid(fieldBoundary.simp, reason = TRUE);
fieldYield.clip <-  intersect(fieldYield.simp, fieldBoundary.simp);
fieldYield@data$ID <- paste(row.names(fieldYield@data)); 
#fieldYield@data$ID <- paste(row.names(fieldYield@data), "1"); 
IDinput <- row.names(fieldYield.simp);
IDcut <- row.names(fieldYield.clip)
fieldYield.clip.df <- data.frame(fieldYield[match(row.names(fieldYield.clip), fieldYield$ID),], row.names=IDcut);
fieldYield.clip.df$ID <- IDcut;
row.names(fieldYield.clip) <- IDcut;
fieldYield.clip2 <- SpatialPointsDataFrame(fieldYield.clip, data=fieldYield.clip.df);
fieldYield <- fieldYield.clip2;
fieldBoundary <- fieldBoundary.simp;

#fieldYield@coords;
#str(fieldYield@data);
#fieldYield@data$DryYield;

#newfieldYield <- aggregate(fieldYield@data$DryYield, by =  list(ID = rep(1, length(fieldYield))), FUN = mean, dissolve = TRUE, areaWeighted = TRUE);


#--------------------## This is where the grid interpolation section starts ##--------------------#; 

###### WORKING! YieldContour Projected ######
projected <- "+proj=utm +zone=13 +ellps=WGS84";
o2 <- cbind(fieldYield@coords, fieldYield@data[which(colnames(fieldYield@data) == dc)]);
colnames(o2) <- c("lon", "lat", "DryYield");
o2 <- as.data.frame(o2);
#o2 <- o2[complete.cases(o2),];

#set DryYield value of 0.000 to mean of non-zero points (doesn't work for quantiles) or delete
zeroList <- which(o2$DryYield <= 0.1);
if(length(zeroList) != 0) {
	o2nz <- o2[-c(zeroList), ];
	o2 <- o2nz;
};
meanYield <- mean(o2$DryYield);
medianYield <- median(o2$DryYield);

#zeroList <- which(o2$DryYield <= 1.0);
#if(length(zeroList) != 0) {
#  for(i in 1:length(zeroList)){
#	 o2[zeroList[i], ]$DryYield = runif(1, min=0.1, max=1.0);
#  };
#};

#o2[which(o2$DryYield == 0.000), ]$DryYield = runif(1, min=0.001, max=0.009); #recode to random value so interval breaks are unique
#o2 <- o2[-c(which(o2$DryYield == 0.000)), ];
#o2 <- o2[c(which(o2$DryYield == 0.000)), ];

nclrYield <- 5; 
#display.brewer.all(n=NULL, type="div", select=NULL, exact.n=TRUE);
plotclrYield <- brewer.pal(nclrYield, cp);
length(plotclrYield);
brksYield <- classIntervals(o2$DryYield, 5, style = "quantile", rtimes = 3, intervalClosure = "right");
brksRounded <- round(brksYield$brks, 3);
brksYield <- classIntervals(o2$DryYield, 5, style = "fixed", fixedBreaks = brksRounded, rtimes = 3, intervalClosure = "right");
colcodeYield <- findColours(brksYield, plotclrYield)
brksYield <- brksYield$brks;


#########
##add 400 (cDim) points of mean intensity (one in each cell) to input file and then append to o2;
#uniformPoints <- expand.grid(seq(min(o2$lon),max(o2$lon),length=20),seq(min(o2$lat),max(o2$lat),length=20));
#uniformPoints$DryYield <- -9999; #add 0.00 so meanYield not on boundary for missing areas;
#uniformPoints <- as.data.frame(uniformPoints);
#colnames(uniformPoints) <- colnames(o2);
##plot(uniformPoints);
#o2 <- rbind(o2, uniformPoints);
##########

coordinates(o2) <- c(1, 2)
proj4string(o2) <- CRS("+proj=longlat")
o2new <- spTransform(o2, CRS(projected))
proj4string(fieldBoundary) <- CRS("+proj=longlat")
fieldBoundaryNew <- spTransform(fieldBoundary, CRS(projected))

fieldArea <- gArea(fieldBoundaryNew)/4046.86;

cXDim = ((max(o2new$lon) - min(o2new$lon)))/2;
cYDim = ((max(o2new$lat) - min(o2new$lat)))/2;
#cDim = 200;
#cSLon = ((max(o2new$lon) - min(o2new$lon)))/cDim;
#cSLat = ((max(o2new$lat) - min(o2new$lat)))/cDim;
#cSLonB = ((max(fieldBoundaryNew@bbox[2,]) - min(fieldBoundaryNew@bbox[2,])))/cDim;
#cSLatB = ((max(fieldBoundaryNew@bbox[1,]) - min(fieldBoundaryNew@bbox[1,])))/cDim;

##below could be factor of 1.000005 or 0.999995 depending on lat/lon or projected coordinate system;
gridT = list(cellcentre.offset=c((min(o2new$lon)),(min(o2new$lat))), cellsize=c(2,2), cells.dim=c(cXDim,cYDim));
#gridT = list(cellcentre.offset=c((min(o2new$lon)),(min(o2new$lat))), cellsize=c(cSLon,cSLat), cells.dim=c(cDim,cDim));
#gridT = list(cellcentre.offset=c(min(o2new$lon),min(o2new$lat)), cellsize=c(cSLon,cSLat), cells.dim=c(cDim,cDim));
grd2 <- GridTopology(gridT$cellcentre.offset, gridT$cellsize, gridT$cells.dim)
SG2 <- SpatialGrid(grd2, proj4string=CRS(projected))
#check to see if grid covers field;
#plot(o2new, col=colcodeYield)
#plot(SG2, col="Blue", add=T)
p2 <- idw(DryYield ~ 1, o2new, newdata=SG2, nmax=100, maxdist=5.6, idp=0.6)

#p2 <- idw(DryYield ~ 1, o2new, newdata=SG2, maxdist=5, nmax=4);
#p2 <- idw(DryYield ~ 1, o2new, newdata=SG2, na.action = na.omit, nmax=4, maxdist=16.5, idp=2); 

#check to see if grid covers field;
#plot(o2new, col=colcodeYield)
#plot(p2, col="Blue", add=T)

### This works for cutting edges! Choose p2B if large empty yield inside of boundary ###
t2 <- p2;
#fullgrid(t2) = FALSE;
fullgrid(t2) = TRUE; #if p2B
#t2.clip = t2[!is.na(overlay(t2, fieldBoundaryNew)),];
t2.clip <- t2;
pim2 <- as.image.SpatialGridDataFrame(t2.clip)
cl2 <- contourLines(pim2)
cL2 <- ContourLines2SLDF(cl2)
proj4string(cL2) <- CRS(projected)
#cL2clip <- gIntersection(fieldBoundary, cL2, byid=TRUE, id=NULL); #this is from rgeos package;

names(attr(colcodeYield, "table"));

t2.clip$Yield <- cut(t2.clip$var1.pred, breaks=brksYield)
r <- raster(t2.clip, layer=3); #3 corresponds to the yield column of the t2.clip@data slot we want
#res(r); #cell grid cellsize from above;
#r.poly <- rasterToPolygons(r, dissolve=TRUE)
r.poly <- polygonizer(r);
r.poly <- unionSpatialPolygons(r.poly, r.poly$DN);
proj4string(r.poly) <- CRS(projected);
#gIsValid(r.poly, reason = TRUE);
r.poly <- gSimplify(r.poly, tol = 0.1, topologyPreserve=FALSE);
r.poly.df <- data.frame(c(1:5));
colnames(r.poly.df) <- c("Yield");
r.poly <- SpatialPolygonsDataFrame(r.poly, data=r.poly.df);
#plot(r.poly, col=plotclrYield, border=T);

#r.poly@data$ADSFLDID = fb;
r.poly@data$SEASON = so;
r.poly@data$c_year = year;
r.poly@data$c_type = type;
r.poly@data <- cbind(SP_ID, r.poly@data);

combinedArea <- gArea(r.poly)/4046.86;
#meanYield <- mean(o2[-c(which(o2$DryYield <= 0.1)), ]$DryYield);
r.poly@data$Field_Area <- fieldArea;
r.poly@data$Field_Combined_Area <- combinedArea;
r.poly@data$Field_Mean_Yield <- meanYield;

COACH <- r.poly@data$Yield
CROP <- r.poly@data$Yield
LABEL <- names(attr(colcodeYield, "table"))
Quintile_Area <- r.poly@data$Yield
Quintile_Mean_Yield <- r.poly@data$Yield
options(scipen=3);

for(i in 1:length(r.poly@data$Yield)){
	if (r.poly@data$Yield[[i]] == 1) {COACH[[i]] = "s2"} else	
	if (r.poly@data$Yield[[i]] == 2) {COACH[[i]] = "s1"} else
	if (r.poly@data$Yield[[i]] == 3) {COACH[[i]] = "normal"} else
	if (r.poly@data$Yield[[i]] == 4) {COACH[[i]] = "a1"} else
	if (r.poly@data$Yield[[i]] == 5) {COACH[[i]] = "a2"} else
	COACH[[i]] = NA;
	CROP[[i]] = type;
	LABEL[[i]] = names(attr(colcodeYield, "table"))[i];
	Quintile_Area[[i]] <- gArea(r.poly[i,])/4046.86;
	Quintile_Mean_Yield[[i]] <- mean(o2[c(which(o2$DryYield > brksYield[i] & o2$DryYield <= brksYield[i+1])), ]$DryYield);
};
r.poly@data$COACH <- COACH;
r.poly@data$LABEL <- LABEL;
r.poly@data$Quintile_Area <- Quintile_Area;
r.poly@data$Quintile_Mean_Yield <- Quintile_Mean_Yield;

r.poly@data$Field_Combined_Area <- format(r.poly@data$Field_Combined_Area, digits=5);
r.poly@data$Field_Area <- format(r.poly@data$Field_Area, digits=5);
r.poly@data$Field_Mean_Yield <- format(r.poly@data$Field_Mean_Yield, digits=5);
r.poly@data$Quintile_Area <- format(r.poly@data$Quintile_Area, digits=4);
r.poly@data$Quintile_Mean_Yield <- format(r.poly@data$Quintile_Mean_Yield, digits=2);
r.poly@data$CombDate <- CombDate;
names(r.poly) <- sub("^Yield$", on, names(r.poly));

r.poly.clip <-  intersect(r.poly, fieldBoundaryNew);
r.poly  <- r.poly.clip;

###
areaBoundary <- gArea(fieldBoundaryNew)/4046.86;
areaPolys <- as.numeric(gArea(r.poly, byid = TRUE)/4046.86);
areaCombined <- sum(areaPolys);
areaBoundary <- round(areaBoundary, 3);
areaPolys <- round(areaPolys, 3);
areaCombined <- round(areaCombined, 3);
r.poly@data$Field_Combined_Area <- areaCombined;
r.poly@data$Quintile_Area <- areaPolys;
r.poly@data$Field_Area <- areaBoundary;
###

#plot(r.poly)
plot(r.poly, col=plotclrYield, border=T);
#plot(r.poly[r.poly$Yield == 2,], col="red")
#gArea(r.poly[which(r.poly$COACH == "s2"),])/4046.86
#plot(r.poly[which(r.poly$COACH == "s2"),]);
#mean(o2[c(which(o2$DryYield > brksYield[1] & o2$DryYield <= brksYield[2])), ]$DryYield);

writePolyShape(r.poly, paste("/Users/tylerpittman/Farm/", op, "Shapefiles/", year, "/", fb, "_", op, "_proj_", year, ".shp", sep=''));


#### Below is to replace boundary attributes, not always needed ###
#fieldBoundaryNew@data$c_year = year;
#fieldBoundaryNew@data$c_type = type;
#fieldBoundaryNew@data$Field_Area <- fieldArea;
#fieldBoundaryNew@data$Field_Area <- format(fieldBoundaryNew@data$Field_Area, digits=5);
#writePolyShape(fieldBoundaryNew, paste("/Users/tylerpittman/Farm/input/2015/boundaries/", fb, "_2015_boundary", ".shp", sep=''));

legendMean <- paste("Mean ", so, ":", sep="");

png(paste("/Users/tylerpittman/Farm/", op, "Shapefiles/", year, "/images/", fb, "_", op, "_proj_", year, ".png", sep=''), res=100, pointsize = 12, units = "px", bg = "transparent");
par(xpd=F, mar=c(0,0,0,1), oma=c(0,0,2,1));
#image(t2.clip, col = plotclrYield, add = FALSE, breaks=brksYield, oldstyle = FALSE, useRaster=TRUE);
#contour(t2.clip, lwd=0.5, col=plotclrYield, levels=brksYield, drawlabels=F, add=T);
#plot(fieldBoundaryNew, lwd=2.0, add=T);
plot(r.poly, col=plotclrYield, border=F);
title(paste(fb, po, "Harvest", year), outer=TRUE);
par(lwd=2); legend("bottomright",  box.lwd=1, legend=rev(names(attr(colcodeYield, "table"))), fill=rev(plotclrYield), cex=0.9, bg="transparent", xjust=0);
par(lwd=2); legend("bottomleft",  box.lwd=1, legend=cbind(c("Field Area:", fieldArea), c("Combined Area:", combinedArea), c(legendMean, meanYield)), cex=0.9, bg="transparent", xjust=0);
#par(lwd=2); legend("bottomright",  box.lwd=1, legend=rev(brksYield[-1]), fill=rev(plotclrYield), cex=0.9, bg="transparent", xjust=0);
dev.off();

#--------------------## This is where the grid interpolation section ends ##--------------------#; 


print(j);
print(fb);
j <- j + 1;

}
mclapply(1:counterField, fieldLoop);
#mclapply(22:counterField, fieldLoop);
#mclapply(21, fieldLoop);
####################### Looping ends here ###############################;
print(k);
k <- k + 1;

}
mclapply(1:counterCombine, combineLoop);
#mclapply(1, combineLoop);
####################### All Looping ends here ###############################;

#stopCluster(cl);
