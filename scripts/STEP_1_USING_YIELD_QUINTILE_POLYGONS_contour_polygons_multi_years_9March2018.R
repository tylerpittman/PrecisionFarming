## Takes farm yield data over all years and fits spatial quintiles
## Tyler Pittman, 9 March 2018
# Rscript --no-save /Users/tylerpittman/Farm/STEP_1_USING_YIELD_QUINTILE_POLYGONS_contour_polygons_multi_years_9March2018.R

Sys.setenv(PATH="/usr/bin:/opt/local/bin:/opt/local/sbin:/Library/Frameworks/GDAL.framework/Programs:/usr/local/mysql/bin:/opt/ImageMagick/bin:/Library/Frameworks/R.framework/Versions/3.2/Resources/bin");
Sys.getenv("PATH");
Sys.which("gdal_polygonize.py");
sessionInfo();

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
#library(imputation); #don't use this anymore, use MICE package instead;
library(xtable);
library(spam);
library(rgdal);
library(parallel);	# Does multicore processing when using mclapply WORKS on Linux!
source("shape2poly.R"); # Reads in shape2poly function and others
source("polygonizer.R"); # Reads in polygonizer function

# Set seeding crops below here in field for each year with given input boundary and combine yield data from previous year;
key <- matrix(c( 
900001, "Canary","Canola","Wheat","Lentils","Wheat", 
900002,"Canary","Lentils","Canary","Lentils","Wheat",
900003,"Wheat","Lentils","Wheat","Lentils","Canola",
900004,"Lentils","Wheat","Lentils","Wheat","Lentils",
900005,"Wheat","Lentils","Canary","Canola","Wheat",
900006,"Lentils","Wheat","Lentils","Wheat","Lentils",
900007,"Wheat","Lentils","Wheat","Lentils","Canola",
900008,"Canola","Wheat","Lentils","Canola","Wheat",
900009,"Canary","Canola","Wheat","Lentils","Canola",
900010,"Wheat","Canola","Wheat","Lentils","Wheat",
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
#"WT1","Wheat","na","Wheat","na","Wheat",
#"WT2","na","Wheat","na","Wheat","na"
), ncol=6, byrow=T);
key <- as.data.frame(key);
colnames(key) <- c("field", "Crop2014", "Crop2015", "Crop2016", "Crop2017", "Crop2018");
key$Crop2014 <- as.character(levels(key$Crop2014))[key$Crop2014];
key$Crop2015 <- as.character(levels(key$Crop2015))[key$Crop2015];
key$Crop2016 <- as.character(levels(key$Crop2016))[key$Crop2016];
key$Crop2017 <- as.character(levels(key$Crop2017))[key$Crop2017];
key$Crop2018 <- as.character(levels(key$Crop2018))[key$Crop2018];

####################### All Looping starts below here ###############################;

####################### Looping starts below here ###############################;
#counterField <- 1;
counterField <- length(key$field);
	fieldLoop <- function(y) {
	j <- y;
	#j <- 10;	
fb = key$field[j];
fy = key$field[j];
#yearStartHarvest = 2014;
if (fb == "900030"){
		yearStartHarvest = 2015;
} else if (fb == "900031"){
		yearStartHarvest = 2015;
} else {
		yearStartHarvest = 2014; 
};
yearCurrent = 2018;
yearBoundary = 2017;
key$Crop <- eval(parse(text=paste("key$Crop", yearCurrent, sep="")));
type = key$Crop[j];

#---(*)---#---(*)---#---(*)---#---(*)---#---(*)---#---(*)---#---(*)---#---(*)---#
#---(*)---#---(*)---#---(*)---#---(*)---#---(*)---#---(*)---#---(*)---#---(*)---#
makeYieldContourYears <- function(yearStartHarvest, yearBoundary, yearCurrent, fy) {

#44.092/36.744=1.199978 is the conversion factor of (canola or canary)/(wheat or lentils);
#1/44.092 = 0.02267985; #1/36.744 = 0.02721533 is tonne/bushels if one wants tonnes of products produced from bushel weights;
canaryW = 2.2;
wheatW = 1;
barleyW = 0.8;
lentilW = 2.2;
canolaW = 1.5;
#canaryW = 1/44.092;
#wheatW = 1/36.744;
#barleyW = 1/45.930;
#lentilW = 1/36.744;
#canolaW = 1/44.092;

#fb="900010b";
#fy="900010b";
#type = "Lentils";
setwd(paste("/Users/tylerpittman/Farm/boundariesSeedingCleanTopcon", yearBoundary, sep=""));
#setwd("/Users/tylerpittman/Farm/boundariesSeedingClean2015");
fieldBoundary <- readShapeSpatial(paste(fb, "", ".shp", sep=""));
fieldBoundary.dbf <- read.dbf(paste(fb, "", ".dbf", sep=""), header=TRUE);
#fieldYield <- readShapeSpatial(paste(fy, "_2015_yield", ".shp", sep=""));
setwd(paste("/Users/tylerpittman/Farm/input/", yearStartHarvest, "/yield", sep=""));
fieldYield <- readOGR(dsn=".", layer=paste(fy, "_", yearStartHarvest, "_yield", sep=""));
fieldYield.dbf <- read.dbf(paste(fy, "_", yearStartHarvest, "_yield", ".dbf", sep=""), header=TRUE);
if (yearStartHarvest == "2014"){
	proj4string(fieldYield) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0");
	fieldYield <- spTransform(fieldYield, CRS("+proj=longlat"));
	fieldYield  <- fieldYield[ , c(which(names(fieldYield) %in% c("Elevatio","DryYield")))];
	kCrop <- eval(parse(text=paste("key$Crop", yearStartHarvest, sep="")));
	kType = kCrop[j];
	if (kType == "Canary"){
		fieldYield$DryYield <- fieldYield$DryYield*canaryW;
	} else if (kType == "Wheat"){
		fieldYield$DryYield <- fieldYield$DryYield*wheatW;
	} else if (kType == "Barley"){
		fieldYield$DryYield <- fieldYield$DryYield*barleyW;
	} else if (kType == "Lentils"){
		fieldYield$DryYield <- fieldYield$DryYield*lentilW; 
	} else {
		fieldYield$DryYield <- fieldYield$DryYield*canolaW;
	};
} else {
	proj4string(fieldYield) <- CRS("+proj=longlat");
	fieldYield  <- fieldYield[ , c(which(names(fieldYield) %in% c("Elevatio","DryYield")))];
	kCrop <- eval(parse(text=paste("key$Crop", yearStartHarvest, sep="")));
	kType = kCrop[j];
	if (kType == "Canary"){
		fieldYield$DryYield <- fieldYield$DryYield*canaryW
	} else if (kType == "Wheat"){
		fieldYield$DryYield <- fieldYield$DryYield*wheatW;
	} else if (kType == "Barley"){
		fieldYield$DryYield <- fieldYield$DryYield*barleyW;
	} else if (kType == "Lentils"){
		fieldYield$DryYield <- fieldYield$DryYield*lentilW; 
	} else {
		fieldYield$DryYield <- fieldYield$DryYield*canolaW;
	};
};
####
#setwd(paste("/Users/tylerpittman/Farm/input/", 2014, "/yield", sep=""));
#fieldYield <- readOGR(dsn=".", layer=paste(900002, "_", 2014, "_yield", sep=""));
#fieldYield.dbf <- read.dbf(paste(900002, "_", 2014, "_yield", ".dbf", sep=""), header=TRUE);
#fieldYield1 <- readOGR(dsn=".", layer=paste(800002, "_", 2014, "_yield", sep=""));
#fieldYield1.dbf <- read.dbf(paste(800002, "_", 2014, "_yield", ".dbf", sep=""), header=TRUE);
#fieldYieldJ <- spRbind(fieldYield, fieldYield1);
#writeSpatialShape(fieldYieldJ, "900002_2014_yield.shp");
#####
year = yearStartHarvest;
#if (yearStartHarvest < (yearCurrent-1)){
if (yearStartHarvest < (yearStartHarvest)){
for (year in ((yearStartHarvest+1):(yearCurrent-1))) {	
	setwd(paste("/Users/tylerpittman/Farm/input/", year, "/yield", sep=""));
	temp.data <- readOGR(dsn=".", layer=paste(fy, "_", year, "_yield", sep=""));
	temp.data.dbf <- read.dbf(paste(fy, "_", year, "_yield", ".dbf", sep=""), header=TRUE);
	
	if (year == "2014"){
		proj4string(temp.data) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0");
		temp.data <- spTransform(fieldYield, CRS("+proj=longlat"));
		temp.data  <- temp.data[ , c(which(names(temp.data) %in% c("Elevatio","DryYield")))];
		kCrop <- eval(parse(text=paste("key$Crop", year, sep="")));
		kType = kCrop[j];
		if (kType == "Canary"){
			temp.data$DryYield <- temp.data$DryYield*canaryW;
		} else if (kType == "Wheat"){
			temp.data$DryYield <- temp.data$DryYield*wheatW;
		} else if (kType == "Barley"){
			temp.data$DryYield <- temp.data$DryYield*barleyW;
		} else if (kType == "Lentils"){
			temp.data$DryYield <- temp.data$DryYield*lentilW; 
		} else {
			temp.data$DryYield <- temp.data$DryYield*canolaW;
		};
		} else {
		proj4string(temp.data) <- CRS("+proj=longlat");
		temp.data  <- temp.data[ , c(which(names(temp.data) %in% c("Elevatio","DryYield")))];
		kCrop <- eval(parse(text=paste("key$Crop", year, sep="")));
		kType = kCrop[j];
		if (kType == "Canary"){
			temp.data$DryYield <- temp.data$DryYield*canaryW;
		} else if (kType == "Wheat"){
			temp.data$DryYield <- temp.data$DryYield*wheatW;
		} else if (kType == "Barley"){
			temp.data$DryYield <- temp.data$DryYield*barleyW;
		} else if (kType == "Lentils"){
			temp.data$DryYield <- temp.data$DryYield*lentilW; 
		} else {
			temp.data$DryYield <- temp.data$DryYield*canolaW;
		};
	};
	fieldYield <- spRbind(fieldYield,temp.data);
}	
};	
#^^^^^#^^^^^#^^^^^#^^^^^#^^^^^#^^^^^#^^^^^#^^^^^#^^^^^#^^^^^#
#fb="900010a";
#fy="900010a";
#fieldYield <- readOGR(dsn=".", layer=paste(fy, "_", year-1, "_yield", sep=""));
#fieldYield.dbf <- read.dbf(paste(fy, "_", year-1, "_yield", ".dbf", sep=""), header=TRUE);
#fb="900010b";
#fy="900010b";
#fieldYield1 <- readOGR(dsn=".", layer=paste(fy, "_", year-1, "_yield", sep=""));
#fieldYield1.dbf <- read.dbf(paste(fy, "_", year-1, "_yield", ".dbf", sep=""), header=TRUE);
#fyJoined<- spRbind(fieldYield, fieldYield1);
#writeSpatialShape(fyJoined, "900010_2015_yield.shp");
#^^^^^#^^^^^#^^^^^#^^^^^#^^^^^#^^^^^#^^^^^#^^^^^#^^^^^#^^^^^#
#fieldYield@coords;
#str(fieldYield@data);
#fieldYield@data$DryYield;

proj4string(fieldBoundary) <- CRS("+proj=longlat");
#fieldYield.simp <- gSimplify(fieldYield, tol = 0.000001, topologyPreserve=FALSE);
#fieldYield.simp <- gSimplify(fieldYield.simp, tol = 0.000001, topologyPreserve=FALSE);
fieldYield.simp <- fieldYield;
fieldBoundary.simp <- gSimplify(fieldBoundary, tol = 0.000001, topologyPreserve=FALSE);
fieldBoundary.simp <- gSimplify(fieldBoundary.simp, tol = 0.000001, topologyPreserve=FALSE);
#gIsValid(fieldYield, reason = TRUE);
#gIsValid(fieldBoundary, reason = TRUE);
#gIsValid(fieldYield.simp, reason = TRUE);
#gIsValid(fieldBoundary.simp, reason = TRUE);
fieldYield.clip <-  intersect(fieldYield.simp, fieldBoundary.simp);
fieldYield@data$ID <- paste(row.names(fieldYield@data)); 
IDinput <- row.names(fieldYield.simp);
IDcut <- row.names(fieldYield.clip)
fieldYield.clip.df <- data.frame(fieldYield[match(row.names(fieldYield.clip), fieldYield$ID),], row.names=IDcut);
fieldYield.clip.df$ID <- IDcut;
row.names(fieldYield.clip) <- IDcut;
fieldYield.clip2 <- SpatialPointsDataFrame(fieldYield.clip, data=fieldYield.clip.df);
fieldYield <- fieldYield.clip2;
fieldBoundary.df <- data.frame(fieldBoundary);
fieldBoundary <- SpatialPolygonsDataFrame(fieldBoundary.simp, data=fieldBoundary.df);

#--------------------## This is where the grid interpolation section starts ##--------------------#; 

###### WORKING! YieldContour Projected ######
projected <- "+proj=utm +zone=13 +ellps=WGS84";
o2 <- cbind(fieldYield@coords, fieldYield@data$DryYield);
colnames(o2) <- c("lon", "lat", "DryYield");
o2 <- as.data.frame(o2);
o2 <- o2[complete.cases(o2),];

#set DryYield value of 0.000 to mean of non-zero points (doesn't work for quantiles) or delete
o2nz <- o2[-c(which(o2$DryYield == 0.000)), ];
meanYield <- mean(o2nz$DryYield)
medianYield <- median(o2nz$DryYield)
#o2[which(o2$DryYield == 0.000), ]$DryYield = meanYield;
o2 <- o2[-c(which(o2$DryYield == 0.000)), ];

########
##add 400 (cDim) points of mean intensity (one in each cell) to input file and then append to o2;
#uniformPoints <- expand.grid(seq(min(o2$lon),max(o2$lon),length=20),seq(min(o2$lat),max(o2$lat),length=20));
#uniformPoints$DryYield <- medianYield + 0.00; #add 0.00 so meanYield not on boundary for missing areas;
#uniformPoints <- as.data.frame(uniformPoints);
#colnames(uniformPoints) <- colnames(o2);
##plot(uniformPoints);
#o2 <- rbind(o2, uniformPoints);
#########

coordinates(o2) <- c(1, 2)
proj4string(o2) <- CRS("+proj=longlat")
o2new <- spTransform(o2, CRS(projected))
proj4string(fieldBoundary) <- CRS("+proj=longlat")
fieldBoundaryNew <- spTransform(fieldBoundary, CRS(projected))

nclrYield <- 5; 
#display.brewer.all(n=NULL, type="div", select=NULL, exact.n=TRUE);
plotclrYield <- brewer.pal(nclrYield,"RdYlGn");
length(plotclrYield);
brksYield <- classIntervals(o2new$DryYield, 5, style = "quantile", rtimes = 3, intervalClosure = c("left", "right"), dataPrecision = NULL);
colcodeYield <- findColours(brksYield, plotclrYield)
brksYield <- brksYield$brks;

fieldArea <- gArea(fieldBoundaryNew)/4046.86;

#cXDim = ((max(o2new$lon) - min(o2new$lon)))/2;
#cYDim = ((max(o2new$lat) - min(o2new$lat)))/2;
## Must use fieldBounaryNew below as there may be missing yield data within boundaries;
cXDim = ((extent(fieldBoundaryNew)@xmax - extent(fieldBoundaryNew)@xmin))/2;
cYDim = ((extent(fieldBoundaryNew)@ymax - extent(fieldBoundaryNew)@ymin))/2;

##below could be factor of 1.000005 or 0.999995 depending on lat/lon or projected coordinate system;
#gridT = list(cellcentre.offset=c((min(o2new$lon)),(min(o2new$lat))), cellsize=c(2,2), cells.dim=c(cXDim,cYDim));
gridT = list(cellcentre.offset=c((extent(fieldBoundaryNew)@xmin),(extent(fieldBoundaryNew)@ymin)), cellsize=c(2,2), cells.dim=c(cXDim,cYDim));
grd2 <- GridTopology(gridT$cellcentre.offset, gridT$cellsize, gridT$cells.dim)
SG2 <- SpatialGrid(grd2, proj4string=CRS(projected))
#check to see if grid covers field;
#plot(o2new, col=colcodeYield)
#plot(SG2, col="Blue", add=T)
p2 <- idw(DryYield ~ 1, o2new, newdata=SG2, nmax=1000, block = c(1000,1000))

########
#add 400 (cDim) points of mean intensity (one in each cell) to input file and then append to o2;
medianYield <- median(p2$var1.pred);
#DryYield <- p2$var1.pred;
#lon <- coordinates(p2)[,1];
#lat <- coordinates(p2)[,2];
DryYield <- o2new$DryYield;
lon <- o2new$lon;
lat <- o2new$lat;
o2n <- cbind(lon, lat, DryYield);
colnames(o2n) <- c("lon", "lat", "DryYield");
o2n <- as.data.frame(o2n);
uniformPoints <- expand.grid(seq(min(o2n$lon),max(o2n$lon),length=80),seq(min(o2n$lat),max(o2n$lat),length=80));
uniformPoints$DryYield <- medianYield + 0.00; #add 0.00 so meanYield not on boundary for missing areas;
uniformPoints <- as.data.frame(uniformPoints);
colnames(uniformPoints) <- colnames(o2n);
#plot(uniformPoints);
o2n <- rbind(o2n, uniformPoints);
colnames(o2n) <- c("lon", "lat", "DryYield");
coordinates(o2n) <- c(1, 2);
proj4string(o2n) <- CRS(projected);
colcodeYieldn <- vector(mode="character", length=length(uniformPoints[,1]));
colcodeYieldn[] <- plotclrYield[3];
colcodeYield2 <- c(colcodeYield, colcodeYieldn);
#plot(o2n, col=colcodeYield2);
p2 <- idw(DryYield ~ 1, o2n, newdata=SG2, nmax=1000, block = c(1000,1000));
#m <- vgm(.59, "Sph", 874, .04);
#p2 <- krige(DryYield~1, o2n, SG2, nmax=100, block = c(1000,1000));
#########


brksYield1 <- classIntervals(p2$var1.pred, 5, style = "quantile", rtimes = 3, intervalClosure = c("left", "right"), dataPrecision = NULL);
colcodeYield1 <- findColours(brksYield1, plotclrYield)
brksYield1 <- brksYield1$brks;
colcodeYield <- colcodeYield1;
brksYield <- brksYield1;
#image(p2);
#check to see if grid covers field;
#plot(o2n, col=colcodeYield2);
#identify(o2new$DryYield, labels=row.names(o2new), plot=TRUE);
#plot(p2, col="Blue", add=T)

### This works for cutting edges! Choose p2B if large empty yield inside of boundary ###
t2 <- p2;
fullgrid(t2) = TRUE; #if p2B
t2.clip <- t2;
pim2 <- as.image.SpatialGridDataFrame(t2.clip)
cl2 <- contourLines(pim2)
cL2 <- ContourLines2SLDF(cl2)
proj4string(cL2) <- CRS(projected)
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
r.poly@data$SEASON = "Seeding Rx";
r.poly@data <- cbind(r.poly@data, fieldBoundary@data);
COACH <- r.poly@data$Yield
CROP <- r.poly@data$Yield
for(i in 1:length(r.poly@data$Yield)){
	if (r.poly@data$Yield[[i]] == 1) {COACH[[i]] = "s2"} else	
	if (r.poly@data$Yield[[i]] == 2) {COACH[[i]] = "s1"} else
	if (r.poly@data$Yield[[i]] == 3) {COACH[[i]] = "normal"} else
	if (r.poly@data$Yield[[i]] == 4) {COACH[[i]] = "a1"} else
	if (r.poly@data$Yield[[i]] == 5) {COACH[[i]] = "a2"} else
	COACH[[i]] = NA;
	CROP[[i]] = type;
};
r.poly@data$COACH <- COACH;
r.poly@data$CROP <- CROP;

##############################
### simplifies polygons ###
r.poly.simp <- gSimplify(r.poly, tol = 0.000001, topologyPreserve=FALSE);
r.poly.simp <- gSimplify(r.poly.simp, tol = 0.000001, topologyPreserve=FALSE);
#gIsValid(r.poly, reason = TRUE);
#gIsValid(r.poly.simp, reason = TRUE);
r.poly@data$ID <- paste(row.names(r.poly@data)); 
IDinput <- row.names(r.poly);
IDcut <- row.names(r.poly.simp);
r.poly.simp.df <- data.frame(r.poly[match(row.names(r.poly.simp), r.poly$ID),], row.names=IDcut);
r.poly.simp.df$ID <- IDcut;
row.names(r.poly.simp) <- IDcut;
r.poly.simp2 <- SpatialPolygonsDataFrame(r.poly.simp, data=r.poly.simp.df);
r.poly <- r.poly.simp2;
##############################

r.poly.clip <-  intersect(r.poly, fieldBoundaryNew);
r.poly  <- r.poly.clip;
r.poly <- r.poly[ , c(which(names(r.poly) %in% c("Yield", "SEASON", "COACH", "CROP", "ID")))];

#plot(r.poly)
#plot(r.poly, col=plotclrYield, border=F)
writePolyShape(r.poly, paste("/Users/tylerpittman/Farm/contour_yields_proj/quintiles", yearCurrent, "fromYears/", fb, "_SeedingRx_proj_", yearCurrent, "from", yearStartHarvest, "Yield", ".shp", sep=''));

eval(parse(text=paste("r.poly", yearStartHarvest, "<-", "r.poly", sep="")));

###### WORKING! YieldContour Lat/Lon ######
r.polyLL <- r.poly;
projected <- "+proj=longlat";
r.polyLLNew <- spTransform(r.polyLL, CRS(projected))
r.polyLLNew <- r.polyLLNew[ , c(which(names(r.polyLLNew) %in% c("Yield", "SEASON", "COACH", "CROP", "ID")))];

#plot(r.polyLLNew)
#plot(r.polyLLNew, col=plotclrYield, border=F)
writePolyShape(r.polyLLNew, paste("/Users/tylerpittman/Farm/contour_yields_ll/quintiles", yearCurrent, "fromYears/", fb, "_SeedingRx_ll_", yearCurrent, "from", yearStartHarvest, "Yield", ".shp", sep=''));

eval(parse(text=paste("r.polyLLNew", yearStartHarvest, "<-", "r.polyLLNew", sep="")));
#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%
#assign('r.polyLLNew2014',r.polyLLNew2014,.GlobalEnv);
eval(parse(text=paste("assign('r.polyLLNew", year, "',r.polyLLNew", year, ",.GlobalEnv)",  sep="")));
assign('r.polyLLNewPoints',o2,.GlobalEnv);

cXDim = ((xmax(r.polyLLNew) - xmin(r.polyLLNew)))/0.00005;
cYDim = ((ymax(r.polyLLNew) - ymin(r.polyLLNew)))/0.00005;
r <- raster(ncol=cXDim, nrow=cYDim);
extent(r) <- extent(r.polyLLNew);
r.polyLLNew_ras <- rasterize(r.polyLLNew, r, "Yield" ,fun='first');
eval(parse(text=paste("r.polyLLNew_ras", yearStartHarvest, "<-", "r.polyLLNew_ras", sep="")));
eval(parse(text=paste("assign('r.polyLLNew_ras", year, "',r.polyLLNew_ras", year, ",.GlobalEnv)",  sep="")));

};
#---(*)---#---(*)---#---(*)---#---(*)---#---(*)---#---(*)---#---(*)---#---(*)---#
#---(*)---#---(*)---#---(*)---#---(*)---#---(*)---#---(*)---#---(*)---#---(*)---#

fileList <- "";
fileList_ras <- "";
uid <- 1;
for (year in ((yearStartHarvest):(yearCurrent-1))) {
	makeYieldContourYears(yearStartHarvest=year, yearBoundary=yearBoundary, yearCurrent=yearCurrent, fy=fb);
	fileList <- paste(fileList, paste("r.polyLLNew", year, ",", sep=""), sep="");
	fileList_ras <- paste(fileList_ras, paste("r.polyLLNew_ras", year, ",", sep=""), sep="");
	eval(parse(text=paste("temp.data", "<-", "r.polyLLNew", year, sep=""))); 
	n <- length(slot(temp.data, "polygons"));
   	temp.data <- spChFIDs(temp.data, as.character(uid:(uid+n-1)));
	temp.data@data$YearHar <- year;
	uid <- uid + n;
	eval(parse(text=paste("r.polyLLNew", year, "<-", "temp.data", sep=""))); 
	
};
fileList <- gsub('.{1}$', '', fileList);
fileList_ras <- gsub('.{1}$', '', fileList_ras);

###
	#n <- length(slot(temp.data, "polygons"));
   	#temp.data <- spChFIDs(temp.data, as.character(uid:(uid+n-1)));
	#temp.data@data$YearHarvest <- year;
	#uid <- uid + n;
  	#fieldPolyNew <- spRbind(fieldPolyNew,temp.data);
###

#fieldYieldQuin <- eval(parse(text=paste("spRbind(", fileList, ")", sep="")));
#fieldYieldQuinJoin <- aggregate(fieldYieldQuin, by='Yield');

#
#o2 <- r.polyLLNewPoints;
#cXDim = ((max(o2$lon) - min(o2$lon)))/0.00005;
#cYDim = ((max(o2$lat) - min(o2$lat)))/0.00005;
#gridT = list(cellcentre.offset=c((min(o2$lon)),(min(o2$lat))), cellsize=c(0.00005,0.00005), cells.dim=c(cXDim,cYDim));
#grd2 <- GridTopology(gridT$cellcentre.offset, gridT$cellsize, gridT$cells.dim);
#SG2 <- SpatialGrid(grd2, proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"));
#r <- SG2;
#

#$$$----#$$$----#$$$----#$$$----#$$$----#$$$----#$$$----#$$$----#$$$----#
#cXDim = ((xmax(r.polyLLNew2014) - xmin(r.polyLLNew2014)))/0.00005;
#cYDim = ((ymax(r.polyLLNew2014) - ymin(r.polyLLNew2014)))/0.00005;
#r <- raster(ncol=cXDim, nrow=cYDim);
#extent(r) <- extent(r.polyLLNew2014);
#r.polyLLNew2014_ras <- rasterize(r.polyLLNew2014, r, "Yield" ,fun='first');
#r.polyLLNew2015_ras <- rasterize(r.polyLLNew2015, r, "Yield" ,fun='first');
#rasStack <- stack(r.polyLLNew2014_ras, r.polyLLNew2015_ras);
#rasMean <- calc(rasStack, fun=function(x){as.integer(mean(x))});
#rasMean <- calc(rasStack, mean);
#rasMean <- na.omit(rasMean);
#
#brksYield1 <- classIntervals(jitter(values(rasMean)), 5, style = "quantile", rtimes = 3, intervalClosure = c("left", "right"), dataPrecision = NULL);
#plotclrYield <- brewer.pal(5,"RdYlGn");
#colcodeYield1 <- findColours(brksYield1, plotclrYield)
#brksYield1 <- brksYield1$brks;
#colcodeYield <- colcodeYield1;
#brksYield <- brksYield1;
#
#t2.clip <- rasMean;
#t2.clip$Yield <- cut(t2.clip@data@values, breaks=brksYield)
#r <- raster(t2.clip, layer=2); #2 corresponds to the yield column of the t2.clip@data slot we want
##res(r); #cell grid cellsize from above;
#r.poly <- polygonizer(r);
#r.poly <- unionSpatialPolygons(r.poly, r.poly$DN);
#proj4string(r.poly) <- CRS("+proj=longlat +ellps=WGS84");
##gIsValid(r.poly, reason = TRUE);
#r.poly <- gSimplify(r.poly, tol = 0.00001, topologyPreserve=FALSE);
#r.poly.df <- data.frame(r.polyLLNew2015@data);
#row.names(r.poly.df) <- r.poly.df$Yield;
#fieldYieldQuinJoin <- SpatialPolygonsDataFrame(r.poly, data=r.poly.df);
#plot(fieldYieldQuinJoin, col=plotclrYield, border=T);
#$$$----#$$$----#$$$----#$$$----#$$$----#$$$----#$$$----#$$$----#$$$----#

lastYearFN <- eval(parse(text=paste((tail(unlist(strsplit(fileList, split=",")), n=1)))));
#tail(unlist(strsplit(fileList, split=",")), n=1) #get the last filename for year;
unlist(strsplit(fileList_ras, split=","));

if (length((yearStartHarvest):(yearCurrent-1)) == 1){
		fieldYieldQuinJoin <- temp.data; 
} else {
		#eval(parse(text=paste(unlist(strsplit(fileList_ras, split=','))[1], "<-", "resample(", fileList_ras, ",method='ngb')", sep="")));
		#r.polyLLNew_ras2014 <-resample(r.polyLLNew_ras2014, r.polyLLNew_ras2015, method="ngb");
		eval(parse(text=paste(unlist(strsplit(fileList_ras, split=','))[1], "<-", "resample(", unlist(strsplit(fileList_ras, split=','))[1], ",", unlist(strsplit(fileList_ras, split=','))[2], ",method='ngb')", sep="")));
		eval(parse(text=paste("rasStack", "<-", "stack(", fileList_ras, ")", sep="")));
		#rasStack <- stack(fileList_ras);
		rasMean <- calc(rasStack, fun=function(x){mean(sqrt(x))}); #this is questionable...;
		#rasMean <- calc(rasStack, sum);
		#rasMean <- overlay(rasStack, fun=mean); 

		brksYield1 <- classIntervals(jitter(values(rasMean)), 5, style = "quantile", rtimes = 3, intervalClosure = c("left", "right"), dataPrecision = NULL);
		plotclrYield <- brewer.pal(5,"RdYlGn");
		colcodeYield1 <- findColours(brksYield1, plotclrYield)
		brksYield1 <- brksYield1$brks;
		colcodeYield <- colcodeYield1;
		brksYield <- brksYield1;

		t2.clip <- rasMean;
		t2.clip$Yield <- cut(t2.clip@data@values, breaks=brksYield)
		#t2.clip$Yield <- cut(t2.clip@data@values, breaks=c(1:5))
		r <- raster(t2.clip, layer=2); #2 corresponds to the yield column of the t2.clip@data slot we want
		#res(r); #cell grid cellsize from above;
		r.poly <- polygonizer(r);
		r.poly <- unionSpatialPolygons(r.poly, r.poly$DN);
		proj4string(r.poly) <- CRS("+proj=longlat +ellps=WGS84");
		#gIsValid(r.poly, reason = TRUE);
		r.poly <- gSimplify(r.poly, tol = 0.00001, topologyPreserve=FALSE);
		r.poly.df <- data.frame(lastYearFN@data);
		row.names(r.poly.df) <- r.poly.df$Yield;
		fieldYieldQuinJoin <- SpatialPolygonsDataFrame(r.poly, data=r.poly.df);
};
fieldYieldQuinJoin <- fieldYieldQuinJoin[ , c(which(names(fieldYieldQuinJoin) %in% c("Yield", "SEASON", "COACH", "CROP")))];
#plot(fieldYieldQuinJoin, col=topo.colors(5, alpha = 1));
#plot(fieldYieldQuinJoin, col=plotclrYield);

writePolyShape(fieldYieldQuinJoin, paste("/Users/tylerpittman/Farm/contour_yields_ll/", fb, "_SeedingRx_ll_", yearCurrent, ".shp", sep=''));

projected <- "+proj=utm +zone=13 +ellps=WGS84";
fieldYieldQuinJoin <- spTransform(fieldYieldQuinJoin, CRS(projected))
writePolyShape(fieldYieldQuinJoin, paste("/Users/tylerpittman/Farm/contour_yields_proj/", fb, "_SeedingRx_ll_", yearCurrent, ".shp", sep=''));
#######################################################################

print(j);
j <- j + 1;

}
mclapply(1:counterField, fieldLoop);
#mclapply(1, fieldLoop);
####################### Looping ends here ###############################;

####################### All Looping ends here ###############################;

#stopCluster(cl);
