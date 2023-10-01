## Takes farm field polygons from Topcon X30 monitor and creates seeding maps for each Tank
## Tyler Pittman, 23 May 2018
# Rscript --no-save /Users/tylerpittman/Farm/STEP_9_SEEDING_MAPS_after_seeding_not_sloughs_2018_23May2018.R

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
library(Hmisc); 	# gives cut2 function
source("shape2poly.R"); # Reads in shape2poly function and others


# Set seeding crops below here in field for each year with given input boundary and combine yield data from previous year;
key <- matrix(c( 
900001,"Wheat","46-0-0-0","8-38-16-0","na","AC Brigade","na","May 15",
#900002,"Wheat","46-0-0-0","8-38-16-0","na","AC Brigade","na","May 15",	
900003,"Canola","46-0-0-0","na","na","21-0-0-24","InVigor L233P","May 9",
900004,"Lentils","CDC Greenland","na","TagTeam","8-38-16-0","na","May 4",
900005,"Wheat","46-0-0-0","8-38-16-0","na","AC Brigade","na","May 16",
900006,"Lentils","CDC Greenland","na","TagTeam","8-38-16-0","na","May 4",
900007,"Canola","46-0-0-0","na","15-30-0-0","21-0-0-24","InVigor L233P","May 9",
900008,"Wheat","46-0-0-0","8-38-16-0","na","AC Brigade","na","May 14",
900009,"Canola","46-0-0-0","na","na","21-0-0-24","InVigor L140P","May 12",
#900010,"Wheat","46-0-0-0","8-38-16-0","na","AC Brigade","na","May 13",
900011,"Wheat","46-0-0-0","8-38-16-0","na","AC Brigade","na","May 14",
900012,"Canola","46-0-0-0","na","na","21-0-0-24","InVigor L140P","May 12",
900013,"Lentils","CDC Greenland","na","TagTeam","8-38-16-0","na","May 8",
900014,"Canola","46-0-0-0","na","na","21-0-0-24","InVigor L233P","May 11",
900015,"Lentils","CDC Greenland","na","TagTeam","8-38-16-0","na","May 7",
900016,"Lentils","CDC Greenland","na","TagTeam","8-38-16-0","na","May 5",
900017,"Wheat","46-0-0-0","8-38-16-0","na","AC Brigade","na","May 17",
900018,"Wheat","46-0-0-0","8-38-16-0","na","AC Brigade","na","May 15",
900020,"Lentils","CDC Greenland","na","TagTeam","8-38-16-0","na","May 7",
900022,"Wheat","46-0-0-0","8-38-16-0","na","AC Brigade","na","May 3",
900024,"Wheat","46-0-0-0","8-38-16-0","na","AC Brigade","na","May 2",
900025,"Wheat","46-0-0-0","8-38-16-0","na","AC Brigade","na","May 2",
900026,"Lentils","CDC Greenland","na","TagTeam","8-38-16-0","na","May 6",
900027,"Lentils","CDC Greenland","na","TagTeam","8-38-16-0","na","May 6",
900028,"Lentils","CDC Greenland","na","TagTeam","8-38-16-0","na","May 18",
##800028,						
900029,"Canola","46-0-0-0","na","na","21-0-0-24","InVigor L233P","May 10",
900030,"Lentils","CDC Greenland","na","TagTeam","8-38-16-0","na","May 5",
900031,"Wheat","46-0-0-0","8-38-16-0","na","AC Brigade","na","May 1"
), ncol=8, byrow=T);
key <- as.data.frame(key);
colnames(key) <- c("field", "Crop", "Tank 1", "Tank 2", "Tank 3", "Tank 4", "Tank 5", "SeedDate");
key$Crop <- as.character(levels(key$Crop))[key$Crop];


####################### All Looping starts below here ###############################;
counterField <- length(key$field);
	fieldLoop <- function(y) {
	j <- y;

####################### Looping starts below here ###############################;
counterTank <- 5;
	tankLoop <- function(z) {
	k <- z;

#j <- 13;
#k <- 5;
#j <- 24;
#k <- 1;
#j <- 6;
#k <- 3;
fb = key$field[j];
type = key$Crop[j];
seedDate = as.character(key$SeedDate[j]);
tn = colnames(key)[k+2];
pr = as.character(key[j, k+2]);
if(isTRUE(pr=="na")) {return(k <- k+1)} #breaks out of mclapply loop if Tank is not used;
year = 2018;
setwd(paste("/Users/tylerpittman/Farm/boundariesSeedingClean", year, sep=""));
#setwd(paste("/Users/tylerpittman/Farm/input/", year, "/boundaries", sep=""));
fieldBoundary <- readShapeSpatial(paste(fb, "", ".shp", sep=""));
fieldBoundary.dbf <- read.dbf(paste(fb, "", ".dbf", sep=""), header=TRUE);
setwd(paste("/Users/tylerpittman/Farm/TractorX30May232018/Clients/Agventure Farms Ltd/Agventure Farms Ltd/",fb, "/CoverageShapefiles/combined/", sep=""));
fieldSeed <- readShapeSpatial(grep(tn, Sys.glob(paste("*", grep(tn, getwd(), ignore.case=TRUE, value=TRUE), "*.shp", sep="")), ignore.case=TRUE, value=TRUE)); 
proj4string(fieldBoundary) <- CRS("+proj=longlat");
proj4string(fieldSeed) <- CRS("+proj=longlat");
#case insensitive i.e. Tank 1 or TANK 1 will load!;
#fieldSeed <- readShapeSpatial(paste(Sys.glob(paste("*", tn, "*.shp", sep="")), sep=""));
#grep(tn, Sys.glob(paste("*", grep(tn, getwd(), ignore.case=TRUE, value=TRUE), "*.shp", sep="")), ignore.case=TRUE, value=TRUE)
#Sys.glob(paste("*", grep(tn, getwd(), ignore.case=TRUE, value=TRUE), "*.shp", sep=""))
#writePolyShape(fieldSeed, "J_AG VENTURE_050415_1209_Tank 5_ lb_ac.shp");
#fieldSeed.dbf <- read.dbf(paste(Sys.glob(paste("*", tn, "*.shp", sep="")), sep=""), header=TRUE);
#fieldSeed <- readShapeSpatial(paste(Sys.glob("*Tank 1*.shp"), sep=""));
#fieldSeed.dbf <- read.dbf(paste(Sys.glob("*Tank 1*.dbf"), sep=""), header=TRUE);
#fieldSeed <- readShapeSpatial(paste(Sys.glob(paste("*", tn, "*.shp", sep="")), sep=""), IDvar="X_lb_ac");


###For canola field 900029, j<-24, k<-1 in 2018 modify;
#fieldSeed@data; 
#fieldSeed@data[which(fieldSeed@data$X_lb_ac == 0), ]$X_lb_ac <- 159.9;
#fieldSeed@data$X_lb_ac <- jitter(fieldSeed@data$X_lb_ac, amount=2);
#writePolyShape(fieldSeed, "/Users/tylerpittman/Farm/TractorX30May232018/Clients/Agventure Farms Ltd/Agventure Farms Ltd/900029/CoverageShapefiles/combined/combined_Tank 1_lb_ac.shp");



#--------*--------*--------* Check for self-intersecting polygons and cuts boundaries --------*--------*--------*--------*#;
##install RFigisGeo
#install.packages("devtools", dep=T);
#require(devtools); 
#install_github("RFigisGeo", "openfigis");

#interP <- gIsValid(fieldSeed, reason = TRUE); #checks for self-intersecting polygons;
#plot(fieldSeed, col = "grey")
#points(-108.105716040001, 50.8053385801112, pch = 3, cex = 3, col = "red")
#axis(1);axis(2);box()

#gArea(fieldSeed)/4046.86

fieldSeed.simp <- gSimplify(fieldSeed, tol = 0.000001, topologyPreserve=FALSE);
#fieldBoundary.simp <- gSimplify(fieldBoundary, tol = 0.000001, topologyPreserve=FALSE);
fieldBoundary.simp <- fieldBoundary;
#fieldSeed.simp <- gSimplify(fieldSeed.simp, tol = 0.000001, topologyPreserve=FALSE);
#fieldBoundary.simp <- gSimplify(fieldBoundary.simp, tol = 0.000001, topologyPreserve=FALSE);
#gIsValid(fieldSeed, reason = TRUE);
#gIsValid(fieldBoundary, reason = TRUE);
#gIsValid(fieldSeed.simp, reason = TRUE);
#gIsValid(fieldBoundary.simp, reason = TRUE);
fieldSeed.clip <-  intersect(fieldSeed.simp, fieldBoundary.simp);
fieldSeed@data$ID <- paste(row.names(fieldSeed@data)); 
#fieldSeed@data$ID <- paste(row.names(fieldSeed@data), row.names(fieldBoundary@data)); 

#plot(fieldSeed);
#plot(fieldBoundary, add=TRUE, col="red");

IDinput <- row.names(fieldSeed.simp);
IDcut <- row.names(fieldSeed.clip)
fieldSeed.clip.df <- data.frame(fieldSeed[match(row.names(fieldSeed.clip), fieldSeed$ID),], row.names=IDcut);
fieldSeed.clip.df$ID <- IDcut;
row.names(fieldSeed.clip) <- IDcut;
fieldSeed.clip2 <- SpatialPolygonsDataFrame(fieldSeed.clip, data=fieldSeed.clip.df);
fieldSeed <- fieldSeed.clip2;
fieldBoundary <- fieldBoundary.simp;


projected <- "+proj=utm +zone=13 +ellps=WGS84";
proj4string(fieldSeed) <- CRS("+proj=longlat")
fieldSeedProj <- spTransform(fieldSeed, CRS(projected));
fieldBoundaryProj <- spTransform(fieldBoundary, CRS(projected));
#seededArea <- gArea(fieldSeedProj)/4046.86;
#sum(sapply(slot(fieldSeed, "polygons"), function(i) slot(i, "area")))/4046.86; #does area right
area_acre <- sapply(slot(fieldSeedProj, "polygons"), function(i) slot(i, "area"))/4046.86;
ibs_poly <- fieldSeedProj@data$X_lb_ac * area_acre;

fieldSeed@data$ibs_poly <- ibs_poly;
fieldSeed@data$area_acre <- area_acre;
#--------*--------*--------*--------*--------*--------*--------*--------*--------*--------*--------*--------*--------*#;


#--------------------## This is where the grid interpolation section starts ##--------------------#; 

###### WORKING! SeedContour Projected ######
polyCentroids <- getSpPPolygonsLabptSlots(fieldSeed);
o2 <- cbind(polyCentroids, fieldSeed@data$X_lb_ac);
colnames(o2) <- c("lon", "lat", "Seed");
o2 <- as.data.frame(o2);
#o2 <- o2[complete.cases(o2),];

##if lat is negative then delete!
#o2nz <- o2[-c(which(o2$lat <= 0.000)), ];
#o2 <- o2nz;
meanSeed <- mean(o2$Seed)
medianSeed <- median(o2$Seed)
##put seed rate of 0.0 to median value, error with Topcon exporting file;
if ( length(which(fieldSeed@data$X_lb_ac == 0.0)) >= 1 ) {
	o2[c(which(o2$Seed == 0.0)), ]$Seed <- medianSeed;
	fieldSeed@data[c(which(fieldSeed@data$X_lb_ac == 0.0)), ]$X_lb_ac <- medianSeed;
};

nclrSeed <- 5; 
#display.brewer.all(n=NULL, type="div", select=NULL, exact.n=TRUE);
plotclrSeed <- brewer.pal(nclrSeed,"PiYG");
length(plotclrSeed);



############
#cXDim = ((xmax(fieldBoundary) - xmin(fieldBoundary)))/0.00005;
#cYDim = ((ymax(fieldBoundary) - ymin(fieldBoundary)))/0.00005;
#r <- raster(ncol=cXDim, nrow=cYDim);
#extent(r) <- extent(fieldBoundary);
#fieldSeed_ras <- rasterize(fieldSeed, r, "X_lb_ac" ,fun='first');

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
#############


brksSeed1 <- classIntervals(fieldSeed@data$X_lb_ac, 5, style = "kmeans", rtimes = 13, intervalClosure = "right");
#brksSeed1 <- classIntervals(jitter(fieldSeed@data$X_lb_ac), 5, style = "quantile", rtimes = 13, intervalClosure = "right");
brksSeed <- brksSeed1;
#brksSeed <- classIntervals(fieldSeed@data$X_lb_ac, 5, style = "quantile", rtimes = 3, intervalClosure = "right");
#brksSeed$brks <- sort(jitter(brksSeed$brks, amount=5));
#brksSeed$brks <- cut2(fieldSeed@data$X_lb_ac, g=5, onlycuts=TRUE);
#brksRounded <- round(brksSeed$brks, 3);
#brksSeed <- classIntervals(fieldSeed@data$X_lb_ac, 5, style = "fixed", fixedBreaks = brksRounded, rtimes = 3, intervalClosure = "right");
brksSeed$brks[1] <- 0;
colcodeSeed <- findColours(brksSeed, plotclrSeed);
brksSeed <- brksSeed$brks;
#plot(fieldSeed, col=colcodeSeed, border=T);

#@#@#@#@#@#@#@#@#@#@ Below works for creating seed polygons #@#@#@#@#@#@#@#@#@#@
##below is raster package way, but leave orphaned holes sometimes
colorOrder <- cbind(plotclrSeed, c(1:length(plotclrSeed)));
colorOrder <- as.data.frame(colorOrder);
colnames(colorOrder) <- c("colcode", "colorOrder");

fieldSeed@data$colcode <- colcodeSeed;
fieldSeed0 <- aggregate(fieldSeed, vars="colcode", sums=list(list(sum, "ibs_poly"), list(sum, "area_acre")), dissolve=TRUE);
#fieldSeed0 <- aggregate(fieldSeed, vars="colcode", sums=list(list(sum, "ibs_poly")), dissolve=TRUE);
IDs <- sapply(slot(fieldSeed0, "polygons"), function(i) slot(i, "ID")); #to retrieve the new IDs;
df <- data.frame(c_year=rep(year, length(IDs)), c_type=rep(type, length(IDs)), row.names=IDs);
fieldSeedSP2  <- spCbind(fieldSeed0, df);
fieldSeedSP2 <- fieldSeedSP2[match(colorOrder$colcode, fieldSeedSP2$colcode),]; #puts the right order for plotting colors;

#plot(fieldSeed, col=fieldSeedSP2@data$colcode, border=T);
#plot(fieldSeedSP2, col=plotclrSeed, border=T);
#plot(fieldSeedSP2, col=fieldSeedSP2@data$colcode, border=T);

#ID <- rep(fb, each=1, length=length(fieldBoundary));
#fieldSeedSP <- unionSpatialPolygons(fieldSeed, colcodeSeed, threshold=0.00000000001, avoidGEOS=TRUE);
#getSpPPolygonsIDSlots(fieldSeedSP);
#sapply(slot(fieldSeedSP, "polygons"), function(i) slot(i, "ID"));
#IDs <- sapply(slot(fieldSeedSP, "polygons"), function(i) slot(i, "ID")); #to retrieve the new IDs;
#df <- data.frame(c_year=rep(year, length(IDs)), c_type=rep(type, length(IDs)), row.names=IDs);
#fieldSeedSP2 <- SpatialPolygonsDataFrame(fieldSeedSP, data=df);
##plot(fieldSeedSP2, col=plotclrSeed, border=T);
##plot(fieldSeedSP2, col=row.names(df), border=T);
#plotclrA <- row.names(df);

#is1 <- unionSpatialPolygons(fieldSeedSP2, as.character(row.names(fieldSeedSP2))) 
#row.names(xtra1) <- plotclrSeed;
#xx1 <- spCbind(xx, xtra1)
#names(xx1)

#outerRings = Filter(function(f){f@ringDir==1},fieldBoundarySP2@polygons[[1]]@Polygons)
#outerBounds = SpatialPolygons(list(Polygons(outerRings,ID=1)))
#plot(outerBounds)
#object.size(as(outerBounds, "SpatialPolygons"));
#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@


coordinates(o2) <- c(1, 2)
proj4string(o2) <- CRS("+proj=longlat")
o2new <- spTransform(o2, CRS(projected))
proj4string(fieldBoundary) <- CRS("+proj=longlat")
fieldBoundaryNew <- spTransform(fieldBoundary, CRS(projected))
proj4string(fieldSeedSP2) <- CRS("+proj=longlat")
fieldSeedNew <- spTransform(fieldSeedSP2, CRS(projected))

fieldArea <- gArea(fieldBoundaryNew)/4046.86;

cDim = 200;
cSLon = ((max(o2new$lon) - min(o2new$lon)))/cDim;
cSLat = ((max(o2new$lat) - min(o2new$lat)))/cDim;
cSLonB = ((max(fieldBoundaryNew@bbox[2,]) - min(fieldBoundaryNew@bbox[2,])))/cDim;
cSLatB = ((max(fieldBoundaryNew@bbox[1,]) - min(fieldBoundaryNew@bbox[1,])))/cDim;


##below could be factor of 1.000005 or 0.999995 depending on lat/lon or projected coordinate system;
gridT = list(cellcentre.offset=c((min(o2new$lon)),(min(o2new$lat))), cellsize=c(cSLon,cSLat), cells.dim=c(cDim,cDim));
grd2 <- GridTopology(gridT$cellcentre.offset, gridT$cellsize, gridT$cells.dim)
SG2 <- SpatialGrid(grd2, proj4string=CRS(projected))
p2 <- idw(Seed ~ 1, o2new, newdata=SG2, nmax=1, maxdist=19, idp=6)
#p2 <- idw(Seed ~ 1, o2new, newdata=SG2, maxdist=5, nmax=4);
#p2 <- idw(Seed ~ 1, o2new, newdata=SG2, na.action = na.omit, nmax=4, maxdist=16.5, idp=2); 

#check to see if grid covers field;
#plot(o2new, col=colcodeSeed)
#plot(p2, col="Blue", add=T)

### This works for cutting edges! Choose p2B if large empty seed inside of boundary ###
t2 <- p2;
#fullgrid(t2) = FALSE;
fullgrid(t2) = TRUE; #if p2B
#t2.clip = t2[!is.na(overlay(t2, fieldBoundaryNew)),]; #Don't need to use this anymore since done above;
pim2 <- as.image.SpatialGridDataFrame(t2)
cl2 <- contourLines(pim2)
cL2 <- ContourLines2SLDF(cl2)
proj4string(cL2) <- CRS(projected)
#cL2clip <- gIntersection(fieldBoundary, cL2, byid=TRUE, id=NULL); #this is from rgeos package;

names(attr(colcodeSeed, "table"));
brksRounded <- round(brksSeed, 3);
brksSeedRounded <- classIntervals(fieldSeed@data$X_lb_ac, 5, style = "fixed", fixedBreaks = brksRounded, rtimes = 3, intervalClosure = "right");
brksSeedRounded$brks[1] <- 0;
colcodeSeedRounded <- findColours(brksSeedRounded, plotclrSeed);


t2$Seed <- cut(t2$var1.pred, right=FALSE, breaks=brksSeed)
summary(t2$Seed);
#t2.clip$Seed <- cut(t2.clip$var1.pred, right=FALSE, breaks=brksSeed)
#summary(t2.clip$Seed); 
( r <- raster(t2, layer=3) ) #3 corresponds to the seed column of the t2.clip@data slot we want
#plot(r);
r.polyCol <- rasterToPolygons(r, dissolve=TRUE)
Seed <- c(1:(length(brksSeed)-1)); #Five quintiles;
r.poly <- fieldSeedNew; #overwrite when using polygons from Topcon X30 as opposed to polygon centroids higher up;
r.poly@data <- cbind(Seed, r.poly@data);
#r.poly@data$ADSFLDID = fb;
r.poly@data$SEASON = paste("Seeding ", tn, sep="");
r.poly@data$c_year = year;
r.poly@data$c_type = type;
#r.poly@data <- cbind(SP_ID, r.poly@data);

#seededArea <- sum(fieldSeedNew@data$area_acre);
#seededArea <- gArea(fieldSeedNew)/4046.86;
#meanSeed <- mean(o2[-c(which(o2$Seed <= 0.1)), ]$Seed);
#meanSeed <- mean(o2$Seed);
#r.poly@data$Field_Area <- fieldArea;
#r.poly@data$Field_Seeded_Area <- seededArea;
#r.poly@data$Field_Mean_Seed <- meanSeed;

SP_ID_1 <- r.poly@data$Seed
COACH <- r.poly@data$Seed
CROP <- r.poly@data$Seed
PRODUCT <- r.poly@data$Seed
LABEL <- names(attr(colcodeSeedRounded, "table"))
Quintile_Area <- r.poly@data$Seed
Quintile_Mean_Seed <- r.poly@data$Seed
options(scipen=3);

for(i in 1:length(r.poly@data$Seed)){
	if (r.poly@data$Seed[[i]] == 1) {COACH[[i]] = "s2"} else	
	if (r.poly@data$Seed[[i]] == 2) {COACH[[i]] = "s1"} else
	if (r.poly@data$Seed[[i]] == 3) {COACH[[i]] = "normal"} else
	if (r.poly@data$Seed[[i]] == 4) {COACH[[i]] = "a1"} else
	if (r.poly@data$Seed[[i]] == 5) {COACH[[i]] = "a2"} else
	COACH[[i]] = NA;
	SP_ID_1[[i]] = as.character(fb);	
	CROP[[i]] = type;
    PRODUCT[[i]] = pr;
	LABEL[[i]] = names(attr(colcodeSeed, "table"))[i];
	#Quintile_Area[[i]] <- gArea(r.poly[i,])/4046.86;
	Quintile_Area[[i]] <- fieldSeedNew@data$area_acre[i];
	#Quintile_Mean_Seed[[i]] <- mean(o2[c(which(o2$Seed > brksSeed[i] & o2$Seed <= brksSeed[i+1])), ]$Seed);
	#Quintile_Mean_Seed[[i]] <- mean(fieldSeed@data[c(which(fieldSeed@data$X_lb_ac > brksSeed[i] & fieldSeed@data$X_lb_ac <= brksSeed[i+1])), ]$X_lb_ac);
	Quintile_Mean_Seed[[i]] <- fieldSeedNew@data$ibs_poly[i]/fieldSeedNew@data$area_acre[i];
};
r.poly@data$SP_ID_1 <- SP_ID_1;
r.poly@data$PRODUCT <- PRODUCT;
r.poly@data$COACH <- COACH;
r.poly@data$LABEL <- LABEL;
r.poly@data$Quintile_Area <- Quintile_Area;
r.poly@data$Quintile_Mean_Seed <- Quintile_Mean_Seed;

r.poly@data$Field_Area <- fieldArea;
r.poly@data$Field_Seeded_Area <- sum(as.numeric(r.poly@data$Quintile_Area));
seededArea <- sum(as.numeric(r.poly@data$Quintile_Area));
r.poly@data$Field_Mean_Seed <- sum(fieldSeedNew@data$ibs_poly)/sum(fieldSeedNew@data$area_acre);

r.poly@data$Field_Seeded_Area <- format(r.poly@data$Field_Seeded_Area, digits=5);
r.poly@data$Field_Area <- format(r.poly@data$Field_Area, digits=5);
r.poly@data$Field_Mean_Seed <- round(r.poly@data$Field_Mean_Seed, 3);
r.poly@data$Quintile_Area <- format(r.poly@data$Quintile_Area, digits=4);
r.poly@data$Quintile_Mean_Seed <- format(r.poly@data$Quintile_Mean_Seed, digits=2);
r.poly@data$SeedDate <- seedDate;

###
areaBoundary <- sapply(slot(fieldBoundaryNew, "polygons"), function(i) slot(i, "area"));
areaBoundary <- areaBoundary/4046.86;
areaPolys <- sapply(slot(r.poly, "polygons"), function(i) slot(i, "area")); #to retrieve the new areas;
areaPolys <- areaPolys/4046.86;
areaSeeded <- sum(areaPolys);
areaBoundary <- round(areaBoundary, 3);
areaPolys <- round(areaPolys, 3);
areaSeeded <- round(areaSeeded, 3);
r.poly@data$Field_Seeded_Area <- areaSeeded;
r.poly@data$Quintile_Area <- areaPolys;
r.poly@data$Field_Area <- areaBoundary;
###


#plot(r.poly)
#plot(r.poly, col=plotclrSeed, border=T);
#gArea(r.poly[which(r.poly$COACH == "s2"),])/4046.86
#gArea(r.poly[which(r.poly$COACH == "s2"),])/4046.86
#plot(r.poly[which(r.poly$COACH == "s2"),]);
#mean(o2[c(which(o2$Seed > brksSeed[1] & o2$Seed <= brksSeed[2])), ]$Seed);

writePolyShape(r.poly, paste("/Users/tylerpittman/Farm/seedShapefiles/", year, "/", fb, "_seed", gsub(" ", "", tn, fixed=T), "_proj_", year, ".shp", sep=""));
##Don't do the below, as it redundant and messes up STEP_8 boundaries;
#writePolyShape(fieldBoundaryProj, paste("/Users/tylerpittman/Farm/input/", year, "/boundaries/", fb, "_", year, "_boundary", ".shp", sep=""));

#### Below is to replace boundary attributes, not always needed ###
#fieldBoundaryNew@data$c_year = year;
#fieldBoundaryNew@data$c_type = type;
#fieldBoundaryNew@data$Field_Area <- fieldArea;
#fieldBoundaryNew@data$Field_Area <- format(fieldBoundaryNew@data$Field_Area, digits=5);
#writePolyShape(fieldBoundaryNew, paste("/Users/tylerpittman/Farm/input/2015/boundaries/", fb, "_2015_boundary", ".shp", sep=''));


png(paste("/Users/tylerpittman/Farm/seedShapefiles/", year, "/images/", fb, "_seed", gsub(" ", "", tn, fixed=T), "_proj_", year, ".png", sep=''), res=100, pointsize = 12, units = "px", bg = "transparent");
par(xpd=F, mar=c(0,0,0,1), oma=c(0,0,2,1));
#image(t2.clip, col = plotclrSeed, add = FALSE, breaks=brksSeed, oldstyle = FALSE, useRaster=TRUE);
#contour(t2.clip, lwd=0.5, col=plotclrSeed, levels=brksSeed, drawlabels=F, add=T);
#plot(fieldBoundaryNew, lwd=2.0, add=T);
plot(r.poly, col=plotclrSeed, border=F);
title(paste(fb, pr, "Seed Quintiles (in lb/acre) \n From", paste("Seeding ", tn, sep=""), year), outer=TRUE);
par(lwd=2); legend("bottomright",  box.lwd=1, legend=rev(names(attr(colcodeSeed, "table"))), fill=rev(plotclrSeed), cex=0.9, bg="transparent", xjust=0);
par(lwd=2); legend("bottomleft",  box.lwd=1, legend=cbind(c("Field Area (acres):", fieldArea), c("Seeded Area (acres):", seededArea), c("Mean Rate (lb/acre):", meanSeed)), cex=0.9, bg="transparent", xjust=0);
#par(lwd=2); legend("bottomright",  box.lwd=1, legend=rev(brksSeed[-1]), fill=rev(plotclrSeed), cex=0.9, bg="transparent", xjust=0);
dev.off();


print(k);
print(pr);
print(fb);
k <- k + 1;

}
mclapply(1:counterTank, tankLoop);
#mclapply(5, tankLoop);
####################### Looping ends here ###############################;
print(j);
j <- j + 1;

}
mclapply(1:counterField, fieldLoop);
#mclapply(25, fieldLoop);
### error happens with field #2....
####################### All Looping ends here ###############################;

#stopCluster(cl);
