## Takes farm field polygons from Topcon X30 monitor and creates boundaries
## Tyler Pittman, 23 May 2018
# Rscript --no-save /Users/tylerpittman/Farm/STEP_7_make_boundaries_after_seeding_23May2018.R

setwd("/Users/tylerpittman/Farm/boundariesSeedingClean2018");

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
library(imputation);
library(xtable);
library(spam);
library(parallel);	# Does multicore processing when using mclapply WORKS on Linux!
source("shape2poly.R"); # Reads in shape2poly function and others


# Set seeding crops below here in field for each year with given input boundary and combine yield data from previous year;
##For some fields, have to go back to Topcon X30, engage all tanks and recreate reports;
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


####################### All Looping starts below here ###############################;

####################### Looping starts below here ###############################;
counterField <- length(key$field);
	fieldLoop <- function(y) {
	j <- y;
	#j <- 9;
fb = key$field[j];
year = 2018;
key$Crop <- eval(parse(text=paste("key$Crop", year, sep="")));
type = key$Crop[j];
#fb="900030";

setwd(paste("/Users/tylerpittman/Farm/TractorX30May232018/Clients/Agventure Farms Ltd/Agventure Farms Ltd/",fb, "/CoverageShapefiles/combined/", sep=""));
fieldBoundary <- readShapeSpatial(paste(Sys.glob("*Tank 1*.shp"), sep=""));
fieldBoundary.dbf <- read.dbf(paste(Sys.glob("*Tank 1*.dbf"), sep=""), header=TRUE);
#plot(fieldBoundary);

#@#@#@#@#@#@#@#@#@#@ Below works for creating field boundaries #@#@#@#@#@#@#@#@#@#@
ID <- rep(fb, each=1, length=length(fieldBoundary));
fieldBoundarySP <- unionSpatialPolygons(fieldBoundary, ID, threshold=0.00000001, avoidGEOS=TRUE);
getSpPPolygonsIDSlots(fieldBoundarySP);
sapply(slot(fieldBoundarySP, "polygons"), function(i) slot(i, "ID"));
#plot(fieldBoundarySP);

IDs <- sapply(slot(fieldBoundarySP, "polygons"), function(i) slot(i, "ID")); #to retrieve the new IDs;
#df <- data.frame(c_year=rep(year, length(IDs)));
df <- data.frame(c_year=rep(year, length(IDs)), c_type=rep(type, length(IDs)), row.names=IDs);
fieldBoundarySP2 <- SpatialPolygonsDataFrame(fieldBoundarySP, data=df);
setwd("/Users/tylerpittman/Farm/boundariesSeedingClean2018");

fieldBoundary <- fieldBoundarySP2;
proj4string(fieldBoundary) <- CRS("+proj=longlat");
fieldBoundary.simp <- gSimplify(fieldBoundary, tol = 0.000001, topologyPreserve=FALSE);
fieldBoundary.simp <- gSimplify(fieldBoundary.simp, tol = 0.000001, topologyPreserve=TRUE);
fieldBoundary.df <- data.frame(fieldBoundary);
fieldBoundary <- SpatialPolygonsDataFrame(fieldBoundary.simp, data=fieldBoundary.df);
fieldBoundarySP2 <- fieldBoundary;

writePolyShape(fieldBoundarySP2, paste(fb, ".shp", sep=""));

#outerRings = Filter(function(f){f@ringDir==1},fieldBoundarySP2@polygons[[1]]@Polygons)
#outerBounds = SpatialPolygons(list(Polygons(outerRings,ID=1)))
#plot(outerBounds)
#object.size(as(outerBounds, "SpatialPolygons"));

#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@


print(j);
print(fb);
j <- j + 1;

}
mclapply(1:counterField, fieldLoop);
#mclapply(31, fieldLoop);
#mclapply(9, fieldLoop); #this field gives problems;
####################### Looping ends here ###############################;

####################### All Looping ends here ###############################;

#stopCluster(cl);
