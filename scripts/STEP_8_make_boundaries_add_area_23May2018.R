## Takes farm field boundaries and adds field area to them for TopoJSON website
## Manually create new boundaries and yield folders under <current year> folder in ~/Farm/input
## Copies boundaries from boundariesSeedingClean<current year> to this new folder
## Tyler Pittman, 23 May 2018
# Rscript --no-save /Users/tylerpittman/Farm/STEP_8_make_boundaries_add_area_23May2018.R

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
library(spBayes); 	## important gives Bayesian spatial modeling functions\
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
	#j <- 1; 
fb = key$field[j];
type = key$Crop[j];
year = 2018;
#fb="900030";

setwd(paste("/Users/tylerpittman/Farm/boundariesSeedingClean", year, sep=""));
fieldBoundary <- readShapeSpatial(paste(fb, ".shp", sep=""));
fieldBoundary.dbf <- read.dbf(paste(fb, ".dbf", sep=""), header=TRUE);
#plot(fieldBoundary);
proj4string(fieldBoundary) <- CRS("+proj=longlat");

fieldBoundary <- spTransform(fieldBoundary, CRS("+proj=utm +zone=13 +ellps=WGS84"));
proj4string(fieldBoundary) <- CRS("+proj=utm +zone=13 +ellps=WGS84");
fieldArea <- gArea(fieldBoundary)/4046.86;
fieldArea <- format(round(fieldArea, 2), nsmall = 2);

r.poly <- fieldBoundary;
r.poly@data$Field_Area <- fieldArea;
r.poly@data$SP_ID_1 <- r.poly@data$SP_ID;

df <- cbind(levels(r.poly@data$SP_ID_1), levels(r.poly@data$SP_ID), r.poly@data$c_year, levels(r.poly@data$c_type), r.poly@data$Field_Area);
df <- data.frame(df);
colnames(df) <- c("SP_ID", "SP_ID_1", "c_year", "c_type", "Field_Area");
row.names(r.poly); #should be 0;
row.names(df) <- 0; #change to 0 from 1;
r.poly <- SpatialPolygonsDataFrame(r.poly, data=df);

r.poly <- r.poly[ , c(which(names(r.poly) %in% c("SP_ID_1", "c_year", "c_type", "Field_Area")))];
setwd(paste("/Users/tylerpittman/Farm/input/", year, "/boundaries", sep=""));
writePolyShape(r.poly, paste(fb, "_", year, "_boundary", ".shp", sep=""));
#writeOGR(fieldBoundary, ".", paste(fb, ".shp", sep=""), driver="ESRI Shapefile");
#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@


print(j);
print(fb);
j <- j + 1;

}
mclapply(1:counterField, fieldLoop);

####################### Looping ends here ###############################;

####################### All Looping ends here ###############################;

#stopCluster(cl);
