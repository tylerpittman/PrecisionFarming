## Takes aggreated farm field yield data from GreenStar 2630 monitor and cuts shapefile information for each field
## Tyler Pittman, 15 September 2018
# Rscript --no-save /Users/tylerpittman/Farm/STEP_10_SPLIT_YIELD_FILES_COMBINE_after_harvest_15September2018.R

## Below renames shapefiles and eliminates unrecognized characters and strings in Ubuntu.
## Use brew to install rename for OSX
##brew install rename
#rename 's/\(//g' *
#rename 's/\)//g' *
#rename 's/agventure-Agventure-//g' *
#rename 's/agventure-//g' *
#rename 's/ //g' *
## Below aggregates shapefiles together.
#for f in *.shp; do ogr2ogr -update -append yields_2018.shp $f -f "ESRI Shapefile" -s_srs EPSG:2956 -t_srs 'EPSG:2956'; done;
## Copy yields_2017.shp to folder /Users/tylerpittman/Farm/input/2017/yield

year = 2018;
setwd(paste("/Users/tylerpittman/Farm/input/", year, "/yield", sep=""));

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
#source("shape2poly.R"); # Reads in shape2poly function and others

# Set seeding crops below here in field for each year with given input boundary and combine yield data from previous year;
key <- matrix(c( 
900001,"Wheat","46-0-0-0","8-38-16-0","na","AC Brigade","na","May 15",
#900002,"Wheat","46-0-0-0","8-38-16-0","na","AC Brigade","na","May 15",	
900003,"Canola","46-0-0-0","na","15-30-0-0","21-0-0-24","InVigor L233P","May 9",
900004,"Lentils","CDC Greenland","na","TagTeam","8-38-16-0","na","May 4",
900005,"Wheat","46-0-0-0","8-38-16-0","na","AC Brigade","na","May 16",
900006,"Lentils","CDC Greenland","na","TagTeam","8-38-16-0","na","May 4",
900007,"Canola","46-0-0-0","na","15-30-0-0","21-0-0-24","InVigor L233P","May 9",
900008,"Wheat","46-0-0-0","8-38-16-0","na","AC Brigade","na","May 14",
900009,"Canola","46-0-0-0","na","15-30-0-0","21-0-0-24","InVigor L140P","May 12",
#900010,"Wheat","46-0-0-0","8-38-16-0","na","AC Brigade","na","May 13",
900011,"Wheat","46-0-0-0","8-38-16-0","na","AC Brigade","na","May 14",
900012,"Canola","46-0-0-0","na","15-30-0-0","21-0-0-24","InVigor L140P","May 12",
900013,"Lentils","CDC Greenland","na","TagTeam","8-38-16-0","na","May 8",
900014,"Canola","46-0-0-0","na","15-30-0-0","21-0-0-24","InVigor L233P","May 11",
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
900029,"Canola","46-0-0-0","na","15-30-0-0","21-0-0-24","InVigor L233P","May 10",
900030,"Lentils","CDC Greenland","na","TagTeam","8-38-16-0","na","May 5",
900031,"Wheat","46-0-0-0","8-38-16-0","na","AC Brigade","na","May 1"
), ncol=8, byrow=T);
key <- as.data.frame(key);
colnames(key) <- c("field", "Crop", "Tank 1", "Tank 2", "Tank 3", "Tank 4", "Tank 5", "SeedDate");
key$Crop <- as.character(levels(key$Crop))[key$Crop];

fieldYield <- readOGR(dsn=".", layer=paste("yields_", year, sep=""));
fieldYield.dbf <- read.dbf(paste("yields_", year, ".dbf", sep=""), header=TRUE);

####################### All Looping starts below here ###############################;
counterField <- length(key$field);
	fieldLoop <- function(y) {
	j <- y;

#j <- 13;
fb = key$field[j];
type = key$Crop[j];
setwd(paste("/Users/tylerpittman/Farm/input/", year, "/boundaries", sep=""));
fieldBoundary <- readShapeSpatial(paste(fb, "", ".shp", sep=""));
fieldBoundary.dbf <- read.dbf(paste(fb, "", ".dbf", sep=""), header=TRUE);

setwd(paste("/Users/tylerpittman/Farm/input/", year, "/yield", sep=""));
fieldYield.clip <-  intersect(fieldYield, fieldBoundary);
writePointsShape(fieldYield.clip, paste(fb, "_", year, "_yield.shp", sep=""));



print(fb);
print(j);
j <- j + 1;

}
mclapply(1:counterField, fieldLoop);
#mclapply(25, fieldLoop);
####################### All Looping ends here ###############################;

#stopCluster(cl);
