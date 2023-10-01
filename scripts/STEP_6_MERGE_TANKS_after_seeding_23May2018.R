## Takes farm field polygons from Topcon X30 monitor and merges seeding maps for each Tank
## Tyler Pittman, 23 May 2018
# Rscript --no-save /Users/tylerpittman/Farm/STEP_6_MERGE_TANKS_after_seeding_23May2018.R

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
library(imputation);
library(xtable);
library(spam);
library(parallel);	# Does multicore processing when using mclapply WORKS on Linux!
library(Hmisc); 	# gives cut2 function
source("shape2poly.R"); # Reads in shape2poly function and others


# Set seeding crops below here in field for each year with given input boundary and combine yield data from previous year;
##For some fields, have to go back to Topcon X30, engage all tanks and recreate reports;
key <- matrix(c( 
900001,"Wheat","46-0-0-0","8-38-16-0","na","AC Brigade","na",
900002,"Wheat","46-0-0-0","8-38-16-0","na","AC Brigade","na",	
900003,"Canola","46-0-0-0","na","15-30-0-0","21-0-0-24","InVigor L233P",
900004,"Lentils","CDC Greenland","na","TagTeam","8-38-16-0","na",
900005,"Wheat","46-0-0-0","8-38-16-0","na","AC Brigade","na",
900006,"Lentils","CDC Greenland","na","TagTeam","8-38-16-0","na",
900007,"Canola","46-0-0-0","na","15-30-0-0","21-0-0-24","InVigor L233P",
900008,"Wheat","46-0-0-0","8-38-16-0","na","AC Brigade","na",
900009,"Canola","46-0-0-0","na","15-30-0-0","21-0-0-24","InVigor L140P",
900010,"Wheat","46-0-0-0","8-38-16-0","na","AC Brigade","na",
900011,"Wheat","46-0-0-0","8-38-16-0","na","AC Brigade","na",
900012,"Canola","46-0-0-0","na","15-30-0-0","21-0-0-24","InVigor L140P",
900013,"Lentils","CDC Greenland","na","TagTeam","8-38-16-0","na",
900014,"Canola","46-0-0-0","na","15-30-0-0","21-0-0-24","InVigor L233P",
900015,"Lentils","CDC Greenland","na","TagTeam","8-38-16-0","na",
900016,"Lentils","CDC Greenland","na","TagTeam","8-38-16-0","na",
900017,"Wheat","46-0-0-0","8-38-16-0","na","AC Brigade","na",
900018,"Wheat","46-0-0-0","8-38-16-0","na","AC Brigade","na",
900020,"Lentils","CDC Greenland","na","TagTeam","8-38-16-0","na",
900022,"Wheat","46-0-0-0","8-38-16-0","na","AC Brigade","na",
900024,"Wheat","46-0-0-0","8-38-16-0","na","AC Brigade","na",
900025,"Wheat","46-0-0-0","8-38-16-0","na","AC Brigade","na",
900026,"Lentils","CDC Greenland","na","TagTeam","8-38-16-0","na",
900027,"Lentils","CDC Greenland","na","TagTeam","8-38-16-0","na",
900028,"Lentils","CDC Greenland","na","TagTeam","8-38-16-0","na",
##800028,						
900029,"Canola","46-0-0-0","na","15-30-0-0","21-0-0-24","InVigor L233P",
900030,"Lentils","CDC Greenland","na","TagTeam","8-38-16-0","na",
900031,"Wheat","46-0-0-0","8-38-16-0","na","AC Brigade","na"
), ncol=7, byrow=T);
key <- as.data.frame(key);
colnames(key) <- c("field", "Crop", "Tank 1", "Tank 2", "Tank 3", "Tank 4", "Tank 5");
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
#j <- 29;
#k <- 1;
#j <- 25;
#k <- 2;
fb = key$field[j];
type = key$Crop[j];
tn = colnames(key)[k+2];
pr = as.character(key[j, k+2]);
if(isTRUE(pr=="na")) {return(k <- k+1)} #breaks out of mclapply loop if Tank is not used;
year = 2017;
setwd("/Users/tylerpittman/Farm/boundariesSeedingCleanTopcon2015");
fieldBoundary <- readShapeSpatial(paste(fb, "", ".shp", sep=""));
fieldBoundary.dbf <- read.dbf(paste(fb, "", ".dbf", sep=""), header=TRUE);
setwd(paste("/Users/tylerpittman/Farm/TractorX30May232018/Clients/Agventure Farms Ltd/Agventure Farms Ltd/",fb, "/CoverageShapefiles/", sep=""));

tankDup <- length(grep(tn, Sys.glob(paste("*", grep(tn, getwd(), ignore.case=TRUE, value=TRUE), "*.shp", sep="")), ignore.case=TRUE, value=TRUE));
joinNames <- NULL;
for (c in 1:tankDup ) {
tmpTN <- paste("fieldSeed", c, sep=""); 
tmpT <- readShapeSpatial(grep(tn, Sys.glob(paste("*", grep(tn, getwd(), ignore.case=TRUE, value=TRUE), "*.shp", sep="")), ignore.case=TRUE, value=TRUE)[c]); #case insensitive i.e. Tank 1 or TANK 1 will load!;
assign(tmpTN, tmpT);
joinNames[c] <- tmpTN; 
}

fieldSeed <- eval(parse(text = paste("rbind(", paste(joinNames, collapse=", "), ", makeUniqueIDs = TRUE)")));
pathDir <- getwd(); 
dir.create("combined", showWarnings = TRUE, recursive = FALSE, mode = "0777");
setwd(paste(pathDir, "/combined", sep=""));
writePolyShape(fieldSeed, paste("combined_", tn, "_lb_ac.shp", sep=""));


print(k);
print(pr);
print(fb);
k <- k + 1;

}
mclapply(1:counterTank, tankLoop);
#mclapply(2:4, tankLoop);
####################### Looping ends here ###############################;
print(j);
j <- j + 1;

}
mclapply(1:counterField, fieldLoop);
#mclapply(25, fieldLoop);
### error happens with field #2....
####################### All Looping ends here ###############################;

#stopCluster(cl);
