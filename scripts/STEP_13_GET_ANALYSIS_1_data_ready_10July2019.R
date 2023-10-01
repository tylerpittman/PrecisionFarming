## Takes farm combine yield data and seeding data from 2015 to 2018
## Tyler Pittman, 1 August 2019
# Rscript --no-save /Users/tylerpittman/Farm/STEP_13_GET_ANALYSIS_1_data_ready_10July2019.R

Sys.setenv(PATH="/usr/bin:/opt/local/bin:/opt/local/sbin:/Library/Frameworks/GDAL.framework/Programs:/usr/local/mysql/bin:/opt/ImageMagick/bin:/Library/Frameworks/R.framework/Versions/3.2/Resources/bin");
Sys.getenv("PATH");
Sys.which("gdal_polygonize.py");
sessionInfo();

setwd("/Users/tylerpittman/Farm/agYieldProject/data");

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
library(rgdal);
library(parallel);	# Does multicore processing when using mclapply WORKS on Linux!
#library(greenbrown); #has many complicated dependencies to install, for equal area parcel division of shapefiles;
#source("shape2poly.R"); # Reads in shape2poly function and others
#source("polygonizer.R"); # Reads in polygonizer function

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

#######################################################################;


#--0--0--0--0--0--0--0--0--0--#;
## Fix .dbf for fields 900025, 900026 and 900006 in 2015, 2016, 2017, 2018 for seedTank4;
library(geojsonio);

### Field 900025;
setwd(paste("/Users/tylerpittman/Farm/seedShapefiles/2015", sep=""));
fixField <- readOGR(dsn=".", layer=paste("900025_seedTank4_proj_2015", sep=""));
#fixFieldtmp <- readOGR(dsn=".", layer=paste("900025_seedTank2_proj_2015", sep=""));
#fixFieldtmp@data;
topoField <- topojson_read("/Users/tylerpittman/Farm/farm/field/900025_seedTank4_proj_2015_ca.json");
ffdf <- as.data.frame(topoField[, c(2:19)]);
ffdf <- ffdf[, c(1:18)];
fixField@data <- ffdf;
writeOGR(fixField, ".", "900025_seedTank4_proj_2015", overwrite_layer=TRUE, driver="ESRI Shapefile");

setwd(paste("/Users/tylerpittman/Farm/seedShapefiles/2016", sep=""));
fixField <- readOGR(dsn=".", layer=paste("900025_seedTank4_proj_2016", sep=""));
#fixFieldtmp <- readOGR(dsn=".", layer=paste("900025_seedTank2_proj_2016", sep=""));
#fixFieldtmp@data;
topoField <- topojson_read("/Users/tylerpittman/Farm/farm/field/900025_seedTank4_proj_2016_ca.json");
ffdf <- as.data.frame(topoField[, c(2:19)]);
ffdf <- ffdf[, c(1:18)];
fixField@data <- ffdf;
writeOGR(fixField, ".", "900025_seedTank4_proj_2016", overwrite_layer=TRUE, driver="ESRI Shapefile");

setwd(paste("/Users/tylerpittman/Farm/seedShapefiles/2017", sep=""));
fixField <- readOGR(dsn=".", layer=paste("900025_seedTank4_proj_2017", sep=""));
#fixFieldtmp <- readOGR(dsn=".", layer=paste("900025_seedTank2_proj_2017", sep=""));
#fixFieldtmp@data;
topoField <- topojson_read("/Users/tylerpittman/Farm/farm/field/900025_seedTank4_proj_2017_ca.json");
ffdf <- as.data.frame(topoField[, c(2:19)]);
ffdf <- ffdf[, c(1:18)];
fixField@data <- ffdf;
writeOGR(fixField, ".", "900025_seedTank4_proj_2017", overwrite_layer=TRUE, driver="ESRI Shapefile");

setwd(paste("/Users/tylerpittman/Farm/seedShapefiles/2018", sep=""));
fixField <- readOGR(dsn=".", layer=paste("900025_seedTank4_proj_2018", sep=""));
#fixFieldtmp <- readOGR(dsn=".", layer=paste("900025_seedTank2_proj_2018", sep=""));
#fixFieldtmp@data;
topoField <- topojson_read("/Users/tylerpittman/Farm/farm/field/900025_seedTank4_proj_2018_ca.json");
ffdf <- as.data.frame(topoField[, c(2:19)]);
ffdf <- ffdf[, c(1:18)];
fixField@data <- ffdf;
writeOGR(fixField, ".", "900025_seedTank4_proj_2018", overwrite_layer=TRUE, driver="ESRI Shapefile");

### Field 900026;
setwd(paste("/Users/tylerpittman/Farm/seedShapefiles/2015", sep=""));
fixField <- readOGR(dsn=".", layer=paste("900026_seedTank4_proj_2015", sep=""));
#fixFieldtmp <- readOGR(dsn=".", layer=paste("900026_seedTank2_proj_2015", sep=""));
#fixFieldtmp@data;
topoField <- topojson_read("/Users/tylerpittman/Farm/farm/field/900026_seedTank4_proj_2015_ca.json");
ffdf <- as.data.frame(topoField[, c(2:19)]);
ffdf <- ffdf[, c(1:18)];
fixField@data <- ffdf;
writeOGR(fixField, ".", "900026_seedTank4_proj_2015", overwrite_layer=TRUE, driver="ESRI Shapefile");

setwd(paste("/Users/tylerpittman/Farm/seedShapefiles/2016", sep=""));
fixField <- readOGR(dsn=".", layer=paste("900026_seedTank4_proj_2016", sep=""));
#fixFieldtmp <- readOGR(dsn=".", layer=paste("900026_seedTank2_proj_2016", sep=""));
#fixFieldtmp@data;
topoField <- topojson_read("/Users/tylerpittman/Farm/farm/field/900026_seedTank4_proj_2016_ca.json");
ffdf <- as.data.frame(topoField[, c(2:19)]);
ffdf <- ffdf[, c(1:18)];
fixField@data <- ffdf;
writeOGR(fixField, ".", "900026_seedTank4_proj_2016", overwrite_layer=TRUE, driver="ESRI Shapefile");

setwd(paste("/Users/tylerpittman/Farm/seedShapefiles/2017", sep=""));
fixField <- readOGR(dsn=".", layer=paste("900026_seedTank4_proj_2017", sep=""));
#fixFieldtmp <- readOGR(dsn=".", layer=paste("900026_seedTank2_proj_2017", sep=""));
#fixFieldtmp@data;
topoField <- topojson_read("/Users/tylerpittman/Farm/farm/field/900026_seedTank4_proj_2017_ca.json");
ffdf <- as.data.frame(topoField[, c(2:19)]);
ffdf <- ffdf[, c(1:18)];
fixField@data <- ffdf;
writeOGR(fixField, ".", "900026_seedTank4_proj_2017", overwrite_layer=TRUE, driver="ESRI Shapefile");

setwd(paste("/Users/tylerpittman/Farm/seedShapefiles/2018", sep=""));
fixField <- readOGR(dsn=".", layer=paste("900026_seedTank4_proj_2018", sep=""));
#fixFieldtmp <- readOGR(dsn=".", layer=paste("900026_seedTank2_proj_2018", sep=""));
#fixFieldtmp@data;
topoField <- topojson_read("/Users/tylerpittman/Farm/farm/field/900026_seedTank4_proj_2018_ca.json");
ffdf <- as.data.frame(topoField[, c(2:19)]);
ffdf <- ffdf[, c(1:18)];
fixField@data <- ffdf;
writeOGR(fixField, ".", "900026_seedTank4_proj_2018", overwrite_layer=TRUE, driver="ESRI Shapefile");

### Field 900025;
setwd(paste("/Users/tylerpittman/Farm/seedShapefiles/2015", sep=""));
fixField <- readOGR(dsn=".", layer=paste("900025_seedTank4_proj_2015", sep=""));
#fixFieldtmp <- readOGR(dsn=".", layer=paste("900025_seedTank2_proj_2015", sep=""));
#fixFieldtmp@data;
topoField <- topojson_read("/Users/tylerpittman/Farm/farm/field/900025_seedTank4_proj_2015_ca.json");
ffdf <- as.data.frame(topoField[, c(2:19)]);
ffdf <- ffdf[, c(1:18)];
fixField@data <- ffdf;
writeOGR(fixField, ".", "900025_seedTank4_proj_2015", overwrite_layer=TRUE, driver="ESRI Shapefile");

setwd(paste("/Users/tylerpittman/Farm/seedShapefiles/2016", sep=""));
fixField <- readOGR(dsn=".", layer=paste("900025_seedTank4_proj_2016", sep=""));
#fixFieldtmp <- readOGR(dsn=".", layer=paste("900025_seedTank2_proj_2016", sep=""));
#fixFieldtmp@data;
topoField <- topojson_read("/Users/tylerpittman/Farm/farm/field/900025_seedTank4_proj_2016_ca.json");
ffdf <- as.data.frame(topoField[, c(2:19)]);
ffdf <- ffdf[, c(1:18)];
fixField@data <- ffdf;
writeOGR(fixField, ".", "900025_seedTank4_proj_2016", overwrite_layer=TRUE, driver="ESRI Shapefile");

setwd(paste("/Users/tylerpittman/Farm/seedShapefiles/2017", sep=""));
fixField <- readOGR(dsn=".", layer=paste("900025_seedTank4_proj_2017", sep=""));
#fixFieldtmp <- readOGR(dsn=".", layer=paste("900025_seedTank2_proj_2017", sep=""));
#fixFieldtmp@data;
topoField <- topojson_read("/Users/tylerpittman/Farm/farm/field/900025_seedTank4_proj_2017_ca.json");
ffdf <- as.data.frame(topoField[, c(2:19)]);
ffdf <- ffdf[, c(1:18)];
fixField@data <- ffdf;
writeOGR(fixField, ".", "900025_seedTank4_proj_2017", overwrite_layer=TRUE, driver="ESRI Shapefile");

setwd(paste("/Users/tylerpittman/Farm/seedShapefiles/2018", sep=""));
fixField <- readOGR(dsn=".", layer=paste("900025_seedTank4_proj_2018", sep=""));
#fixFieldtmp <- readOGR(dsn=".", layer=paste("900025_seedTank2_proj_2018", sep=""));
#fixFieldtmp@data;
topoField <- topojson_read("/Users/tylerpittman/Farm/farm/field/900025_seedTank4_proj_2018_ca.json");
ffdf <- as.data.frame(topoField[, c(2:19)]);
ffdf <- ffdf[, c(1:18)];
fixField@data <- ffdf;
writeOGR(fixField, ".", "900025_seedTank4_proj_2018", overwrite_layer=TRUE, driver="ESRI Shapefile");

### Field 900006;
setwd(paste("/Users/tylerpittman/Farm/seedShapefiles/2015", sep=""));
fixField <- readOGR(dsn=".", layer=paste("900006_seedTank4_proj_2015", sep=""));
#fixFieldtmp <- readOGR(dsn=".", layer=paste("900006_seedTank2_proj_2015", sep=""));
#fixFieldtmp@data;
topoField <- topojson_read("/Users/tylerpittman/Farm/farm/field/900006_seedTank4_proj_2015_ca.json");
ffdf <- as.data.frame(topoField[, c(2:19)]);
ffdf <- ffdf[, c(1:18)];
fixField@data <- ffdf;
writeOGR(fixField, ".", "900006_seedTank4_proj_2015", overwrite_layer=TRUE, driver="ESRI Shapefile");

setwd(paste("/Users/tylerpittman/Farm/seedShapefiles/2016", sep=""));
fixField <- readOGR(dsn=".", layer=paste("900006_seedTank4_proj_2016", sep=""));
#fixFieldtmp <- readOGR(dsn=".", layer=paste("900006_seedTank2_proj_2016", sep=""));
#fixFieldtmp@data;
topoField <- topojson_read("/Users/tylerpittman/Farm/farm/field/900006_seedTank4_proj_2016_ca.json");
ffdf <- as.data.frame(topoField[, c(2:19)]);
ffdf <- ffdf[, c(1:18)];
fixField@data <- ffdf;
writeOGR(fixField, ".", "900006_seedTank4_proj_2016", overwrite_layer=TRUE, driver="ESRI Shapefile");

setwd(paste("/Users/tylerpittman/Farm/seedShapefiles/2017", sep=""));
fixField <- readOGR(dsn=".", layer=paste("900006_seedTank4_proj_2017", sep=""));
#fixFieldtmp <- readOGR(dsn=".", layer=paste("900006_seedTank2_proj_2017", sep=""));
#fixFieldtmp@data;
topoField <- topojson_read("/Users/tylerpittman/Farm/farm/field/900006_seedTank4_proj_2017_ca.json");
ffdf <- as.data.frame(topoField[, c(2:19)]);
ffdf <- ffdf[, c(1:18)];
fixField@data <- ffdf;
writeOGR(fixField, ".", "900006_seedTank4_proj_2017", overwrite_layer=TRUE, driver="ESRI Shapefile");

setwd(paste("/Users/tylerpittman/Farm/seedShapefiles/2018", sep=""));
fixField <- readOGR(dsn=".", layer=paste("900006_seedTank4_proj_2018", sep=""));
#fixFieldtmp <- readOGR(dsn=".", layer=paste("900006_seedTank2_proj_2018", sep=""));
#fixFieldtmp@data;
topoField <- topojson_read("/Users/tylerpittman/Farm/farm/field/900006_seedTank4_proj_2018_ca.json");
ffdf <- as.data.frame(topoField[, c(2:19)]);
ffdf <- ffdf[, c(1:18)];
fixField@data <- ffdf;
writeOGR(fixField, ".", "900006_seedTank4_proj_2018", overwrite_layer=TRUE, driver="ESRI Shapefile");
#--0--0--0--0--0--0--0--0--0--#;


#--1--1--1--1--1--1--1--1--1--#;
setwd(paste("/Users/tylerpittman/Farm/input/2015/yield", sep=""));
fieldYield2015 <- readOGR(dsn=".", layer=paste("yields_2015", sep=""));
fieldYield2015.dbf <- read.dbf(paste("yields_2015", ".dbf", sep=""), header=TRUE);

setwd(paste("/Users/tylerpittman/Farm/input/2016/yield", sep=""));
fieldYield2016 <- readOGR(dsn=".", layer=paste("yields_2016", sep=""));
fieldYield2016.dbf <- read.dbf(paste("yields_2016", ".dbf", sep=""), header=TRUE);

setwd(paste("/Users/tylerpittman/Farm/input/2017/yield", sep=""));
fieldYield2017 <- readOGR(dsn=".", layer=paste("yields_2017", sep=""));
fieldYield2017.dbf <- read.dbf(paste("yields_2017", ".dbf", sep=""), header=TRUE);
## fieldYield2017[which(fieldYield2017$Product0 == "2 Row"), ]; #combined Aug 26 to 27, 2017 is Barley;

setwd(paste("/Users/tylerpittman/Farm/input/2018/yield", sep=""));
fieldYield2018 <- readOGR(dsn=".", layer=paste("yields_2018", sep=""));
fieldYield2018.dbf <- read.dbf(paste("yields_2018", ".dbf", sep=""), header=TRUE);
## fieldYield2018[which(fieldYield2018$Product0 == "0"), ]; #combined Aug 25 to Sep 9, is Canola;

setwd("/Users/tylerpittman/Farm/agYieldProject/data");

###;
fieldYield2015@data$timeFull <- strptime(fieldYield2015@data$Time, format="%m/%d/%Y %I:%M:%S %p");
fieldYield2015@data$Harv_Date <- format(as.Date(fieldYield2015@data$timeFull), format="%B %d");

fieldYield2016@data$timeFull <- strptime(fieldYield2016@data$Time, format="%m/%d/%Y %I:%M:%S %p");
fieldYield2016@data$Harv_Date <- format(as.Date(fieldYield2016@data$timeFull), format="%B %d");

fieldYield2017@data$timeFull <- strptime(fieldYield2017@data$Time, format="%m/%d/%Y %I:%M:%S %p");
fieldYield2017@data$Harv_Date <- format(as.Date(fieldYield2017@data$timeFull), format="%B %d");

fieldYield2018@data$timeFull <- strptime(fieldYield2018@data$Time, format="%m/%d/%Y %I:%M:%S %p");
fieldYield2018@data$Harv_Date <- format(as.Date(fieldYield2018@data$timeFull), format="%B %d");
###;

##colnames(fieldYield2015@data);
##head(fieldYield2015@data, 3);
fieldYield2015sm <- fieldYield2015[, c("Elevatio", "DryYield", "Product0", "Machine", "ProcYear", "Harv_Date")];
writePointsShape(fieldYield2015sm, paste("yields_sm_2015.shp", sep=""));
fieldYield2016sm <- fieldYield2016[, c("Elevatio", "DryYield", "Product0", "Machine", "ProcYear", "Harv_Date")];
writePointsShape(fieldYield2016sm, paste("yields_sm_2016.shp", sep=""));
fieldYield2017sm <- fieldYield2017[, c("Elevatio", "DryYield", "Product0", "Machine", "ProcYear", "Harv_Date")];
writePointsShape(fieldYield2017sm, paste("yields_sm_2017.shp", sep=""));
fieldYield2018sm <- fieldYield2018[, c("Elevatio", "DryYield", "Product0", "Machine", "ProcYear", "Harv_Date")];
writePointsShape(fieldYield2018sm, paste("yields_sm_2018.shp", sep=""));
#--1--1--1--1--1--1--1--1--1--#;

setwd("/Users/tylerpittman/Farm/agYieldProject/data");
fieldYield2015 <- readOGR(dsn=".", layer=paste("yields_sm_2015", sep=""));
fieldYield2016 <- readOGR(dsn=".", layer=paste("yields_sm_2016", sep=""));
fieldYield2017 <- readOGR(dsn=".", layer=paste("yields_sm_2017", sep=""));
fieldYield2018 <- readOGR(dsn=".", layer=paste("yields_sm_2018", sep=""));

yieldData <- rbind(fieldYield2015@data, fieldYield2016@data, fieldYield2017@data, fieldYield2018@data);
save(yieldData, file=paste("yieldData2015to2018", ".RData", sep=""));



#--2--2--2--2--2--2--2--2--2--#;
load(paste("yieldData2015to2018", ".RData", sep=""));
length(yieldData$ProcYear); #4,754,809 data points over four years;
str(yieldData);

levels(yieldData$Product0);
##1] "Brigade"    "Greenland"  "L156H"      "calvi"      "Greenstar"  "2 Row"      "Copeland"   "L140P"      "L252"       "0"          "Canola"     "Green Star"
cha1 <- yieldData[which(yieldData$Product0 == "0"), ];
unique(cha1$ProcYear); #only happened in 2018;

cha2 <- yieldData[which(yieldData$Product0 == "2 Row"), ];
unique(cha2$ProcYear); #only happened in 2017;
unique(cha2$Machine); #only S670 2;

cha3 <- yieldData[which(yieldData$Machine == "2"), ];
unique(cha3$ProcYear); #only happened in 2018;
unique(cha3$Product0); #0 must be Canola?;

cha4 <- yieldData[which(yieldData$Machine == "S670 2"), ];
unique(cha4$ProcYear); #only happened in 2017;
unique(cha4$Product0); #2 Row is the name of Barley?;

yieldData$Product0 <- as.character(yieldData$Product0);
yieldData[which(yieldData$Product0 == "Copeland"), ]$Product0 <- "Barley";
yieldData[which(yieldData$Product0 == "2 Row"), ]$Product0 <- "Barley";
yieldData[which(yieldData$Product0 == "L156H"), ]$Product0 <- "Canola";
yieldData[which(yieldData$Product0 == "L140P"), ]$Product0 <- "Canola";
yieldData[which(yieldData$Product0 == "L252"), ]$Product0 <- "Canola";
yieldData[which(yieldData$Product0 == "0"), ]$Product0 <- "Canola";
yieldData[which(yieldData$Product0 == "Canola"), ]$Product0 <- "Canola";
yieldData[which(yieldData$Product0 == "calvi"), ]$Product0 <- "Canary";
yieldData[which(yieldData$Product0 == "Greenland"), ]$Product0 <- "Lentils";
yieldData[which(yieldData$Product0 == "Greenstar"), ]$Product0 <- "Lentils";
yieldData[which(yieldData$Product0 == "Green Star"), ]$Product0 <- "Lentils";
yieldData[which(yieldData$Product0 == "Brigade"), ]$Product0 <- "Wheat";
yieldData$Product0 <- as.factor(yieldData$Product0);
summary(yieldData$Product0);

yieldData$ProcYear <- as.factor(yieldData$ProcYear);
summary(yieldData$ProcYear);

levels(yieldData$Machine);
##[1] "1H0S670SKC0746032" "1H0S670SPC0747535" "S670 2"            "1"                 "2"    
yieldData$Machine <- as.character(yieldData$Machine);
yieldData[which(yieldData$Machine == "1H0S670SKC0746032"), ]$Machine <- "1";
yieldData[which(yieldData$Machine == "1"), ]$Machine <- "1";
yieldData[which(yieldData$Machine == "1H0S670SPC0747535"), ]$Machine <- "2";
yieldData[which(yieldData$Machine == "S670 2"), ]$Machine <- "2";
yieldData[which(yieldData$Machine == "2"), ]$Machine <- "2";
yieldData$Machine <- as.factor(yieldData$Machine);
summary(yieldData$Machine);

summary(yieldData$DryYield);
#yieldDataConstraint <- yieldData[which(yieldData$DryYield <= 110), ]; #remove points with DryYield above 110 bu/acre;
#yieldDataConstraint <- yieldDataConstraint[which(yieldDataConstraint$DryYield >= 5), ]; #remove points with DryYield below 5 bu/acre;
yieldDataConstraint <- yieldData;
#summary(yieldDataConstraint$DryYield);
#str(yieldDataConstraint);
colnames(yieldDataConstraint) <- c("Elevation", "DryYield", "Product", "Machine", "ProcYear", "Harv_Date", "longitude", "latitude");
save(yieldDataConstraint, file=paste("yieldDataConstraint2015to2018", ".RData", sep=""));
#--2--2--2--2--2--2--2--2--2--#;


#--3--3--3--3--3--3--3--3--3--#;
load(paste("yieldDataConstraint2015to2018", ".RData", sep=""));
str(yieldDataConstraint);

###below loads and merges boundary shapefiles from 2015 fields to assign random effects;
setwd("/Users/tylerpittman/Farm/agYieldProject/data/boundariesSeedingCleanTopcon2015");
shps <- dir("/Users/tylerpittman/Farm/agYieldProject/data/boundariesSeedingCleanTopcon2015", "*.shp");
combinedFields <- do.call(rbind, lapply(shps, rgdal::readOGR));
##plot(combinedFields);
##combinedFields@data;
##str(combinedFields@data);

setwd("/Users/tylerpittman/Farm/agYieldProject/data");
proj4string(combinedFields) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0");
writeOGR(combinedFields, ".", "combinedFields", overwrite_layer=TRUE, driver="ESRI Shapefile");
combinedFields <- readOGR(dsn=".", layer=paste("combinedFields", sep=""));
##str(combinedFields);

coordinates(yieldDataConstraint) <- c("longitude", "latitude");
##str(yieldDataConstraint);
proj4string(yieldDataConstraint) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0");
##str(yieldDataConstraint);

#plot(combinedFields);
#plot(yieldDataConstraint, col="red", add=TRUE);

yieldDataConstraint$Field <- over(yieldDataConstraint, combinedFields)$SP_ID_1; ##does point in polygon calculation;
setwd("/Users/tylerpittman/Farm/agYieldProject/data");
save(yieldDataConstraint, file=paste("yieldDataConstraint2015to2018", ".RData", sep=""));
#--3--3--3--3--3--3--3--3--3--#;


#--4--4--4--4--4--4--4--4--4--#;
##load(paste("yieldDataConstraint2015to2018", ".RData", sep=""));
##str(yieldDataConstraint);
#--4--4--4--4--4--4--4--4--4--#;


#--5--5--5--5--5--5--5--5--5--#;
### ---------- tank1 2015---------- ###;
setwd("/Users/tylerpittman/Farm/seedShapefiles/2015");
shps <- dir("/Users/tylerpittman/Farm/seedShapefiles/2015", "*ank1_proj_2015.shp");
tank1_2015 <- do.call(rbind, lapply(shps, rgdal::readOGR));
#plot(tank1_2015);
#tank1_2015@data;
tank1_2015@data$Quintile_M <- as.numeric(as.character(tank1_2015@data$Quintile_M));
levels(tank1_2015@data$PRODUCT); 
##[1] "41-0-0-4"      "CDC Greenland" "CDC Greenstar";
tank1_2015@data$Seed_Rate <- 0;
tank1_2015@data$Ino_Rate <- 0;
tank1_2015@data$N_Rate <- 0;
tank1_2015@data$P_Rate <- 0;
tank1_2015@data$K_Rate <- 0;
tank1_2015@data$S_Rate <- 0;
tank1_2015@data[which(tank1_2015@data$PRODUCT == "41-0-0-4"), ]$Seed_Rate <- 0*tank1_2015@data[which(tank1_2015@data$PRODUCT == "41-0-0-4"), ]$Quintile_M;
tank1_2015@data[which(tank1_2015@data$PRODUCT == "41-0-0-4"), ]$Ino_Rate <- 0*tank1_2015@data[which(tank1_2015@data$PRODUCT == "41-0-0-4"), ]$Quintile_M;
tank1_2015@data[which(tank1_2015@data$PRODUCT == "41-0-0-4"), ]$N_Rate <- .41*tank1_2015@data[which(tank1_2015@data$PRODUCT == "41-0-0-4"), ]$Quintile_M;
tank1_2015@data[which(tank1_2015@data$PRODUCT == "41-0-0-4"), ]$P_Rate <- 0*tank1_2015@data[which(tank1_2015@data$PRODUCT == "41-0-0-4"), ]$Quintile_M;
tank1_2015@data[which(tank1_2015@data$PRODUCT == "41-0-0-4"), ]$K_Rate <- 0*tank1_2015@data[which(tank1_2015@data$PRODUCT == "41-0-0-4"), ]$Quintile_M;
tank1_2015@data[which(tank1_2015@data$PRODUCT == "41-0-0-4"), ]$S_Rate <- .04*tank1_2015@data[which(tank1_2015@data$PRODUCT == "41-0-0-4"), ]$Quintile_M;
tank1_2015@data[which(tank1_2015@data$PRODUCT == "CDC Greenland"), ]$Seed_Rate <- 1*tank1_2015@data[which(tank1_2015@data$PRODUCT == "CDC Greenland"), ]$Quintile_M;
tank1_2015@data[which(tank1_2015@data$PRODUCT == "CDC Greenland"), ]$Ino_Rate <- 0*tank1_2015@data[which(tank1_2015@data$PRODUCT == "CDC Greenland"), ]$Quintile_M;
tank1_2015@data[which(tank1_2015@data$PRODUCT == "CDC Greenland"), ]$N_Rate <- 0*tank1_2015@data[which(tank1_2015@data$PRODUCT == "CDC Greenland"), ]$Quintile_M;
tank1_2015@data[which(tank1_2015@data$PRODUCT == "CDC Greenland"), ]$P_Rate <- 0*tank1_2015@data[which(tank1_2015@data$PRODUCT == "CDC Greenland"), ]$Quintile_M;
tank1_2015@data[which(tank1_2015@data$PRODUCT == "CDC Greenland"), ]$K_Rate <- 0*tank1_2015@data[which(tank1_2015@data$PRODUCT == "CDC Greenland"), ]$Quintile_M;
tank1_2015@data[which(tank1_2015@data$PRODUCT == "CDC Greenland"), ]$S_Rate <- 0*tank1_2015@data[which(tank1_2015@data$PRODUCT == "CDC Greenland"), ]$Quintile_M;
tank1_2015@data[which(tank1_2015@data$PRODUCT == "CDC Greenstar"), ]$Seed_Rate <- 1*tank1_2015@data[which(tank1_2015@data$PRODUCT == "CDC Greenstar"), ]$Quintile_M;
tank1_2015@data[which(tank1_2015@data$PRODUCT == "CDC Greenstar"), ]$Ino_Rate <- 0*tank1_2015@data[which(tank1_2015@data$PRODUCT == "CDC Greenstar"), ]$Quintile_M;
tank1_2015@data[which(tank1_2015@data$PRODUCT == "CDC Greenstar"), ]$N_Rate <- 0*tank1_2015@data[which(tank1_2015@data$PRODUCT == "CDC Greenstar"), ]$Quintile_M;
tank1_2015@data[which(tank1_2015@data$PRODUCT == "CDC Greenstar"), ]$P_Rate <- 0*tank1_2015@data[which(tank1_2015@data$PRODUCT == "CDC Greenstar"), ]$Quintile_M;
tank1_2015@data[which(tank1_2015@data$PRODUCT == "CDC Greenstar"), ]$K_Rate <- 0*tank1_2015@data[which(tank1_2015@data$PRODUCT == "CDC Greenstar"), ]$Quintile_M;
tank1_2015@data[which(tank1_2015@data$PRODUCT == "CDC Greenstar"), ]$S_Rate <- 0*tank1_2015@data[which(tank1_2015@data$PRODUCT == "CDC Greenstar"), ]$Quintile_M;
setwd("/Users/tylerpittman/Farm/agYieldProject/data/seedShapefiles/2015");
writeOGR(tank1_2015, ".", "tank1_2015", overwrite_layer=TRUE, driver="ESRI Shapefile");
tank1_2015 <- tank1_2015[, c("SP_ID", "SP_ID_1", "SeedDate", "Seed_Rate", "Ino_Rate", "N_Rate", "P_Rate", "K_Rate", "S_Rate")];
projected <- "+proj=longlat";
proj4string(tank1_2015) <- CRS("+proj=utm +zone=12 +ellps=GRS80 +units=m +no_defs");
tank1_2015 <- spTransform(tank1_2015, CRS(projected));
shift.xy <- c(6, 0); #elide longitude 6 degrees to the east;
### update the geometry with elide arguments
tank1_2015 <- elide(tank1_2015, shift = shift.xy);
projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0";
proj4string(tank1_2015) <- CRS(projection);
writeOGR(tank1_2015, ".", "tank1_ll_2015", overwrite_layer=TRUE, driver="ESRI Shapefile");


### ---------- tank2 2015---------- ###;
setwd("/Users/tylerpittman/Farm/seedShapefiles/2015");
shps <- dir("/Users/tylerpittman/Farm/seedShapefiles/2015", "*ank2_proj_2015.shp");
tank2_2015 <- do.call(rbind, lapply(shps, rgdal::readOGR));
#plot(tank2_2015);
#tank2_2015@data;
tank2_2015@data$Quintile_M <- as.numeric(as.character(tank2_2015@data$Quintile_M));
levels(tank2_2015@data$PRODUCT); 
##[1] "21-0-0-24" "8-37-16-0"
tank2_2015@data$Seed_Rate <- 0;
tank2_2015@data$Ino_Rate <- 0;
tank2_2015@data$N_Rate <- 0;
tank2_2015@data$P_Rate <- 0;
tank2_2015@data$K_Rate <- 0;
tank2_2015@data$S_Rate <- 0;
tank2_2015@data[which(tank2_2015@data$PRODUCT == "21-0-0-24"), ]$Seed_Rate <- 0*tank2_2015@data[which(tank2_2015@data$PRODUCT == "21-0-0-24"), ]$Quintile_M;
tank2_2015@data[which(tank2_2015@data$PRODUCT == "21-0-0-24"), ]$Ino_Rate <- 0*tank2_2015@data[which(tank2_2015@data$PRODUCT == "21-0-0-24"), ]$Quintile_M;
tank2_2015@data[which(tank2_2015@data$PRODUCT == "21-0-0-24"), ]$N_Rate <- .21*tank2_2015@data[which(tank2_2015@data$PRODUCT == "21-0-0-24"), ]$Quintile_M;
tank2_2015@data[which(tank2_2015@data$PRODUCT == "21-0-0-24"), ]$P_Rate <- 0*tank2_2015@data[which(tank2_2015@data$PRODUCT == "21-0-0-24"), ]$Quintile_M;
tank2_2015@data[which(tank2_2015@data$PRODUCT == "21-0-0-24"), ]$K_Rate <- 0*tank2_2015@data[which(tank2_2015@data$PRODUCT == "21-0-0-24"), ]$Quintile_M;
tank2_2015@data[which(tank2_2015@data$PRODUCT == "21-0-0-24"), ]$S_Rate <- .24*tank2_2015@data[which(tank2_2015@data$PRODUCT == "21-0-0-24"), ]$Quintile_M;
tank2_2015@data[which(tank2_2015@data$PRODUCT == "8-37-16-0"), ]$Seed_Rate <- 0*tank2_2015@data[which(tank2_2015@data$PRODUCT == "8-37-16-0"), ]$Quintile_M;
tank2_2015@data[which(tank2_2015@data$PRODUCT == "8-37-16-0"), ]$Ino_Rate <- 0*tank2_2015@data[which(tank2_2015@data$PRODUCT == "8-37-16-0"), ]$Quintile_M;
tank2_2015@data[which(tank2_2015@data$PRODUCT == "8-37-16-0"), ]$N_Rate <- .08*tank2_2015@data[which(tank2_2015@data$PRODUCT == "8-37-16-0"), ]$Quintile_M;
tank2_2015@data[which(tank2_2015@data$PRODUCT == "8-37-16-0"), ]$P_Rate <- .37*tank2_2015@data[which(tank2_2015@data$PRODUCT == "8-37-16-0"), ]$Quintile_M;
tank2_2015@data[which(tank2_2015@data$PRODUCT == "8-37-16-0"), ]$K_Rate <- .16*tank2_2015@data[which(tank2_2015@data$PRODUCT == "8-37-16-0"), ]$Quintile_M;
tank2_2015@data[which(tank2_2015@data$PRODUCT == "8-37-16-0"), ]$S_Rate <- 0*tank2_2015@data[which(tank2_2015@data$PRODUCT == "8-37-16-0"), ]$Quintile_M;
setwd("/Users/tylerpittman/Farm/agYieldProject/data/seedShapefiles/2015");
writeOGR(tank2_2015, ".", "tank2_2015", overwrite_layer=TRUE, driver="ESRI Shapefile");
tank2_2015 <- tank2_2015[, c("SP_ID", "SP_ID_1", "SeedDate", "Seed_Rate", "Ino_Rate", "N_Rate", "P_Rate", "K_Rate", "S_Rate")];
projected <- "+proj=longlat";
proj4string(tank2_2015) <- CRS("+proj=utm +zone=12 +ellps=GRS80 +units=m +no_defs");
tank2_2015 <- spTransform(tank2_2015, CRS(projected));
shift.xy <- c(6, 0); #elide longitude 6 degrees to the east;
### update the geometry with elide arguments
tank2_2015 <- elide(tank2_2015, shift = shift.xy);
projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0";
proj4string(tank2_2015) <- CRS(projection);
writeOGR(tank2_2015, ".", "tank2_ll_2015", overwrite_layer=TRUE, driver="ESRI Shapefile");


### ---------- tank3 2015---------- ###;
setwd("/Users/tylerpittman/Farm/seedShapefiles/2015");
shps <- dir("/Users/tylerpittman/Farm/seedShapefiles/2015", "*ank3_proj_2015.shp");
tank3_2015 <- do.call(rbind, lapply(shps, rgdal::readOGR));
#plot(tank3_2015);
#tank3_2015@data;
tank3_2015@data$Quintile_M <- as.numeric(as.character(tank3_2015@data$Quintile_M));
levels(tank3_2015@data$PRODUCT); 
##[1] "TagTeam"
tank3_2015@data$Seed_Rate <- 0;
tank3_2015@data$Ino_Rate <- 0;
tank3_2015@data$N_Rate <- 0;
tank3_2015@data$P_Rate <- 0;
tank3_2015@data$K_Rate <- 0;
tank3_2015@data$S_Rate <- 0;
tank3_2015@data[which(tank3_2015@data$PRODUCT == "TagTeam"), ]$Seed_Rate <- 0*tank3_2015@data[which(tank3_2015@data$PRODUCT == "TagTeam"), ]$Quintile_M;
tank3_2015@data[which(tank3_2015@data$PRODUCT == "TagTeam"), ]$Ino_Rate <- 1*tank3_2015@data[which(tank3_2015@data$PRODUCT == "TagTeam"), ]$Quintile_M;
tank3_2015@data[which(tank3_2015@data$PRODUCT == "TagTeam"), ]$N_Rate <- 0*tank3_2015@data[which(tank3_2015@data$PRODUCT == "TagTeam"), ]$Quintile_M;
tank3_2015@data[which(tank3_2015@data$PRODUCT == "TagTeam"), ]$P_Rate <- 0*tank3_2015@data[which(tank3_2015@data$PRODUCT == "TagTeam"), ]$Quintile_M;
tank3_2015@data[which(tank3_2015@data$PRODUCT == "TagTeam"), ]$K_Rate <- 0*tank3_2015@data[which(tank3_2015@data$PRODUCT == "TagTeam"), ]$Quintile_M;
tank3_2015@data[which(tank3_2015@data$PRODUCT == "TagTeam"), ]$S_Rate <- 0*tank3_2015@data[which(tank3_2015@data$PRODUCT == "TagTeam"), ]$Quintile_M;
setwd("/Users/tylerpittman/Farm/agYieldProject/data/seedShapefiles/2015");
writeOGR(tank3_2015, ".", "tank3_2015", overwrite_layer=TRUE, driver="ESRI Shapefile");
tank3_2015 <- tank3_2015[, c("SP_ID", "SP_ID_1", "SeedDate", "Seed_Rate", "Ino_Rate", "N_Rate", "P_Rate", "K_Rate", "S_Rate")];
projected <- "+proj=longlat";
proj4string(tank3_2015) <- CRS("+proj=utm +zone=12 +ellps=GRS80 +units=m +no_defs");
tank3_2015 <- spTransform(tank3_2015, CRS(projected));
shift.xy <- c(6, 0); #elide longitude 6 degrees to the east;
### update the geometry with elide arguments
tank3_2015 <- elide(tank3_2015, shift = shift.xy);
projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0";
proj4string(tank3_2015) <- CRS(projection);
writeOGR(tank3_2015, ".", "tank3_ll_2015", overwrite_layer=TRUE, driver="ESRI Shapefile");


### ---------- tank4 2015---------- ###;
setwd("/Users/tylerpittman/Farm/seedShapefiles/2015");
shps <- dir("/Users/tylerpittman/Farm/seedShapefiles/2015", "*ank4_proj_2015.shp");
tank4_2015 <- do.call(rbind, lapply(shps, rgdal::readOGR));
#plot(tank4_2015);
#tank4_2015@data;
tank4_2015@data$Quintile_M <- as.numeric(as.character(tank4_2015@data$Quintile_M));
levels(tank4_2015@data$PRODUCT); 
##"11-52-0-0"  "AC Brigade"
tank4_2015@data$Seed_Rate <- 0;
tank4_2015@data$Ino_Rate <- 0;
tank4_2015@data$N_Rate <- 0;
tank4_2015@data$P_Rate <- 0;
tank4_2015@data$K_Rate <- 0;
tank4_2015@data$S_Rate <- 0;
tank4_2015@data[which(tank4_2015@data$PRODUCT == "11-52-0-0"), ]$Seed_Rate <- 0*tank4_2015@data[which(tank4_2015@data$PRODUCT == "11-52-0-0"), ]$Quintile_M;
tank4_2015@data[which(tank4_2015@data$PRODUCT == "11-52-0-0"), ]$Ino_Rate <- 0*tank4_2015@data[which(tank4_2015@data$PRODUCT == "11-52-0-0"), ]$Quintile_M;
tank4_2015@data[which(tank4_2015@data$PRODUCT == "11-52-0-0"), ]$N_Rate <- .11*tank4_2015@data[which(tank4_2015@data$PRODUCT == "11-52-0-0"), ]$Quintile_M;
tank4_2015@data[which(tank4_2015@data$PRODUCT == "11-52-0-0"), ]$P_Rate <- .52*tank4_2015@data[which(tank4_2015@data$PRODUCT == "11-52-0-0"), ]$Quintile_M;
tank4_2015@data[which(tank4_2015@data$PRODUCT == "11-52-0-0"), ]$K_Rate <- 0*tank4_2015@data[which(tank4_2015@data$PRODUCT == "11-52-0-0"), ]$Quintile_M;
tank4_2015@data[which(tank4_2015@data$PRODUCT == "11-52-0-0"), ]$S_Rate <- 0*tank4_2015@data[which(tank4_2015@data$PRODUCT == "11-52-0-0"), ]$Quintile_M;
tank4_2015@data[which(tank4_2015@data$PRODUCT == "AC Brigade"), ]$Seed_Rate <- 1*tank4_2015@data[which(tank4_2015@data$PRODUCT == "AC Brigade"), ]$Quintile_M;
tank4_2015@data[which(tank4_2015@data$PRODUCT == "AC Brigade"), ]$Ino_Rate <- 0*tank4_2015@data[which(tank4_2015@data$PRODUCT == "AC Brigade"), ]$Quintile_M;
tank4_2015@data[which(tank4_2015@data$PRODUCT == "AC Brigade"), ]$N_Rate <- 0*tank4_2015@data[which(tank4_2015@data$PRODUCT == "AC Brigade"), ]$Quintile_M;
tank4_2015@data[which(tank4_2015@data$PRODUCT == "AC Brigade"), ]$P_Rate <- 0*tank4_2015@data[which(tank4_2015@data$PRODUCT == "AC Brigade"), ]$Quintile_M;
tank4_2015@data[which(tank4_2015@data$PRODUCT == "AC Brigade"), ]$K_Rate <- 0*tank4_2015@data[which(tank4_2015@data$PRODUCT == "AC Brigade"), ]$Quintile_M;
tank4_2015@data[which(tank4_2015@data$PRODUCT == "AC Brigade"), ]$S_Rate <- 0*tank4_2015@data[which(tank4_2015@data$PRODUCT == "AC Brigade"), ]$Quintile_M;
setwd("/Users/tylerpittman/Farm/agYieldProject/data/seedShapefiles/2015");
writeOGR(tank4_2015, ".", "tank4_2015", overwrite_layer=TRUE, driver="ESRI Shapefile");
tank4_2015 <- tank4_2015[, c("SP_ID", "SP_ID_1", "SeedDate", "Seed_Rate", "Ino_Rate", "N_Rate", "P_Rate", "K_Rate", "S_Rate")];
projected <- "+proj=longlat";
proj4string(tank4_2015) <- CRS("+proj=utm +zone=12 +ellps=GRS80 +units=m +no_defs");
tank4_2015 <- spTransform(tank4_2015, CRS(projected));
shift.xy <- c(6, 0); #elide longitude 6 degrees to the east;
### update the geometry with elide arguments
tank4_2015 <- elide(tank4_2015, shift = shift.xy);
projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0";
proj4string(tank4_2015) <- CRS(projection);
writeOGR(tank4_2015, ".", "tank4_ll_2015", overwrite_layer=TRUE, driver="ESRI Shapefile");


### ---------- tank5 2015---------- ###;
setwd("/Users/tylerpittman/Farm/seedShapefiles/2015");
shps <- dir("/Users/tylerpittman/Farm/seedShapefiles/2015", "*ank5_proj_2015.shp");
tank5_2015 <- do.call(rbind, lapply(shps, rgdal::readOGR));
#plot(tank5_2015);
#tank5_2015@data;
tank5_2015@data$Quintile_M <- as.numeric(as.character(tank5_2015@data$Quintile_M));
levels(tank5_2015@data$PRODUCT); 
##[1] "L156H"
tank5_2015@data$Seed_Rate <- 0;
tank5_2015@data$Ino_Rate <- 0;
tank5_2015@data$N_Rate <- 0;
tank5_2015@data$P_Rate <- 0;
tank5_2015@data$K_Rate <- 0;
tank5_2015@data$S_Rate <- 0;
tank5_2015@data[which(tank5_2015@data$PRODUCT == "L156H"), ]$Seed_Rate <- 1*tank5_2015@data[which(tank5_2015@data$PRODUCT == "L156H"), ]$Quintile_M;
tank5_2015@data[which(tank5_2015@data$PRODUCT == "L156H"), ]$Ino_Rate <- 0*tank5_2015@data[which(tank5_2015@data$PRODUCT == "L156H"), ]$Quintile_M;
tank5_2015@data[which(tank5_2015@data$PRODUCT == "L156H"), ]$N_Rate <- 0*tank5_2015@data[which(tank5_2015@data$PRODUCT == "L156H"), ]$Quintile_M;
tank5_2015@data[which(tank5_2015@data$PRODUCT == "L156H"), ]$P_Rate <- 0*tank5_2015@data[which(tank5_2015@data$PRODUCT == "L156H"), ]$Quintile_M;
tank5_2015@data[which(tank5_2015@data$PRODUCT == "L156H"), ]$K_Rate <- 0*tank5_2015@data[which(tank5_2015@data$PRODUCT == "L156H"), ]$Quintile_M;
tank5_2015@data[which(tank5_2015@data$PRODUCT == "L156H"), ]$S_Rate <- 0*tank5_2015@data[which(tank5_2015@data$PRODUCT == "L156H"), ]$Quintile_M;
setwd("/Users/tylerpittman/Farm/agYieldProject/data/seedShapefiles/2015");
writeOGR(tank5_2015, ".", "tank5_2015", overwrite_layer=TRUE, driver="ESRI Shapefile");
tank5_2015 <- tank5_2015[, c("SP_ID", "SP_ID_1", "SeedDate", "Seed_Rate", "Ino_Rate", "N_Rate", "P_Rate", "K_Rate", "S_Rate")];
projected <- "+proj=longlat";
proj4string(tank5_2015) <- CRS("+proj=utm +zone=12 +ellps=GRS80 +units=m +no_defs");
tank5_2015 <- spTransform(tank5_2015, CRS(projected));
shift.xy <- c(6, 0); #elide longitude 6 degrees to the east;
### update the geometry with elide arguments
tank5_2015 <- elide(tank5_2015, shift = shift.xy);
projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0";
proj4string(tank5_2015) <- CRS(projection);
writeOGR(tank5_2015, ".", "tank5_ll_2015", overwrite_layer=TRUE, driver="ESRI Shapefile");


### ---------- tank1 2016---------- ###;
setwd("/Users/tylerpittman/Farm/seedShapefiles/2016");
shps <- dir("/Users/tylerpittman/Farm/seedShapefiles/2016", "*ank1_proj_2016.shp");
tank1_2016 <- do.call(rbind, lapply(shps, rgdal::readOGR));
#plot(tank1_2016);
#tank1_2016@data;
tank1_2016@data$Quintile_M <- as.numeric(as.character(tank1_2016@data$Quintile_M));
levels(tank1_2016@data$PRODUCT); 
##[1] "46-0-0-0"      "CDC Greenland" "CDC Greenstar"
tank1_2016@data$Seed_Rate <- 0;
tank1_2016@data$Ino_Rate <- 0;
tank1_2016@data$N_Rate <- 0;
tank1_2016@data$P_Rate <- 0;
tank1_2016@data$K_Rate <- 0;
tank1_2016@data$S_Rate <- 0;
tank1_2016@data[which(tank1_2016@data$PRODUCT == "46-0-0-0"), ]$Seed_Rate <- 0*tank1_2016@data[which(tank1_2016@data$PRODUCT == "46-0-0-0"), ]$Quintile_M;
tank1_2016@data[which(tank1_2016@data$PRODUCT == "46-0-0-0"), ]$Ino_Rate <- 0*tank1_2016@data[which(tank1_2016@data$PRODUCT == "46-0-0-0"), ]$Quintile_M;
tank1_2016@data[which(tank1_2016@data$PRODUCT == "46-0-0-0"), ]$N_Rate <- .46*tank1_2016@data[which(tank1_2016@data$PRODUCT == "46-0-0-0"), ]$Quintile_M;
tank1_2016@data[which(tank1_2016@data$PRODUCT == "46-0-0-0"), ]$P_Rate <- 0*tank1_2016@data[which(tank1_2016@data$PRODUCT == "46-0-0-0"), ]$Quintile_M;
tank1_2016@data[which(tank1_2016@data$PRODUCT == "46-0-0-0"), ]$K_Rate <- 0*tank1_2016@data[which(tank1_2016@data$PRODUCT == "46-0-0-0"), ]$Quintile_M;
tank1_2016@data[which(tank1_2016@data$PRODUCT == "46-0-0-0"), ]$S_Rate <- 0*tank1_2016@data[which(tank1_2016@data$PRODUCT == "46-0-0-0"), ]$Quintile_M;
tank1_2016@data[which(tank1_2016@data$PRODUCT == "CDC Greenland"), ]$Seed_Rate <- 1*tank1_2016@data[which(tank1_2016@data$PRODUCT == "CDC Greenland"), ]$Quintile_M;
tank1_2016@data[which(tank1_2016@data$PRODUCT == "CDC Greenland"), ]$Ino_Rate <- 0*tank1_2016@data[which(tank1_2016@data$PRODUCT == "CDC Greenland"), ]$Quintile_M;
tank1_2016@data[which(tank1_2016@data$PRODUCT == "CDC Greenland"), ]$N_Rate <- 0*tank1_2016@data[which(tank1_2016@data$PRODUCT == "CDC Greenland"), ]$Quintile_M;
tank1_2016@data[which(tank1_2016@data$PRODUCT == "CDC Greenland"), ]$P_Rate <- 0*tank1_2016@data[which(tank1_2016@data$PRODUCT == "CDC Greenland"), ]$Quintile_M;
tank1_2016@data[which(tank1_2016@data$PRODUCT == "CDC Greenland"), ]$K_Rate <- 0*tank1_2016@data[which(tank1_2016@data$PRODUCT == "CDC Greenland"), ]$Quintile_M;
tank1_2016@data[which(tank1_2016@data$PRODUCT == "CDC Greenland"), ]$S_Rate <- 0*tank1_2016@data[which(tank1_2016@data$PRODUCT == "CDC Greenland"), ]$Quintile_M;
tank1_2016@data[which(tank1_2016@data$PRODUCT == "CDC Greenstar"), ]$Seed_Rate <- 1*tank1_2016@data[which(tank1_2016@data$PRODUCT == "CDC Greenstar"), ]$Quintile_M;
tank1_2016@data[which(tank1_2016@data$PRODUCT == "CDC Greenstar"), ]$Ino_Rate <- 0*tank1_2016@data[which(tank1_2016@data$PRODUCT == "CDC Greenstar"), ]$Quintile_M;
tank1_2016@data[which(tank1_2016@data$PRODUCT == "CDC Greenstar"), ]$N_Rate <- 0*tank1_2016@data[which(tank1_2016@data$PRODUCT == "CDC Greenstar"), ]$Quintile_M;
tank1_2016@data[which(tank1_2016@data$PRODUCT == "CDC Greenstar"), ]$P_Rate <- 0*tank1_2016@data[which(tank1_2016@data$PRODUCT == "CDC Greenstar"), ]$Quintile_M;
tank1_2016@data[which(tank1_2016@data$PRODUCT == "CDC Greenstar"), ]$K_Rate <- 0*tank1_2016@data[which(tank1_2016@data$PRODUCT == "CDC Greenstar"), ]$Quintile_M;
tank1_2016@data[which(tank1_2016@data$PRODUCT == "CDC Greenstar"), ]$S_Rate <- 0*tank1_2016@data[which(tank1_2016@data$PRODUCT == "CDC Greenstar"), ]$Quintile_M;
setwd("/Users/tylerpittman/Farm/agYieldProject/data/seedShapefiles/2016");
writeOGR(tank1_2016, ".", "tank1_2016", overwrite_layer=TRUE, driver="ESRI Shapefile");
tank1_2016 <- tank1_2016[, c("SP_ID", "SP_ID_1", "SeedDate", "Seed_Rate", "Ino_Rate", "N_Rate", "P_Rate", "K_Rate", "S_Rate")];
projected <- "+proj=longlat";
proj4string(tank1_2016) <- CRS("+proj=utm +zone=12 +ellps=GRS80 +units=m +no_defs");
tank1_2016 <- spTransform(tank1_2016, CRS(projected));
shift.xy <- c(6, 0); #elide longitude 6 degrees to the east;
### update the geometry with elide arguments
tank1_2016 <- elide(tank1_2016, shift = shift.xy);
projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0";
proj4string(tank1_2016) <- CRS(projection);
writeOGR(tank1_2016, ".", "tank1_ll_2016", overwrite_layer=TRUE, driver="ESRI Shapefile");


### ---------- tank2 2016---------- ###;
setwd("/Users/tylerpittman/Farm/seedShapefiles/2016");
shps <- dir("/Users/tylerpittman/Farm/seedShapefiles/2016", "*ank2_proj_2016.shp");
tank2_2016 <- do.call(rbind, lapply(shps, rgdal::readOGR));
#plot(tank2_2016);
#tank2_2016@data;
tank2_2016@data$Quintile_M <- as.numeric(as.character(tank2_2016@data$Quintile_M));
levels(tank2_2016@data$PRODUCT); 
##[1] "8-38-16-0" "CDC Calvi" 
tank2_2016@data$Seed_Rate <- 0;
tank2_2016@data$Ino_Rate <- 0;
tank2_2016@data$N_Rate <- 0;
tank2_2016@data$P_Rate <- 0;
tank2_2016@data$K_Rate <- 0;
tank2_2016@data$S_Rate <- 0;
tank2_2016@data[which(tank2_2016@data$PRODUCT == "8-38-16-0"), ]$Seed_Rate <- 0*tank2_2016@data[which(tank2_2016@data$PRODUCT == "8-38-16-0"), ]$Quintile_M;
tank2_2016@data[which(tank2_2016@data$PRODUCT == "8-38-16-0"), ]$Ino_Rate <- 0*tank2_2016@data[which(tank2_2016@data$PRODUCT == "8-38-16-0"), ]$Quintile_M;
tank2_2016@data[which(tank2_2016@data$PRODUCT == "8-38-16-0"), ]$N_Rate <- .08*tank2_2016@data[which(tank2_2016@data$PRODUCT == "8-38-16-0"), ]$Quintile_M;
tank2_2016@data[which(tank2_2016@data$PRODUCT == "8-38-16-0"), ]$P_Rate <- .38*tank2_2016@data[which(tank2_2016@data$PRODUCT == "8-38-16-0"), ]$Quintile_M;
tank2_2016@data[which(tank2_2016@data$PRODUCT == "8-38-16-0"), ]$K_Rate <- .16*tank2_2016@data[which(tank2_2016@data$PRODUCT == "8-38-16-0"), ]$Quintile_M;
tank2_2016@data[which(tank2_2016@data$PRODUCT == "8-38-16-0"), ]$S_Rate <- 0*tank2_2016@data[which(tank2_2016@data$PRODUCT == "8-38-16-0"), ]$Quintile_M;
tank2_2016@data[which(tank2_2016@data$PRODUCT == "CDC Calvi"), ]$Seed_Rate <- 1*tank2_2016@data[which(tank2_2016@data$PRODUCT == "CDC Calvi"), ]$Quintile_M;
tank2_2016@data[which(tank2_2016@data$PRODUCT == "CDC Calvi"), ]$Ino_Rate <- 0*tank2_2016@data[which(tank2_2016@data$PRODUCT == "CDC Calvi"), ]$Quintile_M;
tank2_2016@data[which(tank2_2016@data$PRODUCT == "CDC Calvi"), ]$N_Rate <- 0*tank2_2016@data[which(tank2_2016@data$PRODUCT == "CDC Calvi"), ]$Quintile_M;
tank2_2016@data[which(tank2_2016@data$PRODUCT == "CDC Calvi"), ]$P_Rate <- 0*tank2_2016@data[which(tank2_2016@data$PRODUCT == "CDC Calvi"), ]$Quintile_M;
tank2_2016@data[which(tank2_2016@data$PRODUCT == "CDC Calvi"), ]$K_Rate <- 0*tank2_2016@data[which(tank2_2016@data$PRODUCT == "CDC Calvi"), ]$Quintile_M;
tank2_2016@data[which(tank2_2016@data$PRODUCT == "CDC Calvi"), ]$S_Rate <- 0*tank2_2016@data[which(tank2_2016@data$PRODUCT == "CDC Calvi"), ]$Quintile_M;
setwd("/Users/tylerpittman/Farm/agYieldProject/data/seedShapefiles/2016");
writeOGR(tank2_2016, ".", "tank2_2016", overwrite_layer=TRUE, driver="ESRI Shapefile");
tank2_2016 <- tank2_2016[, c("SP_ID", "SP_ID_1", "SeedDate", "Seed_Rate", "Ino_Rate", "N_Rate", "P_Rate", "K_Rate", "S_Rate")];
projected <- "+proj=longlat";
proj4string(tank2_2016) <- CRS("+proj=utm +zone=12 +ellps=GRS80 +units=m +no_defs");
tank2_2016 <- spTransform(tank2_2016, CRS(projected));
shift.xy <- c(6, 0); #elide longitude 6 degrees to the east;
### update the geometry with elide arguments
tank2_2016 <- elide(tank2_2016, shift = shift.xy);
projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0";
proj4string(tank2_2016) <- CRS(projection);
writeOGR(tank2_2016, ".", "tank2_ll_2016", overwrite_layer=TRUE, driver="ESRI Shapefile");


### ---------- tank3 2016---------- ###;
setwd("/Users/tylerpittman/Farm/seedShapefiles/2016");
shps <- dir("/Users/tylerpittman/Farm/seedShapefiles/2016", "*ank3_proj_2016.shp");
tank3_2016 <- do.call(rbind, lapply(shps, rgdal::readOGR));
###### NO INFORMATION FOR INOCULANT RATE IN TANK 3 FOR 2016 SEEDING, USED TANK 5 INSTEAD ######
##plot(tank3_2016);
##tank3_2016@data;
#tank3_2016@data$Quintile_M <- as.numeric(as.character(tank3_2016@data$Quintile_M));
#levels(tank3_2016@data$PRODUCT); 
###[1] "TagTeam"
#tank3_2016@data$Seed_Rate <- 0;
#tank3_2016@data$Ino_Rate <- 0;
#tank3_2016@data$N_Rate <- 0;
#tank3_2016@data$P_Rate <- 0;
#tank3_2016@data$K_Rate <- 0;
#tank3_2016@data$S_Rate <- 0;
#tank3_2016@data[which(tank3_2016@data$PRODUCT == "TagTeam"), ]$Seed_Rate <- 0*tank3_2016@data[which(tank3_2016@data$PRODUCT == "TagTeam"), ]$Quintile_M;
#tank3_2016@data[which(tank3_2016@data$PRODUCT == "TagTeam"), ]$Ino_Rate <- 1*tank3_2016@data[which(tank3_2016@data$PRODUCT == "TagTeam"), ]$Quintile_M;
#tank3_2016@data[which(tank3_2016@data$PRODUCT == "TagTeam"), ]$N_Rate <- 0*tank3_2016@data[which(tank3_2016@data$PRODUCT == "TagTeam"), ]$Quintile_M;
#tank3_2016@data[which(tank3_2016@data$PRODUCT == "TagTeam"), ]$P_Rate <- 0*tank3_2016@data[which(tank3_2016@data$PRODUCT == "TagTeam"), ]$Quintile_M;
#tank3_2016@data[which(tank3_2016@data$PRODUCT == "TagTeam"), ]$K_Rate <- 0*tank3_2016@data[which(tank3_2016@data$PRODUCT == "TagTeam"), ]$Quintile_M;
#tank3_2016@data[which(tank3_2016@data$PRODUCT == "TagTeam"), ]$S_Rate <- 0*tank3_2016@data[which(tank3_2016@data$PRODUCT == "TagTeam"), ]$Quintile_M;
#setwd("/Users/tylerpittman/Farm/agYieldProject/data/seedShapefiles/2016");
#writeOGR(tank3_2016, ".", "tank3_2016", overwrite_layer=TRUE, driver="ESRI Shapefile");
#tank3_2016 <- tank3_2016[, c("SP_ID", "SP_ID_1", "SeedDate", "Seed_Rate", "Ino_Rate", "N_Rate", "P_Rate", "K_Rate", "S_Rate")];
#projected <- "+proj=longlat";
#proj4string(tank3_2016) <- CRS("+proj=utm +zone=12 +ellps=GRS80 +units=m +no_defs");
#tank3_2016 <- spTransform(tank3_2016, CRS(projected));
#shift.xy <- c(6, 0); #elide longitude 6 degrees to the east;
### update the geometry with elide arguments
#tank3_2016 <- elide(tank3_2016, shift = shift.xy);
#projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0";
#proj4string(tank3_2016) <- CRS(projection);
#writeOGR(tank3_2016, ".", "tank3_ll_2016", overwrite_layer=TRUE, driver="ESRI Shapefile");


### ---------- tank4 2016---------- ###;
setwd("/Users/tylerpittman/Farm/seedShapefiles/2016");
shps <- dir("/Users/tylerpittman/Farm/seedShapefiles/2016", "*ank4_proj_2016.shp");
tank4_2016 <- do.call(rbind, lapply(shps, rgdal::readOGR));
#plot(tank4_2016);
#tank4_2016@data;
tank4_2016@data$Quintile_M <- as.numeric(as.character(tank4_2016@data$Quintile_M));
levels(tank4_2016@data$PRODUCT); 
##[1] "AC Brigade" "8-38-16-0" 
tank4_2016@data$Seed_Rate <- 0;
tank4_2016@data$Ino_Rate <- 0;
tank4_2016@data$N_Rate <- 0;
tank4_2016@data$P_Rate <- 0;
tank4_2016@data$K_Rate <- 0;
tank4_2016@data$S_Rate <- 0;
tank4_2016@data[which(tank4_2016@data$PRODUCT == "8-38-16-0"), ]$Seed_Rate <- 0*tank4_2016@data[which(tank4_2016@data$PRODUCT == "8-38-16-0"), ]$Quintile_M;
tank4_2016@data[which(tank4_2016@data$PRODUCT == "8-38-16-0"), ]$Ino_Rate <- 0*tank4_2016@data[which(tank4_2016@data$PRODUCT == "8-38-16-0"), ]$Quintile_M;
tank4_2016@data[which(tank4_2016@data$PRODUCT == "8-38-16-0"), ]$N_Rate <- .08*tank4_2016@data[which(tank4_2016@data$PRODUCT == "8-38-16-0"), ]$Quintile_M;
tank4_2016@data[which(tank4_2016@data$PRODUCT == "8-38-16-0"), ]$P_Rate <- .38*tank4_2016@data[which(tank4_2016@data$PRODUCT == "8-38-16-0"), ]$Quintile_M;
tank4_2016@data[which(tank4_2016@data$PRODUCT == "8-38-16-0"), ]$K_Rate <- .16*tank4_2016@data[which(tank4_2016@data$PRODUCT == "8-38-16-0"), ]$Quintile_M;
tank4_2016@data[which(tank4_2016@data$PRODUCT == "8-38-16-0"), ]$S_Rate <- 0*tank4_2016@data[which(tank4_2016@data$PRODUCT == "8-38-16-0"), ]$Quintile_M;
tank4_2016@data[which(tank4_2016@data$PRODUCT == "AC Brigade"), ]$Seed_Rate <- 1*tank4_2016@data[which(tank4_2016@data$PRODUCT == "AC Brigade"), ]$Quintile_M;
tank4_2016@data[which(tank4_2016@data$PRODUCT == "AC Brigade"), ]$Ino_Rate <- 0*tank4_2016@data[which(tank4_2016@data$PRODUCT == "AC Brigade"), ]$Quintile_M;
tank4_2016@data[which(tank4_2016@data$PRODUCT == "AC Brigade"), ]$N_Rate <- 0*tank4_2016@data[which(tank4_2016@data$PRODUCT == "AC Brigade"), ]$Quintile_M;
tank4_2016@data[which(tank4_2016@data$PRODUCT == "AC Brigade"), ]$P_Rate <- 0*tank4_2016@data[which(tank4_2016@data$PRODUCT == "AC Brigade"), ]$Quintile_M;
tank4_2016@data[which(tank4_2016@data$PRODUCT == "AC Brigade"), ]$K_Rate <- 0*tank4_2016@data[which(tank4_2016@data$PRODUCT == "AC Brigade"), ]$Quintile_M;
tank4_2016@data[which(tank4_2016@data$PRODUCT == "AC Brigade"), ]$S_Rate <- 0*tank4_2016@data[which(tank4_2016@data$PRODUCT == "AC Brigade"), ]$Quintile_M;
setwd("/Users/tylerpittman/Farm/agYieldProject/data/seedShapefiles/2016");
writeOGR(tank4_2016, ".", "tank4_2016", overwrite_layer=TRUE, driver="ESRI Shapefile");
tank4_2016 <- tank4_2016[, c("SP_ID", "SP_ID_1", "SeedDate", "Seed_Rate", "Ino_Rate", "N_Rate", "P_Rate", "K_Rate", "S_Rate")];
projected <- "+proj=longlat";
proj4string(tank4_2016) <- CRS("+proj=utm +zone=12 +ellps=GRS80 +units=m +no_defs");
tank4_2016 <- spTransform(tank4_2016, CRS(projected));
shift.xy <- c(6, 0); #elide longitude 6 degrees to the east;
### update the geometry with elide arguments
tank4_2016 <- elide(tank4_2016, shift = shift.xy);
projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0";
proj4string(tank4_2016) <- CRS(projection);
writeOGR(tank4_2016, ".", "tank4_ll_2016", overwrite_layer=TRUE, driver="ESRI Shapefile");


### ---------- tank5 2016---------- ###;
setwd("/Users/tylerpittman/Farm/seedShapefiles/2016");
shps <- dir("/Users/tylerpittman/Farm/seedShapefiles/2016", "*ank5_proj_2016.shp");
tank5_2016 <- do.call(rbind, lapply(shps, rgdal::readOGR));
#plot(tank5_2016);
#tank5_2016@data;
tank5_2016@data$Quintile_M <- as.numeric(as.character(tank5_2016@data$Quintile_M));
levels(tank5_2016@data$PRODUCT); 
##[1] "TagTeam"
tank5_2016@data$Seed_Rate <- 0;
tank5_2016@data$Ino_Rate <- 0;
tank5_2016@data$N_Rate <- 0;
tank5_2016@data$P_Rate <- 0;
tank5_2016@data$K_Rate <- 0;
tank5_2016@data$S_Rate <- 0;
tank5_2016@data[which(tank5_2016@data$PRODUCT == "TagTeam"), ]$Seed_Rate <- 0*tank5_2016@data[which(tank5_2016@data$PRODUCT == "TagTeam"), ]$Quintile_M;
tank5_2016@data[which(tank5_2016@data$PRODUCT == "TagTeam"), ]$Ino_Rate <- 1*tank5_2016@data[which(tank5_2016@data$PRODUCT == "TagTeam"), ]$Quintile_M;
tank5_2016@data[which(tank5_2016@data$PRODUCT == "TagTeam"), ]$N_Rate <- 0*tank5_2016@data[which(tank5_2016@data$PRODUCT == "TagTeam"), ]$Quintile_M;
tank5_2016@data[which(tank5_2016@data$PRODUCT == "TagTeam"), ]$P_Rate <- 0*tank5_2016@data[which(tank5_2016@data$PRODUCT == "TagTeam"), ]$Quintile_M;
tank5_2016@data[which(tank5_2016@data$PRODUCT == "TagTeam"), ]$K_Rate <- 0*tank5_2016@data[which(tank5_2016@data$PRODUCT == "TagTeam"), ]$Quintile_M;
tank5_2016@data[which(tank5_2016@data$PRODUCT == "TagTeam"), ]$S_Rate <- 0*tank5_2016@data[which(tank5_2016@data$PRODUCT == "TagTeam"), ]$Quintile_M;
setwd("/Users/tylerpittman/Farm/agYieldProject/data/seedShapefiles/2016");
writeOGR(tank5_2016, ".", "tank5_2016", overwrite_layer=TRUE, driver="ESRI Shapefile");
tank5_2016 <- tank5_2016[, c("SP_ID", "SP_ID_1", "SeedDate", "Seed_Rate", "Ino_Rate", "N_Rate", "P_Rate", "K_Rate", "S_Rate")];
projected <- "+proj=longlat";
proj4string(tank5_2016) <- CRS("+proj=utm +zone=12 +ellps=GRS80 +units=m +no_defs");
tank5_2016 <- spTransform(tank5_2016, CRS(projected));
shift.xy <- c(6, 0); #elide longitude 6 degrees to the east;
### update the geometry with elide arguments
tank5_2016 <- elide(tank5_2016, shift = shift.xy);
projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0";
proj4string(tank5_2016) <- CRS(projection);
writeOGR(tank5_2016, ".", "tank5_ll_2016", overwrite_layer=TRUE, driver="ESRI Shapefile");


### ---------- tank1 2017---------- ###;
setwd("/Users/tylerpittman/Farm/seedShapefiles/2017");
shps <- dir("/Users/tylerpittman/Farm/seedShapefiles/2017", "*ank1_proj_2017.shp");
tank1_2017 <- do.call(rbind, lapply(shps, rgdal::readOGR));
#plot(tank1_2017);
#tank1_2017@data;
tank1_2017@data$Quintile_M <- as.numeric(as.character(tank1_2017@data$Quintile_M));
levels(tank1_2017@data$PRODUCT); 
##[1] "CDC Greenland" "46-0-0-0"   
tank1_2017@data$Seed_Rate <- 0;
tank1_2017@data$Ino_Rate <- 0;
tank1_2017@data$N_Rate <- 0;
tank1_2017@data$P_Rate <- 0;
tank1_2017@data$K_Rate <- 0;
tank1_2017@data$S_Rate <- 0;
tank1_2017@data[which(tank1_2017@data$PRODUCT == "46-0-0-0"), ]$Seed_Rate <- 0*tank1_2017@data[which(tank1_2017@data$PRODUCT == "46-0-0-0"), ]$Quintile_M;
tank1_2017@data[which(tank1_2017@data$PRODUCT == "46-0-0-0"), ]$Ino_Rate <- 0*tank1_2017@data[which(tank1_2017@data$PRODUCT == "46-0-0-0"), ]$Quintile_M;
tank1_2017@data[which(tank1_2017@data$PRODUCT == "46-0-0-0"), ]$N_Rate <- .46*tank1_2017@data[which(tank1_2017@data$PRODUCT == "46-0-0-0"), ]$Quintile_M;
tank1_2017@data[which(tank1_2017@data$PRODUCT == "46-0-0-0"), ]$P_Rate <- 0*tank1_2017@data[which(tank1_2017@data$PRODUCT == "46-0-0-0"), ]$Quintile_M;
tank1_2017@data[which(tank1_2017@data$PRODUCT == "46-0-0-0"), ]$K_Rate <- 0*tank1_2017@data[which(tank1_2017@data$PRODUCT == "46-0-0-0"), ]$Quintile_M;
tank1_2017@data[which(tank1_2017@data$PRODUCT == "46-0-0-0"), ]$S_Rate <- 0*tank1_2017@data[which(tank1_2017@data$PRODUCT == "46-0-0-0"), ]$Quintile_M;
tank1_2017@data[which(tank1_2017@data$PRODUCT == "CDC Greenland"), ]$Seed_Rate <- 1*tank1_2017@data[which(tank1_2017@data$PRODUCT == "CDC Greenland"), ]$Quintile_M;
tank1_2017@data[which(tank1_2017@data$PRODUCT == "CDC Greenland"), ]$Ino_Rate <- 0*tank1_2017@data[which(tank1_2017@data$PRODUCT == "CDC Greenland"), ]$Quintile_M;
tank1_2017@data[which(tank1_2017@data$PRODUCT == "CDC Greenland"), ]$N_Rate <- 0*tank1_2017@data[which(tank1_2017@data$PRODUCT == "CDC Greenland"), ]$Quintile_M;
tank1_2017@data[which(tank1_2017@data$PRODUCT == "CDC Greenland"), ]$P_Rate <- 0*tank1_2017@data[which(tank1_2017@data$PRODUCT == "CDC Greenland"), ]$Quintile_M;
tank1_2017@data[which(tank1_2017@data$PRODUCT == "CDC Greenland"), ]$K_Rate <- 0*tank1_2017@data[which(tank1_2017@data$PRODUCT == "CDC Greenland"), ]$Quintile_M;
tank1_2017@data[which(tank1_2017@data$PRODUCT == "CDC Greenland"), ]$S_Rate <- 0*tank1_2017@data[which(tank1_2017@data$PRODUCT == "CDC Greenland"), ]$Quintile_M;
setwd("/Users/tylerpittman/Farm/agYieldProject/data/seedShapefiles/2017");
writeOGR(tank1_2017, ".", "tank1_2017", overwrite_layer=TRUE, driver="ESRI Shapefile");
tank1_2017 <- tank1_2017[, c("SP_ID", "SP_ID_1", "SeedDate", "Seed_Rate", "Ino_Rate", "N_Rate", "P_Rate", "K_Rate", "S_Rate")];
projected <- "+proj=longlat";
proj4string(tank1_2017) <- CRS("+proj=utm +zone=12 +ellps=GRS80 +units=m +no_defs");
tank1_2017 <- spTransform(tank1_2017, CRS(projected));
shift.xy <- c(6, 0); #elide longitude 6 degrees to the east;
### update the geometry with elide arguments
tank1_2017 <- elide(tank1_2017, shift = shift.xy);
projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0";
proj4string(tank1_2017) <- CRS(projection);
writeOGR(tank1_2017, ".", "tank1_ll_2017", overwrite_layer=TRUE, driver="ESRI Shapefile");


### ---------- tank2 2017---------- ###;
setwd("/Users/tylerpittman/Farm/seedShapefiles/2017");
shps <- dir("/Users/tylerpittman/Farm/seedShapefiles/2017", "*ank2_proj_2017.shp");
tank2_2017 <- do.call(rbind, lapply(shps, rgdal::readOGR));
#plot(tank2_2017);
#tank2_2017@data;
tank2_2017@data$Quintile_M <- as.numeric(as.character(tank2_2017@data$Quintile_M));
levels(tank2_2017@data$PRODUCT); 
##[1] "8-38-16-0" "CDC Calvi"
#tank2_2017@data[which(tank2_2017@data$PRODUCT == "na"), ];
tank2_2017@data$Seed_Rate <- 0;
tank2_2017@data$Ino_Rate <- 0;
tank2_2017@data$N_Rate <- 0;
tank2_2017@data$P_Rate <- 0;
tank2_2017@data$K_Rate <- 0;
tank2_2017@data$S_Rate <- 0;
tank2_2017@data[which(tank2_2017@data$PRODUCT == "8-38-16-0"), ]$Seed_Rate <- 0*tank2_2017@data[which(tank2_2017@data$PRODUCT == "8-38-16-0"), ]$Quintile_M;
tank2_2017@data[which(tank2_2017@data$PRODUCT == "8-38-16-0"), ]$Ino_Rate <- 0*tank2_2017@data[which(tank2_2017@data$PRODUCT == "8-38-16-0"), ]$Quintile_M;
tank2_2017@data[which(tank2_2017@data$PRODUCT == "8-38-16-0"), ]$N_Rate <- .08*tank2_2017@data[which(tank2_2017@data$PRODUCT == "8-38-16-0"), ]$Quintile_M;
tank2_2017@data[which(tank2_2017@data$PRODUCT == "8-38-16-0"), ]$P_Rate <- .38*tank2_2017@data[which(tank2_2017@data$PRODUCT == "8-38-16-0"), ]$Quintile_M;
tank2_2017@data[which(tank2_2017@data$PRODUCT == "8-38-16-0"), ]$K_Rate <- .16*tank2_2017@data[which(tank2_2017@data$PRODUCT == "8-38-16-0"), ]$Quintile_M;
tank2_2017@data[which(tank2_2017@data$PRODUCT == "8-38-16-0"), ]$S_Rate <- 0*tank2_2017@data[which(tank2_2017@data$PRODUCT == "8-38-16-0"), ]$Quintile_M;
tank2_2017@data[which(tank2_2017@data$PRODUCT == "CDC Calvi"), ]$Seed_Rate <- 1*tank2_2017@data[which(tank2_2017@data$PRODUCT == "CDC Calvi"), ]$Quintile_M;
tank2_2017@data[which(tank2_2017@data$PRODUCT == "CDC Calvi"), ]$Ino_Rate <- 0*tank2_2017@data[which(tank2_2017@data$PRODUCT == "CDC Calvi"), ]$Quintile_M;
tank2_2017@data[which(tank2_2017@data$PRODUCT == "CDC Calvi"), ]$N_Rate <- 0*tank2_2017@data[which(tank2_2017@data$PRODUCT == "CDC Calvi"), ]$Quintile_M;
tank2_2017@data[which(tank2_2017@data$PRODUCT == "CDC Calvi"), ]$P_Rate <- 0*tank2_2017@data[which(tank2_2017@data$PRODUCT == "CDC Calvi"), ]$Quintile_M;
tank2_2017@data[which(tank2_2017@data$PRODUCT == "CDC Calvi"), ]$K_Rate <- 0*tank2_2017@data[which(tank2_2017@data$PRODUCT == "CDC Calvi"), ]$Quintile_M;
tank2_2017@data[which(tank2_2017@data$PRODUCT == "CDC Calvi"), ]$S_Rate <- 0*tank2_2017@data[which(tank2_2017@data$PRODUCT == "CDC Calvi"), ]$Quintile_M;
setwd("/Users/tylerpittman/Farm/agYieldProject/data/seedShapefiles/2017");
writeOGR(tank2_2017, ".", "tank2_2017", overwrite_layer=TRUE, driver="ESRI Shapefile");
tank2_2017 <- tank2_2017[, c("SP_ID", "SP_ID_1", "SeedDate", "Seed_Rate", "Ino_Rate", "N_Rate", "P_Rate", "K_Rate", "S_Rate")];
projected <- "+proj=longlat";
proj4string(tank2_2017) <- CRS("+proj=utm +zone=12 +ellps=GRS80 +units=m +no_defs");
tank2_2017 <- spTransform(tank2_2017, CRS(projected));
shift.xy <- c(6, 0); #elide longitude 6 degrees to the east;
### update the geometry with elide arguments
tank2_2017 <- elide(tank2_2017, shift = shift.xy);
projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0";
proj4string(tank2_2017) <- CRS(projection);
writeOGR(tank2_2017, ".", "tank2_ll_2017", overwrite_layer=TRUE, driver="ESRI Shapefile");


### ---------- tank3 2017---------- ###;
setwd("/Users/tylerpittman/Farm/seedShapefiles/2017");
shps <- dir("/Users/tylerpittman/Farm/seedShapefiles/2017", "*ank3_proj_2017.shp");
tank3_2017 <- do.call(rbind, lapply(shps, rgdal::readOGR));
#plot(tank3_2017);
#tank3_2017@data;
tank3_2017@data$Quintile_M <- as.numeric(as.character(tank3_2017@data$Quintile_M));
levels(tank3_2017@data$PRODUCT); 
##[1] "TagTeam"
tank3_2017@data$Seed_Rate <- 0;
tank3_2017@data$Ino_Rate <- 0;
tank3_2017@data$N_Rate <- 0;
tank3_2017@data$P_Rate <- 0;
tank3_2017@data$K_Rate <- 0;
tank3_2017@data$S_Rate <- 0;
tank3_2017@data[which(tank3_2017@data$PRODUCT == "TagTeam"), ]$Seed_Rate <- 0*tank3_2017@data[which(tank3_2017@data$PRODUCT == "TagTeam"), ]$Quintile_M;
tank3_2017@data[which(tank3_2017@data$PRODUCT == "TagTeam"), ]$Ino_Rate <- 1*tank3_2017@data[which(tank3_2017@data$PRODUCT == "TagTeam"), ]$Quintile_M;
tank3_2017@data[which(tank3_2017@data$PRODUCT == "TagTeam"), ]$N_Rate <- 0*tank3_2017@data[which(tank3_2017@data$PRODUCT == "TagTeam"), ]$Quintile_M;
tank3_2017@data[which(tank3_2017@data$PRODUCT == "TagTeam"), ]$P_Rate <- 0*tank3_2017@data[which(tank3_2017@data$PRODUCT == "TagTeam"), ]$Quintile_M;
tank3_2017@data[which(tank3_2017@data$PRODUCT == "TagTeam"), ]$K_Rate <- 0*tank3_2017@data[which(tank3_2017@data$PRODUCT == "TagTeam"), ]$Quintile_M;
tank3_2017@data[which(tank3_2017@data$PRODUCT == "TagTeam"), ]$S_Rate <- 0*tank3_2017@data[which(tank3_2017@data$PRODUCT == "TagTeam"), ]$Quintile_M;
setwd("/Users/tylerpittman/Farm/agYieldProject/data/seedShapefiles/2017");
writeOGR(tank3_2017, ".", "tank3_2017", overwrite_layer=TRUE, driver="ESRI Shapefile");
tank3_2017 <- tank3_2017[, c("SP_ID", "SP_ID_1", "SeedDate", "Seed_Rate", "Ino_Rate", "N_Rate", "P_Rate", "K_Rate", "S_Rate")];
projected <- "+proj=longlat";
proj4string(tank3_2017) <- CRS("+proj=utm +zone=12 +ellps=GRS80 +units=m +no_defs");
tank3_2017 <- spTransform(tank3_2017, CRS(projected));
shift.xy <- c(6, 0); #elide longitude 6 degrees to the east;
### update the geometry with elide arguments
tank3_2017 <- elide(tank3_2017, shift = shift.xy);
projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0";
proj4string(tank3_2017) <- CRS(projection);
writeOGR(tank3_2017, ".", "tank3_ll_2017", overwrite_layer=TRUE, driver="ESRI Shapefile");


### ---------- tank4 2017---------- ###;
setwd("/Users/tylerpittman/Farm/seedShapefiles/2017");
shps <- dir("/Users/tylerpittman/Farm/seedShapefiles/2017", "*ank4_proj_2017.shp");
tank4_2017 <- do.call(rbind, lapply(shps, rgdal::readOGR));
#plot(tank4_2017);
#tank4_2017@data;
tank4_2017@data$Quintile_M <- as.numeric(as.character(tank4_2017@data$Quintile_M));
levels(tank4_2017@data$PRODUCT); 
##[1] "8-38-16-0"    "AC Brigade"   "21-0-0-24"    "CDC Copeland"
tank4_2017@data$Seed_Rate <- 0;
tank4_2017@data$Ino_Rate <- 0;
tank4_2017@data$N_Rate <- 0;
tank4_2017@data$P_Rate <- 0;
tank4_2017@data$K_Rate <- 0;
tank4_2017@data$S_Rate <- 0;
tank4_2017@data[which(tank4_2017@data$PRODUCT == "8-38-16-0"), ]$Seed_Rate <- 0*tank4_2017@data[which(tank4_2017@data$PRODUCT == "8-38-16-0"), ]$Quintile_M;
tank4_2017@data[which(tank4_2017@data$PRODUCT == "8-38-16-0"), ]$Ino_Rate <- 0*tank4_2017@data[which(tank4_2017@data$PRODUCT == "8-38-16-0"), ]$Quintile_M;
tank4_2017@data[which(tank4_2017@data$PRODUCT == "8-38-16-0"), ]$N_Rate <- .08*tank4_2017@data[which(tank4_2017@data$PRODUCT == "8-38-16-0"), ]$Quintile_M;
tank4_2017@data[which(tank4_2017@data$PRODUCT == "8-38-16-0"), ]$P_Rate <- .38*tank4_2017@data[which(tank4_2017@data$PRODUCT == "8-38-16-0"), ]$Quintile_M;
tank4_2017@data[which(tank4_2017@data$PRODUCT == "8-38-16-0"), ]$K_Rate <- .16*tank4_2017@data[which(tank4_2017@data$PRODUCT == "8-38-16-0"), ]$Quintile_M;
tank4_2017@data[which(tank4_2017@data$PRODUCT == "8-38-16-0"), ]$S_Rate <- 0*tank4_2017@data[which(tank4_2017@data$PRODUCT == "8-38-16-0"), ]$Quintile_M;
tank4_2017@data[which(tank4_2017@data$PRODUCT == "AC Brigade"), ]$Seed_Rate <- 1*tank4_2017@data[which(tank4_2017@data$PRODUCT == "AC Brigade"), ]$Quintile_M;
tank4_2017@data[which(tank4_2017@data$PRODUCT == "AC Brigade"), ]$Ino_Rate <- 0*tank4_2017@data[which(tank4_2017@data$PRODUCT == "AC Brigade"), ]$Quintile_M;
tank4_2017@data[which(tank4_2017@data$PRODUCT == "AC Brigade"), ]$N_Rate <- 0*tank4_2017@data[which(tank4_2017@data$PRODUCT == "AC Brigade"), ]$Quintile_M;
tank4_2017@data[which(tank4_2017@data$PRODUCT == "AC Brigade"), ]$P_Rate <- 0*tank4_2017@data[which(tank4_2017@data$PRODUCT == "AC Brigade"), ]$Quintile_M;
tank4_2017@data[which(tank4_2017@data$PRODUCT == "AC Brigade"), ]$K_Rate <- 0*tank4_2017@data[which(tank4_2017@data$PRODUCT == "AC Brigade"), ]$Quintile_M;
tank4_2017@data[which(tank4_2017@data$PRODUCT == "AC Brigade"), ]$S_Rate <- 0*tank4_2017@data[which(tank4_2017@data$PRODUCT == "AC Brigade"), ]$Quintile_M;
tank4_2017@data[which(tank4_2017@data$PRODUCT == "21-0-0-24"), ]$Seed_Rate <- 0*tank4_2017@data[which(tank4_2017@data$PRODUCT == "21-0-0-24"), ]$Quintile_M;
tank4_2017@data[which(tank4_2017@data$PRODUCT == "21-0-0-24"), ]$Ino_Rate <- 0*tank4_2017@data[which(tank4_2017@data$PRODUCT == "21-0-0-24"), ]$Quintile_M;
tank4_2017@data[which(tank4_2017@data$PRODUCT == "21-0-0-24"), ]$N_Rate <- .21*tank4_2017@data[which(tank4_2017@data$PRODUCT == "21-0-0-24"), ]$Quintile_M;
tank4_2017@data[which(tank4_2017@data$PRODUCT == "21-0-0-24"), ]$P_Rate <- 0*tank4_2017@data[which(tank4_2017@data$PRODUCT == "21-0-0-24"), ]$Quintile_M;
tank4_2017@data[which(tank4_2017@data$PRODUCT == "21-0-0-24"), ]$K_Rate <- 0*tank4_2017@data[which(tank4_2017@data$PRODUCT == "21-0-0-24"), ]$Quintile_M;
tank4_2017@data[which(tank4_2017@data$PRODUCT == "21-0-0-24"), ]$S_Rate <- .24*tank4_2017@data[which(tank4_2017@data$PRODUCT == "21-0-0-24"), ]$Quintile_M;
tank4_2017@data[which(tank4_2017@data$PRODUCT == "CDC Copeland"), ]$Seed_Rate <- 1*tank4_2017@data[which(tank4_2017@data$PRODUCT == "CDC Copeland"), ]$Quintile_M;
tank4_2017@data[which(tank4_2017@data$PRODUCT == "CDC Copeland"), ]$Ino_Rate <- 0*tank4_2017@data[which(tank4_2017@data$PRODUCT == "CDC Copeland"), ]$Quintile_M;
tank4_2017@data[which(tank4_2017@data$PRODUCT == "CDC Copeland"), ]$N_Rate <- 0*tank4_2017@data[which(tank4_2017@data$PRODUCT == "CDC Copeland"), ]$Quintile_M;
tank4_2017@data[which(tank4_2017@data$PRODUCT == "CDC Copeland"), ]$P_Rate <- 0*tank4_2017@data[which(tank4_2017@data$PRODUCT == "CDC Copeland"), ]$Quintile_M;
tank4_2017@data[which(tank4_2017@data$PRODUCT == "CDC Copeland"), ]$K_Rate <- 0*tank4_2017@data[which(tank4_2017@data$PRODUCT == "CDC Copeland"), ]$Quintile_M;
tank4_2017@data[which(tank4_2017@data$PRODUCT == "CDC Copeland"), ]$S_Rate <- 0*tank4_2017@data[which(tank4_2017@data$PRODUCT == "CDC Copeland"), ]$Quintile_M;
setwd("/Users/tylerpittman/Farm/agYieldProject/data/seedShapefiles/2017");
writeOGR(tank4_2017, ".", "tank4_2017", overwrite_layer=TRUE, driver="ESRI Shapefile");
tank4_2017 <- tank4_2017[, c("SP_ID", "SP_ID_1", "SeedDate", "Seed_Rate", "Ino_Rate", "N_Rate", "P_Rate", "K_Rate", "S_Rate")];
projected <- "+proj=longlat";
proj4string(tank4_2017) <- CRS("+proj=utm +zone=12 +ellps=GRS80 +units=m +no_defs");
tank4_2017 <- spTransform(tank4_2017, CRS(projected));
shift.xy <- c(6, 0); #elide longitude 6 degrees to the east;
### update the geometry with elide arguments
tank4_2017 <- elide(tank4_2017, shift = shift.xy);
projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0";
proj4string(tank4_2017) <- CRS(projection);
writeOGR(tank4_2017, ".", "tank4_ll_2017", overwrite_layer=TRUE, driver="ESRI Shapefile");


### ---------- tank5 2017---------- ###;
setwd("/Users/tylerpittman/Farm/seedShapefiles/2017");
shps <- dir("/Users/tylerpittman/Farm/seedShapefiles/2017", "*ank5_proj_2017.shp");
tank5_2017 <- do.call(rbind, lapply(shps, rgdal::readOGR));
#plot(tank5_2017);
#tank5_2017@data;
tank5_2017@data$Quintile_M <- as.numeric(as.character(tank5_2017@data$Quintile_M));
levels(tank5_2017@data$PRODUCT); 
##[1] "InVigor L140P"
tank5_2017@data$Seed_Rate <- 0;
tank5_2017@data$Ino_Rate <- 0;
tank5_2017@data$N_Rate <- 0;
tank5_2017@data$P_Rate <- 0;
tank5_2017@data$K_Rate <- 0;
tank5_2017@data$S_Rate <- 0;
tank5_2017@data[which(tank5_2017@data$PRODUCT == "InVigor L140P"), ]$Seed_Rate <- 1*tank5_2017@data[which(tank5_2017@data$PRODUCT == "InVigor L140P"), ]$Quintile_M;
tank5_2017@data[which(tank5_2017@data$PRODUCT == "InVigor L140P"), ]$Ino_Rate <- 0*tank5_2017@data[which(tank5_2017@data$PRODUCT == "InVigor L140P"), ]$Quintile_M;
tank5_2017@data[which(tank5_2017@data$PRODUCT == "InVigor L140P"), ]$N_Rate <- 0*tank5_2017@data[which(tank5_2017@data$PRODUCT == "InVigor L140P"), ]$Quintile_M;
tank5_2017@data[which(tank5_2017@data$PRODUCT == "InVigor L140P"), ]$P_Rate <- 0*tank5_2017@data[which(tank5_2017@data$PRODUCT == "InVigor L140P"), ]$Quintile_M;
tank5_2017@data[which(tank5_2017@data$PRODUCT == "InVigor L140P"), ]$K_Rate <- 0*tank5_2017@data[which(tank5_2017@data$PRODUCT == "InVigor L140P"), ]$Quintile_M;
tank5_2017@data[which(tank5_2017@data$PRODUCT == "InVigor L140P"), ]$S_Rate <- 0*tank5_2017@data[which(tank5_2017@data$PRODUCT == "InVigor L140P"), ]$Quintile_M;
setwd("/Users/tylerpittman/Farm/agYieldProject/data/seedShapefiles/2017");
writeOGR(tank5_2017, ".", "tank5_2017", overwrite_layer=TRUE, driver="ESRI Shapefile");
tank5_2017 <- tank5_2017[, c("SP_ID", "SP_ID_1", "SeedDate", "Seed_Rate", "Ino_Rate", "N_Rate", "P_Rate", "K_Rate", "S_Rate")];
projected <- "+proj=longlat";
proj4string(tank5_2017) <- CRS("+proj=utm +zone=12 +ellps=GRS80 +units=m +no_defs");
tank5_2017 <- spTransform(tank5_2017, CRS(projected));
shift.xy <- c(6, 0); #elide longitude 6 degrees to the east;
### update the geometry with elide arguments
tank5_2017 <- elide(tank5_2017, shift = shift.xy);
projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0";
proj4string(tank5_2017) <- CRS(projection);
writeOGR(tank5_2017, ".", "tank5_ll_2017", overwrite_layer=TRUE, driver="ESRI Shapefile");


### ---------- tank1 2018---------- ###;
setwd("/Users/tylerpittman/Farm/seedShapefiles/2018");
shps <- dir("/Users/tylerpittman/Farm/seedShapefiles/2018", "*ank1_proj_2018.shp");
tank1_2018 <- do.call(rbind, lapply(shps, rgdal::readOGR));
#plot(tank1_2018);
#tank1_2018@data;
tank1_2018@data$Quintile_M <- as.numeric(as.character(tank1_2018@data$Quintile_M));
levels(tank1_2018@data$PRODUCT); 
##[1] "46-0-0-0"      "CDC Greenland"   
tank1_2018@data$Seed_Rate <- 0;
tank1_2018@data$Ino_Rate <- 0;
tank1_2018@data$N_Rate <- 0;
tank1_2018@data$P_Rate <- 0;
tank1_2018@data$K_Rate <- 0;
tank1_2018@data$S_Rate <- 0;
tank1_2018@data[which(tank1_2018@data$PRODUCT == "46-0-0-0"), ]$Seed_Rate <- 0*tank1_2018@data[which(tank1_2018@data$PRODUCT == "46-0-0-0"), ]$Quintile_M;
tank1_2018@data[which(tank1_2018@data$PRODUCT == "46-0-0-0"), ]$Ino_Rate <- 0*tank1_2018@data[which(tank1_2018@data$PRODUCT == "46-0-0-0"), ]$Quintile_M;
tank1_2018@data[which(tank1_2018@data$PRODUCT == "46-0-0-0"), ]$N_Rate <- .46*tank1_2018@data[which(tank1_2018@data$PRODUCT == "46-0-0-0"), ]$Quintile_M;
tank1_2018@data[which(tank1_2018@data$PRODUCT == "46-0-0-0"), ]$P_Rate <- 0*tank1_2018@data[which(tank1_2018@data$PRODUCT == "46-0-0-0"), ]$Quintile_M;
tank1_2018@data[which(tank1_2018@data$PRODUCT == "46-0-0-0"), ]$K_Rate <- 0*tank1_2018@data[which(tank1_2018@data$PRODUCT == "46-0-0-0"), ]$Quintile_M;
tank1_2018@data[which(tank1_2018@data$PRODUCT == "46-0-0-0"), ]$S_Rate <- 0*tank1_2018@data[which(tank1_2018@data$PRODUCT == "46-0-0-0"), ]$Quintile_M;
tank1_2018@data[which(tank1_2018@data$PRODUCT == "CDC Greenland"), ]$Seed_Rate <- 1*tank1_2018@data[which(tank1_2018@data$PRODUCT == "CDC Greenland"), ]$Quintile_M;
tank1_2018@data[which(tank1_2018@data$PRODUCT == "CDC Greenland"), ]$Ino_Rate <- 0*tank1_2018@data[which(tank1_2018@data$PRODUCT == "CDC Greenland"), ]$Quintile_M;
tank1_2018@data[which(tank1_2018@data$PRODUCT == "CDC Greenland"), ]$N_Rate <- 0*tank1_2018@data[which(tank1_2018@data$PRODUCT == "CDC Greenland"), ]$Quintile_M;
tank1_2018@data[which(tank1_2018@data$PRODUCT == "CDC Greenland"), ]$P_Rate <- 0*tank1_2018@data[which(tank1_2018@data$PRODUCT == "CDC Greenland"), ]$Quintile_M;
tank1_2018@data[which(tank1_2018@data$PRODUCT == "CDC Greenland"), ]$K_Rate <- 0*tank1_2018@data[which(tank1_2018@data$PRODUCT == "CDC Greenland"), ]$Quintile_M;
tank1_2018@data[which(tank1_2018@data$PRODUCT == "CDC Greenland"), ]$S_Rate <- 0*tank1_2018@data[which(tank1_2018@data$PRODUCT == "CDC Greenland"), ]$Quintile_M;
setwd("/Users/tylerpittman/Farm/agYieldProject/data/seedShapefiles/2018");
writeOGR(tank1_2018, ".", "tank1_2018", overwrite_layer=TRUE, driver="ESRI Shapefile");
tank1_2018 <- tank1_2018[, c("SP_ID", "SP_ID_1", "SeedDate", "Seed_Rate", "Ino_Rate", "N_Rate", "P_Rate", "K_Rate", "S_Rate")];
projected <- "+proj=longlat";
proj4string(tank1_2018) <- CRS("+proj=utm +zone=12 +ellps=GRS80 +units=m +no_defs");
tank1_2018 <- spTransform(tank1_2018, CRS(projected));
shift.xy <- c(6, 0); #elide longitude 6 degrees to the east;
### update the geometry with elide arguments
tank1_2018 <- elide(tank1_2018, shift = shift.xy);
projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0";
proj4string(tank1_2018) <- CRS(projection);
writeOGR(tank1_2018, ".", "tank1_ll_2018", overwrite_layer=TRUE, driver="ESRI Shapefile");


### ---------- tank2 2018---------- ###;
setwd("/Users/tylerpittman/Farm/seedShapefiles/2018");
shps <- dir("/Users/tylerpittman/Farm/seedShapefiles/2018", "*ank2_proj_2018.shp");
tank2_2018 <- do.call(rbind, lapply(shps, rgdal::readOGR));
#plot(tank2_2018);
#tank2_2018@data;
tank2_2018@data$Quintile_M <- as.numeric(as.character(tank2_2018@data$Quintile_M));
levels(tank2_2018@data$PRODUCT); 
##[1] "8-38-16-0"
#tank2_2018@data[which(tank2_2018@data$PRODUCT == "na"), ];
tank2_2018@data$Seed_Rate <- 0;
tank2_2018@data$Ino_Rate <- 0;
tank2_2018@data$N_Rate <- 0;
tank2_2018@data$P_Rate <- 0;
tank2_2018@data$K_Rate <- 0;
tank2_2018@data$S_Rate <- 0;
tank2_2018@data[which(tank2_2018@data$PRODUCT == "8-38-16-0"), ]$Seed_Rate <- 0*tank2_2018@data[which(tank2_2018@data$PRODUCT == "8-38-16-0"), ]$Quintile_M;
tank2_2018@data[which(tank2_2018@data$PRODUCT == "8-38-16-0"), ]$Ino_Rate <- 0*tank2_2018@data[which(tank2_2018@data$PRODUCT == "8-38-16-0"), ]$Quintile_M;
tank2_2018@data[which(tank2_2018@data$PRODUCT == "8-38-16-0"), ]$N_Rate <- .08*tank2_2018@data[which(tank2_2018@data$PRODUCT == "8-38-16-0"), ]$Quintile_M;
tank2_2018@data[which(tank2_2018@data$PRODUCT == "8-38-16-0"), ]$P_Rate <- .38*tank2_2018@data[which(tank2_2018@data$PRODUCT == "8-38-16-0"), ]$Quintile_M;
tank2_2018@data[which(tank2_2018@data$PRODUCT == "8-38-16-0"), ]$K_Rate <- .16*tank2_2018@data[which(tank2_2018@data$PRODUCT == "8-38-16-0"), ]$Quintile_M;
tank2_2018@data[which(tank2_2018@data$PRODUCT == "8-38-16-0"), ]$S_Rate <- 0*tank2_2018@data[which(tank2_2018@data$PRODUCT == "8-38-16-0"), ]$Quintile_M;
setwd("/Users/tylerpittman/Farm/agYieldProject/data/seedShapefiles/2018");
writeOGR(tank2_2018, ".", "tank2_2018", overwrite_layer=TRUE, driver="ESRI Shapefile");
tank2_2018 <- tank2_2018[, c("SP_ID", "SP_ID_1", "SeedDate", "Seed_Rate", "Ino_Rate", "N_Rate", "P_Rate", "K_Rate", "S_Rate")];
projected <- "+proj=longlat";
proj4string(tank2_2018) <- CRS("+proj=utm +zone=12 +ellps=GRS80 +units=m +no_defs");
tank2_2018 <- spTransform(tank2_2018, CRS(projected));
shift.xy <- c(6, 0); #elide longitude 6 degrees to the east;
### update the geometry with elide arguments
tank2_2018 <- elide(tank2_2018, shift = shift.xy);
projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0";
proj4string(tank2_2018) <- CRS(projection);
writeOGR(tank2_2018, ".", "tank2_ll_2018", overwrite_layer=TRUE, driver="ESRI Shapefile");


### ---------- tank3 2018---------- ###;
setwd("/Users/tylerpittman/Farm/seedShapefiles/2018");
shps <- dir("/Users/tylerpittman/Farm/seedShapefiles/2018", "*ank3_proj_2018.shp");
tank3_2018 <- do.call(rbind, lapply(shps, rgdal::readOGR));
#plot(tank3_2018);
#tank3_2018@data;
tank3_2018@data$Quintile_M <- as.numeric(as.character(tank3_2018@data$Quintile_M));
levels(tank3_2018@data$PRODUCT); 
##[1] "TagTeam"  
tank3_2018@data$Seed_Rate <- 0;
tank3_2018@data$Ino_Rate <- 0;
tank3_2018@data$N_Rate <- 0;
tank3_2018@data$P_Rate <- 0;
tank3_2018@data$K_Rate <- 0;
tank3_2018@data$S_Rate <- 0;
tank3_2018@data[which(tank3_2018@data$PRODUCT == "TagTeam"), ]$Seed_Rate <- 0*tank3_2018@data[which(tank3_2018@data$PRODUCT == "TagTeam"), ]$Quintile_M;
tank3_2018@data[which(tank3_2018@data$PRODUCT == "TagTeam"), ]$Ino_Rate <- 1*tank3_2018@data[which(tank3_2018@data$PRODUCT == "TagTeam"), ]$Quintile_M;
tank3_2018@data[which(tank3_2018@data$PRODUCT == "TagTeam"), ]$N_Rate <- 0*tank3_2018@data[which(tank3_2018@data$PRODUCT == "TagTeam"), ]$Quintile_M;
tank3_2018@data[which(tank3_2018@data$PRODUCT == "TagTeam"), ]$P_Rate <- 0*tank3_2018@data[which(tank3_2018@data$PRODUCT == "TagTeam"), ]$Quintile_M;
tank3_2018@data[which(tank3_2018@data$PRODUCT == "TagTeam"), ]$K_Rate <- 0*tank3_2018@data[which(tank3_2018@data$PRODUCT == "TagTeam"), ]$Quintile_M;
tank3_2018@data[which(tank3_2018@data$PRODUCT == "TagTeam"), ]$S_Rate <- 0*tank3_2018@data[which(tank3_2018@data$PRODUCT == "TagTeam"), ]$Quintile_M;
setwd("/Users/tylerpittman/Farm/agYieldProject/data/seedShapefiles/2018");
writeOGR(tank3_2018, ".", "tank3_2018", overwrite_layer=TRUE, driver="ESRI Shapefile");
tank3_2018 <- tank3_2018[, c("SP_ID", "SP_ID_1", "SeedDate", "Seed_Rate", "Ino_Rate", "N_Rate", "P_Rate", "K_Rate", "S_Rate")];
projected <- "+proj=longlat";
proj4string(tank3_2018) <- CRS("+proj=utm +zone=12 +ellps=GRS80 +units=m +no_defs");
tank3_2018 <- spTransform(tank3_2018, CRS(projected));
shift.xy <- c(6, 0); #elide longitude 6 degrees to the east;
### update the geometry with elide arguments
tank3_2018 <- elide(tank3_2018, shift = shift.xy);
projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0";
proj4string(tank3_2018) <- CRS(projection);
writeOGR(tank3_2018, ".", "tank3_ll_2018", overwrite_layer=TRUE, driver="ESRI Shapefile");


### ---------- tank4 2018---------- ###;
setwd("/Users/tylerpittman/Farm/seedShapefiles/2018");
shps <- dir("/Users/tylerpittman/Farm/seedShapefiles/2018", "*ank4_proj_2018.shp");
tank4_2018 <- do.call(rbind, lapply(shps, rgdal::readOGR));
#plot(tank4_2018);
#tank4_2018@data;
tank4_2018@data$Quintile_M <- as.numeric(as.character(tank4_2018@data$Quintile_M));
levels(tank4_2018@data$PRODUCT); 
##[1] "AC Brigade" "21-0-0-24"  "8-38-16-0"
tank4_2018@data$Seed_Rate <- 0;
tank4_2018@data$Ino_Rate <- 0;
tank4_2018@data$N_Rate <- 0;
tank4_2018@data$P_Rate <- 0;
tank4_2018@data$K_Rate <- 0;
tank4_2018@data$S_Rate <- 0;
tank4_2018@data[which(tank4_2018@data$PRODUCT == "8-38-16-0"), ]$Seed_Rate <- 0*tank4_2018@data[which(tank4_2018@data$PRODUCT == "8-38-16-0"), ]$Quintile_M;
tank4_2018@data[which(tank4_2018@data$PRODUCT == "8-38-16-0"), ]$Ino_Rate <- 0*tank4_2018@data[which(tank4_2018@data$PRODUCT == "8-38-16-0"), ]$Quintile_M;
tank4_2018@data[which(tank4_2018@data$PRODUCT == "8-38-16-0"), ]$N_Rate <- .08*tank4_2018@data[which(tank4_2018@data$PRODUCT == "8-38-16-0"), ]$Quintile_M;
tank4_2018@data[which(tank4_2018@data$PRODUCT == "8-38-16-0"), ]$P_Rate <- .38*tank4_2018@data[which(tank4_2018@data$PRODUCT == "8-38-16-0"), ]$Quintile_M;
tank4_2018@data[which(tank4_2018@data$PRODUCT == "8-38-16-0"), ]$K_Rate <- .16*tank4_2018@data[which(tank4_2018@data$PRODUCT == "8-38-16-0"), ]$Quintile_M;
tank4_2018@data[which(tank4_2018@data$PRODUCT == "8-38-16-0"), ]$S_Rate <- 0*tank4_2018@data[which(tank4_2018@data$PRODUCT == "8-38-16-0"), ]$Quintile_M;
tank4_2018@data[which(tank4_2018@data$PRODUCT == "AC Brigade"), ]$Seed_Rate <- 1*tank4_2018@data[which(tank4_2018@data$PRODUCT == "AC Brigade"), ]$Quintile_M;
tank4_2018@data[which(tank4_2018@data$PRODUCT == "AC Brigade"), ]$Ino_Rate <- 0*tank4_2018@data[which(tank4_2018@data$PRODUCT == "AC Brigade"), ]$Quintile_M;
tank4_2018@data[which(tank4_2018@data$PRODUCT == "AC Brigade"), ]$N_Rate <- 0*tank4_2018@data[which(tank4_2018@data$PRODUCT == "AC Brigade"), ]$Quintile_M;
tank4_2018@data[which(tank4_2018@data$PRODUCT == "AC Brigade"), ]$P_Rate <- 0*tank4_2018@data[which(tank4_2018@data$PRODUCT == "AC Brigade"), ]$Quintile_M;
tank4_2018@data[which(tank4_2018@data$PRODUCT == "AC Brigade"), ]$K_Rate <- 0*tank4_2018@data[which(tank4_2018@data$PRODUCT == "AC Brigade"), ]$Quintile_M;
tank4_2018@data[which(tank4_2018@data$PRODUCT == "AC Brigade"), ]$S_Rate <- 0*tank4_2018@data[which(tank4_2018@data$PRODUCT == "AC Brigade"), ]$Quintile_M;
tank4_2018@data[which(tank4_2018@data$PRODUCT == "21-0-0-24"), ]$Seed_Rate <- 0*tank4_2018@data[which(tank4_2018@data$PRODUCT == "21-0-0-24"), ]$Quintile_M;
tank4_2018@data[which(tank4_2018@data$PRODUCT == "21-0-0-24"), ]$Ino_Rate <- 0*tank4_2018@data[which(tank4_2018@data$PRODUCT == "21-0-0-24"), ]$Quintile_M;
tank4_2018@data[which(tank4_2018@data$PRODUCT == "21-0-0-24"), ]$N_Rate <- .21*tank4_2018@data[which(tank4_2018@data$PRODUCT == "21-0-0-24"), ]$Quintile_M;
tank4_2018@data[which(tank4_2018@data$PRODUCT == "21-0-0-24"), ]$P_Rate <- 0*tank4_2018@data[which(tank4_2018@data$PRODUCT == "21-0-0-24"), ]$Quintile_M;
tank4_2018@data[which(tank4_2018@data$PRODUCT == "21-0-0-24"), ]$K_Rate <- 0*tank4_2018@data[which(tank4_2018@data$PRODUCT == "21-0-0-24"), ]$Quintile_M;
tank4_2018@data[which(tank4_2018@data$PRODUCT == "21-0-0-24"), ]$S_Rate <- .24*tank4_2018@data[which(tank4_2018@data$PRODUCT == "21-0-0-24"), ]$Quintile_M;
setwd("/Users/tylerpittman/Farm/agYieldProject/data/seedShapefiles/2018");
writeOGR(tank4_2018, ".", "tank4_2018", overwrite_layer=TRUE, driver="ESRI Shapefile");
tank4_2018 <- tank4_2018[, c("SP_ID", "SP_ID_1", "SeedDate", "Seed_Rate", "Ino_Rate", "N_Rate", "P_Rate", "K_Rate", "S_Rate")];
projected <- "+proj=longlat";
proj4string(tank4_2018) <- CRS("+proj=utm +zone=12 +ellps=GRS80 +units=m +no_defs");
tank4_2018 <- spTransform(tank4_2018, CRS(projected));
shift.xy <- c(6, 0); #elide longitude 6 degrees to the east;
### update the geometry with elide arguments
tank4_2018 <- elide(tank4_2018, shift = shift.xy);
projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0";
proj4string(tank4_2018) <- CRS(projection);
writeOGR(tank4_2018, ".", "tank4_ll_2018", overwrite_layer=TRUE, driver="ESRI Shapefile");


### ---------- tank5 2018---------- ###;
setwd("/Users/tylerpittman/Farm/seedShapefiles/2018");
shps <- dir("/Users/tylerpittman/Farm/seedShapefiles/2018", "*ank5_proj_2018.shp");
tank5_2018 <- do.call(rbind, lapply(shps, rgdal::readOGR));
#plot(tank5_2018);
#tank5_2018@data;
tank5_2018@data$Quintile_M <- as.numeric(as.character(tank5_2018@data$Quintile_M));
levels(tank5_2018@data$PRODUCT); 
##[1] "InVigor L233P" "InVigor L140P"
tank5_2018@data$Seed_Rate <- 0;
tank5_2018@data$Ino_Rate <- 0;
tank5_2018@data$N_Rate <- 0;
tank5_2018@data$P_Rate <- 0;
tank5_2018@data$K_Rate <- 0;
tank5_2018@data$S_Rate <- 0;
tank5_2018@data[which(tank5_2018@data$PRODUCT == "InVigor L233P"), ]$Seed_Rate <- 1*tank5_2018@data[which(tank5_2018@data$PRODUCT == "InVigor L233P"), ]$Quintile_M;
tank5_2018@data[which(tank5_2018@data$PRODUCT == "InVigor L233P"), ]$Ino_Rate <- 0*tank5_2018@data[which(tank5_2018@data$PRODUCT == "InVigor L233P"), ]$Quintile_M;
tank5_2018@data[which(tank5_2018@data$PRODUCT == "InVigor L233P"), ]$N_Rate <- 0*tank5_2018@data[which(tank5_2018@data$PRODUCT == "InVigor L233P"), ]$Quintile_M;
tank5_2018@data[which(tank5_2018@data$PRODUCT == "InVigor L233P"), ]$P_Rate <- 0*tank5_2018@data[which(tank5_2018@data$PRODUCT == "InVigor L233P"), ]$Quintile_M;
tank5_2018@data[which(tank5_2018@data$PRODUCT == "InVigor L233P"), ]$K_Rate <- 0*tank5_2018@data[which(tank5_2018@data$PRODUCT == "InVigor L233P"), ]$Quintile_M;
tank5_2018@data[which(tank5_2018@data$PRODUCT == "InVigor L233P"), ]$S_Rate <- 0*tank5_2018@data[which(tank5_2018@data$PRODUCT == "InVigor L233P"), ]$Quintile_M;
tank5_2018@data[which(tank5_2018@data$PRODUCT == "InVigor L140P"), ]$Seed_Rate <- 1*tank5_2018@data[which(tank5_2018@data$PRODUCT == "InVigor L140P"), ]$Quintile_M;
tank5_2018@data[which(tank5_2018@data$PRODUCT == "InVigor L140P"), ]$Ino_Rate <- 0*tank5_2018@data[which(tank5_2018@data$PRODUCT == "InVigor L140P"), ]$Quintile_M;
tank5_2018@data[which(tank5_2018@data$PRODUCT == "InVigor L140P"), ]$N_Rate <- 0*tank5_2018@data[which(tank5_2018@data$PRODUCT == "InVigor L140P"), ]$Quintile_M;
tank5_2018@data[which(tank5_2018@data$PRODUCT == "InVigor L140P"), ]$P_Rate <- 0*tank5_2018@data[which(tank5_2018@data$PRODUCT == "InVigor L140P"), ]$Quintile_M;
tank5_2018@data[which(tank5_2018@data$PRODUCT == "InVigor L140P"), ]$K_Rate <- 0*tank5_2018@data[which(tank5_2018@data$PRODUCT == "InVigor L140P"), ]$Quintile_M;
tank5_2018@data[which(tank5_2018@data$PRODUCT == "InVigor L140P"), ]$S_Rate <- 0*tank5_2018@data[which(tank5_2018@data$PRODUCT == "InVigor L140P"), ]$Quintile_M;
setwd("/Users/tylerpittman/Farm/agYieldProject/data/seedShapefiles/2018");
writeOGR(tank5_2018, ".", "tank5_2018", overwrite_layer=TRUE, driver="ESRI Shapefile");
tank5_2018 <- tank5_2018[, c("SP_ID", "SP_ID_1", "SeedDate", "Seed_Rate", "Ino_Rate", "N_Rate", "P_Rate", "K_Rate", "S_Rate")];
projected <- "+proj=longlat";
proj4string(tank5_2018) <- CRS("+proj=utm +zone=12 +ellps=GRS80 +units=m +no_defs");
tank5_2018 <- spTransform(tank5_2018, CRS(projected));
shift.xy <- c(6, 0); #elide longitude 6 degrees to the east;
### update the geometry with elide arguments
tank5_2018 <- elide(tank5_2018, shift = shift.xy);
projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0";
proj4string(tank5_2018) <- CRS(projection);
writeOGR(tank5_2018, ".", "tank5_ll_2018", overwrite_layer=TRUE, driver="ESRI Shapefile");
#--5--5--5--5--5--5--5--5--5--#;