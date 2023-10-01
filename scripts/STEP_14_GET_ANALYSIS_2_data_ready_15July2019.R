## Takes farm combine yield data and seeding data from 2015 to 2018
## Tyler Pittman, 1 August 2019
# Rscript --no-save /Users/tylerpittman/Farm/STEP_14_GET_ANALYSIS_2_data_ready_15July2019.R

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
library(readxl);
library(plyr); 		#gives multicolumn frequencies for tables with ddply() and join_all function;
library(dplyr);  	#need dplyr for the bind_rows function 
library(parallel);	# Does multicore processing when using mclapply WORKS on Linux!
#library(greenbrown); #has many complicated dependencies to install, for equal area parcel division of shapefiles;
#source("shape2poly.R"); # Reads in shape2poly function and others
#source("polygonizer.R"); # Reads in polygonizer function

print(citation("sp"), bibtex=T); #bibtex citation for point in polygon over() function from sp package;

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


setwd("/Users/tylerpittman/Farm/agYieldProject/data");
#save(yieldDataConstraint, file=paste("yieldDataConstraint2015to2018", ".RData", sep=""));
#---------------------------------------------------------------------------------------#;


load(paste("yieldDataConstraint2015to2018", ".RData", sep=""));
proj4string(yieldDataConstraint) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0");
str(yieldDataConstraint);
yieldDataConstraint2015 <- yieldDataConstraint[which(yieldDataConstraint$ProcYear == "2015"), ];
yieldDataConstraint2016 <- yieldDataConstraint[which(yieldDataConstraint$ProcYear == "2016"), ];
yieldDataConstraint2017 <- yieldDataConstraint[which(yieldDataConstraint$ProcYear == "2017"), ];
yieldDataConstraint2018 <- yieldDataConstraint[which(yieldDataConstraint$ProcYear == "2018"), ];


setwd("/Users/tylerpittman/Farm/agYieldProject/data/seedShapefiles/2015");
tank1_2015 <- readOGR(dsn=".", layer=paste("tank1_ll_2015", sep=""));
tank2_2015 <- readOGR(dsn=".", layer=paste("tank2_ll_2015", sep=""));
tank3_2015 <- readOGR(dsn=".", layer=paste("tank3_ll_2015", sep=""));
tank4_2015 <- readOGR(dsn=".", layer=paste("tank4_ll_2015", sep=""));
tank5_2015 <- readOGR(dsn=".", layer=paste("tank5_ll_2015", sep=""));

setwd("/Users/tylerpittman/Farm/agYieldProject/data/seedShapefiles/2016");
tank1_2016 <- readOGR(dsn=".", layer=paste("tank1_ll_2016", sep=""));
tank2_2016 <- readOGR(dsn=".", layer=paste("tank2_ll_2016", sep=""));
#tank3_2016 <- readOGR(dsn=".", layer=paste("tank3_ll_2016", sep="")); #No Tank 3 used in 2016;
tank4_2016 <- readOGR(dsn=".", layer=paste("tank4_ll_2016", sep=""));
tank5_2016 <- readOGR(dsn=".", layer=paste("tank5_ll_2016", sep=""));

setwd("/Users/tylerpittman/Farm/agYieldProject/data/seedShapefiles/2017");
tank1_2017 <- readOGR(dsn=".", layer=paste("tank1_ll_2017", sep=""));
tank2_2017 <- readOGR(dsn=".", layer=paste("tank2_ll_2017", sep=""));
tank3_2017 <- readOGR(dsn=".", layer=paste("tank3_ll_2017", sep=""));
tank4_2017 <- readOGR(dsn=".", layer=paste("tank4_ll_2017", sep=""));
tank5_2017 <- readOGR(dsn=".", layer=paste("tank5_ll_2017", sep=""));

setwd("/Users/tylerpittman/Farm/agYieldProject/data/seedShapefiles/2018");
tank1_2018 <- readOGR(dsn=".", layer=paste("tank1_ll_2018", sep=""));
tank2_2018 <- readOGR(dsn=".", layer=paste("tank2_ll_2018", sep=""));
tank3_2018 <- readOGR(dsn=".", layer=paste("tank3_ll_2018", sep=""));
tank4_2018 <- readOGR(dsn=".", layer=paste("tank4_ll_2018", sep=""));
tank5_2018 <- readOGR(dsn=".", layer=paste("tank5_ll_2018", sep=""));


#---$---$---$---$---$---$---$---#;
setwd("/Users/tylerpittman/Farm/agYieldProject/data");
#str(tank1_2015@data);
#proj4string(tank1_2015);
#proj4string(yieldDataConstraint2015);
##ERROR identicalCRS(x, y) is not TRUE
##plot(tank1_2015);
##plot(yieldDataConstraint2015, col="red", add=TRUE);
##colnames(yieldDataConstraint2015@data);



### ------------------------------- BELOW TAKE 14 HOURS;
##---- 2015;
###BELOW does point in polygon calculation;
yieldDataConstraint2015@data$Seed_Date <- over(yieldDataConstraint2015, tank1_2015)$SeedDate;

yieldDataConstraint2015@data$Seed_Rate1 <- over(yieldDataConstraint2015, tank1_2015)$Seed_Rate;
yieldDataConstraint2015@data$Seed_Rate2 <- over(yieldDataConstraint2015, tank2_2015)$Seed_Rate;
yieldDataConstraint2015@data$Seed_Rate3 <- over(yieldDataConstraint2015, tank3_2015)$Seed_Rate;
yieldDataConstraint2015@data$Seed_Rate4 <- over(yieldDataConstraint2015, tank4_2015)$Seed_Rate;
yieldDataConstraint2015@data$Seed_Rate5 <- over(yieldDataConstraint2015, tank5_2015)$Seed_Rate;

yieldDataConstraint2015@data$Ino_Rate1 <- over(yieldDataConstraint2015, tank1_2015)$Ino_Rate;
yieldDataConstraint2015@data$Ino_Rate2 <- over(yieldDataConstraint2015, tank2_2015)$Ino_Rate;
yieldDataConstraint2015@data$Ino_Rate3 <- over(yieldDataConstraint2015, tank3_2015)$Ino_Rate;
yieldDataConstraint2015@data$Ino_Rate4 <- over(yieldDataConstraint2015, tank4_2015)$Ino_Rate;
yieldDataConstraint2015@data$Ino_Rate5 <- over(yieldDataConstraint2015, tank5_2015)$Ino_Rate;

yieldDataConstraint2015@data$N_Rate1 <- over(yieldDataConstraint2015, tank1_2015)$N_Rate;
yieldDataConstraint2015@data$N_Rate2 <- over(yieldDataConstraint2015, tank2_2015)$N_Rate;
yieldDataConstraint2015@data$N_Rate3 <- over(yieldDataConstraint2015, tank3_2015)$N_Rate;
yieldDataConstraint2015@data$N_Rate4 <- over(yieldDataConstraint2015, tank4_2015)$N_Rate;
yieldDataConstraint2015@data$N_Rate5 <- over(yieldDataConstraint2015, tank5_2015)$N_Rate;

yieldDataConstraint2015@data$P_Rate1 <- over(yieldDataConstraint2015, tank1_2015)$P_Rate;
yieldDataConstraint2015@data$P_Rate2 <- over(yieldDataConstraint2015, tank2_2015)$P_Rate;
yieldDataConstraint2015@data$P_Rate3 <- over(yieldDataConstraint2015, tank3_2015)$P_Rate;
yieldDataConstraint2015@data$P_Rate4 <- over(yieldDataConstraint2015, tank4_2015)$P_Rate;
yieldDataConstraint2015@data$P_Rate5 <- over(yieldDataConstraint2015, tank5_2015)$P_Rate;

yieldDataConstraint2015@data$K_Rate1 <- over(yieldDataConstraint2015, tank1_2015)$K_Rate;
yieldDataConstraint2015@data$K_Rate2 <- over(yieldDataConstraint2015, tank2_2015)$K_Rate;
yieldDataConstraint2015@data$K_Rate3 <- over(yieldDataConstraint2015, tank3_2015)$K_Rate;
yieldDataConstraint2015@data$K_Rate4 <- over(yieldDataConstraint2015, tank4_2015)$K_Rate;
yieldDataConstraint2015@data$K_Rate5 <- over(yieldDataConstraint2015, tank5_2015)$K_Rate;

yieldDataConstraint2015@data$S_Rate1 <- over(yieldDataConstraint2015, tank1_2015)$S_Rate;
yieldDataConstraint2015@data$S_Rate2 <- over(yieldDataConstraint2015, tank2_2015)$S_Rate;
yieldDataConstraint2015@data$S_Rate3 <- over(yieldDataConstraint2015, tank3_2015)$S_Rate;
yieldDataConstraint2015@data$S_Rate4 <- over(yieldDataConstraint2015, tank4_2015)$S_Rate;
yieldDataConstraint2015@data$S_Rate5 <- over(yieldDataConstraint2015, tank5_2015)$S_Rate;

yieldDataConstraint2015@data[which(is.na(yieldDataConstraint2015@data$Seed_Rate1)), ]$Seed_Rate1 <- 0;
yieldDataConstraint2015@data[which(is.na(yieldDataConstraint2015@data$Seed_Rate2)), ]$Seed_Rate2 <- 0;
yieldDataConstraint2015@data[which(is.na(yieldDataConstraint2015@data$Seed_Rate3)), ]$Seed_Rate3 <- 0;
yieldDataConstraint2015@data[which(is.na(yieldDataConstraint2015@data$Seed_Rate4)), ]$Seed_Rate4 <- 0;
yieldDataConstraint2015@data[which(is.na(yieldDataConstraint2015@data$Seed_Rate5)), ]$Seed_Rate5 <- 0;

yieldDataConstraint2015@data[which(is.na(yieldDataConstraint2015@data$Ino_Rate1)), ]$Ino_Rate1 <- 0;
yieldDataConstraint2015@data[which(is.na(yieldDataConstraint2015@data$Ino_Rate2)), ]$Ino_Rate2 <- 0;
yieldDataConstraint2015@data[which(is.na(yieldDataConstraint2015@data$Ino_Rate3)), ]$Ino_Rate3 <- 0;
yieldDataConstraint2015@data[which(is.na(yieldDataConstraint2015@data$Ino_Rate4)), ]$Ino_Rate4 <- 0;
yieldDataConstraint2015@data[which(is.na(yieldDataConstraint2015@data$Ino_Rate5)), ]$Ino_Rate5 <- 0;

yieldDataConstraint2015@data[which(is.na(yieldDataConstraint2015@data$N_Rate1)), ]$N_Rate1 <- 0;
yieldDataConstraint2015@data[which(is.na(yieldDataConstraint2015@data$N_Rate2)), ]$N_Rate2 <- 0;
yieldDataConstraint2015@data[which(is.na(yieldDataConstraint2015@data$N_Rate3)), ]$N_Rate3 <- 0;
yieldDataConstraint2015@data[which(is.na(yieldDataConstraint2015@data$N_Rate4)), ]$N_Rate4 <- 0;
yieldDataConstraint2015@data[which(is.na(yieldDataConstraint2015@data$N_Rate5)), ]$N_Rate5 <- 0;

yieldDataConstraint2015@data[which(is.na(yieldDataConstraint2015@data$P_Rate1)), ]$P_Rate1 <- 0;
yieldDataConstraint2015@data[which(is.na(yieldDataConstraint2015@data$P_Rate2)), ]$P_Rate2 <- 0;
yieldDataConstraint2015@data[which(is.na(yieldDataConstraint2015@data$P_Rate3)), ]$P_Rate3 <- 0;
yieldDataConstraint2015@data[which(is.na(yieldDataConstraint2015@data$P_Rate4)), ]$P_Rate4 <- 0;
yieldDataConstraint2015@data[which(is.na(yieldDataConstraint2015@data$P_Rate5)), ]$P_Rate5 <- 0;

yieldDataConstraint2015@data[which(is.na(yieldDataConstraint2015@data$K_Rate1)), ]$K_Rate1 <- 0;
yieldDataConstraint2015@data[which(is.na(yieldDataConstraint2015@data$K_Rate2)), ]$K_Rate2 <- 0;
yieldDataConstraint2015@data[which(is.na(yieldDataConstraint2015@data$K_Rate3)), ]$K_Rate3 <- 0;
yieldDataConstraint2015@data[which(is.na(yieldDataConstraint2015@data$K_Rate4)), ]$K_Rate4 <- 0;
yieldDataConstraint2015@data[which(is.na(yieldDataConstraint2015@data$K_Rate5)), ]$K_Rate5 <- 0;

yieldDataConstraint2015@data[which(is.na(yieldDataConstraint2015@data$S_Rate1)), ]$S_Rate1 <- 0;
yieldDataConstraint2015@data[which(is.na(yieldDataConstraint2015@data$S_Rate2)), ]$S_Rate2 <- 0;
yieldDataConstraint2015@data[which(is.na(yieldDataConstraint2015@data$S_Rate3)), ]$S_Rate3 <- 0;
yieldDataConstraint2015@data[which(is.na(yieldDataConstraint2015@data$S_Rate4)), ]$S_Rate4 <- 0;
yieldDataConstraint2015@data[which(is.na(yieldDataConstraint2015@data$S_Rate5)), ]$S_Rate5 <- 0;

yieldDataConstraint2015@data$Seed_Rate <- yieldDataConstraint2015@data$Seed_Rate1 + yieldDataConstraint2015@data$Seed_Rate2 + yieldDataConstraint2015@data$Seed_Rate3 + yieldDataConstraint2015@data$Seed_Rate4 + yieldDataConstraint2015@data$Seed_Rate5;
yieldDataConstraint2015@data$Ino_Rate <- yieldDataConstraint2015@data$Ino_Rate1 + yieldDataConstraint2015@data$Ino_Rate2 + yieldDataConstraint2015@data$Ino_Rate3 + yieldDataConstraint2015@data$Ino_Rate4 + yieldDataConstraint2015@data$Ino_Rate5;
yieldDataConstraint2015@data$N_Rate <- yieldDataConstraint2015@data$N_Rate1 + yieldDataConstraint2015@data$N_Rate2 + yieldDataConstraint2015@data$N_Rate3 + yieldDataConstraint2015@data$N_Rate4 + yieldDataConstraint2015@data$N_Rate5;
yieldDataConstraint2015@data$P_Rate <- yieldDataConstraint2015@data$P_Rate1 + yieldDataConstraint2015@data$P_Rate2 + yieldDataConstraint2015@data$P_Rate3 + yieldDataConstraint2015@data$P_Rate4 + yieldDataConstraint2015@data$P_Rate5;
yieldDataConstraint2015@data$K_Rate <- yieldDataConstraint2015@data$K_Rate1 + yieldDataConstraint2015@data$K_Rate2 + yieldDataConstraint2015@data$K_Rate3 + yieldDataConstraint2015@data$K_Rate4 + yieldDataConstraint2015@data$K_Rate5;
yieldDataConstraint2015@data$S_Rate <- yieldDataConstraint2015@data$S_Rate1 + yieldDataConstraint2015@data$S_Rate2 + yieldDataConstraint2015@data$S_Rate3 + yieldDataConstraint2015@data$S_Rate4 + yieldDataConstraint2015@data$S_Rate5;

yieldDataConstraint2015 <- yieldDataConstraint2015[, c("Elevation", "DryYield", "Product", "Machine", "ProcYear", "Harv_Date", "Field", "Seed_Date", "Seed_Rate", "Ino_Rate", "N_Rate", "P_Rate", "K_Rate", "S_Rate")];
writePointsShape(yieldDataConstraint2015, paste("yieldDataConstraint_allData_2015.shp", sep=""));
###yieldDataFrameConstraint2015 <- yieldDataConstraint2015@data
###save(yieldDataFrameConstraint2015, file=paste("yieldDataFrameConstraint_allData_2015", ".RData", sep=""), compress="xz");


#---- 2016;
##BELOW does point in polygon calculation, tank 3 wasn't used this year;
yieldDataConstraint2016@data$Seed_Date <- over(yieldDataConstraint2016, tank1_2016)$SeedDate;

yieldDataConstraint2016@data$Seed_Rate1 <- over(yieldDataConstraint2016, tank1_2016)$Seed_Rate;
yieldDataConstraint2016@data$Seed_Rate2 <- over(yieldDataConstraint2016, tank2_2016)$Seed_Rate;
yieldDataConstraint2016@data$Seed_Rate4 <- over(yieldDataConstraint2016, tank4_2016)$Seed_Rate;
yieldDataConstraint2016@data$Seed_Rate5 <- over(yieldDataConstraint2016, tank5_2016)$Seed_Rate;

yieldDataConstraint2016@data$Ino_Rate1 <- over(yieldDataConstraint2016, tank1_2016)$Ino_Rate;
yieldDataConstraint2016@data$Ino_Rate2 <- over(yieldDataConstraint2016, tank2_2016)$Ino_Rate;
yieldDataConstraint2016@data$Ino_Rate4 <- over(yieldDataConstraint2016, tank4_2016)$Ino_Rate;
yieldDataConstraint2016@data$Ino_Rate5 <- over(yieldDataConstraint2016, tank5_2016)$Ino_Rate;

yieldDataConstraint2016@data$N_Rate1 <- over(yieldDataConstraint2016, tank1_2016)$N_Rate;
yieldDataConstraint2016@data$N_Rate2 <- over(yieldDataConstraint2016, tank2_2016)$N_Rate;
yieldDataConstraint2016@data$N_Rate4 <- over(yieldDataConstraint2016, tank4_2016)$N_Rate;
yieldDataConstraint2016@data$N_Rate5 <- over(yieldDataConstraint2016, tank5_2016)$N_Rate;

yieldDataConstraint2016@data$P_Rate1 <- over(yieldDataConstraint2016, tank1_2016)$P_Rate;
yieldDataConstraint2016@data$P_Rate2 <- over(yieldDataConstraint2016, tank2_2016)$P_Rate;
yieldDataConstraint2016@data$P_Rate4 <- over(yieldDataConstraint2016, tank4_2016)$P_Rate;
yieldDataConstraint2016@data$P_Rate5 <- over(yieldDataConstraint2016, tank5_2016)$P_Rate;

yieldDataConstraint2016@data$K_Rate1 <- over(yieldDataConstraint2016, tank1_2016)$K_Rate;
yieldDataConstraint2016@data$K_Rate2 <- over(yieldDataConstraint2016, tank2_2016)$K_Rate;
yieldDataConstraint2016@data$K_Rate4 <- over(yieldDataConstraint2016, tank4_2016)$K_Rate;
yieldDataConstraint2016@data$K_Rate5 <- over(yieldDataConstraint2016, tank5_2016)$K_Rate;

yieldDataConstraint2016@data$S_Rate1 <- over(yieldDataConstraint2016, tank1_2016)$S_Rate;
yieldDataConstraint2016@data$S_Rate2 <- over(yieldDataConstraint2016, tank2_2016)$S_Rate;
yieldDataConstraint2016@data$S_Rate4 <- over(yieldDataConstraint2016, tank4_2016)$S_Rate;
yieldDataConstraint2016@data$S_Rate5 <- over(yieldDataConstraint2016, tank5_2016)$S_Rate;

yieldDataConstraint2016@data[which(is.na(yieldDataConstraint2016@data$Seed_Rate1)), ]$Seed_Rate1 <- 0;
yieldDataConstraint2016@data[which(is.na(yieldDataConstraint2016@data$Seed_Rate2)), ]$Seed_Rate2 <- 0;
yieldDataConstraint2016@data[which(is.na(yieldDataConstraint2016@data$Seed_Rate4)), ]$Seed_Rate4 <- 0;
yieldDataConstraint2016@data[which(is.na(yieldDataConstraint2016@data$Seed_Rate5)), ]$Seed_Rate5 <- 0;

yieldDataConstraint2016@data[which(is.na(yieldDataConstraint2016@data$Ino_Rate1)), ]$Ino_Rate1 <- 0;
yieldDataConstraint2016@data[which(is.na(yieldDataConstraint2016@data$Ino_Rate2)), ]$Ino_Rate2 <- 0;
yieldDataConstraint2016@data[which(is.na(yieldDataConstraint2016@data$Ino_Rate4)), ]$Ino_Rate4 <- 0;
yieldDataConstraint2016@data[which(is.na(yieldDataConstraint2016@data$Ino_Rate5)), ]$Ino_Rate5 <- 0;

yieldDataConstraint2016@data[which(is.na(yieldDataConstraint2016@data$N_Rate1)), ]$N_Rate1 <- 0;
yieldDataConstraint2016@data[which(is.na(yieldDataConstraint2016@data$N_Rate2)), ]$N_Rate2 <- 0;
yieldDataConstraint2016@data[which(is.na(yieldDataConstraint2016@data$N_Rate4)), ]$N_Rate4 <- 0;
yieldDataConstraint2016@data[which(is.na(yieldDataConstraint2016@data$N_Rate5)), ]$N_Rate5 <- 0;

yieldDataConstraint2016@data[which(is.na(yieldDataConstraint2016@data$P_Rate1)), ]$P_Rate1 <- 0;
yieldDataConstraint2016@data[which(is.na(yieldDataConstraint2016@data$P_Rate2)), ]$P_Rate2 <- 0;
yieldDataConstraint2016@data[which(is.na(yieldDataConstraint2016@data$P_Rate4)), ]$P_Rate4 <- 0;
yieldDataConstraint2016@data[which(is.na(yieldDataConstraint2016@data$P_Rate5)), ]$P_Rate5 <- 0;

yieldDataConstraint2016@data[which(is.na(yieldDataConstraint2016@data$K_Rate1)), ]$K_Rate1 <- 0;
yieldDataConstraint2016@data[which(is.na(yieldDataConstraint2016@data$K_Rate2)), ]$K_Rate2 <- 0;
yieldDataConstraint2016@data[which(is.na(yieldDataConstraint2016@data$K_Rate4)), ]$K_Rate4 <- 0;
yieldDataConstraint2016@data[which(is.na(yieldDataConstraint2016@data$K_Rate5)), ]$K_Rate5 <- 0;

yieldDataConstraint2016@data[which(is.na(yieldDataConstraint2016@data$S_Rate1)), ]$S_Rate1 <- 0;
yieldDataConstraint2016@data[which(is.na(yieldDataConstraint2016@data$S_Rate2)), ]$S_Rate2 <- 0;
yieldDataConstraint2016@data[which(is.na(yieldDataConstraint2016@data$S_Rate4)), ]$S_Rate4 <- 0;
yieldDataConstraint2016@data[which(is.na(yieldDataConstraint2016@data$S_Rate5)), ]$S_Rate5 <- 0;

yieldDataConstraint2016@data$Seed_Rate <- yieldDataConstraint2016@data$Seed_Rate1 + yieldDataConstraint2016@data$Seed_Rate2 + yieldDataConstraint2016@data$Seed_Rate4 + yieldDataConstraint2016@data$Seed_Rate5;
yieldDataConstraint2016@data$Ino_Rate <- yieldDataConstraint2016@data$Ino_Rate1 + yieldDataConstraint2016@data$Ino_Rate2 + yieldDataConstraint2016@data$Ino_Rate4 + yieldDataConstraint2016@data$Ino_Rate5;
yieldDataConstraint2016@data$N_Rate <- yieldDataConstraint2016@data$N_Rate1 + yieldDataConstraint2016@data$N_Rate2 + yieldDataConstraint2016@data$N_Rate4 + yieldDataConstraint2016@data$N_Rate5;
yieldDataConstraint2016@data$P_Rate <- yieldDataConstraint2016@data$P_Rate1 + yieldDataConstraint2016@data$P_Rate2 + yieldDataConstraint2016@data$P_Rate4 + yieldDataConstraint2016@data$P_Rate5;
yieldDataConstraint2016@data$K_Rate <- yieldDataConstraint2016@data$K_Rate1 + yieldDataConstraint2016@data$K_Rate2 + yieldDataConstraint2016@data$K_Rate4 + yieldDataConstraint2016@data$K_Rate5;
yieldDataConstraint2016@data$S_Rate <- yieldDataConstraint2016@data$S_Rate1 + yieldDataConstraint2016@data$S_Rate2 + yieldDataConstraint2016@data$S_Rate4 + yieldDataConstraint2016@data$S_Rate5;

yieldDataConstraint2016 <- yieldDataConstraint2016[, c("Elevation", "DryYield", "Product", "Machine", "ProcYear", "Harv_Date", "Field", "Seed_Date", "Seed_Rate", "Ino_Rate", "N_Rate", "P_Rate", "K_Rate", "S_Rate")];
writePointsShape(yieldDataConstraint2016, paste("yieldDataConstraint_allData_2016.shp", sep=""));
##yieldDataFrameConstraint2016 <- yieldDataConstraint2016@data
##save(yieldDataFrameConstraint2016, file=paste("yieldDataFrameConstraint_allData_2016", ".RData", sep=""), compress="xz");


#---- 2017;
##BELOW does point in polygon calculation;
yieldDataConstraint2017@data$Seed_Date <- over(yieldDataConstraint2017, tank1_2017)$SeedDate;

yieldDataConstraint2017@data$Seed_Rate1 <- over(yieldDataConstraint2017, tank1_2017)$Seed_Rate;
yieldDataConstraint2017@data$Seed_Rate2 <- over(yieldDataConstraint2017, tank2_2017)$Seed_Rate;
yieldDataConstraint2017@data$Seed_Rate3 <- over(yieldDataConstraint2017, tank3_2017)$Seed_Rate;
yieldDataConstraint2017@data$Seed_Rate4 <- over(yieldDataConstraint2017, tank4_2017)$Seed_Rate;
yieldDataConstraint2017@data$Seed_Rate5 <- over(yieldDataConstraint2017, tank5_2017)$Seed_Rate;

yieldDataConstraint2017@data$Ino_Rate1 <- over(yieldDataConstraint2017, tank1_2017)$Ino_Rate;
yieldDataConstraint2017@data$Ino_Rate2 <- over(yieldDataConstraint2017, tank2_2017)$Ino_Rate;
yieldDataConstraint2017@data$Ino_Rate3 <- over(yieldDataConstraint2017, tank3_2017)$Ino_Rate;
yieldDataConstraint2017@data$Ino_Rate4 <- over(yieldDataConstraint2017, tank4_2017)$Ino_Rate;
yieldDataConstraint2017@data$Ino_Rate5 <- over(yieldDataConstraint2017, tank5_2017)$Ino_Rate;

yieldDataConstraint2017@data$N_Rate1 <- over(yieldDataConstraint2017, tank1_2017)$N_Rate;
yieldDataConstraint2017@data$N_Rate2 <- over(yieldDataConstraint2017, tank2_2017)$N_Rate;
yieldDataConstraint2017@data$N_Rate3 <- over(yieldDataConstraint2017, tank3_2017)$N_Rate;
yieldDataConstraint2017@data$N_Rate4 <- over(yieldDataConstraint2017, tank4_2017)$N_Rate;
yieldDataConstraint2017@data$N_Rate5 <- over(yieldDataConstraint2017, tank5_2017)$N_Rate;

yieldDataConstraint2017@data$P_Rate1 <- over(yieldDataConstraint2017, tank1_2017)$P_Rate;
yieldDataConstraint2017@data$P_Rate2 <- over(yieldDataConstraint2017, tank2_2017)$P_Rate;
yieldDataConstraint2017@data$P_Rate3 <- over(yieldDataConstraint2017, tank3_2017)$P_Rate;
yieldDataConstraint2017@data$P_Rate4 <- over(yieldDataConstraint2017, tank4_2017)$P_Rate;
yieldDataConstraint2017@data$P_Rate5 <- over(yieldDataConstraint2017, tank5_2017)$P_Rate;

yieldDataConstraint2017@data$K_Rate1 <- over(yieldDataConstraint2017, tank1_2017)$K_Rate;
yieldDataConstraint2017@data$K_Rate2 <- over(yieldDataConstraint2017, tank2_2017)$K_Rate;
yieldDataConstraint2017@data$K_Rate3 <- over(yieldDataConstraint2017, tank3_2017)$K_Rate;
yieldDataConstraint2017@data$K_Rate4 <- over(yieldDataConstraint2017, tank4_2017)$K_Rate;
yieldDataConstraint2017@data$K_Rate5 <- over(yieldDataConstraint2017, tank5_2017)$K_Rate;

yieldDataConstraint2017@data$S_Rate1 <- over(yieldDataConstraint2017, tank1_2017)$S_Rate;
yieldDataConstraint2017@data$S_Rate2 <- over(yieldDataConstraint2017, tank2_2017)$S_Rate;
yieldDataConstraint2017@data$S_Rate3 <- over(yieldDataConstraint2017, tank3_2017)$S_Rate;
yieldDataConstraint2017@data$S_Rate4 <- over(yieldDataConstraint2017, tank4_2017)$S_Rate;
yieldDataConstraint2017@data$S_Rate5 <- over(yieldDataConstraint2017, tank5_2017)$S_Rate;

yieldDataConstraint2017@data[which(is.na(yieldDataConstraint2017@data$Seed_Rate1)), ]$Seed_Rate1 <- 0;
yieldDataConstraint2017@data[which(is.na(yieldDataConstraint2017@data$Seed_Rate2)), ]$Seed_Rate2 <- 0;
yieldDataConstraint2017@data[which(is.na(yieldDataConstraint2017@data$Seed_Rate3)), ]$Seed_Rate3 <- 0;
yieldDataConstraint2017@data[which(is.na(yieldDataConstraint2017@data$Seed_Rate4)), ]$Seed_Rate4 <- 0;
yieldDataConstraint2017@data[which(is.na(yieldDataConstraint2017@data$Seed_Rate5)), ]$Seed_Rate5 <- 0;

yieldDataConstraint2017@data[which(is.na(yieldDataConstraint2017@data$Ino_Rate1)), ]$Ino_Rate1 <- 0;
yieldDataConstraint2017@data[which(is.na(yieldDataConstraint2017@data$Ino_Rate2)), ]$Ino_Rate2 <- 0;
yieldDataConstraint2017@data[which(is.na(yieldDataConstraint2017@data$Ino_Rate3)), ]$Ino_Rate3 <- 0;
yieldDataConstraint2017@data[which(is.na(yieldDataConstraint2017@data$Ino_Rate4)), ]$Ino_Rate4 <- 0;
yieldDataConstraint2017@data[which(is.na(yieldDataConstraint2017@data$Ino_Rate5)), ]$Ino_Rate5 <- 0;

yieldDataConstraint2017@data[which(is.na(yieldDataConstraint2017@data$N_Rate1)), ]$N_Rate1 <- 0;
yieldDataConstraint2017@data[which(is.na(yieldDataConstraint2017@data$N_Rate2)), ]$N_Rate2 <- 0;
yieldDataConstraint2017@data[which(is.na(yieldDataConstraint2017@data$N_Rate3)), ]$N_Rate3 <- 0;
yieldDataConstraint2017@data[which(is.na(yieldDataConstraint2017@data$N_Rate4)), ]$N_Rate4 <- 0;
yieldDataConstraint2017@data[which(is.na(yieldDataConstraint2017@data$N_Rate5)), ]$N_Rate5 <- 0;

yieldDataConstraint2017@data[which(is.na(yieldDataConstraint2017@data$P_Rate1)), ]$P_Rate1 <- 0;
yieldDataConstraint2017@data[which(is.na(yieldDataConstraint2017@data$P_Rate2)), ]$P_Rate2 <- 0;
yieldDataConstraint2017@data[which(is.na(yieldDataConstraint2017@data$P_Rate3)), ]$P_Rate3 <- 0;
yieldDataConstraint2017@data[which(is.na(yieldDataConstraint2017@data$P_Rate4)), ]$P_Rate4 <- 0;
yieldDataConstraint2017@data[which(is.na(yieldDataConstraint2017@data$P_Rate5)), ]$P_Rate5 <- 0;

yieldDataConstraint2017@data[which(is.na(yieldDataConstraint2017@data$K_Rate1)), ]$K_Rate1 <- 0;
yieldDataConstraint2017@data[which(is.na(yieldDataConstraint2017@data$K_Rate2)), ]$K_Rate2 <- 0;
yieldDataConstraint2017@data[which(is.na(yieldDataConstraint2017@data$K_Rate3)), ]$K_Rate3 <- 0;
yieldDataConstraint2017@data[which(is.na(yieldDataConstraint2017@data$K_Rate4)), ]$K_Rate4 <- 0;
yieldDataConstraint2017@data[which(is.na(yieldDataConstraint2017@data$K_Rate5)), ]$K_Rate5 <- 0;

yieldDataConstraint2017@data[which(is.na(yieldDataConstraint2017@data$S_Rate1)), ]$S_Rate1 <- 0;
yieldDataConstraint2017@data[which(is.na(yieldDataConstraint2017@data$S_Rate2)), ]$S_Rate2 <- 0;
yieldDataConstraint2017@data[which(is.na(yieldDataConstraint2017@data$S_Rate3)), ]$S_Rate3 <- 0;
yieldDataConstraint2017@data[which(is.na(yieldDataConstraint2017@data$S_Rate4)), ]$S_Rate4 <- 0;
yieldDataConstraint2017@data[which(is.na(yieldDataConstraint2017@data$S_Rate5)), ]$S_Rate5 <- 0;

yieldDataConstraint2017@data$Seed_Rate <- yieldDataConstraint2017@data$Seed_Rate1 + yieldDataConstraint2017@data$Seed_Rate2 + yieldDataConstraint2017@data$Seed_Rate3 + yieldDataConstraint2017@data$Seed_Rate4 + yieldDataConstraint2017@data$Seed_Rate5;
yieldDataConstraint2017@data$Ino_Rate <- yieldDataConstraint2017@data$Ino_Rate1 + yieldDataConstraint2017@data$Ino_Rate2 + yieldDataConstraint2017@data$Ino_Rate3 + yieldDataConstraint2017@data$Ino_Rate4 + yieldDataConstraint2017@data$Ino_Rate5;
yieldDataConstraint2017@data$N_Rate <- yieldDataConstraint2017@data$N_Rate1 + yieldDataConstraint2017@data$N_Rate2 + yieldDataConstraint2017@data$N_Rate3 + yieldDataConstraint2017@data$N_Rate4 + yieldDataConstraint2017@data$N_Rate5;
yieldDataConstraint2017@data$P_Rate <- yieldDataConstraint2017@data$P_Rate1 + yieldDataConstraint2017@data$P_Rate2 + yieldDataConstraint2017@data$P_Rate3 + yieldDataConstraint2017@data$P_Rate4 + yieldDataConstraint2017@data$P_Rate5;
yieldDataConstraint2017@data$K_Rate <- yieldDataConstraint2017@data$K_Rate1 + yieldDataConstraint2017@data$K_Rate2 + yieldDataConstraint2017@data$K_Rate3 + yieldDataConstraint2017@data$K_Rate4 + yieldDataConstraint2017@data$K_Rate5;
yieldDataConstraint2017@data$S_Rate <- yieldDataConstraint2017@data$S_Rate1 + yieldDataConstraint2017@data$S_Rate2 + yieldDataConstraint2017@data$S_Rate3 + yieldDataConstraint2017@data$S_Rate4 + yieldDataConstraint2017@data$S_Rate5;

yieldDataConstraint2017 <- yieldDataConstraint2017[, c("Elevation", "DryYield", "Product", "Machine", "ProcYear", "Harv_Date", "Field", "Seed_Date", "Seed_Rate", "Ino_Rate", "N_Rate", "P_Rate", "K_Rate", "S_Rate")];
writePointsShape(yieldDataConstraint2017, paste("yieldDataConstraint_allData_2017.shp", sep=""));
##yieldDataFrameConstraint2017 <- yieldDataConstraint2017@data
##save(yieldDataFrameConstraint2017, file=paste("yieldDataFrameConstraint_allData_2017", ".RData", sep=""), compress="xz");


#---- 2018;
##BELOW does point in polygon calculation;
yieldDataConstraint2018@data$Seed_Date <- over(yieldDataConstraint2018, tank1_2018)$SeedDate;

yieldDataConstraint2018@data$Seed_Rate1 <- over(yieldDataConstraint2018, tank1_2018)$Seed_Rate;
yieldDataConstraint2018@data$Seed_Rate2 <- over(yieldDataConstraint2018, tank2_2018)$Seed_Rate;
yieldDataConstraint2018@data$Seed_Rate3 <- over(yieldDataConstraint2018, tank3_2018)$Seed_Rate;
yieldDataConstraint2018@data$Seed_Rate4 <- over(yieldDataConstraint2018, tank4_2018)$Seed_Rate;
yieldDataConstraint2018@data$Seed_Rate5 <- over(yieldDataConstraint2018, tank5_2018)$Seed_Rate;

yieldDataConstraint2018@data$Ino_Rate1 <- over(yieldDataConstraint2018, tank1_2018)$Ino_Rate;
yieldDataConstraint2018@data$Ino_Rate2 <- over(yieldDataConstraint2018, tank2_2018)$Ino_Rate;
yieldDataConstraint2018@data$Ino_Rate3 <- over(yieldDataConstraint2018, tank3_2018)$Ino_Rate;
yieldDataConstraint2018@data$Ino_Rate4 <- over(yieldDataConstraint2018, tank4_2018)$Ino_Rate;
yieldDataConstraint2018@data$Ino_Rate5 <- over(yieldDataConstraint2018, tank5_2018)$Ino_Rate;

yieldDataConstraint2018@data$N_Rate1 <- over(yieldDataConstraint2018, tank1_2018)$N_Rate;
yieldDataConstraint2018@data$N_Rate2 <- over(yieldDataConstraint2018, tank2_2018)$N_Rate;
yieldDataConstraint2018@data$N_Rate3 <- over(yieldDataConstraint2018, tank3_2018)$N_Rate;
yieldDataConstraint2018@data$N_Rate4 <- over(yieldDataConstraint2018, tank4_2018)$N_Rate;
yieldDataConstraint2018@data$N_Rate5 <- over(yieldDataConstraint2018, tank5_2018)$N_Rate;

yieldDataConstraint2018@data$P_Rate1 <- over(yieldDataConstraint2018, tank1_2018)$P_Rate;
yieldDataConstraint2018@data$P_Rate2 <- over(yieldDataConstraint2018, tank2_2018)$P_Rate;
yieldDataConstraint2018@data$P_Rate3 <- over(yieldDataConstraint2018, tank3_2018)$P_Rate;
yieldDataConstraint2018@data$P_Rate4 <- over(yieldDataConstraint2018, tank4_2018)$P_Rate;
yieldDataConstraint2018@data$P_Rate5 <- over(yieldDataConstraint2018, tank5_2018)$P_Rate;

yieldDataConstraint2018@data$K_Rate1 <- over(yieldDataConstraint2018, tank1_2018)$K_Rate;
yieldDataConstraint2018@data$K_Rate2 <- over(yieldDataConstraint2018, tank2_2018)$K_Rate;
yieldDataConstraint2018@data$K_Rate3 <- over(yieldDataConstraint2018, tank3_2018)$K_Rate;
yieldDataConstraint2018@data$K_Rate4 <- over(yieldDataConstraint2018, tank4_2018)$K_Rate;
yieldDataConstraint2018@data$K_Rate5 <- over(yieldDataConstraint2018, tank5_2018)$K_Rate;

yieldDataConstraint2018@data$S_Rate1 <- over(yieldDataConstraint2018, tank1_2018)$S_Rate;
yieldDataConstraint2018@data$S_Rate2 <- over(yieldDataConstraint2018, tank2_2018)$S_Rate;
yieldDataConstraint2018@data$S_Rate3 <- over(yieldDataConstraint2018, tank3_2018)$S_Rate;
yieldDataConstraint2018@data$S_Rate4 <- over(yieldDataConstraint2018, tank4_2018)$S_Rate;
yieldDataConstraint2018@data$S_Rate5 <- over(yieldDataConstraint2018, tank5_2018)$S_Rate;

yieldDataConstraint2018@data[which(is.na(yieldDataConstraint2018@data$Seed_Rate1)), ]$Seed_Rate1 <- 0;
yieldDataConstraint2018@data[which(is.na(yieldDataConstraint2018@data$Seed_Rate2)), ]$Seed_Rate2 <- 0;
yieldDataConstraint2018@data[which(is.na(yieldDataConstraint2018@data$Seed_Rate3)), ]$Seed_Rate3 <- 0;
yieldDataConstraint2018@data[which(is.na(yieldDataConstraint2018@data$Seed_Rate4)), ]$Seed_Rate4 <- 0;
yieldDataConstraint2018@data[which(is.na(yieldDataConstraint2018@data$Seed_Rate5)), ]$Seed_Rate5 <- 0;

yieldDataConstraint2018@data[which(is.na(yieldDataConstraint2018@data$Ino_Rate1)), ]$Ino_Rate1 <- 0;
yieldDataConstraint2018@data[which(is.na(yieldDataConstraint2018@data$Ino_Rate2)), ]$Ino_Rate2 <- 0;
yieldDataConstraint2018@data[which(is.na(yieldDataConstraint2018@data$Ino_Rate3)), ]$Ino_Rate3 <- 0;
yieldDataConstraint2018@data[which(is.na(yieldDataConstraint2018@data$Ino_Rate4)), ]$Ino_Rate4 <- 0;
yieldDataConstraint2018@data[which(is.na(yieldDataConstraint2018@data$Ino_Rate5)), ]$Ino_Rate5 <- 0;

yieldDataConstraint2018@data[which(is.na(yieldDataConstraint2018@data$N_Rate1)), ]$N_Rate1 <- 0;
yieldDataConstraint2018@data[which(is.na(yieldDataConstraint2018@data$N_Rate2)), ]$N_Rate2 <- 0;
yieldDataConstraint2018@data[which(is.na(yieldDataConstraint2018@data$N_Rate3)), ]$N_Rate3 <- 0;
yieldDataConstraint2018@data[which(is.na(yieldDataConstraint2018@data$N_Rate4)), ]$N_Rate4 <- 0;
yieldDataConstraint2018@data[which(is.na(yieldDataConstraint2018@data$N_Rate5)), ]$N_Rate5 <- 0;

yieldDataConstraint2018@data[which(is.na(yieldDataConstraint2018@data$P_Rate1)), ]$P_Rate1 <- 0;
yieldDataConstraint2018@data[which(is.na(yieldDataConstraint2018@data$P_Rate2)), ]$P_Rate2 <- 0;
yieldDataConstraint2018@data[which(is.na(yieldDataConstraint2018@data$P_Rate3)), ]$P_Rate3 <- 0;
yieldDataConstraint2018@data[which(is.na(yieldDataConstraint2018@data$P_Rate4)), ]$P_Rate4 <- 0;
yieldDataConstraint2018@data[which(is.na(yieldDataConstraint2018@data$P_Rate5)), ]$P_Rate5 <- 0;

yieldDataConstraint2018@data[which(is.na(yieldDataConstraint2018@data$K_Rate1)), ]$K_Rate1 <- 0;
yieldDataConstraint2018@data[which(is.na(yieldDataConstraint2018@data$K_Rate2)), ]$K_Rate2 <- 0;
yieldDataConstraint2018@data[which(is.na(yieldDataConstraint2018@data$K_Rate3)), ]$K_Rate3 <- 0;
yieldDataConstraint2018@data[which(is.na(yieldDataConstraint2018@data$K_Rate4)), ]$K_Rate4 <- 0;
yieldDataConstraint2018@data[which(is.na(yieldDataConstraint2018@data$K_Rate5)), ]$K_Rate5 <- 0;

yieldDataConstraint2018@data[which(is.na(yieldDataConstraint2018@data$S_Rate1)), ]$S_Rate1 <- 0;
yieldDataConstraint2018@data[which(is.na(yieldDataConstraint2018@data$S_Rate2)), ]$S_Rate2 <- 0;
yieldDataConstraint2018@data[which(is.na(yieldDataConstraint2018@data$S_Rate3)), ]$S_Rate3 <- 0;
yieldDataConstraint2018@data[which(is.na(yieldDataConstraint2018@data$S_Rate4)), ]$S_Rate4 <- 0;
yieldDataConstraint2018@data[which(is.na(yieldDataConstraint2018@data$S_Rate5)), ]$S_Rate5 <- 0;

yieldDataConstraint2018@data$Seed_Rate <- yieldDataConstraint2018@data$Seed_Rate1 + yieldDataConstraint2018@data$Seed_Rate2 + yieldDataConstraint2018@data$Seed_Rate3 + yieldDataConstraint2018@data$Seed_Rate4 + yieldDataConstraint2018@data$Seed_Rate5;
yieldDataConstraint2018@data$Ino_Rate <- yieldDataConstraint2018@data$Ino_Rate1 + yieldDataConstraint2018@data$Ino_Rate2 + yieldDataConstraint2018@data$Ino_Rate3 + yieldDataConstraint2018@data$Ino_Rate4 + yieldDataConstraint2018@data$Ino_Rate5;
yieldDataConstraint2018@data$N_Rate <- yieldDataConstraint2018@data$N_Rate1 + yieldDataConstraint2018@data$N_Rate2 + yieldDataConstraint2018@data$N_Rate3 + yieldDataConstraint2018@data$N_Rate4 + yieldDataConstraint2018@data$N_Rate5;
yieldDataConstraint2018@data$P_Rate <- yieldDataConstraint2018@data$P_Rate1 + yieldDataConstraint2018@data$P_Rate2 + yieldDataConstraint2018@data$P_Rate3 + yieldDataConstraint2018@data$P_Rate4 + yieldDataConstraint2018@data$P_Rate5;
yieldDataConstraint2018@data$K_Rate <- yieldDataConstraint2018@data$K_Rate1 + yieldDataConstraint2018@data$K_Rate2 + yieldDataConstraint2018@data$K_Rate3 + yieldDataConstraint2018@data$K_Rate4 + yieldDataConstraint2018@data$K_Rate5;
yieldDataConstraint2018@data$S_Rate <- yieldDataConstraint2018@data$S_Rate1 + yieldDataConstraint2018@data$S_Rate2 + yieldDataConstraint2018@data$S_Rate3 + yieldDataConstraint2018@data$S_Rate4 + yieldDataConstraint2018@data$S_Rate5;

yieldDataConstraint2018 <- yieldDataConstraint2018[, c("Elevation", "DryYield", "Product", "Machine", "ProcYear", "Harv_Date", "Field", "Seed_Date", "Seed_Rate", "Ino_Rate", "N_Rate", "P_Rate", "K_Rate", "S_Rate")];
writePointsShape(yieldDataConstraint2018, paste("yieldDataConstraint_allData_2018.shp", sep=""));
##yieldDataFrameConstraint2018 <- yieldDataConstraint2018@data
##save(yieldDataFrameConstraint2018, file=paste("yieldDataFrameConstraint_allData_2018", ".RData", sep=""), compress="xz");
### ------------------------------- ;

#---------------------------------------------------------------------------------------#;
##### ------ ------ ------ ##### ------ ------ ------ #####; 

yieldDataConstraint2015 <- readOGR(dsn=".", layer=paste("yieldDataConstraint_allData_2015", sep=""));
yieldDataConstraint2016 <- readOGR(dsn=".", layer=paste("yieldDataConstraint_allData_2016", sep=""));
yieldDataConstraint2017 <- readOGR(dsn=".", layer=paste("yieldDataConstraint_allData_2017", sep=""));
yieldDataConstraint2018 <- readOGR(dsn=".", layer=paste("yieldDataConstraint_allData_2018", sep=""));

#str(yieldDataConstraint2015@data);
#str(yieldDataConstraint2016@data);
#str(yieldDataConstraint2017@data);
#str(yieldDataConstraint2018@data);
#summary(yieldDataConstraint2015@data$Seed_Rate);
#summary(yieldDataConstraint2016@data$Seed_Rate);
#summary(yieldDataConstraint2017@data$Seed_Rate);
#summary(yieldDataConstraint2018@data$Seed_Rate);

farm2015 <- yieldDataConstraint2015@data[-which(yieldDataConstraint2015@data$Seed_Rate == 0), ];
farm2016 <- yieldDataConstraint2016@data[-which(yieldDataConstraint2016@data$Seed_Rate == 0), ];
farm2017 <- yieldDataConstraint2017@data[-which(yieldDataConstraint2017@data$Seed_Rate == 0), ];
farm2018 <- yieldDataConstraint2018@data[-which(yieldDataConstraint2018@data$Seed_Rate == 0), ];

#str(farm2015);
#str(farm2016);
#str(farm2017);
#str(farm2018);

summary(farm2015$Seed_Rate);
summary(farm2016$Seed_Rate);
summary(farm2017$Seed_Rate);
summary(farm2018$Seed_Rate);

#farm2018$longitude;
#farm2018$latitude;

farm2015$Product <- as.character(farm2015$Product);
farm2015$Machine <- as.character(farm2015$Machine);
farm2015$ProcYear <- as.character(farm2015$ProcYear);
farm2015$Harv_Date <- as.character(farm2015$Harv_Date);
farm2015$Field <- as.character(farm2015$Field);
farm2015$Seed_Date <- as.character(farm2015$Seed_Date);

farm2016$Product <- as.character(farm2016$Product);
farm2016$Machine <- as.character(farm2016$Machine);
farm2016$ProcYear <- as.character(farm2016$ProcYear);
farm2016$Harv_Date <- as.character(farm2016$Harv_Date);
farm2016$Field <- as.character(farm2016$Field);
farm2016$Seed_Date <- as.character(farm2016$Seed_Date);

farm2017$Product <- as.character(farm2017$Product);
farm2017$Machine <- as.character(farm2017$Machine);
farm2017$ProcYear <- as.character(farm2017$ProcYear);
farm2017$Harv_Date <- as.character(farm2017$Harv_Date);
farm2017$Field <- as.character(farm2017$Field);
farm2017$Seed_Date <- as.character(farm2017$Seed_Date);

farm2018$Product <- as.character(farm2018$Product);
farm2018$Machine <- as.character(farm2018$Machine);
farm2018$ProcYear <- as.character(farm2018$ProcYear);
farm2018$Harv_Date <- as.character(farm2018$Harv_Date);
farm2018$Field <- as.character(farm2018$Field);
farm2018$Seed_Date <- as.character(farm2018$Seed_Date);

farmData2015to2018 <- rbind(farm2015, farm2016, farm2017, farm2018);
#str(farmData2015to2018);
farmData2015to2018$Product <- as.factor(farmData2015to2018$Product);
farmData2015to2018$Machine <- as.factor(farmData2015to2018$Machine);
farmData2015to2018$ProcYear <- as.factor(farmData2015to2018$ProcYear);
farmData2015to2018$Harv_Date <- as.factor(farmData2015to2018$Harv_Date);
farmData2015to2018$Field <- as.factor(farmData2015to2018$Field);
farmData2015to2018$Seed_Date <- as.factor(farmData2015to2018$Seed_Date);
str(farmData2015to2018);

setwd("/Users/tylerpittman/Farm/agYieldProject/data");
save(farmData2015to2018, file=paste("farmData2015to2018", ".RData", sep=""), compress="xz");

##### ------ ------ ------ ##### ------ ------ ------ #####; 
#---------------------------------------------------------------------------------------#;
##### ------ ------ ------ ##### ------ ------ ------ #####; 

setwd("/Users/tylerpittman/Farm/agYieldProject/data");
load(paste("farmData2015to2018", ".RData", sep=""));
str(farmData2015to2018);
levels(farmData2015to2018$Field);

setwd("/Users/tylerpittman/Farm");
cultivars <- read_excel("cultivarKeyYear2015to2018.xlsx", sheet = "lt1");
setwd("/Users/tylerpittman/Farm/agYieldProject/data");


###This calculates number of days between seeding and harvest;
farmData2015to2018$Harv_Date <-  as.Date(paste(as.character(farmData2015to2018$Harv_Date), as.character(farmData2015to2018$ProcYear), sep=" "), format="%B %d %Y");
farmData2015to2018$Seed_Date <-  as.Date(paste(as.character(farmData2015to2018$Seed_Date), as.character(farmData2015to2018$ProcYear), sep=" "), format="%B %d %Y");
farmData2015to2018$Mat_Days <- farmData2015to2018$Harv_Date - farmData2015to2018$Seed_Date;
farmData2015to2018$Harv_Date <- format(as.Date(farmData2015to2018$Harv_Date), format="%B %d");
farmData2015to2018$Seed_Date <- format(as.Date(farmData2015to2018$Seed_Date), format="%B %d");
#str(farmData2015to2018);


###This merges cultivar for each product and field based by year;
##levels(farmData2015to2018$Product);
##[1] "Barley"  "Canary"  "Canola"  "Lentils" "Wheat" ;
cultivars <- as.data.frame(cultivars);
farmData2015to2018c <- join_all(list(farmData2015to2018, cultivars), by = c('Field', 'Product', 'ProcYear'), type = "left"); 

farmData2015to2018c$Harv_Date <- as.factor(farmData2015to2018c$Harv_Date);
farmData2015to2018c$Mat_Days <- as.integer(farmData2015to2018c$Mat_Days);
farmData2015to2018c$Seed_Date <- as.factor(farmData2015to2018c$Seed_Date);
farmData2015to2018c$Cultivar <- as.factor(farmData2015to2018c$Cultivar);


###This does data verification checks;
#farmData2015to2018c[which(is.na(farmData2015to2018c$Cultivar)), ]$Field; 
##18,799 observations with missing cultivar. DELETE, these are seeded sloughs going over in neighboring fields;
farmData2015to2018c <- farmData2015to2018c[-which(is.na(farmData2015to2018c$Cultivar)), ]; 

#farmData2015to2018c[which(is.na(farmData2015to2018c$Mat_Days)), ]$Field; 
##74,595 observations with missing Mat_Days. DELETE, these are seeded sloughs going over in neighboring fields;
farmData2015to2018c <- farmData2015to2018c[-which(is.na(farmData2015to2018c$Mat_Days)), ]; 

#farmData2015to2018c[which(is.na(farmData2015to2018c$Ino_Rate)), ]$Field; 
#farmData2015to2018c[which(is.na(farmData2015to2018c$Machine)), ]$Field; 

farmData2015to2018c$Seed_Rate <- jitter(farmData2015to2018c$Seed_Rate);
farmData2015to2018c$N_Rate <- jitter(farmData2015to2018c$N_Rate);
farmData2015to2018c$P_Rate <- jitter(farmData2015to2018c$P_Rate);
farmData2015to2018c$K_Rate <- jitter(farmData2015to2018c$K_Rate);
farmData2015to2018c$S_Rate <- jitter(farmData2015to2018c$S_Rate);
farmData2015to2018c$Ino_Rate <- jitter(farmData2015to2018c$Ino_Rate);

farmData2015to2018c[which(farmData2015to2018c$Product == "Barley"), ]$Ino_Rate <- 0;
farmData2015to2018c[which(farmData2015to2018c$Product == "Canola"), ]$Ino_Rate <- 0;
farmData2015to2018c[which(farmData2015to2018c$Product == "Canary"), ]$Ino_Rate <- 0;
farmData2015to2018c[which(farmData2015to2018c$Product == "Wheat"), ]$Ino_Rate <- 0;



#------- Comment this out 22 October 2019? No, leave it uncommented;
farmData2015to2018c[which(farmData2015to2018c$Product == "Barley"), ]$S_Rate <- 0;
farmData2015to2018c[which(farmData2015to2018c$Product == "Canary"), ]$S_Rate <- 0;
farmData2015to2018c[which(farmData2015to2018c$Product == "Lentils"), ]$S_Rate <- 0;
farmData2015to2018c[which(farmData2015to2018c$Product == "Wheat"), ]$S_Rate <- 0;

farmData2015to2018c[which(farmData2015to2018c$Product == "Canola"), ]$K_Rate <- 0;
#----------------------------------------------------;



str(farmData2015to2018c);
setwd("/Users/tylerpittman/Farm/agYieldProject/data");
save(farmData2015to2018c, file=paste("farmDataComplete2015to2018", ".RData", sep=""), compress="xz");


########## check for 900011 in 2015 #########;
#tank4_2015s <- tank4_2015[which(tank4_2015@data$SP_ID_1 == "900011"), ]; #this is correct for seed rate;
###projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0";
###proj4string(tank4_2015s) <- CRS(projection);
#yieldDataConstraint2015s <- yieldDataConstraint2015[which(yieldDataConstraint2015@data$Field == "900011"), ]; #this is correct for yield;
#plot(tank4_2015s);
#plot(yieldDataConstraint2015s, col="red", add=TRUE);
#yieldDataConstraint2015s@data$Seed_Rate4 <- over(yieldDataConstraint2015s, tank4_2015s)$Seed_Rate;