## Takes farm quintile polygons and creates big map
## Tyler Pittman, 12 March 2018
# Rscript --no-save /Users/tylerpittman/Farm/STEP_4_make_big_map_12March2018.R

## Note that field boundaries aren't exact, can't do this until after harvest 2015

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
library(devEMF); 	# Allows plots to be saved as .emf nicely in Windows
library(ggplot2);
source("shape2poly.R"); # Reads in shape2poly function and others
#library(GISTools);

northarrow <- function(loc,size,bearing=0,cols,cex=1,...) {
	# checking arguments
	if(missing(loc)) stop("loc is missing")
	if(missing(size)) stop("size is missing")
	# default colors are white and black
	if(missing(cols)) cols <- rep(c("white","black"),8)
	# calculating coordinates of polygons
	radii <- rep(size/c(1,4,2,4),4)
	x <- radii[(0:15)+1]*cos((0:15)*pi/8+bearing)+loc[1]
	y <- radii[(0:15)+1]*sin((0:15)*pi/8+bearing)+loc[2]
	# drawing polygons
	for (i in 1:15) {
		x1 <- c(x[i],x[i+1],loc[1])
		y1 <- c(y[i],y[i+1],loc[2])
		polygon(x1,y1,col=cols[i])
	}
	# drawing the last polygon
	polygon(c(x[16],x[1],loc[1]),c(y[16],y[1],loc[2]),col=cols[16])
	# drawing letters
	b <- c("E","N","W","S")
	for (i in 0:3) text((size+par("cxy")[1])*cos(bearing+i*pi/2)+loc[1],
		(size+par("cxy")[2])*sin(bearing+i*pi/2)+loc[2],b[i+1],
		cex=cex)
}

scalebar <- function(loc,length,unit="km",division.cex=.8,...) {
	if(missing(loc)) stop("loc is missing")
	if(missing(length)) stop("length is missing")
	x <- c(0,length/c(4,2,4/3,1),length*1.1)+loc[1]
	y <- c(0,length/(10*3:1))+loc[2]
	cols <- rep(c("black","white"),2)
	for (i in 1:4) rect(x[i],y[1],x[i+1],y[2],col=cols[i])
	for (i in 1:5) segments(x[i],y[2],x[i],y[3])
	#labels <- x[c(1,3)]-loc[1]
	#labels <- append(labels,paste(x[5]-loc[1],unit))
	labels <- c(0, 1)
	labels <- append(labels,paste(2,unit))
	text(x[c(1,3,5)],y[4],labels=labels,adj=.5,cex=division.cex)
}


year = 2018;

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


#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%
j <- 1;

fb = key$field[j];
fy = key$field[j];
key$Crop <- eval(parse(text=paste("key$Crop", year, sep="")));
type = key$Crop[j];
yearStartHarvest = 2014; 
#fb=800028;
#fy=800028;
fieldBoundary <- readShapeSpatial(paste("/Users/tylerpittman/Farm/boundariesSeedingCleanTopcon2017/", fb, ".shp", sep=""));
fieldBoundary.dbf <- read.dbf(paste("/Users/tylerpittman/Farm/boundariesSeedingCleanTopcon2017/", fb, ".dbf", sep=""), header=TRUE);
setwd(paste("/Users/tylerpittman/Farm/contour_yields_ll/", sep=""));
fieldPoly <- readOGR(dsn=".", layer=paste(fy, "_SeedingRx_ll_", year, sep=""));
fieldPoly.dbf <- read.dbf(paste(fy, "_SeedingRx_ll_", year, ".dbf", sep=""), header=TRUE);
fieldPoly <- fieldPoly[ , c(which(names(fieldPoly) %in% c("Yield", "SEASON", "COACH", "CROP", "ID")))];

#####
##Makes field 900010 have the same number of columns (18) as other fields 
#j <- 10;
#fb = key$field[j];
#fy = key$field[j];
#type = key$Crop[j];
#year = 2016;
#fieldPoly1 <- readOGR(dsn=".", layer=paste(fy, "_SeedingRx_ll_", year, sep=""));
#fieldPoly1.dbf <- read.dbf(paste(fy, "_SeedingRx_ll_", year, ".dbf", sep=""), header=TRUE);
##fieldPoly1 <- fieldPoly1[,-c(1, 4)]
#names(fieldPoly1)[names(fieldPoly1)=="SP_ID_1_1"] <- "SP_ID";
#fieldPoly1@data$SP_ID <- c(1:5);
#writePolyShape(fieldPoly1, paste("/Users/tylerpittman/Farm/contour_yields_ll/", fb, "_SeedingRx_ll_", year, ".shp", sep=''));
#####

#fieldYield@coords;
#str(fieldYield@data);
#fieldYield@data$DryYield;

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

#projected <- "+proj=utm +zone=13 +ellps=WGS84";
projected <- "+proj=aea +lat_1=50 +lat_2=51 +lat_0=50.5 +lon_0=-108 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs";
proj4string(fieldBoundary) <- CRS("+proj=longlat");
proj4string(fieldPoly) <- CRS("+proj=longlat");
fieldBoundaryNew <- spTransform(fieldBoundary, CRS(projected));
fieldPolyNew <- spTransform(fieldPoly, CRS(projected));

nclrYield <- 5; 
#display.brewer.all(n=nclrYield, type="div", select=NULL, exact.n=TRUE);
plotclrYieldCanola <- brewer.pal(nclrYield,"PiYG");
plotclrYieldCanary <- brewer.pal(nclrYield,"PuOr");
plotclrYieldWheat <- brewer.pal(nclrYield,"BrBG");
plotclrYieldBarley <- brewer.pal(nclrYield,"PRGn");
plotclrYieldLentils <- brewer.pal(nclrYield,"RdBu");
plotclr <- plotclrYieldWheat;

#plot(fieldPolyNew, col=plotclr);
gCentroid(fieldBoundaryNew)

uid <- 1;
uidB <- 1;
n <- length(slot(fieldPolyNew, "polygons"));
nB <- length(slot(fieldBoundaryNew, "polygons"));
fieldPolyNew <- spChFIDs(fieldPolyNew, as.character(uid:(uid+n-1)));
fieldPolyNew@data$Field <- key$field[j];
fieldBoundaryNew <- spChFIDs(fieldBoundaryNew, as.character(uid:(uidB+nB-1)));
fieldBoundaryNew@data$Field <- key$field[j];
uid <- uid + n;
uidB <- uidB + nB;

for (j in 2:length(key$field)) {
	#j <- 1;
	fb = key$field[j];
	fy = key$field[j];

	temp.dataB <- readShapeSpatial(paste("/Users/tylerpittman/Farm/boundariesSeedingCleanTopcon2017/", fb, ".shp", sep=""));
	temp.dataB.dbf <- read.dbf(paste("/Users/tylerpittman/Farm/boundariesSeedingCleanTopcon2017/", fb, ".dbf", sep=""), header=TRUE);
	setwd(paste("/Users/tylerpittman/Farm/contour_yields_ll/", sep=""));
	temp.data <- readOGR(dsn=".", layer=paste(fy, "_SeedingRx_ll_", year, sep=""));
	temp.data.dbf <- read.dbf(paste(fy, "_SeedingRx_ll_", year, ".dbf", sep=""), header=TRUE);
	temp.data <- temp.data[ , c(which(names(temp.data) %in% c("Yield", "SEASON", "COACH", "CROP", "ID")))];

	projected <- "+proj=aea +lat_1=50 +lat_2=51 +lat_0=50.5 +lon_0=-108 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs";
	proj4string(temp.dataB) <- CRS("+proj=longlat");
	proj4string(temp.data) <- CRS("+proj=longlat");
	temp.dataB <- spTransform(temp.dataB, CRS(projected));
	temp.data <- spTransform(temp.data, CRS(projected));
	
	n <- length(slot(temp.data, "polygons"));
	nB <- length(slot(temp.dataB, "polygons"));
   	temp.data <- spChFIDs(temp.data, as.character(uid:(uid+n-1)));
	temp.data@data$Field <- key$field[j];
	temp.dataB <- spChFIDs(temp.dataB, as.character(uidB:(uidB+nB-1)));
	temp.dataB@data$Field <- key$field[j];
	uid <- uid + n;
	uidB <- uidB + nB;
  	fieldPolyNew <- spRbind(fieldPolyNew,temp.data);
	fieldBoundaryNew <- spRbind(fieldBoundaryNew,temp.dataB);

	if (key$Crop[j] == "Canary"){
		plotclr.tmp = plotclrYieldCanary; 
	} else if (key$Crop[j] == "Wheat"){
		plotclr.tmp = plotclrYieldWheat;
	} else if (key$Crop[j] == "Barley"){
		plotclr.tmp = plotclrYieldBarley;
	} else if (key$Crop[j] == "Lentils"){
		plotclr.tmp = plotclrYieldLentils; 
	} else {
		plotclr.tmp = plotclrYieldCanola; 
	};
	plotclr <- cbind(plotclr, plotclr.tmp);
}
plotclr <- as.vector(plotclr);
bounde <- extent(fieldPolyNew);
xmin <- bounde@xmin;
xmax <- bounde@xmax;
ymin <- bounde@ymin;
ymax <- bounde@ymax;
#plot(fieldPolyNew, col=plotclr, border=F);
#plot.window( xlim=c(xmin, xmax), ylim=c(ymin, ymax) );

par(mar=c(0,0,0,0), oma=c(1,0,5,0));
plot(fieldPolyNew, col=plotclr, border=F, xlim=c(xmin, xmax), ylim=c(ymin, ymax));
plot(fieldBoundaryNew, add=TRUE);
shadowtext(coordinates(fieldBoundaryNew), labels=fieldBoundaryNew@data$Field, bg="white", col="black", cex=0.50);
#text(coordinates(fieldBoundaryNew), labels=fieldBoundaryNew@data$ADSFLDID);
#text(gCentroid(fieldBoundaryNew)@coords, labels=fieldBoundaryNew@data$ADSFLDID);
#legend("bottomright", paste("Bushel Yield Weighting - Canary: ", canaryW, ", Wheat: ", wheatW, ", Barley: ", barleyW, ", Lentils: ", lentilW, ", Canola: ", canolaW, sep=""), cex = 0.75, bty="n");
legend("bottomright", paste("Aggregation based on combining yield quintiles from each harvest with median imputation.", sep=""), cex = 0.75, bty="n");
title(paste("Rx Mid Row Bander Seeding Rates, Spring ", year, " (in Lbs/Acre) from Harvest ", yearStartHarvest, " to ", year-1, " Yield", sep=""), outer=TRUE);
legend("bottomleft", legend=c("120", "122.5", "125", "127.5", "130", "170", "172.5", "175", "177.5", "180", "105", "107.5", "110", "112.5", "115", "35", "37.5", "40", "42.5", "45", "175", "177.5", "180", "182.5", "185"), fill=c(plotclrYieldCanary, plotclrYieldWheat, plotclrYieldBarley, plotclrYieldLentils, plotclrYieldCanola), title="Canary      Wheat       Barley       Lentils       Canola\n 46-0-0-0    46-0-0-0    46-0-0-0     8-38-16-0     46-0-0-0", cex=1.0, bg="white", bty="n", ncol=5);
box(lty = '1373', col = 'black');
length(fieldBoundaryNew@polygons)

## this spreads out labels so no overlap on plot, uses y coordinate! ##
tmp.y <- spread.labs(coordinates(fieldBoundaryNew)[,2], 175, maxiter=1000, min=0, max=100);
spreadLabel <- cbind(coordinates(fieldBoundaryNew)[,1], tmp.y);

par(mar=c(0,0,0,0), oma=c(1,0,5,0));
plot(fieldPolyNew, col=plotclr, border=F, xlim=c(xmin, xmax), ylim=c(ymin, ymax));
plot(fieldBoundaryNew, add=TRUE);
shadowtext(spreadLabel, labels=fieldBoundaryNew@data$Field, bg="white", col="black", cex=0.50);
#legend("bottomright", paste("Bushel Yield Weighting - Canary: ", canaryW, ", Wheat: ", wheatW, ", Barley: ", barleyW, ", Lentils: ", lentilW, ", Canola: ", canolaW, sep=""), cex = 0.75, bty="n");
legend("bottomright", paste("Aggregation based on combining yield quintiles from each harvest with median imputation.", sep=""), cex = 0.75, bty="n");
title(paste("Rx Mid Row Bander Seeding Rates, Spring ", year, " (in Lbs/Acre) from Harvest ", yearStartHarvest, " to ", year-1, " Yield", sep=""), outer=TRUE);
legend("bottomleft", legend=c("120", "122.5", "125", "127.5", "130", "170", "172.5", "175", "177.5", "180", "105", "107.5", "110", "112.5", "115", "35", "37.5", "40", "42.5", "45", "175", "177.5", "180", "182.5", "185"), fill=c(plotclrYieldCanary, plotclrYieldWheat, plotclrYieldBarley, plotclrYieldLentils, plotclrYieldCanola), title="Canary      Wheat       Barley       Lentils       Canola\n 46-0-0-0    46-0-0-0    46-0-0-0     8-38-16-0     46-0-0-0", cex=1.0, bg="white", bty="n", ncol=5);
box(lty = '1373', col = 'black');

################

svg(paste("/Users/tylerpittman/Farm/", "SeedingRx_Midrow_", year, ".svg", sep=''),width=10.5,height=8, family="times",)
par(mar=c(0,0,0,0), oma=c(1,0,5,0));
plot(fieldPolyNew, col=plotclr, border=F, xlim=c(xmin, xmax), ylim=c(ymin, ymax));
plot(fieldBoundaryNew, add=TRUE);
##map.scale(-1500,40500,4000,"Kilometers",4,1,sfcol='red');
#map.scale(-1500,40500,3218.69,"Miles",2,1,sfcol='red');
#north.arrow(-1500,41000,75,lab="N",col='black');
northarrow(c(-1500,40500),500,bearing=0,cex=1);
scalebar(loc=c(-3000,39000),length=3218.69,unit="Miles",division.cex=.8); #make sure to set right values in scalebar() function above for two miles;
shadowtext(spreadLabel, labels=fieldBoundaryNew@data$Field, bg="white", col="black", cex=0.50);
#text(coordinates(fieldBoundaryNew), labels=fieldBoundaryNew@data$ADSFLDID);
#text(gCentroid(fieldBoundaryNew)@coords, labels=fieldBoundaryNew@data$ADSFLDID);
#legend("bottomright", paste("Bushel Yield Weighting - Canary: ", canaryW, ", Wheat: ", wheatW, ", Barley: ", barleyW, ", Lentils: ", lentilW, ", Canola: ", canolaW, sep=""), cex = 0.75, bty="n");
legend("bottomright", paste("Aggregation based on combining yield quintiles from each harvest with median imputation.", sep=""), cex = 0.75, bty="n");
title(paste("Rx Mid Row Bander Seeding Rates, Spring ", year, " (in Lbs/Acre) from Harvest ", yearStartHarvest, " to ", year-1, " Yield", sep=""), outer=TRUE);
legend("bottomleft", legend=c("120", "122.5", "125", "127.5", "130", "170", "172.5", "175", "177.5", "180", "105", "107.5", "110", "112.5", "115", "35", "37.5", "40", "42.5", "45", "175", "177.5", "180", "182.5", "185"), fill=c(plotclrYieldCanary, plotclrYieldWheat, plotclrYieldBarley, plotclrYieldLentils, plotclrYieldCanola), title="Canary      Wheat       Barley       Lentils       Canola\n 46-0-0-0    46-0-0-0    46-0-0-0     8-38-16-0     46-0-0-0", cex=1.0, bg="white", bty="n", ncol=5);
box(lty = '1373', col = 'black');
dev.off()


png(paste("/Users/tylerpittman/Farm/", "SeedingRx_Midrow_", year, ".png", sep=''),width=10000,height=8000)
par(mar=c(0,0,0,0), oma=c(1,0,5,0));
plot(fieldPolyNew, col=plotclr, border=F, xlim=c(xmin, xmax), ylim=c(ymin, ymax));
plot(fieldBoundaryNew, add=TRUE);
northarrow(c(-1500,40500),500,bearing=0,cex=1);
scalebar(loc=c(-3000,39000),length=3218.69,unit="Miles",division.cex=.8); #make sure to set right values in scalebar() function above for two miles;
shadowtext(spreadLabel, labels=fieldBoundaryNew@data$Field, bg="white", col="black", cex=0.25);
#text(coordinates(fieldBoundaryNew), labels=fieldBoundaryNew@data$ADSFLDID);
#text(gCentroid(fieldBoundaryNew)@coords, labels=fieldBoundaryNew@data$ADSFLDID);
#legend("bottomright", paste("Bushel Yield Weighting - Canary: ", canaryW, ", Wheat: ", wheatW, ", Barley: ", barleyW, ", Lentils: ", lentilW, ", Canola: ", canolaW, sep=""), cex = 0.75, bty="n");
legend("bottomright", paste("Aggregation based on combining yield quintiles from each harvest with median imputation.", sep=""), cex = 0.75, bty="n");
title(paste("Rx Mid Row Bander Seeding Rates, Spring ", year, " (in Lbs/Acre) from Harvest ", yearStartHarvest, " to ", year-1, " Yield", sep=""), outer=TRUE);
legend("bottomleft", legend=c("120", "122.5", "125", "127.5", "130", "170", "172.5", "175", "177.5", "180", "105", "107.5", "110", "112.5", "115", "35", "37.5", "40", "42.5", "45", "175", "177.5", "180", "182.5", "185"), fill=c(plotclrYieldCanary, plotclrYieldWheat, plotclrYieldBarley, plotclrYieldLentils, plotclrYieldCanola), title="Canary      Wheat       Barley       Lentils       Canola\n 46-0-0-0    46-0-0-0    46-0-0-0     8-38-16-0     46-0-0-0", cex=1.0, bg="white", bty="n", ncol=5);
box(lty = '1373', col = 'black');
dev.off()

################################
## This is temporary fix to check 900030 yield from Harvest 2014 with boundary
#fieldBoundary3 <- readShapeSpatial(paste("/Users/tylerpittman/Farm/input/", "900030", "_2015_boundary", ".shp", sep=""));
#fieldBoundary3.dbf <- read.dbf(paste("/Users/tylerpittman/Farm/input/", "900030", "_2015_boundary", ".dbf", sep=""), header=TRUE);
#setwd("/Users/tylerpittman/Farm/contour_yields_ll");
#fieldPoly3 <- readOGR(dsn=".", layer=paste("900030", "_SeedingRx_proj_2015", sep=""));
#fieldPoly3.dbf <- read.dbf(paste("900030", "_SeedingRx_proj_2015", ".dbf", sep=""), header=TRUE);
#
#projected <- "+proj=longlat";
#proj4string(fieldPoly3) <- CRS("+proj=utm +zone=13 +ellps=WGS84");
#fieldPoly4 <- spTransform(fieldPoly3, CRS(projected));
#
#plot(fieldPoly4, col=plotclrYieldWheat, border=F);
#
#writePolyShape(fieldPoly4, paste("/Users/tylerpittman/Farm/contour_yields_ll/", "900030", "_SeedingRx_ll_2015", ".shp", sep=''));

#fieldPoly3 <- readOGR(dsn=".", layer=paste("900030", "_SeedingRx_ll_2015", sep=""));
#fieldPoly3.dbf <- read.dbf(paste("900030", "_SeedingRx_ll_2015", ".dbf", sep=""), header=TRUE);
#head(fieldPoly3@data,3);
#fieldPoly6 <- fieldPoly3[,-(2)];
#head(fieldPoly6@data,3);
#writePolyShape(fieldPoly6, paste("/Users/tylerpittman/Farm/contour_yields_ll/", "900030", "_SeedingRx_ll_2015", ".shp", sep=''));
#fieldPoly7 <- readOGR(dsn=".", layer=paste("900030", "_SeedingRx_ll_2015", sep=""));
################################


#p <- ggplot() +
#    geom_polygon(data = fieldPolyNew, aes(x = long, y = lat, group = group),
#        fill = NA, color = "black", size = 0.25) +
#    coord_map()
#ggsave(p, file = "/Users/tylerpittman/Farm/map3.png", width = 5, height = 4.5, type = "cairo-png")
#
#
#
#county <- readOGR(dsn = ".", layer = "gz_2010_13_060_00_500k")
#county <- fortify(county, region="COUNTY")
#
#p <- ggplot() +
#    geom_polygon(data = plotData, aes(x = long, y = lat, group = group,
#        fill = percent)) +
#    geom_polygon(data = county, aes(x = long, y = lat, group = group),
#        fill = NA, color = "black", size = 0.25) +
#    coord_map()
#ggsave(p, file = "/Users/tylerpittman/Farm/map3.png", width = 5, height = 4.5, type = "cairo-png")
#

#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%#@%

