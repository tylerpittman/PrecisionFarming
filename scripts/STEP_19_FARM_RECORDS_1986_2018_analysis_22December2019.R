## Analysis of farm records and rainfall 1986 to 2018 
## Tyler Pittman, 22 December 2019
# Rscript --no-save /Users/tylerpittman/Farm/STEP_19_FARM_RECORDS_1986_2018_analysis_22December2019.R

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
library(stargazer); 	#use to create nice summary table
library(tidyr);
###library(reshape2);	#gives melt() function
library(postMCMCglmm); #need for ranef() function, may not install on cedar;
library(reporttools);	#gives tableContinous() summary function;
library(gridExtra);		#gives grid.arrange() to stack ggplots;
library(doBy);			#gives summaryBy() function for descriptive statistics;
library(splitstackshape); #gives stratified() function for stratified random sample;
library(readxl);
library(parallel);	# Does multicore processing when using mclapply WORKS on Linux!
#library(greenbrown); #has many complicated dependencies to install, for equal area parcel division of shapefiles;
library(SDMTools);  #gives scalebar() function;
library(PCICt);		#gives growing.season.length() function;
#source("shape2poly.R"); # Reads in shape2poly function and others
#source("polygonizer.R"); # Reads in polygonizer function
source("http://www.math.mcmaster.ca/bolker/R/misc/legendx.R"); # custom source allows box.cex() function to be used in legend() for plot!;



#-----------------------------------------!!!!!-----------------------------------------#;
#-----------------------------------------!!!!!-----------------------------------------#;
#-------- Does field record map for reference --------#;

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

scalebars <- function(loc,length,unit="km",division.cex=.8,...) {
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


fieldBoundary <- readShapeSpatial(paste("/Users/tylerpittman/Farm/agYieldProject/", "fieldBoundaries4269.shp", sep=""));
fieldBoundary.dbf <- read.dbf(paste("/Users/tylerpittman/Farm/agYieldProject/", "fieldBoundaries4269.dbf", sep=""), header=TRUE);

projected <- "+proj=aea +lat_1=50 +lat_2=51 +lat_0=50.5 +lon_0=-108 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs";
proj4string(fieldBoundary) <- CRS("+proj=longlat +datum=NAD83");
fieldBoundaryNew <- spTransform(fieldBoundary, CRS(projected));

bounde <- extent(fieldBoundaryNew);
xmin <- bounde@xmin;
xmax <- bounde@xmax;
ymin <- bounde@ymin;
ymax <- bounde@ymax;

## this spreads out labels so no overlap on plot, uses y and x coordinates! ##
#tmp.y <- spread.labs(coordinates(fieldBoundaryNew)[,2], 75, maxiter=1000, min=-Inf, max=Inf);
#tmp.x <- spread.labs(coordinates(fieldBoundaryNew)[,1], 75, maxiter=1000, min=-Inf, max=Inf);
#spreadLabel <- cbind(tmp.x, tmp.y);
#spreadLabel1 <- cbind(tmp.x, tmp.y-100);
#spreadLabel2 <- cbind(tmp.x, tmp.y+100);

projection(fieldBoundaryNew); #units are meters;
plot(fieldBoundaryNew, border=T, xlim=c(xmin, xmax), ylim=c(ymin, ymax));
labPosition <- polygonsLabel(fieldBoundaryNew, labels=sub(".*0{2,}", "", fieldBoundaryNew@data$SP_ID_1), method="inpolygon", gridpoints=60, polypart="largest", cex=0.9, doPlot=TRUE);
dev.off()

#tmp.y <- spread.labs(labPosition[,2], mindiff=75, stepsize=1/10, maxiter=1000, min=-Inf, max=Inf);
#spreadLabelNM <- cbind(labPosition[,1], tmp.y);

## manual label spacing ##;
customLabs <- cbind(labPosition[,1], labPosition[,2], sub(".*0{2,}", "", fieldBoundaryNew@data$SP_ID_1));
#customLabs <- cbind(labPosition[,1], tmp.y, sub(".*0{2,}", "", fieldBoundaryNew@data$SP_ID_1));

key <- matrix(c( 
"-6500","34453.49954443","1b", 
"-7616.64179282","36753.94847462","21", 
"-6797.12456845","36687.84602203","19", 
"-8027.20587613","36735.38037611","20", 
"-8520.2470076","39800","17a",
"-7202.25036497","36813.09094767","18",
"-8929.33050734","40200","16a",
"-14946.7858769","35959.24115435","15", 
"-14219.89073122","36700","13a",
"-15346.34036346","35413.35535664","14b",
"-14374.26791539","37595.84759138","11", 
"-14638.46627995","37000","12a",
"-11794.33273255","36700","9a",
"-8125.74457399","34600","7a", 
"-15001.68292238","37546.83941591","10b",
"-6201.11312053","34950","5b", 
"-12190.53418218","37000","8a", 
"-5177.36434342","35172.10361052","6b", 
"-6500","33643.14088482","3b", 
"-6500","33249.10256493","4b", 
"-7100","34234.31657129","2a", 
"-6500","34025.38484509","2b", 
"-7100","34626.12615784","1a", 
"-7100","33423.07825264","4a", 
"-7100","33838.39313313","3a", 
"-6200.55743","35400","5a", 
"-7707.62905757","34600","7c", 
"-7913.86116682","34200","7b", 
"-7509.54281884","34200","7d", 
"-8723.31452092","39400","16b",
"-8316.00725717","39200","17b",
"-11992.60151124","35800","8b",
"-11587.87151196","36200","9b", 
"-14027.1823015","36200","13b",
"-14423.43660903","35800","12b",
"-15340.04210106","36566.43583032","14a",
"-15357.63351346","37800","10a",
"-5581.21201278","35376.2271084","6a" 
), ncol=3, byrow=T);
key <- as.data.frame(key);
colnames(key) <- c("x", "y", "field");
key$x <- as.integer(levels(key$x))[key$x];
key$y <- as.integer(levels(key$y))[key$y];
labPosition <- cbind(key$x, key$y);

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("/Users/tylerpittman/Farm/agYieldProject/", "FieldRecord_Boundaries_Manuscript_1986to2018", ".pdf", sep=''), bg="transparent", width=7, height=5, pointsize=12, family="FreeSans");
par(mar=c(0,0,0,0), oma=c(0,0,0,0));
plot(fieldBoundaryNew, border=T, xlim=c(xmin, xmax), ylim=c(ymin, ymax));
northarrow(c(-6000,39300),500,bearing=0,cex=1.1);
scalebars(loc=c(-15000,34200),length=3218.69,unit="mi",division.cex=1.0); #make sure to set right values in scalebar() function above for two miles;
Scalebar(x=-15000, y=33500, distance=4000, unit="km", scale=0.001, t.cex=1.0);
#shadowtext(spreadLabel, labels=sub(".*0{2,}", "", fieldBoundaryNew@data$SP_ID_1), bg="white", col="black", cex=0.9);
shadowtext(labPosition, labels=sub(".*0{2,}", "", fieldBoundaryNew@data$SP_ID_1), bg="white", col="black", cex=0.9);
#shadowtext(spreadLabelNM, labels=sub(".*0{2,}", "", fieldBoundaryNew@data$SP_ID_1), bg="white", col="black", cex=0.9);
##polygonsLabel(fieldBoundaryNew, labels=sub(".*0{2,}", "", fieldBoundaryNew@data$SP_ID_1), method="random", gridpoints=10, polypart="largest", cex=0.9);
#shadowtext(spreadLabel1, labels=fieldBoundaryNew@data$AreaRecA, bg="white", col="blue", cex=0.25);
#shadowtext(spreadLabel2, labels=fieldBoundaryNew@data$AreaAcr, bg="white", col="red", cex=0.25);
#legend(x=-16000, y=33300, paste("Recorded field area (acres) in blue, actual field area in red.", sep=""), cex = 0.8, bty="n");
#box(lty = '1373', col = 'black');
dev.off()

labPosition1 <- cbind(labPosition[,1], labPosition[,2]-150);
labPosition2 <- cbind(labPosition[,1], labPosition[,2]+150);

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("/Users/tylerpittman/Farm/agYieldProject/", "FieldRecord_Boundaries_1986to2018", ".pdf", sep=''), bg="transparent", width=7, height=5, pointsize=12, family="FreeSans");
par(mar=c(0,0,0,0), oma=c(0,0,0,0));
plot(fieldBoundaryNew, border=T, xlim=c(xmin, xmax), ylim=c(ymin, ymax));
northarrow(c(-6000,39300),500,bearing=0,cex=1.1);
scalebars(loc=c(-15000,34200),length=3218.69,unit="mi",division.cex=1.0); #make sure to set right values in scalebar() function above for two miles;
Scalebar(x=-15000, y=33500, distance=4000, unit="km", scale=0.001, t.cex=1.0);
shadowtext(labPosition, labels=sub(".*0{2,}", "", fieldBoundaryNew@data$SP_ID_1), bg="white", col="black", cex=0.9);
shadowtext(labPosition1, labels=fieldBoundaryNew@data$AreaRecA, bg="white", col="blue", cex=0.4);
shadowtext(labPosition2, labels=fieldBoundaryNew@data$AreaAcr, bg="white", col="red", cex=0.4);
legend(x=-16000, y=33300, paste("Recorded field area (acres) in blue, actual field area in red.", sep=""), cex = 0.8, bty="n");
#box(lty = '1373', col = 'black');
dev.off()
#-----------------------------------------!!!!!-----------------------------------------#;
#-----------------------------------------!!!!!-----------------------------------------#;




#-----------------------------------------!!!!!-----------------------------------------#;
#-----------------------------------------!!!!!-----------------------------------------#;
#-------- Rainfall and temperature data for farm and swift Current 1986 to 2018 --------#;

##---##---##---##---##---##---##---##---##---##---##---##---##;
### This loads Weather Office data;
##http://climate.weather.gc.ca/climate_data/daily_data_e.html?hlyRange=2011-06-27%7C2019-07-14&dlyRange=2011-06-30%7C2019-07-13&mlyRange=%7C&StationID=48976&Prov=SK&urlExtension=_e.html&searchType=stnName&optLimit=yearRange&StartYear=1840&EndYear=2019&selRowPerPage=25&Line=1&searchMethod=contains&Month=7&Day=16&txtStationName=Swift+current&timeframe=2&Year=2019

################################################################################;
# This for loop is important for open weather .csvs for different sites;
yearK <- as.vector(as.numeric(1986):as.numeric(2018));
output <- matrix(NA, 1, length(c("date", "month", "day", "maxTemp", "minTemp", "meanTemp", "totalPrecip", "maxWinDir", "maxWindSpeed", "Station")));
output <- as.data.frame(output);
colnames(output) <- c("date", "month", "day", "maxTemp", "minTemp", "meanTemp", "totalPrecip", "maxWinDir", "maxWindSpeed", "Station");
output <- output[-1,];
for(i in 1:as.numeric(2018-1986+1)){
	if (yearK[i] >= 2018){
		tmp <- read.csv(paste("sc", yearK[i], ".csv", sep=""), fileEncoding="latin1", header=F, skip=26);
	} else {
		tmp <- read.csv(paste("sc", yearK[i], ".csv", sep=""), fileEncoding="latin1", header=F, skip=25);
	};
	tmp <- tmp[ , c(1,3,4,6,8,10,20,24,26)];
	colnames(tmp) <- c("date", "month", "day", "maxTemp", "minTemp", "meanTemp", "totalPrecip", "maxWinDir", "maxWindSpeed");
	if (yearK[i] >= 2012){
		tmp$Station <- "71142";
	} else if (yearK[i] == 1992){
		tmp$Station <- "71870";
	} else {
		tmp$Station <- "71146";
	};
	output <- rbind(output, tmp);
	print(yearK[i]);
};
output$date <- as.Date(output$date, "%Y-%m-%d");
output$month <- as.integer(output$month);
output$day <- as.integer(output$day);
output$maxTemp <- as.numeric(output$maxTemp);
output$minTemp <- as.numeric(output$minTemp);
output$meanTemp <- as.numeric(output$meanTemp);
output$totalPrecip <- as.numeric(output$totalPrecip);
output$maxWinDir <- as.integer(output$maxWinDir);
output$maxWindSpeed <- as.integer(output$maxWindSpeed);
################################################################################;

str(output);
#output$totalPrecip;
#output$date;
#output$Station;

wo0 <- output;
wo0$Station <- as.factor(wo0$Station);
wo0$Month <- months(wo0$date);
wo0$Year <- format(wo0$date, format="%Y");
wo0$DryTempH <- wo0$maxTemp;
wo0$DryTempL <- wo0$minTemp;
wo0$DryTempA <- wo0$meanTemp;
wo0$PrecipSum <- wo0$totalPrecip;
str(wo0);

write.table(wo0, file = "SwiftCurrentWeather_1986to2018.csv", row.names=FALSE, na="", col.names=TRUE, sep=",");


### This loads farm rainfall data;
farmRF <- read_excel("rainfallFarm1986to2018.xlsx", sheet="RF");
farmRF <- as.data.frame(farmRF);
farmRF$Date <-  as.Date(paste(as.character(farmRF$Date), sep=""), format="%Y-%m-%d");
farmRF$PrecipSum <- farmRF$RF_inches * 25.4;
farmRF$Station <- "Farm";
farmRF$Station <- as.factor(farmRF$Station);
farmRF$Month <- months(farmRF$Date);
farmRF$Year <- format(farmRF$Date, format="%Y");
str(farmRF);
######;

maxTempM <- aggregate(DryTempH ~ Month + Year + Station, wo0, max);
minTempM <- aggregate(DryTempL ~ Month + Year + Station, wo0, min);
avgTempM <- aggregate(DryTempA ~ Month + Year + Station, wo0, mean);
sumPrecipM <- aggregate(PrecipSum ~ Month + Year + Station, wo0, sum);



#### Farm rainfall for each year 2015 to 2018;
sumPrecipM <- aggregate(PrecipSum ~ Month + Year, farmRF, sum);
wo4 <- join_all(list(maxTempM, minTempM, avgTempM, sumPrecipM), by = c('Month', 'Year'), type = "left"); 
annTemp <- mean(wo0$DryTempA, na.rm=TRUE);
annTemp <- as.data.frame(annTemp);
colnames(annTemp) <- c("AnnTemp");

##below takes maxTempM, minTempM and avgTempM from above section;
wo4$AnnTemp <- annTemp$AnnTemp;
wo4$Date <- as.Date(paste(as.character(wo4$Month), "/01/", as.character(wo4$Year), sep=""), format="%B/%d/%Y");
wo4 <- wo4[order(wo4$Date), ]; #sorts by Date for below cumsum;
wo4[which(is.na(wo4$PrecipSum)), ]$PrecipSum <- 0;
wo4$cummPrecip <- ave(wo4$PrecipSum, list(wo4$Year), FUN=cumsum); 
wo4[which(wo4$Month == "December"), ]$cummPrecip <- 0;
wo4[which(wo4$cummPrecip == 0), ]$cummPrecip <- NA;
str(wo4);

summary(wo4$DryTempH); #max is 39.8 C;
summary(wo4$DryTempL); #low is -41.0 C;
###-------------------------###;


###---------------------------------##;
## Calculates length of frost-free period by year;

weatherFarm <- wo0;
#str(weatherFarm);
weatherFarm$Year;
weatherFarm[which(is.na(weatherFarm$DryTempL)), ]$Year;

weatherFarm$Frost_Free <- "Yes";
weatherFarm[which(weatherFarm$DryTempL < 0), ]$Frost_Free <- "No";
#weatherFarm$Frost_Free;

growing_season1 <- weatherFarm %>% 
    group_by(run = data.table::rleid(Frost_Free), Year, Frost_Free) %>% 
    summarise(count = n())
growing_season2 <- growing_season1 %>%
	group_by(Year, Frost_Free) %>% 
	filter(Frost_Free == "Yes") %>%
	top_n(n=1, wt = count)

growing_season2$count;
growing_season2$Year;	
mean(growing_season2$count);
	
#weatherFarm_check <- weatherFarm[which(weatherFarm$Year == 1987), ]; #137 days;
###---------------------------------##;





#-----------------------------------------!!!!!-----------------------------------------#;
#-----------------------------------------!!!!!-----------------------------------------#;



#-----------------------------------------!!!!!-----------------------------------------#;
#-----------------------------------------!!!!!-----------------------------------------#;
#-------------------------- Farm records from 1986 to 2018 -----------------------------#;

### This loads farm field record data;
farmFieldRecords <- read_excel("fieldRecords1986to2018.xlsx", sheet="fieldSeedHarvRecord");
farmFieldRecords <- as.data.frame(farmFieldRecords);
colnames(farmFieldRecords);
farmFieldRecords$SeedDate <-  as.Date(paste(as.character(farmFieldRecords$SeedDate), sep=""), format="%m/%d/%Y");
farmFieldRecords$HarvDate <-  as.Date(paste(as.character(farmFieldRecords$HarvDate), sep=""), format="%m/%d/%Y");
farmFieldRecords$SwaDesDate <-  as.Date(paste(as.character(farmFieldRecords$SwaDesDate), sep=""), format="%m/%d/%Y");

#### Convert to metric units (from lbs/acre to kg/ha);
farmFieldRecords$SeedRate <- farmFieldRecords$SeedRate * 1.12085;
farmFieldRecords$NRatio <- farmFieldRecords$NRatio * 1.12085;
farmFieldRecords$PRatio <- farmFieldRecords$PRatio * 1.12085;
farmFieldRecords$KRatio <- farmFieldRecords$KRatio * 1.12085;
farmFieldRecords$SRatio <- farmFieldRecords$SRatio * 1.12085;
## from bu/acre to kg/ha;
farmFieldRecords[which(farmFieldRecords$Crop == "Wheat"), ]$Yield <- farmFieldRecords[which(farmFieldRecords$Crop == "Wheat"), ]$Yield * 67.25;
farmFieldRecords[which(farmFieldRecords$Crop == "Lentils"), ]$Yield <- farmFieldRecords[which(farmFieldRecords$Crop == "Lentils"), ]$Yield * 67.25;
farmFieldRecords[which(farmFieldRecords$Crop == "Barley"), ]$Yield <- farmFieldRecords[which(farmFieldRecords$Crop == "Barley"), ]$Yield * 53.80;
farmFieldRecords[which(farmFieldRecords$Crop == "Canola"), ]$Yield <- farmFieldRecords[which(farmFieldRecords$Crop == "Canola"), ]$Yield * 56.04;
farmFieldRecords[which(farmFieldRecords$Crop == "Canary"), ]$Yield <- farmFieldRecords[which(farmFieldRecords$Crop == "Canary"), ]$Yield * 56.04;
farmFieldRecords[which(farmFieldRecords$Crop == "Chickpea"), ]$Yield <- farmFieldRecords[which(farmFieldRecords$Crop == "Chickpea"), ]$Yield * 67.25;
farmFieldRecords[which(farmFieldRecords$Crop == "Field Pea"), ]$Yield <- farmFieldRecords[which(farmFieldRecords$Crop == "Field Pea"), ]$Yield * 67.25;
######;


### This part does summation of tonnes of different products over the years ###;
farmFieldSize <- read_excel("fieldRecords1986to2018.xlsx", sheet="fieldSizeRecord");
farmFieldSize <- as.data.frame(farmFieldSize);
farmFieldSize$Hectares <- farmFieldSize$Acres * 0.404686;
######;

farmFieldRecordsHectares <- join_all(list(farmFieldRecords, farmFieldSize), by = c('Field'), type = "left"); 
farmFieldRecordsHectares$YieldTonne <- farmFieldRecordsHectares$Yield * farmFieldRecordsHectares$Hectares / 1000;
farmFieldRecordsHectares$sumTonnesYear <- ave(farmFieldRecordsHectares$YieldTonne, list(farmFieldRecordsHectares$Year), FUN=sum); 
farmFieldRecordsHectares$sumTonnesCrop <- ave(farmFieldRecordsHectares$YieldTonne, list(farmFieldRecordsHectares$Crop), FUN=sum); 
#-----------------------------------------!!!!!-----------------------------------------#;
#-----------------------------------------!!!!!-----------------------------------------#;

farmFieldRecordsHectares[which(is.na(farmFieldRecordsHectares$SwaDesc)), ]$SwaDesc <- "SC";

str(farmFieldRecordsHectares); #770 field records over the years;
farmFieldRecordsHectares$Year <- as.factor(farmFieldRecordsHectares$Year);
farmFieldRecordsHectares$Field <- as.factor(farmFieldRecordsHectares$Field);
farmFieldRecordsHectares$Crop <- as.factor(farmFieldRecordsHectares$Crop);
farmFieldRecordsHectares$Cultivar <- as.factor(farmFieldRecordsHectares$Cultivar);
farmFieldRecordsHectares$SwaDesc <- as.factor(farmFieldRecordsHectares$SwaDesc);
farmFieldRecordsHectares$InocType <- as.factor(farmFieldRecordsHectares$InocType);
farmFieldRecordsHectares$Hail <- as.factor(farmFieldRecordsHectares$Hail);
farmFieldRecordsHectares$Flood <- as.factor(farmFieldRecordsHectares$Flood);
farmFieldRecordsHectares$Fungicide <- as.factor(farmFieldRecordsHectares$Fungicide);
farmFieldRecordsHectares$Pests <- as.factor(farmFieldRecordsHectares$Pests);
farmFieldRecordsHectares$Copper <- as.factor(farmFieldRecordsHectares$Copper);
farmFieldRecordsHectares$Drought <- as.factor(farmFieldRecordsHectares$Drought);
farmFieldRecordsHectares$TreatedSeed <- as.factor(farmFieldRecordsHectares$TreatedSeed);
farmFieldRecordsHectares$Terminated <- as.factor(farmFieldRecordsHectares$Terminated);
levels(farmFieldRecordsHectares$Field);
levels(farmFieldRecordsHectares$Crop);

farmFieldRecordsHectares$SwaDesc <- relevel(farmFieldRecordsHectares$SwaDesc, ref="SC");

###This checks for lentil crops that were straight combined, should be zero;
farmFieldRecordsHectares[which(farmFieldRecordsHectares$SwaDesc == "SC" & farmFieldRecordsHectares$Crop == "Lentils"), ];

which(is.na(farmFieldRecordsHectares$Crop)); #All yield points have a product;
summary(farmFieldRecordsHectares$Yield);

###This calculates number of days between seeding and harvest;
farmFieldRecordsHectares$Mat_Date <- farmFieldRecordsHectares$SwaDesDate;
farmFieldRecordsHectares[which(is.na(farmFieldRecordsHectares$Mat_Date)), ]$Mat_Date <- farmFieldRecordsHectares[which(is.na(farmFieldRecordsHectares$Mat_Date)), ]$HarvDate;

farmFieldRecordsHectares$Mat_Days <- farmFieldRecordsHectares$Mat_Date - farmFieldRecordsHectares$SeedDate;
farmFieldRecordsHectares$Mat_Days <- as.integer(farmFieldRecordsHectares$Mat_Days);

#--^-^--^-^--^-^--^-^--^-^--^-^--^-^--^-^--^-^--^-^--#;
### Cross-Tabulation section;
#library(descr);
#crosstab <- crosstab(farmFieldRecordsHectares$Field, farmFieldRecordsHectares$Crop, weight = NULL);
#str(crosstab);
#write.csv(crosstab, file="tmp.csv");

#crosstab <- table(farmFieldRecordsHectares$Field, farmFieldRecordsHectares$Crop);
#write.csv(crosstab, file="tmp.csv");

#library(questionr);
#crosstab <- wtd.table(farmFieldRecordsHectares$Field, farmFieldRecordsHectares$Crop, weights=farmFieldRecordsHectares$Yield);
#write.csv(crosstab, file="tmp.csv");
#--^-^--^-^--^-^--^-^--^-^--^-^--^-^--^-^--^-^--^-^--#;
#-----------------------------------------------------------------------------------------------------------------------------#;


#--------------------------------------------------------------------------------------------------------------------------------#;
#--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--#;
#--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--#;
############### This is the data wrangling section ##################;
## Merges farm field records to weather data;
str(farmFieldRecordsHectares);
str(wo0);
wo0$Year <- as.factor(wo0$Year);
#SeedDate
#HarvDate

farmFieldsWeather <- join_all(list(farmFieldRecordsHectares, wo0), by = c('Year'), type = "left"); 

#farmFieldsWeather$SeedDate; 
#farmFieldsWeather$date; 
##Check:
#farmFieldsWeather[which((farmFieldsWeather$Field == "3a") & (farmFieldsWeather$Year == 2010)), ]; #is complete;

##https://rdrr.io/cran/dplyr/man/filter.html
mydataMtmp1 <- farmFieldsWeather %>% group_by(Year, Field) %>% filter(date >= mean(SeedDate, na.rm = FALSE));
#mydataMtmp1$date;
farmFieldsWeatherGS <- mydataMtmp1 %>% group_by(Year, Field) %>% filter(date <= mean(Mat_Date, na.rm = FALSE));
#farmFieldsWeatherGS$date;
#str(farmFieldsWeatherGS);

## Have to use custom function for below to handle na values;
##https://stackoverflow.com/questions/31753005/using-ave-in-r-without-na-values
farmFieldsWeatherGS$RF_sum <- ave(farmFieldsWeatherGS$PrecipSum, list(farmFieldsWeatherGS$Year, farmFieldsWeatherGS$Field), FUN=function(x) sum(x, na.rm=T)); 
farmFieldsWeatherGS$Ave_Temp <- ave(farmFieldsWeatherGS$meanTemp, list(farmFieldsWeatherGS$Year, farmFieldsWeatherGS$Field), FUN=function(x) mean(x, na.rm=T)); 
farmFieldsWeatherGS$Min_Temp <- ave(farmFieldsWeatherGS$minTemp, list(farmFieldsWeatherGS$Year, farmFieldsWeatherGS$Field), FUN=function(x) min(x, na.rm=T));
farmFieldsWeatherGS$Max_Temp <- ave(farmFieldsWeatherGS$maxTemp, list(farmFieldsWeatherGS$Year, farmFieldsWeatherGS$Field), FUN=function(x) max(x, na.rm=T));  

#################;
tmpCheck <- farmFieldsWeatherGS[which((farmFieldsWeatherGS$Field == "3a") & (farmFieldsWeatherGS$Year == 2010)), ]; #is complete;
tmpCheck <- as.data.frame(tmpCheck); #is complete, one NA in temp on last day of growing season;
#################;

farmRecordsWeatherGS <- farmFieldsWeatherGS %>% group_by(Year, Field) %>% filter(date == mean(SeedDate, na.rm = FALSE));
farmRecordsWeatherGS <- as.data.frame(farmRecordsWeatherGS);

farmRecordsWeatherGS <- farmRecordsWeatherGS[, c("Year", "Field", "Crop", "Cultivar", "SeedRate", "SeedDate", "HarvDate", "Mat_Date", "SwaDesc", "SwaDesDate", "Yield", "NRatio", "PRatio", "KRatio", "SRatio", "InocType", "Hail", "Flood", "Fungicide", "Pests", "Copper", "Drought", "TreatedSeed", "Terminated", "Acres", "Hectares", "YieldTonne", "sumTonnesYear", "sumTonnesCrop", "Mat_Days", "RF_sum", "Ave_Temp", "Min_Temp", "Max_Temp")];
farmRecordsWeatherGSpart <- farmRecordsWeatherGS[, c("Year", "Field", "RF_sum", "Ave_Temp", "Min_Temp", "Max_Temp")];


##Join back and keep summerfallow field information;
joinComplete <- join_all(list(farmFieldRecordsHectares, farmRecordsWeatherGSpart), by = c('Year', 'Field'), type = "left"); 
farmRecordsWeatherGS <- joinComplete;

farmRecordsWeatherGS$SeedRate <- round(farmRecordsWeatherGS$SeedRate, digits=1);
farmRecordsWeatherGS$Yield <- round(farmRecordsWeatherGS$Yield, digits=1);
farmRecordsWeatherGS$NRatio <- round(farmRecordsWeatherGS$NRatio, digits=1);
farmRecordsWeatherGS$PRatio <- round(farmRecordsWeatherGS$PRatio, digits=1);
farmRecordsWeatherGS$KRatio <- round(farmRecordsWeatherGS$KRatio, digits=1);
farmRecordsWeatherGS$SRatio <- round(farmRecordsWeatherGS$SRatio, digits=1);
farmRecordsWeatherGS$Acres <- round(farmRecordsWeatherGS$Acres, digits=2);
farmRecordsWeatherGS$Hectares <- round(farmRecordsWeatherGS$Hectares, digits=2);
farmRecordsWeatherGS$YieldTonne <- round(farmRecordsWeatherGS$YieldTonne, digits=1);
farmRecordsWeatherGS$sumTonnesYear <- round(farmRecordsWeatherGS$sumTonnesYear, digits=1);
farmRecordsWeatherGS$sumTonnesCrop <- round(farmRecordsWeatherGS$sumTonnesCrop, digits=1);
farmRecordsWeatherGS$Ave_Temp <- round(farmRecordsWeatherGS$Ave_Temp, digits=1);

##sort, remove NAs and write to .csv file;
farmRecordsWeatherGS <- farmRecordsWeatherGS[order(farmRecordsWeatherGS$Year, farmRecordsWeatherGS$Field), ];
#farmRecordsWeatherGS[which(farmRecordsWeatherGS$SeedRate == 0), ]$SeedRate <- NA;
#farmRecordsWeatherGS[which(farmRecordsWeatherGS$NRatio == 0), ]$NRatio <- NA;
#farmRecordsWeatherGS[which(farmRecordsWeatherGS$PRatio == 0), ]$PRatio <- NA;
#farmRecordsWeatherGS[which(farmRecordsWeatherGS$KRatio == 0), ]$KRatio <- NA;
#farmRecordsWeatherGS[which(farmRecordsWeatherGS$SRatio == 0), ]$SRatio <- NA;
#farmRecordsWeatherGS[which(farmRecordsWeatherGS$Yield == 0), ]$Yield <- NA;
write.table(farmRecordsWeatherGS, file="farmFieldRecordsWeather_1986to2018.csv", row.names=FALSE, na="", col.names=TRUE, sep=",");
#--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--#;
#--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--(%-%)--#;
#--------------------------------------------------------------------------------------------------------------------------------#;


#-()-()-()-()-()-()-()-()-()-()-()-()-()-()-()-#;
#-()-()-()-()-()-()-()-()-()-()-()-()-()-()-()-#;
### Determine the previous crop;
#str(farmRecordsWeatherGS);
smallRec <- farmRecordsWeatherGS[, c("Year", "Field", "Crop")];
#t(smallRec);

##https://dplyr.tidyverse.org/reference/lead-lag.html
#right <- mutate(smallRec, prev=lag(Crop, order_by=Field))
#arrange(right, Field)

##This does the correct way of getting previous year value for crop!;
cropPY <- smallRec %>%
   group_by(Field) %>%
   mutate(prevCrop=lag(Crop, order_by=Field)) %>%
   as.data.frame()

cropPY <- cropPY[, c("Year", "Field", "prevCrop")];

tmpJoin <- join_all(list(farmRecordsWeatherGS, cropPY), by=c("Year", "Field"), type = "left");
farmRecordsWeatherGS <- tmpJoin;

write.table(farmRecordsWeatherGS, file="farmFieldRecordsWeather_1986to2018.csv", row.names=FALSE, na="", col.names=TRUE, sep=",");
#-()-()-()-()-()-()-()-()-()-()-()-()-()-()-()-#;
#-()-()-()-()-()-()-()-()-()-()-()-()-()-()-()-#;


#str(farmRecordsWeatherGS);

priorFarm1 <- list(
  R=list(V=1, nu = 0),
  G=list(G1=list(V        = diag(1),
                 nu       = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),
         G2=list(V        = diag(1),
                 nu       = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),
         G3=list(V        = diag(1),
                 nu       = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),
		 G4=list(V        = diag(1),
                 nu       = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2)));
                 
priorFarm1_int <- list(
  R=list(V=1, nu = 0),
  G=list(G1=list(V        = diag(1),
                 nu       = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),
         G2=list(V        = diag(1),
                 nu       = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),
         G3=list(V        = diag(1),
                 nu       = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),
		 G4=list(V        = diag(1),
                 nu       = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),
         G5=list(V 	      = diag(2), 
         		 nu 	  = 1,
         		 alpha.mu = rep(0, 1),
                alpha.V   = diag(1)*(0)^2))
         );                 
                 
priorFarm2 <- list(
  R=list(V=1, nu = 0),
  G=list(G1=list(V        = diag(1),
                 nu       = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),
         G2=list(V        = diag(1),
                 nu       = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),
		 G3=list(V        = diag(1),
                 nu       = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2)));
                 
priorFarm3 <- list(
  R=list(V=1, nu = 0),
  G=list(G1=list(V        = diag(1),
                 nu       = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),
		 G2=list(V        = diag(1),
                 nu       = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2)));
                 
priorFarm4 <- list(
  R=list(V=1, nu = 0),
  G=list(G1=list(V        = diag(1),
                 nu       = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2)));

priorFarm5 <- list(
  R=list(V=1, nu = 0));
  
#########;
### Below does Table 1 cultivar types:
farmRecordsWeatherGS1 <- farmRecordsWeatherGS[which(farmRecordsWeatherGS$Crop == "Barley"), ];
unique(farmRecordsWeatherGS1$Cultivar);

farmRecordsWeatherGS1 <- farmRecordsWeatherGS[which(farmRecordsWeatherGS$Crop == "Canary"), ];
unique(farmRecordsWeatherGS1$Cultivar);

farmRecordsWeatherGS1 <- farmRecordsWeatherGS[which(farmRecordsWeatherGS$Crop == "Canola"), ];
unique(farmRecordsWeatherGS1$Cultivar);

farmRecordsWeatherGS1 <- farmRecordsWeatherGS[which(farmRecordsWeatherGS$Crop == "Chickpea"), ];
unique(farmRecordsWeatherGS1$Cultivar);

farmRecordsWeatherGS1 <- farmRecordsWeatherGS[which(farmRecordsWeatherGS$Crop == "Field Pea"), ];
unique(farmRecordsWeatherGS1$Cultivar);

farmRecordsWeatherGS1 <- farmRecordsWeatherGS[which(farmRecordsWeatherGS$Crop == "Lentils"), ];
unique(farmRecordsWeatherGS1$Cultivar);

farmRecordsWeatherGS1 <- farmRecordsWeatherGS[which(farmRecordsWeatherGS$Crop == "Wheat"), ];
unique(farmRecordsWeatherGS1$Cultivar);
#########;

#---\/--\/--\/--\/--\/--\/--\/--\/--\/--\/--\/--\/--\/--\/--\/--\/--\/--\/---#;
### 2. Area Description section summary information;
str(farmRecordsWeatherGS);
summary(farmRecordsWeatherGS$Ave_Temp);

summary(farmRecordsWeatherGS$Min_Temp);
farmRecordsWeatherGS[which(farmRecordsWeatherGS$Min_Temp == -10.0), ]; #year 2016;

summary(farmRecordsWeatherGS$Max_Temp);
farmRecordsWeatherGS[which(farmRecordsWeatherGS$Max_Temp == 39.8), ]; #year 2018;

summary(farmRecordsWeatherGS$RF_sum);
farmRecordsWeatherGS[which(farmRecordsWeatherGS$RF_sum == 46.9), ]; #year 2018;
farmRecordsWeatherGS[which(farmRecordsWeatherGS$RF_sum == 515.2), ]; #year 2010;

length(unique(farmRecordsWeatherGS$Field)); #38 disjoined fields;

farmRecordsWeatherGS[!duplicated(farmRecordsWeatherGS$Field), ]$Hectares; #returns field size of first occurrence (by year) of each field; 

sum(farmRecordsWeatherGS[which(farmRecordsWeatherGS$Year == 1986), ]$Hectares); #256.94 cultivated hectares in 1986;
sum(farmRecordsWeatherGS[which(farmRecordsWeatherGS$Year == 2018), ]$Hectares); #1276.65 cultivated hectares in 2018;
sum(farmRecordsWeatherGS[!duplicated(farmRecordsWeatherGS$Field), ]$Hectares); #1276.65 cultivated hectares in 2018;
#---\/--\/--\/--\/--\/--\/--\/--\/--\/--\/--\/--\/--\/--\/--\/--\/--\/--\/---#;

####-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----####;
####-----@-----@-----@-----@--Does for all crops --@-----@-----@-----@-----@-----####;
#Statistical modeling part;
str(farmRecordsWeatherGS);

##remove field records with missing values for Yield (i.e. summerfallow);
farmRecordsWeatherGS <- farmRecordsWeatherGS[which(!is.na(farmRecordsWeatherGS$RF_sum)), ];
                 
length(farmRecordsWeatherGS$SeedRate);
summary(farmRecordsWeatherGS$Yield);
summary(farmRecordsWeatherGS$SeedRate);
summary(farmRecordsWeatherGS$NRatio);
summary(farmRecordsWeatherGS$PRatio);
summary(farmRecordsWeatherGS$KRatio);
summary(farmRecordsWeatherGS$SRatio);
summary(farmRecordsWeatherGS$InocType);
summary(farmRecordsWeatherGS$Mat_Days);
summary(farmRecordsWeatherGS$SwaDesc);
summary(farmRecordsWeatherGS$Drought);
summary(farmRecordsWeatherGS$Flood);
summary(farmRecordsWeatherGS$Hail);
summary(farmRecordsWeatherGS$Copper);
summary(farmRecordsWeatherGS$Fungicide);
summary(farmRecordsWeatherGS$Pests);
summary(farmRecordsWeatherGS$TreatedSeed);
summary(farmRecordsWeatherGS$Ave_Temp);
summary(farmRecordsWeatherGS$Max_Temp);
summary(farmRecordsWeatherGS$Min_Temp);
summary(farmRecordsWeatherGS$RF_sum);
length(which(summary(farmRecordsWeatherGS$Cultivar) > 0));
length(which(summary(farmRecordsWeatherGS$Field) > 0));
length(which(summary(farmRecordsWeatherGS$Year) > 0));
length(which(summary(farmRecordsWeatherGS$prevCrop) > 0));
##below sd() doesn't work with NA values;
sd(farmRecordsWeatherGS$SeedRate);
sd(farmRecordsWeatherGS$NRatio);
sd(farmRecordsWeatherGS$PRatio);
sd(farmRecordsWeatherGS$KRatio);
sd(farmRecordsWeatherGS$SRatio);
sd(farmRecordsWeatherGS$InocType);
sd(farmRecordsWeatherGS$Mat_Days);
sd(farmRecordsWeatherGS$SwaDesc);
sd(farmRecordsWeatherGS$Drought);
sd(farmRecordsWeatherGS$Flood);
sd(farmRecordsWeatherGS$Hail);
sd(farmRecordsWeatherGS$Copper);
sd(farmRecordsWeatherGS$Fungicide);
sd(farmRecordsWeatherGS$Pests);
sd(farmRecordsWeatherGS$TreatedSeed);
sd(farmRecordsWeatherGS$Ave_Temp);
sd(farmRecordsWeatherGS$Max_Temp);
sd(farmRecordsWeatherGS$Min_Temp);
sd(farmRecordsWeatherGS$RF_sum);

options(scipen=3); #suppresses scientific notation;

#model1 <- MCMCglmm(Yield ~ Crop + SeedRate + NRatio + PRatio + KRatio + SRatio + Mat_Days + Drought + Flood + Hail + Copper + Fungicide + Pests + TreatedSeed + Ave_Temp + Max_Temp + Min_Temp + RF_sum, random=~Cultivar + Field + Year + prevCrop, family="gaussian", data=farmRecordsWeatherGS, saveX=TRUE, verbose=TRUE, burnin=1388889, nitt=13888889, thin=6250, pr=TRUE, prior=priorFarm1);
#model1 <- MCMCglmm(Yield ~ Crop + SeedRate + NRatio + PRatio + KRatio + SRatio + Mat_Days + Drought + Flood + Hail + Copper + Fungicide + Pests + TreatedSeed + Ave_Temp + Max_Temp + Min_Temp + RF_sum, random=~Cultivar + Field + Year + prevCrop, family="gaussian", data=farmRecordsWeatherGS, saveX=TRUE, verbose=TRUE, burnin=500, nitt=2000, thin=25, pr=TRUE, prior=priorFarm1);
###model1 <- MCMCglmm(Yield ~ Crop + SeedRate + NRatio + PRatio + KRatio + SRatio + InocType + Mat_Days + Drought + Flood + Hail + Copper + Fungicide + Pests + TreatedSeed + Ave_Temp + Max_Temp + Min_Temp + RF_sum, random=~Cultivar + Field + Year + prevCrop, family="gaussian", data=farmRecordsWeatherGS, saveX=TRUE, verbose=TRUE, burnin=500, nitt=2000, thin=25, pr=TRUE, prior=priorFarm1);


#save(model1, file=paste("farmRecordsRainfall_summary_Model_1.RData", sep=""), compress="xz");
load(paste("farmRecordsRainfall_summary_Model_1.RData", sep=""));
#load(paste("STR_model_published/farmRecordsRainfall_summary_Model_1.RData", sep=""));

#https://stats.stackexchange.com/questions/307021/why-does-mcmcglmm-result-in-small-effective-sample-sizes-for-logistic-regression/310291
#autocorr.diag(model1$VCV); #use this to choose thinning interval, autocorrelation less than 0.01 ideally;

fe1 <- summary(model1);
#str(fe1$solutions);
fe1$solutions[,1] <- fe1$solutions[,1];
fe1$solutions[,2] <- fe1$solutions[,2];
fe1$solutions[,3] <- fe1$solutions[,3];

options(scipen=3); #suppresses scientific notation;
writeLines(capture.output(fe1), paste("farmRecordsRainfall_Model_1.txt", sep=""));
writeLines(capture.output(summary(model1)), paste("farmRecordsRainfall_Model_1.txt", sep=""));

### DO THIS FOR DIAGNOSTICS OF HOW CHAIN FIT with mcmcglmm;
#http://www.maths.bath.ac.uk/~jjf23/mixchange/split.html

#colnames(model1$Sol); #1:24 are explanatory variables;
EVs1 <- model1$Sol[, 1:24];
intCultivars1 <- model1$Sol[, 25:55];
intFields1 <- model1$Sol[, 56:92];
intProcYears1 <- model1$Sol[, 93:125];
intPrevCrops1 <- model1$Sol[, 126:133];

length(posterior.mode(intFields1)); #38;  #mean intercept;

### https://github.com/tmalsburg/MCMCglmm-intro
###

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farmRecordsRainfall_model1_VCV_farm.pdf", sep=''), bg="transparent", width=10, height=8, pointsize=12, family="FreeSans");
par(mfrow=c(length(colnames(model1$VCV)),2), mar=c(2,2,1,0));
plot(model1$VCV, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farmRecordsRainfall_model1_sol_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(EVs1)),2), mar=c(2,2,1,0));
plot(EVs1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farmRecordsRainfall_model1_Cultivars_farm.pdf", sep=''), bg="transparent", width=11, height=30, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intCultivars1)),2), mar=c(2,2,1,0));
plot(intCultivars1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farmRecordsRainfall_model1_Fields_farm.pdf", sep=''), bg="transparent", width=11, height=30, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intFields1)),2), mar=c(2,2,1,0));
plot(intFields1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farmRecordsRainfall_model1_ProcYear_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intProcYears1)),2), mar=c(2,2,1,0));
plot(intProcYears1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farmRecordsRainfall_model1_PrevCrop_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intPrevCrops1)),2), mar=c(2,2,1,0));
plot(intPrevCrops1, auto.layout=F);
dev.off()

colnames(model1$VCV);
ICC_Cultivar <- model1$VCV[, 1]/(rowSums(model1$VCV)); 
ICC_Field <- model1$VCV[, 2]/(rowSums(model1$VCV)); 
ICC_ProcYear <- model1$VCV[, 3]/(rowSums(model1$VCV)); 
ICC_PrevCrop <- model1$VCV[, 4]/(rowSums(model1$VCV)); 
dft1.1 <- cbind(ICC = posterior.mode(ICC_Cultivar), CI = HPDinterval(ICC_Cultivar));
dft2.1 <- cbind(ICC = posterior.mode(ICC_Field), CI = HPDinterval(ICC_Field));
dft3.1 <- cbind(ICC = posterior.mode(ICC_ProcYear), CI = HPDinterval(ICC_ProcYear));
dft4.1 <- cbind(ICC = posterior.mode(ICC_PrevCrop), CI = HPDinterval(ICC_PrevCrop));
dft1 <- rbind(dft1.1, dft2.1, dft3.1, dft4.1);
row.names(dft1) <- c("ICC_Cultivar", "ICC_Field", "ICC_ProcYear", "ICC_PrevCrop");
options(scipen=3); #suppresses scientific notation;
writeLines(capture.output(dft1), paste("farmRecordsRainfall_ICC_Model_1.txt", sep=""));

#---;
####-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----####;
####-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----####;


#### IGNORE BELOW SECTION, NOT USED IN TABLE 3 BUT IS IN TABLE 2 ####;
####-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----####;
####-----@-----@-----@-----@--- Does for barley ---@-----@-----@-----@-----@-----####;
farmRecordsWeatherGS1 <- farmRecordsWeatherGS[which(farmRecordsWeatherGS$Crop == "Barley"), ];
str(farmRecordsWeatherGS1);

##remove field records with missing values for Yield (i.e. summerfallow);
farmRecordsWeatherGS1 <- farmRecordsWeatherGS1[which(!is.na(farmRecordsWeatherGS1$RF_sum)), ];
                 
length(farmRecordsWeatherGS1$SeedRate);
summary(farmRecordsWeatherGS1$Yield);
summary(farmRecordsWeatherGS1$SeedRate);
summary(farmRecordsWeatherGS1$NRatio);
summary(farmRecordsWeatherGS1$PRatio);
summary(farmRecordsWeatherGS1$KRatio);
summary(farmRecordsWeatherGS1$SRatio);
summary(farmRecordsWeatherGS1$InocType);
summary(farmRecordsWeatherGS1$Mat_Days);
summary(farmRecordsWeatherGS1$SwaDesc);
summary(farmRecordsWeatherGS1$Drought);
summary(farmRecordsWeatherGS1$Flood);
summary(farmRecordsWeatherGS1$Hail);
summary(farmRecordsWeatherGS1$Copper);
summary(farmRecordsWeatherGS1$Fungicide);
summary(farmRecordsWeatherGS1$Pests);
summary(farmRecordsWeatherGS1$TreatedSeed);
summary(farmRecordsWeatherGS1$Ave_Temp);
summary(farmRecordsWeatherGS1$Max_Temp);
summary(farmRecordsWeatherGS1$Min_Temp);
summary(farmRecordsWeatherGS1$RF_sum);
length(which(summary(farmRecordsWeatherGS1$Cultivar) > 0));
length(which(summary(farmRecordsWeatherGS1$Field) > 0));
length(which(summary(farmRecordsWeatherGS1$Year) > 0));
length(which(summary(farmRecordsWeatherGS1$prevCrop) > 0));
##below sd() doesn't work with NA values;
sd(farmRecordsWeatherGS1$Yield);
sd(farmRecordsWeatherGS1$SeedRate);
sd(farmRecordsWeatherGS1$NRatio);
sd(farmRecordsWeatherGS1$PRatio);
sd(farmRecordsWeatherGS1$KRatio);
sd(farmRecordsWeatherGS1$SRatio);
sd(farmRecordsWeatherGS1$InocType);
sd(farmRecordsWeatherGS1$Mat_Days);
sd(farmRecordsWeatherGS1$SwaDesc);
sd(farmRecordsWeatherGS1$Drought);
sd(farmRecordsWeatherGS1$Flood);
sd(farmRecordsWeatherGS1$Hail);
sd(farmRecordsWeatherGS1$Copper);
sd(farmRecordsWeatherGS1$Fungicide);
sd(farmRecordsWeatherGS1$Pests);
sd(farmRecordsWeatherGS1$TreatedSeed);
sd(farmRecordsWeatherGS1$Ave_Temp);
sd(farmRecordsWeatherGS1$Max_Temp);
sd(farmRecordsWeatherGS1$Min_Temp);
sd(farmRecordsWeatherGS1$RF_sum);

options(scipen=3); #suppresses scientific notation;

#model1 <- MCMCglmm(Yield ~ Ave_Temp + RF_sum, random=~Year + prevCrop, family="gaussian", data=farmRecordsWeatherGS1, saveX=TRUE, verbose=TRUE, burnin=1388889, nitt=13888889, thin=6250, pr=TRUE, prior=priorFarm3); #2000 samples;
#model1 <- MCMCglmm(Yield ~ Ave_Temp + RF_sum, random=~Year + prevCrop, family="gaussian", data=farmRecordsWeatherGS1, saveX=TRUE, verbose=TRUE, burnin=500, nitt=2000, thin=25, pr=TRUE, prior=priorFarm3); #2000 samples;
###model1 <- MCMCglmm(Yield ~ SeedRate + NRatio + PRatio + KRatio + SRatio + InocType + Mat_Days + Drought + Flood + Hail + Copper + Fungicide + Pests + TreatedSeed + Ave_Temp + Max_Temp + Min_Temp + RF_sum, random=~Cultivar + Field + Year + prevCrop, family="gaussian", data=farmRecordsWeatherGS1, saveX=TRUE, verbose=TRUE, burnin=500, nitt=2000, thin=25, pr=TRUE, prior=priorFarm1);


#save(model1, file=paste("farmRecordsRainfall_summary_barley_Model_1.RData", sep=""), compress="xz");
load(paste("farmRecordsRainfall_summary_barley_Model_1.RData", sep=""));
#load(paste("STR_model_published/farmRecordsRainfall_summary_barley_Model_1.RData", sep=""));

#https://stats.stackexchange.com/questions/307021/why-does-mcmcglmm-result-in-small-effective-sample-sizes-for-logistic-regression/310291
#autocorr.diag(model1$VCV); #use this to choose thinning interval, autocorrelation less than 0.01 ideally;

fe1 <- summary(model1);
#str(fe1$solutions);
fe1$solutions[,1] <- fe1$solutions[,1];
fe1$solutions[,2] <- fe1$solutions[,2];
fe1$solutions[,3] <- fe1$solutions[,3];

options(scipen=3); #suppresses scientific notation;
writeLines(capture.output(fe1), paste("farmRecordsRainfall_barley_Model_1.txt", sep=""));
writeLines(capture.output(summary(model1)), paste("farmRecordsRainfall_barley_Model_1.txt", sep=""));

### DO THIS FOR DIAGNOSTICS OF HOW CHAIN FIT with mcmcglmm;
#http://www.maths.bath.ac.uk/~jjf23/mixchange/split.html

#colnames(model1$Sol); #1:3 are explanatory variables;
EVs1 <- model1$Sol[, 1:3];
#intCultivars1 <- model1$Sol[, 25:55];
#intFields1 <- model1$Sol[, 56:93];
intProcYears1 <- model1$Sol[, 4:6];
intPrevCrops1 <- model1$Sol[, 7:9];

length(posterior.mode(intFields1)); #38;  #mean intercept;

### https://github.com/tmalsburg/MCMCglmm-intro
###

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farmRecordsRainfall_barley_model1_VCV_farm.pdf", sep=''), bg="transparent", width=10, height=8, pointsize=12, family="FreeSans");
par(mfrow=c(length(colnames(model1$VCV)),2), mar=c(2,2,1,0));
plot(model1$VCV, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farmRecordsRainfall_barley_model1_sol_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(EVs1)),2), mar=c(2,2,1,0));
plot(EVs1, auto.layout=F);
dev.off()

#pdf.options(encoding='CP1250');
#cairo_pdf(file=paste("farmRecordsRainfall_barley_model1_Cultivars_farm.pdf", sep=''), bg="transparent", width=11, height=30, pointsize=6, family="FreeSans");
#par(mfrow=c(length(colnames(intCultivars1)),2), mar=c(2,2,1,0));
#plot(intCultivars1, auto.layout=F);
#dev.off()

#pdf.options(encoding='CP1250');
#cairo_pdf(file=paste("farmRecordsRainfall_barley_model1_Fields_farm.pdf", sep=''), bg="transparent", width=11, height=30, pointsize=6, family="FreeSans");
#par(mfrow=c(length(colnames(intFields1)),2), mar=c(2,2,1,0));
#plot(intFields1, auto.layout=F);
#dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farmRecordsRainfall_barley_model1_ProcYear_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intProcYears1)),2), mar=c(2,2,1,0));
plot(intProcYears1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farmRecordsRainfall_barley_model1_PrevCrop_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intPrevCrops1)),2), mar=c(2,2,1,0));
plot(intPrevCrops1, auto.layout=F);
dev.off()

colnames(model1$VCV);
#ICC_Cultivar <- model1$VCV[, 1]/(rowSums(model1$VCV)); 
#ICC_Field <- model1$VCV[, 2]/(rowSums(model1$VCV)); 
ICC_ProcYear <- model1$VCV[, 1]/(rowSums(model1$VCV)); 
ICC_PrevCrop <- model1$VCV[, 2]/(rowSums(model1$VCV)); 
#dft1.1 <- cbind(ICC = posterior.mode(ICC_Cultivar), CI = HPDinterval(ICC_Cultivar));
#dft2.1 <- cbind(ICC = posterior.mode(ICC_Field), CI = HPDinterval(ICC_Field));
dft3.1 <- cbind(ICC = posterior.mode(ICC_ProcYear), CI = HPDinterval(ICC_ProcYear));
dft4.1 <- cbind(ICC = posterior.mode(ICC_PrevCrop), CI = HPDinterval(ICC_PrevCrop));
dft1 <- rbind(dft3.1, dft4.1);
row.names(dft1) <- c("ICC_ProcYear", "ICC_PrevCrop");
options(scipen=3); #suppresses scientific notation;
writeLines(capture.output(dft1), paste("farmRecordsRainfall_ICC_barley_Model_1.txt", sep=""));

#---;
####-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----####;
####-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----####;


####-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----####;
####-----@-----@-----@-----@--- Does for canary ---@-----@-----@-----@-----@-----####;
farmRecordsWeatherGS1 <- farmRecordsWeatherGS[which(farmRecordsWeatherGS$Crop == "Canary"), ];
str(farmRecordsWeatherGS1);

##remove field records with missing values for Yield (i.e. summerfallow);
farmRecordsWeatherGS1 <- farmRecordsWeatherGS1[which(!is.na(farmRecordsWeatherGS1$RF_sum)), ];
                 
length(farmRecordsWeatherGS1$SeedRate);
summary(farmRecordsWeatherGS1$Yield);
summary(farmRecordsWeatherGS1$SeedRate);
summary(farmRecordsWeatherGS1$NRatio);
summary(farmRecordsWeatherGS1$PRatio);
summary(farmRecordsWeatherGS1$KRatio);
summary(farmRecordsWeatherGS1$SRatio);
summary(farmRecordsWeatherGS1$InocType);
summary(farmRecordsWeatherGS1$Mat_Days);
summary(farmRecordsWeatherGS1$SwaDesc);
summary(farmRecordsWeatherGS1$Drought);
summary(farmRecordsWeatherGS1$Flood);
summary(farmRecordsWeatherGS1$Hail);
summary(farmRecordsWeatherGS1$Copper);
summary(farmRecordsWeatherGS1$Fungicide);
summary(farmRecordsWeatherGS1$Pests);
summary(farmRecordsWeatherGS1$TreatedSeed);
summary(farmRecordsWeatherGS1$Ave_Temp);
summary(farmRecordsWeatherGS1$Max_Temp);
summary(farmRecordsWeatherGS1$Min_Temp);
summary(farmRecordsWeatherGS1$RF_sum);
length(which(summary(farmRecordsWeatherGS1$Cultivar) > 0));
length(which(summary(farmRecordsWeatherGS1$Field) > 0));
length(which(summary(farmRecordsWeatherGS1$Year) > 0));
length(which(summary(farmRecordsWeatherGS1$prevCrop) > 0));
##below sd() doesn't work with NA values;
sd(farmRecordsWeatherGS1$Yield);
sd(farmRecordsWeatherGS1$SeedRate);
sd(farmRecordsWeatherGS1$NRatio);
sd(farmRecordsWeatherGS1$PRatio);
sd(farmRecordsWeatherGS1$KRatio);
sd(farmRecordsWeatherGS1$SRatio);
sd(farmRecordsWeatherGS1$InocType);
sd(farmRecordsWeatherGS1$Mat_Days);
sd(farmRecordsWeatherGS1$SwaDesc);
sd(farmRecordsWeatherGS1$Drought);
sd(farmRecordsWeatherGS1$Flood);
sd(farmRecordsWeatherGS1$Hail);
sd(farmRecordsWeatherGS1$Copper);
sd(farmRecordsWeatherGS1$Fungicide);
sd(farmRecordsWeatherGS1$Pests);
sd(farmRecordsWeatherGS1$TreatedSeed);
sd(farmRecordsWeatherGS1$Ave_Temp);
sd(farmRecordsWeatherGS1$Max_Temp);
sd(farmRecordsWeatherGS1$Min_Temp);
sd(farmRecordsWeatherGS1$RF_sum);

options(scipen=3); #suppresses scientific notation;

#model1 <- MCMCglmm(Yield ~ SeedRate + NRatio + PRatio + KRatio + Mat_Days + SwaDesc + Fungicide + Pests + Ave_Temp + RF_sum + RF_sum*Ave_Temp, random=~Cultivar + Field + Year + prevCrop, family="gaussian", data=farmRecordsWeatherGS1, saveX=TRUE, verbose=TRUE, burnin=1388889, nitt=13888889, thin=6250, pr=TRUE, prior=priorFarm1); #2000 samples;
#model1 <- MCMCglmm(Yield ~ SeedRate + NRatio + PRatio + KRatio + Mat_Days + SwaDesc + Fungicide + Pests + Ave_Temp + RF_sum + Ave_Temp*RF_sum, random=~Cultivar + Field + Year + prevCrop, family="gaussian", data=farmRecordsWeatherGS1, saveX=TRUE, verbose=TRUE, burnin=500, nitt=2000, thin=25, pr=TRUE, prior=priorFarm1); #2000 samples;
#model1 <- MCMCglmm(Yield ~ SeedRate + NRatio + PRatio + KRatio + Mat_Days + SwaDesc + Fungicide + Pests + Ave_Temp + RF_sum + Ave_Temp*RF_sum, random=~Cultivar + Field + Year + prevCrop, family="gaussian", data=farmRecordsWeatherGS1, saveX=TRUE, verbose=TRUE, burnin=1388889, nitt=13888889, thin=6250, pr=TRUE, prior=priorFarm1); #2000 samples;

#model1 <- MCMCglmm(Yield ~ SeedRate + NRatio + PRatio + KRatio + Mat_Days + Drought + Flood + Fungicide + Pests + TreatedSeed + Ave_Temp + Max_Temp + Min_Temp + RF_sum, random=~Cultivar + Field + Year + prevCrop, family="gaussian", data=farmRecordsWeatherGS1, saveX=TRUE, verbose=TRUE, burnin=1388889, nitt=13888889, thin=6250, pr=TRUE, prior=priorFarm1); #2000 samples;
#model1 <- MCMCglmm(Yield ~ SeedRate + NRatio + PRatio + KRatio + Mat_Days + Drought + Flood + Fungicide + Pests + TreatedSeed + Ave_Temp + Max_Temp + Min_Temp + RF_sum, random=~Cultivar + Field + Year + prevCrop, family="gaussian", data=farmRecordsWeatherGS1, saveX=TRUE, verbose=TRUE, burnin=500, nitt=2000, thin=25, pr=TRUE, prior=priorFarm1); #2000 samples;
###model1 <- MCMCglmm(Yield ~ SeedRate + NRatio + PRatio + KRatio + SRatio + InocType + Mat_Days + Drought + Flood + Hail + Copper + Fungicide + Pests + TreatedSeed + Ave_Temp + Max_Temp + Min_Temp + RF_sum, random=~Cultivar + Field + Year + prevCrop, family="gaussian", data=farmRecordsWeatherGS1, saveX=TRUE, verbose=TRUE, burnin=500, nitt=2000, thin=25, pr=TRUE, prior=priorFarm1);


#save(model1, file=paste("farmRecordsRainfall_summary_canary_Model_1.RData", sep=""), compress="xz");
load(paste("farmRecordsRainfall_summary_canary_Model_1.RData", sep=""));
#load(paste("STR_model_published/farmRecordsRainfall_summary_canary_Model_1.RData", sep=""));

#https://stats.stackexchange.com/questions/307021/why-does-mcmcglmm-result-in-small-effective-sample-sizes-for-logistic-regression/310291
#autocorr.diag(model1$VCV); #use this to choose thinning interval, autocorrelation less than 0.01 ideally;

fe1 <- summary(model1);
#str(fe1$solutions);
fe1$solutions[,1] <- fe1$solutions[,1];
fe1$solutions[,2] <- fe1$solutions[,2];
fe1$solutions[,3] <- fe1$solutions[,3];

options(scipen=3); #suppresses scientific notation;
writeLines(capture.output(fe1), paste("farmRecordsRainfall_canary_Model_1.txt", sep=""));
writeLines(capture.output(summary(model1)), paste("farmRecordsRainfall_canary_Model_1.txt", sep=""));

### DO THIS FOR DIAGNOSTICS OF HOW CHAIN FIT with mcmcglmm;
#http://www.maths.bath.ac.uk/~jjf23/mixchange/split.html

#colnames(model1$Sol); #1:12 are explanatory variables;
EVs1 <- model1$Sol[, 1:12];
intCultivars1 <- model1$Sol[, 13:15];
intFields1 <- model1$Sol[, 16:45];
intProcYears1 <- model1$Sol[, 46:68];
intPrevCrops1 <- model1$Sol[, 69:73];

length(posterior.mode(intFields1)); #38;  #mean intercept;

### https://github.com/tmalsburg/MCMCglmm-intro
###

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farmRecordsRainfall_canary_model1_VCV_farm.pdf", sep=''), bg="transparent", width=10, height=8, pointsize=12, family="FreeSans");
par(mfrow=c(length(colnames(model1$VCV)),2), mar=c(2,2,1,0));
plot(model1$VCV, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farmRecordsRainfall_canary_model1_sol_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(EVs1)),2), mar=c(2,2,1,0));
plot(EVs1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farmRecordsRainfall_canary_model1_Cultivars_farm.pdf", sep=''), bg="transparent", width=11, height=30, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intCultivars1)),2), mar=c(2,2,1,0));
plot(intCultivars1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farmRecordsRainfall_canary_model1_Fields_farm.pdf", sep=''), bg="transparent", width=11, height=30, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intFields1)),2), mar=c(2,2,1,0));
plot(intFields1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farmRecordsRainfall_canary_model1_ProcYear_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intProcYears1)),2), mar=c(2,2,1,0));
plot(intProcYears1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farmRecordsRainfall_canary_model1_PrevCrop_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intPrevCrops1)),2), mar=c(2,2,1,0));
plot(intPrevCrops1, auto.layout=F);
dev.off()

colnames(model1$VCV);
ICC_Cultivar <- model1$VCV[, 1]/(rowSums(model1$VCV)); 
ICC_Field <- model1$VCV[, 2]/(rowSums(model1$VCV)); 
ICC_ProcYear <- model1$VCV[, 3]/(rowSums(model1$VCV)); 
ICC_PrevCrop <- model1$VCV[, 4]/(rowSums(model1$VCV)); 
dft1.1 <- cbind(ICC = posterior.mode(ICC_Cultivar), CI = HPDinterval(ICC_Cultivar));
dft2.1 <- cbind(ICC = posterior.mode(ICC_Field), CI = HPDinterval(ICC_Field));
dft3.1 <- cbind(ICC = posterior.mode(ICC_ProcYear), CI = HPDinterval(ICC_ProcYear));
dft4.1 <- cbind(ICC = posterior.mode(ICC_PrevCrop), CI = HPDinterval(ICC_PrevCrop));
dft1 <- rbind(dft1.1, dft2.1, dft3.1, dft4.1);
row.names(dft1) <- c("ICC_Cultivar", "ICC_Field", "ICC_ProcYear", "ICC_PrevCrop");
options(scipen=3); #suppresses scientific notation;
writeLines(capture.output(dft1), paste("farmRecordsRainfall_ICC_canary_Model_1.txt", sep=""));

#---;
####-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----####;
####-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----####;


####-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----####;
####-----@-----@-----@-----@--- Does for canola ---@-----@-----@-----@-----@-----####;
farmRecordsWeatherGS1 <- farmRecordsWeatherGS[which(farmRecordsWeatherGS$Crop == "Canola"), ];
str(farmRecordsWeatherGS1);

##remove field records with missing values for Yield (i.e. summerfallow);
farmRecordsWeatherGS1 <- farmRecordsWeatherGS1[which(!is.na(farmRecordsWeatherGS1$RF_sum)), ];
                 
length(farmRecordsWeatherGS1$SeedRate);
summary(farmRecordsWeatherGS1$Yield);
summary(farmRecordsWeatherGS1$SeedRate);
summary(farmRecordsWeatherGS1$NRatio);
summary(farmRecordsWeatherGS1$PRatio);
summary(farmRecordsWeatherGS1$KRatio);
summary(farmRecordsWeatherGS1$SRatio);
summary(farmRecordsWeatherGS1$InocType);
summary(farmRecordsWeatherGS1$Mat_Days);
summary(farmRecordsWeatherGS1$SwaDesc);
summary(farmRecordsWeatherGS1$Drought);
summary(farmRecordsWeatherGS1$Flood);
summary(farmRecordsWeatherGS1$Hail);
summary(farmRecordsWeatherGS1$Copper);
summary(farmRecordsWeatherGS1$Fungicide);
summary(farmRecordsWeatherGS1$Pests);
summary(farmRecordsWeatherGS1$TreatedSeed);
summary(farmRecordsWeatherGS1$Ave_Temp);
summary(farmRecordsWeatherGS1$Max_Temp);
summary(farmRecordsWeatherGS1$Min_Temp);
summary(farmRecordsWeatherGS1$RF_sum);
length(which(summary(farmRecordsWeatherGS1$Cultivar) > 0));
length(which(summary(farmRecordsWeatherGS1$Field) > 0));
length(which(summary(farmRecordsWeatherGS1$Year) > 0));
length(which(summary(farmRecordsWeatherGS1$prevCrop) > 0));
##below sd() doesn't work with NA values;
sd(farmRecordsWeatherGS1$Yield);
sd(farmRecordsWeatherGS1$SeedRate);
sd(farmRecordsWeatherGS1$NRatio);
sd(farmRecordsWeatherGS1$PRatio);
sd(farmRecordsWeatherGS1$KRatio);
sd(farmRecordsWeatherGS1$SRatio);
sd(farmRecordsWeatherGS1$InocType);
sd(farmRecordsWeatherGS1$Mat_Days);
sd(farmRecordsWeatherGS1$SwaDesc);
sd(farmRecordsWeatherGS1$Drought);
sd(farmRecordsWeatherGS1$Flood);
sd(farmRecordsWeatherGS1$Hail);
sd(farmRecordsWeatherGS1$Copper);
sd(farmRecordsWeatherGS1$Fungicide);
sd(farmRecordsWeatherGS1$Pests);
sd(farmRecordsWeatherGS1$TreatedSeed);
sd(farmRecordsWeatherGS1$Ave_Temp);
sd(farmRecordsWeatherGS1$Max_Temp);
sd(farmRecordsWeatherGS1$Min_Temp);
sd(farmRecordsWeatherGS1$RF_sum);

options(scipen=3); #suppresses scientific notation;

#model1 <- MCMCglmm(Yield ~ SeedRate + NRatio + PRatio + SRatio + Mat_Days + SwaDesc + Fungicide + Pests + Ave_Temp + RF_sum + RF_sum*Ave_Temp, family="gaussian", data=farmRecordsWeatherGS1, saveX=TRUE, verbose=TRUE, burnin=1388889, nitt=13888889, thin=6250, pr=TRUE, prior=priorFarm5); #2000 samples;
#model1 <- MCMCglmm(Yield ~ SeedRate + NRatio + PRatio + SRatio + Mat_Days + SwaDesc + Fungicide + Pests + Ave_Temp + RF_sum + Ave_Temp*RF_sum, family="gaussian", data=farmRecordsWeatherGS1, saveX=TRUE, verbose=TRUE, burnin=500, nitt=2000, thin=25, pr=TRUE, prior=priorFarm5); #2000 samples;
#model1 <- MCMCglmm(Yield ~ SeedRate + NRatio + PRatio + SRatio + Mat_Days + SwaDesc + Fungicide + Pests + Ave_Temp + RF_sum + Ave_Temp*RF_sum, random=~Cultivar + Year + prevCrop, family="gaussian", data=farmRecordsWeatherGS1, saveX=TRUE, verbose=TRUE, burnin=1388889, nitt=13888889, thin=6250, pr=TRUE, prior=priorFarm2); #2000 samples;

#model1 <- MCMCglmm(Yield ~ SeedRate + NRatio + PRatio + SRatio + Mat_Days + Fungicide + Pests + Ave_Temp + Max_Temp + Min_Temp + RF_sum, random=~Cultivar + Year + prevCrop, family="gaussian", data=farmRecordsWeatherGS1, saveX=TRUE, verbose=TRUE, burnin=1388889, nitt=13888889, thin=6250, pr=TRUE, prior=priorFarm2); #2000 samples;
#model1 <- MCMCglmm(Yield ~ SeedRate + NRatio + PRatio + SRatio + Mat_Days + Fungicide + Pests + Ave_Temp + Max_Temp + Min_Temp + RF_sum, random=~Cultivar + Year + prevCrop, family="gaussian", data=farmRecordsWeatherGS1, saveX=TRUE, verbose=TRUE, burnin=500, nitt=2000, thin=25, pr=TRUE, prior=priorFarm2); #2000 samples;
###model1 <- MCMCglmm(Yield ~ SeedRate + NRatio + PRatio + KRatio + SRatio + InocType + Mat_Days + Drought + Flood + Hail + Copper + Fungicide + Pests + TreatedSeed + Ave_Temp + Max_Temp + Min_Temp + RF_sum, random=~Cultivar + Field + Year + prevCrop, family="gaussian", data=farmRecordsWeatherGS1, saveX=TRUE, verbose=TRUE, burnin=500, nitt=2000, thin=25, pr=TRUE, prior=priorFarm1);


#save(model1, file=paste("farmRecordsRainfall_summary_canola_Model_1.RData", sep=""), compress="xz");
load(paste("farmRecordsRainfall_summary_canola_Model_1.RData", sep=""));
#load(paste("STR_model_published/farmRecordsRainfall_summary_canola_Model_1.RData", sep=""));

#https://stats.stackexchange.com/questions/307021/why-does-mcmcglmm-result-in-small-effective-sample-sizes-for-logistic-regression/310291
#autocorr.diag(model1$VCV); #use this to choose thinning interval, autocorrelation less than 0.01 ideally;

fe1 <- summary(model1);
#str(fe1$solutions);
fe1$solutions[,1] <- fe1$solutions[,1];
fe1$solutions[,2] <- fe1$solutions[,2];
fe1$solutions[,3] <- fe1$solutions[,3];

options(scipen=3); #suppresses scientific notation;
writeLines(capture.output(fe1), paste("farmRecordsRainfall_canola_Model_1.txt", sep=""));
writeLines(capture.output(summary(model1)), paste("farmRecordsRainfall_canola_Model_1.txt", sep=""));

### DO THIS FOR DIAGNOSTICS OF HOW CHAIN FIT with mcmcglmm;
#http://www.maths.bath.ac.uk/~jjf23/mixchange/split.html

#colnames(model1$Sol); #1:12 are explanatory variables;
EVs1 <- model1$Sol[, 1:12];
#intCultivars1 <- model1$Sol[, 13:17];
#intFields1 <- model1$Sol[, 19:48];
#intProcYears1 <- model1$Sol[, 13:19];
#intPrevCrops1 <- model1$Sol[, 25:29];

length(posterior.mode(intFields1)); #38;  #mean intercept;

### https://github.com/tmalsburg/MCMCglmm-intro
###

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farmRecordsRainfall_canola_model1_VCV_farm.pdf", sep=''), bg="transparent", width=10, height=8, pointsize=12, family="FreeSans");
par(mfrow=c(length(colnames(model1$VCV)),2), mar=c(2,2,1,0));
plot(model1$VCV, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farmRecordsRainfall_canola_model1_sol_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(EVs1)),2), mar=c(2,2,1,0));
plot(EVs1, auto.layout=F);
dev.off()

#pdf.options(encoding='CP1250');
#cairo_pdf(file=paste("farmRecordsRainfall_canola_model1_Cultivars_farm.pdf", sep=''), bg="transparent", width=11, height=30, pointsize=6, family="FreeSans");
#par(mfrow=c(length(colnames(intCultivars1)),2), mar=c(2,2,1,0));
#plot(intCultivars1, auto.layout=F);
#dev.off()

#pdf.options(encoding='CP1250');
#cairo_pdf(file=paste("farmRecordsRainfall_canola_model1_Fields_farm.pdf", sep=''), bg="transparent", width=11, height=30, pointsize=6, family="FreeSans");
#par(mfrow=c(length(colnames(intFields1)),2), mar=c(2,2,1,0));
#plot(intFields1, auto.layout=F);
#dev.off()

#pdf.options(encoding='CP1250');
#cairo_pdf(file=paste("farmRecordsRainfall_canola_model1_ProcYear_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
#par(mfrow=c(length(colnames(intProcYears1)),2), mar=c(2,2,1,0));
#plot(intProcYears1, auto.layout=F);
#dev.off()

#pdf.options(encoding='CP1250');
#cairo_pdf(file=paste("farmRecordsRainfall_canola_model1_PrevCrop_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
#par(mfrow=c(length(colnames(intPrevCrops1)),2), mar=c(2,2,1,0));
#plot(intPrevCrops1, auto.layout=F);
#dev.off()

#colnames(model1$VCV);
#ICC_Cultivar <- model1$VCV[, 1]/(rowSums(model1$VCV)); 
#ICC_Field <- model1$VCV[, 2]/(rowSums(model1$VCV)); 
#ICC_ProcYear <- model1$VCV[, 1]/(rowSums(model1$VCV)); 
#ICC_PrevCrop <- model1$VCV[, 3]/(rowSums(model1$VCV)); 
#dft1.1 <- cbind(ICC = posterior.mode(ICC_Cultivar), CI = HPDinterval(ICC_Cultivar));
#dft2.1 <- cbind(ICC = posterior.mode(ICC_Field), CI = HPDinterval(ICC_Field));
#dft3.1 <- cbind(ICC = posterior.mode(ICC_ProcYear), CI = HPDinterval(ICC_ProcYear));
#dft4.1 <- cbind(ICC = posterior.mode(ICC_PrevCrop), CI = HPDinterval(ICC_PrevCrop));
#dft1 <- rbind(dft3.1);
#row.names(dft1) <- c("ICC_ProcYear");
#options(scipen=3); #suppresses scientific notation;
#writeLines(capture.output(dft1), paste("farmRecordsRainfall_ICC_canola_Model_1.txt", sep=""));

#---;
####-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----####;
####-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----####;


#### IGNORE BELOW SECTION, NOT USED IN TABLE 3 BUT IS IN TABLE 2 ####;
####-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----####;
####-----@-----@-----@-----@--- Does for chickpea ---@-----@-----@-----@-----@-----####;
farmRecordsWeatherGS1 <- farmRecordsWeatherGS[which(farmRecordsWeatherGS$Crop == "Chickpea"), ];
str(farmRecordsWeatherGS1);

##remove field records with missing values for Yield (i.e. summerfallow);
farmRecordsWeatherGS1 <- farmRecordsWeatherGS1[which(!is.na(farmRecordsWeatherGS1$RF_sum)), ];
                 
length(farmRecordsWeatherGS1$SeedRate);
summary(farmRecordsWeatherGS1$Yield);
summary(farmRecordsWeatherGS1$SeedRate);
summary(farmRecordsWeatherGS1$NRatio);
summary(farmRecordsWeatherGS1$PRatio);
summary(farmRecordsWeatherGS1$KRatio);
summary(farmRecordsWeatherGS1$SRatio);
summary(farmRecordsWeatherGS1$InocType);
summary(farmRecordsWeatherGS1$Mat_Days);
summary(farmRecordsWeatherGS1$SwaDesc);
summary(farmRecordsWeatherGS1$Drought);
summary(farmRecordsWeatherGS1$Flood);
summary(farmRecordsWeatherGS1$Hail);
summary(farmRecordsWeatherGS1$Copper);
summary(farmRecordsWeatherGS1$Fungicide);
summary(farmRecordsWeatherGS1$Pests);
summary(farmRecordsWeatherGS1$TreatedSeed);
summary(farmRecordsWeatherGS1$Ave_Temp);
summary(farmRecordsWeatherGS1$Max_Temp);
summary(farmRecordsWeatherGS1$Min_Temp);
summary(farmRecordsWeatherGS1$RF_sum);
length(which(summary(farmRecordsWeatherGS1$Cultivar) > 0));
length(which(summary(farmRecordsWeatherGS1$Field) > 0));
length(which(summary(farmRecordsWeatherGS1$Year) > 0));
length(which(summary(farmRecordsWeatherGS1$prevCrop) > 0));
##below sd() doesn't work with NA values;
sd(farmRecordsWeatherGS1$Yield);
sd(farmRecordsWeatherGS1$SeedRate);
sd(farmRecordsWeatherGS1$NRatio);
sd(farmRecordsWeatherGS1$PRatio);
sd(farmRecordsWeatherGS1$KRatio);
sd(farmRecordsWeatherGS1$SRatio);
sd(farmRecordsWeatherGS1$InocType);
sd(farmRecordsWeatherGS1$Mat_Days);
sd(farmRecordsWeatherGS1$SwaDesc);
sd(farmRecordsWeatherGS1$Drought);
sd(farmRecordsWeatherGS1$Flood);
sd(farmRecordsWeatherGS1$Hail);
sd(farmRecordsWeatherGS1$Copper);
sd(farmRecordsWeatherGS1$Fungicide);
sd(farmRecordsWeatherGS1$Pests);
sd(farmRecordsWeatherGS1$TreatedSeed);
sd(farmRecordsWeatherGS1$Ave_Temp);
sd(farmRecordsWeatherGS1$Max_Temp);
sd(farmRecordsWeatherGS1$Min_Temp);
sd(farmRecordsWeatherGS1$RF_sum);

options(scipen=3); #suppresses scientific notation;

#model1 <- MCMCglmm(Yield ~ Ave_Temp + RF_sum, random=~Year + prevCrop, family="gaussian", data=farmRecordsWeatherGS1, saveX=TRUE, verbose=TRUE, burnin=1388889, nitt=13888889, thin=6250, pr=TRUE, prior=priorFarm3); #2000 samples;
#model1 <- MCMCglmm(Yield ~ Ave_Temp + RF_sum, random=~Year + prevCrop, family="gaussian", data=farmRecordsWeatherGS1, saveX=TRUE, verbose=TRUE, burnin=500, nitt=2000, thin=25, pr=TRUE, prior=priorFarm3); #2000 samples;
###model1 <- MCMCglmm(Yield ~ SeedRate + NRatio + PRatio + KRatio + SRatio + InocType + Mat_Days + Drought + Flood + Hail + Copper + Fungicide + Pests + TreatedSeed + Ave_Temp + Max_Temp + Min_Temp + RF_sum, random=~Cultivar + Field + Year + prevCrop, family="gaussian", data=farmRecordsWeatherGS1, saveX=TRUE, verbose=TRUE, burnin=500, nitt=2000, thin=25, pr=TRUE, prior=priorFarm1);


#save(model1, file=paste("farmRecordsRainfall_summary_chickpea_Model_1.RData", sep=""), compress="xz");
load(paste("farmRecordsRainfall_summary_chickpea_Model_1.RData", sep=""));
#load(paste("STR_model_published/farmRecordsRainfall_summary_chickpea_Model_1.RData", sep=""));

#https://stats.stackexchange.com/questions/307021/why-does-mcmcglmm-result-in-small-effective-sample-sizes-for-logistic-regression/310291
#autocorr.diag(model1$VCV); #use this to choose thinning interval, autocorrelation less than 0.01 ideally;

fe1 <- summary(model1);
#str(fe1$solutions);
fe1$solutions[,1] <- fe1$solutions[,1];
fe1$solutions[,2] <- fe1$solutions[,2];
fe1$solutions[,3] <- fe1$solutions[,3];

options(scipen=3); #suppresses scientific notation;
writeLines(capture.output(fe1), paste("farmRecordsRainfall_chickpea_Model_1.txt", sep=""));
writeLines(capture.output(summary(model1)), paste("farmRecordsRainfall_chickpea_Model_1.txt", sep=""));

### DO THIS FOR DIAGNOSTICS OF HOW CHAIN FIT with mcmcglmm;
#http://www.maths.bath.ac.uk/~jjf23/mixchange/split.html

#colnames(model1$Sol); #1:3 are explanatory variables;
EVs1 <- model1$Sol[, 1:3];
#intCultivars1 <- model1$Sol[, 13:18];
#intFields1 <- model1$Sol[, 19:48];
intProcYears1 <- model1$Sol[, 4:9];
intPrevCrops1 <- model1$Sol[, 10:12];

length(posterior.mode(intFields1)); #38;  #mean intercept;

### https://github.com/tmalsburg/MCMCglmm-intro
###

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farmRecordsRainfall_chickpea_model1_VCV_farm.pdf", sep=''), bg="transparent", width=10, height=8, pointsize=12, family="FreeSans");
par(mfrow=c(length(colnames(model1$VCV)),2), mar=c(2,2,1,0));
plot(model1$VCV, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farmRecordsRainfall_chickpea_model1_sol_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(EVs1)),2), mar=c(2,2,1,0));
plot(EVs1, auto.layout=F);
dev.off()

#pdf.options(encoding='CP1250');
#cairo_pdf(file=paste("farmRecordsRainfall_chickpea_model1_Cultivars_farm.pdf", sep=''), bg="transparent", width=11, height=30, pointsize=6, family="FreeSans");
#par(mfrow=c(length(colnames(intCultivars1)),2), mar=c(2,2,1,0));
#plot(intCultivars1, auto.layout=F);
#dev.off()

#pdf.options(encoding='CP1250');
#cairo_pdf(file=paste("farmRecordsRainfall_chickpea_model1_Fields_farm.pdf", sep=''), bg="transparent", width=11, height=30, pointsize=6, family="FreeSans");
#par(mfrow=c(length(colnames(intFields1)),2), mar=c(2,2,1,0));
#plot(intFields1, auto.layout=F);
#dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farmRecordsRainfall_chickpea_model1_ProcYear_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intProcYears1)),2), mar=c(2,2,1,0));
plot(intProcYears1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farmRecordsRainfall_chickpea_model1_PrevCrop_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intPrevCrops1)),2), mar=c(2,2,1,0));
plot(intPrevCrops1, auto.layout=F);
dev.off()

colnames(model1$VCV);
#ICC_Cultivar <- model1$VCV[, 1]/(rowSums(model1$VCV)); 
#ICC_Field <- model1$VCV[, 2]/(rowSums(model1$VCV)); 
ICC_ProcYear <- model1$VCV[, 1]/(rowSums(model1$VCV)); 
ICC_PrevCrop <- model1$VCV[, 2]/(rowSums(model1$VCV)); 
#dft1.1 <- cbind(ICC = posterior.mode(ICC_Cultivar), CI = HPDinterval(ICC_Cultivar));
#dft2.1 <- cbind(ICC = posterior.mode(ICC_Field), CI = HPDinterval(ICC_Field));
dft3.1 <- cbind(ICC = posterior.mode(ICC_ProcYear), CI = HPDinterval(ICC_ProcYear));
dft4.1 <- cbind(ICC = posterior.mode(ICC_PrevCrop), CI = HPDinterval(ICC_PrevCrop));
dft1 <- rbind(dft3.1, dft4.1);
row.names(dft1) <- c("ICC_ProcYear", "ICC_PrevCrop");
options(scipen=3); #suppresses scientific notation;
writeLines(capture.output(dft1), paste("farmRecordsRainfall_ICC_chickpea_Model_1.txt", sep=""));

#---;
####-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----####;
####-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----####;


#### IGNORE BELOW SECTION, NOT USED IN TABLE 3 BUT IS IN TABLE 2 ####;
####-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----####;
####-----@-----@-----@-----@--- Does for field pea ---@-----@-----@-----@-----@-----####;
farmRecordsWeatherGS1 <- farmRecordsWeatherGS[which(farmRecordsWeatherGS$Crop == "Field Pea"), ];
str(farmRecordsWeatherGS1);

##remove field records with missing values for Yield (i.e. summerfallow);
farmRecordsWeatherGS1 <- farmRecordsWeatherGS1[which(!is.na(farmRecordsWeatherGS1$RF_sum)), ];
                 
length(farmRecordsWeatherGS1$SeedRate);
summary(farmRecordsWeatherGS1$Yield);
summary(farmRecordsWeatherGS1$SeedRate);
summary(farmRecordsWeatherGS1$NRatio);
summary(farmRecordsWeatherGS1$PRatio);
summary(farmRecordsWeatherGS1$KRatio);
summary(farmRecordsWeatherGS1$SRatio);
summary(farmRecordsWeatherGS1$InocType);
summary(farmRecordsWeatherGS1$Mat_Days);
summary(farmRecordsWeatherGS1$SwaDesc);
summary(farmRecordsWeatherGS1$Drought);
summary(farmRecordsWeatherGS1$Flood);
summary(farmRecordsWeatherGS1$Hail);
summary(farmRecordsWeatherGS1$Copper);
summary(farmRecordsWeatherGS1$Fungicide);
summary(farmRecordsWeatherGS1$Pests);
summary(farmRecordsWeatherGS1$TreatedSeed);
summary(farmRecordsWeatherGS1$Ave_Temp);
summary(farmRecordsWeatherGS1$Max_Temp);
summary(farmRecordsWeatherGS1$Min_Temp);
summary(farmRecordsWeatherGS1$RF_sum);
length(which(summary(farmRecordsWeatherGS1$Cultivar) > 0));
length(which(summary(farmRecordsWeatherGS1$Field) > 0));
length(which(summary(farmRecordsWeatherGS1$Year) > 0));
length(which(summary(farmRecordsWeatherGS1$prevCrop) > 0));
##below sd() doesn't work with NA values;
sd(farmRecordsWeatherGS1$Yield);
sd(farmRecordsWeatherGS1$SeedRate);
sd(farmRecordsWeatherGS1$NRatio);
sd(farmRecordsWeatherGS1$PRatio);
sd(farmRecordsWeatherGS1$KRatio);
sd(farmRecordsWeatherGS1$SRatio);
sd(farmRecordsWeatherGS1$InocType);
sd(farmRecordsWeatherGS1$Mat_Days);
sd(farmRecordsWeatherGS1$SwaDesc);
sd(farmRecordsWeatherGS1$Drought);
sd(farmRecordsWeatherGS1$Flood);
sd(farmRecordsWeatherGS1$Hail);
sd(farmRecordsWeatherGS1$Copper);
sd(farmRecordsWeatherGS1$Fungicide);
sd(farmRecordsWeatherGS1$Pests);
sd(farmRecordsWeatherGS1$TreatedSeed);
sd(farmRecordsWeatherGS1$Ave_Temp);
sd(farmRecordsWeatherGS1$Max_Temp);
sd(farmRecordsWeatherGS1$Min_Temp);
sd(farmRecordsWeatherGS1$RF_sum);

options(scipen=3); #suppresses scientific notation;

#model1 <- MCMCglmm(Yield ~ Ave_Temp + RF_sum, random=~Year + prevCrop, family="gaussian", data=farmRecordsWeatherGS1, saveX=TRUE, verbose=TRUE, burnin=1388889, nitt=13888889, thin=6250, pr=TRUE, prior=priorFarm3); #2000 samples;
#model1 <- MCMCglmm(Yield ~ Ave_Temp + RF_sum, random=~Year + prevCrop, family="gaussian", data=farmRecordsWeatherGS1, saveX=TRUE, verbose=TRUE, burnin=500, nitt=2000, thin=25, pr=TRUE, prior=priorFarm3); #2000 samples;
###model1 <- MCMCglmm(Yield ~ SeedRate + NRatio + PRatio + KRatio + SRatio + InocType + Mat_Days + Drought + Flood + Hail + Copper + Fungicide + Pests + TreatedSeed + Ave_Temp + Max_Temp + Min_Temp + RF_sum, random=~Cultivar + Field + Year + prevCrop, family="gaussian", data=farmRecordsWeatherGS1, saveX=TRUE, verbose=TRUE, burnin=500, nitt=2000, thin=25, pr=TRUE, prior=priorFarm1);


#save(model1, file=paste("farmRecordsRainfall_summary_fieldpea_Model_1.RData", sep=""), compress="xz");
load(paste("farmRecordsRainfall_summary_fieldpea_Model_1.RData", sep=""));
#load(paste("STR_model_published/farmRecordsRainfall_summary_fieldpea_Model_1.RData", sep=""));

#https://stats.stackexchange.com/questions/307021/why-does-mcmcglmm-result-in-small-effective-sample-sizes-for-logistic-regression/310291
#autocorr.diag(model1$VCV); #use this to choose thinning interval, autocorrelation less than 0.01 ideally;

fe1 <- summary(model1);
#str(fe1$solutions);
fe1$solutions[,1] <- fe1$solutions[,1];
fe1$solutions[,2] <- fe1$solutions[,2];
fe1$solutions[,3] <- fe1$solutions[,3];

options(scipen=3); #suppresses scientific notation;
writeLines(capture.output(fe1), paste("farmRecordsRainfall_fieldpea_Model_1.txt", sep=""));
writeLines(capture.output(summary(model1)), paste("farmRecordsRainfall_fieldpea_Model_1.txt", sep=""));

### DO THIS FOR DIAGNOSTICS OF HOW CHAIN FIT with mcmcglmm;
#http://www.maths.bath.ac.uk/~jjf23/mixchange/split.html

#colnames(model1$Sol); #1:3 are explanatory variables;
EVs1 <- model1$Sol[, 1:3];
#intCultivars1 <- model1$Sol[, 13:18];
#intFields1 <- model1$Sol[, 19:48];
intProcYears1 <- model1$Sol[, 4:5];
intPrevCrops1 <- model1$Sol[, 6:9];

length(posterior.mode(intFields1)); #38;  #mean intercept;

### https://github.com/tmalsburg/MCMCglmm-intro
###

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farmRecordsRainfall_fieldpea_model1_VCV_farm.pdf", sep=''), bg="transparent", width=10, height=8, pointsize=12, family="FreeSans");
par(mfrow=c(length(colnames(model1$VCV)),2), mar=c(2,2,1,0));
plot(model1$VCV, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farmRecordsRainfall_fieldpea_model1_sol_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(EVs1)),2), mar=c(2,2,1,0));
plot(EVs1, auto.layout=F);
dev.off()

#pdf.options(encoding='CP1250');
#cairo_pdf(file=paste("farmRecordsRainfall_fieldpea_model1_Cultivars_farm.pdf", sep=''), bg="transparent", width=11, height=30, pointsize=6, family="FreeSans");
#par(mfrow=c(length(colnames(intCultivars1)),2), mar=c(2,2,1,0));
#plot(intCultivars1, auto.layout=F);
#dev.off()

#pdf.options(encoding='CP1250');
#cairo_pdf(file=paste("farmRecordsRainfall_fieldpea_model1_Fields_farm.pdf", sep=''), bg="transparent", width=11, height=30, pointsize=6, family="FreeSans");
#par(mfrow=c(length(colnames(intFields1)),2), mar=c(2,2,1,0));
#plot(intFields1, auto.layout=F);
#dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farmRecordsRainfall_fieldpea_model1_ProcYear_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intProcYears1)),2), mar=c(2,2,1,0));
plot(intProcYears1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farmRecordsRainfall_fieldpea_model1_PrevCrop_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intPrevCrops1)),2), mar=c(2,2,1,0));
plot(intPrevCrops1, auto.layout=F);
dev.off()

colnames(model1$VCV);
#ICC_Cultivar <- model1$VCV[, 1]/(rowSums(model1$VCV)); 
#ICC_Field <- model1$VCV[, 2]/(rowSums(model1$VCV)); 
ICC_ProcYear <- model1$VCV[, 1]/(rowSums(model1$VCV)); 
ICC_PrevCrop <- model1$VCV[, 2]/(rowSums(model1$VCV)); 
#dft1.1 <- cbind(ICC = posterior.mode(ICC_Cultivar), CI = HPDinterval(ICC_Cultivar));
#dft2.1 <- cbind(ICC = posterior.mode(ICC_Field), CI = HPDinterval(ICC_Field));
dft3.1 <- cbind(ICC = posterior.mode(ICC_ProcYear), CI = HPDinterval(ICC_ProcYear));
dft4.1 <- cbind(ICC = posterior.mode(ICC_PrevCrop), CI = HPDinterval(ICC_PrevCrop));
dft1 <- rbind(dft3.1, dft4.1);
row.names(dft1) <- c("ICC_ProcYear", "ICC_PrevCrop");
options(scipen=3); #suppresses scientific notation;
writeLines(capture.output(dft1), paste("farmRecordsRainfall_ICC_fieldpea_Model_1.txt", sep=""));

#---;
####-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----####;
####-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----####;


####-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----####;
####-----@-----@-----@-----@--- Does for lentils ---@-----@-----@-----@-----@-----####;
farmRecordsWeatherGS1 <- farmRecordsWeatherGS[which(farmRecordsWeatherGS$Crop == "Lentils"), ];
str(farmRecordsWeatherGS1);

##remove field records with missing values for Yield (i.e. summerfallow);
farmRecordsWeatherGS1 <- farmRecordsWeatherGS1[which(!is.na(farmRecordsWeatherGS1$RF_sum)), ];
                 
length(farmRecordsWeatherGS1$SeedRate);
summary(farmRecordsWeatherGS1$Yield);
summary(farmRecordsWeatherGS1$SeedRate);
summary(farmRecordsWeatherGS1$NRatio);
summary(farmRecordsWeatherGS1$PRatio);
summary(farmRecordsWeatherGS1$KRatio);
summary(farmRecordsWeatherGS1$SRatio);
summary(farmRecordsWeatherGS1$InocType);
summary(farmRecordsWeatherGS1$Mat_Days);
summary(farmRecordsWeatherGS1$SwaDesc);
summary(farmRecordsWeatherGS1$Drought);
summary(farmRecordsWeatherGS1$Flood);
summary(farmRecordsWeatherGS1$Hail);
summary(farmRecordsWeatherGS1$Copper);
summary(farmRecordsWeatherGS1$Fungicide);
summary(farmRecordsWeatherGS1$Pests);
summary(farmRecordsWeatherGS1$TreatedSeed);
summary(farmRecordsWeatherGS1$Ave_Temp);
summary(farmRecordsWeatherGS1$Max_Temp);
summary(farmRecordsWeatherGS1$Min_Temp);
summary(farmRecordsWeatherGS1$RF_sum);
length(which(summary(farmRecordsWeatherGS1$Cultivar) > 0));
length(which(summary(farmRecordsWeatherGS1$Field) > 0));
length(which(summary(farmRecordsWeatherGS1$Year) > 0));
length(which(summary(farmRecordsWeatherGS1$prevCrop) > 0));
##below sd() doesn't work with NA values;
sd(farmRecordsWeatherGS1$Yield);
sd(farmRecordsWeatherGS1$SeedRate);
sd(farmRecordsWeatherGS1$NRatio);
sd(farmRecordsWeatherGS1$PRatio);
sd(farmRecordsWeatherGS1$KRatio);
sd(farmRecordsWeatherGS1$SRatio);
sd(farmRecordsWeatherGS1$InocType);
sd(farmRecordsWeatherGS1$Mat_Days);
sd(farmRecordsWeatherGS1$SwaDesc);
sd(farmRecordsWeatherGS1$Drought);
sd(farmRecordsWeatherGS1$Flood);
sd(farmRecordsWeatherGS1$Hail);
sd(farmRecordsWeatherGS1$Copper);
sd(farmRecordsWeatherGS1$Fungicide);
sd(farmRecordsWeatherGS1$Pests);
sd(farmRecordsWeatherGS1$TreatedSeed);
sd(farmRecordsWeatherGS1$Ave_Temp);
sd(farmRecordsWeatherGS1$Max_Temp);
sd(farmRecordsWeatherGS1$Min_Temp);
sd(farmRecordsWeatherGS1$RF_sum);

options(scipen=3); #suppresses scientific notation;

#model1 <- MCMCglmm(Yield ~ SeedRate + NRatio + PRatio + KRatio + SRatio + InocType + Mat_Days + SwaDesc + Fungicide + Pests + TreatedSeed + Ave_Temp + RF_sum + RF_sum*Ave_Temp, random=~Cultivar + Field + Year + prevCrop, family="gaussian", data=farmRecordsWeatherGS1, saveX=TRUE, verbose=TRUE, burnin=1388889, nitt=13888889, thin=6250, pr=TRUE, prior=priorFarm1); #2000 samples;
#model1 <- MCMCglmm(Yield ~ SeedRate + NRatio + PRatio + KRatio + SRatio + InocType + Mat_Days + SwaDesc + Fungicide + Pests + TreatedSeed + Ave_Temp + RF_sum + Ave_Temp*RF_sum, random=~Cultivar + Field + Year + prevCrop, family="gaussian", data=farmRecordsWeatherGS1, saveX=TRUE, verbose=TRUE, burnin=500, nitt=2000, thin=25, pr=TRUE, prior=priorFarm1); #2000 samples;
#model1 <- MCMCglmm(Yield ~ SeedRate + NRatio + PRatio + KRatio + SRatio + InocType + Mat_Days + SwaDesc + Fungicide + Pests + TreatedSeed + Ave_Temp + RF_sum + Ave_Temp*RF_sum, random=~Cultivar + Field + Year + prevCrop, family="gaussian", data=farmRecordsWeatherGS1, saveX=TRUE, verbose=TRUE, burnin=1388889, nitt=13888889, thin=6250, pr=TRUE, prior=priorFarm1); #2000 samples;

#model1 <- MCMCglmm(Yield ~ SeedRate + NRatio + PRatio + KRatio + SRatio + InocType + Mat_Days + Drought + Flood + Hail + Fungicide + Pests + TreatedSeed + Ave_Temp + Max_Temp + Min_Temp + RF_sum, random=~Cultivar + Field + Year + prevCrop, family="gaussian", data=farmRecordsWeatherGS1, saveX=TRUE, verbose=TRUE, burnin=1388889, nitt=13888889, thin=6250, pr=TRUE, prior=priorFarm1); #2000 samples;
#model1 <- MCMCglmm(Yield ~ SeedRate + NRatio + PRatio + KRatio + SRatio + InocType + Mat_Days + Drought + Flood + Hail + Fungicide + Pests + TreatedSeed + Ave_Temp + Max_Temp + Min_Temp + RF_sum, random=~Cultivar + Field + Year + prevCrop, family="gaussian", data=farmRecordsWeatherGS1, saveX=TRUE, verbose=TRUE, burnin=500, nitt=2000, thin=25, pr=TRUE, prior=priorFarm1); #2000 samples;
###model1 <- MCMCglmm(Yield ~ SeedRate + NRatio + PRatio + KRatio + SRatio + InocType + Mat_Days + Drought + Flood + Hail + Copper + Fungicide + Pests + TreatedSeed + Ave_Temp + Max_Temp + Min_Temp + RF_sum, random=~Cultivar + Field + Year + prevCrop, family="gaussian", data=farmRecordsWeatherGS1, saveX=TRUE, verbose=TRUE, burnin=500, nitt=2000, thin=25, pr=TRUE, prior=priorFarm1);


#save(model1, file=paste("farmRecordsRainfall_summary_lentils_Model_1.RData", sep=""), compress="xz");
load(paste("farmRecordsRainfall_summary_lentils_Model_1.RData", sep=""));
#load(paste("STR_model_published/farmRecordsRainfall_summary_lentils_Model_1.RData", sep=""));

#https://stats.stackexchange.com/questions/307021/why-does-mcmcglmm-result-in-small-effective-sample-sizes-for-logistic-regression/310291
#autocorr.diag(model1$VCV); #use this to choose thinning interval, autocorrelation less than 0.01 ideally;

fe1 <- summary(model1);
#str(fe1$solutions);
fe1$solutions[,1] <- fe1$solutions[,1];
fe1$solutions[,2] <- fe1$solutions[,2];
fe1$solutions[,3] <- fe1$solutions[,3];

options(scipen=3); #suppresses scientific notation;
writeLines(capture.output(fe1), paste("farmRecordsRainfall_lentils_Model_1.txt", sep=""));
writeLines(capture.output(summary(model1)), paste("farmRecordsRainfall_lentils_Model_1.txt", sep=""));

### DO THIS FOR DIAGNOSTICS OF HOW CHAIN FIT with mcmcglmm;
#http://www.maths.bath.ac.uk/~jjf23/mixchange/split.html

#colnames(model1$Sol); #1:15 are explanatory variables;
EVs1 <- model1$Sol[, 1:15];
intCultivars1 <- model1$Sol[, 16:22];
intFields1 <- model1$Sol[, 23:60];
intProcYears1 <- model1$Sol[, 61:84];
intPrevCrops1 <- model1$Sol[, 85:87];

length(posterior.mode(intFields1)); #38;  #mean intercept;

### https://github.com/tmalsburg/MCMCglmm-intro
###

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farmRecordsRainfall_lentils_model1_VCV_farm.pdf", sep=''), bg="transparent", width=10, height=8, pointsize=12, family="FreeSans");
par(mfrow=c(length(colnames(model1$VCV)),2), mar=c(2,2,1,0));
plot(model1$VCV, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farmRecordsRainfall_lentils_model1_sol_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(EVs1)),2), mar=c(2,2,1,0));
plot(EVs1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farmRecordsRainfall_lentils_model1_Cultivars_farm.pdf", sep=''), bg="transparent", width=11, height=30, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intCultivars1)),2), mar=c(2,2,1,0));
plot(intCultivars1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farmRecordsRainfall_lentils_model1_Fields_farm.pdf", sep=''), bg="transparent", width=11, height=30, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intFields1)),2), mar=c(2,2,1,0));
plot(intFields1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farmRecordsRainfall_lentils_model1_ProcYear_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intProcYears1)),2), mar=c(2,2,1,0));
plot(intProcYears1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farmRecordsRainfall_lentils_model1_PrevCrop_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intPrevCrops1)),2), mar=c(2,2,1,0));
plot(intPrevCrops1, auto.layout=F);
dev.off()

colnames(model1$VCV);
ICC_Cultivar <- model1$VCV[, 1]/(rowSums(model1$VCV)); 
ICC_Field <- model1$VCV[, 2]/(rowSums(model1$VCV)); 
ICC_ProcYear <- model1$VCV[, 3]/(rowSums(model1$VCV)); 
ICC_PrevCrop <- model1$VCV[, 4]/(rowSums(model1$VCV)); 
dft1.1 <- cbind(ICC = posterior.mode(ICC_Cultivar), CI = HPDinterval(ICC_Cultivar));
dft2.1 <- cbind(ICC = posterior.mode(ICC_Field), CI = HPDinterval(ICC_Field));
dft3.1 <- cbind(ICC = posterior.mode(ICC_ProcYear), CI = HPDinterval(ICC_ProcYear));
dft4.1 <- cbind(ICC = posterior.mode(ICC_PrevCrop), CI = HPDinterval(ICC_PrevCrop));
dft1 <- rbind(dft1.1, dft2.1, dft3.1, dft4.1);
row.names(dft1) <- c("ICC_Cultivar", "ICC_Field", "ICC_ProcYear", "ICC_PrevCrop");
options(scipen=3); #suppresses scientific notation;
writeLines(capture.output(dft1), paste("farmRecordsRainfall_ICC_lentils_Model_1.txt", sep=""));

#---;
####-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----####;
####-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----####;


####-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----####;
####-----@-----@-----@-----@--- Does for wheat ---@-----@-----@-----@-----@-----####;
farmRecordsWeatherGS1 <- farmRecordsWeatherGS[which(farmRecordsWeatherGS$Crop == "Wheat"), ];
str(farmRecordsWeatherGS1);

##remove field records with missing values for Yield (i.e. summerfallow);
farmRecordsWeatherGS1 <- farmRecordsWeatherGS1[which(!is.na(farmRecordsWeatherGS1$RF_sum)), ];
                 
length(farmRecordsWeatherGS1$SeedRate);
summary(farmRecordsWeatherGS1$Yield);
summary(farmRecordsWeatherGS1$SeedRate);
summary(farmRecordsWeatherGS1$NRatio);
summary(farmRecordsWeatherGS1$PRatio);
summary(farmRecordsWeatherGS1$KRatio);
summary(farmRecordsWeatherGS1$SRatio);
summary(farmRecordsWeatherGS1$InocType);
summary(farmRecordsWeatherGS1$Mat_Days);
summary(farmRecordsWeatherGS1$SwaDesc);
summary(farmRecordsWeatherGS1$Drought);
summary(farmRecordsWeatherGS1$Flood);
summary(farmRecordsWeatherGS1$Hail);
summary(farmRecordsWeatherGS1$Copper);
summary(farmRecordsWeatherGS1$Fungicide);
summary(farmRecordsWeatherGS1$Pests);
summary(farmRecordsWeatherGS1$TreatedSeed);
summary(farmRecordsWeatherGS1$Ave_Temp);
summary(farmRecordsWeatherGS1$Max_Temp);
summary(farmRecordsWeatherGS1$Min_Temp);
summary(farmRecordsWeatherGS1$RF_sum);
length(which(summary(farmRecordsWeatherGS1$Cultivar) > 0));
length(which(summary(farmRecordsWeatherGS1$Field) > 0));
length(which(summary(farmRecordsWeatherGS1$Year) > 0));
length(which(summary(farmRecordsWeatherGS1$prevCrop) > 0));
##below sd() doesn't work with NA values;
sd(farmRecordsWeatherGS1$Yield);
sd(farmRecordsWeatherGS1$SeedRate);
sd(farmRecordsWeatherGS1$NRatio);
sd(farmRecordsWeatherGS1$PRatio);
sd(farmRecordsWeatherGS1$KRatio);
sd(farmRecordsWeatherGS1$SRatio);
sd(farmRecordsWeatherGS1$InocType);
sd(farmRecordsWeatherGS1$Mat_Days);
sd(farmRecordsWeatherGS1$SwaDesc);
sd(farmRecordsWeatherGS1$Drought);
sd(farmRecordsWeatherGS1$Flood);
sd(farmRecordsWeatherGS1$Hail);
sd(farmRecordsWeatherGS1$Copper);
sd(farmRecordsWeatherGS1$Fungicide);
sd(farmRecordsWeatherGS1$Pests);
sd(farmRecordsWeatherGS1$TreatedSeed);
sd(farmRecordsWeatherGS1$Ave_Temp);
sd(farmRecordsWeatherGS1$Max_Temp);
sd(farmRecordsWeatherGS1$Min_Temp);
sd(farmRecordsWeatherGS1$RF_sum);

options(scipen=3); #suppresses scientific notation;

#model1 <- MCMCglmm(Yield ~ SeedRate + NRatio + PRatio + KRatio + SRatio + Mat_Days + SwaDesc + Copper + Fungicide + Pests + TreatedSeed + Ave_Temp + RF_sum + RF_sum*Ave_Temp, random=~Cultivar + Field + Year + prevCrop, family="gaussian", data=farmRecordsWeatherGS1, saveX=TRUE, verbose=TRUE, burnin=1388889, nitt=13888889, thin=6250, pr=TRUE, prior=priorFarm1); #2000 samples;
#model1 <- MCMCglmm(Yield ~ SeedRate + NRatio + PRatio + KRatio + SRatio + Mat_Days + SwaDesc + Copper + Fungicide + Pests + TreatedSeed + Ave_Temp + RF_sum + Ave_Temp*RF_sum, random=~Cultivar + Field + Year + prevCrop, family="gaussian", data=farmRecordsWeatherGS1, saveX=TRUE, verbose=TRUE, burnin=500, nitt=2000, thin=25, pr=TRUE, prior=priorFarm1); #2000 samples;
#model1 <- MCMCglmm(Yield ~ SeedRate + NRatio + PRatio + KRatio + SRatio + Mat_Days + SwaDesc + Copper + Fungicide + Pests + TreatedSeed + Ave_Temp + RF_sum + Ave_Temp*RF_sum, random=~Cultivar + Field + Year + prevCrop, family="gaussian", data=farmRecordsWeatherGS1, saveX=TRUE, verbose=TRUE, burnin=1388889, nitt=13888889, thin=6250, pr=TRUE, prior=priorFarm1); #2000 samples;

#model1 <- MCMCglmm(Yield ~ SeedRate + NRatio + PRatio + KRatio + SRatio + Mat_Days + Drought + Hail + Copper + Fungicide + Pests + TreatedSeed + Ave_Temp + Max_Temp + Min_Temp + RF_sum, random=~Cultivar + Field + Year + prevCrop, family="gaussian", data=farmRecordsWeatherGS1, saveX=TRUE, verbose=TRUE, burnin=1388889, nitt=13888889, thin=6250, pr=TRUE, prior=priorFarm1); #2000 samples;
#model1 <- MCMCglmm(Yield ~ SeedRate + NRatio + PRatio + KRatio + SRatio + Mat_Days + Drought + Hail + Copper + Fungicide + Pests + TreatedSeed + Ave_Temp + Max_Temp + Min_Temp + RF_sum, random=~Cultivar + Field + Year + prevCrop, family="gaussian", data=farmRecordsWeatherGS1, saveX=TRUE, verbose=TRUE, burnin=500, nitt=2000, thin=25, pr=TRUE, prior=priorFarm1); #2000 samples;
###model1 <- MCMCglmm(Yield ~ SeedRate + NRatio + PRatio + KRatio + SRatio + InocType + Mat_Days + Drought + Flood + Hail + Copper + Fungicide + Pests + TreatedSeed + Ave_Temp + Max_Temp + Min_Temp + RF_sum, random=~Cultivar + Field + Year + prevCrop, family="gaussian", data=farmRecordsWeatherGS1, saveX=TRUE, verbose=TRUE, burnin=500, nitt=2000, thin=25, pr=TRUE, prior=priorFarm1);


#save(model1, file=paste("farmRecordsRainfall_summary_wheat_Model_1.RData", sep=""), compress="xz");
load(paste("farmRecordsRainfall_summary_wheat_Model_1.RData", sep=""));
#load(paste("STR_model_published/farmRecordsRainfall_summary_wheat_Model_1.RData", sep=""));

#https://stats.stackexchange.com/questions/307021/why-does-mcmcglmm-result-in-small-effective-sample-sizes-for-logistic-regression/310291
#autocorr.diag(model1$VCV); #use this to choose thinning interval, autocorrelation less than 0.01 ideally;

fe1 <- summary(model1);
#str(fe1$solutions);
fe1$solutions[,1] <- fe1$solutions[,1];
fe1$solutions[,2] <- fe1$solutions[,2];
fe1$solutions[,3] <- fe1$solutions[,3];

options(scipen=3); #suppresses scientific notation;
writeLines(capture.output(fe1), paste("farmRecordsRainfall_wheat_Model_1.txt", sep=""));
writeLines(capture.output(summary(model1)), paste("farmRecordsRainfall_wheat_Model_1.txt", sep=""));

### DO THIS FOR DIAGNOSTICS OF HOW CHAIN FIT with mcmcglmm;
#http://www.maths.bath.ac.uk/~jjf23/mixchange/split.html

#colnames(model1$Sol); #1:16 are explanatory variables;
EVs1 <- model1$Sol[, 1:16];
intCultivars1 <- model1$Sol[, 17:24];
intFields1 <- model1$Sol[, 25:62];
intProcYears1 <- model1$Sol[, 63:95];
intPrevCrops1 <- model1$Sol[, 96:101];

length(posterior.mode(intFields1)); #38;  #mean intercept;

### https://github.com/tmalsburg/MCMCglmm-intro
###

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farmRecordsRainfall_wheat_model1_VCV_farm.pdf", sep=''), bg="transparent", width=10, height=8, pointsize=12, family="FreeSans");
par(mfrow=c(length(colnames(model1$VCV)),2), mar=c(2,2,1,0));
plot(model1$VCV, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farmRecordsRainfall_wheat_model1_sol_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(EVs1)),2), mar=c(2,2,1,0));
plot(EVs1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farmRecordsRainfall_wheat_model1_Cultivars_farm.pdf", sep=''), bg="transparent", width=11, height=30, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intCultivars1)),2), mar=c(2,2,1,0));
plot(intCultivars1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farmRecordsRainfall_wheat_model1_Fields_farm.pdf", sep=''), bg="transparent", width=11, height=30, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intFields1)),2), mar=c(2,2,1,0));
plot(intFields1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farmRecordsRainfall_wheat_model1_ProcYear_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intProcYears1)),2), mar=c(2,2,1,0));
plot(intProcYears1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farmRecordsRainfall_wheat_model1_PrevCrop_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intPrevCrops1)),2), mar=c(2,2,1,0));
plot(intPrevCrops1, auto.layout=F);
dev.off()

colnames(model1$VCV);
ICC_Cultivar <- model1$VCV[, 1]/(rowSums(model1$VCV)); 
ICC_Field <- model1$VCV[, 2]/(rowSums(model1$VCV)); 
ICC_ProcYear <- model1$VCV[, 3]/(rowSums(model1$VCV)); 
ICC_PrevCrop <- model1$VCV[, 4]/(rowSums(model1$VCV)); 
dft1.1 <- cbind(ICC = posterior.mode(ICC_Cultivar), CI = HPDinterval(ICC_Cultivar));
dft2.1 <- cbind(ICC = posterior.mode(ICC_Field), CI = HPDinterval(ICC_Field));
dft3.1 <- cbind(ICC = posterior.mode(ICC_ProcYear), CI = HPDinterval(ICC_ProcYear));
dft4.1 <- cbind(ICC = posterior.mode(ICC_PrevCrop), CI = HPDinterval(ICC_PrevCrop));
dft1 <- rbind(dft1.1, dft2.1, dft3.1, dft4.1);
row.names(dft1) <- c("ICC_Cultivar", "ICC_Field", "ICC_ProcYear", "ICC_PrevCrop");
options(scipen=3); #suppresses scientific notation;
writeLines(capture.output(dft1), paste("farmRecordsRainfall_ICC_wheat_Model_1.txt", sep=""));

#---;
####-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----####;
####-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----####;

