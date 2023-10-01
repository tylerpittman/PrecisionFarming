## Creates a nice plot of farm weather data (ILACADEN2 and ILACADEN3) [Jul 2018 to Jun 2019]; 
## Swift Current weather data, and farm rainfall data from Jan 2015 to Dec 2018;
## Tyler Pittman, 22 July 2019
# Rscript --no-save /Users/tylerpittman/Farm/STEP_15_farm_weather_data_paper_15July2019.R

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
library(reshape2);
library(ggplot2);
library(plyr); 		#gives multicolumn frequencies for tables with ddply() and join_all function;
library(dplyr);  	#need dplyr for the bind_rows function 
library(cobs);		# fits nicer splines that can be strictly increasing;
library(parallel);	# Does multicore processing when using mclapply WORKS on Linux!
#library(greenbrown); #has many complicated dependencies to install, for equal area parcel division of shapefiles;
#source("shape2poly.R"); # Reads in shape2poly function and others
#source("polygonizer.R"); # Reads in polygonizer function
source("http://www.math.mcmaster.ca/bolker/R/misc/legendx.R"); # custom source allows box.cex() function to be used in legend() for plot!;

#weatherFarmJuly2018toJune2019

ILACADEN2 <- read_excel("weatherFarmJuly2018toJune2019.xlsx", sheet = "ILACADEN2");
ILACADEN3 <- read_excel("weatherFarmJuly2018toJune2019.xlsx", sheet = "ILACADEN3");

ILACADEN2 <- as.data.frame(ILACADEN2);
ILACADEN3 <- as.data.frame(ILACADEN3);

str(ILACADEN2);
str(ILACADEN3);

ILACADEN2$Date <-  as.Date(paste(as.character(ILACADEN2$Date), sep=""), format="%m/%d/%Y");
ILACADEN3$Date <-  as.Date(paste(as.character(ILACADEN3$Date), sep=""), format="%m/%d/%Y");

##gsub("[a-zA-Z ]", "", ILACADEN2$DryTempH);
ILACADEN2$DryTempH <- as.numeric(gsub("[a-zA-Z ]", "", ILACADEN2$DryTempH));
ILACADEN2$DryTempA <- as.numeric(gsub("[a-zA-Z ]", "", ILACADEN2$DryTempA));
ILACADEN2$DryTempL <- as.numeric(gsub("[a-zA-Z ]", "", ILACADEN2$DryTempL));
ILACADEN2$DewPointH <- as.numeric(gsub("[a-zA-Z ]", "", ILACADEN2$DewPointH));
ILACADEN2$DewPointA <- as.numeric(gsub("[a-zA-Z ]", "", ILACADEN2$DewPointA));
ILACADEN2$DewPointL <- as.numeric(gsub("[a-zA-Z ]", "", ILACADEN2$DewPointL));
ILACADEN2$WindH <- as.numeric(gsub("[a-zA-Z /]", "", ILACADEN2$WindH));
ILACADEN2$WindA <- as.numeric(gsub("[a-zA-Z /]", "", ILACADEN2$WindA));
ILACADEN2$WindL <- as.numeric(gsub("[a-zA-Z /]", "", ILACADEN2$WindL));
ILACADEN2$PressureH <- as.numeric(gsub("[a-zA-Z ]", "", ILACADEN2$PressureH));
ILACADEN2$PressureL <- as.numeric(gsub("[a-zA-Z ]", "", ILACADEN2$PressureL));
ILACADEN2$PrecipSum <- as.numeric(gsub("[a-zA-Z ]", "", ILACADEN2$PrecipSum));

ILACADEN3$DryTempH <- as.numeric(gsub("[a-zA-Z ]", "", ILACADEN3$DryTempH));
ILACADEN3$DryTempA <- as.numeric(gsub("[a-zA-Z ]", "", ILACADEN3$DryTempA));
ILACADEN3$DryTempL <- as.numeric(gsub("[a-zA-Z ]", "", ILACADEN3$DryTempL));
ILACADEN3$DewPointH <- as.numeric(gsub("[a-zA-Z ]", "", ILACADEN3$DewPointH));
ILACADEN3$DewPointA <- as.numeric(gsub("[a-zA-Z ]", "", ILACADEN3$DewPointA));
ILACADEN3$DewPointL <- as.numeric(gsub("[a-zA-Z ]", "", ILACADEN3$DewPointL));
ILACADEN3$WindH <- as.numeric(gsub("[a-zA-Z /]", "", ILACADEN3$WindH));
ILACADEN3$WindA <- as.numeric(gsub("[a-zA-Z /]", "", ILACADEN3$WindA));
ILACADEN3$WindL <- as.numeric(gsub("[a-zA-Z /]", "", ILACADEN3$WindL));
ILACADEN3$PressureH <- as.numeric(gsub("[a-zA-Z ]", "", ILACADEN3$PressureH));
ILACADEN3$PressureL <- as.numeric(gsub("[a-zA-Z ]", "", ILACADEN3$PressureL));
ILACADEN3$PrecipSum <- as.numeric(gsub("[a-zA-Z ]", "", ILACADEN3$PrecipSum));

ILACADEN2$Station <- "ILACADEN2";
ILACADEN3$Station <- "ILACADEN3";

#str(ILACADEN2);
#str(ILACADEN3);

weather0 <- rbind(ILACADEN2, ILACADEN3);
weather0$Station <- as.factor(weather0$Station);
str(weather0);

weather0$Month <- months(weather0$Date);
weather0$Year <- format(weather0$Date, format="%Y");
maxTempM <- aggregate(DryTempH ~ Month + Year + Station, weather0, max);
minTempM <- aggregate(DryTempL ~ Month + Year + Station, weather0, min);
avgTempM <- aggregate(DryTempA ~ Month + Year + Station, weather0, mean);
sumPrecipM <- aggregate(PrecipSum ~ Month + Year + Station, weather0, sum);

weather1 <- join_all(list(maxTempM, minTempM, avgTempM, sumPrecipM), by = c('Month', 'Year', 'Station'), type = "left"); 

annTemp <- aggregate(DryTempA ~ Station, weather0, mean);
colnames(annTemp) <- c("Station", "AnnTemp");

weather1 <- join_all(list(weather1, annTemp), by = c('Station'), type = "left"); 
weather1$Date <- as.Date(paste(as.character(weather1$Month), "/01/", as.character(weather1$Year), sep=""), format="%B/%d/%Y");
weather1 <- weather1[order(weather1$Station, weather1$Date), ]; #sorts by Station and Date for below cumsum;
weather1$cummPrecip <- ave(weather1$PrecipSum, list(weather1$Station), FUN=cumsum); 



##---##---##---##---##---##---##---##---##---##---##---##---##;
### This loads Weather Office data;
##http://climate.weather.gc.ca/climate_data/daily_data_e.html?hlyRange=2011-06-27%7C2019-07-14&dlyRange=2011-06-30%7C2019-07-13&mlyRange=%7C&StationID=48976&Prov=SK&urlExtension=_e.html&searchType=stnName&optLimit=yearRange&StartYear=1840&EndYear=2019&selRowPerPage=25&Line=1&searchMethod=contains&Month=7&Day=16&txtStationName=Swift+current&timeframe=2&Year=2019

sc2015 <- read.csv(paste("sc2015.csv", sep=""), fileEncoding="latin1", header=F, skip=25);
sc2015 <- sc2015[ , c(1,3,4,6,8,10,20,24,26)];
colnames(sc2015) <- c("date", "month", "day", "maxTemp", "minTemp", "meanTemp", "totalPrecip", "maxWinDir", "maxWindSpeed");
sc2015$date <- as.Date(sc2015$date, "%Y-%m-%d");
sc2015$month <- as.integer(sc2015$month);
sc2015$day <- as.integer(sc2015$day);
sc2015$maxTemp <- as.numeric(sc2015$maxTemp);
sc2015$minTemp <- as.numeric(sc2015$minTemp);
sc2015$meanTemp <- as.numeric(sc2015$meanTemp);
sc2015$totalPrecip <- as.numeric(sc2015$totalPrecip);
sc2015$maxWinDir <- as.integer(sc2015$maxWinDir);
sc2015$maxWindSpeed <- as.integer(sc2015$maxWindSpeed);

sc2016 <- read.csv(paste("sc2016.csv", sep=""), fileEncoding="latin1", header=F, skip=25);
sc2016 <- sc2016[ , c(1,3,4,6,8,10,20,24,26)];
colnames(sc2016) <- c("date", "month", "day", "maxTemp", "minTemp", "meanTemp", "totalPrecip", "maxWinDir", "maxWindSpeed");
sc2016$date <- as.Date(sc2016$date, "%Y-%m-%d");
sc2016$month <- as.integer(sc2016$month);
sc2016$day <- as.integer(sc2016$day);
sc2016$maxTemp <- as.numeric(sc2016$maxTemp);
sc2016$minTemp <- as.numeric(sc2016$minTemp);
sc2016$meanTemp <- as.numeric(sc2016$meanTemp);
sc2016$totalPrecip <- as.numeric(sc2016$totalPrecip);
sc2016$maxWinDir <- as.integer(sc2016$maxWinDir);
sc2016$maxWindSpeed <- as.integer(sc2016$maxWindSpeed);

sc2017 <- read.csv(paste("sc2017.csv", sep=""), fileEncoding="latin1", header=F, skip=25);
sc2017 <- sc2017[ , c(1,3,4,6,8,10,20,24,26)];
colnames(sc2017) <- c("date", "month", "day", "maxTemp", "minTemp", "meanTemp", "totalPrecip", "maxWinDir", "maxWindSpeed");
sc2017$date <- as.Date(sc2017$date, "%Y-%m-%d");
sc2017$month <- as.integer(sc2017$month);
sc2017$day <- as.integer(sc2017$day);
sc2017$maxTemp <- as.numeric(sc2017$maxTemp);
sc2017$minTemp <- as.numeric(sc2017$minTemp);
sc2017$meanTemp <- as.numeric(sc2017$meanTemp);
sc2017$totalPrecip <- as.numeric(sc2017$totalPrecip);
sc2017$maxWinDir <- as.integer(sc2017$maxWinDir);
sc2017$maxWindSpeed <- as.integer(sc2017$maxWindSpeed);

sc2018 <- read.csv(paste("sc2018.csv", sep=""), fileEncoding="latin1", header=F, skip=26);
sc2018 <- sc2018[ , c(1,3,4,6,8,10,20,24,26)];
colnames(sc2018) <- c("date", "month", "day", "maxTemp", "minTemp", "meanTemp", "totalPrecip", "maxWinDir", "maxWindSpeed");
sc2018$date <- as.Date(sc2018$date, "%Y-%m-%d");
sc2018$month <- as.integer(sc2018$month);
sc2018$day <- as.integer(sc2018$day);
sc2018$maxTemp <- as.numeric(sc2018$maxTemp);
sc2018$minTemp <- as.numeric(sc2018$minTemp);
sc2018$meanTemp <- as.numeric(sc2018$meanTemp);
sc2018$totalPrecip <- as.numeric(sc2018$totalPrecip);
sc2018$maxWinDir <- as.integer(sc2018$maxWinDir);
sc2018$maxWindSpeed <- as.integer(sc2018$maxWindSpeed);

wo0 <- rbind(sc2015, sc2016, sc2017, sc2018);
wo0$Station <- "71142";
wo0$Station <- as.factor(wo0$Station);
wo0$Month <- months(wo0$date);
wo0$Year <- format(wo0$date, format="%Y");
wo0$DryTempH <- wo0$maxTemp;
wo0$DryTempL <- wo0$minTemp;
wo0$DryTempA <- wo0$meanTemp;
wo0$PrecipSum <- wo0$totalPrecip;
str(wo0);

### This loads farm rainfall data;
farmRF <- read_excel("rainfallFarm2015to2018.xlsx", sheet="RF");
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

wo1 <- join_all(list(maxTempM, minTempM, avgTempM, sumPrecipM), by = c('Month', 'Year', 'Station'), type = "left"); 

annTemp <- aggregate(DryTempA ~ Station, wo0, mean);
colnames(annTemp) <- c("Station", "AnnTemp");

wo1 <- join_all(list(wo1, annTemp), by = c('Station'), type = "left"); 
wo1$Date <- as.Date(paste(as.character(wo1$Month), "/01/", as.character(wo1$Year), sep=""), format="%B/%d/%Y");
wo1 <- wo1[order(wo1$Station, wo1$Date), ]; #sorts by Station and Date for below cumsum;
#wo1$cummPrecip <- ave(wo1$PrecipSum, list(wo1$Station), FUN=cumsum); 
wo1$cummPrecip <- ave(wo1$PrecipSum, list(wo1$Station, wo1$Year), FUN=cumsum); 
str(wo1);



#### Farm rainfall for each year 2015 to 2018;
sumPrecipM <- aggregate(PrecipSum ~ Month + Year, farmRF, sum);
sumPrecipM$Station <- "71142";
wo4 <- join_all(list(maxTempM, minTempM, avgTempM, sumPrecipM), by = c('Month', 'Year', 'Station'), type = "left"); 
annTemp <- aggregate(DryTempA ~ Station, wo0, mean);
colnames(annTemp) <- c("Station", "AnnTemp");

##below takes maxTempM, minTempM and avgTempM from above section;
wo4 <- join_all(list(wo4, annTemp), by = c('Station'), type = "left"); 
wo4$Date <- as.Date(paste(as.character(wo4$Month), "/01/", as.character(wo4$Year), sep=""), format="%B/%d/%Y");
wo4 <- wo4[order(wo4$Station, wo4$Date), ]; #sorts by Station and Date for below cumsum;
wo4[which(is.na(wo4$PrecipSum)), ]$PrecipSum <- 0;
wo4$cummPrecip <- ave(wo4$PrecipSum, list(wo4$Station, wo4$Year), FUN=cumsum); 
wo4[which(wo4$Month == "December"), ]$cummPrecip <- 0;
wo4[which(wo4$cummPrecip == 0), ]$cummPrecip <- NA;
str(wo4);

###-------------------------###;

maxTempM <- aggregate(DryTempH ~ Month + Station, wo0, max);
minTempM <- aggregate(DryTempL ~ Month + Station, wo0, min);
avgTempM <- aggregate(DryTempA ~ Month + Station, wo0, mean);
sumPrecipM0 <- aggregate(PrecipSum ~ Month + Year + Station, wo0, sum);
sumPrecipM <- aggregate(PrecipSum ~ Month + Station, sumPrecipM0, mean);

wo2 <- join_all(list(maxTempM, minTempM, avgTempM, sumPrecipM), by = c('Month', 'Station'), type = "left"); 

annTemp <- aggregate(DryTempA ~ Station, wo0, mean);
colnames(annTemp) <- c("Station", "AnnTemp");

wo2 <- join_all(list(wo2, annTemp), by = c('Station'), type = "left"); 
wo2$Date <- as.Date(paste(as.character(wo2$Month), "/01/", "1900", sep=""), format="%B/%d/%Y");
wo2 <- wo2[order(wo2$Station, wo2$Date), ]; #sorts by Station and Date for below cumsum;
wo2$cummPrecip <- ave(wo2$PrecipSum, list(wo2$Station), FUN=cumsum);
str(wo2);


#### Farm rainfall average for each month 2015 to 2018;
sumPrecipM0 <- aggregate(PrecipSum ~ Month + Year + Station, farmRF, sum);
sumPrecipM <- aggregate(PrecipSum ~ Month, sumPrecipM0, mean);

##below takes maxTempM, minTempM and avgTempM from above section;
wo3 <- join_all(list(maxTempM, minTempM, avgTempM, sumPrecipM), by = c('Month'), type = "left"); 
wo3$Date <- as.Date(paste(as.character(wo3$Month), "/01/", "1900", sep=""), format="%B/%d/%Y");
wo3 <- wo3[order(wo3$Station, wo3$Date), ]; #sorts by Station and Date for below cumsum;
wo3[which(is.na(wo3$PrecipSum)), ]$PrecipSum <- 0;
wo3$cummPrecip <- ave(wo3$PrecipSum, list(wo3$Station), FUN=cumsum);
wo3[which(wo3$Month == "December"), ]$cummPrecip <- 0;
wo3[which(wo3$cummPrecip == 0), ]$cummPrecip <- NA;
str(wo3);

summary(wo3$DryTempH); #max is 39.8 C;
summary(wo3$DryTempL); #low is -36.8 C;
##---##---##---##---##---##---##---##---##---##---##---##---##;


#---------------#---------------#---------------#---------------#---------------#;    
###########;
### below for ggplot doesn't work well for climatogram;
#https://stackoverflow.com/questions/3099219/ggplot-with-2-y-axes-on-each-side-and-different-scales    
#https://stackoverflow.com/questions/6142944/how-can-i-plot-with-2-different-y-axes

weather1a <- weather1[which(weather1$Station == "ILACADEN2"), ];
str(weather1a);

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farm_weather_ILACADEN2_July2018toJune2019", ".pdf", sep=''), bg="transparent", width=8, height=6, pointsize=12, family="FreeSans");
par(mar=c(5, 5, 4, 6) + 0.1);
plot(weather1a$Date, weather1a$cummPrecip, pch=2, axes=FALSE, ylim=c(0, ceiling(max(weather1a$cummPrecip))+15), xlab="", ylab="", type="l", main=" ", lwd=2, col=NA, yaxs="i");
smoothingSpline <- smooth.spline(weather1a$Date, weather1a$cummPrecip, spar=0.25);
loessSpline <- loess(cummPrecip ~ as.numeric(Date), data=weather1a, span=0.50);
lines(smoothingSpline, lwd=2);
axis(2, ylim=c(0, 1), col="black", las=1, at=c(seq(from=0,to=1000,by=25)));  ## las=1 makes horizontal labels;
mtext(expression(paste("Monthly total (RF" ["mtot"], ") and mean monthly", sep="")), side=2, line=3.5);
mtext(expression(paste("cumulative (RF" ["mmc"], ") rainfall (mm)", sep="")), side=2, line=2.5);
box();
par(new=TRUE);
plot(weather1a$Date, weather1a$AnnTemp, lty=2,  xlab="", ylab="", ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), axes=FALSE, type="l", col=NA, yaxs="i");
abline(h=unique(weather1a$AnnTemp), col="black", lty=2);
mtext("Monthly min, max, average temperature and", side=4, col="black", line=2.5); 
mtext(expression(paste("average annual (Avg T" ["ann"], ") temperature (°C)", sep="")), side=4, col="black", line=3.5); 
axis(4, ylim=c(floor(min(weather1a$DryTempL)), ceiling(max(weather1a$DryTempH))), at=c(seq(from=-50,to=50,by=10)), col.axis="black", las=1, bg="white", col="black");
mtext("Months", side=1, col="black", line=2.5);  
par(new=TRUE);
plot(weather1a$Date, weather1a$PrecipSum, main=" ", xlab="", ylab="", axes=FALSE, ylim=c(0, ceiling(max(weather1a$cummPrecip))+15), col="grey", type="h", lend=1, lwd=35, yaxs="i");
#plot(weather1a$Date, weather1a$PrecipSum, main=" ", xlab="", ylab="", axes=FALSE, ylim=c(0, ceiling(max(weather1a$cummPrecip))), col="grey50", type="h", lend=1, lwd=35, yaxs="i");
#lines(weather1a$Date, weather1a$PrecipSum-0.75, type="h", lwd=35-2, col="grey", lend=1);
par(new=TRUE);
plot(weather1a$Date, weather1a$DryTempH, pch=24,  xlab="", ylab="", ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), axes=FALSE, type="p", bg="white", col="black", yaxs="i");
par(new=TRUE);
plot(weather1a$Date, weather1a$DryTempA, pch=22,  xlab="", ylab="", ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), axes=FALSE, type="p", bg="grey50", col="black", yaxs="i");
par(new=TRUE);
plot(weather1a$Date, weather1a$DryTempL, pch=21, ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), xaxt='n', yaxt='n', ann=FALSE, type="p", bg="white", col="black", yaxs="i");
axis(1, weather1a$Date, format(weather1a$Date, "%b"), cex.axis=0.7, col="black", col.axis="black");
legend("bottom", c(expression("RF"[mtot]), expression("RF"[mmc]), "Min T", "Max T", "Avg T", expression("Avg T"[ann])), fill=c("grey", NA, NA, NA, NA, NA), border=c("grey", NA, NA, NA, NA, NA), box.cex=c(2.5, 1.0), pt.cex=c(1, 1), pch=c(NA, NA, 21, 24, 22, NA), pt.bg=c(NA, NA, "white", "white", "grey50", NA), col=c(NA, "black", "black", "black", "black", "black"), lty=c(NA, 1, NA, NA, NA, 2), lwd=c(NA, 2, NA, NA, NA, 1), inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="o", cex=0.6, x.intersp=0.1);
dev.off(); 



weather1a <- weather1[which(weather1$Station == "ILACADEN3"), ];
str(weather1a);

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farm_weather_ILACADEN3_July2018toJune2019", ".pdf", sep=''), bg="transparent", width=8, height=6, pointsize=12, family="FreeSans");
par(mar=c(5, 5, 4, 6) + 0.1);
plot(weather1a$Date, weather1a$cummPrecip, pch=2, axes=FALSE, ylim=c(0, ceiling(max(weather1a$cummPrecip))+15), xlab="", ylab="", type="l", main=" ", lwd=2, col=NA, yaxs="i");
smoothingSpline <- smooth.spline(weather1a$Date, weather1a$cummPrecip, spar=0.25);
loessSpline <- loess(cummPrecip ~ as.numeric(Date), data=weather1a, span=0.50);
lines(smoothingSpline, lwd=2);
axis(2, ylim=c(0, 1), col="black", las=1, at=c(seq(from=0,to=1000,by=25)));  ## las=1 makes horizontal labels;
mtext(expression(paste("Monthly total (RF" ["mtot"], ") and mean monthly", sep="")), side=2, line=3.5);
mtext(expression(paste("cumulative (RF" ["mmc"], ") rainfall (mm)", sep="")), side=2, line=2.5);
box();
par(new=TRUE);
plot(weather1a$Date, weather1a$AnnTemp, lty=2,  xlab="", ylab="", ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), axes=FALSE, type="l", col=NA, yaxs="i");
abline(h=unique(weather1a$AnnTemp), col="black", lty=2);
mtext("Monthly min, max, average temperature and", side=4, col="black", line=2.5); 
mtext(expression(paste("average annual (Avg T" ["ann"], ") temperature (°C)", sep="")), side=4, col="black", line=3.5); 
axis(4, ylim=c(floor(min(weather1a$DryTempL)), ceiling(max(weather1a$DryTempH))), at=c(seq(from=-50,to=50,by=10)), col.axis="black", las=1, bg="white", col="black");
mtext("Months", side=1, col="black", line=2.5);  
par(new=TRUE);
plot(weather1a$Date, weather1a$PrecipSum, main=" ", xlab="", ylab="", axes=FALSE, ylim=c(0, ceiling(max(weather1a$cummPrecip))+15), col="grey", type="h", lend=1, lwd=35, yaxs="i");
#plot(weather1a$Date, weather1a$PrecipSum, main=" ", xlab="", ylab="", axes=FALSE, ylim=c(0, ceiling(max(weather1a$cummPrecip))), col="grey50", type="h", lend=1, lwd=35, yaxs="i");
#lines(weather1a$Date, weather1a$PrecipSum-0.75, type="h", lwd=35-2, col="grey", lend=1);
par(new=TRUE);
plot(weather1a$Date, weather1a$DryTempH, pch=24,  xlab="", ylab="", ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), axes=FALSE, type="p", bg="white", col="black", yaxs="i");
par(new=TRUE);
plot(weather1a$Date, weather1a$DryTempA, pch=22,  xlab="", ylab="", ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), axes=FALSE, type="p", bg="grey50", col="black", yaxs="i");
par(new=TRUE);
plot(weather1a$Date, weather1a$DryTempL, pch=21, ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), xaxt='n', yaxt='n', ann=FALSE, type="p", bg="white", col="black", yaxs="i");
axis(1, weather1a$Date, format(weather1a$Date, "%b"), cex.axis=0.7, col="black", col.axis="black");
legend("bottom", c(expression("RF"[mtot]), expression("RF"[mmc]), "Min T", "Max T", "Avg T", expression("Avg T"[ann])), fill=c("grey", NA, NA, NA, NA, NA), border=c("grey", NA, NA, NA, NA, NA), box.cex=c(2.5, 1.0), pt.cex=c(1, 1), pch=c(NA, NA, 21, 24, 22, NA), pt.bg=c(NA, NA, "white", "white", "grey50", NA), col=c(NA, "black", "black", "black", "black", "black"), lty=c(NA, 1, NA, NA, NA, 2), lwd=c(NA, 2, NA, NA, NA, 1), inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="o", cex=0.6, x.intersp=0.1);
dev.off(); 



#---------------#---------------#---------------#---------------#---------------#;    
weather1a <- wo1[which(wo1$Station == "71142"), ];
str(weather1a);

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farm_weather_71142_Jan2015toDec2018", ".pdf", sep=''), bg="transparent", width=8, height=6, pointsize=12, family="FreeSans");
par(mar=c(5, 5, 4, 6) + 0.1);
plot(weather1a$Date, weather1a$cummPrecip, pch=2, axes=FALSE, ylim=c(0, ceiling(max(weather1a$cummPrecip))+15), xlab="", ylab="", type="l", main=" ", lwd=2, col=NA, yaxs="i");
smoothingSpline <- smooth.spline(weather1a$Date, weather1a$cummPrecip, spar=0.25);
loessSpline <- loess(cummPrecip ~ as.numeric(Date), data=weather1a, span=0.50);
lines(smoothingSpline, lwd=2);
axis(2, ylim=c(0, 1), col="black", las=1, at=c(seq(from=0,to=1000,by=25)));  ## las=1 makes horizontal labels;
mtext(expression(paste("Monthly total (RF" ["mtot"], ") and mean monthly", sep="")), side=2, line=3.5);
mtext(expression(paste("cumulative (RF" ["mmc"], ") rainfall (mm)", sep="")), side=2, line=2.5);
box();
par(new=TRUE);
plot(weather1a$Date, weather1a$AnnTemp, lty=2,  xlab="", ylab="", ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), axes=FALSE, type="l", col=NA, yaxs="i");
abline(h=unique(weather1a$AnnTemp), col="black", lty=2);
mtext("Monthly min, max, average temperature and", side=4, col="black", line=2.5); 
mtext(expression(paste("average annual (Avg T" ["ann"], ") temperature (°C)", sep="")), side=4, col="black", line=3.5); 
axis(4, ylim=c(floor(min(weather1a$DryTempL)), ceiling(max(weather1a$DryTempH))), at=c(seq(from=-50,to=50,by=10)), col.axis="black", las=1, bg="white", col="black");
mtext("Months", side=1, col="black", line=2.5);  
par(new=TRUE);
plot(weather1a$Date, weather1a$PrecipSum, main=" ", xlab="", ylab="", axes=FALSE, ylim=c(0, ceiling(max(weather1a$cummPrecip))+15), col="grey", type="h", lend=1, lwd=35, yaxs="i");
#plot(weather1a$Date, weather1a$PrecipSum, main=" ", xlab="", ylab="", axes=FALSE, ylim=c(0, ceiling(max(weather1a$cummPrecip))), col="grey50", type="h", lend=1, lwd=35, yaxs="i");
#lines(weather1a$Date, weather1a$PrecipSum-0.75, type="h", lwd=35-2, col="grey", lend=1);
par(new=TRUE);
plot(weather1a$Date, weather1a$DryTempH, pch=24,  xlab="", ylab="", ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), axes=FALSE, type="p", bg="white", col="black", yaxs="i");
par(new=TRUE);
plot(weather1a$Date, weather1a$DryTempA, pch=22,  xlab="", ylab="", ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), axes=FALSE, type="p", bg="grey50", col="black", yaxs="i");
par(new=TRUE);
plot(weather1a$Date, weather1a$DryTempL, pch=21, ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), xaxt='n', yaxt='n', ann=FALSE, type="p", bg="white", col="black", yaxs="i");
axis(1, weather1a$Date, format(weather1a$Date, "%b"), cex.axis=0.7, col="black", col.axis="black");
legend("bottom", c(expression("RF"[mtot]), expression("RF"[mmc]), "Min T", "Max T", "Avg T", expression("Avg T"[ann])), fill=c("grey", NA, NA, NA, NA, NA), border=c("grey", NA, NA, NA, NA, NA), box.cex=c(2.5, 1.0), pt.cex=c(1, 1), pch=c(NA, NA, 21, 24, 22, NA), pt.bg=c(NA, NA, "white", "white", "grey50", NA), col=c(NA, "black", "black", "black", "black", "black"), lty=c(NA, 1, NA, NA, NA, 2), lwd=c(NA, 2, NA, NA, NA, 1), inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="o", cex=0.6, x.intersp=0.1);
dev.off(); 


weather1a <- wo2[which(wo2$Station == "71142"), ];
str(weather1a);

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farm_weather_71142_AvgJan2015toDec2018", ".pdf", sep=''), bg="transparent", width=8, height=6, pointsize=12, family="FreeSans");
par(mar=c(5, 5, 4, 6) + 0.1);
plot(weather1a$Date, weather1a$cummPrecip, pch=2, axes=FALSE, ylim=c(0, ceiling(max(weather1a$cummPrecip))+15), xlab="", ylab="", type="l", main=" ", lwd=2, col=NA, yaxs="i");
smoothingSpline <- smooth.spline(weather1a$Date, weather1a$cummPrecip, spar=0.25);
loessSpline <- loess(cummPrecip ~ as.numeric(Date), data=weather1a, span=0.50);
lines(smoothingSpline, lwd=2);
axis(2, ylim=c(0, 1), col="black", las=1, at=c(seq(from=0,to=1000,by=25)));  ## las=1 makes horizontal labels;
mtext(expression(paste("Monthly total (RF" ["mtot"], ") and mean monthly", sep="")), side=2, line=3.5);
mtext(expression(paste("cumulative (RF" ["mmc"], ") rainfall (mm)", sep="")), side=2, line=2.5);
box();
par(new=TRUE);
plot(weather1a$Date, weather1a$AnnTemp, lty=2,  xlab="", ylab="", ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), axes=FALSE, type="l", col=NA, yaxs="i");
abline(h=unique(weather1a$AnnTemp), col="black", lty=2);
mtext("Monthly min, max, average temperature and", side=4, col="black", line=2.5); 
mtext(expression(paste("average annual (Avg T" ["ann"], ") temperature (°C)", sep="")), side=4, col="black", line=3.5); 
axis(4, ylim=c(floor(min(weather1a$DryTempL)), ceiling(max(weather1a$DryTempH))), at=c(seq(from=-50,to=50,by=10)), col.axis="black", las=1, bg="white", col="black");
mtext("Months", side=1, col="black", line=2.5);  
par(new=TRUE);
plot(weather1a$Date, weather1a$PrecipSum, main=" ", xlab="", ylab="", axes=FALSE, ylim=c(0, ceiling(max(weather1a$cummPrecip))+15), col="grey", type="h", lend=1, lwd=35, yaxs="i");
#plot(weather1a$Date, weather1a$PrecipSum, main=" ", xlab="", ylab="", axes=FALSE, ylim=c(0, ceiling(max(weather1a$cummPrecip))), col="grey50", type="h", lend=1, lwd=35, yaxs="i");
#lines(weather1a$Date, weather1a$PrecipSum-0.75, type="h", lwd=35-2, col="grey", lend=1);
par(new=TRUE);
plot(weather1a$Date, weather1a$DryTempH, pch=24,  xlab="", ylab="", ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), axes=FALSE, type="p", bg="white", col="black", yaxs="i");
par(new=TRUE);
plot(weather1a$Date, weather1a$DryTempA, pch=22,  xlab="", ylab="", ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), axes=FALSE, type="p", bg="grey50", col="black", yaxs="i");
par(new=TRUE);
plot(weather1a$Date, weather1a$DryTempL, pch=21, ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), xaxt='n', yaxt='n', ann=FALSE, type="p", bg="white", col="black", yaxs="i");
axis(1, weather1a$Date, format(weather1a$Date, "%b"), cex.axis=0.7, col="black", col.axis="black");
legend("bottom", c(expression("RF"[mtot]), expression("RF"[mmc]), "Min T", "Max T", "Avg T", expression("Avg T"[ann])), fill=c("grey", NA, NA, NA, NA, NA), border=c("grey", NA, NA, NA, NA, NA), box.cex=c(2.5, 1.0), pt.cex=c(1, 1), pch=c(NA, NA, 21, 24, 22, NA), pt.bg=c(NA, NA, "white", "white", "grey50", NA), col=c(NA, "black", "black", "black", "black", "black"), lty=c(NA, 1, NA, NA, NA, 2), lwd=c(NA, 2, NA, NA, NA, 1), inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="o", cex=0.6, x.intersp=0.1);
dev.off(); 


weather1a <- wo3[which(wo3$Station == "71142"), ];
str(weather1a);

##ceiling(max(weather1a[!is.na(weather1a$cummPrecip), ]$cummPrecip))+15;
##smoothingSpline <- smooth.spline(weather1a[!is.na(weather1a$cummPrecip), ]$Date, weather1a[!is.na(weather1a$cummPrecip), ]$cummPrecip, spar=0.25);

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farm_weather_71142_FarmRFAvgJan2015toDec2018", ".pdf", sep=''), bg="transparent", width=8, height=6, pointsize=12, family="FreeSans");
par(mar=c(5, 5, 4, 6) + 0.1);
plot(weather1a$Date, weather1a$cummPrecip, pch=2, axes=FALSE, ylim=c(0, ceiling(max(weather1a[!is.na(weather1a$cummPrecip), ]$cummPrecip))+15), xlab="", ylab="", type="l", main=" ", lwd=2, col=NA, yaxs="i");
smoothingSpline <- smooth.spline(weather1a[!is.na(weather1a$cummPrecip), ]$Date, weather1a[!is.na(weather1a$cummPrecip), ]$cummPrecip, spar=0.25);
lines(smoothingSpline, lwd=2);
axis(2, ylim=c(0, 1), col="black", las=1, at=c(seq(from=0,to=1000,by=25)));  ## las=1 makes horizontal labels;
mtext(expression(paste("Monthly total (RF" ["mtot"], ") and mean monthly", sep="")), side=2, line=3.5);
mtext(expression(paste("cumulative (RF" ["mmc"], ") rainfall (mm)", sep="")), side=2, line=2.5);
box();
par(new=TRUE);
plot(weather1a$Date, weather1a$AnnTemp, lty=2,  xlab="", ylab="", ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), axes=FALSE, type="l", col=NA, yaxs="i");
abline(h=unique(weather1a$AnnTemp), col="black", lty=2);
mtext("Monthly min, max, average temperature and", side=4, col="black", line=2.5); 
mtext(expression(paste("average annual (Avg T" ["ann"], ") temperature (°C)", sep="")), side=4, col="black", line=3.5); 
axis(4, ylim=c(floor(min(weather1a$DryTempL)), ceiling(max(weather1a$DryTempH))), at=c(seq(from=-50,to=50,by=10)), col.axis="black", las=1, bg="white", col="black");
mtext("Months", side=1, col="black", line=2.5);  
par(new=TRUE);
plot(weather1a$Date, weather1a$PrecipSum, main=" ", xlab="", ylab="", axes=FALSE, ylim=c(0, ceiling(max(weather1a[!is.na(weather1a$cummPrecip), ]$cummPrecip))+15), col="grey", type="h", lend=1, lwd=35, yaxs="i");
#plot(weather1a$Date, weather1a$PrecipSum, main=" ", xlab="", ylab="", axes=FALSE, ylim=c(0, ceiling(max(weather1a$cummPrecip))), col="grey50", type="h", lend=1, lwd=35, yaxs="i");
#lines(weather1a$Date, weather1a$PrecipSum-0.75, type="h", lwd=35-2, col="grey", lend=1);
par(new=TRUE);
plot(weather1a$Date, weather1a$DryTempH, pch=24,  xlab="", ylab="", ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), axes=FALSE, type="p", bg="white", col="black", yaxs="i");
par(new=TRUE);
plot(weather1a$Date, weather1a$DryTempA, pch=22,  xlab="", ylab="", ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), axes=FALSE, type="p", bg="grey50", col="black", yaxs="i");
par(new=TRUE);
plot(weather1a$Date, weather1a$DryTempL, pch=21, ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), xaxt='n', yaxt='n', ann=FALSE, type="p", bg="white", col="black", yaxs="i");
axis(1, weather1a$Date, format(weather1a$Date, "%b"), cex.axis=0.7, col="black", col.axis="black");
legend("bottom", c(expression("RF"[mtot]), expression("RF"[mmc]), "Min T", "Max T", "Avg T", expression("Avg T"[ann])), fill=c("grey", NA, NA, NA, NA, NA), border=c("grey", NA, NA, NA, NA, NA), box.cex=c(2.5, 1.0), pt.cex=c(1, 1), pch=c(NA, NA, 21, 24, 22, NA), pt.bg=c(NA, NA, "white", "white", "grey50", NA), col=c(NA, "black", "black", "black", "black", "black"), lty=c(NA, 1, NA, NA, NA, 2), lwd=c(NA, 2, NA, NA, NA, 1), inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="o", cex=0.6, x.intersp=0.1);
dev.off(); 



weather1a <- wo4[which(wo4$Station == "71142"), ];
str(weather1a);


pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farm_weather_71142_FarmRFJan2015toDec2018", ".pdf", sep=''), bg="transparent", width=8, height=6, pointsize=12, family="FreeSans");
par(mar=c(4, 5, 2, 5) + 0.0);
plot(weather1a$Date, weather1a$cummPrecip, pch=2, axes=FALSE, ylim=c(0, ceiling(max(weather1a[!is.na(weather1a$cummPrecip), ]$cummPrecip))+15), xlab="", ylab="", type="l", main=" ", lwd=2, col=NA, yaxs="i");
#smoothingSpline <- smooth.spline(weather1a[!is.na(weather1a$cummPrecip), ]$Date, weather1a[!is.na(weather1a$cummPrecip), ]$cummPrecip, spar=0.25);
#lines(smoothingSpline, lwd=2);
smoothingSpline1 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2015), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2015), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline2 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2016), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2016), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline3 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2017), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2017), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline4 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2018), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2018), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
lines(smoothingSpline1, lwd=2);
lines(smoothingSpline2, lwd=2);
lines(smoothingSpline3, lwd=2);
lines(smoothingSpline4, lwd=2);
axis(2, ylim=c(0, 1), col="black", las=1, at=c(seq(from=0,to=1000,by=25)));  ## las=1 makes horizontal labels;
mtext(expression(paste("Monthly total (RF" ["mtot"], ") and mean monthly", sep="")), side=2, line=3.5);
mtext(expression(paste("cumulative (RF" ["mmc"], ") rainfall (mm)", sep="")), side=2, line=2.5);
box();
par(new=TRUE);
plot(weather1a$Date, weather1a$AnnTemp, lty=2,  xlab="", ylab="", ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), axes=FALSE, type="l", col=NA, yaxs="i");
abline(h=unique(weather1a$AnnTemp), col="black", lty=2);
mtext("Monthly min, max, average temperature and", side=4, col="black", line=2.5); 
mtext(expression(paste("average annual (Avg T" ["ann"], ") temperature (°C)", sep="")), side=4, col="black", line=3.5); 
axis(4, ylim=c(floor(min(weather1a$DryTempL)), ceiling(max(weather1a$DryTempH))), at=c(seq(from=-50,to=50,by=10)), col.axis="black", las=1, bg="white", col="black");
mtext("Months", side=1, col="black", line=2.5);  
par(new=TRUE);
plot(weather1a$Date, weather1a$PrecipSum, main=" ", xlab="", ylab="", axes=FALSE, ylim=c(0, ceiling(max(weather1a[!is.na(weather1a$cummPrecip), ]$cummPrecip))+15), col="grey", type="h", lend=1, lwd=10, yaxs="i");
#plot(weather1a$Date, weather1a$PrecipSum, main=" ", xlab="", ylab="", axes=FALSE, ylim=c(0, ceiling(max(weather1a$cummPrecip))), col="grey50", type="h", lend=1, lwd=35, yaxs="i");
#lines(weather1a$Date, weather1a$PrecipSum-0.75, type="h", lwd=35-2, col="grey", lend=1);
par(new=TRUE);
plot(weather1a$Date, weather1a$DryTempH, pch=24,  xlab="", ylab="", ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), axes=FALSE, type="p", bg="white", col="black", yaxs="i");
par(new=TRUE);
plot(weather1a$Date, weather1a$DryTempA, pch=22,  xlab="", ylab="", ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), axes=FALSE, type="p", bg="grey50", col="black", yaxs="i");
par(new=TRUE);
plot(weather1a$Date, weather1a$DryTempL, pch=21, ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), xaxt='n', yaxt='n', ann=FALSE, type="p", bg="white", col="black", yaxs="i");
at <- weather1a$Date;
at1 <- at[rep(c(FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE), 4)];
at2 <- at[rep(c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE), 4)];
at3 <- at[rep(c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE), 4)];
axis(side=1, at=at1, format(at1, "%b\n"), labels=FALSE, tck=-0.01, cex.axis=0.7, col="black", col.axis="black"); 
axis(side=1, at=at2, format(at2, "%b\n"), tck=-0.02, cex.axis=0.75, col="black", col.axis="black"); 
axis(side=1, at=at3, format(at3, "%b\n%Y"), tck=-0.02, cex.axis=0.75, col="black", col.axis="black");
#axis(1, weather1a$Date, format(weather1a$Date, "%b"), cex.axis=0.7, col="black", col.axis="black");
legend("bottom", c(expression("RF"[mtot]), expression("RF"[mmc]), "Min T", "Max T", "Avg T", expression("Avg T"[ann])), fill=c("grey", NA, NA, NA, NA, NA), border=c("grey", NA, NA, NA, NA, NA), box.cex=c(1.2, 0.8), pt.cex=c(1, 1), pch=c(NA, NA, 21, 24, 22, NA), pt.bg=c(NA, NA, "white", "white", "grey50", NA), col=c(NA, "black", "black", "black", "black", "black"), lty=c(NA, 1, NA, NA, NA, 2), lwd=c(NA, 2, NA, NA, NA, 1), inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="o", cex=0.75, x.intersp=0.1);
dev.off();   