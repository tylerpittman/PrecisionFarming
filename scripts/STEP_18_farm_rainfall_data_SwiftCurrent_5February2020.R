## Swift Current weather data, and farm rainfall data from Jan 1986 to Dec 2018;
## Tyler Pittman, 5 February 2020
# Rscript --no-save /Users/tylerpittman/Farm/STEP_18_farm_rainfall_data_SwiftCurrent_5February2020.R

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
output$totalPrecip;
output$date;
output$Station;

wo0 <- output;
wo0$Station <- as.factor(wo0$Station);
wo0$Month <- months(wo0$date);
wo0$Year <- format(wo0$date, format="%Y");
wo0$DryTempH <- wo0$maxTemp;
wo0$DryTempL <- wo0$minTemp;
wo0$DryTempA <- wo0$meanTemp;
wo0$PrecipSum <- wo0$totalPrecip;
str(wo0);


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

### Examine to find warmest July based on DryTempA ###;
#wo4[which(wo4$Month %in% c("July")), ];



#---------------#---------------#---------------#---------------#---------------#;    


weather1a <- wo4;
str(weather1a);


pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farm_weather_71142_FarmRFJan1986toDec2018", ".pdf", sep=''), bg="transparent", width=50, height=6, pointsize=12, family="FreeSans");
par(mar=c(4, 5, 2, 5) + 0.0);
plot(weather1a$Date, weather1a$cummPrecip, pch=2, axes=FALSE, ylim=c(0, ceiling(max(weather1a[!is.na(weather1a$cummPrecip), ]$cummPrecip))+15), xlab="", ylab="", type="l", main=" ", lwd=2, col=NA, yaxs="i");
smoothingSpline1 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1986), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1986), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline2 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1987), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1987), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline3 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1988), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1988), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline4 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1989), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1989), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline5 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1990), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1990), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline6 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1991), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1991), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline7 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1992), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1992), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline8 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1993), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1993), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline9 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1994), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1994), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline10 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1995), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1995), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline11 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1996), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1996), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline12 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1997), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1997), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline13 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1998), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1998), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline14 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1999), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1999), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline15 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2000), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2000), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline16 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2001), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2001), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline17 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2002), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2002), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline18 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2003), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2003), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline19 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2004), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2004), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline20 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2005), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2005), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline21 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2006), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2006), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline22 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2007), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2007), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline23 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2008), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2008), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline24 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2009), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2009), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline25 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2010), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2010), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline26 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2011), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2011), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline27 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2012), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2012), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline28 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2013), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2013), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline29 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2014), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2014), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline30 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2015), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2015), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline31 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2016), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2016), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline32 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2017), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2017), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline33 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2018), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2018), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
lines(smoothingSpline1, lwd=2);
lines(smoothingSpline2, lwd=2);
lines(smoothingSpline3, lwd=2);
lines(smoothingSpline4, lwd=2);
lines(smoothingSpline5, lwd=2);
lines(smoothingSpline6, lwd=2);
lines(smoothingSpline7, lwd=2);
lines(smoothingSpline8, lwd=2);
lines(smoothingSpline9, lwd=2);
lines(smoothingSpline10, lwd=2);
lines(smoothingSpline11, lwd=2);
lines(smoothingSpline12, lwd=2);
lines(smoothingSpline13, lwd=2);
lines(smoothingSpline14, lwd=2);
lines(smoothingSpline15, lwd=2);
lines(smoothingSpline16, lwd=2);
lines(smoothingSpline17, lwd=2);
lines(smoothingSpline18, lwd=2);
lines(smoothingSpline19, lwd=2);
lines(smoothingSpline20, lwd=2);
lines(smoothingSpline21, lwd=2);
lines(smoothingSpline22, lwd=2);
lines(smoothingSpline23, lwd=2);
lines(smoothingSpline24, lwd=2);
lines(smoothingSpline25, lwd=2);
lines(smoothingSpline26, lwd=2);
lines(smoothingSpline27, lwd=2);
lines(smoothingSpline28, lwd=2);
lines(smoothingSpline29, lwd=2);
lines(smoothingSpline30, lwd=2);
lines(smoothingSpline31, lwd=2);
lines(smoothingSpline32, lwd=2);
lines(smoothingSpline33, lwd=2);
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
par(new=TRUE);
plot(weather1a$Date, weather1a$DryTempH, pch=24,  xlab="", ylab="", ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), axes=FALSE, type="p", bg="white", col="black", yaxs="i");
par(new=TRUE);
plot(weather1a$Date, weather1a$DryTempA, pch=22,  xlab="", ylab="", ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), axes=FALSE, type="p", bg="grey50", col="black", yaxs="i");
par(new=TRUE);
plot(weather1a$Date, weather1a$DryTempL, pch=21, ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), xaxt='n', yaxt='n', ann=FALSE, type="p", bg="white", col="black", yaxs="i");
at <- weather1a$Date;
#length(at);
at1 <- at[rep(c(FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE), 33)];
at2 <- at[rep(c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE), 33)];
at3 <- at[rep(c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE), 33)];
axis(side=1, at=at1, format(at1, "%b\n"), labels=FALSE, tck=-0.01, cex.axis=0.7, col="black", col.axis="black"); 
axis(side=1, at=at2, format(at2, "%b\n"), tck=-0.02, cex.axis=0.75, col="black", col.axis="black"); 
axis(side=1, at=at3, format(at3, "%b\n%Y"), tck=-0.02, cex.axis=0.75, col="black", col.axis="black");
#axis(1, weather1a$Date, format(weather1a$Date, "%b"), cex.axis=0.7, col="black", col.axis="black");
legend("bottom", c(expression("RF"[mtot]), expression("RF"[mmc]), "Min T", "Max T", "Avg T", expression("Avg T"[ann])), fill=c("grey", NA, NA, NA, NA, NA), border=c("grey", NA, NA, NA, NA, NA), box.cex=c(1.2, 0.8), pt.cex=c(1, 1), pch=c(NA, NA, 21, 24, 22, NA), pt.bg=c(NA, NA, "white", "white", "grey50", NA), col=c(NA, "black", "black", "black", "black", "black"), lty=c(NA, 1, NA, NA, NA, 2), lwd=c(NA, 2, NA, NA, NA, 1), inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="o", cex=0.75, x.intersp=0.1);
dev.off();   

#----------------------------------------------------#;

weather1a <- wo4[which(wo4$Year %in% c(1986:1996)), ];
annTemp <- mean(weather1a$DryTempA, na.rm=TRUE);
annTemp <- as.data.frame(annTemp);
colnames(annTemp) <- c("AnnTemp");

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farm_weather_1_FarmRFJan1986toDec1996", ".pdf", sep=''), bg="transparent", width=11, height=6, pointsize=12, family="FreeSans");
par(mar=c(4, 5, 2, 5) + 0.0);
plot(weather1a$Date, weather1a$cummPrecip, pch=2, axes=FALSE, ylim=c(0, 525), xlab="", ylab="", type="l", main=" ", lwd=2, col=NA, yaxs="i");
smoothingSpline1 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1986), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1986), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline2 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1987), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1987), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline3 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1988), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1988), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline4 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1989), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1989), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline5 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1990), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1990), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline6 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1991), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1991), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline7 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1992), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1992), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline8 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1993), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1993), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline9 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1994), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1994), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline10 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1995), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1995), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline11 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1996), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1996), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
lines(smoothingSpline1, lwd=2);
lines(smoothingSpline2, lwd=2);
lines(smoothingSpline3, lwd=2);
lines(smoothingSpline4, lwd=2);
lines(smoothingSpline5, lwd=2);
lines(smoothingSpline6, lwd=2);
lines(smoothingSpline7, lwd=2);
lines(smoothingSpline8, lwd=2);
lines(smoothingSpline9, lwd=2);
lines(smoothingSpline10, lwd=2);
lines(smoothingSpline11, lwd=2);
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
plot(weather1a$Date, weather1a$PrecipSum, main=" ", xlab="", ylab="", axes=FALSE, ylim=c(0, ceiling(max(weather1a[!is.na(weather1a$cummPrecip), ]$cummPrecip))+15), col="grey", type="h", lend=1, lwd=5, yaxs="i");
par(new=TRUE);
plot(weather1a$Date, weather1a$DryTempH, pch=24,  xlab="", ylab="", ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), axes=FALSE, type="p", bg="white", col="black", yaxs="i");
par(new=TRUE);
plot(weather1a$Date, weather1a$DryTempA, pch=22,  xlab="", ylab="", ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), axes=FALSE, type="p", bg="grey50", col="black", yaxs="i");
par(new=TRUE);
plot(weather1a$Date, weather1a$DryTempL, pch=21, ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), xaxt='n', yaxt='n', ann=FALSE, type="p", bg="white", col="black", yaxs="i");
at <- weather1a$Date;
#length(at);
at1 <- at[rep(c(FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE), 11)];
at2 <- at[rep(c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE), 11)];
at3 <- at[rep(c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE), 11)];
axis(side=1, at=at1, paste(substring(format(at1, "%b\n"), 1, 1), "\n", sep=""), labels=FALSE, tck=-0.01, cex.axis=0.7, col="black", col.axis="black"); 
axis(side=1, at=at2, paste(substring(format(at2, "%b\n"), 1, 1), "\n", sep=""), tck=-0.02, cex.axis=0.75, col="black", col.axis="black"); 
axis(side=1, at=at3, paste(substring(format(at3, "%b\n%Y"), 1, 1), "\n", format(at3, "%Y"), sep=""), tck=-0.02, cex.axis=0.75, col="black", col.axis="black");
#axis(1, weather1a$Date, format(weather1a$Date, "%b"), cex.axis=0.7, col="black", col.axis="black");
legend("bottom", c(expression("RF"[mtot]), expression("RF"[mmc]), "Min T", "Max T", "Avg T", expression("Avg T"[ann])), fill=c("grey", NA, NA, NA, NA, NA), border=c("grey", NA, NA, NA, NA, NA), box.cex=c(1.2, 0.8), pt.cex=c(1, 1), pch=c(NA, NA, 21, 24, 22, NA), pt.bg=c(NA, NA, "white", "white", "grey50", NA), col=c(NA, "black", "black", "black", "black", "black"), lty=c(NA, 1, NA, NA, NA, 2), lwd=c(NA, 2, NA, NA, NA, 1), inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="o", cex=0.75, x.intersp=0.1);
dev.off();  


weather1a <- wo4[which(wo4$Year %in% c(1997:2007)), ];
annTemp <- mean(weather1a$DryTempA, na.rm=TRUE);
annTemp <- as.data.frame(annTemp);
colnames(annTemp) <- c("AnnTemp");

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farm_weather_2_FarmRFJan1997toDec2007", ".pdf", sep=''), bg="transparent", width=11, height=6, pointsize=12, family="FreeSans");
par(mar=c(4, 5, 2, 5) + 0.0);
plot(weather1a$Date, weather1a$cummPrecip, pch=2, axes=FALSE, ylim=c(0, 525), xlab="", ylab="", type="l", main=" ", lwd=2, col=NA, yaxs="i");
smoothingSpline12 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1997), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1997), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline13 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1998), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1998), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline14 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1999), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1999), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline15 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2000), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2000), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline16 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2001), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2001), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline17 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2002), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2002), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline18 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2003), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2003), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline19 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2004), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2004), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline20 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2005), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2005), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline21 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2006), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2006), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline22 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2007), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2007), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
lines(smoothingSpline12, lwd=2);
lines(smoothingSpline13, lwd=2);
lines(smoothingSpline14, lwd=2);
lines(smoothingSpline15, lwd=2);
lines(smoothingSpline16, lwd=2);
lines(smoothingSpline17, lwd=2);
lines(smoothingSpline18, lwd=2);
lines(smoothingSpline19, lwd=2);
lines(smoothingSpline20, lwd=2);
lines(smoothingSpline21, lwd=2);
lines(smoothingSpline22, lwd=2);
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
plot(weather1a$Date, weather1a$PrecipSum, main=" ", xlab="", ylab="", axes=FALSE, ylim=c(0, ceiling(max(weather1a[!is.na(weather1a$cummPrecip), ]$cummPrecip))+15), col="grey", type="h", lend=1, lwd=5, yaxs="i");
par(new=TRUE);
plot(weather1a$Date, weather1a$DryTempH, pch=24,  xlab="", ylab="", ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), axes=FALSE, type="p", bg="white", col="black", yaxs="i");
par(new=TRUE);
plot(weather1a$Date, weather1a$DryTempA, pch=22,  xlab="", ylab="", ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), axes=FALSE, type="p", bg="grey50", col="black", yaxs="i");
par(new=TRUE);
plot(weather1a$Date, weather1a$DryTempL, pch=21, ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), xaxt='n', yaxt='n', ann=FALSE, type="p", bg="white", col="black", yaxs="i");
at <- weather1a$Date;
#length(at);
at1 <- at[rep(c(FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE), 11)];
at2 <- at[rep(c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE), 11)];
at3 <- at[rep(c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE), 11)];
axis(side=1, at=at1, paste(substring(format(at1, "%b\n"), 1, 1), "\n", sep=""), labels=FALSE, tck=-0.01, cex.axis=0.7, col="black", col.axis="black"); 
axis(side=1, at=at2, paste(substring(format(at2, "%b\n"), 1, 1), "\n", sep=""), tck=-0.02, cex.axis=0.75, col="black", col.axis="black"); 
axis(side=1, at=at3, paste(substring(format(at3, "%b\n%Y"), 1, 1), "\n", format(at3, "%Y"), sep=""), tck=-0.02, cex.axis=0.75, col="black", col.axis="black");
#axis(1, weather1a$Date, format(weather1a$Date, "%b"), cex.axis=0.7, col="black", col.axis="black");
legend("bottom", c(expression("RF"[mtot]), expression("RF"[mmc]), "Min T", "Max T", "Avg T", expression("Avg T"[ann])), fill=c("grey", NA, NA, NA, NA, NA), border=c("grey", NA, NA, NA, NA, NA), box.cex=c(1.2, 0.8), pt.cex=c(1, 1), pch=c(NA, NA, 21, 24, 22, NA), pt.bg=c(NA, NA, "white", "white", "grey50", NA), col=c(NA, "black", "black", "black", "black", "black"), lty=c(NA, 1, NA, NA, NA, 2), lwd=c(NA, 2, NA, NA, NA, 1), inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="o", cex=0.75, x.intersp=0.1);
dev.off(); 


weather1a <- wo4[which(wo4$Year %in% c(2008:2018)), ];
annTemp <- mean(weather1a$DryTempA, na.rm=TRUE);
annTemp <- as.data.frame(annTemp);
colnames(annTemp) <- c("AnnTemp");

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farm_weather_3_FarmRFJan2008toDec2018", ".pdf", sep=''), bg="transparent", width=11, height=6, pointsize=12, family="FreeSans");
par(mar=c(4, 5, 2, 5) + 0.0);
plot(weather1a$Date, weather1a$cummPrecip, pch=2, axes=FALSE, ylim=c(0, 525), xlab="", ylab="", type="l", main=" ", lwd=2, col=NA, yaxs="i");
smoothingSpline23 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2008), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2008), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline24 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2009), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2009), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline25 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2010), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2010), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline26 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2011), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2011), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline27 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2012), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2012), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline28 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2013), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2013), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline29 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2014), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2014), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline30 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2015), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2015), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline31 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2016), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2016), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline32 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2017), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2017), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline33 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2018), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 2018), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
lines(smoothingSpline23, lwd=2);
lines(smoothingSpline24, lwd=2);
lines(smoothingSpline25, lwd=2);
lines(smoothingSpline26, lwd=2);
lines(smoothingSpline27, lwd=2);
lines(smoothingSpline28, lwd=2);
lines(smoothingSpline29, lwd=2);
lines(smoothingSpline30, lwd=2);
lines(smoothingSpline31, lwd=2);
lines(smoothingSpline32, lwd=2);
lines(smoothingSpline33, lwd=2);
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
plot(weather1a$Date, weather1a$PrecipSum, main=" ", xlab="", ylab="", axes=FALSE, ylim=c(0, ceiling(max(weather1a[!is.na(weather1a$cummPrecip), ]$cummPrecip))+15), col="grey", type="h", lend=1, lwd=5, yaxs="i");
par(new=TRUE);
plot(weather1a$Date, weather1a$DryTempH, pch=24,  xlab="", ylab="", ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), axes=FALSE, type="p", bg="white", col="black", yaxs="i");
par(new=TRUE);
plot(weather1a$Date, weather1a$DryTempA, pch=22,  xlab="", ylab="", ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), axes=FALSE, type="p", bg="grey50", col="black", yaxs="i");
par(new=TRUE);
plot(weather1a$Date, weather1a$DryTempL, pch=21, ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), xaxt='n', yaxt='n', ann=FALSE, type="p", bg="white", col="black", yaxs="i");
at <- weather1a$Date;
#length(at);
at1 <- at[rep(c(FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE), 11)];
at2 <- at[rep(c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE), 11)];
at3 <- at[rep(c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE), 11)];
axis(side=1, at=at1, paste(substring(format(at1, "%b\n"), 1, 1), "\n", sep=""), labels=FALSE, tck=-0.01, cex.axis=0.7, col="black", col.axis="black"); 
axis(side=1, at=at2, paste(substring(format(at2, "%b\n"), 1, 1), "\n", sep=""), tck=-0.02, cex.axis=0.75, col="black", col.axis="black"); 
axis(side=1, at=at3, paste(substring(format(at3, "%b\n%Y"), 1, 1), "\n", format(at3, "%Y"), sep=""), tck=-0.02, cex.axis=0.75, col="black", col.axis="black");
#axis(1, weather1a$Date, format(weather1a$Date, "%b"), cex.axis=0.7, col="black", col.axis="black");
legend("bottom", c(expression("RF"[mtot]), expression("RF"[mmc]), "Min T", "Max T", "Avg T", expression("Avg T"[ann])), fill=c("grey", NA, NA, NA, NA, NA), border=c("grey", NA, NA, NA, NA, NA), box.cex=c(1.2, 0.8), pt.cex=c(1, 1), pch=c(NA, NA, 21, 24, 22, NA), pt.bg=c(NA, NA, "white", "white", "grey50", NA), col=c(NA, "black", "black", "black", "black", "black"), lty=c(NA, 1, NA, NA, NA, 2), lwd=c(NA, 2, NA, NA, NA, 1), inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="o", cex=0.75, x.intersp=0.1);
dev.off(); 



#-----------------------------------------!!!!!-----------------------------------------#;
#-----------------------------------------!!!!!-----------------------------------------#;
#-----------------------------------------!!!!!-----------------------------------------#;
#-----------------------------------------!!!!!-----------------------------------------#;
annTemp <- mean(wo4$DryTempA, na.rm=TRUE);
annTemp <- as.data.frame(annTemp);
colnames(annTemp) <- c("AnnTemp");
weather1a <- wo4[which(wo4$Year %in% c(1986:1996)), ];
weather2a <- wo4[which(wo4$Year %in% c(1997:2007)), ];
weather3a <- wo4[which(wo4$Year %in% c(2008:2018)), ];


#width=7, height=9;
pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farm_weather_BigPlot_FarmRFJan1986toDec1996", ".pdf", sep=''), bg="transparent", width=7, height=9, pointsize=12, family="FreeSans");
par(mfrow=c(3,1), mar=c(0, 0, 2, 0) + 0.0, oma=c(3.6, 3.8, 0, 3.5));
plot(weather1a$Date, weather1a$cummPrecip, pch=2, axes=FALSE, ylim=c(0, 525), xlab="", ylab="", type="l", main=" ", lwd=2, col=NA, yaxs="i");
smoothingSpline1 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1986), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1986), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline2 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1987), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1987), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline3 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1988), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1988), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline4 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1989), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1989), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline5 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1990), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1990), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline6 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1991), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1991), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline7 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1992), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1992), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline8 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1993), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1993), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline9 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1994), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1994), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline10 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1995), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1995), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline11 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1996), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1996), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
lines(smoothingSpline1, lwd=2);
lines(smoothingSpline2, lwd=2);
lines(smoothingSpline3, lwd=2);
lines(smoothingSpline4, lwd=2);
lines(smoothingSpline5, lwd=2);
lines(smoothingSpline6, lwd=2);
lines(smoothingSpline7, lwd=2);
lines(smoothingSpline8, lwd=2);
lines(smoothingSpline9, lwd=2);
lines(smoothingSpline10, lwd=2);
lines(smoothingSpline11, lwd=2);
axis(2, ylim=c(0, 1), col="black", las=1, at=c(seq(from=0,to=1000,by=25)));  ## las=1 makes horizontal labels;
#mtext(expression(paste("Monthly total (RF" ["mtot"], ") and mean monthly", sep="")), side=2, line=3.5);
#mtext(expression(paste("cumulative (RF" ["mmc"], ") rainfall (mm)", sep="")), side=2, line=2.5);
box();
par(new=TRUE);
plot(weather1a$Date, weather1a$AnnTemp, lty=2,  xlab="", ylab="", ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), axes=FALSE, type="l", col=NA, yaxs="i");
abline(h=unique(weather1a$AnnTemp), col="black", lty=2);
#mtext("Monthly min, max, average temperature and", side=4, col="black", line=2.5); 
#mtext(expression(paste("average annual (Avg T" ["ann"], ") temperature (°C)", sep="")), side=4, col="black", line=3.5); 
axis(4, ylim=c(floor(min(weather1a$DryTempL)), ceiling(max(weather1a$DryTempH))), at=c(seq(from=-50,to=50,by=10)), col.axis="black", las=1, bg="white", col="black");
#mtext("Months", side=1, col="black", line=2.5);  
par(new=TRUE);
plot(weather1a$Date, weather1a$PrecipSum, main=" ", xlab="", ylab="", axes=FALSE, ylim=c(0, ceiling(max(weather1a[!is.na(weather1a$cummPrecip), ]$cummPrecip))+15), col="grey", type="h", lend=1, lwd=3.5, yaxs="i");
par(new=TRUE);
plot(weather1a$Date, weather1a$DryTempH, pch=24,  xlab="", ylab="", ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), axes=FALSE, type="p", bg="white", col="black", yaxs="i");
par(new=TRUE);
plot(weather1a$Date, weather1a$DryTempA, pch=22,  xlab="", ylab="", ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), axes=FALSE, type="p", bg="grey50", col="black", yaxs="i");
par(new=TRUE);
plot(weather1a$Date, weather1a$DryTempL, pch=21, ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), xaxt='n', yaxt='n', ann=FALSE, type="p", bg="white", col="black", yaxs="i");
at <- weather1a$Date;
#length(at);
at1 <- at[rep(c(FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE), 11)];
at2 <- at[rep(c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE), 11)];
at3 <- at[rep(c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE), 11)];
axis(side=1, at=at1, paste(substring(format(at1, "%b\n"), 1, 1), "\n", sep=""), labels=FALSE, tck=-0.005, cex.axis=0.81, col="black", col.axis="black"); 
axis(side=1, at=at2, paste(substring(format(at2, "%b\n"), 1, 1), "\n", sep=""), tck=-0.015, cex.axis=0.81, col="black", col.axis="black"); 
axis(side=1, at=at3, paste(substring(format(at3, "%b\n%Y"), 1, 1), "\n", format(at3, "%Y"), sep=""), tck=-0.015, cex.axis=0.81, col="black", col.axis="black");
legend("bottom", c(expression("RF"[mtot]), expression("RF"[mmc]), "Min T", "Max T", "Avg T", expression("Avg T"[ann])), fill=c("grey", NA, NA, NA, NA, NA), border=c("grey", NA, NA, NA, NA, NA), box.cex=c(1.2, 0.8), pt.cex=c(1, 1), pch=c(NA, NA, 21, 24, 22, NA), pt.bg=c(NA, NA, "white", "white", "grey50", NA), col=c(NA, "black", "black", "black", "black", "black"), lty=c(NA, 1, NA, NA, NA, 2), lwd=c(NA, 2, NA, NA, NA, 1), inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="o", cex=1.0, x.intersp=0.1);


plot(weather2a$Date, weather2a$cummPrecip, pch=2, axes=FALSE, ylim=c(0, 525), xlab="", ylab="", type="l", main=" ", lwd=2, col=NA, yaxs="i");
smoothingSpline12 <- cobs(as.numeric(weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 1997), ]$Date), weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 1997), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline13 <- cobs(as.numeric(weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 1998), ]$Date), weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 1998), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline14 <- cobs(as.numeric(weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 1999), ]$Date), weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 1999), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline15 <- cobs(as.numeric(weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 2000), ]$Date), weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 2000), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline16 <- cobs(as.numeric(weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 2001), ]$Date), weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 2001), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline17 <- cobs(as.numeric(weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 2002), ]$Date), weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 2002), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline18 <- cobs(as.numeric(weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 2003), ]$Date), weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 2003), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline19 <- cobs(as.numeric(weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 2004), ]$Date), weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 2004), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline20 <- cobs(as.numeric(weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 2005), ]$Date), weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 2005), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline21 <- cobs(as.numeric(weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 2006), ]$Date), weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 2006), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline22 <- cobs(as.numeric(weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 2007), ]$Date), weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 2007), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
lines(smoothingSpline12, lwd=2);
lines(smoothingSpline13, lwd=2);
lines(smoothingSpline14, lwd=2);
lines(smoothingSpline15, lwd=2);
lines(smoothingSpline16, lwd=2);
lines(smoothingSpline17, lwd=2);
lines(smoothingSpline18, lwd=2);
lines(smoothingSpline19, lwd=2);
lines(smoothingSpline20, lwd=2);
lines(smoothingSpline21, lwd=2);
lines(smoothingSpline22, lwd=2);
axis(2, ylim=c(0, 1), col="black", las=1, at=c(seq(from=0,to=1000,by=25)));  ## las=1 makes horizontal labels;
#mtext(expression(paste("Monthly total (RF" ["mtot"], ") and mean monthly", sep="")), side=2, line=3.5);
#mtext(expression(paste("cumulative (RF" ["mmc"], ") rainfall (mm)", sep="")), side=2, line=2.5);
#box();
par(new=TRUE);
plot(weather2a$Date, weather2a$AnnTemp, lty=2,  xlab="", ylab="", ylim=c(floor(min(weather2a$DryTempL))-4, ceiling(max(weather2a$DryTempH))+2), axes=FALSE, type="l", col=NA, yaxs="i");
abline(h=unique(weather2a$AnnTemp), col="black", lty=2);
#mtext("Monthly min, max, average temperature and", side=4, col="black", line=2.5); 
#mtext(expression(paste("average annual (Avg T" ["ann"], ") temperature (°C)", sep="")), side=4, col="black", line=3.5); 
axis(4, ylim=c(floor(min(weather2a$DryTempL)), ceiling(max(weather2a$DryTempH))), at=c(seq(from=-50,to=50,by=10)), col.axis="black", las=1, bg="white", col="black");
#mtext("Months", side=1, col="black", line=2.5);  
par(new=TRUE);
plot(weather2a$Date, weather2a$PrecipSum, main=" ", xlab="", ylab="", axes=FALSE, ylim=c(0, ceiling(max(weather2a[!is.na(weather2a$cummPrecip), ]$cummPrecip))+15), col="grey", type="h", lend=1, lwd=3.5, yaxs="i");
par(new=TRUE);
plot(weather2a$Date, weather2a$DryTempH, pch=24,  xlab="", ylab="", ylim=c(floor(min(weather2a$DryTempL))-4, ceiling(max(weather2a$DryTempH))+2), axes=FALSE, type="p", bg="white", col="black", yaxs="i");
par(new=TRUE);
plot(weather2a$Date, weather2a$DryTempA, pch=22,  xlab="", ylab="", ylim=c(floor(min(weather2a$DryTempL))-4, ceiling(max(weather2a$DryTempH))+2), axes=FALSE, type="p", bg="grey50", col="black", yaxs="i");
par(new=TRUE);
plot(weather2a$Date, weather2a$DryTempL, pch=21, ylim=c(floor(min(weather2a$DryTempL))-4, ceiling(max(weather2a$DryTempH))+2), xaxt='n', yaxt='n', ann=FALSE, type="p", bg="white", col="black", yaxs="i");
at <- weather2a$Date;
#length(at);
at1 <- at[rep(c(FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE), 11)];
at2 <- at[rep(c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE), 11)];
at3 <- at[rep(c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE), 11)];
axis(side=1, at=at1, paste(substring(format(at1, "%b\n"), 1, 1), "\n", sep=""), labels=FALSE, tck=-0.005, cex.axis=0.81, col="black", col.axis="black"); 
axis(side=1, at=at2, paste(substring(format(at2, "%b\n"), 1, 1), "\n", sep=""), tck=-0.015, cex.axis=0.81, col="black", col.axis="black"); 
axis(side=1, at=at3, paste(substring(format(at3, "%b\n%Y"), 1, 1), "\n", format(at3, "%Y"), sep=""), tck=-0.015, cex.axis=0.81, col="black", col.axis="black");


plot(weather3a$Date, weather3a$cummPrecip, pch=2, axes=FALSE, ylim=c(0, 525), xlab="", ylab="", type="l", main=" ", lwd=2, col=NA, yaxs="i");
smoothingSpline23 <- cobs(as.numeric(weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2008), ]$Date), weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2008), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline24 <- cobs(as.numeric(weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2009), ]$Date), weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2009), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline25 <- cobs(as.numeric(weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2010), ]$Date), weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2010), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline26 <- cobs(as.numeric(weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2011), ]$Date), weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2011), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline27 <- cobs(as.numeric(weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2012), ]$Date), weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2012), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline28 <- cobs(as.numeric(weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2013), ]$Date), weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2013), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline29 <- cobs(as.numeric(weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2014), ]$Date), weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2014), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline30 <- cobs(as.numeric(weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2015), ]$Date), weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2015), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline31 <- cobs(as.numeric(weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2016), ]$Date), weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2016), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline32 <- cobs(as.numeric(weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2017), ]$Date), weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2017), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline33 <- cobs(as.numeric(weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2018), ]$Date), weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2018), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
lines(smoothingSpline23, lwd=2);
lines(smoothingSpline24, lwd=2);
lines(smoothingSpline25, lwd=2);
lines(smoothingSpline26, lwd=2);
lines(smoothingSpline27, lwd=2);
lines(smoothingSpline28, lwd=2);
lines(smoothingSpline29, lwd=2);
lines(smoothingSpline30, lwd=2);
lines(smoothingSpline31, lwd=2);
lines(smoothingSpline32, lwd=2);
lines(smoothingSpline33, lwd=2);
axis(2, ylim=c(0, 1), col="black", las=1, at=c(seq(from=0,to=1000,by=25)));  ## las=1 makes horizontal labels;
#mtext(expression(paste("Monthly total (RF" ["mtot"], ") and mean monthly", sep="")), side=2, line=3.5);
#mtext(expression(paste("cumulative (RF" ["mmc"], ") rainfall (mm)", sep="")), side=2, line=2.5);
#box();
par(new=TRUE);
plot(weather3a$Date, weather3a$AnnTemp, lty=2,  xlab="", ylab="", ylim=c(floor(min(weather3a$DryTempL))-4, ceiling(max(weather3a$DryTempH))+2), axes=FALSE, type="l", col=NA, yaxs="i");
abline(h=unique(weather3a$AnnTemp), col="black", lty=2);
#mtext("Monthly min, max, average temperature and", side=4, col="black", line=2.5); 
#mtext(expression(paste("average annual (Avg T" ["ann"], ") temperature (°C)", sep="")), side=4, col="black", line=3.5); 
axis(4, ylim=c(floor(min(weather3a$DryTempL)), ceiling(max(weather3a$DryTempH))), at=c(seq(from=-50,to=50,by=10)), col.axis="black", las=1, bg="white", col="black");
par(new=TRUE);
plot(weather3a$Date, weather3a$PrecipSum, main=" ", xlab="", ylab="", axes=FALSE, ylim=c(0, ceiling(max(weather3a[!is.na(weather3a$cummPrecip), ]$cummPrecip))+15), col="grey", type="h", lend=1, lwd=3.5, yaxs="i");
par(new=TRUE);
plot(weather3a$Date, weather3a$DryTempH, pch=24,  xlab="", ylab="", ylim=c(floor(min(weather3a$DryTempL))-4, ceiling(max(weather3a$DryTempH))+2), axes=FALSE, type="p", bg="white", col="black", yaxs="i");
par(new=TRUE);
plot(weather3a$Date, weather3a$DryTempA, pch=22,  xlab="", ylab="", ylim=c(floor(min(weather3a$DryTempL))-4, ceiling(max(weather3a$DryTempH))+2), axes=FALSE, type="p", bg="grey50", col="black", yaxs="i");
par(new=TRUE);
plot(weather3a$Date, weather3a$DryTempL, pch=21, ylim=c(floor(min(weather3a$DryTempL))-4, ceiling(max(weather3a$DryTempH))+2), xaxt='n', yaxt='n', ann=FALSE, type="p", bg="white", col="black", yaxs="i");
at <- weather3a$Date;
#length(at);
at1 <- at[rep(c(FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE), 11)];
at2 <- at[rep(c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE), 11)];
at3 <- at[rep(c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE), 11)];
axis(side=1, at=at1, paste(substring(format(at1, "%b\n"), 1, 1), "\n", sep=""), labels=FALSE, tck=-0.005, cex.axis=0.81, col="black", col.axis="black"); 
axis(side=1, at=at2, paste(substring(format(at2, "%b\n"), 1, 1), "\n", sep=""), tck=-0.015, cex.axis=0.81, col="black", col.axis="black"); 
axis(side=1, at=at3, paste(substring(format(at3, "%b\n%Y"), 1, 1), "\n", format(at3, "%Y"), sep=""), tck=-0.015, cex.axis=0.81, col="black", col.axis="black");
#axis(1, weather3a$Date, format(weather3a$Date, "%b"), cex.axis=0.7, col="black", col.axis="black");
#legend("bottom", c(expression("RF"[mtot]), expression("RF"[mmc]), "Min T", "Max T", "Avg T", expression("Avg T"[ann])), fill=c("grey", NA, NA, NA, NA, NA), border=c("grey", NA, NA, NA, NA, NA), box.cex=c(1.2, 0.8), pt.cex=c(1, 1), pch=c(NA, NA, 21, 24, 22, NA), pt.bg=c(NA, NA, "white", "white", "grey50", NA), col=c(NA, "black", "black", "black", "black", "black"), lty=c(NA, 1, NA, NA, NA, 2), lwd=c(NA, 2, NA, NA, NA, 1), inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="o", cex=0.75, x.intersp=0.1);
mtext(expression(paste("Monthly min, max, average temperature and average annual (Avg T" ["ann"], ") temperature (°C)", sep="")), outer=TRUE, side=4, col="black", cex=0.9, line=2.7);
mtext(expression(paste("Monthly total (RF" ["mtot"], ") and mean monthly cumulative (RF" ["mmc"], ") rainfall (mm)", sep="")), outer=TRUE, side=2, cex=0.9, line=2.5);
mtext("Month by Year", side=1, col="black", cex=0.9, line=2.5);  
dev.off(); 

#-----------------------------------------!!!!!-----------------------------------------#;
#-----------------------------------------!!!!!-----------------------------------------#;
#-----------------------------------------!!!!!-----------------------------------------#;
#-----------------------------------------!!!!!-----------------------------------------#;


#-----------------------------------------!!!!!-----------------------------------------#;
#-----------------------------------------!!!!!-----------------------------------------#;
#-----------------------------------------!!!!!-----------------------------------------#;
#-----------------------------------------!!!!!-----------------------------------------#;

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

gs <- cbind(as.numeric(growing_season2$count), as.numeric(growing_season2$Year));
gs <- as.data.frame(gs);
colnames(gs) <- c("FF_Duration", "Year");

wo4_tmp <- join_all(list(wo4, gs), by = c('Year'), type = "left"); 
wo4 <- wo4_tmp;
#weatherFarm_check <- weatherFarm[which(weatherFarm$Year == 1987), ]; #137 days;
###---------------------------------##;


annTemp <- mean(wo4$DryTempA, na.rm=TRUE);
annTemp <- as.data.frame(annTemp);
colnames(annTemp) <- c("AnnTemp");
weather1a <- wo4[which(wo4$Year %in% c(1986:1996)), ];
weather2a <- wo4[which(wo4$Year %in% c(1997:2007)), ];
weather3a <- wo4[which(wo4$Year %in% c(2008:2018)), ];


#width=7, height=9;
pdf.options(encoding='CP1250');
cairo_pdf(file=paste("farm_weather_BigPlot_FF_FarmRFJan1986toDec1996", ".pdf", sep=''), bg="transparent", width=7, height=9, pointsize=12, family="FreeSans");
par(mfrow=c(3,1), mar=c(0, 0, 2, 0) + 0.0, oma=c(3.6, 3.8, 0, 3.5));
plot(weather1a$Date, weather1a$cummPrecip, pch=2, axes=FALSE, ylim=c(0, 525), xlab="", ylab="", type="l", main=" ", lwd=2, col=NA, yaxs="i");
smoothingSpline1 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1986), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1986), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline2 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1987), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1987), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline3 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1988), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1988), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline4 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1989), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1989), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline5 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1990), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1990), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline6 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1991), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1991), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline7 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1992), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1992), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline8 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1993), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1993), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline9 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1994), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1994), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline10 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1995), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1995), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline11 <- cobs(as.numeric(weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1996), ]$Date), weather1a[which(weather1a$cummPrecip >= 0 & weather1a$Year == 1996), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
lines(smoothingSpline1, lwd=2);
lines(smoothingSpline2, lwd=2);
lines(smoothingSpline3, lwd=2);
lines(smoothingSpline4, lwd=2);
lines(smoothingSpline5, lwd=2);
lines(smoothingSpline6, lwd=2);
lines(smoothingSpline7, lwd=2);
lines(smoothingSpline8, lwd=2);
lines(smoothingSpline9, lwd=2);
lines(smoothingSpline10, lwd=2);
lines(smoothingSpline11, lwd=2);
axis(2, ylim=c(0, 1), col="black", las=1, at=c(seq(from=0,to=1000,by=25)));  ## las=1 makes horizontal labels;
#mtext(expression(paste("Monthly total (RF" ["mtot"], ") and mean monthly", sep="")), side=2, line=3.5);
#mtext(expression(paste("cumulative (RF" ["mmc"], ") rainfall (mm)", sep="")), side=2, line=2.5);
box();
par(new=TRUE);
plot(weather1a$Date, weather1a$AnnTemp, lty=2,  xlab="", ylab="", ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), axes=FALSE, type="l", col=NA, yaxs="i");
abline(h=unique(weather1a$AnnTemp), col="black", lty=2);
#mtext("Monthly min, max, average temperature and", side=4, col="black", line=2.5); 
#mtext(expression(paste("average annual (Avg T" ["ann"], ") temperature (°C)", sep="")), side=4, col="black", line=3.5); 
axis(4, ylim=c(floor(min(weather1a$DryTempL)), ceiling(max(weather1a$DryTempH))), at=c(seq(from=-50,to=50,by=10)), col.axis="black", las=1, bg="white", col="black");
#mtext("Months", side=1, col="black", line=2.5);  
par(new=TRUE);
plot(weather1a$Date, weather1a$PrecipSum, main=" ", xlab="", ylab="", axes=FALSE, ylim=c(0, ceiling(max(weather1a[!is.na(weather1a$cummPrecip), ]$cummPrecip))+15), col="grey", type="h", lend=1, lwd=3.5, yaxs="i");
par(new=TRUE);
plot(weather1a$Date, weather1a$DryTempH, pch=24,  xlab="", ylab="", ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), axes=FALSE, type="p", bg="white", col="black", yaxs="i");
par(new=TRUE);
plot(weather1a$Date, weather1a$DryTempA, pch=22,  xlab="", ylab="", ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), axes=FALSE, type="p", bg="grey50", col="black", yaxs="i");
par(new=TRUE);
plot(weather1a$Date, weather1a$DryTempL, pch=21, ylim=c(floor(min(weather1a$DryTempL))-4, ceiling(max(weather1a$DryTempH))+2), xaxt='n', yaxt='n', ann=FALSE, type="p", bg="white", col="black", yaxs="i");
at <- weather1a$Date;
ff <- weather1a$FF_Duration;
#length(at);
at1 <- at[rep(c(FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE), 11)];
at2 <- at[rep(c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE), 11)];
at3 <- at[rep(c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE), 11)];
ff3 <- ff[rep(c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE), 11)];
axis(side=1, at=at1, paste(substring(format(at1, "%b\n"), 1, 1), "\n", sep=""), labels=FALSE, tck=-0.005, cex.axis=0.81, col="black", col.axis="black"); 
axis(side=1, at=at2, paste(substring(format(at2, "%b\n"), 1, 1), "\n", sep=""), tck=-0.015, cex.axis=0.81, col="black", col.axis="black"); 
axis(side=1, at=at3, paste(substring(format(at3, "%b\n%Y"), 1, 1), "\n", format(at3, "%Y"), ": ", ff3, "", sep=""), tck=-0.015, cex.axis=0.81, col="black", col.axis="black");
legend("bottom", c(expression("RF"[mtot]), expression("RF"[mmc]), "Min T", "Max T", "Avg T", expression("Avg T"[ann])), fill=c("grey", NA, NA, NA, NA, NA), border=c("grey", NA, NA, NA, NA, NA), box.cex=c(1.2, 0.8), pt.cex=c(1, 1), pch=c(NA, NA, 21, 24, 22, NA), pt.bg=c(NA, NA, "white", "white", "grey50", NA), col=c(NA, "black", "black", "black", "black", "black"), lty=c(NA, 1, NA, NA, NA, 2), lwd=c(NA, 2, NA, NA, NA, 1), inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="o", cex=1.0, x.intersp=0.1);


plot(weather2a$Date, weather2a$cummPrecip, pch=2, axes=FALSE, ylim=c(0, 525), xlab="", ylab="", type="l", main=" ", lwd=2, col=NA, yaxs="i");
smoothingSpline12 <- cobs(as.numeric(weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 1997), ]$Date), weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 1997), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline13 <- cobs(as.numeric(weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 1998), ]$Date), weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 1998), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline14 <- cobs(as.numeric(weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 1999), ]$Date), weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 1999), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline15 <- cobs(as.numeric(weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 2000), ]$Date), weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 2000), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline16 <- cobs(as.numeric(weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 2001), ]$Date), weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 2001), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline17 <- cobs(as.numeric(weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 2002), ]$Date), weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 2002), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline18 <- cobs(as.numeric(weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 2003), ]$Date), weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 2003), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline19 <- cobs(as.numeric(weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 2004), ]$Date), weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 2004), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline20 <- cobs(as.numeric(weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 2005), ]$Date), weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 2005), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline21 <- cobs(as.numeric(weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 2006), ]$Date), weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 2006), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline22 <- cobs(as.numeric(weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 2007), ]$Date), weather2a[which(weather2a$cummPrecip >= 0 & weather2a$Year == 2007), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
lines(smoothingSpline12, lwd=2);
lines(smoothingSpline13, lwd=2);
lines(smoothingSpline14, lwd=2);
lines(smoothingSpline15, lwd=2);
lines(smoothingSpline16, lwd=2);
lines(smoothingSpline17, lwd=2);
lines(smoothingSpline18, lwd=2);
lines(smoothingSpline19, lwd=2);
lines(smoothingSpline20, lwd=2);
lines(smoothingSpline21, lwd=2);
lines(smoothingSpline22, lwd=2);
axis(2, ylim=c(0, 1), col="black", las=1, at=c(seq(from=0,to=1000,by=25)));  ## las=1 makes horizontal labels;
#mtext(expression(paste("Monthly total (RF" ["mtot"], ") and mean monthly", sep="")), side=2, line=3.5);
#mtext(expression(paste("cumulative (RF" ["mmc"], ") rainfall (mm)", sep="")), side=2, line=2.5);
#box();
par(new=TRUE);
plot(weather2a$Date, weather2a$AnnTemp, lty=2,  xlab="", ylab="", ylim=c(floor(min(weather2a$DryTempL))-4, ceiling(max(weather2a$DryTempH))+2), axes=FALSE, type="l", col=NA, yaxs="i");
abline(h=unique(weather2a$AnnTemp), col="black", lty=2);
#mtext("Monthly min, max, average temperature and", side=4, col="black", line=2.5); 
#mtext(expression(paste("average annual (Avg T" ["ann"], ") temperature (°C)", sep="")), side=4, col="black", line=3.5); 
axis(4, ylim=c(floor(min(weather2a$DryTempL)), ceiling(max(weather2a$DryTempH))), at=c(seq(from=-50,to=50,by=10)), col.axis="black", las=1, bg="white", col="black");
#mtext("Months", side=1, col="black", line=2.5);  
par(new=TRUE);
plot(weather2a$Date, weather2a$PrecipSum, main=" ", xlab="", ylab="", axes=FALSE, ylim=c(0, ceiling(max(weather2a[!is.na(weather2a$cummPrecip), ]$cummPrecip))+15), col="grey", type="h", lend=1, lwd=3.5, yaxs="i");
par(new=TRUE);
plot(weather2a$Date, weather2a$DryTempH, pch=24,  xlab="", ylab="", ylim=c(floor(min(weather2a$DryTempL))-4, ceiling(max(weather2a$DryTempH))+2), axes=FALSE, type="p", bg="white", col="black", yaxs="i");
par(new=TRUE);
plot(weather2a$Date, weather2a$DryTempA, pch=22,  xlab="", ylab="", ylim=c(floor(min(weather2a$DryTempL))-4, ceiling(max(weather2a$DryTempH))+2), axes=FALSE, type="p", bg="grey50", col="black", yaxs="i");
par(new=TRUE);
plot(weather2a$Date, weather2a$DryTempL, pch=21, ylim=c(floor(min(weather2a$DryTempL))-4, ceiling(max(weather2a$DryTempH))+2), xaxt='n', yaxt='n', ann=FALSE, type="p", bg="white", col="black", yaxs="i");
at <- weather2a$Date;
ff <- weather2a$FF_Duration;
#length(at);
at1 <- at[rep(c(FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE), 11)];
at2 <- at[rep(c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE), 11)];
at3 <- at[rep(c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE), 11)];
ff3 <- ff[rep(c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE), 11)];
axis(side=1, at=at1, paste(substring(format(at1, "%b\n"), 1, 1), "\n", sep=""), labels=FALSE, tck=-0.005, cex.axis=0.81, col="black", col.axis="black"); 
axis(side=1, at=at2, paste(substring(format(at2, "%b\n"), 1, 1), "\n", sep=""), tck=-0.015, cex.axis=0.81, col="black", col.axis="black"); 
axis(side=1, at=at3, paste(substring(format(at3, "%b\n%Y"), 1, 1), "\n", format(at3, "%Y"), ": ", ff3, "", sep=""), tck=-0.015, cex.axis=0.81, col="black", col.axis="black");


plot(weather3a$Date, weather3a$cummPrecip, pch=2, axes=FALSE, ylim=c(0, 525), xlab="", ylab="", type="l", main=" ", lwd=2, col=NA, yaxs="i");
smoothingSpline23 <- cobs(as.numeric(weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2008), ]$Date), weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2008), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline24 <- cobs(as.numeric(weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2009), ]$Date), weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2009), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline25 <- cobs(as.numeric(weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2010), ]$Date), weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2010), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline26 <- cobs(as.numeric(weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2011), ]$Date), weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2011), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline27 <- cobs(as.numeric(weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2012), ]$Date), weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2012), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline28 <- cobs(as.numeric(weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2013), ]$Date), weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2013), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline29 <- cobs(as.numeric(weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2014), ]$Date), weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2014), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline30 <- cobs(as.numeric(weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2015), ]$Date), weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2015), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline31 <- cobs(as.numeric(weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2016), ]$Date), weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2016), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline32 <- cobs(as.numeric(weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2017), ]$Date), weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2017), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
smoothingSpline33 <- cobs(as.numeric(weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2018), ]$Date), weather3a[which(weather3a$cummPrecip >= 0 & weather3a$Year == 2018), ]$cummPrecip, constraint="increase", lambda=0, degree=1, tau=0.5);
lines(smoothingSpline23, lwd=2);
lines(smoothingSpline24, lwd=2);
lines(smoothingSpline25, lwd=2);
lines(smoothingSpline26, lwd=2);
lines(smoothingSpline27, lwd=2);
lines(smoothingSpline28, lwd=2);
lines(smoothingSpline29, lwd=2);
lines(smoothingSpline30, lwd=2);
lines(smoothingSpline31, lwd=2);
lines(smoothingSpline32, lwd=2);
lines(smoothingSpline33, lwd=2);
axis(2, ylim=c(0, 1), col="black", las=1, at=c(seq(from=0,to=1000,by=25)));  ## las=1 makes horizontal labels;
#mtext(expression(paste("Monthly total (RF" ["mtot"], ") and mean monthly", sep="")), side=2, line=3.5);
#mtext(expression(paste("cumulative (RF" ["mmc"], ") rainfall (mm)", sep="")), side=2, line=2.5);
#box();
par(new=TRUE);
plot(weather3a$Date, weather3a$AnnTemp, lty=2,  xlab="", ylab="", ylim=c(floor(min(weather3a$DryTempL))-4, ceiling(max(weather3a$DryTempH))+2), axes=FALSE, type="l", col=NA, yaxs="i");
abline(h=unique(weather3a$AnnTemp), col="black", lty=2);
#mtext("Monthly min, max, average temperature and", side=4, col="black", line=2.5); 
#mtext(expression(paste("average annual (Avg T" ["ann"], ") temperature (°C)", sep="")), side=4, col="black", line=3.5); 
axis(4, ylim=c(floor(min(weather3a$DryTempL)), ceiling(max(weather3a$DryTempH))), at=c(seq(from=-50,to=50,by=10)), col.axis="black", las=1, bg="white", col="black");
par(new=TRUE);
plot(weather3a$Date, weather3a$PrecipSum, main=" ", xlab="", ylab="", axes=FALSE, ylim=c(0, ceiling(max(weather3a[!is.na(weather3a$cummPrecip), ]$cummPrecip))+15), col="grey", type="h", lend=1, lwd=3.5, yaxs="i");
par(new=TRUE);
plot(weather3a$Date, weather3a$DryTempH, pch=24,  xlab="", ylab="", ylim=c(floor(min(weather3a$DryTempL))-4, ceiling(max(weather3a$DryTempH))+2), axes=FALSE, type="p", bg="white", col="black", yaxs="i");
par(new=TRUE);
plot(weather3a$Date, weather3a$DryTempA, pch=22,  xlab="", ylab="", ylim=c(floor(min(weather3a$DryTempL))-4, ceiling(max(weather3a$DryTempH))+2), axes=FALSE, type="p", bg="grey50", col="black", yaxs="i");
par(new=TRUE);
plot(weather3a$Date, weather3a$DryTempL, pch=21, ylim=c(floor(min(weather3a$DryTempL))-4, ceiling(max(weather3a$DryTempH))+2), xaxt='n', yaxt='n', ann=FALSE, type="p", bg="white", col="black", yaxs="i");
at <- weather3a$Date;
ff <- weather3a$FF_Duration;
#length(at);
at1 <- at[rep(c(FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE), 11)];
at2 <- at[rep(c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE), 11)];
at3 <- at[rep(c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE), 11)];
ff3 <- ff[rep(c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE), 11)];
axis(side=1, at=at1, paste(substring(format(at1, "%b\n"), 1, 1), "\n", sep=""), labels=FALSE, tck=-0.005, cex.axis=0.81, col="black", col.axis="black"); 
axis(side=1, at=at2, paste(substring(format(at2, "%b\n"), 1, 1), "\n", sep=""), tck=-0.015, cex.axis=0.81, col="black", col.axis="black"); 
axis(side=1, at=at3, paste(substring(format(at3, "%b\n%Y"), 1, 1), "\n", format(at3, "%Y"), ": ", ff3, "", sep=""), tck=-0.015, cex.axis=0.81, col="black", col.axis="black");
#axis(1, weather3a$Date, format(weather3a$Date, "%b"), cex.axis=0.7, col="black", col.axis="black");
#legend("bottom", c(expression("RF"[mtot]), expression("RF"[mmc]), "Min T", "Max T", "Avg T", expression("Avg T"[ann])), fill=c("grey", NA, NA, NA, NA, NA), border=c("grey", NA, NA, NA, NA, NA), box.cex=c(1.2, 0.8), pt.cex=c(1, 1), pch=c(NA, NA, 21, 24, 22, NA), pt.bg=c(NA, NA, "white", "white", "grey50", NA), col=c(NA, "black", "black", "black", "black", "black"), lty=c(NA, 1, NA, NA, NA, 2), lwd=c(NA, 2, NA, NA, NA, 1), inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="o", cex=0.75, x.intersp=0.1);
mtext(expression(paste("Monthly min, max, average temperature and average annual (Avg T" ["ann"], ") temperature (°C)", sep="")), outer=TRUE, side=4, col="black", cex=0.9, line=2.7);
mtext(expression(paste("Monthly total (RF" ["mtot"], ") and mean monthly cumulative (RF" ["mmc"], ") rainfall (mm)", sep="")), outer=TRUE, side=2, cex=0.9, line=2.5);
mtext("Month by Year: frost-free days", side=1, col="black", cex=0.9, line=2.5);  
dev.off(); 

#-----------------------------------------!!!!!-----------------------------------------#;
#-----------------------------------------!!!!!-----------------------------------------#;
#-----------------------------------------!!!!!-----------------------------------------#;
#-----------------------------------------!!!!!-----------------------------------------#;
