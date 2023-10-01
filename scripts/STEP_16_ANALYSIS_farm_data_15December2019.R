## Analysis of combine yield data and seeding data from 2015 to 2018
## Adds in record of elemental sulfur added to fields for 2016 and 2017 growing years (optional)
## Tyler Pittman, 15 December 2019
# Rscript --no-save /Users/tylerpittman/Farm/STEP_16_ANALYSIS_farm_data_15December2019.R

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
#source("shape2poly.R"); # Reads in shape2poly function and others
#source("polygonizer.R"); # Reads in polygonizer function

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
), ncol=6, byrow=T);
key <- as.data.frame(key);
colnames(key) <- c("field", "Crop2014", "Crop2015", "Crop2016", "Crop2017", "Crop2018");
key$Crop2014 <- as.character(levels(key$Crop2014))[key$Crop2014];
key$Crop2015 <- as.character(levels(key$Crop2015))[key$Crop2015];
key$Crop2016 <- as.character(levels(key$Crop2016))[key$Crop2016];
key$Crop2017 <- as.character(levels(key$Crop2017))[key$Crop2017];
key$Crop2018 <- as.character(levels(key$Crop2018))[key$Crop2018];


#---------------------------------------------------------------------------------------#;


setwd("/Users/tylerpittman/Farm/agYieldProject/data");
load(paste("farmDataComplete2015to2018", ".RData", sep=""));
str(farmData2015to2018c); #4,268,459 obs.;
levels(farmData2015to2018c$Field);
levels(farmData2015to2018c$Product);
which(is.na(farmData2015to2018c$Product)); #All yield points have a product;
summary(farmData2015to2018c$Mat_Days);
summary(farmData2015to2018c$Elevation);


#---------------------------------------------------------------------------------------#;
#---------------------------------------------------------------------------------------#;
### This part was discovered in Dad's field records ###;
## Adds in record of elemental sulfur added to fields for 2016 and 2017 growing years (11 October 2019);
## (2016) 80 lbs/acre sulphur spread previous October to fields: 1, 3, 7, 9, 12, 14, 17;
## (2017) 80 lbs/acre sulphur spread previous October to fields: 6, 8, 15, 16, 18, 21;
### 80 lbs/acre = 89.6681 kg/ha, but this conversion is done in a section below this ###;
### elemental sulphur is 0-0-0-90 to 99, assume 90% sulphur or 72 lbs/acre;

#farmData2015to2018c[which(farmData2015to2018c$ProcYear == 2016 & farmData2015to2018c$Field == 900001), ]$S_Rate = farmData2015to2018c[which(farmData2015to2018c$ProcYear == 2016 & farmData2015to2018c$Field == 900001), ]$S_Rate + 72; 
#farmData2015to2018c[which(farmData2015to2018c$ProcYear == 2016 & farmData2015to2018c$Field == 900003), ]$S_Rate = farmData2015to2018c[which(farmData2015to2018c$ProcYear == 2016 & farmData2015to2018c$Field == 900003), ]$S_Rate + 72; 
#farmData2015to2018c[which(farmData2015to2018c$ProcYear == 2016 & farmData2015to2018c$Field == 900007), ]$S_Rate = farmData2015to2018c[which(farmData2015to2018c$ProcYear == 2016 & farmData2015to2018c$Field == 900007), ]$S_Rate + 72; 
#farmData2015to2018c[which(farmData2015to2018c$ProcYear == 2016 & farmData2015to2018c$Field == 900009), ]$S_Rate = farmData2015to2018c[which(farmData2015to2018c$ProcYear == 2016 & farmData2015to2018c$Field == 900009), ]$S_Rate + 72; 
#farmData2015to2018c[which(farmData2015to2018c$ProcYear == 2016 & farmData2015to2018c$Field == 900012), ]$S_Rate = farmData2015to2018c[which(farmData2015to2018c$ProcYear == 2016 & farmData2015to2018c$Field == 900012), ]$S_Rate + 72; 
#farmData2015to2018c[which(farmData2015to2018c$ProcYear == 2016 & farmData2015to2018c$Field == 900014), ]$S_Rate = farmData2015to2018c[which(farmData2015to2018c$ProcYear == 2016 & farmData2015to2018c$Field == 900014), ]$S_Rate + 72; 
#farmData2015to2018c[which(farmData2015to2018c$ProcYear == 2016 & farmData2015to2018c$Field == 900017), ]$S_Rate = farmData2015to2018c[which(farmData2015to2018c$ProcYear == 2016 & farmData2015to2018c$Field == 900017), ]$S_Rate + 72; 
#
#farmData2015to2018c[which(farmData2015to2018c$ProcYear == 2017 & farmData2015to2018c$Field == 900006), ]$S_Rate = farmData2015to2018c[which(farmData2015to2018c$ProcYear == 2017 & farmData2015to2018c$Field == 900006), ]$S_Rate + 72; 
#farmData2015to2018c[which(farmData2015to2018c$ProcYear == 2017 & farmData2015to2018c$Field == 900008), ]$S_Rate = farmData2015to2018c[which(farmData2015to2018c$ProcYear == 2017 & farmData2015to2018c$Field == 900008), ]$S_Rate + 72; 
#farmData2015to2018c[which(farmData2015to2018c$ProcYear == 2017 & farmData2015to2018c$Field == 900015), ]$S_Rate = farmData2015to2018c[which(farmData2015to2018c$ProcYear == 2017 & farmData2015to2018c$Field == 900015), ]$S_Rate + 72; 
#farmData2015to2018c[which(farmData2015to2018c$ProcYear == 2017 & farmData2015to2018c$Field == 900016), ]$S_Rate = farmData2015to2018c[which(farmData2015to2018c$ProcYear == 2017 & farmData2015to2018c$Field == 900016), ]$S_Rate + 72; 
#farmData2015to2018c[which(farmData2015to2018c$ProcYear == 2017 & farmData2015to2018c$Field == 900018), ]$S_Rate = farmData2015to2018c[which(farmData2015to2018c$ProcYear == 2017 & farmData2015to2018c$Field == 900018), ]$S_Rate + 72; 
#farmData2015to2018c[which(farmData2015to2018c$ProcYear == 2017 & farmData2015to2018c$Field == 900021), ]$S_Rate = farmData2015to2018c[which(farmData2015to2018c$ProcYear == 2017 & farmData2015to2018c$Field == 900021), ]$S_Rate + 72; 

farmData2015to2018c$S_Rate <- jitter(farmData2015to2018c$S_Rate);
farmData2015to2018c[which(farmData2015to2018c$S_Rate <= 0), ]$S_Rate <- 0;
farmData2015to2018c[which(farmData2015to2018c$P_Rate <= 0), ]$P_Rate <- 0;
farmData2015to2018c[which(farmData2015to2018c$N_Rate <= 0), ]$N_Rate <- 0;
farmData2015to2018c[which(farmData2015to2018c$K_Rate <= 0), ]$K_Rate <- 0;
#---------------------------------------------------------------------------------------#;
#---------------------------------------------------------------------------------------#;


farmData2015to2018c[1,]$latitude; #5 digits;
farmData2015to2018c[1,]$longitude; #4 digits;


#-()-()-()-()-()-()-()-()-()-()-()-()-()-()-()-#;
### Determine the previous crop, link to 2014 information for observations in ProcYear 2015;
cropFY <- read_excel("cropField2014to2018.xlsx", sheet = "crop");
cropFY <- as.data.frame(cropFY);
str(cropFY);

farmData2015to2018c$F <- sub(".*0{2,}", "", farmData2015to2018c$Field);
farmData2015to2018c$F <- as.factor(farmData2015to2018c$F);
str(farmData2015to2018c);

farmData2015c <- farmData2015to2018c[which(farmData2015to2018c$ProcYear == 2015), ];
farmData2016c <- farmData2015to2018c[which(farmData2015to2018c$ProcYear == 2016), ];
farmData2017c <- farmData2015to2018c[which(farmData2015to2018c$ProcYear == 2017), ];
farmData2018c <- farmData2015to2018c[which(farmData2015to2018c$ProcYear == 2018), ];

cropFY2015 <- cropFY[,c(1,2)];
cropFY2016 <- cropFY[,c(1,3)];
cropFY2017 <- cropFY[,c(1,4)];
cropFY2018 <- cropFY[,c(1,5)];
colnames(cropFY2015) <- c("F", "PrevCrop");
colnames(cropFY2016) <- c("F", "PrevCrop");
colnames(cropFY2017) <- c("F", "PrevCrop");
colnames(cropFY2018) <- c("F", "PrevCrop");

tmp2015 <- join_all(list(farmData2015c, cropFY2015), by = c('F'), type = "left");
tmp2016 <- join_all(list(farmData2016c, cropFY2016), by = c('F'), type = "left");
tmp2017 <- join_all(list(farmData2017c, cropFY2017), by = c('F'), type = "left");
tmp2018 <- join_all(list(farmData2018c, cropFY2018), by = c('F'), type = "left");

tmp2016[which(tmp2016$F == "10" & tmp2016$longitude <= -108.21505), ]$PrevCrop = "Lentils";
tmp2016[which(tmp2016$F == "10" & tmp2016$longitude >= -108.21505), ]$PrevCrop = "Canola";
tmp2018[which(tmp2018$F == "14" & tmp2018$latitude <= 50.8233), ]$PrevCrop = "Lentils";
tmp2018[which(tmp2018$F == "14" & tmp2018$latitude >= 50.8233), ]$PrevCrop = "Barley";

farmData2015to2018c <- rbind(tmp2015, tmp2016, tmp2017, tmp2018);
farmData2015to2018c$PrevCrop <- as.factor(farmData2015to2018c$PrevCrop);
levels(farmData2015to2018c$PrevCrop);
#-()-()-()-()-()-()-()-()-()-()-()-()-()-()-()-#;


#-----------------------------------------------------------------------------------------------------------------------------#;
#-----&-----&-----&-----&-----&-----#;
#### Convert to metric units;
farmData2015to2018c$Seed_Rate <- farmData2015to2018c$Seed_Rate * 1.12085;
farmData2015to2018c$Ino_Rate <- farmData2015to2018c$Ino_Rate * 1.12085;
farmData2015to2018c$N_Rate <- farmData2015to2018c$N_Rate * 1.12085;
farmData2015to2018c$P_Rate <- farmData2015to2018c$P_Rate * 1.12085;
farmData2015to2018c$K_Rate <- farmData2015to2018c$K_Rate * 1.12085;
farmData2015to2018c$S_Rate <- farmData2015to2018c$S_Rate * 1.12085;

farmData2015to2018c[which(farmData2015to2018c$Product == "Wheat"), ]$DryYield <- farmData2015to2018c[which(farmData2015to2018c$Product == "Wheat"), ]$DryYield * 67.25;
farmData2015to2018c[which(farmData2015to2018c$Product == "Lentils"), ]$DryYield <- farmData2015to2018c[which(farmData2015to2018c$Product == "Lentils"), ]$DryYield * 67.25;
farmData2015to2018c[which(farmData2015to2018c$Product == "Barley"), ]$DryYield <- farmData2015to2018c[which(farmData2015to2018c$Product == "Barley"), ]$DryYield * 53.80;
farmData2015to2018c[which(farmData2015to2018c$Product == "Canola"), ]$DryYield <- farmData2015to2018c[which(farmData2015to2018c$Product == "Canola"), ]$DryYield * 56.04;
farmData2015to2018c[which(farmData2015to2018c$Product == "Canary"), ]$DryYield <- farmData2015to2018c[which(farmData2015to2018c$Product == "Canary"), ]$DryYield * 56.04;

farmData2015to2018c$Elevation <- farmData2015to2018c$Elevation * 0.3048;
#-----&-----&-----&-----&-----&-----#;

###############################################################;
#### This determine mean elevation for 20 quarters we own #####;

#unique(farmData2015to2018c$Field);
elevation <- farmData2015to2018c[which(farmData2015to2018c$Field %in% c(900001:900021)), ]$Elevation;
summary(elevation);
###############################################################;


#----------- Data cleaning section ----------------#;
str(farmData2015to2018c); #4,268,459 obs.;

#df1 = farmData2015to2018c %>%
#  group_by(Field, Product, ProcYear) %>%
#  filter(!(abs(DryYield - mean(DryYield)) > 3*sd(DryYield))) ;
#str(df1); #4,064,516 obs. removing outliers outside 3 standard deviations of mean ;
#summary(df1$DryYield); #mean is 34.39;
#  	# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#	#0.00   17.97   32.76   34.39   50.95  212.84 
#sd(df1$DryYield); #sd is 23.5546;

df1 = farmData2015to2018c %>%
  group_by(Field, Product, ProcYear) %>%
  filter(!(abs(DryYield - mean(DryYield)) > 2*sd(DryYield))) ;
str(df1); #4,064,516 obs. removing outliers outside 2 standard deviations of mean ;
summary(df1$DryYield); #mean is 35.53;
	#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   #0.00   20.16   33.84   35.53   51.71  167.04 
sd(df1$DryYield); #sd is 22.8374;


farmData2015to2018c <- df1;
farmData2015to2018c <- farmData2015to2018c[which(farmData2015to2018c$DryYield >= 5.6), ]; #remove points with DryYield below 5.6 kg/ha;
str(farmData2015to2018c); #3,627,023 obs.;
#--------------------------------------------------#;


setwd("/Users/tylerpittman/Farm/agYieldProject/data");
save(farmData2015to2018c, file=paste("farmDataCompleteMetric2015to2018", ".RData", sep=""), compress="xz");
#load(paste("farmDataCompleteMetric2015to2018", ".RData", sep=""));

str(farmData2015to2018c); #3,627,023 obs.;
levels(farmData2015to2018c$Field);
levels(farmData2015to2018c$Product);
which(is.na(farmData2015to2018c$Product)); #All yield points have a product;
summary(farmData2015to2018c$Mat_Days);

########## check for 900011 in 2015 #########;
#tank4_2015s <- tank4_2015[which(tank4_2015@data$SP_ID_1 == "900011"), ]; #this is correct for seed rate;
###projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0";
###proj4string(tank4_2015s) <- CRS(projection);
#yieldDataConstraint2015s <- yieldDataConstraint2015[which(yieldDataConstraint2015@data$Field == "900011"), ]; #this is correct for yield;
#plot(tank4_2015s);
#plot(yieldDataConstraint2015s, col="red", add=TRUE);
#yieldDataConstraint2015s@data$Seed_Rate4 <- over(yieldDataConstraint2015s, tank4_2015s)$Seed_Rate;

summary(farmData2015to2018c$DryYield); #mean is 2577.890 kg/ha;
sd(farmData2015to2018c$DryYield); #sd is 1343.525 kg/ha;


### Detect and remove outliers by group; #Field, Product and ProcYear;
###https://stackoverflow.com/questions/28687515/search-for-and-remove-outliers-from-a-dataframe-grouped-by-a-variable
#df1 = farmData2015to2018c %>%
#  group_by(Field, Product, ProcYear) %>%
#  filter(!(abs(DryYield - mean(DryYield)) > 2*sd(DryYield))) %>%
#  summarise_each(funs(mean), DryYield);


summary(farmData2015to2018c$Elevation); #mean is 649.9 meters;
sd(farmData2015to2018c$Elevation); #sd is 6.684401 meters;

summary(farmData2015to2018c$Product);
legend1_text <- levels(farmData2015to2018c$Product);
twoCol <- c("grey", "white");

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("DryYield_product_boxplot.pdf", sep=''), bg="transparent", width=11, height=6, pointsize=12, family="FreeSans");
par(family="FreeSans", mfrow=c(1,1), mar=c(4.5, 4.1, 0.0, 2.1), oma=c(1,0,4,0), xpd=TRUE);
boxplot(DryYield ~ Product, data=farmData2015to2018c, ylab="Dry Yield (kg/ha)", xaxt="n", xlab="", col=twoCol);
text(x=seq_along(legend1_text), y=par("usr")[3]-1, srt=45, adj=1, labels=legend1_text, xpd=TRUE, cex=0.9);
options(scipen=999);
mtext(xpd=NA, outer=TRUE, cex=1.5, line=0.5, text=paste("Boxplot of Dry Yield by \n Crop Product"));
#mtext(text=paste("For nursing homes with organizations as owners.", "\nProcessing date:", processingDate), cex=1.0, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
dev.off();

summary(farmData2015to2018c$ProcYear);
legend1_text <- levels(farmData2015to2018c$ProcYear);
twoCol <- c("grey", "white");

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("DryYield_year_boxplot.pdf", sep=''), bg="transparent", width=11, height=6, pointsize=12, family="FreeSans");
par(family="FreeSans", mfrow=c(1,1), mar=c(4.5, 4.1, 0.0, 2.1), oma=c(1,0,4,0), xpd=TRUE);
boxplot(DryYield ~ ProcYear, data=farmData2015to2018c, ylab="Dry Yield (kg/ha)", xaxt="n", xlab="", col=twoCol);
text(x=seq_along(legend1_text), y=par("usr")[3]-1, srt=45, adj=1, labels=legend1_text, xpd=TRUE, cex=0.9);
options(scipen=999);
mtext(xpd=NA, outer=TRUE, cex=1.5, line=0.5, text=paste("Boxplot of Dry Yield by \n Year"));
#mtext(text=paste("For nursing homes with organizations as owners.", "\nProcessing date:", processingDate), cex=1.0, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
dev.off();


#### Below works nicely for categorical variables, used for Excel input ####
#-------;
neat.table <- function(x, name){
  xx <- data.frame(x)
  names(xx) <- c("Value", "Count")
  xx$Fraction <- with(xx, Count/sum(Count))
  data.frame(Variable = name, xx)
}

x <- lapply(farmData2015to2018c[, c("Product", "Machine", "ProcYear", "Field", "Cultivar")], table)
do.call(rbind, lapply(seq_along(x), function(i)neat.table(x[i], names(x[i]))))

mydata2015 <- farmData2015to2018c[which(farmData2015to2018c$ProcYear == "2015"), ];
x <- lapply(mydata2015[, c("Product", "Machine", "ProcYear", "Field", "Cultivar")], table)
do.call(rbind, lapply(seq_along(x), function(i)neat.table(x[i], names(x[i]))))

mydata2016 <- farmData2015to2018c[which(farmData2015to2018c$ProcYear == "2016"), ];
x <- lapply(mydata2016[, c("Product", "Machine", "ProcYear", "Field", "Cultivar")], table)
do.call(rbind, lapply(seq_along(x), function(i)neat.table(x[i], names(x[i]))))

mydata2017 <- farmData2015to2018c[which(farmData2015to2018c$ProcYear == "2017"), ];
x <- lapply(mydata2017[, c("Product", "Machine", "ProcYear", "Field", "Cultivar")], table)
do.call(rbind, lapply(seq_along(x), function(i)neat.table(x[i], names(x[i]))))

mydata2018 <- farmData2015to2018c[which(farmData2015to2018c$ProcYear == "2018"), ];
x <- lapply(mydata2018[, c("Product", "Machine", "ProcYear", "Field", "Cultivar")], table)
do.call(rbind, lapply(seq_along(x), function(i)neat.table(x[i], names(x[i]))))

#####;

mydata2 <- farmData2015to2018c[, c("Elevation", "DryYield", "Seed_Rate", "Ino_Rate", "N_Rate", "P_Rate", "K_Rate", "S_Rate", "Mat_Days", "Product", "ProcYear")];
#colnames(mydata2);
str(mydata2);
#mydata2;

##This part takes 2 minutes to run;  
mydata2 %>%
    group_by(Product, ProcYear) %>%
    mutate(id = 1:n()) %>%
    ungroup() %>%
    gather(temp, val, Elevation, DryYield, Seed_Rate, Ino_Rate, N_Rate, P_Rate, K_Rate, S_Rate, Mat_Days) %>%
    unite(temp1, temp, Product, ProcYear, sep = '_') %>%
    spread(temp1, val) %>%
    select(-id) %>%
    as.data.frame() %>%
    stargazer(digits=2, type='text')


vars <- farmData2015to2018c[, c("Elevation", "DryYield", "Seed_Rate", "Ino_Rate", "N_Rate", "P_Rate", "K_Rate", "S_Rate", "Mat_Days")];
group <- interaction(as.factor(farmData2015to2018c$Product), as.factor(farmData2015to2018c$ProcYear)); #use this code for multiplying factor levels to use for grouping in tables;
##levels(group);
#result <- tableContinuous(vars=vars, group=group, prec=2, cap="Farming Characteristics by Crop and Year", lab="tab1_descr_stat", stats=c("mean", "s", "n"), longtable=TRUE);

### display default statistics, only use a subset of observations, grouped analysis
#result <- tableContinuous(vars=vars, group=group, prec=2, cap="Farming Characteristics by Crop and Year", lab="tab1_descr_stat", stats=c("mean", "s", "n"), longtable=TRUE);
##cat(gsub("\\\\hline\n[^\n]+& all &[^\n]+\n", "", result)); #get rid of all category;
##cat(gsub("\\hline", "", result)); #get rid of \hline;


#look at page 64 of Hadfield's Course Notes to specify random effects;
#us variance structure, of idh variance structure?;
citation();
citation("MCMCglmm");
citation("igraph");
packageVersion("MCMCglmm");
packageVersion("igraph");

#https://gkhajduk.github.io/2017-10-25-cleanMCMCglmm/
clean.MCMC <- function(x) {
    sols <- summary(x)$solutions  ## pull out relevant info from model summary
    Gcovs <- summary(x)$Gcovariances
    Rcovs <- summary(x)$Rcovariances
    fixed <- data.frame(row.names(sols), sols, row.names = NULL)  ## convert to dataframes with the row.names as the first col
    random <- data.frame(row.names(Gcovs), Gcovs, row.names = NULL)
    residual <- data.frame(row.names(Rcovs), Rcovs, row.names = NULL)
    names(fixed)[names(fixed) == "row.names.sols."] <- "variable"  ## change the columns names to variable, so they all match
    names(random)[names(random) == "row.names.Gcovs."] <- "variable"
    names(residual)[names(residual) == "row.names.Rcovs."] <- "variable"
    fixed$effect <- "fixed"  ## add ID column for type of effect (fixed, random, residual)
    random$effect <- "random"
    residual$effect <- "residual"
    modelTerms <- as.data.frame(bind_rows(fixed, random, residual))  # merge it all together
}

getName.MCMC <- function(x) deparse(substitute(x))  # add the model name

#---------------#;
##farmData2015to2018c$lonlat <- round(farmData2015to2018c$longitude, digits=3) * round(farmData2015to2018c$latitude, digits=4);
##farmData2015to2018c$lonlat <- paste(farmData2015to2018c$longitude, "c", farmData2015to2018c$latitude, sep="");
farmData2015to2018c$lonlat <- paste(round(farmData2015to2018c$longitude, digits=4), "c", round(farmData2015to2018c$latitude, digits=4), sep="");
#length(farmData2015to2018c$lonlat);
#length(unique(farmData2015to2018c$lonlat)); #217,898 unique lonlat points out of 3,586,250;


#-------!!!!!!!!-------!!!!!!!!-------!!!!!!!!-------!!!!!!!!-------!!!!!!!!-------#;
#####;
##This does stratified random sampling using field, product and year;
set.seed(1); #this is done for stratified() function;
#load(paste("farmDataSample", ".RData", sep=""));

### 13 August 2019;
### Do this possibly for a second paper, using repeat fertilizer observations for a given lonlat point over the four years;
#farmDataSample <- stratified(farmData2015to2018c, c("Field"), 1000); #122,031 observations;
###;

farmDataSample <- stratified(farmData2015to2018c, c("Field", "Product", "ProcYear"), 1000); #122,031 observations;
save(farmDataSample, file=paste("farmDataSample.RData", sep=""), compress="xz");
#-------!!!!!!!!-------!!!!!!!!-------!!!!!!!!-------!!!!!!!!-------!!!!!!!!-------#;



#head(farmDataSample);
table(farmDataSample$Field, farmDataSample$Product, farmDataSample$ProcYear);
#####;

summary(farmDataSample$lonlat);
length(unique(farmDataSample$lonlat)); #88,234 unique lonlat points for sample, #217,898 unique lonlat points for full data;
str(farmDataSample);

#---------------#;
priorFarm1 <- list(
  R=list(V=1, nu = 0),
  G=list(G1=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),
         G2=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),   
         G3=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),     
         G4=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2)));
                 
priorFarm2 <- list(
  R=list(V=1, nu = 0),
  G=list(G1=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),
         G2=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(10, 1),
                 alpha.V  = diag(1)*(5)^2),   
         G3=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(35, 1),
                 alpha.V  = diag(1)*(5)^2),     
         G4=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2)));
                 
priorFarm3 <- list(
  R=list(V=1, nu = 0),
  G=list(G1=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),
         G2=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(25, 1),
                 alpha.V  = diag(1)*(4)^2),     
         G3=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2)));

priorFarm4 <- list(
  R=list(V=1, nu = 0),
  G=list(G1=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),
         G2=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),  
         G3=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),
         G4=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),
         G5=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(5, 1),
                 alpha.V  = diag(1)*(10)^2),
         G6=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),
         G7=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2)));


## look at McDonald 2016 manuscript for setting priors;                 
priorFarm5 <- list(
  R=list(V=1, nu = 0),
  G=list(G1=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),
         G2=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(10)^2),  
         G3=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(10)^2),
         G4=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),
         G5=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(10)^2),
         G6=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),
         G7=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2)));


prior.ex <- list(G = list(G1 = list(V = 1, nu = 2, alpha.mu = 0, alpha.V = 10000),
                      G2 = list(V = 1, nu = 2, alpha.mu = 0, alpha.V = 10000),
                      G3 = list(V = 1, nu = 2, alpha.mu = 0, alpha.V = 10000),
                      G4 = list(V = 1, nu = 2, alpha.mu = 0, alpha.V = 10000),
                      G5 = list(V = 1, nu = 2, alpha.mu = 0, alpha.V = 10000),
                      G6 = list(V = 1, nu = 2, alpha.mu = 0, alpha.V = 10000),
                      G7 = list(V = 1, nu = 2, alpha.mu = 0, alpha.V = 10000)),
                      R = list(V=1, fix=1));

prior.m5 <- list(
  R=list(V=1, n=1, fix=1),
  G=list(G1=list(V        = diag(2),
                 n        = 2,
                 alpha.mu = rep(0, 2),
                 alpha.V  = diag(2)*25^2),
         G2=list(V        = diag(2),
                 n        = 2,
                 alpha.mu = rep(0, 2),
                 alpha.V  = diag(2)*25^2),
         G3=list(V        = diag(2),
                 n        = 2,
                 alpha.mu = rep(0, 2),
                 alpha.V  = diag(2)*25^2),
         G4=list(V        = diag(2),
                 n        = 2,
                 alpha.mu = rep(0, 2),
                 alpha.V  = diag(2)*25^2),
         G5=list(V        = diag(2),
                 n        = 2,
                 alpha.mu = rep(0, 2),
                 alpha.V  = diag(2)*25^2),  
         G6=list(V        = diag(2),
                 n        = 2,
                 alpha.mu = rep(0, 2),
                 alpha.V  = diag(2)*25^2),      
         G7=list(V        = diag(2),
                 n        = 2,
                 alpha.mu = rep(0, 2),
                 alpha.V  = diag(1)*25^2)));


priorFarmSel <- list(
  R=list(V=1, nu = 0),
  G=list(G1=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),
         G2=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2)));

priorFarmSel1 <- list(
  R=list(V=1, nu = 0),
  G=list(G1=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),
         G4=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2)));
                 
priorFarmSel2 <- list(
  R=list(V=1, nu = 0),
  G=list(G1=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),
         G2=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),
         G3=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),
         G4=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),
         G5=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2)));

priorsProper1 <- list(
  B=list(B1=list(mu=0, V=10^10),
  		 B2=list(mu=0, V=10^10),
  		 B3=list(mu=0, V=10^10),
		 B4=list(mu=0, V=10^10)),
  R=list(V=1, nu = 0.001),
  G=list(G1=list(V        = diag(1),
                 nu        = 0.001,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(10)^2),
         G2=list(V        = diag(1),
                 nu        = 0.001,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(10)^2),
         G3=list(V        = diag(1),
                 nu        = 0.001,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(10)^2),
         G4=list(V        = diag(1),
                 nu        = 0.001,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(10)^2),
         G5=list(V        = diag(1),
                 nu        = 0.001,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(10)^2),
         G6=list(V        = diag(1),
                 nu        = 0.001,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(10)^2),
         G7=list(V        = diag(1),
                 nu        = 0.001,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(10)^2)));

priorsProper1 <- list(
  R=list(V=1, nu = 0.001),
  G=list(G1=list(V        = diag(1),
                 nu        = 0.001,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(10)^2),
         G2=list(V        = diag(1),
                 nu        = 0.001,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(10)^2),
         G3=list(V        = diag(1),
                 nu        = 0.001,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(10)^2),
         G4=list(V        = diag(1),
                 nu        = 0.001,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(10)^2),
         G5=list(V        = diag(1),
                 nu        = 0.001,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(10)^2),
         G6=list(V        = diag(1),
                 nu        = 0.001,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(10)^2),
         G7=list(V        = diag(1),
                 nu        = 0.001,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(10)^2)));
                 
priorsProper <- list(
  R=list(V=1, nu = 0.001),
  G=list(G1=list(V        = diag(1),
                 nu        = 0.001,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),
         G2=list(V        = diag(1),
                 nu        = 0.001,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),
         G3=list(V        = diag(1),
                 nu        = 0.001,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),
         G4=list(V        = diag(1),
                 nu        = 0.001,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),
         G5=list(V        = diag(1),
                 nu        = 0.001,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),
         G6=list(V        = diag(1),
                 nu        = 0.001,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),
         G7=list(V        = diag(1),
                 nu        = 0.001,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2)));

priorFarmSel3 <- list(
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
         G5=list(V        = diag(1),
                 nu       = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),
         G6=list(V        = diag(1),
                 nu       = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),
         G7=list(V        = diag(1),
                 nu       = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2)));

##how to specify fixed effect as random slope;            
##https://stat.ethz.ch/pipermail/r-sig-mixed-models/2013q1/019896.html
###https://stat.ethz.ch/pipermail/r-sig-mixed-models/2011q2/016366.html


####-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----####;
####-----@-----@-----@-----@--- Does for all crops-@-----@-----@-----@-----@-----####;

#model1 <- MCMCglmm(DryYield ~ Seed_Rate + Ino_Rate + N_Rate +  P_Rate + K_Rate + S_Rate + Product + ProcYear + Machine, random=~Field + lonlat, family="gaussian", data=farmDataSample, saveX=TRUE, verbose=TRUE, burnin=55556, nitt=555600, thin=250, pr=TRUE, prior=priorFarmSel1);
##model1 <- MCMCglmm(DryYield ~ Seed_Rate + Ino_Rate + N_Rate +  P_Rate + K_Rate + S_Rate + Product + ProcYear + Machine, random=~Field + lonlat, family="gaussian", data=farmDataSample, saveX=TRUE, verbose=TRUE, burnin=500, nitt=2000, thin=25, pr=TRUE, prior=priorFarmSel1);

##model1 <- MCMCglmm(DryYield ~ Seed_Rate + Ino_Rate + N_Rate +  P_Rate + K_Rate + S_Rate + Product + Mat_Days + Machine + Elevation, random=~Field + ProcYear + lonlat, family="gaussian", data=farmDataSample, saveX=TRUE, verbose=TRUE, burnin=55556, nitt=555600, thin=250, pr=TRUE, prior=priorFarmSel1);
##model1 <- MCMCglmm(DryYield ~ Seed_Rate + Ino_Rate + N_Rate +  P_Rate + K_Rate + S_Rate + Product + Mat_Days + Machine + Elevation, random=~Field + ProcYear + lonlat, family="gaussian", data=farmDataSample, saveX=TRUE, verbose=TRUE, burnin=277778, nitt=2777800, thin=1250, pr=TRUE, prior=priorFarmSel1);
#3model1 <- MCMCglmm(DryYield ~ Seed_Rate + Ino_Rate + N_Rate +  P_Rate + K_Rate + S_Rate + Product + Mat_Days + Machine + Elevation, random=~Field + ProcYear + lonlat, family="gaussian", data=farmDataSample, saveX=TRUE, verbose=TRUE, burnin=500, nitt=2000, thin=25, pr=TRUE, prior=priorFarmSel1);

##model1 <- MCMCglmm(DryYield ~ Seed_Rate + Ino_Rate + N_Rate +  P_Rate + K_Rate + S_Rate + ProcYear + Product + Mat_Days + Machine + Elevation, random=~Field + lonlat, family="gaussian", data=farmDataSample, saveX=TRUE, verbose=TRUE, burnin=13890, nitt=138900, thin=25, pr=TRUE, prior=priorFarmSel);
##model1 <- MCMCglmm(DryYield ~ Seed_Rate + Ino_Rate + N_Rate +  P_Rate + K_Rate + S_Rate + ProcYear + Product + Mat_Days + Machine + Elevation, random=~Field + lonlat, family="gaussian", data=farmDataSample, saveX=TRUE, verbose=TRUE, burnin=500, nitt=2000, thin=25, pr=TRUE, prior=priorFarmSel);

##model1 <- MCMCglmm(DryYield ~ Seed_Rate + Ino_Rate + N_Rate +  P_Rate + K_Rate + S_Rate, random=~Field + ProcYear + Product + Mat_Days + Machine + Elevation + lonlat, family="gaussian", data=farmDataSample, saveX=TRUE, verbose=TRUE, burnin=55556, nitt=555600, thin=250, pr=TRUE, prior=priorFarm4);
##model1 <- MCMCglmm(DryYield ~ Seed_Rate + Ino_Rate + N_Rate +  P_Rate + K_Rate + S_Rate, random=~Field + ProcYear + Product + Mat_Days + Machine + Elevation + lonlat, family="gaussian", data=farmDataSample, saveX=TRUE, verbose=TRUE, burnin=500, nitt=2000, thin=25, pr=TRUE, prior=prior.m5);

##model1 <- MCMCglmm(DryYield ~ Product + Machine + Cultivar + Mat_Days + Seed_Rate + Ino_Rate + N_Rate +  P_Rate + K_Rate + S_Rate + Elevation, random=~Field + ProcYear + lonlat, family="gaussian", data=farmDataSample, saveX=TRUE, verbose=TRUE, burnin=30000, nitt=280000, thin=50, pr=TRUE, prior=priorFarm2);
##model1 <- MCMCglmm(DryYield ~ Product + Machine + Cultivar + Mat_Days + Seed_Rate + Ino_Rate + N_Rate +  P_Rate + K_Rate + S_Rate + Elevation, random=~Field + ProcYear + lonlat, family="gaussian", data=farmDataSample, saveX=TRUE, verbose=TRUE, burnin=30000, nitt=280000, thin=50, pr=TRUE, prior=priorFarm1);
##model1 <- MCMCglmm(DryYield ~ Product + Machine + Cultivar + Mat_Days + Seed_Rate + Ino_Rate + N_Rate +  P_Rate + K_Rate + S_Rate + Elevation, random=~Field + ProcYear + lonlat, family="gaussian", data=farmDataSample, saveX=TRUE, verbose=TRUE, burnin=500, nitt=2000, thin=5, pr=TRUE, prior=priorFarm1);



#### https://github.com/tmalsburg/MCMCglmm-intro
#### Try this mc.cores option for Gelman-Rubin criterion to check convergence;
#set.seed(1)
#m6 <- mclapply(1:4, function(i) {
#	model1 <- MCMCglmm(DryYield ~ Seed_Rate + Ino_Rate + N_Rate +  P_Rate + K_Rate + S_Rate + ProcYear + Product + Mat_Days + Machine + Elevation, random=~Field + lonlat, family="gaussian", data=farmDataSample, saveX=TRUE, verbose=TRUE, burnin=500, nitt=2000, thin=25, pr=TRUE, prior=priorFarmSel);
#}, mc.cores=4)
#
#m6 <- lapply(m6, function(m) m$Sol)
#m6 <- do.call(mcmc.list, m6)
#
#save(m6, file=paste("farm_summary_Model_1.RData", sep=""), compress="xz");


#setwd("/Users/tylerpittman/Farm/agYieldProject/data/STR_model_published");
#save(model1, file=paste("farm_summary_Model_1.RData", sep=""), compress="xz");
#setwd("/Users/tylerpittman/Farm/agYieldProject/data");
load(paste("STR_model_published/farm_summary_Model_1.RData", sep=""));

#https://stats.stackexchange.com/questions/307021/why-does-mcmcglmm-result-in-small-effective-sample-sizes-for-logistic-regression/310291
#autocorr.diag(model1$VCV); #use this to choose thinning interval, autocorrelation less than 0.01 ideally;

fe1 <- summary(model1);
#str(fe1$solutions);
fe1$solutions[,1] <- fe1$solutions[,1];
fe1$solutions[,2] <- fe1$solutions[,2];
fe1$solutions[,3] <- fe1$solutions[,3];

writeLines(capture.output(fe1), paste("farm_Model_1.txt", sep=""));
writeLines(capture.output(summary(model1)), paste("farm_summary_Model_1.txt", sep=""));

### DO THIS FOR DIAGNOSTICS OF HOW CHAIN FIT with mcmcglmm;
#http://www.maths.bath.ac.uk/~jjf23/mixchange/split.html

#colnames(model1$Sol); #1:15 are explanatory variables;
EVs1 <- model1$Sol[, 1:15];
intFields1 <- model1$Sol[, 16:46];
#intProcYears1 <- model1$Sol[, 49:52];

length(posterior.mode(intFields1)); #31;  #mean intercept;

### https://github.com/tmalsburg/MCMCglmm-intro
###

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("model1_VCV_farm.pdf", sep=''), bg="transparent", width=10, height=8, pointsize=12, family="FreeSans");
par(mfrow=c(length(colnames(model1$VCV)),2), mar=c(2,2,1,0));
plot(model1$VCV, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("model1_sol_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(EVs1)),2), mar=c(2,2,1,0));
plot(EVs1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("model1_Fields_farm.pdf", sep=''), bg="transparent", width=11, height=30, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intFields1)),2), mar=c(2,2,1,0));
plot(intFields1, auto.layout=F);
dev.off()

#pdf.options(encoding='CP1250');
#cairo_pdf(file=paste("model1_ProcYear_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
#par(mfrow=c(length(colnames(intProcYears1)),2), mar=c(2,2,1,0));
#plot(intProcYears1, auto.layout=F);
#dev.off()

#---;

####-----^*^-----^*^-----^*^-----^*^-----^*^-----^*^-----^*^-----^*^-----####;
####-----^*^-----^*^----- summary statistics-----^*^-----^*^-----^*^-----####;
farmDataSample1 <- farmDataSample[which(farmDataSample$Product == "Barley"), ];
length(farmDataSample1$Seed_Rate);
summary(farmDataSample1$DryYield);
summary(farmDataSample1$Seed_Rate);
summary(farmDataSample1$N_Rate);
summary(farmDataSample1$P_Rate);
summary(farmDataSample1$K_Rate);
summary(farmDataSample1$S_Rate);
summary(farmDataSample1$Ino_Rate);
sd(farmDataSample1$DryYield);
sd(farmDataSample1$Seed_Rate);
sd(farmDataSample1$N_Rate);
sd(farmDataSample1$P_Rate);
sd(farmDataSample1$K_Rate);
sd(farmDataSample1$S_Rate);
sd(farmDataSample1$Ino_Rate);

farmDataSample1 <- farmDataSample[which(farmDataSample$Product == "Canary"), ];
length(farmDataSample1$Seed_Rate);
summary(farmDataSample1$DryYield);
summary(farmDataSample1$Seed_Rate);
summary(farmDataSample1$N_Rate);
summary(farmDataSample1$P_Rate);
summary(farmDataSample1$K_Rate);
summary(farmDataSample1$S_Rate);
summary(farmDataSample1$Ino_Rate);
sd(farmDataSample1$DryYield);
sd(farmDataSample1$Seed_Rate);
sd(farmDataSample1$N_Rate);
sd(farmDataSample1$P_Rate);
sd(farmDataSample1$K_Rate);
sd(farmDataSample1$S_Rate);
sd(farmDataSample1$Ino_Rate);

farmDataSample1 <- farmDataSample[which(farmDataSample$Product == "Canola"), ];
length(farmDataSample1$Seed_Rate);
summary(farmDataSample1$DryYield);
summary(farmDataSample1$Seed_Rate);
summary(farmDataSample1$N_Rate);
summary(farmDataSample1$P_Rate);
summary(farmDataSample1$K_Rate);
summary(farmDataSample1$S_Rate);
summary(farmDataSample1$Ino_Rate);
sd(farmDataSample1$DryYield);
sd(farmDataSample1$Seed_Rate);
sd(farmDataSample1$N_Rate);
sd(farmDataSample1$P_Rate);
sd(farmDataSample1$K_Rate);
sd(farmDataSample1$S_Rate);
sd(farmDataSample1$Ino_Rate);

farmDataSample1 <- farmDataSample[which(farmDataSample$Product == "Lentils"), ];
length(farmDataSample1$Seed_Rate);
summary(farmDataSample1$DryYield);
summary(farmDataSample1$Seed_Rate);
summary(farmDataSample1$N_Rate);
summary(farmDataSample1$P_Rate);
summary(farmDataSample1$K_Rate);
summary(farmDataSample1$S_Rate);
summary(farmDataSample1$Ino_Rate);
sd(farmDataSample1$DryYield);
sd(farmDataSample1$Seed_Rate);
sd(farmDataSample1$N_Rate);
sd(farmDataSample1$P_Rate);
sd(farmDataSample1$K_Rate);
sd(farmDataSample1$S_Rate);
sd(farmDataSample1$Ino_Rate);

farmDataSample1 <- farmDataSample[which(farmDataSample$Product == "Wheat"), ];
length(farmDataSample1$Seed_Rate);
summary(farmDataSample1$DryYield);
summary(farmDataSample1$Seed_Rate);
summary(farmDataSample1$N_Rate);
summary(farmDataSample1$P_Rate);
summary(farmDataSample1$K_Rate);
summary(farmDataSample1$S_Rate);
summary(farmDataSample1$Ino_Rate);
sd(farmDataSample1$DryYield);
sd(farmDataSample1$Seed_Rate);
sd(farmDataSample1$N_Rate);
sd(farmDataSample1$P_Rate);
sd(farmDataSample1$K_Rate);
sd(farmDataSample1$S_Rate);
sd(farmDataSample1$Ino_Rate);


####-----^*^-----^*^-----^*^-----^*^-----^*^-----^*^-----^*^-----^*^-----####;
####-----^*^-----^*^-----^*^-----^*^-----^*^-----^*^-----^*^-----^*^-----####;




####-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----####;
####-----@-----@-----@-----@--- Does for barley ---@-----@-----@-----@-----@-----####;

farmDataSample1 <- farmDataSample[which(farmDataSample$Product == "Barley"), ];
str(farmDataSample1);

#model1 <- MCMCglmm(DryYield ~ Seed_Rate + N_Rate + P_Rate + K_Rate, random=~Field + ProcYear + Machine + Cultivar + PrevCrop + Field:lonlat + Field:Elevation, family="gaussian", data=farmDataSample1, saveX=TRUE, verbose=TRUE, burnin=27778, nitt=277800, thin=125, pr=TRUE, prior=priorFarmSel3);
#model1 <- MCMCglmm(DryYield ~ Seed_Rate + N_Rate + P_Rate + K_Rate, random=~Field + ProcYear + Machine + Cultivar + PrevCrop + Field:lonlat + Field:Elevation, family="gaussian", data=farmDataSample1, saveX=TRUE, verbose=TRUE, burnin=500, nitt=2000, thin=25, pr=TRUE, prior=priorFarmSel3);

#save(model1, file=paste("farm_summary_barley_Model_1.RData", sep=""), compress="xz");
#load(paste("farm_summary_barley_Model_1.RData", sep=""));
load(paste("STR_model_published/farm_summary_barley_Model_1.RData", sep=""));

#https://stats.stackexchange.com/questions/307021/why-does-mcmcglmm-result-in-small-effective-sample-sizes-for-logistic-regression/310291
#autocorr.diag(model1$VCV); #use this to choose thinning interval, autocorrelation less than 0.01 ideally;

fe1 <- summary(model1);
#str(fe1$solutions);
fe1$solutions[,1] <- fe1$solutions[,1];
fe1$solutions[,2] <- fe1$solutions[,2];
fe1$solutions[,3] <- fe1$solutions[,3];

writeLines(capture.output(fe1), paste("farm_barley_Model_1.txt", sep=""));
writeLines(capture.output(summary(model1)), paste("farm_summary_barley_Model_1.txt", sep=""));

### DO THIS FOR DIAGNOSTICS OF HOW CHAIN FIT with mcmcglmm;
#http://www.maths.bath.ac.uk/~jjf23/mixchange/split.html

#colnames(model1$Sol); #1:5 are explanatory variables;
EVs1 <- model1$Sol[, 1:5];
intFields1 <- model1$Sol[, 6:8];
intProcYears1 <- model1$Sol[, 8:9];
intMachine1 <- model1$Sol[, 10:11];
intCultivar1 <- model1$Sol[, 12:13];
intPrevCrop1 <- model1$Sol[, 13:14];

length(posterior.mode(intFields1)); #3;  #mean intercept;

### https://github.com/tmalsburg/MCMCglmm-intro
###

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("barley_model1_VCV_farm.pdf", sep=''), bg="transparent", width=10, height=8, pointsize=12, family="FreeSans");
par(mfrow=c(length(colnames(model1$VCV)),2), mar=c(2,2,1,0));
plot(model1$VCV, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("barley_model1_sol_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(EVs1)),2), mar=c(2,2,1,0));
plot(EVs1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("barley_model1_Fields_farm.pdf", sep=''), bg="transparent", width=11, height=30, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intFields1)),2), mar=c(2,2,1,0));
plot(intFields1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("barley_model1_ProcYear_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intProcYears1)),2), mar=c(2,2,1,0));
plot(intProcYears1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("barley_model1_Machine_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intMachine1)),2), mar=c(2,2,1,0));
plot(intMachine1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("barley_model1_Cultivar_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intCultivar1)),2), mar=c(2,2,1,0));
plot(intCultivar1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("barley_model1_PrevCrop_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intPrevCrop1)),2), mar=c(2,2,1,0));
plot(intPrevCrop1, auto.layout=F);
dev.off()

colnames(model1$VCV);
ICC_Field <- model1$VCV[, 1]/(rowSums(model1$VCV)); 
ICC_ProcYear <- model1$VCV[, 2]/(rowSums(model1$VCV)); 
ICC_Machine <- model1$VCV[, 3]/(rowSums(model1$VCV)); 
ICC_Cultivar <- model1$VCV[, 4]/(rowSums(model1$VCV)); 
ICC_PrevCrop <- model1$VCV[, 5]/(rowSums(model1$VCV)); 
ICC_lonlat <- model1$VCV[, 6]/(rowSums(model1$VCV)); 
ICC_elevation <- model1$VCV[, 7]/(rowSums(model1$VCV)); 
dft1.1 <- cbind(ICC = posterior.mode(ICC_Field), CI = HPDinterval(ICC_Field));
dft2.1 <- cbind(ICC = posterior.mode(ICC_ProcYear), CI = HPDinterval(ICC_ProcYear));
dft3.1 <- cbind(ICC = posterior.mode(ICC_Machine), CI = HPDinterval(ICC_Machine));
dft4.1 <- cbind(ICC = posterior.mode(ICC_Cultivar), CI = HPDinterval(ICC_Cultivar));
dft5.1 <- cbind(ICC = posterior.mode(ICC_PrevCrop), CI = HPDinterval(ICC_PrevCrop));
dft6.1 <- cbind(ICC = posterior.mode(ICC_lonlat), CI = HPDinterval(ICC_lonlat));
dft7.1 <- cbind(ICC = posterior.mode(ICC_elevation), CI = HPDinterval(ICC_elevation));
dft1 <- rbind(dft1.1, dft2.1, dft3.1, dft4.1, dft5.1, dft6.1, dft7.1);
row.names(dft1) <- c("ICC_Field", "ICC_ProcYear", "ICC_Machine", "ICC_Cultivar", "ICC_PrevCrop", "ICC_lonlat", "ICC_elevation");
writeLines(capture.output(dft1), paste("barley_ICC_Model_1.txt", sep=""));

#---;

####-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----####;
####-----@-----@-----@-----@--- Does for canary ---@-----@-----@-----@-----@-----####;

farmDataSample1 <- farmDataSample[which(farmDataSample$Product == "Canary"), ];
str(farmDataSample1);

#model1 <- MCMCglmm(DryYield ~ Seed_Rate + N_Rate + P_Rate + K_Rate, random=~Field + ProcYear + Machine + Cultivar + PrevCrop + Field:lonlat + Field:Elevation, family="gaussian", data=farmDataSample1, saveX=TRUE, verbose=TRUE, burnin=27778, nitt=277800, thin=125, pr=TRUE, prior=priorFarmSel3);
#model1 <- MCMCglmm(DryYield ~ Seed_Rate + N_Rate + P_Rate + K_Rate, random=~Field + ProcYear + Machine + Cultivar + PrevCrop + Field:lonlat + Field:Elevation, family="gaussian", data=farmDataSample1, saveX=TRUE, verbose=TRUE, burnin=500, nitt=2000, thin=25, pr=TRUE, prior=priorFarmSel3);

#save(model1, file=paste("farm_summary_canary_Model_1.RData", sep=""), compress="xz");
#load(paste("farm_summary_canary_Model_1.RData", sep=""));
load(paste("STR_model_published/farm_summary_canary_Model_1.RData", sep=""));

#https://stats.stackexchange.com/questions/307021/why-does-mcmcglmm-result-in-small-effective-sample-sizes-for-logistic-regression/310291
#autocorr.diag(model1$VCV); #use this to choose thinning interval, autocorrelation less than 0.01 ideally;

fe1 <- summary(model1);
#str(fe1$solutions);
fe1$solutions[,1] <- fe1$solutions[,1];
fe1$solutions[,2] <- fe1$solutions[,2];
fe1$solutions[,3] <- fe1$solutions[,3];

writeLines(capture.output(fe1), paste("farm_canary_Model_1.txt", sep=""));
writeLines(capture.output(summary(model1)), paste("farm_summary_canary_Model_1.txt", sep=""));

### DO THIS FOR DIAGNOSTICS OF HOW CHAIN FIT with mcmcglmm;
#http://www.maths.bath.ac.uk/~jjf23/mixchange/split.html

#colnames(model1$Sol); #1:5 are explanatory variables;
EVs1 <- model1$Sol[, 1:5];
intFields1 <- model1$Sol[, 6:11];
intProcYears1 <- model1$Sol[, 12:13];
intMachine1 <- model1$Sol[, 14:15];
intCultivar1 <- model1$Sol[, 16:17];
intPrevCrop1 <- model1$Sol[, 17:18];

length(posterior.mode(intFields1)); #6;  #mean intercept;

### https://github.com/tmalsburg/MCMCglmm-intro
###

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("canary_model1_VCV_farm.pdf", sep=''), bg="transparent", width=10, height=8, pointsize=12, family="FreeSans");
par(mfrow=c(length(colnames(model1$VCV)),2), mar=c(2,2,1,0));
plot(model1$VCV, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("canary_model1_sol_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(EVs1)),2), mar=c(2,2,1,0));
plot(EVs1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("canary_model1_Fields_farm.pdf", sep=''), bg="transparent", width=11, height=30, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intFields1)),2), mar=c(2,2,1,0));
plot(intFields1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("canary_model1_ProcYear_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intProcYears1)),2), mar=c(2,2,1,0));
plot(intProcYears1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("canary_model1_Machine_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intMachine1)),2), mar=c(2,2,1,0));
plot(intMachine1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("canary_model1_Cultivar_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intCultivar1)),2), mar=c(2,2,1,0));
plot(intCultivar1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("canary_model1_PrevCrop_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intPrevCrop1)),2), mar=c(2,2,1,0));
plot(intPrevCrop1, auto.layout=F);
dev.off()

colnames(model1$VCV);
ICC_Field <- model1$VCV[, 1]/(rowSums(model1$VCV)); 
ICC_ProcYear <- model1$VCV[, 2]/(rowSums(model1$VCV)); 
ICC_Machine <- model1$VCV[, 3]/(rowSums(model1$VCV)); 
ICC_Cultivar <- model1$VCV[, 4]/(rowSums(model1$VCV)); 
ICC_PrevCrop <- model1$VCV[, 5]/(rowSums(model1$VCV)); 
ICC_lonlat <- model1$VCV[, 6]/(rowSums(model1$VCV)); 
ICC_elevation <- model1$VCV[, 7]/(rowSums(model1$VCV)); 
dft1.1 <- cbind(ICC = posterior.mode(ICC_Field), CI = HPDinterval(ICC_Field));
dft2.1 <- cbind(ICC = posterior.mode(ICC_ProcYear), CI = HPDinterval(ICC_ProcYear));
dft3.1 <- cbind(ICC = posterior.mode(ICC_Machine), CI = HPDinterval(ICC_Machine));
dft4.1 <- cbind(ICC = posterior.mode(ICC_Cultivar), CI = HPDinterval(ICC_Cultivar));
dft5.1 <- cbind(ICC = posterior.mode(ICC_PrevCrop), CI = HPDinterval(ICC_PrevCrop));
dft6.1 <- cbind(ICC = posterior.mode(ICC_lonlat), CI = HPDinterval(ICC_lonlat));
dft7.1 <- cbind(ICC = posterior.mode(ICC_elevation), CI = HPDinterval(ICC_elevation));
dft1 <- rbind(dft1.1, dft2.1, dft3.1, dft4.1, dft5.1, dft6.1, dft7.1);
row.names(dft1) <- c("ICC_Field", "ICC_ProcYear", "ICC_Machine", "ICC_Cultivar", "ICC_PrevCrop", "ICC_lonlat", "ICC_elevation");
writeLines(capture.output(dft1), paste("canary_ICC_Model_1.txt", sep=""));

#---;

####-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----####;
####-----@-----@-----@-----@--- Does for canola ---@-----@-----@-----@-----@-----####;

farmDataSample1 <- farmDataSample[which(farmDataSample$Product == "Canola"), ];
str(farmDataSample1);

#model1 <- MCMCglmm(DryYield ~ Seed_Rate + N_Rate +  P_Rate + S_Rate, random=~Field + ProcYear + Machine + Cultivar + PrevCrop + Field:lonlat + Field:Elevation, family="gaussian", data=farmDataSample1, saveX=TRUE, verbose=TRUE, burnin=27778, nitt=277800, thin=125, pr=TRUE, prior=priorFarmSel3);
#model1 <- MCMCglmm(DryYield ~ Seed_Rate + N_Rate +  P_Rate + S_Rate, random=~Field + ProcYear + Machine + Cultivar + PrevCrop + Field:lonlat + Field:Elevation, family="gaussian", data=farmDataSample1, saveX=TRUE, verbose=TRUE, burnin=500, nitt=2000, thin=25, pr=TRUE, prior=priorFarmSel3);

#save(model1, file=paste("farm_summary_canola_Model_1.RData", sep=""), compress="xz");
#load(paste("farm_summary_canola_Model_1.RData", sep=""));
load(paste("STR_model_published/farm_summary_canola_Model_1.RData", sep=""));

#https://stats.stackexchange.com/questions/307021/why-does-mcmcglmm-result-in-small-effective-sample-sizes-for-logistic-regression/310291
#autocorr.diag(model1$VCV); #use this to choose thinning interval, autocorrelation less than 0.01 ideally;

fe1 <- summary(model1);
#str(fe1$solutions);
fe1$solutions[,1] <- fe1$solutions[,1];
fe1$solutions[,2] <- fe1$solutions[,2];
fe1$solutions[,3] <- fe1$solutions[,3];

writeLines(capture.output(fe1), paste("farm_canola_Model_1.txt", sep=""));
writeLines(capture.output(summary(model1)), paste("farm_summary_canola_Model_1.txt", sep=""));

### DO THIS FOR DIAGNOSTICS OF HOW CHAIN FIT with mcmcglmm;
#http://www.maths.bath.ac.uk/~jjf23/mixchange/split.html

#colnames(model1$Sol); #1:5 are explanatory variables;
EVs1 <- model1$Sol[, 1:5];
intFields1 <- model1$Sol[, 6:21];
intProcYears1 <- model1$Sol[, 22:24];
intMachine1 <- model1$Sol[, 25:26];
intCultivar1 <- model1$Sol[, 27:29];
intPrevCrop1 <- model1$Sol[, 30:33];

length(posterior.mode(intFields1)); #16;  #mean intercept;

### https://github.com/tmalsburg/MCMCglmm-intro
###

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("canola_model1_VCV_farm.pdf", sep=''), bg="transparent", width=10, height=8, pointsize=12, family="FreeSans");
par(mfrow=c(length(colnames(model1$VCV)),2), mar=c(2,2,1,0));
plot(model1$VCV, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("canola_model1_sol_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(EVs1)),2), mar=c(2,2,1,0));
plot(EVs1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("canola_model1_Fields_farm.pdf", sep=''), bg="transparent", width=11, height=30, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intFields1)),2), mar=c(2,2,1,0));
plot(intFields1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("canola_model1_ProcYear_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intProcYears1)),2), mar=c(2,2,1,0));
plot(intProcYears1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("canola_model1_Machine_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intMachine1)),2), mar=c(2,2,1,0));
plot(intMachine1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("canola_model1_Cultivar_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intCultivar1)),2), mar=c(2,2,1,0));
plot(intCultivar1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("canola_model1_PrevCrop_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intPrevCrop1)),2), mar=c(2,2,1,0));
plot(intPrevCrop1, auto.layout=F);
dev.off()

colnames(model1$VCV);
ICC_Field <- model1$VCV[, 1]/(rowSums(model1$VCV)); 
ICC_ProcYear <- model1$VCV[, 2]/(rowSums(model1$VCV)); 
ICC_Machine <- model1$VCV[, 3]/(rowSums(model1$VCV)); 
ICC_Cultivar <- model1$VCV[, 4]/(rowSums(model1$VCV)); 
ICC_PrevCrop <- model1$VCV[, 5]/(rowSums(model1$VCV)); 
ICC_lonlat <- model1$VCV[, 6]/(rowSums(model1$VCV)); 
ICC_elevation <- model1$VCV[, 7]/(rowSums(model1$VCV)); 
dft1.1 <- cbind(ICC = posterior.mode(ICC_Field), CI = HPDinterval(ICC_Field));
dft2.1 <- cbind(ICC = posterior.mode(ICC_ProcYear), CI = HPDinterval(ICC_ProcYear));
dft3.1 <- cbind(ICC = posterior.mode(ICC_Machine), CI = HPDinterval(ICC_Machine));
dft4.1 <- cbind(ICC = posterior.mode(ICC_Cultivar), CI = HPDinterval(ICC_Cultivar));
dft5.1 <- cbind(ICC = posterior.mode(ICC_PrevCrop), CI = HPDinterval(ICC_PrevCrop));
dft6.1 <- cbind(ICC = posterior.mode(ICC_lonlat), CI = HPDinterval(ICC_lonlat));
dft7.1 <- cbind(ICC = posterior.mode(ICC_elevation), CI = HPDinterval(ICC_elevation));
dft1 <- rbind(dft1.1, dft2.1, dft3.1, dft4.1, dft5.1, dft6.1, dft7.1);
row.names(dft1) <- c("ICC_Field", "ICC_ProcYear", "ICC_Machine", "ICC_Cultivar", "ICC_PrevCrop", "ICC_lonlat", "ICC_elevation");
writeLines(capture.output(dft1), paste("canola_ICC_Model_1.txt", sep=""));

#---;

####-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----####;
####-----@-----@-----@-----@--- Does for lentils ---@-----@-----@-----@-----@-----####;

farmDataSample1 <- farmDataSample[which(farmDataSample$Product == "Lentils"), ];
str(farmDataSample1);


#model1 <- MCMCglmm(DryYield ~ Seed_Rate + Ino_Rate + N_Rate + P_Rate + K_Rate, random=~Field + ProcYear + Machine + Cultivar + PrevCrop + Field:lonlat + Field:Elevation, family="gaussian", data=farmDataSample1, saveX=TRUE, verbose=TRUE, burnin=27778, nitt=277800, thin=125, pr=TRUE, prior=priorFarmSel3);
#model1 <- MCMCglmm(DryYield ~ Seed_Rate + Ino_Rate + N_Rate + P_Rate + K_Rate, random=~Field + ProcYear + Machine + Cultivar + PrevCrop + Field:lonlat + Field:Elevation, family="gaussian", data=farmDataSample1, saveX=TRUE, verbose=TRUE, burnin=500, nitt=2000, thin=25, pr=TRUE, prior=priorFarmSel3);

#save(model1, file=paste("farm_summary_lentils_Model_1.RData", sep=""), compress="xz");
#load(paste("farm_summary_lentils_Model_1.RData", sep=""));
load(paste("STR_model_published/farm_summary_lentils_Model_1.RData", sep=""));

#https://stats.stackexchange.com/questions/307021/why-does-mcmcglmm-result-in-small-effective-sample-sizes-for-logistic-regression/310291
#autocorr.diag(model1$VCV); #use this to choose thinning interval, autocorrelation less than 0.01 ideally;

fe1 <- summary(model1);
#str(fe1$solutions);
fe1$solutions[,1] <- fe1$solutions[,1];
fe1$solutions[,2] <- fe1$solutions[,2];
fe1$solutions[,3] <- fe1$solutions[,3];

writeLines(capture.output(fe1), paste("farm_lentils_Model_1.txt", sep=""));
writeLines(capture.output(summary(model1)), paste("farm_summary_lentils_Model_1.txt", sep=""));

### DO THIS FOR DIAGNOSTICS OF HOW CHAIN FIT with mcmcglmm;
#http://www.maths.bath.ac.uk/~jjf23/mixchange/split.html

#colnames(model1$Sol); #1:7 are explanatory variables;
EVs1 <- model1$Sol[, 1:6];
intFields1 <- model1$Sol[, 7:36];
intProcYears1 <- model1$Sol[, 37:40];
intMachine1 <- model1$Sol[, 41:42];
intCultivar1 <- model1$Sol[, 43:44];
intPrevCrop1 <- model1$Sol[, 45:47];

length(posterior.mode(intFields1)); #30;  #mean intercept;

### https://github.com/tmalsburg/MCMCglmm-intro
###

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("lentils_model1_VCV_farm.pdf", sep=''), bg="transparent", width=10, height=8, pointsize=12, family="FreeSans");
par(mfrow=c(length(colnames(model1$VCV)),2), mar=c(2,2,1,0));
plot(model1$VCV, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("lentils_model1_sol_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(EVs1)),2), mar=c(2,2,1,0));
plot(EVs1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("lentils_model1_Fields_farm.pdf", sep=''), bg="transparent", width=11, height=30, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intFields1)),2), mar=c(2,2,1,0));
plot(intFields1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("lentils_model1_ProcYear_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intProcYears1)),2), mar=c(2,2,1,0));
plot(intProcYears1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("lentils_model1_Machine_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intMachine1)),2), mar=c(2,2,1,0));
plot(intMachine1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("lentils_model1_Cultivar_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intCultivar1)),2), mar=c(2,2,1,0));
plot(intCultivar1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("lentils_model1_PrevCrop_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intPrevCrop1)),2), mar=c(2,2,1,0));
plot(intPrevCrop1, auto.layout=F);
dev.off()

colnames(model1$VCV);
ICC_Field <- model1$VCV[, 1]/(rowSums(model1$VCV)); 
ICC_ProcYear <- model1$VCV[, 2]/(rowSums(model1$VCV)); 
ICC_Machine <- model1$VCV[, 3]/(rowSums(model1$VCV)); 
ICC_Cultivar <- model1$VCV[, 4]/(rowSums(model1$VCV)); 
ICC_PrevCrop <- model1$VCV[, 5]/(rowSums(model1$VCV)); 
ICC_lonlat <- model1$VCV[, 6]/(rowSums(model1$VCV)); 
ICC_elevation <- model1$VCV[, 7]/(rowSums(model1$VCV)); 
dft1.1 <- cbind(ICC = posterior.mode(ICC_Field), CI = HPDinterval(ICC_Field));
dft2.1 <- cbind(ICC = posterior.mode(ICC_ProcYear), CI = HPDinterval(ICC_ProcYear));
dft3.1 <- cbind(ICC = posterior.mode(ICC_Machine), CI = HPDinterval(ICC_Machine));
dft4.1 <- cbind(ICC = posterior.mode(ICC_Cultivar), CI = HPDinterval(ICC_Cultivar));
dft5.1 <- cbind(ICC = posterior.mode(ICC_PrevCrop), CI = HPDinterval(ICC_PrevCrop));
dft6.1 <- cbind(ICC = posterior.mode(ICC_lonlat), CI = HPDinterval(ICC_lonlat));
dft7.1 <- cbind(ICC = posterior.mode(ICC_elevation), CI = HPDinterval(ICC_elevation));
dft1 <- rbind(dft1.1, dft2.1, dft3.1, dft4.1, dft5.1, dft6.1, dft7.1);
row.names(dft1) <- c("ICC_Field", "ICC_ProcYear", "ICC_Machine", "ICC_Cultivar", "ICC_PrevCrop", "ICC_lonlat", "ICC_elevation");
writeLines(capture.output(dft1), paste("lentils_ICC_Model_1.txt", sep=""));

#---;

####-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----@-----####;
####-----@-----@-----@-----@--- Does for wheat ---@-----@-----@-----@-----@-----####;

farmDataSample1 <- farmDataSample[which(farmDataSample$Product == "Wheat"), ];
str(farmDataSample1);

#model1 <- MCMCglmm(DryYield ~ Seed_Rate + N_Rate +  P_Rate + K_Rate, random=~Field + ProcYear + Machine + Cultivar + PrevCrop + Field:lonlat + Field:Elevation, family="gaussian", data=farmDataSample1, saveX=TRUE, verbose=TRUE, burnin=27778, nitt=277800, thin=125, pr=TRUE, prior=priorFarmSel3);
#model1 <- MCMCglmm(DryYield ~ Seed_Rate + N_Rate +  P_Rate + K_Rate, random=~Field + ProcYear + Machine + Cultivar + PrevCrop + Field:lonlat + Field:Elevation, family="gaussian", data=farmDataSample1, saveX=TRUE, verbose=TRUE, burnin=500, nitt=2000, thin=25, pr=TRUE, prior=priorFarmSel3);

#save(model1, file=paste("farm_summary_wheat_Model_1.RData", sep=""), compress="xz");
#load(paste("farm_summary_wheat_Model_1.RData", sep=""));
load(paste("STR_model_published/farm_summary_wheat_Model_1.RData", sep=""));


#https://stats.stackexchange.com/questions/307021/why-does-mcmcglmm-result-in-small-effective-sample-sizes-for-logistic-regression/310291
#autocorr.diag(model1$VCV); #use this to choose thinning interval, autocorrelation less than 0.01 ideally;

fe1 <- summary(model1);
#str(fe1$solutions);
fe1$solutions[,1] <- fe1$solutions[,1];
fe1$solutions[,2] <- fe1$solutions[,2];
fe1$solutions[,3] <- fe1$solutions[,3];

writeLines(capture.output(fe1), paste("farm_wheat_Model_1.txt", sep=""));
writeLines(capture.output(summary(model1)), paste("farm_summary_wheat_Model_1.txt", sep=""));

### DO THIS FOR DIAGNOSTICS OF HOW CHAIN FIT with mcmcglmm;
#http://www.maths.bath.ac.uk/~jjf23/mixchange/split.html

#colnames(model1$Sol); #1:5 are explanatory variables;
EVs1 <- model1$Sol[, 1:5];
intFields1 <- model1$Sol[, 6:36];
intProcYears1 <- model1$Sol[, 37:40];
intMachine1 <- model1$Sol[, 41:42];
intCultivar1 <- model1$Sol[, 43:44];
intPrevCrop1 <- model1$Sol[, 44:46];

length(posterior.mode(intFields1)); #31;  #mean intercept;

### https://github.com/tmalsburg/MCMCglmm-intro
###

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("wheat_model1_VCV_farm.pdf", sep=''), bg="transparent", width=10, height=8, pointsize=12, family="FreeSans");
par(mfrow=c(length(colnames(model1$VCV)),2), mar=c(2,2,1,0));
plot(model1$VCV, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("wheat_model1_sol_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(EVs1)),2), mar=c(2,2,1,0));
plot(EVs1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("wheat_model1_Fields_farm.pdf", sep=''), bg="transparent", width=11, height=30, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intFields1)),2), mar=c(2,2,1,0));
plot(intFields1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("wheat_model1_ProcYear_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intProcYears1)),2), mar=c(2,2,1,0));
plot(intProcYears1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("wheat_model1_Machine_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intMachine1)),2), mar=c(2,2,1,0));
plot(intMachine1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("wheat_model1_Cultivar_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intCultivar1)),2), mar=c(2,2,1,0));
plot(intCultivar1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("wheat_model1_PrevCrop_farm.pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intPrevCrop1)),2), mar=c(2,2,1,0));
plot(intPrevCrop1, auto.layout=F);
dev.off()

colnames(model1$VCV);
ICC_Field <- model1$VCV[, 1]/(rowSums(model1$VCV)); 
ICC_ProcYear <- model1$VCV[, 2]/(rowSums(model1$VCV)); 
ICC_Machine <- model1$VCV[, 3]/(rowSums(model1$VCV)); 
ICC_Cultivar <- model1$VCV[, 4]/(rowSums(model1$VCV)); 
ICC_PrevCrop <- model1$VCV[, 5]/(rowSums(model1$VCV)); 
ICC_lonlat <- model1$VCV[, 6]/(rowSums(model1$VCV)); 
ICC_elevation <- model1$VCV[, 7]/(rowSums(model1$VCV)); 
dft1.1 <- cbind(ICC = posterior.mode(ICC_Field), CI = HPDinterval(ICC_Field));
dft2.1 <- cbind(ICC = posterior.mode(ICC_ProcYear), CI = HPDinterval(ICC_ProcYear));
dft3.1 <- cbind(ICC = posterior.mode(ICC_Machine), CI = HPDinterval(ICC_Machine));
dft4.1 <- cbind(ICC = posterior.mode(ICC_Cultivar), CI = HPDinterval(ICC_Cultivar));
dft5.1 <- cbind(ICC = posterior.mode(ICC_PrevCrop), CI = HPDinterval(ICC_PrevCrop));
dft6.1 <- cbind(ICC = posterior.mode(ICC_lonlat), CI = HPDinterval(ICC_lonlat));
dft7.1 <- cbind(ICC = posterior.mode(ICC_elevation), CI = HPDinterval(ICC_elevation));
dft1 <- rbind(dft1.1, dft2.1, dft3.1, dft4.1, dft5.1, dft6.1, dft7.1);
row.names(dft1) <- c("ICC_Field", "ICC_ProcYear", "ICC_Machine", "ICC_Cultivar", "ICC_PrevCrop", "ICC_lonlat", "ICC_elevation");
writeLines(capture.output(dft1), paste("wheat_ICC_Model_1.txt", sep=""));

#---;