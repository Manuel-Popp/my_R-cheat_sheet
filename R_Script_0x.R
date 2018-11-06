# Useful commands in R ~ using RStudio R version 3.3.2 (2016-10-31)
########################################################
# ~ Contents ~
#  0 - Installing multiple packages if required
#  1 - Setting the working directory; finding the script location and adapting the wd
#  2 - downloading data from a website
#  3 - creating tables and textfiles a) .csv; b) .csv2; c) .rds d) .txt
#  4 - reading tables a) .csv; b) .csv2; c) .rds; d) .txt; e) .ods
#  5 - clear workspace
#  6 - creating Pie plot a) w/ ggplot2; b) w/ base R
#  7 - creating Barplot
#  8 - creating Histogr?m a) w/ base R; b) w/ ggplot
#  9 - calculating mean and median
# 10 - calculating the standard deviation
# 11 - calculating the coefficient of variation
# 12 - calculating quartile distance
# 13 - calculating the Range
# 14 - How to deal with NA values
# 15 - Stacked Area Plot w/ggplot2
# 16 - mapping plots/ plot on map w/ base R
# 17 - ggplot2 +additionals
# 18 - load and edit a map w/ ggmap
# 19 - create grob from a map or plot
# 20 - combining multiple plots on one layer
# 21 - t-test w/ R
# 22 - correlation (Spearman, Bravais-Pearson (Pearson's r), Kendall)
# 23 - correlation: strength and significance
# 24 - create maps
# 25 - combine data
########################################################

########################################################
#0. Installing multiple packages if required
########################################################
packages <- c("classInt",
              "rgdal",
              "raster"
              )
for(i in 1:NROW(packages)){
  if (!require(packages[i], character.only= TRUE)) {
    install.packages(packages[i])
    library(packages[i], character.only= TRUE)
  }
}

########################################################
#1. set working directory and creade folder surroundings
########################################################
if (!require("rstudioapi")) install.packages("rstudioapi")
library(rstudioapi)    
file_dir <- rstudioapi::getActiveDocumentContext()$path
r_script <- gsub(basename(file_dir), "", file_dir)

setwd(r_script)
setwd("..")
wd <- getwd()

dir.create(paste(wd, "/data", sep=""))

# 1 a) moving script to file
#   dir.create(paste(wd, "/r_script", sep=""))
#   file.rename(paste(r_script, "/R_Script_1.R", sep=""),
#            paste(wd, "/r_script/R_Script_1.R", sep=""))
#   @me: How tell RStudio that the file has been moved to avoid the "File moved or deleted, wanna close window now?" notification?


########################################################
#2. downloading data
########################################################
#   unfortunately, this does not work for ILIAS. (For dropbox it can be accomplished by simply setting dl=0 to dl=1 when the file is not password protected)
... <- download.file("https://www.dropbox.com/s/i5y2qn0ylr52dlf/data_example.csv?dl=1",
                        paste(wd, "/data", "/data_example.csv", sep="")
              )
#alternative w/ password protected site:
if (!require("httr")) install.packages("httr")
... <- GET("https://ilias.studium.kit.edu/goto.php?target=file_772770_download&client_id=produktiv", 
       authenticate("user", "psswrd"))

########################################################
#3. creating tables a) .csv b) .csv2 c) .rds d) .txt
########################################################

path_1 <- paste(
  wd, "/data/examples", sep=""
)
dir.create(path_1)

# a) .csv
write.csv2(
  example_data, file= paste(
    path_1, "/example_1.csv", sep=""
    )
  )

# b) .csv2
write.table(
  example_data, file= paste(
    path_1, "/example_2.csv", sep=""
  )
)

# c) .rds
saveRDS(
  example_data, file= paste(
    path_1, "/example_3.rds", sep=""
  )
)

# d) .txt
write.table(
  citation, paste(
    wd, "/info", "/citation.txt", sep=""
    ), sep=""
  )

# e) .ods need package: "readODF"
read_ods(
  paste(
    wd, "/data", "/name.ods", sep=""
  )
)
########################################################
#4. reading tables a) .csv; b) .csv2; c) .rds
########################################################
#a) .csv
... <- read.table(
  paste(wd, "/data/example/example_1.csv", sep=""
  ), sep=";", header= TRUE
)

#b) .csv2 (notice: ending is just .csv though)
... <- read.table(
  paste(wd, "/data/example/example_2.csv", sep=""
  ), sep="", header= TRUE
)

#c) .rds
... <- readRDS(
  paste(wd, "/data/example/example_3.rds", sep=""
  )
)
########################################################
#5. clear workspace (I saved wd here though)
########################################################
setwd(wd)
rm(list=ls())
wd <- getwd()
########################################################
#6. create Pie plot a) w/ ggplot2; b) w/ base R
########################################################
#a) ggplot2 Pie Plot
example_data <- read.table(
  paste(wd, "/data/data_example.csv", sep=""
  ), sep=",", header= TRUE
)

if (!require("ggplot2")) install.packages("ggplot2")
if (!require("reshape2")) install.packages("reshape2")
if (!require("data.table")) install.packages("data.table")
temp <- example_data[, c(7, 8)]
agg <- aggregate(temp, by= list(temp$formation), FUN=mean, na.rm=TRUE)
agg <- agg[, c(1, 3)]
names(agg) <- c("formation", "number of species")

library(reshape2)
expl_melted <- melt(agg, id.vars= "formation")

library(data.table)

setnames(expl_melted, old = names(expl_melted), new = c("formation", "variable", "number of species"))

library(ggplot)
pie <- ggplot(expl_melted, aes(x= "", fill= formation)) +
  geom_bar(width = 1) +
  coord_polar("y", start=0) +
  labs(title = "habitats of the species found in the example project", fill= "") +
  xlab("") +
  ylab("") +
  labs(fill = "habitat") +
  scale_fill_manual(values = c("forestgreen", "lightgreen", "chocolate4"))
pie

#b) base R Pie Plot
pie <- pie(table(example_data$formation), col = c("forestgreen", "lightgreen", "chocolate4"))
########################################################
#7. create Barplot
########################################################
bp <- barplot(table(example_data$formation), col = c("forestgreen", "lightgreen", "chocolate4"))

########################################################
#8. create Histogram a) w/ base R; b) w/ ggplot
########################################################
#a) w/ base R
histo <- hist(example_data[,4],
     main="temperature histogram", 
     xlab="temperature [?C]", 
     border="black", 
     col="firebrick1")
# or
h <- hist(example_data$temperature, plot=FALSE, breaks=10)
h$counts <- cumsum(h$counts)
h$density <- cumsum(h$density)

plot(h, freq=TRUE, main="(Cumulative) histogram of Temperature", col="gray", border="white")
box()

# alternative
plot(h$mids, h$density, type="b", pch=16, xlab="Temperature (°C)", ylab="Cumm. density")

#b) w/ ggplot
#   Here dat= data table, Aug= values for August, ggdat[,2]= dat$Aug, ggdat[,3]= variable
#   (in this example there are 2 histograms plotted on one another to compare past and current temperatures in August)

if (!require("ggplot2")) install.packages("ggplot2")
if (!require("reshape2")) install.packages("reshape2")

ggmin <- round(min(dat$Aug), digits= 1)
ggmax <- round(max(dat$Aug), digits= 1)
ggbrks <- seq(ggmin, ggmax, len= ((ggmax-ggmin)*5-.5)) #writes a sequence of numbers each differing by 0.5 from min to max
gg_dat <- dat[,c(1,9, 20)] #combines the columns "Year", "Aug" and "past/current" from the original data table

ggplot(gg_dat, aes(gg_dat[,2], fill= gg_dat[,3])) + geom_density(alpha= 0.4) +
  scale_fill_manual(values = c("forestgreen", "firebrick1")) +
  labs(title= "Histogram: Temperature August", fill = "Timeline") +
  scale_x_continuous(name= "Temperature", breaks = ggbrks, limits=c(ggmin, ggmax)) +
  ylab("Frequency") +
  theme_minimal() +
  theme(plot.title= element_text(hjust= 0.5, size= rel(1.6)),
        axis.title= element_text(size= rel(1.3)),
        legend.title= element_text(size = rel((1.2))),
        legend.text= element_text(size = rel((1.1))),
        axis.text = element_text(size = rel((1.1)))
  )

########################################################
#9. Calculate mean and median
########################################################
mean_temp <- mean(example_data$temperature)
median_temp <- median(example_data$temperature)

########################################################
#10. Calculate the standard deviation
########################################################
standard_dev <- sd(example_data$temperature, na.rm= TRUE)

########################################################
#11 Calculate the coefficient of variation
########################################################
coeff_var <- standard_dev/mean_temp*100

########################################################
#12 Calculate quartile distance
########################################################
quart_dist <- IQR(example_data$temperature)

########################################################
#13 Calculate the Range
########################################################
range <-range(example_data$temperature)

########################################################
#14 How to deal with NA values
########################################################
example_data$temperature[c(4, 7)] <- NA

### It is not possible to calculate with NA values! This is why I "#" them out - I was said we must not create errors in our R script
#mean(dat$temperature)
#mean(dat$temperature, na.rm=T)

### Find NAs
#example_data$temperature == NA ### Does not help!
is.na(example_data$temperature)

### Exclude NA values from vector
example_data$temperature[!is.na(example_data$temperature)]

### Exclude all rows with NAs! Sometimes this makes sense, sometimes not.
dat1 <- na.omit(example_data)
summary(dat1)
nrow(example_data)
nrow(dat1)

########################################################
#15 Stacked Area Plot w/ggplot2
########################################################
require(ggplot2)
library(jpeg)

jpeg(filename= "past_disturbance.jpg")
plot_1 <- ggplot(data= pastPlot_melted, aes_(x= pastPlot_melted$year, y= pastPlot_melted$total,
                                             fill= pastPlot_melted$disturbance)) +
  geom_area(alpha=0.8, mapping = NULL, stat = "identity",
            position = "stack", na.rm = FALSE, show.legend = NA, inherit.aes = TRUE
  ) + labs(title = "overall past disturbances", fill= "") +
  xlab("Year") + ylab("disturbance damage [ha timber/Year]") +
  scale_fill_manual(values= c("deepskyblue", "black", "red"))
print(plot_1) #send data to jpeg device
dev.off()
print(plot_1) #show plot

setwd("..")

########################################################
#16 mapping plots/ plot on map w/ base R
########################################################
require(rworldmap)
require(mapproj)
require(TeachingDemos)

grid.newpage()
setwd(paste(wd,"/plot", sep=""))

plot(getMap(), xlim = c(-5, 35), ylim = c(30, 70), asp = 1)
for(i in 1:8){
  subplot(
    barplot(
      height = unlist(Example_agg[i, 2:4], use.names= FALSE), width = c(.5, .5, .5),
      axes=FALSE,
      col= c("blue", "black", "red"), ylim= range(Example_agg[,2:4])),
    x=example_agg[i, 5], y=example_agg[i, 6], vadj = 0, size=c(.5, 1)
  )
}
addMapLegendBoxes(colourVector = c("blue", "black", "red"), x = "topleft", cutVector = c("variable_1", "variable_2", "variable_3"),
                  cex = 1, pt.cex = 2)
plot_expl <- recordPlot()
jpeg(filename = "example.jpg", quality = 100, width = mapsize[1], height = mapsize[2])
print(plot_expl)
dev.off()

########################################################
#17 editing in ggplot2
########################################################
#basic
ggplot(data, aes_(x= X_Axis_label)) +
  theme_light() +
  geom_bar()

#split the data into several plots by category
+
  facet_wrap(~ category)

#split the data into several plots by two dimensions (by two categories)
+
  facet_wrap(category_1 ~ category_2)

#use fill to divide bars into two colours by percentage for category
ggplot(data, aes(x= x_axis_label, fill= category))

# add a line
+ 
  geom_line(alpha=0.8, mapping = NULL, stat = "identity",
            position = "identity", na.rm = FALSE, show.legend = NA, inherit.aes = TRUE)
# add points (to a line)
+
  geom_point(size = 3, inherit.aes = TRUE,
             position = "identity", na.rm = FALSE)
# edit labs titles
+
  labs(title = "Title")
+
  xlab("x_value")
+
  ylab("variable_value")
# edit colors
+
  scale_colour_manual(values= c("deepskyblue", "black", "red"))
# remove legend title
+
  theme(legend.title=element_blank())

#histogram, binwidth= x provides the option for scaling the x axis resolution by grouping
#by intervals of the size x
ggplot(data, aes(x= x_axis_label)) +
  geom_histogram(binwidth= 1)

########################################################
#18 load a map w/ ggmap
########################################################
if (!require("ggmap")) install.packages("ggmap")

mapsize <- c(640, round(as.numeric(640*dev.size("cm")[2]/dev.size("cm")[1])))
map <- get_googlemap(center = c(lon = 7.4446, lat = 46.9479), zoom = 4,
                     maptype = "satellite",
                     size = mapsize, 
                     style = c("feature:administrative.country|element:labels|visibility:off",
                               "feature:administrative.city|element:labels|visibility:off")
)     ##,"feature:administrative.country|element:geometry.fill|lightness:90" doesn't work!!
ggmap(map, extent = "device")

########################################################
#19 create grob from map or plot
########################################################
library(gridExtra)
library(EBImage)
library(RGraphics)

#create grob from "back_map"
grid.newpage()
mapsize <- c(640, round(as.numeric(640*dev.size("cm")[2]/dev.size("cm")[1])))
mapratio <- mapsize[1]/mapsize[2]
gob <- rasterGrob(unclass(back_map), width = 1, height = mapratio, just = "center")
grid.draw(gob)
#   create white shade to lay over background map because google's map style option "lightness" seems to work for
#   elements only, not for the map itself
shade <- ggplotGrob(ggplot() + theme(panel.background = element_rect(fill = alpha("white", 0.5), color = "white")))

#create grob from "plot_combined"
g <- ggplotGrob(plot_combined)

########################################################
#20 combining plots
########################################################
if (!require("cowplot")) install.packages("cowplot")
library(cowplot)

Theme <- theme(axis.text=element_text(size=6), axis.title=element_text(size=6,face="bold"),
               plot.title =element_text(size=8, face='bold'), legend.position="none")

plot_combined <- plot_grid(plot_1 + Theme, plot_2 + Theme, get_legend(plot_1),
                           plot_3 + Theme, plot_4 + Theme, plot_5 + Theme,
                           plot_6 + Theme, plot_7 + Theme, plot_8 + Theme, nrow=3, ncol=3)
print(plot_combined)

#for standard R Plot
par(mfrow=c(1,2))

########################################################
#21 t-test w/ R
########################################################

#single variable t-test
t.test(y~x)                   #where y is numeric and x is a binary factor
                              #(x can be the character vector defining the two groups being compared e.g "m"/"w")
                              #Used to test whether the mean is actually around a given value. e.g. "Is the
                              #price I paied for my car really the average market price for this car?"
                              #Calculates s*? which is an estimated sigma? (=variance) based on the sample values.
#paired t-test
t.test(y1, y2)                #where y1 and y2 are numeric, independent 2-group t-test
                              #Test of two independent groups. e.g. "Do men pay as much as women when
                              #shopping for shoes?". Tests whether the mean values differ between the groups
                              #or how probable is it, that the average values of the respective groups differ.

t.test(y1, y2, paired= TRUE)  #where y1 & y2 are numeric and the objects are linked to each other
                              #e.g. testing a group of people for their alcohol level before and after
                              #a party. Each person is tested twice.
                              #Tests how big the probability is, that there is a difference between the first
                              #and the second test (with one parameter changed between the 1st and 2nd test).

#   The alpha level applied to a t-test represents the probability of rejecting the null hypothesis
#   (= assumtion, that both, y1 and y2 are from the same basic population), though it is true. When
#   the p-value, which represents the possibility that the null hypothesis is true and the differences
#   between y1 and y2 are based on random chance, is smaller than the alpha level, the null hypothesis
#   can be regarded as falsified and there is strong evidence, that the two samples y1 and y2 differ.
#   As opposite to the alpha level, the beta level represents the probability, that the null hypothesis
#   is not rejected, however it is false. A common alpha level to work with is 0.05.

#   When running a t-test, the following assumptions are made:
#   The data is collected from a representative, randomly selected portion of the total polulation
#   and the scale of measurement applied follows a metric or ordinal scale.
#   Furthermore, it is necessary that the data shows a normal distribution and that the sample is of
#   a reasonalbe size, so the distribution should approach a normal bell-shaped curve. Also, when
#   doing a two sample t-test, the variance of the samples must be homogenous, which means the sd of
#   the respective groups must be approximately eaqual.

########################################################
#22 Correlation (Spearman, Bravais-Pearson, Kendall)
########################################################
spearman <- cor(dat_$x, dat_$y, method= "spearman")
#   [-1,1]; 1= total positive linear correlation; 0= no linear correlation; -1= total negative linear correlation
#   Data must be at least ordinal and the scores on one variable must be monotonically related to the other variable
pearson <- cor(dat$x, dat$y, method= "pearson")
#   [-1,1]; 1= total positive linear correlation; 0= no linear correlation; ???1= total negative linear correlation
#   Variables have to be continuous; both variables should be normally distributed; linear relationship suggested;
#   homoscedasticity of the data (data is equally distributed about the regression line)
#
#   I          . .        I            .: '.'
#   I        . .          I          . .. ...
#   I      . .            I       . . . ..' .
#   I    . .              I   . . ..  
#   I  . .                I  .        
#   I. .                  I.
#   I._________________   I _________________
#   variation const        variation increasing= no homoscedasticity
#
kendall <- cor(dat$x, dat$y, method= "kendall")

########################################################
#23 Correlation: strength and significance
########################################################
#Test for association between paired samples, using one of Pearson's product moment correlation coefficient,
#Kendall's \(\tau\) or Spearman's \(\rho\)
cor.test(dat$x, dat$y, alternative = c("two.sided", "less", "greater"),
         method = c("pearson", "kendall", "spearman"),
         exact = NULL, conf.level = 0.95)

#Gives out Pearson's r as well as the p-value, the t-value and whether the null hypothesis or the alternative hypothesis
#is the more probable one, depending on the significance level (conv.level) S= 1-alpha with alpha = max. accepted proba-
#bility of a type 1 error. The smaller alpha, the higher is the confidence level and the higher the conv. level, the more
#likely is a type 2 error and the less likely is a type 1 error. That means the lower alpha the less likely one will
#consider the null hypothesis false when it is true. A bigger alpha makes it more probable that the alternative hypothesis
#is regarded as the more probable one despite it actually is the wrong one.

#devtools
if(!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("dkahle/ggmap", ref = "tidyup")

########################################################
#24 Creating maps
########################################################
#administrative boundaries
getData('GADM', country = 'DE', level = 1)
# Level decides the scale (0 = Country, 1 = Bundesland, 2 = Kreis, 3 = Stadt, 4 = Siedlung)

# DEM (elevation data)
getData('SRTM', lon=5, lat=45)
#or
getData('alt', country = "DE")

### examples: note on which levels country, states, cities etc are
mycolors <- c("blue", "red", "green", "black")
shp <- getData("GADM", country = 'DE', level = 4, path = "../")
BaWu <- shp[shp@data$NAME_1 == "Baden-Württemberg",]
ka <- BaWu[c(
  which(BaWu@data$NAME_4 == "Bruchsal"))
  ,]
plot(ka)
points(latitude ~ longitude, data = sites, col = mycolors[sites$vegclass], pch = 16, asp = 1)
legend(
  "topleft", legend = levels(
    as.factor(sites$vegclass)
  ), col = mycolors, pch=16, title = "vegetation class"
)


plot_01 <- ggplot(ka, aes(x = long, y = lat)) +
  geom_path(aes(group = group)) +
  geom_point(data = sites, aes(x = sites$longitude, y = sites$latitude, col = sites$vegclass)) +
  labs(x = "longitude", y = "latitude") +
  theme(legend.position = "topright")

plot_01

# plot DEM instead of administrative boundaries
DEM <- getData("alt", country = "DE")
plot(
  DEM, xlim = c(
    min(sites$longitude - .05, na.rm = TRUE), max(sites$longitude + .05, na.rm = TRUE)
  ),
  ylim = c(
    min(sites$latitude - .01, na.rm = TRUE), max(sites$latitude + .01, na.rm = TRUE)
  )
)
points(latitude ~ longitude, data = sites, col = mycolors[sites$vegclass], pch = 16, asp = 1)

# using ggplot
dem <- rasterToPoints(DEM)
df <-  data.frame(dem)
colnames(df) = c("lon", "lat", "alt")
dfx <- subset(df, df$lon >= min(sites$longitude - .05, na.rm = TRUE) &
                df$lon <= max(sites$longitude + .05, na.rm = TRUE) &
                df$lat >= min(sites$latitude - .01, na.rm = TRUE) &
                df$lat <= max(sites$latitude + .01, na.rm = TRUE)
)

plot_01 <- ggplot(dfx, aes(x = lon, y = lat)) +
  geom_raster(aes(fill = alt)) +
  scale_fill_gradient(low = "black", high = "white") +
  geom_point(data = sites, aes(x = sites$longitude, y = sites$latitude, col = sites$vegclass))

plot_01

########################################################
#25 combine data
########################################################
# When combining two sets of data, it is important whether the entries are in the same order in both sources.
# the match() command tells us at which place in the second data set the entry with a certain name is stored.
# It gives out the row in which the entries must be ordered to make both sources fit together.
m <- match(df1$ID, df2$ID)
df1$additional_data <- df2[m]
df1[!is.na(df1$additional_data),] # remove the rows where no data was added, hence the rows that exist only in one data frame
