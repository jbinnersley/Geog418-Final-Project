#Libraries
#install.packages("spgwr")
#install.packages("spatstat")
#install.packages("tmap")
#install.packages("gstat")
#install.packages("sf")
#install.packages("raster")
#install.packages("rgdal")
#install.packages("e1071")
#install.packages("spdep")

library(spgwr)
library(spatstat)
library(tmap)
library(gstat)
library(sf)
library(raster)
library(rgdal)
library(e1071)
library(spdep)

#Set working directory
dir <- "C:/Users/JesseB/Desktop/Geography Courses/geog 418/Labs/Final/Jesse Binnersley/Census"
setwd(dir)

#Reading in particulate matter dataset
#Read in PM2.5 data:
pm2.5 <- readOGR(".", "Pm25Sample") 
pm2.5 <- spTransform(pm2.5, CRS("+init=epsg:26910"))

#Reading in dissemination tract and income data
#Read in census income data:
income <- read.csv("./Income.csv")  
#Select only ID and Income columns:
colnames(income) <- c("DAUID", "Income") 
#Read in dissemination tract shapefile:
census.tracts <- readOGR(".", "BC_DA") 
#Merge income and dissemination data:
income.tracts <- merge(census.tracts,income, by = "DAUID") 
#Determine the number of columns in the dataframe:
nrow(income.tracts)
#Remove NA values:
income.tracts <- income.tracts[!is.na(income.tracts$Income),]
#Reproject the data:
income.tracts <- spTransform(income.tracts, CRS("+init=epsg:26910"))

#Create choropleth map of income:

map_Income <- tm_shape(income.tracts) +
  tm_polygons(col = "Income",
              title = "Median Income",
              style = "jenks",
              palette = "viridis", n = 6) +
  tm_legend(legend.position = c("LEFT", "BOTTOM"))

map_Income


#Descriptive Stats

meanIncome <- mean(income.tracts$Income, na.rm = TRUE)
sdIncome <- sd(income.tracts$Income, na.rm = TRUE)
modePop <- as.numeric(names(sort(table(income.tracts$Income), decreasing = TRUE))[1])
#Median
medPop <- median(income.tracts$Income, na.rm = TRUE)
#Skewness
skewPop <- skewness(income.tracts$Income, na.rm = TRUE)[1]
#Kurtosis
kurtPop <- kurtosis(income.tracts$Income, na.rm = TRUE)[1]
#CoV
CoVPop <- (sdIncome / meanIncome) * 100
#TotalArea <- sum(df$Income)
#Normal distribution test
normIncome_PVAL <- shapiro.test(income.tracts$Income)$p.value



##Global Morans I Test Rooks Case
income.nb <- poly2nb(income.tracts)
income.net <- nb2lines(income.nb, coords=coordinates(income.tracts))
crs(income.net) <- crs(income.tracts)

income.lw <- nb2listw(income.nb, zero.policy = TRUE, style = "W")
print.listw(income.lw, zero.policy = TRUE)

mi <- moran.test(income.tracts$Income, income.lw, zero.policy = TRUE)
mi

moran.range <- function(lw) {
  wmat <- listw2mat(lw)
  return(range(eigen((wmat + t(wmat))/2)$values))
}
moran.range(income.lw)

mI <- mi$estimate[[1]]
eI <- mi$estimate[[2]]
var <- mi$estimate[[3]]

z <- (mI-eI)/sqrt(var)

#Local Morans I
lisa.test <- localmoran(income.tracts$Income, income.lw)

income.tracts$Ii <- lisa.test[,1]
income.tracts$E.Ii<- lisa.test[,2]
income.tracts$Var.Ii<- lisa.test[,3]
income.tracts$Z.Ii<- lisa.test[,4]
income.tracts$P<- lisa.test[,5]

map_LISA <- tm_shape(income.tracts) + 
  tm_polygons(col = "Z.Ii", 
              title = "Local Moran's I", 
              style = "jenks", 
              palette = "viridis", n = 6) 


map_LISA


#Spatial Interpolation

#Create a grid called grd to use in your interpolation
# Create an empty grid where n is the total number of cells
grd <- as.data.frame(spsample(pm2.5, "regular", n=50000))
names(grd)       <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
# Create SpatialPixel object:
gridded(grd)     <- TRUE  
# Create SpatialGrid object:
fullgrid(grd)    <- TRUE  
#Reproject the grid:
proj4string(grd) <- proj4string(income.tracts)



P.idw <- gstat::idw(PM25 ~ 1, pm2.5, newdata=grd, idp=2)
r       <- raster(P.idw)
r.m     <- mask(r, income.tracts)

surfaceMap <- tm_shape(r.m) + 
  tm_raster(n=5,palette = "viridis",
            title="Predicted pollution \n(PM2.5)") + 
  tm_shape(pm2.5) + tm_dots(size=0.05) +
  tm_legend(legend.outside=TRUE)
surfaceMap

income.tracts$Pm2.5 <- round(extract(r, income.tracts, fun = mean)[,1], 5)

######Linear Regression##########
#Let's say your dataset with both PM2.5 and Income 
#are stored in a dataset called income.tracts.
#Plot income and PM2.5 from the income.tracts dataset you created
plot(income.tracts$Income~income.tracts$Pm2.5)

#Notice that there are a lot of 0's in this dataset. If you decide to remove them, use the following line:
income.tracts.no0 <-  income.tracts[which(income.tracts$Pm2.5 > 0), ]

#Now plot the data again
plot(income.tracts.no0$Income~income.tracts.no0$Pm2.5)

#Perform a linear regression on the two variables. You should decide which one is dependent.
lm.model <- lm(income.tracts.no0$Income~income.tracts.no0$Pm2.5)
#Add the regression model to the plot you created
plot(income.tracts.no0$Income~income.tracts.no0$Pm2.5)
abline(lm.model, col = "red")
#Get the summary of the results
summary(lm.model)

#add the fitted values to your spatialpolygon dataframe
income.tracts.no0$predictlm <- lm.model$fitted.values

#You want to determine if the model residuals are spatially clustered. 
#add the residuals to your spatialpolygon dataframe
income.tracts.no0$residuals <- residuals.lm(lm.model)

#Observe the result to make sure it looks correct
head(income.tracts.no0)

#Now, create choropleth map of residuals
map_resid <- tm_shape(income.tracts.no0) +
  tm_polygons(col = "residuals",
              title = "Median Income",
              style = "jenks",
              palette = "viridis", n = 6,
              midpoint = NA) 

map_resid

##Global Morans I Test Rooks Case
residuals.nb <- poly2nb(income.tracts.no0)
residuals.net <- nb2lines(residuals.nb, coords=coordinates(income.tracts.no0))
crs(residuals.net) <- crs(income.tracts.no0)

residuals.lw <- nb2listw(residuals.nb, zero.policy = TRUE, style = "W")
print.listw(residuals.lw, zero.policy = TRUE)

mi <- moran.test(income.tracts.no0$residuals, residuals.lw, zero.policy = TRUE)
mi

moran.range <- function(lw) {
  wmat <- listw2mat(lw)
  return(range(eigen((wmat + t(wmat))/2)$values))
}
moran.range(residuals.lw)

mI <- mi$estimate[[1]]
eI <- mi$estimate[[2]]
var <- mi$estimate[[3]]

z <- (mI-eI)/sqrt(var)


#Point Pattern Analysis (lab2)


kma <- pm2.5
kma$x <- coordinates(kma)[,1]
kma$y <- coordinates(kma)[,2]


#check for and remove duplicated points
#first, finds zero distance among points to see if there are any duplicates
zd <- zerodist(kma)

#if there are duplicates, remove them
kma <- remove.duplicates(kma)

#create an "extent" object which can be used to create the observation window for spatstat
kma.ext <- as.matrix(extent(kma)) 

#observation window
window <- as.owin(list(xrange = kma.ext[1,], yrange = kma.ext[2,]))

#create ppp oject from spatstat
kma.ppp <- ppp(x = kma$x, y = kma$y, window = window)



#####
##QUADRAT ANALYSIS
##First, determine the number of qusdrats 
quads <- 10

qcount <- quadratcount(kma.ppp, nx = quads, ny = quads)


plot(kma.ppp, pch = "+", cex = 0.5)
plot(qcount, add = T, col = "red", title = "TITLE")

qcount.df <- as.data.frame(qcount)

##Second, count the number of quadrats with a distinct number of points.
qcount.df <- plyr::count(qcount.df,'Freq')

##Change the column names so that x=number of points and f=frequency of quadrats with x point.
colnames(qcount.df) <- c("x","f")

sum.f.x2 <- sum(qcount.df$x * qcount.df$f * qcount.df$x)
M <- quads*quads
N <- sum(qcount.df$f * qcount.df$x)
sum.fx.2 <- sum(qcount.df$x*qcount.df$f) * sum(qcount.df$x*qcount.df$f)

VAR <- (sum.f.x2-(sum.fx.2/M))/(M-1)
MEAN <- N/M
VMR <- VAR/MEAN

##Finally, perform the test statistic to test for the existence of a random spatial pattern.
chi.square = VMR * (M - 1)
p = 1 - pchisq(chi.square, (M - 1))


###############
income.tracts.no0.coords <- sp::coordinates(income.tracts.no0)
#Observe the result:
head(income.tracts.no0.coords)
#Now add the coordinates back to the spatialpolygondataframe
income.tracts.no0$X <- income.tracts.no0.coords[,1]
income.tracts.no0$Y <- income.tracts.no0.coords[,2]

###Determine the bandwidth for GWR: this will take a while
GWRbandwidth <- gwr.sel(income.tracts.no0$Income~income.tracts.no0$Pm2.5, 
                        data=income.tracts.no0, coords=cbind(income.tracts.no0$X,income.tracts.no0$Y),adapt=T) 

###Perform GWR on the two variables with the bandwidth determined above
###This will take a looooooong while
gwr.model = gwr(income.tracts.no0$Income~income.tracts.no0$Pm2.5, 
                data=income.tracts.no0, coords=cbind(income.tracts.no0$X,income.tracts.no0$Y), 
                adapt=GWRbandwidth, hatmatrix=TRUE, se.fit=TRUE) 

#Print the results of the model
gwr.model

#Look at the results in detail
results<-as.data.frame(gwr.model$SDF)
head(results)

#Now for the magic. Let's add our local r-square values to the map
income.tracts.no0$localr <- results$localR2

#Create choropleth map of r-square values
map_r2 <- tm_shape(income.tracts.no0) +
  tm_polygons(col = "localr",
              title = "R2 values",
              style = "jenks",
              palette = "viridis", n = 6)
map_r2

#Time for more magic. Let's map the coefficients
income.tracts.no0$coeff <- results$income.tracts.no0.Pm2.5
#Create choropleth map of the coefficients
map_coef <- tm_shape(income.tracts.no0) +
  tm_polygons(col = "coeff",
              title = "Coefficients",
              style = "jenks",
              palette = "viridis", n = 6,
              midpoint = NA)
map_coef
