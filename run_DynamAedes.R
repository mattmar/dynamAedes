##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
## R code developed by: 
## Daniele Da Re, Matteo Marcantonio, Diego Montecino
## ---------------------------------------------------##
## Latest update: September 2020
## ---------------------------------------------------##
## Scope: Run DynamAedes.R using i) a simulated 2-year 
## long spatially and temporally correlated temperature  
## time series; ii) simulated road segments; iii) 100 
## eggs as introduction propagule.
## ---------------------------------------------------##
## DOI: 
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
### Load required packages
library(raster)
library(sp)
library(gstat)
library(spatstat)
library(maptools)
library(rgeos)
library(eesim)

### Prepare input data
##(1) Create lattice arena with some spatial autocorrelation
# 10 km squared arena with 250 m resolution (1600 total cells)
gridDim <- 10000
res <- 250
xy <- expand.grid(x=seq(1,gridDim,res), y=seq(1,gridDim,res))
## Create spatial autocorrelation using semivariogram then simple kriging
# Variogram model, with defined sill (value that the semivarion attains at the range) and range (distance of 0 spatial correlation)
varioMod <- vgm(psill=0.005, range=100, model='Exp') # psill = partial sill = (sill-nugget)
# Set up an additional variable from simple kriging
zDummy <- gstat(formula=z~1, locations = ~x+y, dummy=TRUE, 
	beta=1, model=varioMod, nmax=1)
# Generate 2 randomly autocorrelated predictor data fields
set.seed(3)
xyz <- predict(zDummy, newdata=xy, nsim=2)
# Generate an autocorrelated response variable:
# Deterministic on predictors, with normal error
e <- rnorm(nrow(xyz), 18,sd=5)
xyz$resp <- .7*xyz$sim1 + 1.4*xyz$sim2 + e
test <- xyz
# Insert the two variables in the RasterLayer object
gridded(test) <- ~x+y
r <- raster(test[1])
raster::plot(r)
# Parameters for autocorrelation
autocorr_factor <- values(r)
df <- data.frame("id"=1:nrow(test), coordinates(r))
bbox <- as(extent(r), "SpatialPolygons")

##(2) Simulate 2-year temperature data with seasonal trend -----
ndays = 365*2 #[days]
set.seed(123)
sim_temp <- create_sims(n_reps = 1, n = ndays, central = 18, sd = 2,
	exposure_type = "continuous", exposure_trend = "cos1",
	exposure_amp = -.3, average_outcome = 12,
	outcome_trend = "cos1", outcome_amp = 0.8, 
	rr = 1.0005)
# Check temperature distribution and what would be the "average" temporal trends for all pixels
hist(sim_temp[[1]]$x)
plot(1:ndays,sim_temp[[1]]$x*mean(sapply(autocorr_factor,max)))
#Merge spatial autocorrelation to simulated temperatures
mat <- mclapply(1:ncell(r), function(x) {
	d_t <- sim_temp[[1]]$x*autocorr_factor[[x]]
	return(d_t)
},mc.cores=2)
mat <- do.call(rbind,mat)
# Check distribution of autocorrelated temperatures
hist(mat)
# Format temperature data
names(mat) <- paste0("d_", 1:ndays)
df_temp <- cbind(df, mat)

##(3) Simulate an arbitrary road segment for medium-range dispersal 
set.seed(123)
pts <- spsample(bbox, 5, type="random")
roads <- spLines(pts)
# Check simulated segment
raster::plot(r)
raster::plot(roads, add=T)
# Create buffer around the road segment
buff=buffer(roads, width=100)
# Check grid, road segment and buffer
raster::plot(r)
raster::plot(buff, add=T)
raster::plot(roads, add=T, col="red")
# Create a distance matrix for cells inside the buffer along roads
df_sp=df
coordinates(df_sp)=~x+y
df_sp=intersect(df_sp,buff)
# Check selected cells
raster::plot(r)
raster::plot(buff, add=T)
raster::plot(df_sp, add=T, col="red")
# Compute Euclidean distance between each cell in the buffer
dist_matrix <- as.matrix(dist(coordinates(df_sp)))

### Run the model with the simulated input data
## Define temperature matrix
# IMPORTANT: temperature must be multiplied by 1000
w <- df_temp[,-c(1:3)]*1000 
## Define matrix of coordinates for each cell in the grid
cc <- df_temp[,c("x","y")]
## Assign colnames from rownames to the distance matrix.
# IMPORTANT: The distance matrix has to be rounded to the thousands 
colnames(dist_matrix) <- row.names(dist_matrix)
dist_matrix <- apply(dist_matrix,2,function(x) round(x/1000,1)*1000) 
hist(dist_matrix, xlab="Distance (meters)")
## Define cells into which introduce propagules
intro.vector <- as.numeric(row.names(dist_matrix)) 
## Define the day of introduction
str = 121
## Define the end-day of life cycle
endr = 365*2; 
## Define the number of eggs to be introduced
ie = 100; 
## Define number of iterations
it = 5
## Define the number of parallel processes (for sequential itarations set nc=1)
nc = 5
## Set folder where the *.RDS output will be saved
outfolder <- ('./output')
## Set working directory: 
# IMPORTANT: This must be the directory where `DynamAedes.R` is placed
setwd('./')
## Source the function
source('DynamAedes.R')

### Run the model
simout <- DynamAedes(temps.matrix=w, cells.coords=cc, road.dist.matrix=dist_matrix, startd=str,endd=endr, n.clusters=nc, cluster.type="SOCK",iter=it,intro.cells=intro.vector,intro.eggs=ie, compressed.output=TRUE,country="es",suffix=paste(outfolder,"/DynamAedes_testrun_dayintro_",str,"_end",endr,"_niters",it,"_neggs",ie,sep=""))