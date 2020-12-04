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
library(parallel)

### Prepare input data
##(1) Create lattice arena with some spatial autocorrelation
# 2 km squared arena with 250 m resolution (1600 total cells)
gridDim <- 2000
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
set.seed(123); sim_temp <- create_sims(n_reps = 1, n = ndays, central = 11.9, sd = 4,
	exposure_type = "continuous", exposure_trend = "cos1",
	exposure_amp = -0.9, average_outcome = 1,
	outcome_trend = "cos1", outcome_amp = -0.5, 
	rr = 1.005,start.date="2018-01-01")

## T series similar to Barcelona ##
# sim_temp <- create_sims(n_reps = 1, n = ndays, central = 18.2, sd = 2.5, exposure_type = "continuous", exposure_trend = "cos1", exposure_amp = -0.3, average_outcome = 25, outcome_trend = "cos1", outcome_amp = -1, rr = 1.0005,start.date="2019-01-01")
####################################

# Check temperature distribution and what would be the "average" temporal trends for all pixels
hist(sim_temp[[1]]$x)
plot(sim_temp[[1]]$date,sim_temp[[1]]$x)
plot(1:ndays,sim_temp[[1]]$x*mean(sapply(autocorr_factor,max)))
sim_temp[[1]]$month<-format(sim_temp[[1]]$date,"%m")
aggregate(x~month,data=sim_temp[[1]],"mean")

# Compare with some real data
temp <-  readRDS("/home/matteo/own_data/PoD/topics/aedes_genmod/data/t_fresno.RDS")
temp <- temp[temp$yearmoda>="2018-01-01"&temp$yearmoda<="2019-12-31",]
plot(temp$yearmoda,temp$temp, col="red")
abline(v=as.Date("2018-01-01")+121)

sim_temp[[1]]$x<-temp$temp

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
w <- as.matrix(df_temp[1,-c(1:3)]*1000)
#w <- as.matrix(df_temp[,-c(1:3)]*1000)

storage.mode(w) <- "integer"

## Define matrix of coordinates for each cell in the grid
cc <- as.matrix(df_temp[,c("x","y")])
storage.mode(cc) <- "integer"
## Assign colnames from rownames to the distance matrix.
# IMPORTANT: The distance matrix has to be rounded to the thousands 
colnames(dist_matrix) <- row.names(dist_matrix)
dist_matrix <- apply(dist_matrix,2,function(x) round(x/1000,1)*1000) 
#hist(dist_matrix, xlab="Distance (meters)")
storage.mode(dist_matrix) <- "integer"
## Define cells into which introduce propagules
intro.vector <- type.convert(as.numeric(row.names(dist_matrix)))
## Define the day of the year for introduction; day 1 is the first day in the temperature time series
str = as.integer(strftime(as.Date("2018-05-15"),"%j"))
## Define the last day of the simulation; day 1 is the first day in the temperature time series
endr = as.integer(as.Date("2020-07-21") - as.Date("2018-01-01"))
## Define the number of eggs to be introduced
ie = 100; 
## Define number of iterations
it = 5
## Define the number of parallel processes (for sequential itarations set nc=1)
nc = 5
## Set folder where the *.RDS output will be saved
outfolder <- ('./output')
if(!dir.exists(outfolder)){dir.create(outfolder)}
## Set working directory: 
# IMPORTANT: This must be the directory where `DynamAedes.R` is placed
setwd('/home/matteo/own_data/PoD/topics/aedes_genmod/')
#setwd('/home/matteo/GitHub/euaeae/')

## Source the function
source('dynamAedes.r')

species="aegypti"; dispersal=FALSE; temps.matrix=w; cells.coords=cc; lat=0; long=0; road.dist.matrix=dist_matrix; intro.year=2020; startd=1; endd=10; n.clusters=1; cluster.type="SOCK"; iter=1;  intro.cells=NULL; intro.adults=0; intro.immatures=10;  intro.eggs=10; sparse.output=FALSE; compressed.output=FALSE;suffix="dynamAedes"; country="it"

### Run the model
aeg.s <- dynamAedes(species="aegypti", scale="lc", temps.matrix=w, cells.coords=cc, road.dist.matrix=dist_matrix, startd=str, endd=175, n.clusters=8, cluster.type="SOCK",iter=20,intro.cells=intro.vector,intro.eggs=ie, compressed.output=TRUE,country="es",suffix=paste(outfolder,"/dynamAedes_aeg_testrun_dayintro_",str,"_end",endr,"_niters",it,"_neggs",ie,sep=""))

aeg.loc <- dynamAedes(species="aegypti", scale="ws", temps.matrix=w, startd=str, endd=175, n.clusters=8,iter=20, intro.eggs=ie, country="es",suffix=paste(outfolder,"/dynamAedes_aeg_testrun_dayintro_",str,"_end",endr,"_niters",it,"_neggs",ie,sep=""))

albo.s <- dynamAedes(species="albopictus", temps.matrix=w, cells.coords=cc, road.dist.matrix=dist_matrix, startd=str,endd=endr, n.clusters=nc, cluster.type="SOCK",iter=it,intro.cells=intro.vector,intro.eggs=ie, compressed.output=TRUE,country="es",suffix=paste(outfolder,"/dynamAedes_albo_testrun_dayintro_",str,"_end",endr,"_niters",it,"_neggs",ie,sep=""),lat=44.3,long=8.9,intro.year=2018)

kore.s <- dynamAedes(species="koreicus", temps.matrix=w, cells.coords=cc, road.dist.matrix=dist_matrix, startd=str,endd=endr, n.clusters=nc, cluster.type="SOCK",iter=it,intro.cells=intro.vector,intro.eggs=ie, compressed.output=TRUE,country="es",suffix=paste(outfolder,"/dynamAedes_kore_testrun_dayintro_",str,"_end",endr,"_niters",it,"_neggs",ie,sep=""))

#Coords: Bolzano: 46.5,11.3; Amsterdam: 51.6,4.5; Genova: 44.3,8.9; Barcelona: 41,3. Rio: -22.9,-43.4

s <- aeg.s
days <- max(sapply(s, function(x) length(x)))
## Define temperature dates, assuming the two years are 2019 and 2020
tdates <- seq(as.Date("2018-01-01")+str,as.Date("2018-01-01")+endr,by="day")

#%%%%%%%%%%%%%%%%%%%##
### Derive probability of a successfull introduction at the end of the simulated period
pofe <- sum(unlist(mclapply(s, function(x) {
	days <- 365*2
	pe <- length(which(sum(x[days][[1]])>0)) 
	gc()
	return(pe)
},mc.cores=8))) / length(s)

### Print the probability
cat("\n### p of successfull introduction is:",pofe," ###\n")

#%%%%%%%%%%%%%%%%%%%#
## Derive abundance 95% CI for each life stage and in each day
# Define function
dabu95ci <- function(outl=NA,st=1,cores=1,days=0){
	out <- apply(do.call(rbind.data.frame,mclapply(outl, function(x) {
		lapply(1:days, function(y) {
			if(y<=length(x)) {sum(x[[y]][st,],na.rm=T)} else {NA}})},mc.cores=cores)),2,quantile,probs=c(0.25,0.50,0.75),na.rm=T);
	colnames(out)<-NULL
	outo <- rbind.data.frame(out,
		stage=rep(st,nrow(out)),
		day=as.factor(1:ncol(out)))
	return(t(outo))
}
# Apply function and format dataset for ggplot

albo.s <- aeg.s
aeg.s <- aeg.s
kore.s <- aeg.s

dabu_df <- rbind.data.frame(
	rbind.data.frame(
		dabu95ci(albo.s,1,days=max(sapply(albo.s, function(x) length(x)))),dabu95ci(albo.s,2,days=max(sapply(albo.s, function(x) length(x)))),dabu95ci(albo.s,3,days=max(sapply(albo.s, function(x) length(x)))),dabu95ci(albo.s,4,days=max(sapply(albo.s, function(x) length(x))))),
	rbind.data.frame(
		dabu95ci(aeg.s,1,days=max(sapply(aeg.s, function(x) length(x)))),dabu95ci(aeg.s,2,days=max(sapply(aeg.s, function(x) length(x)))),dabu95ci(aeg.s,3,days=max(sapply(aeg.s, function(x) length(x)))), dabu95ci(aeg.s,4,days=max(sapply(aeg.s, function(x) length(x))))),
	rbind.data.frame(
		dabu95ci(kore.s,1,days=max(sapply(kore.s, function(x) length(x)))),dabu95ci(kore.s,2,days=max(sapply(kore.s, function(x) length(x)))),dabu95ci(kore.s,3,days=max(sapply(kore.s, function(x) length(x)))), dabu95ci(kore.s,4,days=max(sapply(kore.s, function(x) length(x)))))
	)

dabu_df$stage <- as.factor(dabu_df$stage)
levels(dabu_df$stage) <- c("Egg","Immature","Adult","Diapause egg")

dabu_df$date <- as.Date(origin=as.Date("2018-01-01")+str,dabu_df$day)
dabu_df$sp <- c(
	rep("albopictus",max(sapply(albo.s, function(x) length(x)))*4),
	rep("aegypti",max(sapply(aeg.s, function(x) length(x)))*4),
	rep("koreicus",max(sapply(kore.s, function(x) length(x)))*4)
	)

# Plot simple abundance time trend
ggplot(dabu_df, aes(y=`50%`,x=date,group=stage,col=stage)) + 
geom_ribbon(aes(ymin=`25%`,ymax=`75%`,fill=stage),col="white",alpha=0.1,outline.type="full") +
geom_line() +
xlab("Date") +
facet_wrap(sp~stage,scale="free_y",ncol=4) +
theme(legend.pos="top") +
ggtitle("Genova - Interquartile abundance per stage") +
scale_x_date(date_breaks = "1 months", date_labels = "%b") +
scale_y_continuous(trans = "log10")

#ggsave("~/dynamaedes_test_genova.png",dpi=400)

sapply(aeg.s,length)

asim <- as.data.frame(t(sapply(aeg.s[[1]],rowSums)))[,1:3]

tdates <- seq(as.Date("2018-01-01")+str,as.Date("2018-01-01")+nrow(asim)+str-1,by="day")
names(asim) <- c("e","i","a")


asim$date <- as.character(tdates)
asim <- melt(asim,id.vars="date")
asim$date <- as.Date(asim$date)

ggplot(asim, aes(x=date,y=value, group=variable, col=variable)) +
geom_line(lwd=1) + 
#facet_wrap(~variable, scale='free_y') +
coord_cartesian(ylim=c(0,2000)) +
#scale_x_date(date_breaks = "2 months", date_labels = "%b-%y")
scale_x_date(date_breaks = "1 months", date_labels = "%b-%y",limits=c(as.Date("2018-11-30"),as.Date("2019-03-31")))


###
