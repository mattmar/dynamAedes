#Simulate three scale dataset
library(raster)
library(sp)
library(gstat)
library(spatstat)
library(maptools)
library(rgeos)
library(parallel)
library(eesim)

#---- 1 create lattice arena ----
gridDim <- 40 # 10000m/250 m = 40 columns and rows
xy <- expand.grid(x=1:gridDim, y=1:gridDim)
varioMod <- vgm(psill=0.005, range=100, model='Exp') # psill = partial sill = (sill-nugget)
# Set up an additional variable from simple kriging
zDummy <- gstat(formula=z~1, 
                locations = ~x+y, 
                dummy=TRUE,
                beta=1, 
                model=varioMod, 
                nmax=1)
# Generate a randomly autocorrelated predictor data field
set.seed(123)
xyz <- predict(zDummy, newdata=xy, nsim=1)
utm32N <- "+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
r <- raster(nrow=40, ncol=40, crs=utm32N, ext=extent(0,10000, 0,10000))
values(r)=xyz$sim1
plot(r)

df <- data.frame("id"=1:nrow(xyz), coordinates(r))
bbox <- as(extent(r), "SpatialPolygons")

# Store Parameters for autocorrelation
autocorr_factor <- values(r)


#---- 2 simulate temperatures ----
ndays = 365*3 #length of the time series in days
set.seed(123)
sim_temp <- create_sims(n_reps = 1, 
                        n = ndays, 
                        central = 18, 
                        sd = 2, 
                        exposure_type = "continuous", 
                        exposure_trend = "cos1", exposure_amp = -.3, 
                        average_outcome = 12,
                        outcome_trend = "cos1",
                        outcome_amp = 0.8, 
                        rr = 1.0005)

mat <- mclapply(1:ncell(r), function(x) {
  d_t <- sim_temp[[1]]$x*autocorr_factor[[x]]
  return(d_t)
},mc.cores=1) #mc.cores=1 in Windows OS, which does not support mclapply function

mat <- do.call(rbind,mat)

# Format temperature data for local and regional scale
names(mat) <- paste0("d_", 1:ndays)
df_temp <- cbind(df, mat)

#Float numbers in the temperature matrix would slow the computational speed, thus we first multiply them by 1000 and then transform them in integer numbers.
df_temp_loc<- sapply(df_temp[,-c(1:3)], function(x) as.integer(x*1000))

#We can now define a two-column matrix of coordinates to identify each cell in the lattice grid.
cc <- df_temp[,c("x","y")]


# Format temperature data for weather station scale
df_temp_ws=sim_temp[[1]]
df_temp_ws= as.data.frame((t(df_temp_ws[, "x"])))*1000
names_col=paste0("d_", 1:ndays)
names(df_temp_ws)=names_col

#---- 3 simulate roads network ----
set.seed(123)
pts <- spsample(bbox, 5, type="random")
roads <- spLines(pts)
buff <- buffer(roads, width=100)
crs(buff) <- crs(r)
# Check grid, road segment and buffer
raster::plot(r)
raster::plot(buff, add=T)
raster::plot(roads, add=T, col="red")
df_sp=df
coordinates(df_sp)=~x+y
df_sp=intersect(df_sp,buff)

#create dist matrix and format it 
dist_matrix <- as.matrix(dist(coordinates(df_sp)))
colnames(dist_matrix) <- row.names(dist_matrix)
dist_matrix <- apply(dist_matrix,2,function(x) round(x/1000,1)*1000) 

#---- 4 model settings ----
#model parameters
nc=3 #number of cluster
cl=3
it=5 #number of iterations
ie=100 #number of introduced eggs 
str=150 #day into
endr=150+365+50 #day end simulations
outfolder <- ('./output')
## Set working directory: 
# IMPORTANT: This must be the directory where `DynamAedes.R` is placed
setwd('./')

library(parallel)
library(doParallel)
#run the model: weather station scale
simout=dynamAedes(species="aegypti", 
                  scale="ws", 
                  ihwv=1, 
                  verbose = FALSE, 
                  temps.matrix=df_temp_ws, 
                  startd=str, endd=endr, 
                  n.clusters=cl,
                  cluster.type="SOCK", #in Windows OS
                  iter=it, 
                  intro.eggs=ie, 
                  country="it",
                  suffix=paste(outfolder,"/aegypti_","ws", "_dynamAedes_testrun_", "dayintro_",str,"_to_","_enday_",endr,"_niters",it,"_neggs",ie,sep=""))


#run the model: local scale
simout=dynamAedes(species="aegypti", 
                  scale="lc", 
                  cells.coords=cc,
                  road.dist.matrix=dist_matrix,
                  ihwv=10, 
                  verbose = FALSE, 
                  temps.matrix=df_temp_loc, 
                  startd=str, endd=endr, 
                  n.clusters=cl,
                  cluster.type="SOCK", #in Windows OS
                  iter=it, 
                  intro.eggs=ie, 
                  country="it",
                  suffix=paste(outfolder,"/aegypti_","lc", "_dynamAedes_testrun_", "dayintro_",str,"_to_","_enday_",endr,"_niters",it,"_neggs",ie,sep=""))

#run the model: regional scale
simout=dynamAedes(species="aegypti", 
                  scale="rg", 
                  cells.coords=cc,
                  ihwv=100, 
                  verbose = FALSE, 
                  temps.matrix=df_temp_loc, 
                  startd=str, endd=endr, 
                  n.clusters=cl,
                  cluster.type="SOCK", #in Windows OS
                  iter=it, 
                  intro.eggs=ie, 
                  country="it",
                  suffix=paste(outfolder,"/aegypti_","rg", "_dynamAedes_testrun_", "dayintro_",str,"_to_","_enday_",endr,"_niters",it,"_neggs",ie,sep=""))

#---- 5 Analize outputs ----
library(dynamAedes)
source("new_functions_script.R")

#load output
infile=list.files("output/", pattern="rg", full.names = T)
simout=readRDS(infile)

length(simout)#simulations
length(simout[[1]])#days
str(simout[[1]][5])

## 1 get prob succesfull intro
#1.1 overall
my_out=prob_succ_intro(input_sim=simout, eval.date=c(150), st=1) #come lunghezza massima del biennio considerato
my_out

#1.2 spati1al
my_out=prob_succ_intro_spatial(input_sim=simout, eval.date=c(150, 300), coords=cc)
plot(my_out)


## 2 Derive Overall abundance CI for each life stage and in each day

#2.1 Overall
# Derive maxium simulation length in days from run_DynamAedes.R output
days <- max(sapply(simout, function(x) length(x)))

dabu_df <- rbind.data.frame(indiv_abund_ci(simout, st=1,days=days, breaks=c(0.25,0.50,0.75)),
                            indiv_abund_ci(simout, st=2,days=days, breaks=c(0.25,0.50,0.75)),
                            indiv_abund_ci(simout, st=3,days=days, breaks=c(0.25,0.50,0.75) ),
                            indiv_abund_ci(simout, st=4,days=days, breaks=c(0.25,0.50,0.75)))
dabu_df$stage <- as.factor(dabu_df$stage)
levels(dabu_df$stage) <- c("Egg", "Immature", "Adult", "Egg_dia")
head(dabu_df)


##2.2 SPATIAL Derive abundance CI for each life stage and in each day
#output:mappe per due giorni e per uno stadio biologico
myout=indiv_abund_ci_spatial(input_sim=simout, coords = cc, eval.date=c(150, 300), st=1,  breaks=c(0.25,0.50,0.75)) 
plot(myout$d_150)



##3 Derive median and interquartile number of invaded cells per day
icel <- rbind.data.frame(invaded_cells_ci(simout,days=days, breaks=c(0.25,0.50,0.75)))
names(icel) <- c("cells_ic25%","cells_ic50%","cells_ic75%","day")
head(icel)

#get female abundace estimation
res = 250
dfp <- merge(dabu_df[dabu_df$stage%in%"Adult",], icel, by="day")
# dfp[,c(6:8)] = (dfp[,c(6:8)]*res)/40
dfp[,c(6:8)] = dfp[,c(6:8)]*(res^2)/10^6
head(dfp)

plot(dfp$day, dfp$`cells_ic50%`)

##3 spatial rates of mechanistic traits
pippo=get_rates_spatial(coords=cc, 
                        daily_temp =df_temp[,c(5,10,100,250)], 
                        species="aegypti", 
                        rate_fun = ".e.hatch_rate.f", 
                        spatial=FALSE, 
                        rate=FALSE)

pippo=get_rates_spatial(coords=cc, 
                        daily_temp =df_temp[,c(5,10,100,250)], 
                        species="aegypti", 
                        rate_fun = ".e.hatch_rate.f", 
                        spatial=TRUE, 
                        rate=FALSE)
plot(pippo[[1]])

#MISSING: Distance (km) from cell of introduction
