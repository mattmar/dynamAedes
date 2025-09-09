#' Life cycle simulation of \emph{Aedes} mosquitoes
#'
#' Function to simulate population dynamics of \emph{Aedes} mosquitoes
#' @param species character. Select what species to model: \code{"aegypti"}, \code{"albopictus"}, \code{"japonicus"}, \code{"koreicus"}. Default \code{species = "aegypti"}. 
#' @param intro.eggs positive integer. number of introduced eggs, default \code{intro.eggs = 0}.
#' @param intro.deggs positive integer. number of introduced diapause eggs, default \code{intro.deggs = 100}.
#' @param intro.adults positive integer. number of introduced adults, default \code{intro.adults = 0}. 
#' @param intro.juveniles positive integer. number of introduced juveniles, default \code{intro.juveniles = 0}.
#' @param scale character. Define the model spatial scale: punctual/weather station "ws", local "lc", or regional "rg". Active and passive dispersal is enabled only for \code{scale = "lc"}. Default \code{scale = "ws"}.
#' @param intro.cells positive integer. One or more cells (id) where to introduce the population at local ("lc") scale. If intro.cells=NULL, then a random cell is used for introduction; If intro.cells is a vector of cell ids then a cell is drawn at random from the vector (with repetition) for introduction in each model iteration. 
#' @param jhwv positive integer. Juvenile-habitat water volume, define the volume (L) of water habitat presents in each spatial unit (parametrised with data retrieved from \doi{10.1111/1365-2664.12620}). Default \code{lhwv = 1}.
#' @param startd Character date (ISO format "%Y-%m-%d"). Date of start of simulations.
#' @param endd Character date (ISO format "%Y-%m-%d"). Date of end of simulation. It can be \code{NA}; then it will be derived using the number of columns in \code{temps.matrix}.
#' @param iter positive integer. Define the number of model iterations. 
#' @param temps.matrix matrix. A matrix of daily (average) temperatures (in degrees **Celsius degree x 1000**) used to fit the life cycle rates. This matrix must be organised with the daily temperature observations as columns and the geographic position of the i-grid cell as rows. \bold{Importantly}, the first column must match \code{startd} date.
#' @param cells.coords matrix. A matrix reporting the spatial coordinates of the temperature observations.
#' @param coords.proj4 string. Proj4 string of cell coordinates used for the calculation of photoperiod.
#' @param lat numeric. Latitude value of the area of interested used to derive the photoperiod (and thus the diapause eggs allocation function). 
#' @param long numeric. Longitude value of the area of interested used to derive the photoperiod (and thus the diapause eggs allocation function)
#' @param road.dist.matrix matrix. when \code{scale = "lc"}, defines the matrix containing the distances (in meters) between grid cells intersecting the road network for the mosquito passive dispersal process. colnames(road.dist.matrix) must correspond to the same cell in temps.matrix (may be an issue if the coordinates used for creating road.dist.matrix have been subset to match road segments).
#' @param avgpdisp optional. when \code{scale = "lc"}, define the average car trip distance (in meters) for the mosquito passive dispersal process. The value can be set by the users (positive numeric), or the estimates made by \href{https://publications.jrc.ec.europa.eu/repository/handle/JRC77079}{Pasaoglu et al. 2012}) for the following European countries: France "fra", Germany "deu", Italy "ita", Poland "pol", Spain "esp" and the United Kingdom "uk". The average passive dispersal distance must be smaller than the maximum distance in **road.dist.matrix**.
#' @param pDispersal boolean. if TRUE (default) then the model considers passive dispersal.
#' @param cellsize positive integer. When \code{scale = "lc"}, defines the minimal distance for the active dispersal kernel and should match the spatial resolution of temps.matrix to avoid inconsistencies. Default cellsize = 250.
#' @param maxadisp positive integer. When \code{scale = "lc"}, defines the maximum daily active dispersal, default maxadisp = 600.
#' @param dispbins positive integer. When scale = "lc", defines the resolution of the active dispersal kernel, default dispbins = 10.
#' @param n.clusters positive integer. Defines the number of parallel processes.
#' @param cluster.type character. Defines the type of cluster, default "PSOCK".
#' @param compressed.output logical. Default TRUE, if FALSE provide abundance for each model's subcompartiment; if FALSE abundances are summed per compartment.
#' @param suffix character. Model output suffix for output RDS. 
#' @param verbose integer. if 1 then an overview of population dynamics is printed in the console, if 2 more information on population dynamics are printed out. Default is 0 (silent).
#' @param seeding logical, default \code{FALSE}, if \code{seeding=TRUE} a fixed seed is applied for result reproducibility.  
#' @return Matrix or a list of matrices containing, for each iteration, the number of individuals in each life stage per day (and for each grid cell of the study area if scale="lc" or "rg"). If the argument compressed.output=FALSE (default TRUE), the model returns the daily number of individuals in each life stage sub-compartment.	
#' @example inst/examples/dynamAedes.m.R
#' @seealso Beta regression functions were taken from the R package \code{aomisc}, which may be available at \url{https://github.com/OnofriAndreaPG/aomisc}.
#' @author Matteo Marcantonio \email{marcantoniomatteo@gmail.com}, Daniele Da Re \email{daniele.dare@uclouvain.be}
#' @export

dynamAedes.m <- function(species=NULL, intro.eggs=0, intro.deggs=0, intro.adults=0, intro.juveniles=0, scale=NULL, intro.cells=NULL, jhwv=NULL, temps.matrix=NULL, startd=1, endd=NA, cells.coords=NULL, coords.proj4=NA, lat=NA, long=NA, road.dist.matrix=NULL, avgpdisp=NA, pDispersal=TRUE, iter=1, n.clusters=1, cluster.type="PSOCK", compressed.output=TRUE,	suffix=NA, cellsize=250, maxadisp=600, dispbins=10, verbose=0, seeding=FALSE) {

#%%%%%%%%%%%%%%%%%%%#
### Initial checks
# Check species names (abbreviation allowed: at least two characters)
if(nchar(species) < 2) {
	stop("Please provide an abbreviation that is at least two characters long.")
}
full_species <- c("aegypti", "albopictus", "koreicus", "japonicus")
matches <- grep(paste0("^", species), full_species)
if(length(matches) == 1) {
	species <- full_species[matches]
	} else if(length(matches) > 1) {
		stop("Ambiguous abbreviation provided, please be more specific, exiting...")
		} else {
			stop("Mosquito species not supported, exiting...")
		}

#Safer version of sample
.resample <- function(x, ...) x[sample.int(length(x), ...)]

# Check dayspan as well as date format function
check_date_format <- function(date) {
	if (!grepl("^\\d{4}-\\d{2}-\\d{2}$", date)) {
		stop("Dates in the wrong format: change them to '%Y-%m-%d'.")
	}
}

# Check start date format
check_date_format(as.character(startd))

# Compute vector of dates for eggs introductions
myDates <- seq.Date(from = as.Date(startd), to = as.Date(endd), by="day")

# Determine dayspan
if (is.na(endd)) {
	dayspan <- ncol(temps.matrix) - 1
	} else {
  # Check end date format
  check_date_format(as.character(endd))
  dayspan <- as.integer(as.Date(endd) - as.Date(startd))
}

# Set dayspan length
cells.coords.photo <- cells.coords
if( dayspan>ncol(temps.matrix) ) {
	stop("You're trying to run the model for more days than columns in 'temps.matrix', exiting..." ) 
}
## Set maxlat for checks
if(scale!="ws") {
	maxlat <- abs(max(cells.coords[,2]))
	} else {
		maxlat <- lat
	}
	if( species!="aegypti" & scale=="rg" & maxlat>90 ) {
		if( is.na(coords.proj4) ) {
			stop("No proj4 string for input coordinates. Please set 'coords.proj4' option.")
			} else {
				cells.coords.photo <- terra::crds(terra::project(terra::vect(x=cells.coords, type="points", atts=NULL, crs=coords.proj4), "EPSG:4326"), df=TRUE)
			}
		}
		if( maxlat>90 )	{
			projected=TRUE } else {
				projected=FALSE
			}

### Preamble: declare variables and prepare the parallel environment for the life cycle ###
legind <- 0
## Define globally the average distance of a trip by car (km)
if( is.na(avgpdisp) ) {
	car.avg.trip <- 22.04 * 1000
	} else if( avgpdisp=="ita" ) {
		car.avg.trip <- 18.99 * 1000
		} else if( avgpdisp=="deu" ) {
			car.avg.trip <- 23.31 * 1000
			} else if( avgpdisp=="esp" ) {
				car.avg.trip <- 28.97 * 1000
				} else if( avgpdisp=="fra" ) {
					car.avg.trip <- 19.29 * 1000
					} else if( avgpdisp=="pol" ) {
						car.avg.trip <- 22.65 * 1000
						} else if( avgpdisp=="uk" ) {
							car.avg.trip <- 19.03 * 1000
							} else if ( is.numeric(avgpdisp) ) {
								car.avg.trip <- avgpdisp
								} else (stop("avgpdisp not supported yet..."))
## Derive daylength for laying of diapausing eggs in albopictus/koreicus/japonicus
if( species!="aegypti" ){
	doy <- as.numeric(format(seq(as.POSIXct(startd), as.POSIXct(as.Date(startd)+dayspan), by='day'), "%j"))
	if( scale=="rg" ) {
		photo.matrix <- lapply(doy, function(x){geosphere::daylength(lat=cells.coords.photo[,2], doy=x)})
		photo.matrix <- do.call(cbind,photo.matrix)
		} else if( !is.na(lat)&!is.na(long) ) {
			dl <- daylength(lat,doy)
			photo.matrix <- matrix(dl, nrow=1)
			} else (stop("Something's wrong with scale or lat and long"))
			} else {
				dl <- rep(24,(dayspan))
			}
if (species == "aegypti") {
  photo.matrix <- matrix(0, nrow = nrow(temps.matrix), ncol = dayspan)
  dl <- rep(24, dayspan)
}
## Set dispersal according to scale
dispersal <- if(scale=="lc"){TRUE}else if(scale=="rg"|scale=="ws"){FALSE}else{stop("Wrong scale. Exiting...")}
## Define the `margin` for apply
mrg <- if(scale=="ws"){2}else if(scale=="lc"|scale=="rg"){1}
if( !dispersal ) message("\n ### Model without dispersal ### \n") 
## Define the type of cluster computing environment 
if(cluster.type=="PSOCK") {
	cl <- parallel::makeCluster(spec=n.clusters, type=cluster.type, nnodes=n.clusters, outfile="")
	} else(message("The only supported `cluster.type` is SOCK"))
## Register the environment 
doParallel::registerDoParallel(cl, cores=n.clusters)
foreach::registerDoSEQ()
options(warn = 2)  # Turn warnings into errors
if(seeding) parallel::clusterEvalQ(cl, set.seed(2021))
## Define space dimensionality into which simulations occour
space <- nrow(temps.matrix)
## Set a progress bar
message("##########################################\n## Life cycle iterations have begun... ##\n##########################################")
pb <- txtProgressBar(char = "%", min = 0, max = iter, style = 3)
### End of preamble ###
#%%%%%%%%%%%%%%%%%%%%%#
### Iterations: start parallelised introduction "iteration" ###
rs <- foreach::foreach( iteration=1:iter, .packages="foreach", .export = ls(globalenv()) ) %dopar% {
## Condition to satisfy to stop the life cycle: sum(pop) == 0, in case the day before extinction has happened
stopit <- FALSE
## Update progress bar
legind <- legind+n.clusters
## Vector of propagules to initiate the life cycle
# If intro.cells is a vector of cells then sample a value for each iteration
if( scale=="lc" ) {
	if( nrow(temps.matrix)<2 | !exists("road.dist.matrix") | !exists("cells.coords") ) {
		stop("If scale='lc' then temps.matrix|road.dist.matrix|cells.coords) must exist and nrows must be > 1")
	}
	if( length(intro.cells)>=1 ) {
		intro.cell <- intro.cells[sample(1:length(intro.cells),1)]
		} else {intro.cell <- NA}
# if intro.cell is not NA than use intro.cell, else sample at random a cell along roads (column of road.dist.matrix)
# Eggs
if( intro.eggs!=0 ) {
	e.intro.n <- rep(0,space)
	if( !is.na(intro.cell) ) {
		# e.intro.n[intro.cell] <- intro.eggs #DDR
	  e.intro.n[intro.cell] <- c(intro.eggs, rep(0, length = length(myDates)-1))
		} else {
			# e.intro.n[sample(as.integer(colnames(road.dist.matrix)),1)] <- intro.eggs #DDR
		  e.intro.n[sample(as.integer(colnames(road.dist.matrix)),1)] <- c(intro.eggs, rep(0, length = length(myDates)-1))
		}
		} else e.intro.n <- intro.eggs
		if( intro.deggs!=0 ) {
# Diapause eggs
d.intro.n <- rep(0,space)
if( !is.na(intro.cell) ) {
	d.intro.n[intro.cell] <- intro.deggs
	} else {
		d.intro.n[sample(as.integer(colnames(road.dist.matrix)),1)] <- intro.deggs
	}
	} else d.intro.n <- intro.deggs
# Immatures
if( intro.juveniles!=0 ) {
	i.intro.n <- rep(0,space)
	if( !is.na(intro.cell) ) {
		i.intro.n[intro.cell] <- intro.juveniles
		} else {
			i.intro.n[sample(as.integer(colnames(road.dist.matrix)),1)] <- intro.juveniles
		} 
		} else i.intro.n <- intro.juveniles
# Adults
if( intro.adults!=0 ) {
	a.intro.n <- rep(0,space)
	if( !is.na(intro.cell) ) {
		a.intro.n[intro.cell] <- intro.adults
		} else {
			a.intro.n[sample(as.integer(colnames(road.dist.matrix)),1)] <- intro.adults
		}
		} else a.intro.n <- intro.adults
		} else if( scale=="ws" ) {
			if( nrow(temps.matrix)>1 ) {
				stop( "if scale='lc' then nrow(temps.matrix) must be 1" )
				} else {
				  if(is.list(intro.eggs)){ # multiple introductions
				    e.intro.n <- rep(0, length = length(myDates))
				    count.mIntro <- which(myDates %in% as.Date(names(intro.eggs[[iteration]])))
				    e.intro.n[count.mIntro]<- as.numeric(intro.eggs[[iteration]])}
				  else{# single introduction on the first day
				    e.intro.n <- c(intro.eggs, rep(0, length = length(myDates)-1))
				  }
					d.intro.n <- intro.deggs; i.intro.n <- intro.juveniles; a.intro.n <- intro.adults; road.dist.matrix <- as.data.frame(c(0,0)); names(road.dist.matrix) <- 1
				}
				} else if( scale=="rg" ) {
					if( nrow(temps.matrix)<1 ) {
						stop( "if scale='rg' then nrow(temps.matrix) must be > 1" )
						} else {
						  if(is.list(intro.eggs)){ # multiple introductions
						    e.intro.n <- rep(0, length = length(myDates))
						    count.mIntro <- which(myDates %in% as.Date(names(intro.eggs[[iteration]])))
						    e.intro.n[count.mIntro]<- as.numeric(intro.eggs[[iteration]])
						  }else{ # single introduction on the first day
						    e.intro.n <- c(intro.eggs, rep(0, length = length(myDates)-1))
						  }
              d.intro.n <- intro.deggs; i.intro.n <- intro.juveniles; a.intro.n <- intro.adults; road.dist.matrix <- as.data.frame(c(0,0)); names(road.dist.matrix) <- 1
						}
						} else stop("Wrong scale.")
### Day cycle: Start sequential "day" life cycle into the "iteration" loop ###
## Define counter which serves as an introduction benchmark
if( exists("counter") ) {
	rm(counter)
}
setTxtProgressBar(pb, legind)
foreach(day = 2:dayspan, .combine=c, .export = ls(globalenv())) %do% {
	while( !stopit ) {
		if( !exists("counter") ) {
# Index for dynamic array eggs
de <- if(species=="koreicus"|species=="japonicus") 2 else 1
# Index for dynamic array juveniles
dj <- if(species=="koreicus"|species=="japonicus") 2 else 1
# Index for dynamic array adults
da <- if(species=="koreicus"|species=="japonicus") 2.25 else 1

# Define objects required to store data during a day
counter <- 0; i.temp.v <- 0; d.temp.v <- 0; e.temp.v <- 0; a.egg.n <- 0; a.new.n <- 0; a.degg.n <- 0 
p.life.a <- array(0,c(4,nrow(temps.matrix),6*de), dimnames = list(c("egg", "juvenile", "adult", "diapause_egg"), NULL, paste0("sc",1:(6*de)))); #DDR
# # Determine number of subcompartments per life stage
# n_egg_sub      <- 6 * de         # corrected from 4 to 6
# n_juv_sub      <- 6 * dj
# n_adult_sub    <- 5              # fixed
# n_diapause_sub <- 6 * de         # matching egg structure
# 
# # Max dimension across all life stages
# n_total_sub <- max(n_egg_sub, n_juv_sub, n_adult_sub, n_diapause_sub)
# 
# # Allocate p.life.a
# p.life.a <- array(
#   0L,
#   dim = c(4, nrow(temps.matrix), n_total_sub),
#   dimnames = list(
#     c("egg", "juvenile", "adult", "diapause_egg"),
#     NULL,
#     paste0("sc", seq_len(n_total_sub))
#   )
# )
storage.mode(p.life.a) <- "integer"; outl <- list()
} else counter <- append(counter,day)
### Header:
## Gonotrophic cycle
## Derive daily rate for gonotrophic cycle, i.e. blood meal to oviposition, then transform rate in daily probabiltiy to terminate the gonotrophic cycle.
a.gono.p <- .a.gono_rate.f(temps.matrix[,day]/1000, species)
## Oviposition rate
## Derive oviposition rate, i.e. number of eggs laid per female per day.
a.batc.n <- .a.ovi_rate.f(temps.matrix[,day]/1000, species)
## Adult survival
## Derive daily adult female survival rate
a.surv.p <- .a.surv_rate.f(temps.matrix[,day]/1000, species)
## Immature emergence
## Derive daily immature emergence rate
i.emer.p <- .i.emer_rate.f(temps.matrix[,day]/1000, species)
## Immature survival
## Derive daily immature survival rate
i.mort_rate.v <- -log(.i.surv_rate.f(temps.matrix[,day]/1000, species))
# Set allocation of diapause/non-diapause eggs
if(length(counter) != 1) {
    if(scale == "rg" && species != "aegypti") {
        if(any(photo.matrix[, day - 1] > photo.matrix[, day])) {
            e.diap.p <- .e.dia_rate.f(photo.matrix[, day], species)
        } else {
            e.diap.p <- rep(0, ncol(photo.matrix))
        }
    } else {
        e.diap.p <- if(dl[day - 1] > dl[day]) {
            .e.dia_rate.f(dl[day], species)
        } else {
            0
        } 
    }
} else {
    e.diap.p <- ifelse(scale == "rg", rep(0, ncol(photo.matrix)), 0)
}

## Derive daily egg hatching rate
e.hatc.p <- .e.hatch_rate.f(temps.matrix[,day]/1000, species)
## Derive daily egg survival rate
e.surv.p <- .e.surv_rate.f(temps.matrix[,day]/1000, species)
d.surv.p <- if( species!="aegypti" ) {.d.surv_rate.f(temps.matrix[,day]/1000, species)} else {0}
# Binned (10m) (Log-normal) probability density for active dispersal up to 600m 
if( dispersal ) {f.adis.p <- .a.a_disp.f(sp=species, max.a.disp=maxadisp, disp.bins=dispbins)}
## Gamma probability density of long passive dispersal (from DOI: 10.2790/7028); from 0 to maximum distance of road segments with 1000 m resolution.
if( dispersal ) {f.pdis.p <- dgamma(seq(1,max(road.dist.matrix,na.rm=T),1000),shape=car.avg.trip/(10000/car.avg.trip), scale=10000/car.avg.trip)}
### Events in the (`E`) egg compartment
## `E` has eight sub-compartment: 1:7 for eggs 1-7 days old that can only die or survive, 8 for eggs older than 7 days that can die/survive/hatch
## Binomial draw to find numbers of eggs that die or survive
p.life.a[1,,2:(4*de)] <- apply(t(p.life.a[1,,1:(4*de-1)]),MARGIN=mrg,function(x) rbinom(n=space, size=x, prob=e.surv.p))
if(species!="aegypti") {p.life.a[4,,2:(4*de)] <- apply(t(p.life.a[4,,1:(4*de-1)]),MARGIN=mrg,function(x) rbinom(size=x,n=space,prob=d.surv.p))} else {p.life.a[4,,2:(4*de)] <- 0}
## Introduce eggs if day==1; introduction happens in E sub-compartment 8 as it can be assumed that eggs are most likely to be introduced in an advanced stage of development 
p.life.a[1,,(4*de)] <- if( length(counter)==1 ) {
	e.intro.n[1]
	} else {
	  p.life.a[1,,(4*de)] + e.intro.n[day]
	}
# Diapause eggs
p.life.a[4,,(4*de)] <- if( length(counter)==1 ) {
	d.intro.n
	} else p.life.a[4,,(4*de)]
# Add eggs laid by females the day before (t-1) stored in a.egg.n (end of the day)
p.life.a[1,,1] <- a.egg.n
if(species!="aegypti") {p.life.a[4,,c(1*de)] <- a.degg.n}
# Add eggs that did not hatch yesterday to egg that today are ready to hatch
p.life.a[1,,c(4*de)] <- p.life.a[1,,c(4*de)] + e.temp.v
if(species!="aegypti") {p.life.a[4,,c(4*de)] <- p.life.a[4,,c(4*de)] + d.temp.v}
# Binomial draw to find numbers of eggs 8-d+ old that hatch today (per cell)
e.hatc.n <- rbinom(length(1:space), p.life.a[1,,c(4*de)], prob=e.hatc.p)
if( species=="albopictus" ) {
	if( any(photo.matrix[,day]>photo.matrix[,day-1]) & any(photo.matrix[,day]>11.44) ) {
		d.hatc.n <- rep(0,space)
		ddays <- which(photo.matrix[,day]>11.44)
		d.hatc.n[ddays] <- rbinom(length(ddays), p.life.a[4,ddays,c(4*de)], prob=e.hatc.p)
		} else {d.hatc.n <- 0}
		}else if( species=="koreicus"|species=="japonicus" ) {
			if( photo.matrix[,day]>photo.matrix[,day-1] & any(photo.matrix[,day]>10.71) ) {
				d.hatc.n <- rep(0,space)
				ddays <- which(photo.matrix[,day]>10.71)
				d.hatc.n[ddays] <- rbinom(length(ddays), p.life.a[4,ddays,c(4*de)], prob=e.hatc.p)
				} else {d.hatc.n <- 0}
				} else {d.hatc.n <- 0}
# Remove hatched eggs from eggs 8d+ old
e.temp.v <- p.life.a[1,,(4*de)] - e.hatc.n
if( species!="aegypti" ) {d.temp.v <- p.life.a[4,,c(4*de)] - if(species!="aegypti") {d.hatc.n} else {0}}
# Apply fixed mortality to non hatched 8d+ old eggs
e.temp.v <- rbinom(length(1:space), e.temp.v, prob=0.99)
d.temp.v <- rbinom(length(1:space), d.temp.v, prob=0.99)
### Events in the (`I`) immature compartment
## `I` has 6 sub-compartments representing days from hatching; an immature can survive/die for the first 5 days after hatching, from the 5th day on, it can survive/die and `emerge`.
## Derive mortality rate due to density and add to mortality rate due to temperature sum and derive probability of survival in each cell.
imm.v <- if( scale=="ws" ) {
	sum(p.life.a[2,,2:(6*dj)])
	} else rowSums(p.life.a[2,,2:(6*dj)])
## Derive density-dependent mortality,*2 is to report densities at 1L (original model is for a 2L water habitat.) / jhwv transform density to new liter/cell habitat volume
i.ddmort_rate.v <- exp(.i.ddmort_rate.f(list(i.dens.v=(imm.v*2)/jhwv)))
i.surv.p <- 1-(1-exp(-(i.mort_rate.v + i.ddmort_rate.v)))
## Binomial draw to find numbers of immature that die or survive-and-move to the next compartment
p.life.a[2,,2:(6*dj)] <- apply(t(p.life.a[2,,1:(6*dj-1)]), MARGIN=mrg, FUN=function(x) rbinom(size=x, n=space, prob=i.surv.p))
## Introduce `I` if day==1; introduction happens in `I` sub-compartment 6 
p.life.a[2,,(6*dj)] <- if( length(counter)==1 ) {
	i.intro.n
	} else p.life.a[2,,(6*dj)]
## Add immatures hatched the same day
p.life.a[2,,1] <- e.hatc.n + if(species!="aegypti") {d.hatc.n} else {0}
## Add immatures that did not emerge yesterday to immatures that today are ready to emerge
p.life.a[2,,(6*dj)] <- p.life.a[2,,(6*dj)] + i.temp.v
## Find numbers of immature 5d+ old that emerge before applying mortality (applied as newly emerged adults today)
i.emer.n <- rbinom(length(1:space), p.life.a[2,,(6*dj)], prob=i.emer.p)
## Remove emerged immatures from immatures 5d+ old
i.temp.v <- p.life.a[2,,(6*dj)] - i.emer.n
## Apply mortality to non emerged 5d+ old immatures
i.temp.v <- sapply(1:space, function(x){rbinom(1,i.temp.v[x],prob=i.surv.p[x])})
### Events in the (`A`) adult compartment
## `A` has 5 sub-compartments representing: adults in day 1 and 2 of oviposition [2:3]; 2d+ old adults host-seeking and non ovipositing [4]; 2d+ old blod-fed adults which are not yet laying [1]; 1d old adults, non-laying and non-dispersing [5].
## Introduce blood-fed females if day is 1
p.life.a[3,,1] <- if( length(counter)==1 ) {
	a.intro.n
	} else p.life.a[3,,1]
## Binomial random draw to find newly emerged females (removing males adult from newly emerged adults)
a.new.n <- rbinom(space,i.emer.n, prob=0.5)
## Add blood-fed females which today matured eggs to females with matured eggs from yesterday 
n.ovir.a <- rbinom(space, p.life.a[3,,1], prob=a.gono.p)
p.life.a[3,,2] <- n.ovir.a
## Remove females which today matured eggs from the host-fed compartment
p.life.a[3,,1] <- p.life.a[3,,1] - n.ovir.a
## Find number of eggs laid today by ovipositing females
if( species!="aegypti" & any(e.diap.p!=0) & sum(p.life.a[3,,])>0 ) {
	if(verbose==2) print("Laying diapausing eggs")
## Total number of eggs laid per cell
a.tegg.n <- sapply(1:space, function(x) sum(rpois(sum(p.life.a[3,x,2:3]), a.batc.n[x])))
## Proportion of diapausing and normal eggs
if(scale=="rg") {
	a.degg.n <- sapply(1:space, function(x){rbinom(1,a.tegg.n[x],prob=(e.diap.p[x]))})
	} else {
	a.degg.n <- sapply(1:space, function(x){rbinom(1,a.tegg.n[x],prob=(e.diap.p))})
	}
	a.egg.n <- a.tegg.n-a.degg.n
	} else {
		a.egg.n <- sapply(1:space, function(x) sum(rpois(sum(p.life.a[3,x,2:3]), a.batc.n[x])))
		a.degg.n <- 0
	}
## Find number of adult females surviving today
p.life.a[3,,1:5] <- apply(t(p.life.a[3,,1:5]),MARGIN=mrg,FUN=function(x) rbinom(size=x,n=space,prob=a.surv.p))
## Short-distance active dispersal
if( dispersal ) {
# It happens only if there is any host-seeking female [4], which are the only actively dispersing
if( any(which(p.life.a[3,,4]>0)) ) {
# Find cells (origin) which have at least 1 host-seeking dispersing female
f.ocel.v <- which(p.life.a[3,,4]>0)
# Find how many host-seeking females move and at what distance (e.g., 0-600 m - every 10 m)
f.adis.v <- sapply(f.ocel.v, function(x) rmultinom(1,p.life.a[3,x,4],f.adis.p))
# Find a landing cell for each distance at which each set of adult females from the same origin disperses to.
for (i in 1:length(f.ocel.v)) {
# Get cell of origin from cell index 
e <- f.ocel.v[i]
# Build a distance matrix between cell of origin and all other cells in the grid; round with option -1 returns the dispersal distance to the closest ten to match the resolution of the dispersal kernel. 
cell.dist.matrix <- terra::distance(matrix(cells.coords[e,],ncol=2),cells.coords,lonlat=!projected)
subcells <- which(cell.dist.matrix<=maxadisp)
cell.dist.matrix <- round(as.matrix(cell.dist.matrix[subcells],ncol=1,as.table=T), -1)
row.names(cell.dist.matrix) <- subcells
# Set matrix columns as the number of the row in the grid
# names(cell.dist.matrix) <- 1:nrow(cells.coords)
# Find set of dispersing cells for each `origin` based on the selected dispersing distances; any mosquito that disperses <250 m stays in the cell of origin. Note: *10-10 is to transform the cell id (1-60) to 100-600 to match the range of the dispersal kernel and set the first dispersal distance to 0 (not 1).
a.plan.l <- sapply(which(f.adis.v[,i]<max(cell.dist.matrix))*10-10, function(x) {
	if( x<cellsize ) {
		x <- as.numeric()
		} else {
			x<-cell.dist.matrix[as.numeric(which(cell.dist.matrix==x)),1]
			}; return(x)
			})
# Set dispersing distances which do not exist in the distance matrix (cell.dist.matrix) as NULL; for example distances > 600, maximum of the dispersal kernel
a.plan.l <- lapply(a.plan.l, function(x) if( (length(x)>0) & (is.null(names(x)) | any(is.na(names(x)))) ) {
	x=numeric()
	} else x)
# Randomly choose a landing cell for each dispersing distance, thus the direction of dispersal is completely random within each distance class (non directional or anisotrophic dispersal). Set cells not selected for landing as -999
a.land.v <- sapply(a.plan.l, function(x) {if( length(x)>0 ) {
	.resample(x,1,replace=FALSE)
	} else -999})
toret <- as.integer(which(a.land.v>=0))
# If there is any dispersal distance farther than 0 (alias >250 m, see above) for a cell of origin, then proceed with moving dispersing females from origin to landing cells
if( length(f.adis.v[which(f.adis.v[,i]<max(cell.dist.matrix)),i][toret])>0 ) {
# Remove dispersing individuals from `origin`
p.life.a[3,e,4] <- p.life.a[3,e,4] - sum(f.adis.v[which(f.adis.v[,i]<max(cell.dist.matrix)),i][toret])
# Add dispersing individuals in the chosen landing cell for each distance
p.life.a[3,as.integer(names(a.land.v[toret])),4] <- p.life.a[3,as.integer(names(a.land.v[toret])),4] + f.adis.v[which(f.adis.v[,i]<max(cell.dist.matrix)),i][toret]
} else if(verbose==2) print("Actively dispersing females stay in the cell of origin. Jumping to the next cell of origin...")
}
} else {if(verbose==2) print("No active dispersing females today...")}
## Medium-distance passive dispersal
# It happens only if a cell with at least 1 female touches a road segement; thus the order of colnames(road.dist.matrix)  must be the same as the row names in `p.life.a`
if( any(which(p.life.a[3,,]>0)%in%colnames(road.dist.matrix)) & pDispersal ) {
# Extract cells (`origin`) which contain medium-distance dispersing females
f.opac.n <- unique(which(p.life.a[3,,]>0,arr.ind=T)[,1])[which(unique(which(p.life.a[3,,]>0,arr.ind=T)[,1])%in%colnames(road.dist.matrix))]
# Binomial draw to select fraction of females which are moved by a car (from DOI: 10.1038/s41598-017-12652-5) in each of the selected `origin`
f.pdis.n <- lapply(f.opac.n, function(x) sapply(p.life.a[3,x,], function(y) rbinom(1,y,prob=0.0051)))
# Find distance along roads at which females from each `origin` will move (each row*1000 is a distance category)
f.mdis.n <- lapply(f.pdis.n, function(x) sapply(x, function(y) rmultinom(1,y,f.pdis.p)))
# Select only distances>0
landing_d <- lapply(f.mdis.n, function(x) which(rowSums(x)>0))
# Move dispersing females to a landing cells and remove them from `origin` cells
if( length(landing_d)>0 ) {
	for ( op in 1:length(f.opac.n) ) {
		for ( lp in 1:length(landing_d[op]) ) {
# Select landing cells along road corresponing to the dispersing distance and subset them at random to a single landing cell for each `origin`
ll <- as.integer(colnames(road.dist.matrix)[which(road.dist.matrix[,which(colnames(road.dist.matrix)%in%f.opac.n[op])]%in%(landing_d[[op]]*1000))])
if( length(ll)>0 ){
	ll <- sample(ll,1,replace=T)
# Remove long dispersing females from cell of `origin`
p.life.a[3,f.opac.n[op],] <- p.life.a[3,f.opac.n[op],] - f.mdis.n[[op]][landing_d[[op]][lp],]
# Add long dispersing females to cells of `landing`
p.life.a[3,ll,] <- p.life.a[3,ll,] + f.mdis.n[[op]][landing_d[[op]][lp],]
}else if(verbose==2) print("Passively dispersing females stay in the cell of origin. Jumping to the next cell of origin...")
}
}
} else if(verbose==2) print("No medium-dispersing females today...")
}
}
## Print information on population structure today
if( verbose==1 ) {
	message("\n", as.Date(startd)+day, ". Day ",length(counter),"-- of iteration ",iteration," has ended. Population is e: ",sum(p.life.a[1,,])," i: ",sum(p.life.a[2,,])," a: " ,sum(p.life.a[3,,]), " d: ",sum(p.life.a[4,,]), " eh: ", sum(e.hatc.n+if(species!="aegypti") {d.hatc.n} else {0}), " el: ",sum(a.egg.n), " \n", "Max t: ",max(temps.matrix[,day]/1000), " \n", "Max hp: ",round(max(e.hatc.p),3), " \n")
	} else if( verbose==2 ) {
		message("\n", as.Date(startd)+day, ". Day ",length(counter),"-- of iteration ",iteration," has ended. Population is:\ne:",sum(p.life.a[1,,])," i: ",sum(p.life.a[2,,])," a: " ,sum(p.life.a[3,,]), " a-hs: " ,sum(p.life.a[3,,4]), " d: ",sum(p.life.a[4,,]), " eh: ", sum(e.hatc.n+if(species!="aegypti") {d.hatc.n} else {0}), " el: ",sum(a.egg.n), " \n", "Max t: ",max(temps.matrix[,day]/1000), " \n")
	}
## Verify exinction
stopit <- sum(p.life.a,na.rm=TRUE)==0
if ( stopit & verbose>=0 ) message("Iteration ", iteration, " is extinct...\n")
## Some (unnecessary?) garbage cleaning
gc()
## if TRUE arrays are compressed (by summing) in matrices so that information on sub-compartements is lost.
## Pre-end-of-day housekeeping
p.life.a[1,,(4*de)] <- p.life.a[1,,(4*de)] - e.hatc.n # Remove new hatched eggs from matrix
# Compress or not? Anyway, make daily output
p.life.aout <- if( compressed.output ) {
	p.life.a_out <- apply(p.life.a, MARGIN=c(1, 2), sum)
	} else { 
		p.life.a_out <- p.life.a
	}
## End-of-day housekeeping to prepare next day
# Make a new host-seeking compartment and slide blood-fed and ovipositing female status to prepare the life cycle for tomorrow (t+1); remove hatached eggs from embryonated eggs
p.life.a[3,,1] <- p.life.a[3,,4] # Host-seeking to blood-fed
p.life.a[3,,4] <- p.life.a[3,,3] + p.life.a[3,,5] # end gono + new to host seeking
p.life.a[3,,5] <- a.new.n # new vector into matrix
p.life.a[3,,3] <- p.life.a[3,,2] # ovi d1 to ovi d2
p.life.a[3,,2] <- 0 # clean ovi d1 (d1 can only last for one day)
# Return the daily summary
return(list(p.life.aout))
} #end of stopif condition
}
}
if( compressed.output ) {
	attributes(rs) <- list(compressed=TRUE)
	} else attributes(rs) <- list(compressed=FALSE)
### Complete final tasks, then return data and exit:
if( !is.na(suffix) ) {
	message(paste("\n\n\nAll iterations concluded. Saving the output to: ",suffix,".RDS\n\n\n",sep=""))
	saveRDS(rs, paste(suffix,".RDS",sep=""))
	} else message("\n\n\n########################################\n## Iterations concluded.##\n########################################\n\n\n")
# Close cluster
parallel::stopCluster(cl)

#Return results as a dynamAedesClass
return(
	new("dynamAedesClass", 
		species=species, 
		scale=ifelse(scale=="ws","Weather Station", ifelse(scale=="lc","Local","Regional")), 
		start_date=startd, end_date=ifelse(is.na(endd), "1990-01-01", endd), 
		n_iterations=iter, 
		stage_intro=ifelse(any(sapply(intro.eggs, function(x){x!=0})), "egg", ifelse(intro.juveniles!=0, "juvenile", ifelse(intro.adults!=0, "adult", "Diapause egg"))),
		n_intro=ifelse(any(sapply(intro.eggs, function(x){x!=0})), intro.eggs[[1]], ifelse(intro.juveniles!=0, intro.juveniles, ifelse(intro.adults!=0,intro.adults, intro.deggs))),
		# stage_intro=ifelse(intro.eggs!=0, "egg", ifelse(intro.juveniles!=0, "juvnile", ifelse(intro.adults!=0, "adult", "Diapause egg"))), #DDR
		# n_intro=ifelse(intro.eggs!=0, intro.eggs, ifelse(intro.juveniles!=0, intro.juveniles, ifelse(intro.adults!=0,intro.adults, intro.deggs))), #DDR
		coordinates=if(!scale=="ws") {matrix(cells.coords, ncol=2)} else{matrix(c(long, lat), byrow=TRUE, nrow=1)},
		compressed_output=compressed.output,
		jhwv=jhwv, 
		simulation=rs,
		dispersal=if(scale=="lc") {list(
			pDispersal=pDispersal,
			avg_p_disp_distance=avgpdisp,
			max_a_disp_distance=maxadisp,
			bins_a_disp_distance=dispbins,
			cellSize_a_disp_distance=cellsize
			)} else {as.list("Simulations without dispersal")}
		)
	)
}
