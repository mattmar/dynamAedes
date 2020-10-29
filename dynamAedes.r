##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
## R (v4.0.2+) code developed by: 
## Matteo Marcantonio, Daniele Da Re, Diego Montecino
## ---------------------------------------------------##
## Latest updated and tested: September 2020
## ---------------------------------------------------##
## Scope: Population dynamical model for `Aedes aegypti`
## ---------------------------------------------------##
## Basic features: Process-based and time-discrete,
## 3 life-stages and sub-stages, temperature-dependent 
## rates and active/passive dispersal based on 
## probability kernels
## ---------------------------------------------------##
## Basic nomenclature: .a.=adult, .i.=immature, .e.=egg
## .p.=probability, .r.=rate, .n.=number, .v.=vector 
## .m.=matrix, .a.=array, .f.=function
## ---------------------------------------------------##
## DOI: 
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
dynamAedes <- function(species="aegypti",temps.matrix=NULL,cells.coords=NULL,
	road.dist.matrix=NULL,startd=1,endd=10,n.clusters=1,
	cluster.type="SOCK",iter=1,intro.cells=NULL,intro.adults=0,
	intro.immatures=0,intro.eggs=0,sparse.output=FALSE,
	compressed.output=FALSE,suffix="DynamAedes",country="it") {
	#%%%%%%%%%%%%%%%%%%%#
	### Preamble: declare variables and prepare the parallel environment for the life cycle ###
	## Install all required packages if not already installed
	source("./of/libraries.r")
	## Define globally a "safer version of "sample" function
	resample <- function(x, ...) x[sample.int(length(x), ...)]
	## Define globally the average distance of a trip by car: data taken (from DOI: 10.2790/7028)
	if( country=="it" ) {
		car.avg.trip <- 18.43
	} else if( country=="nl" ) {
		car.avg.trip <- 23.14
	} else if( country=="es" ) {
		car.avg.trip <- 28.14
	} else ( car.avg.trip=23.24 )
	## Export variables to the global environment
	sapply(c("libraries","resample","cluster.type","car.avg.trip","suffix"), function(x) {assign(x,get(x),envir=.GlobalEnv)})
	## Load required packages
	suppressPackageStartupMessages(libraries(c("foreach","doSNOW","actuar","fields","slam","Matrix","epiR")))
	## Define the type of cluster computing environment 
	if(cluster.type=="SOCK" || cluster.type=="FORK") {
		cl <- makeCluster(n.clusters,type=cluster.type, outfile="",useXDR=FALSE,methods=FALSE,output="")
	} else(message("Default cluster.type is SOCK"))
 	## Register the environment and export newly defined variables and packages to the global environment 
	doSNOW::registerDoSNOW(cl)
	parallel::clusterExport(cl=cl, varlist=c("libraries","resample","car.avg.trip","suffix")) # This loads functions in each child R process
	parallel::clusterCall(cl=cl, function() libraries(c("foreach","slam","epiR"))) # This loads packages in each child R process
  	## Define space dimensionality into which simulations occour
	space <- nrow(temps.matrix)
	## Define time dimensionality into which simulations occour
	days <- ncol(temps.matrix) #n of day to simulate, equal to the number of column of the temperature matrix
	### End of preamble ###
	#%%%%%%%%%%%%%%%%%%%%%#
	### Iterations: start parallelised introduction "iteration" ###
	rs <- foreach( iteration=1:iter ) %dopar% {
		## Condition to satisfy to stop the life cycle: sum(pop) == 0, in case the day before extinction has happened
		stopit <- FALSE
		## Vector of propagules to initiate the life cycle
		# If intro.cells is a vector of cells than sample a value for each iteration
		if( length(intro.cells)>1 ) {
			intro.cell <- sample(intro.cells,1)
		} else {intro.cell <- NA}
		# if intro.cell is not NA than use intro.cell, else sample at random a cell along roads (column of road.dist.matrix)
		# Eggs
		if( intro.eggs!=0 ) {
			e.intro.n <- rep(0,space)
			if( !is.na(intro.cell) ) {
				e.intro.n[intro.cell] <- intro.eggs
			} else {
				e.intro.n[sample(as.integer(colnames(road.dist.matrix)),1)] <- intro.eggs
			}
		}
		# Immatures
		if( intro.immatures!=0 ) {
			i.intro.n <- rep(0,space)
			if( !is.na(intro.cell) ) {
				i.intro.n[intro.cell] <- intro.immatures
			} else {
				i.intro.n[sample(as.integer(colnames(road.dist.matrix)),1)] <- intro.immatures
			} 
		} else i.intro.n <- intro.immatures
		# Adults
		if( intro.adults!=0 ) {
			a.intro.n <- rep(0,space)
			if( !is.na(intro.cell) ) {
				a.intro.n[intro.cell] <- intro.adults
			} else {
				a.intro.n[sample(as.integer(colnames(road.dist.matrix)),1)] <- intro.adults
			}
		} else a.intro.n <- intro.adults

		### Day cycle: Start sequential "day" life cycle into the "iteration" loop ###
		## Define counter which serves as an introduction benchmark
		if( exists("counter") ) {
			rm(counter)
		}
		foreach(day = startd:endd, .combine=c) %do% {
			if( !stopit ) {
				if( !exists("counter") ) {
					# Define objects required to store data during a day
					counter <- 0; i.surv.m <- matrix(0,ncol=5,nrow=2); i.temp.v <- 0; e.temp.v <- 0; a.egg.n <- 0; p.life.a <- array(0,c(3,nrow(temps.matrix),6)); outl <- list()
				} else counter <- append(counter,day)

				### Header: load functions for life cycle parameters. (/1000 because T*1000)
				source(paste("./lcf/",species,".r",sep=""))
				## Gonotrophic cycle
				## Derive daily rate for gonotrophic cycle, i.e. blood meal to oviposition, then transform rate in daily probabiltiy to terminate the gonotrophic cycle.
				a.gono.p <- 1-exp(-(a.gono_rate.f(temps.matrix[,day]/1000)))
				## Oviposition rate
				## Derive oviposition rate, i.e. number of eggs laid per female per day.
				a.batc.n <- a.ovi_rate.f(temps.matrix[,day]/1000)
				## Adult survival
				## Derive daily adult female survival rate and transform rate in daily probabiltiy to survive.
				a.surv.p <- 1-(1-exp(-a.mort_rate.f(temps.matrix[,day]/1000)))
				## Immature survival
				## Derive daily immature survival rate, then transform rate in daily probabiltiy to survive.
				i.surv.p <- 1-(1-exp(-i.mort_rate.f(temps.matrix[,day]/1000)))
				## Immature emergence
				## Derive daily immature emergence rate then transform rate in daily probabiltiy to emerge.
				i.emer.p <- 1-exp(-i.emer_rate.f(temps.matrix[,day]/1000))
				## Derive daily egg survival rate then transform rate in daily probabiltiy to survive.
				e.surv.p <- 1-exp(-e.surv_rate.f(temps.matrix[,day]/1000))
				## Derive daily egg hatching rate then transform rate in daily probabiltiy to hatch.
				e.hatc.p <- 1-exp(-e.hatch_rate.f(temps.matrix[,day]/1000))
				## Gamma probability density of long passive dispersal (from DOI: 10.2790/7028); from 0 to maximum distance of road segments with resolution of 1000 m.
				f.pdis.p <- dgamma(seq(1,max(road.dist.matrix,na.rm=T),1000),shape=car.avg.trip/(10000/car.avg.trip), scale=10000/car.avg.trip)
				
				### Events in the (`E`) egg compartment
				## `E` has four sub-compartment: 1:3 for eggs 1-3 days old that can only die or survive, 4 for eggs older than 3 days that can die/survive/hatch
				## Binomial draw to find numbers of eggs that die or survive
				p.life.a[1,,2:4] <- apply(t(p.life.a[1,,1:3]),MARGIN=1,function(x) rbinom(size=x,n=space,p=e.surv.p))
				## Introduce eggs if day==1; introduction happens in E sub-compartment 4 as it can be assumed that eggs are most likely to be introduced in an advanced stage of development 
				p.life.a[1,,4] <- if(length(counter)==1) {
					e.intro.n
				} else p.life.a[1,,4]
				# Add eggs laid by females the day before (t-1) stored in a.egg.n (end of the day)
				p.life.a[1,,1] <- a.egg.n
				# Add eggs that did not hatch yesterday to egg that today are ready to hatch
				p.life.a[1,,4] <- p.life.a[1,,4] + e.temp.v
				# Binomial draw to find numbers of eggs 4-day+ old that hatch today
				e.hatc.n <- rbinom(length(1:space), p.life.a[1,,4], e.hatc.p)
				# Remove hatched eggs from eggs 4d+ old
				e.temp.v <- p.life.a[1,,4] - e.hatc.n
				# Apply mortality to non hatched 4d+ old eggs
				e.temp.v <- rbinom(length(1:space), e.temp.v, e.surv.p)
				
				### Events in the (`I`) immature compartment
				## `I` has 6 sub-compartments representing days from hatching; an immature can survive/die for the first 5 days after hatching, from the 5th day on, it can survive/die and `emerge`.
				## Binomial draw to find numbers of immature that die or survive-and-move to the next compartment
				p.life.a[2,,2:6] <- apply(t(p.life.a[2,,1:5]), MARGIN = 1, FUN= function(x) rbinom(size=x, n=space, p=i.surv.p))
				## Introduce `I` if day==1; introduction happens in `I` sub-compartment 6 
				p.life.a[2,,6] <- if( length(counter)==1 ) {
					i.intro.n
				} else p.life.a[2,,6]
				## Add immatures hatched the same day
				p.life.a[2,,1] <- e.hatc.n
				## Add immatures that did not emerge yesterday to immatures that today are ready to emerge
				p.life.a[2,,6] <- p.life.a[2,,6] + i.temp.v
				## Find numbers of immature 5d+ old that emerge before applying mortality (applied as newly emerged adults today)
				i.emer.n <- rbinom(length(1:space), p.life.a[2,,6], i.emer.p)
				## Remove emerged immatures from immatures 5d+ old
				i.temp.v <- p.life.a[2,,6] - i.emer.n
				## Apply mortality to non emerged 5d+ old immatures
				i.temp.v <- sapply(1:space, function(x){rbinom(1,i.temp.v[x],i.surv.p[x])})
				
				### Events in the (`A`) adult compartment
				## `A` has 5 sub-compartments representing: adults in day 1 and 2 of oviposition [2:3]; 2d+ old adults host-seeking and non ovipositing [4]; 2d+ old blod-fed adults which are not yet laying [1]; 1d old adults, non-laying and non-dispersing [5].
				## Introduce blood-fed females if day is 1
				p.life.a[3,,1] <- if( length(counter)==1 ) {
					a.intro.n
				} else p.life.a[3,,1]
				## Binomial random draw to remove males adult from newly emerged adults
				p.life.a[3,,5] <- rbinom(space,i.emer.n, 0.5)
				## Add blood-fed females which today matured eggs to females with matured eggs from yesteday 
				n.ovir.a <- rbinom(space, p.life.a[3,,1], a.gono.p)
				p.life.a[3,,2] <- n.ovir.a
				## Remove females which today matured eggs from the host-fed compartment
				p.life.a[3,,1] <- p.life.a[3,,1] - n.ovir.a
				## Find number of eggs laid today by ovipositing females
				a.egg.n <- sapply(1:space, function(x) sum(rpois(sum(p.life.a[3,x,2:3]), a.batc.n[x])))
				## Find number of adult females surviving today
				p.life.a[3,,1:5] <- apply(t(p.life.a[3,,1:5]),MARGIN=1,FUN=function(x) rbinom(size=x,n=space,p=a.surv.p))

				## Short-distance active dispersal
				# It happens only if there is any host-seeking female [4], which are the only actively dispersing
				if( any(which(p.life.a[3,,4]>0)) ) {
					# Find cells (origin) which have at least 1 host-seeking dispersing female
					f.ocel.v <- which(p.life.a[3,,4]>0)
					# Find how many host-seeking females move and at what distance (0-600 m - every 10 m)
					f.adis.v <- sapply(f.ocel.v, function(x) rmultinom(1,p.life.a[3,x,4],f.adis.p))
					# Find a landing cell for each distance at which each set of adult females from the same origin disperses to.
					for (i in 1:length(f.ocel.v)) {
						# Get cell of origin from cell index
						e <- f.ocel.v[i]
						# Build a distance matrix between cell of origin and all other cells in the grid; round with option -1 returns the dispersal distance to the closest ten to match the resolution of the dispersal kernel. 
						cell.dist.matrix <- sapply(fields::rdist(cells.coords[e,],cells.coords,compact=T),function(x) round(x,-1))
						# Set matrix columns as the number of the row in the grid
						names(cell.dist.matrix) <- 1:nrow(cells.coords)
						# Find set of dispersing cells for each `origin` based on the selected dispersing distances; any mosquito that disperses <250 m stays in the cell of origin. Note: *10-10 is to transform the cell id (1-60) to 100-600 to match the range of the dispersal kernel and set the first dispersal distance to 0 (not 1).
						a.plan.l <- sapply(which(f.adis.v[,i]<max(cell.dist.matrix))*10-10, function(x) {
							if( x<250 ) {
								x <- as.numeric()
							} else {
								x<-cell.dist.matrix[as.numeric(which(cell.dist.matrix==x))]
							}; return(x)
						})
						# Set dispersing distances which do not exist in the distance matrix (cell.dist.matrix) as NULL; for example distances > 600, maximum of the dispersal kernel
						a.plan.l <- lapply(a.plan.l, function(x) if( (length(x)>0) & (is.null(names(x)) | any(is.na(names(x)))) ) {
							x=numeric()
						} else x)
						# Randomly choose a landing cell for each dispersing distance, thus the direction of dispersal is completely random within each distance class (non directional or anisotrophic dispersal). Set cells not selected for landing as -999
						a.land.v <- sapply(a.plan.l, function(x) {if( length(x)>0 ) {
							resample(x,1,replace=FALSE)
						} else -999})
						toret <- as.integer(which(a.land.v>=0))
						# If there is any dispersal distance farther than 0 (alias >250 m, see above) for a cell of origin, then proceed with moving dispersing females from origin to landing cells
						if( length(f.adis.v[which(f.adis.v[,i]<max(cell.dist.matrix)),i][toret])>0 ) {
							# Remove dispersing individuals from `origin`
							p.life.a[3,e,4] <- p.life.a[3,e,4] - sum(f.adis.v[which(f.adis.v[,i]<max(cell.dist.matrix)),i][toret])
							# Add dispersing individuals in the chosen landing cell for each distance
							p.life.a[3,as.integer(names(a.land.v[toret])),4] <- p.life.a[3,as.integer(names(a.land.v[toret])),4] + f.adis.v[which(f.adis.v[,i]<max(cell.dist.matrix)),i][toret]
						} else print("Actively dispersing females stay in the cell of origin. Jumping to the next cell of origin...")
					}
				}else print("No active dispersing females today...")

				## Medium-distance passive dispersal
				# It happens only if a cell with at least 1 female touches a road segement; thus  the order of colnames(road.dist.matrix)  must be the same of `p.life.a`
				if( any(which(p.life.a[3,,]>0)%in%colnames(road.dist.matrix)) ) {
					# Extract cells (`origin`) which contain medium-distance dispersing females
					f.opac.n <- unique(which(p.life.a[3,,]>0,arr.ind=T)[,1])[which(unique(which(p.life.a[3,,]>0,arr.ind=T)[,1])%in%colnames(road.dist.matrix))]
					# Binomial draw to select fraction of females which are moved by a car (from DOI: 10.1038/s41598-017-12652-5) in each of the selected `origin`
					f.pdis.n <- lapply(f.opac.n, function(x) sapply(p.life.a[3,x,], function(y) rbinom(1,y,0.0051)))
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
								}else(print("Passively dispersing females stay in the cell of origin. Jumping to the next cell of origin..."))
							}
						}
					} else print("No medium-dispersing females today...")
				}
				## End-of-day housekeeping
				# Make a new host-seeking compartment and slide blood-fed and ovipositing female status to prepare the life cycle for tomorrow (t+1)
				p.life.a[3,,1] <- p.life.a[3,,4]
				p.life.a[3,,4] <- p.life.a[3,,3] + p.life.a[3,,5]
				p.life.a[3,,3] <- p.life.a[3,,2]
				p.life.a[3,,2] <- 0
				# Print information on population structure today
				message("\nday ",length(counter),"-- of iteration ",iteration," has ended. Population is e: ",sum(p.life.a[1,,])," i: ", sum(p.life.a[2,,])," a: " ,sum(p.life.a[3,,]), " eh: ", sum(e.hatc.n), " el: ",sum(a.egg.n), " \n")
				# Condition for exinction
				stopit <- sum(p.life.a)==0
				# Some (unnecessary?) garbage cleaning
				gc()
				# if TRUE arrays are compressed (by summing) in matrices so that information on sub-compartements is irreparabily lost.
				if( compressed.output ) {
					p.life.aout <- apply(p.life.a, MARGIN=c(1, 2), sum)
				}
				# If TRUE a sparse array is returned (save memory but complex to process)
				if(sparse.output) return(list(as.simple_sparse_array(p.life.a))) else return(list(p.life.aout))
			}else{
				message("Extinct")
			} #end of stopif condition
		}
	}
	### Complete final tasks, then return data and exit:
	cat(paste("\nIterations concluded. Saving the output to: ",suffix,".RDS\n\n\n",sep=""))
	saveRDS(rs, paste(suffix,".RDS",sep=""))
	# Close MPI clusters
	if(cluster.type == "MPI") {
		message("MPI")
		mpi.quit() 
	}
	# Close cluster
	stopCluster(cl)
	return(rs)
}