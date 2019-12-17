## Demographic model for an organism with three life stages (egg, immature [instars and pupa] and adult).
# Spatially explicit; T-dependent daily survival and development rates.
# Short active and long passive dispersal along road networks
#.a=adult, i.=immature, e.=egg
#.p=probability, .r=rate, .n=number, .v=vector, .m=matrix, .a=array, .f=function

zanzinv <- function(temps.matrix=NULL,cells.coords=NULL,road.dist.matrix=NULL,startd=1,endd=10,n.clusters=1,cluster.type="SOCK",iter=1,intro.cell=NULL,intro.adults=0,intro.immatures=0,intro.eggs=0,sparse.output=FALSE,compressed.output=FALSE) {

	### Preamble: define variables for the model ###
	## Export variables in the global environment
	source("./other/libraries.r")
	## Define globally a "safer version of "sample" function
	resample <- function(x, ...) x[sample.int(length(x), ...)]
	sapply(c("libraries","resample","cluster.type"), function(x) {assign(x,get(x),envir= .GlobalEnv)})
	## Load packages
	suppressPackageStartupMessages(libraries(c("foreach","doSNOW","Rmpi","actuar","fields","slam","Matrix")))
	## Type of cluster
	if(cluster.type=="SOCK" || cluster.type=="FORK") {
		cl <- makeCluster(n.clusters,type=cluster.type, outfile="",useXDR=FALSE,methods=FALSE,output="")
	} else if(cluster.type=="MPI") {
		cl <- snow::makeMPIcluster(n.clusters,outfile="",useXDR=FALSE,methods=FALSE,output="")
	}
 	## Start the parallelized loop over iter
	doSNOW::registerDoSNOW(cl)
	parallel::clusterExport(cl=cl, varlist=c("libraries", "resample"))
	parallel::clusterCall(cl=cl, function() libraries(c("foreach","slam")))
 	## Space
	space <- nrow(temps.matrix)
	## Time
	days <- ncol(temps.matrix) #n of day to simulate
	## Vector of propagules to initiate the life cycle
	if(intro.eggs!=0) {e.intro.n <- rep(0,space); if(!is.na(intro.cell)) e.intro.n[intro.cell] <- intro.eggs else  e.intro.n[sample(as.integer(colnames(road.dist.matrix)),1)] <- intro.eggs} else e.intro.n <- rep(0,space)
	if(intro.immatures!=0) {i.intro.n <- rep(0,space); if(!is.na(intro.cell)) i.intro.n[intro.cell] <- intro.immatures else i.intro.n[sample(as.integer(colnames(road.dist.matrix)),1)] <- intro.immatures} else i.intro.n <- rep(0,space)
	if(intro.adults!=0) {a.intro.n <- rep(0,space); if(!is.na(intro.cell)) a.intro.n[intro.cell] <- intro.adults else a.intro.n[sample(as.integer(colnames(road.dist.matrix)),1)] <- intro.adults} else a.intro.n <- rep(0,space)
	
	### Parallelized iterations of the life cycle
	rs <- foreach(iteration=1:iter) %dopar% {
		stopit<-FALSE
		### Life cycle ###
		if( exists("counter") ) rm(counter)
			foreach(day = startd:endd, .combine=c) %do% {
				#print("start of the day")
				if(!stopit) {
					if( !exists("counter") ) {
						counter <- 0; i.surv.m <- matrix(0,ncol=5,nrow=2); i.temp.v <- 0; e.temp.v <- 0; a.egg.n <- 0; p.life.a <- array(0,c(3,nrow(temps.matrix),5)); outl <- list()
					}else counter <- append(counter,day)
					### Header: load temperature dependent functions ###
					## Gonotrophic cycle
					source("./lc/a.gono_rate.f.r")
					## Derive daily rate for gonotrophic cycle,i.e. blood meal to oviposition. 
					a.gono.r <- a.gono_rate.f(temps.matrix[,day]/1000)
					## Transform rate in daily probabiltiy to terminate the gonotrophic cycle.
					a.gono.p <- 1-exp(-(a.gono.r))
					## Oviposition rate
					source("./lc/a.ovi_rate.f.r")
					## Derive oviposition rate, i.e., number of eggs laid per female per day 
					a.batc.n <- a.ovi_rate.f(temps.matrix[,day]/1000)
					## Adult survival
					source("./lc/a.surv_rate.f.r")
					## Derive daily adult female survival rate
					a.surv.r <- a.surv_rate.f(temps.matrix[,day]/1000)
					## Transform rate in daily probabiltiy to survive.
					a.surv.p <- exp(-a.surv.r)
					## Add difference between lab and field survival (from Brady et al. 2014)
					a.surv.p <- ifelse(a.surv.p<0,0,a.surv.p)
					## Immature survival
					source("./lc/i.surv_rate.f.r")
					## Derive daily immature survival rate
					i.surv.r <- i.surv_rate.f(temps.matrix[,day]/1000)
					## Transform rate in daily probabiltiy to survive.
					i.surv.p <- exp(-i.surv.r)
					## Immature emergence
					source("./lc/i.emer_rate.f.r")
					## Derive daily immature emergence rate
					i.emer.r <- i.emer_rate.f(temps.matrix[,day]/1000)
					## Transform rate in daily probabiltiy to survive.
					i.emer.p <- exp(-i.emer.r)
					## Probability of short active dispersal (from Marcantonio et al 2019)
					f.adis.p <- dlnorm(seq(0,600,10),meanlog=4.95,sdlog=0.66)
					## Probability of long passive dispersal (from Pasaoglu et al. 2012)
					f.pdis.p <- rgamma(seq(1,max(road.dist.matrix,na.rm=T),1000),shape=16/(10000/16), scale=10000/16)
					## Probability of egg survival
					e.surv.p <- 0.99
					## Probability of hatching (from Soares-Pinheiro et al. 2015)
					e.hatc.p <- exp(-rgamma(length(temps.matrix[,day]/1000),shape=(7.6/7.4)^2,rate=(7.6/7.4^2)))

					## Egg compartment ##
					# E has four sub-compartment: 1:3 for eggs 1-3 days old that can only die or survive, 4 can die/survive/hatch
					# Binomial draw to find numbers of eggs that die or survive
					p.life.a[1,,2:4] <- t(sapply(1:space, function(x) {sapply(p.life.a[1,x,1:3],rbinom,n=1,p=e.surv.p)}))
					# Introduce if day is 1
					p.life.a[1,,4] <- if(length(counter)==1) e.intro.n else p.life.a[1,,4]
					# Add eggs laid by females the day before
					p.life.a[1,,1] <- a.egg.n
					# Add eggs that not hatched yesterday to egg that today are ready to hatch
					p.life.a[1,,4] <- p.life.a[1,,4] + e.temp.v
					# Random binomial draw to find numbers of eggs 4d+ old that don't hatch or hatch
					e.hatc.n <- sapply(1:space, function(x){rbinom(1,p.life.a[1,x,4],e.hatc.p[x])})
					# Remove hatched eggs from eggs 4d+ old
					e.temp.v <- p.life.a[1,,4] - e.hatc.n
					# Apply mortality to non hatched 4d+ old eggs
					e.temp.v <- sapply(1:space, function(x){rbinom(1,e.temp.v[x],e.surv.p)})

					## Immature compartment ##
					# I has five sub-compartments representing days from hatching; an immature can survive/die for the first four days after hatching but from the fifth day on, it can survive/die/emerge. 
					# Binomial draw to find numbers of immature that die or survive passing to the next compartment
					p.life.a[2,,2:5] <- t(sapply(1:space, function(x) {sapply(p.life.a[2,x,1:4],rbinom,n=1,p=i.surv.p[x])}))
					# Introduce if day is 1
					p.life.a[2,,5] <- if( length(counter)==1 ) i.intro.n else p.life.a[2,,5]
					# Add immatures hatched the day before
					p.life.a[2,,1] <- e.hatc.n
					# Add immatures that did not emerge yesterday to immature that today are ready to emerge
					p.life.a[2,,5] <- p.life.a[2,,5] + i.temp.v
					# Random binomial draw to find numbers of immature 5d+ old that not emerge or emerge
					i.emer.n <- sapply(1:space, function(x){rbinom(1,p.life.a[2,x,5], c(i.emer.p[x]))})
					# Remove emerged immatures from imamtures 5d+ old
					i.temp.v <- p.life.a[2,,5] - i.emer.n
					# Apply mortality to non emerged 5d+ old immatures
					i.temp.v <- sapply(1:space, function(x){rbinom(1,i.temp.v[x],i.surv.p[x])})

					## Adult compartment ##
					# A has five sub-compartments representing: adults in 1,2 of oviposition[2:3]; 2d+ old adults host-seeking and non ovipositing[4]; 2d+ old adults which digested the blood meal but not yet laying[1]; 1d old adults, non-laying and non-dispersing [5].
					# Introduce ovipositing females if day is 1
					p.life.a[3,,1] <- if( length(counter)==1 ) a.intro.n else p.life.a[3,,1]
					# Remove males adult from newly emerged adults
					p.life.a[3,,5] <- sapply(1:space,function(x) rbinom(1,i.emer.n[x],0.5))
					# Add to females ready to oviposit the number of females which pass from host-seeking to ovipositing
					n.ovir.a <- sapply(1:space, function(x) rbinom(1, p.life.a[3,x,4], a.gono.p[x]))
					p.life.a[3,,1] <- p.life.a[3,,1] + n.ovir.a
					# Remove females which stop host-seeking from the host-seeking compartment
					p.life.a[3,,4] <- p.life.a[3,,4] - n.ovir.a
						# Find number of eggs laid by ovipositing females
					a.egg.n <- sapply(1:space, function(x) sum(rpois(sum(p.life.a[3,x,2:3])[which(sum(p.life.a[3,x,2:3])>0)], a.batc.n[x])))
					# Find number of adult females surviving
					p.life.a[3,,1:4] <- t(sapply(1:space, function(x) {sapply(p.life.a[3,x,1:4],rbinom,n=1,p=a.surv.p[x])}))
					#print(c(sum(a.intro.n), sum(p.life.a[3,,1]), a.surv.p[28622],day))

					## Short-distance active dispersal
					# It happens only if there is any host-seeking female, which are the only actively dispersing
						if( any(which(p.life.a[3,,4]>0)) ) {
						# Find cells which have at least 1 host-seeking dispersing adult (cells of origin)
						f.ocel.v <- which(p.life.a[3,,4]>0)
						# Find how many host-seeking females move and at what distance (0:600 m)
						f.adis.v <- sapply(f.ocel.v, function(x) rmultinom(1,p.life.a[3,x,4],f.adis.p))
						#print(p.life.a[3,,1:5][unique(which(p.life.a[3,,1:5]>0,arr.ind=T)[,1]),unique(which(p.life.a[3,,1:5]>0,arr.ind=T)[,2])])
						#print("just before active-dispersal")
						# Find a landing cell for each distance at which a set of adult females disperse 
						for (i in 1:length(f.ocel.v)) {
							#print(c("In short-dispersal",length(f.ocel.v)))
							# Cell of origin from cell index
							e <- f.ocel.v[i]
							# Build a distance matrix between cell of origin and all other cells in the system
							cell.dist.matrix <- sapply(fields::rdist(cells.coords[e,],cells.coords,compact=T),function(x) round(x,-1))
							names(cell.dist.matrix) <- 1:nrow(cells.coords)
							# Find set of dispersing cells for each cell of origin based on the selected dispersing distance; any mosquito that disperses less than 250 meters stays in the cell of origin
							a.plan.l <- sapply(which(f.adis.v[,i]<max(cell.dist.matrix))*10-10, function(x) {
								if(x<250) {x<-as.numeric()} else {x<-cell.dist.matrix[as.numeric(which(cell.dist.matrix==x))]}; return(x)})
							# Set dispersing distances which do not exist in the distance matrix cell.dist.matrix as NULL
							a.plan.l <- lapply(a.plan.l, function(x) if( (length(x)>0) & (is.null(names(x))|any(is.na(names(x)))) ) x=numeric() else x)
							# Randomly choose one of the landing cells for each dispersing distance
							a.land.v <- sapply(a.plan.l, function(x) {if(length(x)>0) resample(x,1,replace=FALSE) else -999})
							#print(a.land.v)
							toret <- as.integer(which(a.land.v>=0))
							#print("just before adding active-adults to cells")
							#print(length(f.adis.v[which(f.adis.v[,i]<max(cell.dist.matrix)),i][toret])>0)
							if( length(f.adis.v[which(f.adis.v[,i]<max(cell.dist.matrix)),i][toret])>0 ) {
								#print(c("toreret is ", e, toret, a.land.v[toret]))
								# Remove dispersing individuals from cell of origin
								p.life.a[3,e,4] <- p.life.a[3,e,4] - sum(f.adis.v[which(f.adis.v[,i]<max(cell.dist.matrix)),i][toret])
								#print("then")
								# Add dispersing individuals in the chosen landing cell for each distance
								p.life.a[3,as.integer(names(a.land.v[toret])),4] <- p.life.a[3,as.integer(names(a.land.v[toret])),4] + f.adis.v[which(f.adis.v[,i]<max(cell.dist.matrix)),i][toret]
							#print(p.life.a[3,,1:5][unique(which(p.life.a[3,,1:5]>0,arr.ind=T)[,1]),unique(which(p.life.a[3,,1:5]>0,arr.ind=T)[,2])])
							#if(any(which(p.life.a[3,,1:5]<0))) print(which(p.life.a[3,,1:5]<0,arr.ind=T),cat("\n\n\n")); stop("negative")
							}else print("Jumping to the next cell of origin")
						}
						#print("out from short dispersal") 
					}else print("No short-dispersal")
					#print("just before passive-dispersal")

					## Long-distance passive dispersal
					# It happens only if any cell with at least one adult is in proximity of roads; therefore it is key that colnames(road.dist.matrix) and the order of cells correspond to p.life.a
					if( any(which(p.life.a[3,,]>0)%in%colnames(road.dist.matrix)) ) {
						# Extract cells whose contain long-distance dispersing adults
						f.opac.n <- unique(which(p.life.a[3,,]>0,arr.ind=T)[,1])[which(unique(which(p.life.a[3,,]>0,arr.ind=T)[,1])%in%colnames(road.dist.matrix))]
						# Select only 0.0001 adults in those cells, meaning that, on average, 1 adult on 10000 is moved by a car
						f.pdis.n <- lapply(f.opac.n, function(x) sapply(p.life.a[3,x,], function(y) rbinom(1,y,0.0001)))
						# Disperse adults at ld distance along roads, each row*1000 is a distance category
						f.mdis.n <- lapply(f.pdis.n, function(x) sapply(x, function(y) rmultinom(1,y,f.pdis.p)))
						# Select drawn distances
						landing_d <- lapply(f.mdis.n, function(x) which(rowSums(x)>0))
						# Place dispersing adults in landing cells and remove them from origin cells
						if( length(landing_d)>0 ) {
							for ( op in 1:length(f.opac.n) ) {
								for ( lp in 1:length(landing_d[op]) ) {
								# Select landing cells along road corresponing to the dispersing distance and subset at random to a single landing cell 
									ll <- as.integer(colnames(road.dist.matrix)[which(road.dist.matrix[,which(colnames(road.dist.matrix)%in%f.opac.n[op])]%in%(landing_d[[op]]*1000))])
									if( length(ll)>0 ){
										ll <- sample(ll,1,replace=T)
									# Remove long dispersing mosquitoes from cell of origin
										p.life.a[3,f.opac.n[op],] <- p.life.a[3,f.opac.n[op],] - f.mdis.n[[op]][landing_d[[op]][lp],]
									# Add long dispersing mosquitoes to cells of landing
										p.life.a[3,ll,] <- p.life.a[3,ll,] + f.mdis.n[[op]][landing_d[[op]][lp],]
									#print(p.life.a[3,f.opac.n[op],])
									}else(print("No long-dispersal"))
								}
							}
						}
					}
					# Make a new host-seeking compartment and slide ovipositing female status for new day
					p.life.a[3,,4] <- p.life.a[3,,3] + p.life.a[3,,4] + p.life.a[3,,5]
					p.life.a[3,,2:3] <- p.life.a[3,,1:2]
					p.life.a[3,,1] <- 0
					message("\nday ",length(counter)," has ended. Population is ",(sum(p.life.a)-sum(i.emer.n))," individuals \n")
					# Condition for exinction
					stopit <- sum(p.life.a)==0
					gc()
					# Condition for 2D output
					if(compressed.output) p.life.aout <- apply(p.life.a, MARGIN=c(1, 2), sum)
					if(sparse.output) return(list(as.simple_sparse_array(p.life.a))) else return(list(p.life.aout))
				}else{message("Extinct")} #end of stopif condition
			} #end of days
			#print("out of day")
		} #end of iterations
		print("out of iteration")
			if(cluster.type == "MPI") {
			message("MPI")
			saveRDS(rs,"~/output_zanzinv.RDS")
			mpi.quit() #Close mpi clusters
		}
		stopCluster(cl)
		return(rs)
	} #end fo function