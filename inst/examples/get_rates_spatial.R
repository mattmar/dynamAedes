\dontrun{
# Make a toy temperature time series and coordinates
w <- matrix(as.integer(runif(250,0,45)*1000),ncol=10)
cc <- matrix(c(rep(seq(1,5),5),rev(rep(seq(1,5),each=5))),ncol=2)

# Plot adult daily survival rate
plot(get_rates_spatial(coords = cc,
	temps.matrix = w,
	species = "aegypti",
	rate_fun = ".a.surv_rate.f",
	spatial = TRUE,
	rate = TRUE,
	n.clusters = 1
	))
}