# dynamAedes

## Run the model

###*Punctual* (ws) simulations: no spatial dynamics Requirements: vector of temperatures (\*1000)
aeg.ws <- dynamAedes(species="aegypti", scale="ws", temps.matrix=w, startd=135, endd=175, n.clusters=2, iter=10, intro.eggs=100, ihwv=1).
###*Local* (lc) simulations: passive+active dispersal, cells of the grid are interconnected. Requirements: grid of temperature data (\*1000), distance matrix, coordinates.
aeg.lc <- dynamAedes(species="aegypti", scale="lc", temps.matrix=wsp, cells.coords=cc, road.dist.matrix=dist_matrix, startd=135, endd=175, n.clusters=2, iter=10, intro.adults=100, ihwv=1)
###*Regional* (rg) simulations: spatial dynamics but no dispersal, cells in the grid are indipendent units. Requires: grid of temperature data (\*1000).
aeg.rg <- dynamAedes(species="aegypti", scale="rg", temps.matrix=wsp, startd=135, endd=175, n.clusters=2,iter=10, intro.eggs=100, ihwv=1)

## Processing the output

### Get the probability (proportion) of iterations which have a viable population on the specified date(s)
psucc <- psi(aeg.lc, eval_date=1:100, st=0) #Arena-wide matrix
psucc_sp <- psi_sp(aeg.lc, coords=cc, eval_date=100) #Pixel-by-pixel raster (only for spatial models)
### Get population abundace CI for the specified stage and date(s)
abu <- adci(aeg.lc, eval_date=1:100, stage=1) #Arena-wide matrix
abu_sp <- adci_sp(aeg.lc, coords=cc, eval_date=100, stage=1) #Pixel-by-pixel raster
### Get the number of invaded cells CI (any stage) for the specified date(s)
inv_cells <- icci(aeg.lc, eval_date=1:100) #Arena-wide matrix
### Get the dispersal distance CI for the specified date(s)
disp <- dici(aeg.lc, coords=cc, eval_date=1:100, space=FALSE)
### Get invaded cells for the specified date(s) per iteration
disp <- dici(aeg.lc, coords=cc, eval_date=100, space=TRUE) #Stack of binary rasters
### Get physiological rates in each pixel given a temperature matrix (row=cells; column=days) 
rates <- get_rates_spatial(coords=cc,daily_temp=(wsp)[,200:210]/1000, species="aegypti", rate_fun=".a.ovi_rate.f",spatial=TRUE, np=8)

