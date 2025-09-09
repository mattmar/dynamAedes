#' Get life-history rate based on temperature or photoperiod
#'
#' Computes one of the internal temperature- or photoperiod-dependent life-history rates 
#' for a given mosquito species. Supports numeric vectors and [terra::SpatRaster()] inputs.
#'
#' @param temp Numeric vector or [terra::SpatRaster()] with temperature values (°C).
#' @param photoperiod Numeric vector or [terra::SpatRaster()] with photoperiod values (hours).
#' @param species Character. Species name. One of: `"aegypti"`, `"albopictus"`, `"koreicus"`, `"japonicus"`.
#' @param rate Character. One of:
#'   - `"a.surv"`: Adult daily survival rate  
#'   - `"a.gono"`: Gonotrophic cycle rate (1/duration)  
#'   - `"a.ovi"`: Oviposition rate (eggs per female/day)  
#'   - `"i.surv"`: Immature daily survival rate  
#'   - `"i.emer"`: Immature daily emergence rate  
#'   - `"i.ddmort"`: Immature density-dependent mortality (based on density, not temperature)  
#'   - `"e.hatch"`: Hatching rate of eggs  
#'   - `"e.surv"`: Egg daily survival rate  
#'   - `"d.surv"`: Diapausing egg daily survival rate  
#'   - `"e.diap"`: Diapause allocation (based on photoperiod)
#'   
#' @author Matteo Marcantonio \email{marcantoniomatteo@gmail.com}, Daniele Da Re \email{dare.daniele@gmail.com}
#' @return A numeric vector (if input is numeric) or [terra::SpatRaster()] if raster input.
#'
#' @import terra
#' @export
#'
#' @examples
#' \dontrun{
#' # Example 1: Vector of temperatures
#' temps <- seq(-5, 40, by = 1)
#' x <- get_rate(temp = temps, species = "aegypti", rate = "a.surv")
#'
#' # Example 2: SpatRaster temperature
#' r <- terra::rast(nrows = 10, ncols = 10, nlyrs = 1)
#' values(r) <- runif(ncell(r), 10, 35)
#' r_rate <- get_rate(temp = r, species = "albopictus", rate = "a.gono")
#' plot(r_rate)
#'
#' # Example 3: Spatiotemporal SpatRaster for photoperiod
#' photo_r <- terra::rast(nrows = 5, ncols = 5, nlyrs = 3)
#' values(photo_r) <- rep(c(13, 12, 11), each = ncell(photo_r))
#' time(photo_r) <- as.Date("2020-06-01") + 0:2
#' diap_stack <- get_rate(photoperiod = photo_r, species = "albopictus", rate = "e.diap")
#' plot(diap_stack)
#' }
get_rate <- function(temp = NULL, photoperiod = NULL, species, rate) {
  # Map internal rate functions
  rate_function <- switch(rate,
                          "a.surv" = .a.surv_rate.f,
                          "a.gono" = .a.gono_rate.f,
                          "a.ovi"  = .a.ovi_rate.f,
                          "i.surv" = .i.surv_rate.f,
                          "i.emer" = .i.emer_rate.f,
                          "i.ddmort" = .i.ddmort_rate.f,
                          "e.hatch" = .e.hatch_rate.f,
                          "e.surv" = .e.surv_rate.f,
                          "d.surv" = .d.surv_rate.f,
                          "e.diap" = .e.dia_rate.f,
                          stop("Unsupported rate type.")
  )
  
  if (rate == "e.diap") {
    if (is.null(photoperiod)) stop("For 'e.diap', provide a photoperiod input.")
    
    if (inherits(photoperiod, "SpatRaster")) {
      r <- terra::app(photoperiod, fun = function(x) rate_function(x, species))
      names(r) <- paste0(rate, "_", species, "_", seq_len(nlyr(r)))
      return(r)
    } else if (is.numeric(photoperiod)) {
      return(rate_function(photoperiod, species))
    } else {
      stop("Unsupported photoperiod input format.")
    }
  }
  
  if (is.null(temp)) stop("Temperature input is required for selected rate.")
  
  if (inherits(temp, "SpatRaster")) {
    r <- terra::app(temp, fun = function(x) rate_function(x, species))
    names(r) <- paste0(rate, "_", species, "_", seq_len(nlyr(r)))
    return(r)
  } else if (is.numeric(temp)) {
    return(rate_function(temp, species))
  } else {
    stop("Unsupported temperature input format.")
  }
}