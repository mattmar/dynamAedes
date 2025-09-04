#' S4 class representing the output of dynamAedes.m
#' @slot species Character. The simulated species.
#' @slot scale Character. The scale of the simulation.
#' @slot start_date Character. The introduction date.
#' @slot end_date Character. The end date of the simulation.
#' @slot n_iterations Numeric. The number of iterations.
#' @slot stage_intro Character. The introduced mosquito stage.
#' @slot n_intro numeric. The number of propagules introduced.
#' @slot coordinates Matrix. The coordinates of each cell.
#' @slot compressed_output Logical. If output is at stage or substage level.
#' @slot jhwv Numeric. The volume of water in the system.

setClass("dynamAedesClass", slots=list(
	species="character", 
	scale="character", 
	start_date="character", 
	end_date="character", 
	n_iterations="numeric", 
	stage_intro="character", 
	n_intro="numeric", 
	coordinates="matrix", 
	compressed_output="logical", 
	jhwv = "numeric",
	dispersal = "list",
	simulation = "list")
)

#' Summary method for dynamAedesClass

#' Provides a summary of simulations based on the dynamAedesClass.
#' 
#' @param object An object of class \code{dynamAedesClass}.
#'
#' @return A character vector with the summary details of the simulation.
#'
#' @examples
#' \dontrun{
#' summary(sim)
#' }
#'
#' @aliases summary,dynamAedesClass-summary
#' @export

# create own method to print? a summary of the iterations
setMethod("summary", "dynamAedesClass", function(object) {
  cat("Summary of dynamAedes simulations:\n")
  cat("----------------------------------\n")
  
  # Just to format strings properly
  formatString <- function(label, value) {
    sprintf("%-25s %s", paste0(label,":"), value)
  }
  species <- paste("Aedes", object@species)
  
  cat(formatString("Species", species), "\n")
  cat(formatString("Scale", toupper(object@scale)), "\n")
  if(object@scale=="lc") {
  	cat(formatString("  Avg passive dispersal", object@dispersal$avg_p_disp_distance), "\n")
  	cat(formatString("  Max active dispersal", object@dispersal$max_a_disp_distance), "\n")
  }
  cat(formatString("Start Date", object@start_date), "\n")
  cat(formatString("End Date", object@end_date), "\n")
  cat(formatString("Number of Iterations", object@n_iterations), "\n")
  cat(formatString("Introduced Stage", object@stage_intro), "\n")
  cat(formatString("Number Introduced", object@n_intro), "\n")
  cat(formatString("Is Output Compressed?", ifelse(object@compressed_output, "Yes", "No")), "\n")
  cat(formatString("Water in the System", object@jhwv), "L", "\n")
  cat(formatString("Min days with population", max(object)), "\n")
  cat(formatString("Max days with population", min(object)), "\n")

  cat("\n")    
  # invisible(object)
})

#' Max method for dynamAedesClass
#' Provides the max number of days with at least one
#' propagule in the system (any stage) along iterations.
#' @param x An object of class \code{dynamAedesClass}.
#' @param na.rm logic.
#' 
#' @return An integer.
#'
#' @examples
#' \dontrun{
#' max(sim)
#' }
#' 
#' @aliases max,dynamAedesClass-max
#' @export

setMethod("max", "dynamAedesClass", 
          function(x, na.rm = FALSE) {
            max_days <- max(sapply(x@simulation, length), na.rm = na.rm)
            return(max_days)
          })

#' Min method for dynamAedesClass

#' Provides the min number of days with at least one
#' propagule in the system (any stage) along iterations.
#' @param x An object of class \code{dynamAedesClass}.
#' @param na.rm logic.
#'
#' @return An integer.
#'
#' @examples
#' \dontrun{
#' min(sim)
#' }
#' 
#' @aliases min,dynamAedesClass-min
#' @export

setMethod("min", "dynamAedesClass", 
          function(x, na.rm = FALSE) {
            min_days <- min(sapply(x@simulation, length), na.rm = na.rm)
            return(min_days)
          })