#' S4 class representing the output of dynamAedes.m
#'
#' @slot species Character. The simulated species.
#' @slot scale Character. The scale of the simulation.
#' @slot start_date Character. The introduction date.
#' @slot end_date Character. The end date of the simulation.
#' @slot n_iterations Numeric. The number of iterations.
#' @slot stage_intro Character. The introduced mosquito stage.
#' @slot n_intro Numeric. The number of propagules introduced.
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