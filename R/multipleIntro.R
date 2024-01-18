#' Multiple introductions of Aedes sp. propagules
#' 
#' @description
#' `multipleIntro` allows testing scenarios of multiple introduction of invasive mosquitoes. 
#' 
#' @details
#' For a given iteration, the function returns the number of eggs introduced and the dates at which they are introduced. 
#' The function allows to introduce mosquito eggs multiple times for each model iteration, determining how many eggs are introduced and when, based on four possible scenarios: 
#' - \code{"ff"}, i.e. "fixed fixed": egg introductions at fixed intervals in time, and fixed number of eggs introduced each time
#' - \code{"fr"}, i.e. "fixed random": egg introductions at fixed intervals in time, but random number of eggs introduced each time
#' - \code{"rf"}, i.e. "random fixed": egg introductions at random intervals in time, but fixed number of eggs introduced each time
#' - \code{"rr"}, i.e. "random random": egg introductions at random intervals in time, and random number of eggs introduced each time
#' 
#' The function allows introducing two different mosquito life history stages: egg and female adult. By introducing adults, the number of eggs produced by a single female is sampled from a normal distribution
#' with \eqn{\mu = 77} and \eqn{\sigma = 15}, being based on the findings presented in Aida et al. (2011)\url{https://doi.org/10.1016/S2221-1691(11)60103-2}
#' The \code{"fr"} and \code{"rr"} scenarios currently support only the species Aedes albopictus (Aida et al.,2011).
#'
#' @param startd character. The date of the start of simulations, also corresponding to the first introduction event of mosquito propagules.
#' @param endd character. The date of the end of simulations.
#' @param n.intro numeric. The number of times propagules will be introduced throughout the simulations.
#' @param intro.eggs numeric. The average number of eggs introduced per single introduction event. Note: specify either \code{intro.eggs} or \code{intro.adults}. If both are specified, \code{intro.eggs} will be used as default.
#' @param intro.adults numeric. The average number of female mosquitoes introduced per single introduction event. Note: specify either \code{intro.eggs} or \code{intro.adults}. If both are specified, \code{intro.eggs} will be used as default.
#' @param sd_eggs numeric. The standard deviation around the average number of eggs introduced by \code{intro.eggs} or produced by \code{intro.adults} per single introduction event.
#' @param scenario character. The number of eggs and frequency of their introduction scenario. The letter word defines the frequency of egg introduction, which is either at "fixed" or "random" intervals, thus "f" of "r". The second word defines the number of eggs introduced, which is either "fixed" or "random",  thus "f" of "r". The latter corresponds to a number randomly sampled from a normal distribution with a mean defined by either \code{intro.eggs} or \code{intro.adults}. Default \code{scenario = "ff"}.
#' @param iter positive integer. The number of model iterations.
#' 
#' @returns A list of vectors containing the number of eggs introduced each time and the date at which said eggs are introduced. The length of the list is equal to the number of iterations.
#' 
#' @examples
#' \dontrun{out <- multipleIntro(startd = "2022-06-10",
#'                        endd = "2022-08-10",
#'                        n.intro = 4,
#'                        intro.eggs = 100,
#'                        sd_eggs = 15,
#'                        scenario = "ff",
#'                        iter = 100)
#' }
#' @author Emmanuelle Kern \email{emmanuelle.kern22@imperial.ac.uk}, Ilaria Dorigatti  \email{i.dorigatti@imperial.ac.uk}
#' @export

multipleIntro <- function(startd = NULL, endd = NULL, n.intro = NULL, intro.eggs = NULL, intro.adults = NULL, sd_eggs = NULL, scenario = "ff", iter = NULL) {
  
  # create matrix with dates as columns and iterations as rows
  dates <- seq(as.Date(startd),as.Date(endd), by=1)   
  n_col <- length(dates)
  matrix_temp <- matrix(NA, nrow = iter, ncol = n_col)     
  colnames(matrix_temp) <- as.character(dates)
  
  # introducing intro.adults adult females or intro.eggs number of eggs
  if(is.null(intro.eggs)) {
    eggs_per_female <- round(stats::rnorm(n = intro.adults, mean = 77, sd = 15)) # from Aida et al (2011). Later on we might add here a switch for different species
    mean_eggs <- sum(c(eggs_per_female))
  } else {
    mean_eggs <- intro.eggs
  }
  
  ### FIXED INTRODUCTIONS OF EGGS OVER TIME
  # scenario 1: fixed distribution (egg number)
  if(scenario == "ff") {
    for(i in 1:iter){
      col_selected <- round(seq(from = 1, to = ncol(matrix_temp), length.out = n.intro))
      matrix_temp[i, col_selected] <- mean_eggs 
    }
  }
  
  # scenario 2: random distribution (egg number)
  if(scenario == "fr") {
    for(i in 1:iter){
      col_selected <- round(seq(from = 1, to = ncol(matrix_temp), length.out = n.intro))
      matrix_temp[i, col_selected] <- abs(round(stats::rnorm(n = n.intro, mean = mean_eggs, sd = sd_eggs),0))
    }
  }
  
  ### RANDOM INTRODUCTIONS OF EGGS OVER TIME
  if(scenario == "rf" | scenario == "rr") {
    # introduce eggs on first day (so that simulations can start)
    matrix_temp[,1] <- abs(round(stats::rnorm(n = iter, mean = mean_eggs, sd = sd_eggs),0))
    num_col <- list(n_col)
    
    # scenario 3: fixed distribution (egg number)
    if(scenario == "rf") {
      for(i in 1:iter){
        col_selected <- sample(2:ncol(matrix_temp), size = n.intro-1, replace = FALSE)
        matrix_temp[i, col_selected] <- mean_eggs
      }
    }
    
    # scenario 4: random distribution (egg number)
    else {
      for(i in 1:iter){
        col_selected <- sample(2:ncol(matrix_temp), size = n.intro-1, replace = FALSE)      
        matrix_temp[i, col_selected] <- abs(round(stats::rnorm(n = n.intro-1, mean = mean_eggs, sd = sd_eggs),0)) 
      }
    }
  }
  
  # simplify output: for each iteration, show only the number of eggs & the date of introduction as a dataframe
  final_matrix <- apply(X = matrix_temp, MARGIN = 1, FUN = function(x){data.frame(t(x[!is.na(x)]))})
  # format output
  final_matrix <- lapply(final_matrix, function(x){ temp <- unlist(x); 
  names(temp) <- as.Date(gsub("X","", names(temp)), "%Y.%m.%d"); 
  return(temp)})
  # return output
  return(final_matrix)
}