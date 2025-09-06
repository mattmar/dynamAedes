#' Summaries of mosquito abundance
#'
#' Summaries of mosquito abundance at each life stage and sub-stage for each day. 
#' @param input_sim matrix. \code{dynamAedes} compressed or uncompressed output matrix. 
#' @param stage character. "Eggs", "Juveniles", "Adults", or "DiapauseEggs" or any shorter attempt longer than 3 letters.
#' @param sub_stage character. For uncompressed outputs only, defines the substage of interest. Please see cheat-sheet table. 
#' @param eval_date positive integer. Define the day to evaluate from the first day of introduction. Note that this can be particularly demanding in the case of spatial outputs.
#' @param breaks numeric vector. Quantile breaks, default the first, the second and the third quartile: \code{c(0.25,0.5,0.75)}.
#' @param n.clusters integer. Number of parallel processes.
#' @param type character. The type of output. Set "O" to force overall time-only summary over spatial. Default is "N" (normal).

#' @return Returns a data frame or a raster with the summary of mosquito abundance at each life stage for each day. 
#' @author Matteo Marcantonio \email{marcantoniomatteo@gmail.com}, Daniele Da Re \email{dare.daniele@gmail.com}
#' @export

adci <- function(input_sim = NULL, stage = NULL, sub_stage = NULL, 
                     breaks = c(0.25, 0.5, 0.75), eval_date = NULL, n.clusters = 1, 
                     type = "N") {
  
  # Preliminary Checks
  if (!is.numeric(eval_date)) stop("eval_date not defined, exiting...")
  if (is.null(input_sim)) stop("Provide input simulations")
  
  # Translate compartments, a bit wordy but it should grep widely.  
  get_comp <- function(stage) {
    stage <- tolower(stage)
    if (grepl("^eg", stage)) return(1)
    if (grepl("^juv", stage)) return(2)
    if (grepl("^ad", stage)) return(3)
    if (grepl("^dia", stage)) return(4)
    return(NULL)    
  }
  
  comp <- get_comp(stage)
  stage <- c("Egg", "Juvenile", "Adult", "DiapauseEgg")[comp]
  
  coords <- input_sim@coordinates
  dd <- max(sapply(input_sim@simulation, function(x) length(x))) # Max days in iterations
  
  # Non-spatial Simulation
  if (input_sim@scale == "ws") {
    message("Weather station scale input simulation provided.") 
    
    if (input_sim@compressed_output) {
      message("Compressed simulation provided")
      message(paste("Computing abundance for ", stage, " estimated over the temporal profile"))
      
      all.matrix <- lapply(1:input_sim@n_iterations, function(x) {
        sapply(eval_date, function(y) {
          unlist(input_sim@simulation[[x]][y])[comp]
        })
      })
      all.matrix <- lapply(all.matrix, function(row) {
        length_diff = length(eval_date) - length(row)
        if (length_diff > 0) {
          c(row, rep(NA, length_diff))
        } else {
          row
        }
      })
      all.matrix <- do.call(cbind, all.matrix)
      all.matrix <- data.frame(t(apply(all.matrix, 1, function(x) {
        quantile(x, probs=breaks, na.rm=T)
      })), day=eval_date, stage=stage)
      
    } else {
      message("Uncompressed simulation provided")
      
      # Check if sub_stage is NULL and handle accordingly.
      if (is.null(sub_stage)) {
        sub_stage <- AedeslifeHistoryList$codesheet[!is.na(AedeslifeHistoryList$codesheet[[stage]]), stage]
        message(paste("Computing", stage, "abundance..." ))
        
        compute_values <- function(x) {
          vals <- sapply(eval_date, function(y) {
            sum(do.call(rbind, input_sim@simulation[[x]][y][1])[sub_stage],na.rm=TRUE)
          })
          c(vals, rep(NA, max(0, length(eval_date) - length(vals))))
        }
        all.matrix <- do.call(cbind, parallel::mclapply(1:length(input_sim@simulation), mc.cores=n.clusters, compute_values))
        all.matrix <- data.frame(t(apply(all.matrix, 1, function(x) quantile(x, probs=breaks, na.rm=T))), day=eval_date, stage=stage)
        
      } else {
        tabCoord <- which(apply(AedeslifeHistoryList$speciesheet, 1:2, function(x) grepl(sub_stage, x)), arr.ind = TRUE)[1,]
        message(paste("Computing sub-compartment", sub_stage, names(AedeslifeHistoryList$speciesheet)[tabCoord[2]], "abundance..." ))
        
        all.matrix <- do.call(cbind, parallel::mclapply(1:input_sim@n_iterations, mc.cores=n.clusters, function(y) {
          sapply(eval_date, function(x) sum(input_sim@simulation[[y]][x][[1]][tabCoord[2],,tabCoord[1]]))
        }))
        all.matrix <- data.frame(t(apply(all.matrix, 1, function(x) quantile(x, probs=breaks, na.rm=TRUE))), day=eval_date, stage=stage, sub_stage = sub_stage)
      }
    }
    
  } else {
    # Spatial Simulation
    message(toupper(paste(input_sim@scale)), " scale input simulation provided") 
    
    if (input_sim@compressed_output) {
      message("Compressed simulation provided")
      if (nrow(coords) <= 1 || is.null(coords) || type == "O") {
        message("Computing overall median temporal abundance estimate.")
        # Here it derives the sum of all the individuals of a stage in a date
        all.matrix <- do.call(cbind, parallel::mclapply(1:input_sim@n_iterations, mc.cores=n.clusters, function(x) {
          sapply(eval_date, function(y) sum(input_sim@simulation[[x]][y][[1]][comp,]))
        }))
        # Then it calculates quantiles of the sum across iterations
        all.matrix <- data.frame(t(apply(all.matrix, 1, function(x) {
          quantile(x, probs=breaks, na.rm=TRUE)
        })), day=eval_date, stage=stage)
        
      } else {
        message("Computing spatio-temporal abundance estimate.")
        
        all.matrix <- parallel::mclapply(1:input_sim@n_iterations, mc.cores= n.clusters, function(x) {
          sapply(eval_date, function(y) {
            # result <- do.call(cbind, input_sim@simulation[[x]][y][1])[comp, , drop = FALSE]
            # as.matrix(result)
            temp_result <- do.call(cbind, input_sim@simulation[[x]][y][1])
            
            # If temp_result is NULL, create a matrix with the correct number of columns
            if (is.null(temp_result)) {
              num_cols <- ncol(do.call(cbind, input_sim@simulation[[1]][eval_date[1]][1])) # Get number of columns from a valid entry
              return(matrix(0, nrow = length(comp), ncol = num_cols))  
            }
            result <- temp_result[comp, , drop = FALSE]
            as.matrix(result)
          })
        })
        
        
        all.matrix <- abind::abind(all.matrix, along=3)
        
        all.matrix <- lapply(breaks, function(x) {
          quant_res <- apply(all.matrix, c(1,2), quantile, probs=x, na.rm=T)
          rast_obj <- rast(data.frame(coords, quant_res), type="xyz")
          setNames(rast_obj, paste0("day", eval_date))
        })
        names(all.matrix) <- paste("q",breaks,sep="_")
      }
    } else {
      message("Uncompressed simulation provided")
      
      if (is.null(sub_stage) && type=="O") {
        message(paste("Computing non-spatial", stage, "abundance..."))
        
        all.matrix <- do.call(cbind, parallel::mclapply(1:input_sim@n_iterations, mc.cores=n.clusters, function(y) {
          sapply(eval_date, function(x) sum(input_sim@simulation[[y]][x][[1]][comp,,]))
        }))
        all.matrix <- data.frame(t(apply(all.matrix, 1, function(x) quantile(x, probs=breaks, na.rm=TRUE))),
                                 day=eval_date, stage=stage)
        
      } else if (is.null(sub_stage) && type=="N")
      {
        message("Computing spatio-temporal abundance estimate at stage level.")
        all.matrix <- parallel::mclapply(1:input_sim@n_iterations, mc.cores= n.clusters, function(x) {
          sapply(eval_date, function(y) {
            result <- rowSums(input_sim@simulation[[x]][y][1][[1]][comp, , ])
            as.matrix(result)
          })
        })
        
        all.matrix <- abind::abind(all.matrix, along=3)
        
        all.matrix <- lapply(breaks, function(x) {
          quant_res <- apply(all.matrix, c(1,2), quantile, probs=x, na.rm=T)
          rast_obj <- rast(data.frame(coords, quant_res), type="xyz")
          setNames(rast_obj, paste0("day", eval_date))
        })
        names(all.matrix) <- paste("q",breaks,sep="_")
      } else {
        if (nrow(coords) <= 1 || is.null(coords) || type == "O") {
          tabCoord <- which(sapply(AedeslifeHistoryList$speciesheet, grepl, pattern=sub_stage), arr.ind=TRUE)[1,]
          message(paste("Computing non-spatial sub-compartment", sub_stage, names(AedeslifeHistoryList$speciesheet)[tabCoord[2]], "abundance..."))
          
          all.matrix <- do.call(cbind, parallel::mclapply(1:input_sim@n_iterations, mc.cores=n.clusters, function(y) {
            sapply(eval_date, function(x) sum(input_sim@simulation[[y]][x][[1]][tabCoord[2],,tabCoord[1]]))
          }))
          all.matrix <- data.frame(t(apply(all.matrix, 1, function(x) quantile(x, probs=breaks, na.rm=TRUE))), day=eval_date, stage=stage, sub_stage = sub_stage)
        } 
        else {
          message("Computing spatio-temporal abundance estimate at substage level.")
          tabCoord <- which(sapply(AedeslifeHistoryList$speciesheet, grepl, pattern=sub_stage), arr.ind=TRUE)[1,]
          message(paste("Computing sub-compartment", sub_stage, names(AedeslifeHistoryList$speciesheet)[tabCoord[2]], "abundance..."))
          
          # Compute all matrix based on the specific sub-stage
          all.matrix <- parallel::mclapply(1:input_sim@n_iterations, mc.cores=n.clusters, function(y) {
            sapply(eval_date, function(x) {
              result <- input_sim@simulation[[y]][x][[1]][tabCoord[2],,tabCoord[1]]
              as.matrix(result)
            })
          })
          all.matrix <- abind::abind(all.matrix, along=3)
          
          all.matrix.raster <- lapply(breaks, function(x) {
            quant_res <- apply(all.matrix, c(1,2), quantile, probs=x, na.rm=T)
            rast_obj <- rast(data.frame(coords, quant_res), type="xyz")
            setNames(rast_obj, paste0("day", eval_date))
          })
          names(all.matrix.raster) <- paste(sub_stage, "q", breaks, sep="_")
          all.matrix <- all.matrix.raster
        }
      }
    }
  }
  return(all.matrix)
}