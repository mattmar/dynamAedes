#Install (if not installed) and load required packages 
libraries <- function(pkg) {
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if ( length(new.pkg) ) {
        cat("\n        Installing required packages...")
        install.packages(new.pkg, dependencies = TRUE)
        cat("\n        Loading required packages...")
        sapply(pkg, function(x) suppressPackageStartupMessages(sapply(pkg, function(x) suppressPackageStartupMessages(library(x, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE,verbose=FALSE)))))
    } else {
        cat("\n        Loading required packages...")
        suppressPackageStartupMessages(sapply(pkg, function(x) suppressPackageStartupMessages(library(x, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE,verbose=FALSE))))
    }
}
