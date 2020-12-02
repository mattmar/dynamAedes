#Install (if not installed) and load required packages 
libraries <- function(pkg,vvv=verbose) {
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if ( length(new.pkg) ) {
        if( vvv ) {cat("\n        Installing required packages...")}
        install.packages(new.pkg, dependencies = TRUE)
        if( vvv ) {cat("\n        Loading required packages...")}
        sapply(pkg, function(x) suppressPackageStartupMessages(sapply(pkg, function(x) suppressPackageStartupMessages(library(x, character.only=TRUE,quietly=TRUE,warn.conflicts=FALSE,verbose=FALSE)))))
    } else {
        if( vvv ) {cat("\n        Loading required packages...")}
        suppressPackageStartupMessages(sapply(pkg, function(x) suppressPackageStartupMessages(library(x, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE,verbose=FALSE))))
    }
}
