# dynamAedes 2.1.1
* Adding rgdal in the description

# dynamAedes 2.1.0
* Adding website with examples for the three scales (pkgdown)
* Coordinates now internally transformed (if not in lat/long) for the derivation of the photoperiod
* Adding option **coords.proj4** to specify coordinate projection (in proj4 format)
* Renaming dynamAedes as dynamAedes

# dynamAedes 2.0.3
* At regional scale now diapause derived pixel by pixel using the photoperiod of input coordinates.
* You can now introduce diapause eggs using dynamAedes argument
 "intro.degg"

# dynamAedes 2.0.2
* Minor bug fix for idci and icci functions
* Change option name for passive dispersal; option to use a numeric vector for average passive dispersal
* Reshape vignette (simplified) for packaging
* Change vignette builder to knitr
* Add an initial test for reproducibility
* Add code example for dynamAedes

# dynamAedes 2.0.1
* All functions reshaped for publication
* Physiological rates updates
* Vignette updated

# dynamAedes 2.0.0
* Enhanced physiological rates for koreicus/japonicus
* Proportion of diapausing eggs dependent on photoperiod
* Diapausing eggs with fixed survival
* Diapausing eggs are laid only if photoperiod is decreasing and hatch only if photoperiod is rising. 
* Increased minimum duration of egg and juvenile stages to match species biology

# dynamAedes 1.0.1
* Functions to process the output are now integrated
* Adding dispersal kernel for albopictus
* README drafted
* Package now compiles with CRAN standards

# dynamAedes 1.0.0
First release which includes functions to:
* Simulate life cycle