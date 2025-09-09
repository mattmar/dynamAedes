# dynamAedes 3.0.0
* Fixes bug on dybamAedes.m at the regional scale when the species is Aedes aegypti.
* Updates adci function
* Implements Multiple introduction events through the multipleIntro function
* Adds the get_rate function to compute mosquito temperature/photoperiod dependent rates

# dynamAedes 2.2.9
Substitutes legend.pos with legend.position in vignettes.

# dynamAedes 2.2.8
* Adds reference file for life history traits `AedeslifeHistoryList`.
* Adds function `spreader`.
* Adds vignettes 4 (Uncompressed output) and 5 (spreader).
* Fixes minor issues in `adci`.
* check summary print for punctual and local scale: done

# dynamAedes 2.2.7
* Adapts all auxiliary functions to the `dynamAedesClass` output.
* Adds `adci` function that allows for substages statistics and merges space/non-space outputs.
* Adds `summary` and `max` methods for class `dynamAedesClass`.
* Removes dplyr from vignettes and packages.

# dynamAedes 2.2.6
* Introduces output object as a dynamAedesClass
* Allows for species abbreviation at least 2 characters long (e.g., "ae")

# dynamAedes 2.2.5
* Removes definitely rgdal from package
* Clears NAMESPACE

# dynamAedes 2.2.0
* Substitutes `fields::rdist` with terra::distance which solves a bug in active dispersal
* Adds a switch for active dispersal only `pDispersal=FALSE`
* Removes warnings changing Beta function with polynomial for egg survival in Aedes aegypti
* Raised 0 hatching rate to 15°C for aegypti
* Other minor cosmetic changes
* Modify verbose option for debugging, now an integer with values 0 (silent), 1 or 2.

# dynamAedes 2.1.2
* Transitioning from 'rgdal' and 'sp' to 'terra'
* Transitioning from 'insol' to 'geosphere'
* Removing block-loading of 'tidyverse' (speeding up loading?)
* Helper functions not exported any more

# dynamAedes 2.1.1
* Adding rgdal in the description

# dynamAedes 2.1.0
* Adding website with examples for the three scales (pkgdown)
* Coordinates now internally transformed (if not in lat/long) for the derivation of the photoperiod
* Adding option **coords.proj4** to specify coordinate projection (in proj4 format)
* Renaming dynamAedes as dynamAedes.m

# dynamAedes 2.0.3
* At regional scale now diapause derived pixel by pixel using the photoperiod of input coordinates.
* You can now introduce diapause eggs using dynamAedes argument "intro.degg"

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