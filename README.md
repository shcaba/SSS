## Simple Stock Synthesis 

Simple Stock Synthesis (SSS) is an assessment method for application to data-limited stocks.  
*Code available on this website comes with no warranty or guarantee of accuracy.*

## Installation

```S
install.packages("devtools")
library(devtools)

devtools::install_github("shcaba/SSS", build_vignettes = TRUE)
```

Note: devtools may give this message: "*WARNING: Rtools is required to build R packages, but is not currently installed.*" However, Rtools is NOT required for installing sss via devtools, so ignore the warning.

Once you have installed the sss package, it can be loaded using:

```S
library(sss)
```
A vignette is not yet available but is being created and will be available using:

```S
vignette(sss)
```
