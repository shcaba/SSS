## Simple Stock Synthesis 

Simple Stock Synthesis (SSS) is an assessment method for application to data-limited stocks that estimates catch limits.
It is an age-structured version of other catch-only methods such as DBSRA and CMSY.

The "background literature" folder contains papers on the methods origin and development. 

An GUI version of SSS is available using the SS-DL tool (https://github.com/shcaba/SS-DL-tool).
The SS-DL tool incorporates the SSS approach into a suite of methods ranging from data-limited to data-rich methods all in one framework/package.

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
