## ISEC workshop: Combining different data types in a single model

This repo contains information and materials for the ISEC workshop 'Combining different data types in a single model' 29th June 2018.
See the (draft) program for the day [here](https://docs.google.com/document/d/1j_arfcNETRIf7pmsJ3fZxuWrI4XWUzdJJjyMIh5bIkc/edit?usp=sharing).
The rest of this page contains information and materials about the three workshop tutorials.

---

### *Rahel Sollman* - **Combining spatial capture-recapture, telemetry, and point count data**

Rahel will talk about how to integrate abundance and individual location data into spatial capture-recapture models to estimate animal abundance.
Rahel will walk us through how to carry out these analyses using the `oSCR` package ([Github](https://github.com/jaroyle/oSCR)).

#### installation

You can download the most recent version of `oSCR` from GitHub, using the `devtools` package (which is on CRAN).

```r
# install devtools if you don't already have it installed, with:
# install.packages("devtools")
devtools::install_github("jaroyle/oSCR")
```
If you're a Windows user and this the first time you've tried installing a package from GitHub, you may need to install the [RTools software](https://cran.r-project.org/bin/windows/Rtools/) (note: this is not an R package), which provides the tools you need to build packages from source code.

#### examples

[ materials to be added ]

---

### *Brett McClintock* - **Combining animal movement and environmental data with hidden Markov models**

Brett will be talking about using Hidden Markov Models to combine movement data with other data streams to infer behavioural states of animals.
Brett and Theo Michelot developed the `momentuHMM` package ([CRAN](https://CRAN.R-project.org/package=momentuHMM), [GitHub](https://github.com/bmcclintock/momentuHMM)) to make it easier to do these types of analyses.
A paper introducing momentuHMM is freely available at Methods in Ecology and Evolution ([pdf](https://besjournals.onlinelibrary.wiley.com/doi/epdf/10.1111/2041-210X.12995)).

#### installation

To get started, you will need to load R and install the latest version of the `momentuHMM` package from CRAN.

```r
install.packages("momentuHMM")
```

To run the main example, you'll also need to install the setRNG package from CRAN.

```r
install.packages("setRNG")
```

#### examples

The main example Brett will work through is an analysis of Northern Fur Seal foraging trips.
You can download the code for the Northern Fur Seal example analysis [here](https://raw.githubusercontent.com/bmcclintock/momentuHMM/master/vignettes/nfsExample.R)
(you can either copy-paste the code into a file, or go to `save page as ...` in your browser)

There are a couple of more advanced example analyses that you can download and play with, and that Brett might have time to talk about:
 - Modelling turtles with position uncertainty ([link](https://raw.githubusercontent.com/bmcclintock/momentuHMM/master/vignettes/turtleExample.R))
 - Modelling Northern Fulmars with position uncertainty ([link](https://raw.githubusercontent.com/bmcclintock/momentuHMM/master/vignettes/northernFulmarExample.R)) 

---

### *Bob O’Hara* - **Combining presence-only, presence-absence and abundance data with point-process models**

Bob will talk about combining different types of distribution data together in a single point process model of a species' distribution.
Bob has developed a package `PointedSDMS` ([GitHub](https://github.com/oharar/PointedSDMs)) to fit these models using the [INLA](http://www.r-inla.org/) inference software.

#### installation

You can install `PointedSDMs` from GitHub:

```r
devtools::install_github("oharar/PointedSDMs")
```

You'll also need to install INLA from the INLA team's own repository, using the following code:

```r
install.packages("INLA", repos = "https://inla.r-inla-download.org/R/stable"))
```

#### examples

[ materials to be added ]
