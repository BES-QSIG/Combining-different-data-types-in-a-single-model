## ISEC workshop: Combining different data types in a single model

This repo contains information and materials for the ISEC workshop 'Combining different data types in a single model' 29th June 2018.
See the (draft) program for the day [here](https://docs.google.com/document/d/1j_arfcNETRIf7pmsJ3fZxuWrI4XWUzdJJjyMIh5bIkc/edit?usp=sharing).
The rest of this page contains information and materials about the three workshop tutorials.

---

### *Rahel Sollman* - **Combining spatial capture-recapture, telemetry, and point count data**

Rahel will talk about how to integrate abundance and individual location data into spatial capture-recapture models to estimate animal abundance.
You can download Rahel's slides [here](https://github.com/BES-QSIG/Combining-different-data-types-in-a-single-model/blob/master/presentations/rahel_sollman_scr_telemetry.pdf).
[This paper](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecs2.2203) by Linden *et al.* is another good reference.
Rahel will walk us through how to carry out these analyses using the `oSCR` package ([Github](https://github.com/jaroyle/oSCR)).

#### installation

You can download the most recent version of `oSCR` from GitHub, using the `devtools` package (which is on CRAN), and you'll also need to manually install the `car` package from CRAN.

```r
# install devtools if you don't already have it installed, with:
# install.packages("devtools")
devtools::install_github("jaroyle/oSCR")
install.packages("car")
```
If you're a Windows user and this the first time you've tried installing a package from GitHub, you may need to install the [RTools software](https://cran.r-project.org/bin/windows/Rtools/) (note: this is not an R package), which provides the tools you need to build packages from source code.

Rahel's example also requires the `scrbook` package for some plotting functions. `scrbook` is hosted on [the developers' own website](https://sites.google.com/site/spatialcapturerecapture/scrbook-r-package), and you should also be able to install it with the following commmands:
```r
ext <- ifelse(.Platform$OS.type == "windows", ".zip", ".tar.gz")
scrbook_url <- paste0("https://sites.google.com/site/spatialcapturerecapture/scrbook-r-package/scrbook_0.28-2",
                      ext)
install.packages(scrbook_url, repos = NULL)
```

#### examples

Rahel's example will use synthetic data representing species counts from a spatial capture-recapture design, combined with data on the positions of individuals from telemetry.
You can **download the code for the example analysis [here](https://raw.githubusercontent.com/BES-QSIG/Combining-different-data-types-in-a-single-model/master/scripts/sollmann_workshop_script.R)**
(you can either copy-paste the code into a file, or go to `save page as ...` in your browser)

---

### *Brett McClintock* - **Combining animal movement and environmental data with hidden Markov models**

Brett will be talking about using Hidden Markov Models to combine movement data with other data streams to infer behavioural states of animals.
You can download Brett's slides [here](https://github.com/BES-QSIG/Combining-different-data-types-in-a-single-model/blob/master/presentations/brett_mcclintock_hmm.pdf).
Brett and Theo Michelot developed the `momentuHMM` package ([CRAN](https://CRAN.R-project.org/package=momentuHMM), [GitHub](https://github.com/bmcclintock/momentuHMM)) to make it easier to do these types of analyses.
A paper introducing momentuHMM is freely available at Methods in Ecology and Evolution ([pdf](https://besjournals.onlinelibrary.wiley.com/doi/epdf/10.1111/2041-210X.12995)).

#### installation

To get started, you will need to load R and install the `momentuHMM` package from CRAN.

```r
install.packages("momentuHMM")
```

To run the main example, you'll also need to install the setRNG package from CRAN.

```r
install.packages("setRNG")
```

#### examples

The main example Brett will work through is an analysis of Northern Fur Seal foraging trips with position uncertainty.
You can **download the code for the Northern Fur Seal example analysis [here](https://raw.githubusercontent.com/bmcclintock/momentuHMM/master/vignettes/nfsExample.R)**
(you can either copy-paste the code into a file, or go to `save page as ...` in your browser)

There are a couple of more advanced example analyses that you can download and play with, and that Brett might have time to talk about:
 - Modelling Loggerhead Turtle movements relative to ocean currents with position uncertainty ([link](https://raw.githubusercontent.com/bmcclintock/momentuHMM/master/vignettes/turtleExample.R))
 - Modelling Northern Fulmar movements relative to fishing vessels ([link](https://raw.githubusercontent.com/bmcclintock/momentuHMM/master/vignettes/northernFulmarExample.R)) 

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

You'll need the following packages from CRAN to run the example too:
```r
install.packages(c("spatstat", "RColorBrewer", "mapview"))
```

#### examples

Bob will work through a species distribution modelling analysis of New Zealand's Solitary Tinamou, using distribution data from multiple sources.
This analysis is provided as a vignette in the package, and **you can download the code for the Solitary Tinamou analysis [here](https://raw.githubusercontent.com/goldingn/PointedSDMs/master/vignettes/SolitaryTinomou.R)** 


