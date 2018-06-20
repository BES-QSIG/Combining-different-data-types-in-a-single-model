## ISEC workshop: Combining different data types in a single model

This repo contains information and materials for the ISEC workshop 'Combining different data types in a single model' 29th June 2018


## program
Here's a link to the draft program for the day: [link](https://docs.google.com/document/d/1j_arfcNETRIf7pmsJ3fZxuWrI4XWUzdJJjyMIh5bIkc/edit?usp=sharing)

## materials

### *Rahel Sollman* - **Combining spatial capture-recapture, telemetry, and point count data**

<materials to be added>

### *Brett McClintock* - **Combining animal movement and environmental data with hidden Markov models**

Brett will be talking about combining movement data with other data streams to infer behavioural states of animals.
Brett and Theo Michelot developed the `momentuHMM` package ([CRAN](https://CRAN.R-project.org/package=momentuHMM), [GitHub](https://github.com/bmcclintock/momentuHMM)) to make it easier to do these types of analyses.
A paper introducing momentuHMM is freely available at Methods in Ecology and Evolution [pdf](https://besjournals.onlinelibrary.wiley.com/doi/epdf/10.1111/2041-210X.12995)

The main example Brett will work through is an analysis of Northern Fur Seal foraging trips.
To get started, you will need to load R and install the `momentuHMM` package from CRAN:

```r
install.packages("momentuHMM")
```

You can download the code for the Northern Fur Seal example analysis [here](https://raw.githubusercontent.com/bmcclintock/momentuHMM/master/vignettes/nfsExample.R)
(you can either copy-paste the code into a file, or go to `save page as ...` in your browser)

There are a couple of more advanced example analyses that you can download and play with, and that Brett might have time to talk about:
 - Turtles [link](https://raw.githubusercontent.com/bmcclintock/momentuHMM/master/vignettes/turtleExample.R)
 - Northern Fulmars [link](https://raw.githubusercontent.com/bmcclintock/momentuHMM/master/vignettes/northernFulmarExample.R) 

### *Bob O’Hara* - **Combining presence-only, presence-absence and abundance data with point-process models**

<materials to be added>
