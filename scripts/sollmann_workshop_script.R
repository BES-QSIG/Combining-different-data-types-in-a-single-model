###This is an R script that demonstrates how to integrate telemetry data with SCR models using oSCR.
###It is based on code developed by A. Royle and D. Linden (https://rdrr.io/github/jaroyle/oSCR/man/telemetry.html)
###The original ML code is a modification from the Supplement in 
###Royle, J. Andrew, Richard B. Chandler, Catherine C. Sun, and Angela K. Fuller. "Integrating resource selection information with spatial capture-recapture." 
###Methods in Ecology and Evolution 4, no. 6 (2013): 520-530.

##had to install these packages separately

#install.packages("devtools")
#install.packages("car")
#install.packages("FNN")

library(raster)
library(scrbook)
library(FNN)
library(devtools)

#install_github("jaroyle/oSCR")
library(oSCR)

#################################################################################################
###### PART 1 - Telemetry in sparse data situations ##############################################

### discrete state-space; objective: estimate density
N<-15  ##true abundance, small population

##dimensions and discrete cells of state space S
max.x<-max.y<-25
min.x<-min.y<-1
gr<-expand.grid(min.x:max.x,min.y:max.y)
G<-dim(gr)[1]

##distance between cells within S
Dmat<-as.matrix(dist(gr))

###place 49 sampling devices in the middle
# Make a trap array
X<- cbind( sort(rep( seq(7,19,2),7)), rep( seq(7,19,2),7))
J<-nrow(X)

# This maps the trap locations to the raster cells
raster.point<-rep(NA,nrow(X))
for(j in 1:nrow(X)){ 
  raster.point[j]<- (1:dim(gr)[1])[ (X[j,1]==gr[,1]) & (X[j,2] == gr[,2])]
}

sigma<-1.5   ##half-normal scale parameter, here in units of grid cell resolution
lam0<-0.1    ##baseline trap encounter rate
K<-5         ##number of sampling occasions


###generate activity centers
set.seed(2078)
s<-sample(1:G, N)


###plot everything
plot(c(1,25), c(1,25), type="n", xlab="X", ylab="Y")
points(gr[s,1], gr[s,2], pch=19, col="grey") #activity centers
points(X, pch=19) #detectors


#### generate SCR data from hair snare type detector (max 1 detection per trap and occasion)
### for analysis in oSCR, data need to be 3D: individual by trap by occasion

##distance from activity centers to all traps
dist<-e2dist(gr[s,], X)

det<-array(NA, c(N,J,K))

for (i in 1:N){
  
  for (j in 1:J){
    #loglam<-log(lam0)-dist[i,j]^2/(2*sigma^2)
    #p<-1-exp(-exp(loglam))
    lam<-lam0 * exp(-dist[i,j]^2/(2*sigma^2))  ##expected detection rate
    p<-1-exp(-lam)                            ##cloglog link - converts rate to probability
    det[i,j,]<-rbinom(K, 1, p)                ## generate detections based on Binomial model
  }
}

nobs<-sum(apply(det,1,sum)>0)  ##number detected
##if nobs=0, model will not run

###subset detection matrix to detected animals only
who<-which(apply(det,1,sum)>0)
obs<-det[who,,]

##check at how many traps animals were captured
##no spatial recaptures --> model will not run
apply(obs, 1, function(x){y<-apply(x, 1,sum)
return(sum(y>0))})

# nt<-NULL
# for (i in 1:9){
#   sub<-obs[i,,]
#   y<-apply(sub, 1, sum)
#   nt[i]<-sum(y>0)
# }

###generate telemetry data, for small number of individuals only

Ntel<-3
n.locs<-10

###pick individuals in the center of S so as to not truncate locations
###In real life analysis, you'd increase S to include all locations
poss.tel<- gr[s,1]>(0+4*sigma) & gr[s,1]<(25-4*sigma) & gr[s,2]>(0+4*sigma) & gr[s,2]<(25-4*sigma)
tel.guys<-sort(sample(1:N,Ntel, prob=as.numeric(poss.tel)))
sid<-s[tel.guys]      #activity center grid cells for telemetered guys
cap.tel <- tel.guys   #which row in original capture history for each telemetry guy
stel<-gr[sid,]        ##activity center coordinates for telemetered guys


###generate telemetry locations
tel.locs<-matrix(NA, Ntel*n.locs, 2)

for (i in 1:Ntel){
  tel.locs[(i*n.locs - (n.locs-1)) : (i*n.locs),1]<-rnorm(n.locs, stel[i,1], sigma)
  tel.locs[(i*n.locs - (n.locs-1)) : (i*n.locs), 2]<-rnorm(n.locs, stel[i,2], sigma)
  ##equivalent to bivariate normal distribution of locations around activity center
}

### add these to plot

for (i in 1:Ntel){
  points(tel.locs[(i*n.locs - (n.locs-1)) : (i*n.locs),], pch=19, cex=0.7, col=i+1)
  points(stel[i,1], stel[i,2], pch=8,  col=i+1)
}


##########################################################################################
####### format and bundle data for analysis in oSCR ######################################

##convert telemetry location matrix to data frame
teldat<-data.frame(animalid=rep(1:Ntel, each=n.locs), 
                   X= tel.locs[,1],  
                   Y= tel.locs[,2])

## create state space data frame (for oSCR)
ssDF <- data.frame(X=gr[,1],Y=gr[,2])

##process telemetry data to obtain number of telemetry locations in each grid cell of the
##state-space (for each individual) - use oSCR function telemetry.processor
telem<-telemetry.processor(list(ssDF),list(teldat) )
str(telem)


### this is a little clunky but the way oSCR is set up, we need telemetry 
### information of detected individuals first, then of non-detected individuals

##extract locations per grid cell
nfreq=telem$nfreq[[1]]

##identify which telemetered animal is in obs - "who" is index of which simulated animal was
## detected 
bb<-pmatch(tel.guys, who)
## if bb is all NA (no telemetered animal was detected, telemetry data needs to be set up without
## cap.tel)
if (sum(!is.na(bb))==0) 
  print("No telemetered animal detected! Build telemetry object without cap.tel")

##change order of telemetry frequencies so that animals in obs come first
new.freq<-rbind(nfreq[!is.na(bb),], nfreq[is.na(bb),])

##bb also gives you the new position of telemetry guys in obs (excludig NA)
cap.tel<-bb[!is.na(bb)]

##with that, we can compile telemetry data
telemetry<-list(fixfreq=list(new.freq), cap.tel=list(cap.tel))

# Create the scrFrame 
sftel <- make.scrFrame(caphist = list(obs),
                       traps = list(data.frame(X=X[,1],Y=X[,2])), 
                       telemetry = telemetry)
##if you have no spatial recaptures, a warning message will appear
##there will always be a warning (if you use cap.tel) about the order of telemetry information

##model without telemetry
fit1 <- oSCR.fit(scrFrame=sftel,ssDF=list(ssDF),DorN="D",encmod="CLOG",
                 trimS=5*sigma,
                 model=list(D~1,p0~1,sigma~1))

###parameter estimates on link scale
fit1$outStats
##sigma on log-scale
##p0 on log scale for CLOG model
##d0 density on log scale

###back-transform estimates 
predS <- get.real(fit1,type = "sig")
predp <- get.real(fit1,type = "det")  #this takes a minute
predp[[1]][[1]][1,] ##baseline detection, constant in all pixels
predD <- get.real(fit1,type = "dens")  #this takes a minute
predD[[1]][1,] ##density per pixel, constant in all pixels 

####predict density across state-space
dsurf<-predict.oSCR(fit1, sftel, list(ssDF), override.trim = T)

plot(dsurf$r[[1]])

##

####### fit the same model with telemetry data, considering both data
####### sets as independent

fit2 <- oSCR.fit(scrFrame=sftel,ssDF=list(ssDF),DorN="D",encmod="CLOG",
                 trimS=5*sigma,
                 telemetry="ind",
                 model=list(D~1,p0~1,sigma~1))
fit2$outStats
fit1$outStats
##estimates of detection parameters improved, density slightly improved

dsurf2<-predict.oSCR(fit2, sftel, list(ssDF), override.trim = T)
##note the output of estimated N in S

plot(dsurf2$r[[1]])
## estimates of activity centers improved


#######fit the same model with telemetry data, considering both data
# ####### sets as dependent
# fit3 <- oSCR.fit(scrFrame=sftel,ssDF=list(ssDF),DorN="D",encmod="CLOG",
#                  trimS=5*sigma,
#                  telemetry="dep",
#                  model=list(D~1,p0~1,sigma~1))
# fit3$outStats




###########################################################################################
#### analyze simulated data with spatial covariate ########################################

source("DataSimFun.R")
## This function simulates SCR and telemetry data as above with the specified input on a 
## discrete state-space of size 37 by 37, with 49 traps, spaced by 2 units, in the middle.
## It also simulates a spatially correlated covariate.
## Note that the state space is adequate for sigma =2 but would be too small for larger 
## values of sigma

alpha0 <- -1 #log(lam0)
sigma<- 2
alpha2<- 1   # This is the effect of the covariate on resource selection
Ntel<-4      # number of individuals with telemetry devices
Nfixes<-20   # number of telemetry fixes per individual
N<- 100      # population size
K<-3         # number of occasions

set.seed(1234)

dat<-generateData(alpha0, sigma, alpha2, Ntel, Nfixes, N, K)

##simple plot of covariate elevation
par(mar=c(5,4,2,5))
spatial.plot(dat$gr,dat$z,cx=3, col="red")

##add trap locations to plot
points(dat$X, pch=19)

##add telemetry locations to plot
points(dat$gr[dat$locs[1,]>0, ], pch=19, col="purple")
points(dat$gr[dat$locs[2,]>0, ], pch=19, col="red")
points(dat$gr[dat$locs[3,]>0, ], pch=19, col="blue")
points(dat$gr[dat$locs[4,]>0, ], pch=19, col="darkgrey")

##look at SCR capture locations
spiderplot(dat$y, dat$X)

###### Compile data for analysis
ntraps<-dim(dat$X)[1]
trapCovs <- list(z=list(data.frame(z=matrix(dat$z[dat$raster.point],ntraps,K))))
trapCovs <- make.trapCovs(trapCovs)

# Set up the ssDF and the rsfDF
ssDF <- rsfDF <- data.frame(X=dat$gr[,1],Y=dat$gr[,2],z=dat$z)

# temeletry data
telemetry <- list(fixfreq=list(dat$locs), cap.tel=list(dat$cap.tel))

# Create the scrFrame 
sftel <- make.scrFrame(caphist = list(dat$y),
                       traps = list(data.frame(X=dat$X[,1],Y=dat$X[,2])), 
                       trapCovs = trapCovs,
                       telemetry = telemetry)
##Note: you can also not provide trapCovs; they will be extracted automatically
##      if you provide an rsfDF instead
##The simulation function returns telemetry data in the right order


# fit the SCR model with NO telemetry integration, no z
# these models will take a few minutes to run (2 mins on my laptop)

fit0z <- oSCR.fit(scrFrame=sftel,ssDF=list(ssDF),DorN="D",encmod="CLOG",
                 trimS=7,
                 model=list(D~1,p0~1,sigma~1))

# fit the SCR model with NO telemetry integration (z from traps only)
fit1z <- oSCR.fit(scrFrame=sftel,ssDF=list(ssDF),DorN="D",encmod="CLOG",
                 trimS=7,
                 model=list(D~1,p0~z,sigma~1))

# Next we use the telemetry information to inform about both beta.z (RSF=TRUE) and
# also 'sigma', assuming independence between data (captures vs. collars)
fit2z <- oSCR.fit(scrFrame=sftel,ssDF=list(ssDF),DorN="D",encmod="CLOG",
                 rsfDF=list(rsfDF),RSF=TRUE,telemetry="ind",
                 trimS = 7,
                 model=list(D~1,p0~z,sigma~1))

# Now we use the telemetry information again (RSF = TRUE) but assume dependence
# between data since some collared guys were captured
fit3z <- oSCR.fit(scrFrame=sftel,ssDF=list(ssDF),DorN="D",encmod="CLOG",
                 rsfDF=list(rsfDF),RSF=TRUE,telemetry="dep",
                 trimS = 7,
                 model=list(D~1,p0~z,sigma~1))


# Here we fit the SCR model with RSF = FALSE, which only uses the
# telemetry data to inform about 'sigma' NOT the RSF parameters
fit4z <- oSCR.fit(scrFrame=sftel,ssDF=list(ssDF),DorN="D",encmod="CLOG",
                 rsfDF=list(rsfDF),RSF=FALSE,telemetry="ind",
                 trimS = 7,
                 model=list(D~1,p0~z,sigma~1,path~1))

##no telemetry, no estimation of effect of z
fit0z$outStats

##no telemetry
fit1z$outStats

##telemetry for sigma and space use
fit2z$outStats

##as above, dependent data sets
fit3z$outStats

##telemetry for sigma
fit4z$outStats
                       
###backtransformed estimates and density plot from fit3          
predS <- get.real(fit3z,type = "sig")
predS
predp <- get.real(fit3z,type = "det")  #this takes a minute
predp[[1]][[1]][1:3,] ##baseline detection, varies across pixels
predD <- get.real(fit3z,type = "dens")  #this takes a minute
predD[[1]][1,] ##density per pixel

####predict density across state-space
dsurfz<-predict.oSCR(fit3z, sftel, list(rsfDF), override.trim = T)

plot(dsurfz$r[[1]])
       
#### for model comparison
fit3z$AIC


 


##############################################################################################
############ NY bear example, within oSCR ####################################################
### based on telemetry() helpfile
### models take a little too long for this workshop 

## load data
data("nybears")
attach(nybears)
str(nybears)

# Note: 8/21/2016 some of the data objects in this data object are not correct and
#  will be updated at a later date, but the elements used here are correct.
colnames(teldata2)<- c("ind","X","Y")
colnames(ssgrid)<- c("X","Y")

ntraps<- length(trap2raster)
# nybears$traplocs is wrong... gotta do this:
traplocs <- ssgrid[trap2raster,]
ntraps <- nrow(traplocs)

##simple plot of covariate elevation
par(mar=c(5,4,2,5))
spatial.plot(ssgrid,elevation,cx=2, col="red")

##add trap locations to plot
points(traplocs, pch=19)

##add telemetry locations to plot
points(ssgrid[ntel[1,]>0, ], pch=19)
points(ssgrid[ntel[2,]>0, ], pch=19, col="red")
points(ssgrid[ntel[3,]>0, ], pch=19, col="blue")

##extract elevation at each trap
trapCovs <- list(z=list(data.frame(z=matrix(elevation[trap2raster],ntraps,K))))
trapCovs <- make.trapCovs(trapCovs)

##make state space data frame, RSF data frame (identical)
ssDF <- rsfDF <- data.frame(ssgrid,z=elevation)

# distribute the binomial captures (n x ntraps) among the K surveys
# could simply simulate a binary array of captures (n x ntraps x K)
# only possible when there is no variation in p0 among surveys

y.arr <- array(0,dim=c(nrow(y2d),ntraps,K))
for(i in 1:nrow(y2d)){
  for(j in 1:ntraps){
    which.K <- sample(1:K,y2d[i,j])
    y.arr[i,j,which.K] <- 1
  }
}

# create telemetry list using precalculated fix frequencies 
# order ntel so that captured animals come first
# cap.tel indicates location in the observation data
telemetry.bears <- list(fixfreq=list(ntel[3:1,]),cap.tel=list(c(14, 26)))

###compile data in scrFrame
sftel.bears <- make.scrFrame(caphist = list(y.arr),
                       traps = list(traplocs), 
                       trapCovs = trapCovs,
                       telemetry = telemetry.bears)

######################## fit a series of model #####################################
###fit model, no telemetry information, no trap covariate
fit0 <- oSCR.fit(scrFrame=sftel.bears,ssDF=list(ssDF),DorN="D",encmod="CLOG",
                model=list(D~1,p0~1,sigma~1))

###fit model, no telemetry information, with trap covariate
fit1 <- oSCR.fit(scrFrame=sftel.bears,ssDF=list(ssDF),DorN="D",encmod="CLOG",
                 model=list(D~1,p0~z,sigma~1))

###fit model, telemetry only informs sigma
fit2 <- oSCR.fit(scrFrame=sftel.bears,ssDF=list(ssDF),DorN="D",encmod="CLOG",
                rsfDF=list(rsfDF),RSF=FALSE,telemetry="ind",
                model=list(D~1,p0~z,sigma~1))

###fit model, telemetry informs RSF and sigma, data sets are independent
fit3 <- oSCR.fit(scrFrame=sftel.bears,ssDF=list(ssDF),DorN="D",encmod="CLOG",
                rsfDF=list(rsfDF),RSF=TRUE,telemetry="ind",
                model=list(D~1,p0~z,sigma~1))

###fit model, telemetry informs RSF and sigma, dependency between data sets
fit4 <- oSCR.fit(scrFrame=sftel.bears,ssDF=list(ssDF),DorN="D",encmod="CLOG",
                rsfDF=list(rsfDF),RSF=TRUE,telemetry="dep",
                model=list(D~1,p0~z,sigma~1))


##############################################################################################















