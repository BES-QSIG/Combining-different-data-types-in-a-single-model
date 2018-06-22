##### write data simulation function

##### This function generates telemetry and SCR data on a discrete state-space
##### of size 37 by 37, with 7 traps, spaced by 2 units, in the middle
##### The state-space has a spatially correlated covariate, z, that influences
##### space use. SCR data simulated under a cloglog model
##### The function takes as input:
## alpha0: log(baseline detection rate)
## sigma: movement/scale parameter of detection function
## alpha2: effect of covariate on space use/detection
## Ntel: number of individuals with telemetry data
## Nfixes: number of independent telemetry locations per collared individual
## N: total population size
## K: number of sampling occasions
##### The function returns:
## z: covariate, one value for each grid cell of S
## cap.tel: The index of which rows in the detection array the telemetered 
##         individuals belong to
## locs: Ntel by grid cell matrix with telemetry frequency 
## y: observation array, n by traps by K
## X: trap matrix
## raster.point: grid cell each trap locations falls in

generateData<-function(alpha0, sigma, alpha2, Ntel, Nfixes, N, K){

nx<-ny<-37
gr<-expand.grid(1:nx,1:ny)
gr<-as.matrix(gr)
Dmat<-as.matrix(dist(gr))
V<-exp(-Dmat/5)
# Below change to matrix multiplication
z<-as.vector( crossprod( t(chol(V)), rnorm(dim(gr)[1]) ) )

# Simulate activity centers of all N individuals in the population
Sid<- sample(1:dim(gr)[1],N,replace=TRUE)
# and coordinates
S<-gr[Sid,]

poss.tel<- gr[Sid,1]>(1+4*sigma) & gr[Sid,1]<(nx-4*sigma) &
          gr[Sid,2]>(1+4*sigma) & gr[Sid,2]<(ny-4*sigma)
tel.guys<-sort(sample(1:N,Ntel, prob=as.numeric(poss.tel)))
sid<-Sid[tel.guys]
stel<-gr[sid,]

n<-matrix(NA,nrow=Ntel,ncol=dim(gr)[1])

# for each telemetered guy simulate a number of fixes.
# note that nfix = 0 for most of the landscape pixels

lammat<-matrix(NA,nrow=Ntel,ncol=dim(gr)[1])
for(i in 1:Ntel){
  d<- Dmat[sid[i],]
  lam<- exp(1 - (1/(2*sigma*sigma))*d*d + alpha2* z)
  n[i,]<-rmultinom(1,Nfixes,lam/sum(lam))
  lammat[i,]<-lam
}

# Make a trap array
X<- cbind( sort(rep( seq(7,31,4),7)), rep( seq(7,31,4),7))
ntraps<-nrow(X)

raster.point<-rep(NA,nrow(X))
# This just maps the trap locations to the raster cells
for(j in 1:nrow(X)){ 
  raster.point[j]<- (1:dim(gr)[1])[ (X[j,1]==gr[,1]) & (X[j,2] == gr[,2])]
}

# Hazard model is used. This seems the most sensible. 
D<- e2dist(S,X) ## N x ntraps
Zmat<- matrix(z[raster.point],nrow=N,ncol=ntraps,byrow=TRUE) # note make dims the same
loglam<- alpha0 -(1/(2*sigma*sigma))*D*D + alpha2*Zmat
p<- 1-exp(-exp(loglam))
# Now simulate SCR data
y<-array(NA,c(N,ntraps, K))
for(i in 1:N){
  for (k in 1:K){
  y[i,,k]<- rbinom(ntraps,1, p[i,])
  }
}

# Subset data to captured individuals
cap<-which(apply(y,1,sum)>0)
y<-y[cap,,]

bb<-pmatch(tel.guys, cap)
## if bb is all NA (no telemetered animal was detected, telemetry data needs to be set up without
## cap.tel)

##change order of telemetry frequencies so that animals in obs come first
new.freq<-rbind(n[!is.na(bb),], n[is.na(bb),])

##bb also gives you the new position of telemetry guys in obs (excludig NA)
cap.tel<-bb[!is.na(bb)]
  
return(list(z=z, cap.tel=cap.tel, locs=new.freq, y=y, X=X, raster.point=raster.point, gr=gr))
     
} 




