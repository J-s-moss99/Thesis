
################################
## 2. von Bertalanffy growth ##
################################


load(file="recapture_data.RData")

# ginc: growth increment (lrec-lrel)
recaptures.rjc$ginc <- recaptures.rjc$lrec - recaptures.rjc$lrel

recaptures.rjc$ginc[recaptures.rjc$ginc<0]
# very small negative growth => measurement error?


######################################################
# exploratory data analysis ##########################
######################################################


library(lattice)

# recapture vs. release length

lmin <- min(c(min(recaptures.rjc$lrel),min(recaptures.rjc$lrec)))
lmax <- max(c(max(recaptures.rjc$lrel),max(recaptures.rjc$lrec)))
xyplot(lrec~lrel,
       groups=experiment,
       xlab=expression(l[rel]),
       ylab=expression(l[rec]),
       recaptures.rjc,type='p',xlim=c(lmin,lmax),ylim=c(lmin,lmax))

# release length and time-at-liberty

xyplot(tlib~lrel,
       groups=experiment,
       xlab=expression(l[rel]),
       ylab=expression(t[lib]),
       recaptures.rjc,type='p')

# empirical growth rate: (lrec-lrel)/tlib  in terms of length-at-release

xyplot(((lrec-lrel)/tlib)~lrel,
       groups=experiment,
       xlab=expression(l[rel]),
       ylab=expression(k[emp]),
       recaptures.rjc,type='p')

# summarise the data

lapply(recaptures.rjc[,c("lrel","lrec","tlib","ginc")],quantile,c(0.025,0.5,0.975))





######################################################
# fit von Bertalanffy (Linf and k NOT t0!) to data ###
######################################################

# predicted growth increment function

growth.inc <- function(l,tau,Linf,k) return((Linf-l)*(1-exp(-k*tau)))

# normal likelihood (to deal with negative increments)
# standard deviation proportional to length-at-release

mod1.loglkhd <- function(theta) {
  
  Linf <- theta[1]
  k <- theta[2]
  phi <- theta[3]
  sigmal <- phi*recaptures.rjc$lrel
  
  # model-predicted growth increments
  
  del.hat <- growth.inc(recaptures.rjc$lrel,recaptures.rjc$tlib,Linf,k)
  
  # log-likelihood
  
  logl <- sum(dnorm(recaptures.rjc$ginc,del.hat,sigmal,TRUE))
  
  return(-logl)
}

# use nlminb() to estimate MLE for Linf, t and phi

initial.guess <- c(110,0.03,0.04)
res.mod1      <- nlminb(initial.guess,mod1.loglkhd,control=list(trace=1))

# did it converge (is convergence flag zero?)

res.mod1[['convergence']]




######################################################
# summary plots ######################################
######################################################

# observed versus predicted length-at-recapture
par(mfrow=c(1,1))
lrec.pred <- recaptures.rjc$lrel+growth.inc(recaptures.rjc$lrel,recaptures.rjc$tlib,res.mod1$par[1],res.mod1$par[2])
ymin <- min(c(min(recaptures.rjc$lrec)),min(lrec.pred))
ymax <- max(c(max(recaptures.rjc$lrec)),max(lrec.pred))
plot(recaptures.rjc$lrel,recaptures.rjc$lrec,ylim=c(ymin,ymax),xlab=expression(l[rel]),ylab=expression(l[rec]),type='p')
points(recaptures.rjc$lrel,lrec.pred,pch=2,col='magenta')
legend('bottomright',pch=c(1,2),col=c("black","magenta"),legend=c("Observed","Predicted"),bty='n')

# fish can't shrink altough there are some observations

# residuals (of growth increment NOT length-at-recapture)

resid.mod1 <- recaptures.rjc$ginc-growth.inc(recaptures.rjc$lrel,recaptures.rjc$tlib,res.mod1$par[1],res.mod1$par[2])

# any pattern in residuals with length-at-release (LOESS smoother)?
# (plotting tip: order data frame in ascending order for length-at-release)

tmp.df <- data.frame(lrel=recaptures.rjc$lrel,res=resid.mod1)
tmp.df <- tmp.df[with(tmp.df, order(lrel)),]
resid.patt <- loess(res~lrel,tmp.df)

plot(tmp.df$lrel,tmp.df$res,xlab=expression(l[rel]),ylab='Residual',type='p')
lines(tmp.df$lrel,resid.patt$fitted,lty=1,col='magenta')
abline(h=0,lty=2)

# are the residual normally distributed?

qqnorm(resid.mod1)
qqline(resid.mod1)


######################################################
# alternative likelihood (Cauchy) ####################
######################################################

# We have fat tails!
# Normal distro for growth increments looks dodgy
# Let's see what difference a Cauchy distribution makes
# Has fatter tails to deal with over-dispersed data

# Cauchy likelihood
# scale parameter (determining variability) proportional to length-at-release

mod2.loglkhd <- function(theta) {
  
  Linf <- theta[1]
  k <- theta[2]
  phi <- exp(theta[3])
  gaml <- phi*recaptures.rjc$lrel
  
  # model-predicted growth increments
  
  del.hat <- growth.inc(recaptures.rjc$lrel,recaptures.rjc$tlib,Linf,k)
  
  # log-likelihood
  
  logl <- sum(dcauchy(recaptures.rjc$ginc,del.hat,gaml,TRUE))
  
  return(-logl)
}

# use nlminb() to estimate MLE for Linf, t and phi

initial.guess <- c(70,0.3,log(0.05))
res.mod2 <- nlminb(initial.guess,mod2.loglkhd,control=list(trace=1))

# did it converge (is convergence flag zero?)

res.mod2[['convergence']]

# residuals

resid.mod2 <- recaptures.rjc$ginc-growth.inc(recaptures.rjc$lrel,recaptures.rjc$tlib,res.mod2$par[1],res.mod2$par[2])


qqnorm(resid.mod2)
qqline(resid.mod2)

######################################################
# what difference does the alternative distro make? ##
######################################################

# difference in parameter estimates

res.mod1$par
res.mod2$par

# is the fit "better" (simple option: look at residual SD)

sd(resid.mod1)
sd(resid.mod2)

######################################################
# generating transition matrix #######################
######################################################

# define the length partition (first lower length class; last upper)
# no rules to define bins, but look to your (catch) data, but be shore that you cover all sizes (L infinity)
Linf <- res.mod1$par[1]
k    <- res.mod1$par[2]
tau  <- 365                # growht increment per day
l0   <- 15 
l    <- l0  
for(t in 1 :15){
  l_old <- l0
  incr <- growth.inc(l0,tau,Linf,k)
  l <- c(l,incr)
  l0 <- l_old + incr
}

plot(1:16,cumsum(l),type="l",col="red")

Linf <- res.mod2$par[1]
k    <- res.mod2$par[2]
tau  <- 365                # growht increment per day
l0   <- 15 
l    <- l0  
for(t in 1 :15){
  l_old <- l0
  incr <- growth.inc(l0,tau,Linf,k)
  l <- c(l,incr)
  l0 <- l_old + incr
}
lines(1:16,cumsum(l),col="blue")


lpar <- c(25,36,45,53,60,71,75,78,80,130)
Nbins <- length(lpar)-1

# define time-period: 1 year in this case
tau <- 1

# using normal likelihood
# size transition matrix, estimate probability to be in a next bin, given you current length
Gamma.normal <- matrix(nrow=Nbins,ncol=Nbins)
for(i in 1:Nbins) {
  for(j in 1:Nbins) {
    # reference length: mid-point of partition interval
    lref <- mean(lpar[i],lpar[i+1])
    lu <- lpar[j+1]
    ll <- lpar[j]
    mul <- lref+growth.inc(lref,tau,res.mod1$par[1],res.mod1$par[2])  # expected new length aft tau
    sigmal <- lref*res.mod1$par[3]
    Gamma.normal[i,j] <- pnorm(lu,mul,sigmal,lower.tail=TRUE)-pnorm(ll,mul,sigmal,lower.tail=TRUE)  # take the part in which you are interested
  }
}
Gamma.normal
rowSums(Gamma.normal)



# using Cauchy likelihood

Gamma.cauchy <- matrix(nrow=Nbins,ncol=Nbins)
for(i in 1:Nbins) {
  for(j in 1:Nbins) {
    # reference length: mid-point of partition interval
    lref <- mean(lpar[i],lpar[i+1])
    lu <- lpar[j+1]
    ll <- lpar[j]
    mul <- lref+growth.inc(lref,tau,res.mod2$par[1],res.mod2$par[2])
    gammal <- lref*exp(res.mod2$par[3])
    Gamma.cauchy[i,j] <- pcauchy(lu,mul,gammal,lower.tail=TRUE)-pcauchy(ll,mul,gammal,lower.tail=TRUE)
  }
}

#######################################################################
# use (normalised) Frobenius norm to compare similarity/dispersion in transition matrices
#######################################################################
# value range from 0 to 1
# Closer to 1 means little spread in transition matrix
# Closer to 0 means length classes get smeared out much more

Frob.norm <- sqrt(sum(diag(Gamma.normal%*%t(Gamma.normal))))/sqrt(Nbins)
Frob.cauchy <- sqrt(sum(diag(Gamma.cauchy%*%t(Gamma.cauchy))))/sqrt(Nbins)

# The Cauchy is fat-tailed (i.e. larger deviations from growth curve especially as length increases) do we see this effect in Frobenius norms?

###################################
# select the optimal length bin => values as close to 1 on the off-diagonal elements

optimal_length_bin_selection <- function(min_l, pars){
  lpar <- min_l+exp(pars)
  Nbins <- length(lpar)-1
  # using normal likelihood
  # size transition matrix, estimate probability to be in a next bin, given you current length
  Gamma.normal <- matrix(nrow=Nbins,ncol=Nbins)
  for(i in 1:Nbins) {
    for(j in 1:Nbins) {
      # reference length: mid-point of partition interval
      lref <- mean(c(lpar[i],lpar[i+1]))
      lu <- lpar[j+1]
      ll <- lpar[j]
      mul <- lref+growth.inc(lref,tau,res.mod1$par[1],res.mod1$par[2])  # expected new length aft tau
      sigmal <- lref*res.mod1$par[3]
      Gamma.normal[i,j] <- pnorm(lu,mul,sigmal,lower.tail=TRUE)-pnorm(ll,mul,sigmal,lower.tail=TRUE)  # take the part in which you are interested
    }
  }
  return(sum((diag(Gamma.normal[-Nbins,-1])-rep(1,(Nbins-1)))^2))
}

pars <- log(10*c(1:9))
res  <- nlminb(pars,optimal_length_bin_selection,min_l=15)
res

lpar <- 15+exp(res$par)

lpar <- lpar[-1]

Nbins <- length(lpar)-1

# define time-period: 1 year in this case
tau <- 1

# using normal likelihood
# size transition matrix, estimate probability to be in a next bin, given you current length
Gamma.normal <- matrix(nrow=Nbins,ncol=Nbins)
for(i in 1:Nbins) {
  for(j in 1:Nbins) {
    # reference length: mid-point of partition interval
    lref <- mean(c(lpar[i],lpar[i+1]))
    lu <- lpar[j+1]
    ll <- lpar[j]
    mul <- lref+growth.inc(lref,tau,res.mod1$par[1],res.mod1$par[2])  # expected new length aft tau
    sigmal <- lref*res.mod1$par[3]
    Gamma.normal[i,j] <- pnorm(lu,mul,sigmal,lower.tail=TRUE)-pnorm(ll,mul,sigmal,lower.tail=TRUE)  # take the part in which you are interested
  }
}
Gamma.normal
rowSums(Gamma.normal)
Gamma.normal[8,8] <- Gamma.normal[8,8] +  (1-rowSums(Gamma.normal)[8])
rowSums(Gamma.normal)






