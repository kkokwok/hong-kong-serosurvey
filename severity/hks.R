# Copyright Steven Riley (sr@stevenriley.net)

# This file is part of the library idsource.

# idsource is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This work is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this idsource.  If not, see <http://www.gnu.org/licenses/>.

# Function to generate a dummy serosurvey protocol
hks.gen.dum.prot <- function(N=400,ave.t1=50,ave.t2=220,ave.t3=420,ave.t4=600,pmean=15) {
  rtn <- data.frame(	UID=1:N,
                     t1=ave.t1-ceiling(pmean/2) + rpois(N,pmean),
                     t2=ave.t2-ceiling(pmean/2) + rpois(N,pmean),
                     t3=ave.t3-ceiling(pmean/2) + rpois(N,pmean),
                     t4=ave.t4-ceiling(pmean/2) + rpois(N,pmean))
  rtn
}

# Function to simulate a vector of N from an arbitrary density of infections
# This needs to be changed to allow some people not to be infected
# Needs an extra parameter
hks.sim.incidence.days <- function(N,vecInc,overallp=1) {
  nsteps <- length(vecInc)
  cuminc <- vector(mode="numeric",length=nsteps)
  cuminc[1] <- vecInc[1]
  for (i in 2:nsteps) cuminc[i] <- cuminc[i-1]+vecInc[i]
  totalinc <- cuminc[nsteps]
  rtn <- vector(mode="numeric",length=N)
  for (i in 1:N) {
    if (runif(1) < overallp) {
      inf <- ceiling(runif(1)*totalinc)
      day <- 1
      while (cuminc[day] < inf) day <- day + 1
      rtn[i] <- day
    } else {
      rtn[i] < -1
    }
  }
  rtn
}

# Another simulation function to take a set of dt and number and then
# return a random deviate drawn according to the temporal distribution
# Probability between starts is linearly interpolated
# Can use the function above to draw the months first then choose ranomly from there.
hks.sim.incidence.months <- function(N,vecStarts,vecEnds,vecInc) {
  months <- hks.sim.incidence.days(N,vecInc)
  rtn <- months
  for (i in 1:N) {
    month <- months[i]
    rtn[i] <- vecStarts[month] + (vecEnds[month] - vecStarts[month]) * runif(1)
  }
  rtn
}

# Simulate HI readings given day of blood draws, times of infection.
# Assume everyone starts from zero
# Takes a matrix of times of samples and number of people, vector of the times of infection, single parameter for boosting and singel parameter for waning
# Assumes all are zero at time t=0
hks.sim.study <- function(matTimes,vecInfs,t_min,t_max,fnP,t_r=21,t_year=365.25) {

  # Declare required auxilliary variables and make some checks
  tmpdim <- dim(matTimes)
  noPeople <- tmpdim[1]
  noSample <- tmpdim[2]
  vecY0 <- rep(0,noPeople)
  matP <- read.csv(fnP,row.names=1)
  shp <- 10

  # Setup the return variable
  rtn <- matTimes

  # Start main loop
  for (i in 1:noPeople) {

    # Assign person specific variables
    y0 	<- vecY0[i]
    ti 	<- vecInfs[i]
    b 	<- matP["boost","val"]
    w 	<- matP["wane","val"]

    # Start time loop
    for (j in 1:noSample) {

      # Assign time specific variables
      t_meas <- matTimes[i,j]
      y <- hks.time.func.1(t_meas,y0,ti,b,w,t_r,t_year)

      # Simulate from the modified poisson
      # tmp <- 0
      # while (tmp < 1) tmp <- rgamma(1,shape=shp,scale=(y+1)/shp)
      # tmp <- tmp - 1
      # rtn[i,j] <- floor(tmp)
      # rtn[i,j] <- srg.mod.gamma.rng(1,y,shp)
      rtn[i,j] <- srg.titre.rng(1,y)

    }

  }

  # Return in the same format as the actual data so that swapping out real and sim is easy
  chain <- list(vt=vecInfs,y0=vecY0,ps=matP)
  list(mt=matTimes,mh=rtn,chain=chain)

}


hks.test.foldrise <- function(vect1,vect2,vecinf,td=0) {
  npeople <- length(vect1)
  rtn <- vector(length=npeople)
  for (i in 1:npeople) {
    if (((vecinf[i]-td) > vect1[i]) & ((vecinf[i]-td) <= vect2[i])) rtn[i] <- TRUE
    else rtn[i] <- FALSE
  }
  rtn
}

# Function to simulate the model that we will use to test the inference methods for
# multiannual serological data
hks.simulate.infection <- function(
  dumdat = hks.gen.dum.prot(
    N=400,
    ave.t1=50,
    ave.t2=250,
    ave.t3=150,
    ave.t4=400,
    pmean=15
  ),
  pairdefs = list(c("t1","t2"),c("t1","t3"),c("t1","t4")),
  single_tmax = 100,
  t_rise = 0,
  tmax = 360*2,
  overallp = 0.6,
  t01=100,
  R01=1.6,
  Tg1=2.6,
  t02=300,
  R02=2.5,
  Tg2=2.6,
  p1=0.5,
  p2=0.3,
  waveBounds = c(1,250,400)
) {

  # Set up some auxillary variables
  npeople <- dim(dumdat)[1]

  # Draw infection events from the waves
  # Simulate the two epidemics for the period of time required
  mod1 <- tex.seir.basic(De=1.48,Tg=Tg1,R0=R01,I0=1,dt=1,noTimeSteps=single_tmax,deterministic=TRUE)
  mod2 <- tex.seir.basic(De=1.48,Tg=Tg2,R0=R02,I0=1,dt=1,noTimeSteps=single_tmax,deterministic=TRUE)
  siminc <- vector(mode="numeric",length=tmax)
  siminc[t01:(t01+100)] <- siminc[t01:(t01+100)] + p1*mod1$inf_inc
  siminc[t02:(t02+100)] <- siminc[t02:(t02+100)] + p2*mod2$inf_inc
  notps <- length(siminc)
  nowaves <- length(waveBounds)-1
  startWaves <- waveBounds[1:nowaves]
  endWaves <- waveBounds[2:(nowaves+1)]
  cuminc <- vector(mode="numeric",length=length(siminc))
  cuminc[1] <- siminc[1]
  for (i in 2:notps) cuminc[i] <- cuminc[i-1]+siminc[i]

  waveinc <- array(0,dim=c(length(siminc),nowaves))

  for (i in 1:nowaves) {
    waveinc[startWaves[i]:endWaves[i],i] <- cuminc[startWaves[i]:endWaves[i]]-cuminc[startWaves[i]]
    waveinc[(endWaves[i]+1):length(siminc),i] <- cuminc[endWaves[i]] - cuminc[startWaves[i]]
  }

  rtn <- dumdat

  rtn <- cbind(rtn,ti=hks.sim.incidence.days(npeople,siminc,overallp=overallp))

  rtn <- cbind(rtn, 	inf.t12 = hks.test.foldrise(rtn$t1,rtn$t2,rtn$ti,td=t_rise),
               inf.t13 = hks.test.foldrise(rtn$t1,rtn$t3,rtn$ti,td=t_rise),
               inf.t14 = hks.test.foldrise(rtn$t1,rtn$t4,rtn$ti,td=t_rise))

  wavevars 	<- data.frame(index=1:npeople)
  nopairs 	<- length(pairdefs)

  for (pair in pairdefs) {
    tvar1 <- pair[1]
    tvar2 <- pair[2]
    for (i in 1:nowaves) {
      tmpvec <- vector(mode="numeric",length=npeople)
      maxwave <- max(waveinc[,i])
      for (j in 1:npeople) {
        if (rtn[j,tvar2]-t_rise < 1) stop("lkasdla")
        if (rtn[j,tvar1]-t_rise < 1) stop("lkasdla")
        propinc <- (waveinc[rtn[j,tvar2]-t_rise,i] - waveinc[rtn[j,tvar1]-t_rise,i]) / maxwave
        tmpvec[j] <- propinc
      }
      wavevars <- cbind(wavevars,tmpvec)
      nocols <- dim(wavevars)[2]
      names(wavevars)[nocols] <- paste(tvar1,".",tvar2,".w",i,sep="")
    }
  }

  # Up to here - this doesn't look right
  # plot(rtn$t3,wavevars$t1.t3.w1)

  wavevars$index = NULL

  rtn <- cbind(rtn,wavevars)

  # Output data in a similar form to our table
  list(simdat=rtn,vecinc=siminc)

}

# A simple likelihood model will probably work
hks.testlike <- function(prob,vecPropCAR,vecOutcome) {
  atomlnlike <- function(prob,car,out) {
    if (out) rtn <- log(car*prob)
    else rtn <- log(1-car*prob)
    rtn
  }
  value <- 0
  for (i in 1:length(vecPropCAR)) value <- value + atomlnlike(prob,vecPropCAR[i],vecOutcome[i])
  value
}

# Function to take the proportion of a single wave and add it to a second wave
hks.pe.cis <- function(y,x1) {

  ci_aux <- function(prob,vecPropCAR,vecOutcome,maxlnlike=0){hks.testlike(prob,vecPropCAR,vecOutcome)-maxlnlike+1.96}

  lb <- 0.00000001
  opt_pest <- optimize(
    f=hks.testlike,
    interval=c(lb,1.0-lb),
    vecPropCAR=x1,
    vecOutcome=y1,
    tol=1e-10,
    maximum=TRUE
  )

  pe 		<- opt_pest$maximum
  maxml 	<- opt_pest$objective

  lb 		<- uniroot(ci_aux,c(0,opt_pest$maximum),vecPropCAR=x1,vecOutcome=y1,maxlnlike=opt_pest$objective)
  ub 		<- uniroot(ci_aux,c(opt_pest$maximum,1),vecPropCAR=x1,vecOutcome=y1,maxlnlike=opt_pest$objective)

  list(pe=pe,lb=lb$root,ub=ub$root)

}

# Demonstrate that the obvioius logistic regression doesn't work
hks.logistic.fails <- function(z) {

  y1 <- z$simdat[,"inf.t12"]
  y2 <- z$simdat[,"inf.t14"]
  x1 <- z$simdat[,"t1.t2.w1"]
  x2 <- z$simdat[,"t1.t4.w1"]
  x3 <- 1-exp(-x1)
  x4 <- 1-exp(-x2)
  x5 <- z$simdat[,"t1.t4.w2"]

  mod1 <- glm(y1 ~ x1 +0, family=binomial(link="logit"))
  mod2 <- glm(y2 ~ x2 +0, family=binomial(link="logit"))
  mod3 <- glm(y1 ~ x1 +0, family=poisson)
  mod4 <- glm(y2 ~ x2 +0, family=poisson)

  print(summary(mod1))
  print(summary(mod2))
  print(summary(mod3))
  print(summary(mod4))

  # Lines below demonstrate that the linear model doen't work
  # The last line should be much higher because 54% are infected with
  # less than full exposure ??
  print(sum(y1)/length(y1))
  print(mean(x1))
  print(1/(1+exp(-mod1$coefficients[1])))
  print(sum(y2)/length(y2))
  print(mean(x2))
  print(1/(1+exp(-mod2$coefficients[1])))

  hist(x2)

}

# Demonstrate that the bernoilli likleihood seems to work
# Demonstrate that the benoilli likelihood works well
hks.bernoulli.works <- function() {
  car <- c(rep(0.5,100),rep(1,100))
  out <- c(rep(1,25),rep(0,75),rep(1,50),rep(0,50))
  delta <- 0.0001
  prob <- 0.500
  a <- hks.testlike(prob,car,out)
  b <- hks.testlike(prob-delta,car,out)
  c <- hks.testlike(prob+delta,car,out)
  print(c(a-b,a-c))
}

hks.contact.fig.age.specific <- function(df,
                                         outstem="/Users/sriley/Dropbox/shares/on_contact/sr_hk_contacts/plos_med_HK_contact/") {

  df$age.decile <- ifelse(df$age > 50,6,(df$age %/% 10 )  + 1)
  df$location.mid <- df$location.min+(df$location.max-df$location.min) / 2

  # Utilities for the chart
  xlabloc <- 1.5+(0:5)*2
  xlabs <- c("2-9","10-19","20-29","30-39","40-49","50+")

  # Make figure illustrate contacts
  pdf(file=paste(outstem,"Fig_4A_main_contacts_raw.pdf"))
  ymax=50
  boxplot(
    contact.time.gte.10.min~ffold*age.decile,
    data=df,
    col=c("green","red"),
    frame.plot=TRUE,
    axes=FALSE,
    xlab="Age decile",
    ylab="Number of contacts greater than or equal to 10 minutes",
    ylim=c(0,ymax))
  axis(2)
  axis(1,las=1,at=xlabloc,labels=xlabs)
  mtext("Figure 4, part A")
  legend(1.5, ymax, c("Not infected", "Infected"),fill = c("green", "red"))
  dev.off()

  # Make figure illustrate contacts
  pdf(file=paste(outstem,"Fig_4A_inset_contacts_raw.pdf"))
  boxplot(
    contact.time.gte.10.min~ffold*age.decile,
    data=df,
    col=c("green","red"),
    frame.plot=TRUE,
    axes=FALSE,
    xlab="",
    ylab="")
  axis(2)
  axis(1,las=1,at=xlabloc,labels=rep("",length(xlabloc)))
  mtext("Figure 4, part A, Inset")
  dev.off()

  # Make figure illustrate locations
  pdf(file=paste(outstem,"Fig_4B_locations_raw.pdf"))
  boxplot(location.mid~ffold*age.decile,data=df,col=c("green","red"),frame.plot=TRUE,axes=FALSE,xlab="Age decile",ylab="Median location")
  axis(2)
  axis(1,las=1,at=xlabloc,labels=xlabs)
  mtext("Figure 4, part B")
  dev.off()

}

hks.contact.fig.location.contact.compare <- function(df,
                                                     outstem="/Users/sriley/Dropbox/shares/on_contact/sr_hk_contacts/plos_med_HK_contact/") {

  age.decile <- ifelse(df$age > 70,8,(df$age %/% 10 )  + 1)
  location.mid <- df$location.min+(df$location.max-df$location.min) / 2
  contacts.mid <- df$contact.all.total.min+(df$contact.all.total.max-df$contact.all.total.min) / 2

  # Utilities for the chart
  xlabloc <- 1:8
  xlabs <- c("2-9","10-19","20-29","30-39","40-49","50-59","60-69","70+")

  # Make figure illustrate correlation
  pdf(file=paste(outstem,"Fig_2A_locations_raw.pdf"))
  boxplot(location.mid~age.decile,frame.plot=TRUE,axes=FALSE,xlab="Age deciles",ylab="Locations (mid of max to min)")
  axis(2)
  axis(1,las=1,at=xlabloc,labels=xlabs)
  mtext("Figure 2, part A")
  dev.off()

  # Make figure illustrate correlation
  pdf(file=paste(outstem,"Fig_2B_contacts_main_raw.pdf"))
  boxplot(
    contacts.mid~age.decile,
    frame.plot=TRUE,
    axes=FALSE,
    xlab="Age deciles",
    ylab="Contacts (mid of max to min)",
    ylim=c(0,60))
  axis(2)
  axis(1,las=1,at=xlabloc,labels=xlabs)
  mtext("Figure 2, part B, main")
  dev.off()

  # Make figure illustrate correlation
  pdf(file=paste(outstem,"Fig_2B_contacts_inset_raw.pdf"))
  boxplot(
    contacts.mid~age.decile,
    frame.plot=TRUE,
    axes=FALSE)
  axis(2)
  mtext("Figure 2, part B, inset")
  dev.off()

  # Make figure illustrate correlation
  pdf(file=paste(outstem,"Fig_2C_correlation_raw.pdf"))
  plot(
    location.mid,
    contacts.mid,
    log="y",
    xlab="Locations (mid of max to min)",
    ylab="Contacts (mid of max to min)",
    pch=20,
    cex=0.5)
  tmpxvals <- 1+13*(0:100)/100
  points(tmpxvals,tmpxvals,type="l")
  text(12,10,"y=x")
  mtext("Figure 2, Part C")
  dev.off()

}


# Calculate some summary statistics mentioned in the text
hks.contact.text.items <- function(df) {

  # Useful temporary variables
  noParts <- dim(df)[1]
  mincontacts <- sum(df$contact.all.total.min)
  maxcontacts <- sum(df$contact.all.total.max)
  midcontacts <- df$contact.all.total.min+(df$contact.all.total.max-df$contact.all.total.min)/2
  minlocations <- sum(df$location.min)
  maxlocations <- sum(df$location.max)
  midlocations <- df$location.min + (df$location.max-df$location.min)/2

  # Total numbers of contacts
  print(paste("Total minimum number of contacts was",mincontacts))
  print(paste("Total maximum number of contacts was",maxcontacts))
  print(paste("Range was",mincontacts/noParts,maxcontacts/noParts),digits=3)

  # Total numbers of locations
  print(paste("Total minimum number of locations was",minlocations))
  print(paste("Total maximum number of contacts was",maxlocations))
  print(paste("Range was",minlocations/noParts,maxlocations/noParts),digits=3)

  # Correlation between contacts and locations
  print(paste("cor(midcontacts,midlocations,method=\"pearson\"))",cor(midcontacts,midlocations,method="pearson")),digits=3)
  print(paste("cor(midcontacts-midlocations,midlocations),method=\"pearson\")",cor(midcontacts-midlocations,midlocations),method="pearson"),digits=3)
  print(paste("cor(log(midcontacts),midlocations,method=\"pearson\"))",cor(log(midcontacts),midlocations,method="pearson")),digits=3)
  print(paste("cor(log(midcontacts-midlocations-min(midcontacts-midlocations)+1),midlocations,method=\"pearson\"))",cor(log(midcontacts-midlocations-min(midcontacts-midlocations)+1),midlocations,method="pearson")),digits=3)

  print(paste("cor(midcontacts,midlocations,method=\"spearman\"))",cor(midcontacts,midlocations,method="spearman")),digits=3)
  print(paste("cor(midcontacts-midlocations,midlocations),method=\"spearman\")",cor(midcontacts-midlocations,midlocations),method="spearman"),digits=3)
  print(paste("cor(log(midcontacts),midlocations,method=\"spearman\"))",cor(log(midcontacts),midlocations,method="spearman")),digits=3)
  print(paste("cor(log(midcontacts-midlocations-min(midcontacts-midlocations)+1),midlocations,method=\"spearman\"))",cor(log(midcontacts-midlocations-min(midcontacts-midlocations)+1),midlocations,method="spearman")),digits=3)

  print(paste("min(midcontacts-midlocations)",min(midcontacts-midlocations)),digits=3)

  # Some related plots
  pdf(file="hks.contact.text.items_output.pdf")
  hist(midcontacts-midlocations,breaks=-4:500)
  hist(midcontacts-midlocations,breaks=-4:500,xlim=c(-4,20))
  dev.off()

}

hks.four.p.sim <- function(
  vecPs=c(0.1,0.2,0.3,0.4),
  npeeps=10,
  sizepop=10000,
  alpha=0.9,
  phi=0.1,
  sens=0.99,
  spec=0.99) {

  # Set up output variables
  rtnsero <- data.frame(
    y13=rep(-1,npeeps),y12=rep(-1,npeeps),y23=rep(-1,npeeps),
    y35=rep(-1,npeeps),y34=rep(-1,npeeps),y45=rep(-1,npeeps)
  )
  rtninf <- matrix(nrow=npeeps,ncol=4)
  rtnns <- vector(length=4,mode="numeric")

  # Generate the ns
  rtnns[1] <- rbinom(1,sizepop,phi*vecPs[1])
  rtnns[2] <- rbinom(1,sizepop-rtnns[1],phi*vecPs[2])
  tmpdoubles <- rbinom(1,rtnns[1]+rtnns[2],phi*alpha*vecPs[3])
  rtnns[3] <- rbinom(1,sizepop-rtnns[1]-rtnns[2],phi*vecPs[3]) + tmpdoubles
  rtnns[4] <- rbinom(1,sizepop-rtnns[1]-rtnns[2]-rtnns[3],vecPs[4]) + rbinom(1,rtnns[1]+rtnns[2]-tmpdoubles,phi*alpha*vecPs[3])

  # Simulate the first period
  rtninf[,1] <- rbinom(npeeps,1,vecPs[1])

  # Simulate the second period
  rtninf[,2] <- rbinom(npeeps,1-rtninf[,1],vecPs[2])

  # Simulate the third period
  tmpC <- rbinom(npeeps,rowSums(rtninf[,1:2]),alpha*vecPs[3])
  rtninf[,3] <- rbinom(npeeps,1-rowSums(rtninf[,1:2]),vecPs[3]) + tmpC

  # Simulate the fourth period
  rtninf[,4] <- rbinom(npeeps,1-apply(rtninf[,1:3],1,max),vecPs[4]) + rbinom(npeeps,rowSums(rtninf[,1:2])-tmpC,alpha*vecPs[4])

  # Set up start periods and end periods for assay results
  starts <- c(1,1,2,3,3,4)
  ends <- c(2,1,2,4,3,4)
  noperiods <- length(starts)

  # Generate assay results
  # This must be consistent with rtnseronaming
  # Consider automating the naming of rtnsero
  for (i in 1:noperiods) {
    rtnsero[,i] <- rbinom(npeeps,rowSums(cbind(rep(0,npeeps),rtninf[,starts[i]:ends[i]])),spec) +
      rbinom(npeeps,1-rowSums(cbind(rep(0,npeeps),rtninf[,starts[i]:ends[i]])),1-sens)
  }

  # Return output variables
  list(sero=rtnsero,inf=rtninf,ns=rtnns)

}

# The double likelihood looks right, but the quadruple might not
# Needs switching to log and quadruple testing with more
hks.four.p.like <- function(ytab,p1=0.1,p2=0.2,p3=0.3,p4=0.4,alpha=0.9) {

  nopeeps <- dim(ytab)[1]

  y13 <- ytab[,"y13"]
  y12 <- ytab[,"y12"]
  y23 <- ytab[,"y23"]
  y35 <- ytab[,"y35"]
  y34 <- ytab[,"y34"]
  y45 <- ytab[,"y45"]

  rtndouble <- (p1+p2)^y13 * (1-p1-p2)^(1-y13) * (p3+p4)^(y35*(1-y13)) * (1-p3-p4)^((1-y35)*(1-y13))
  rtndouble <- rtndouble * (alpha*(p3+p4))^(y35*y13) * (1-alpha*(p3+p4))^((1-y35)*y13)

  rtnquad <- p1^y12 * (1-p1)^(1-y12) * p2^(y23*(1-y12)) * (1-p2)^((1-y23)*(1-y12)) * p3^(y34*(1-y12-y23))
  rtnquad <- rtnquad * (1-p3)^((1-y34)*(1-y12-y23)) * (alpha*p3)^(y34*(y12+y23)) * (1-alpha*p3)^((1-y34)*(y12+y23))
  rtnquad <- rtnquad * p4^(y45*(1-y34)*(1-y12-y23)) * (1-p4)^((1-y45)*(1-y34)*(1-y12-y23))
  rtnquad <- rtnquad * (alpha*p4)^(y45*(1-y34)*(y12+y23)) * (1-alpha*p4)^((1-y45)*(1-y34)*(y12+y23))

  list(pPair=sum(log(rtndouble)),pQuad=sum(log(rtnquad)))

}

# Make a survival matrix from the data
hks.aux.survive <- function(vecInc,phi=0.001,N=7e7,waveBounds = c(250,400)) {
  nwaves <- length(waveBounds)
  lengthvec <- length(vecInc)
  nodays <- waveBounds[nwaves]-1
  if (nodays > lengthvec) stop("laksdjalkjd")
  rtn <- matrix(nrow=nodays,ncol=nwaves)
  waveTots <- vector(length=nwaves,mode="numeric")

  # Record wave totals
  completeWaveBounds <- c(1,waveBounds)
  for (i in 1:(nwaves)) {
    waveTots[i] <- sum(vecInc[completeWaveBounds[i]:completeWaveBounds[i+1]])
  }

  curday <- 1
  curwave <- 1
  startnextwave <- waveBounds[1]
  while (curday <= nodays) {
    if (curday == waveBounds[curwave]) curwave <- curwave + 1
    for (i in 1:nwaves) {
      if (curday < completeWaveBounds[i]) rtn[curday,i] <- 0
      else if (curday >= completeWaveBounds[i] & curday < completeWaveBounds[i+1]) {
        if (waveTots[i] < 1e-10) rtn[curday,i] <- 0
        else rtn[curday,i] <- vecInc[curday]/waveTots[i]
      }
      else rtn[curday,i] <- 0
    }
    curday <- curday + 1
  }
  for (i in 2:nodays) rtn[i,] <- rtn[i-1,] + rtn[i,]
  rtn
}

# Function to load and cleant he QMH data
hks.read.clean.QMH <- function(filename) {
  require("date")
  rtn <- read.csv(filename)
  rtn[,"week_ending"] <- as.date(as.character(rtn[,"week_ending"]),order="dmy")
  rtn
}

# Function to load and cleant he QMH data
hks.read.clean.chp <- function(filename) {
  require("date")
  rtn <- read.csv(filename)
  rtn[,"month_starting"] <- as.date(as.character(rtn[,"month_starting"]),order="dmy")
  rtn[,"month_ending"] <- as.date(as.character(rtn[,"month_ending"]),order="dmy")
  rtn$TotalA <- rowSums(rtn[,c("sH3N2","sH1N1","pH1N1","H5N1","H9N2")])
  rtn
}

# Function to load and cleant gp weekly data with strains
hks.read.clean.gpw <- function(filename) {
  require("date")
  rtn <- read.csv(filename)

  rtn[,"date.wk.start"] <- as.date(as.character(rtn[,"date.wk.start"]),order="dmy")
  rtn[,"date.wk.end"] <- as.date(as.character(rtn[,"date.wk.end"]),order="dmy")

  rtn$gp.h1n1 <- rtn$A.H1N1pdm / rtn$n.samps * rtn$GP.per.1000
  rtn$gp.h3n2 <- rtn$A.H3N2 / rtn$n.samps * rtn$GP.per.1000
  rtn$gp.b <- rtn$B / rtn$n.samps * rtn$GP.per.1000
  rtn$gp.A <- rtn$gp.h1n1 + rtn$gp.h3n2

  rtn
}

# Read and clean the study data from the multiple waves
hks.read.clean.study <- function(filename) {

  # Load up the raw data
  rtn <- read.csv(filename, stringsAsFactors=FALSE)

  # Define all the quads
  rtn$is.quad <- rtn$blood_1=="Yes" & rtn$blood_2=="Yes" & rtn$blood_3=="Yes" & (rtn$blood_4 %in% c("Yes","yes","YES"))

  # Correct the dates of the bloods
  # These now handled by the lines on wrote below

# rtn$blood_1_date<- as.date(paste(substr(rtn$blood_1_date,start=1,stop=4),'/',substr(rtn$blood_1_date,start=6,stop=7),'/',substr(rtn$blood_1_date,start=9,stop=10)),order="ymd")
# rtn$blood_2_date<- as.date(paste(substr(rtn$blood_2_date,start=1,stop=4),'/',substr(rtn$blood_2_date,start=6,stop=7),'/',substr(rtn$blood_2_date,start=9,stop=10)),order="ymd")
# rtn$blood_3_date<- as.date(paste(substr(rtn$blood_3_date,start=1,stop=4),'/',substr(rtn$blood_3_date,start=6,stop=7),'/',substr(rtn$blood_3_date,start=9,stop=10)),order="ymd")
# rtn$blood_4_date<- as.date(paste(substr(rtn$blood_4_date,start=1,stop=4),'/',substr(rtn$blood_4_date,start=6,stop=7),'/',substr(rtn$blood_4_date,start=9,stop=10)),order="ymd")



     #rtn$blood_1_date <-Sys.Date()
     #rtn$blood_2_date <-Sys.Date()
     #rtn$blood_3_date <-Sys.Date()
     #rtn$blood_4_date <-Sys.Date()

#rtn$blood_1_date<- format( rtn$blood_1_date, format="%d-%b-%Y")
#rtn$blood_2_date<-format( rtn$blood_2_date, format="%d-%b-%Y")
#rtn$blood_3_date<-format( rtn$blood_3_date, format="%d-%b-%Y")
#rtn$blood_4_date<-format( rtn$blood_4_date, format="%d-%b-%Y")

  #rtn$blood_1_date <- as.date(as.character(rtn$blood_1_date),order="dmy")
  #rtn$blood_2_date <- as.date(as.character(rtn$blood_2_date),order="dmy")
  #rtn$blood_3_date <- as.date(as.character(rtn$blood_3_date),order="dmy")
  #rtn$blood_4_date <- as.date(as.character(rtn$blood_4_date),order="dmy")

rtn$blood_1_date <- as.Date(rtn$blood_1_date)
rtn$blood_2_date <- as.Date(rtn$blood_2_date)
rtn$blood_3_date <- as.Date(rtn$blood_3_date)
rtn$blood_4_date <- as.Date(rtn$blood_4_date)

  # hist(c(rtn$blood_1_date,rtn$blood_2_date,rtn$blood_3_date,rtn$blood_4_date),breaks=17800+(0:40)*35)
  rtn$pT1 <- log((rtn$H1N1.T1 / 5),2)
  rtn$pT2 <- log((rtn$H1N1.T2 / 5),2)
  rtn$pT3 <- log((rtn$H1N1.T3 / 5),2)
  rtn$pT4 <- log((rtn$H1N1.T4 / 5),2)

  rtn$sT1 <- log((rtn$H3N2.T1 / 5),2)
  rtn$sT2 <- log((rtn$H3N2.T2 / 5),2)
  rtn$sT3 <- log((rtn$H3N2.T3 / 5),2)
  rtn$sT4 <- log((rtn$H3N2.T4 / 5),2)

  # Add soem age groups
  rtn$ag_dec <- cut(rtn$Age_1,breaks=0:8*10,labels=FALSE)
  rtn$ag1 <- cut(rtn$Age_1,breaks=c(0,19,45,65,110),labels=FALSE)
  rtn$ag2 <- cut(rtn$Age_1,breaks=c(0,22,48.1,56.7,110),labels=FALSE)

  rtn$ff_r12_ph1n1 <- rtn$pT2 - rtn$pT1 > 1.9 & rtn$H1N1.T2 > 20
  rtn$ff_r13_ph1n1 <- rtn$pT3 - rtn$pT1 > 1.9 & rtn$H1N1.T3 > 20
  rtn$ff_r14_ph1n1 <- rtn$pT4 - rtn$pT1 > 1.9 & rtn$H1N1.T4 > 20
  rtn$ff_r23_ph1n1 <- rtn$pT3 - rtn$pT2 > 1.9 & rtn$H1N1.T3 > 20
  rtn$ff_r24_ph1n1 <- rtn$pT4 - rtn$pT2 > 1.9 & rtn$H1N1.T4 > 20
  rtn$ff_r34_ph1n1 <- rtn$pT4 - rtn$pT3 > 1.9 & rtn$H1N1.T4 > 20

  rtn$ff_r12_sh3n2 <- rtn$sT2 - rtn$sT1 > 1.9 & rtn$H3N2.T2 > 20
  rtn$ff_r13_sh3n2 <- rtn$sT3 - rtn$sT1 > 1.9 & rtn$H3N2.T3 > 20
  rtn$ff_r14_sh3n2 <- rtn$sT4 - rtn$sT1 > 1.9 & rtn$H3N2.T4 > 20
  rtn$ff_r23_sh3n2 <- rtn$sT3 - rtn$sT2 > 1.9 & rtn$H3N2.T3 > 20
  rtn$ff_r24_sh3n2 <- rtn$sT4 - rtn$sT2 > 1.9 & rtn$H3N2.T4 > 20
  rtn$ff_r34_sh3n2 <- rtn$sT4 - rtn$sT3 > 1.9 & rtn$H3N2.T4 > 20

  rtn$ff_any_pH1N1 <- rtn$ff_r12_ph1n1 |
    rtn$ff_r13_ph1n1 |
    rtn$ff_r14_ph1n1 |
    rtn$ff_r23_ph1n1 |
    rtn$ff_r24_ph1n1 |
    rtn$ff_r34_ph1n1

  rtn$ff_any_sH3N2 <- rtn$ff_r12_sh3n2 |
    rtn$ff_r13_sh3n2 |
    rtn$ff_r14_sh3n2 |
    rtn$ff_r23_sh3n2 |
    rtn$ff_r24_sh3n2 |
    rtn$ff_r34_sh3n2

  rtn$ff_r12_both <- rtn$ff_r12_ph1n1 & rtn$ff_r12_sh3n2
  rtn$ff_r13_both <- rtn$ff_r13_ph1n1 & rtn$ff_r13_sh3n2
  rtn$ff_r14_both <- rtn$ff_r14_ph1n1 & rtn$ff_r14_sh3n2
  rtn$ff_r23_both <- rtn$ff_r23_ph1n1 & rtn$ff_r23_sh3n2
  rtn$ff_r24_both <- rtn$ff_r24_ph1n1 & rtn$ff_r24_sh3n2
  rtn$ff_r34_both <- rtn$ff_r34_ph1n1 & rtn$ff_r34_sh3n2

  # On's correction for the blood 4 data entry
  rtn$blood_4[rtn$blood_4=="no"] <- "No"
  rtn$blood_4[rtn$blood_4=="NO"] <- "No"
  rtn$blood_4[rtn$blood_4=="no "] <- "No"
  rtn$blood_4[rtn$blood_4=="No "] <- "No"
  rtn$blood_4[rtn$blood_4=="YES"] <- "Yes"
  rtn$blood_4[rtn$blood_4=="yes" ] <- "Yes"

  # Steven's number of blood samples field
  blood_mat <- as.matrix(cbind(
    rtn$blood_1=="Yes",rtn$blood_2=="Yes",rtn$blood_3=="Yes",rtn$blood_4=="Yes"))
  rtn$blood_count <- apply(blood_mat,1,sum,na.rm=TRUE)

  # Need to break up the effects below
  rtn$nomissing <- !is.na(rtn$ff_r23_sh3n2) &
    !is.na(rtn$blood_3_date) & !is.na(rtn$Age_1)

  rtn$qDat <- rtn$is.quad==TRUE &
    !is.na(rtn$is.quad) & rtn$nomissing

rtn$vac0910_h1n1_cleaned <- rep(0,dim(rtn)[1])

for (i in 1: dim(rtn)[1]) {
  if (rtn$vac0910_h1n1_cleaned.r3[i]=="Yes" | rtn$mono.r2[i]=="Yes")
   rtn$vac0910_h1n1_cleaned[i]="Yes"

  if (rtn$vac0910_h1n1_cleaned.r3[i]=="No" &  rtn$mono.r2[i]=="No")
   rtn$vac0910_h1n1_cleaned[i]="No"

  if (rtn$vac0910_h1n1_cleaned.r3[i]=="Missing" & rtn$mono.r2[i]=="No")
   rtn$vac0910_h1n1_cleaned[i]="Missing"

   if (rtn$vac0910_h1n1_cleaned.r3[i]=="No" & rtn$mono.r2[i]=="Missing")
    rtn$vac0910_h1n1_cleaned[i]="No"

  if (rtn$vac0910_h1n1_cleaned.r3[i]=="Missing" & rtn$mono.r2[i]=="Missing")
    rtn$vac0910_h1n1_cleaned[i]="Missing"

}

rtn$vac0910_cleaned <- rep(0,dim(rtn)[1])

for (i in 1: dim(rtn)[1]) {
  if (rtn$vac0910_cleaned.r3[i]=="Yes" | rtn$seasonal.r2[i]=="Yes")
   rtn$vac0910_cleaned[i]="Yes"

  if (rtn$vac0910_cleaned.r3[i]=="No"  & rtn$seasonal.r2[i]=="No")
   rtn$vac0910_cleaned[i]="No"

  if (rtn$vac0910_cleaned.r3[i]=="Missing" & rtn$seasonal.r2[i]=="No")
    rtn$vac0910_cleaned[i]="Missing"

  if (rtn$vac0910_cleaned.r3[i]=="No" & rtn$seasonal.r2[i]=="Missing")
    rtn$vac0910_cleaned[i]="No"

  if (rtn$vac0910_cleaned.r3[i]=="Missing" & rtn$seasonal.r2[i]=="Missing")
    rtn$vac0910_cleaned[i]="Missing"

}

rtn$vac1011_cleaned <- rep(0,dim(rtn)[1])

for (i in 1: dim(rtn)[1]) {
  if (rtn$vac1011.r3[i]=="Yes" | rtn$vac1011.r4[i]=="Yes")
   rtn$vac1011_cleaned[i]="Yes"

  if (rtn$vac1011.r3[i]=="No" & rtn$vac1011.r4[i]=="No")
   rtn$vac1011_cleaned[i]="No"

  if (rtn$vac1011.r3[i]=="Missing" & rtn$vac1011.r4[i]=="No")
    rtn$vac1011_cleaned[i]="No"

  if (rtn$vac1011.r3[i]=="No" & rtn$vac1011.r4[i]=="Missing")
    rtn$vac1011_cleaned[i]="Missing"

  if (rtn$vac1011.r3[i]=="Missing" & rtn$vac1011.r4[i]=="Missing")
    rtn$vac1011_cleaned[i]="Missing"
}

  rtn$novac <- rtn$vac1011_cleaned=="No" & rtn$vac0910_h1n1_cleaned=="No" & rtn$vac0910_cleaned=="No"

  rtn$vac <- rtn$vac1011_cleaned=="Yes" | rtn$vac0910_h1n1_cleaned=="Yes" | rtn$vac0910_cleaned=="Yes"

  # Make some other vaccine -based subgroups
  rtn$vac0910Only <- rtn$vac0910_h1n1_cleaned=="No" & rtn$vac0910_cleaned=="Yes" & !(rtn$novac)
  rtn$vac0910PanOnly <- rtn$vac0910_h1n1_cleaned=="Yes" & rtn$vac0910_cleaned=="No" & !(rtn$novac)

  # rtn$blood_4 <- droplevels(rtn$blood_4)
  rtn

}

hks.plot.waves.aux <- function(
  auxD,
  stuD,
  filename,
  cwidth=10,
  cheight=15) {

  # Breaks to use
   #breaks <- 17800+0:200*7
  breaks <- 14147+0:200*7 # as.numeric(as.Date("2008-09-25"))
  #date_labels <- c("01-Jan-2009","01-Jul-2009","01-Jan-2010","01-Jul-2010","01-Jan-2011","01-Jul-2011","31-Dec-2011")
  date_labels <- c("2009-01-01","2009-07-01","2010-01-01","2010-07-01","2011-01-01","2011-07-01","2011-12-31")
  #date_vals <- as.date(date_labels)
  date_vals <- as.Date(date_labels)
  #date_range <- c(min(date_vals),max(date_vals))
  date_range <- as.numeric(c(min(date_vals),max(date_vals)))

  # Cycle through different $blood_4_date to see the recorded wave dates
  pdf(filename, width=cwidth/cm(1), height=cheight/cm(1))
  par(mai=c(0,0,0,0))
  par(fig=c(0.15,0.85,0.25,0.95))
  # mgp and pos here
  #hist(	c(stuD$blood_1_date,stuD$blood_2_date,stuD$blood_3_date,stuD$blood_4_date),
  hist(  as.numeric(c(stuD$blood_1_date,stuD$blood_2_date,stuD$blood_3_date,stuD$blood_4_date)),
        breaks=breaks,lwd=1,
        main="",
        xlab="Date (17800 = 25-Sep-2008) ",
        xlim=date_range,
        axes=FALSE)
  axis(1,at=date_vals,labels=date_labels,las=3)
  axis(2)
  #hist(	c(stuD$blood_1_date),breaks=breaks,add=TRUE,col="grey",border="grey",lwd=30)
  #hist(	c(stuD$blood_2_date),breaks=breaks,add=TRUE,col="grey",border="grey",lwd=30)
  #hist(	c(stuD$blood_3_date),breaks=breaks,add=TRUE,col="grey",border="grey",lwd=30)
  #hist(	c(stuD$blood_4_date),breaks=breaks,add=TRUE,border="grey",col="grey",lwd=30)
  hist( as.numeric(c(stuD$blood_1_date)),breaks=breaks,add=TRUE,col="grey",border="grey",lwd=30)
  hist(	as.numeric(c(stuD$blood_2_date)),breaks=breaks,add=TRUE,col="grey",border="grey",lwd=30)
  hist(	as.numeric(c(stuD$blood_3_date)),breaks=breaks,add=TRUE,col="grey",border="grey",lwd=30)
  hist(	as.numeric(c(stuD$blood_4_date)),breaks=breaks,add=TRUE,border="grey",col="grey",lwd=30)

  par(new=TRUE)
  auxD$month_starting <- auxD$month_starting - 3653 # V added
  auxD$month_ending <- auxD$month_ending - 3653  # V added

  plot(auxD$month_starting+(auxD$month_ending - auxD$month_starting)/2,
       auxD$TotalA,
       type="l",
       col="black",
       lwd=3,
       xlim=date_range,
       axes = FALSE)
  axis(4)
  points(auxD$month_starting+(auxD$month_ending - auxD$month_starting)/2,auxD$sH3N2,type="l",col="green",lwd=2)
  points(auxD$month_starting+(auxD$month_ending - auxD$month_starting)/2,auxD$pH1N1,type="l",col="red",lwd=2)

  #abline(v=as.date("10-Apr-2009"))
  #abline(v=as.date("15-May-2010"))
  #abline(v=as.date("01-Dec-2010"))
  abline(v=as.numeric(as.Date("2009-04-10")))
  abline(v=as.numeric(as.Date("2010-05-15")))
  abline(v=as.numeric(as.Date("2010-12-01")))


  dev.off()

}

# A simple likelihood model for two waves of infection
# Figure out how to vectorize this calculation to improve speed
hks.testlike.2w <- function(prob,vecPropCAR1,vecPropCAR2,vecOutcome) {
  atomlnlike <- function(vecProb,car1,car2,out) {
    if (out) rtn <- log(car1*prob[1](1-car2*prob[2])+car2*prob[2](1-car1*prob[1]))
    else rtn <- log((1-car1*prob[1])*(1-car2*prob[2]))
    rtn
  }
  value <- 0
  for (i in 1:length(vecPropCAR)) value <- value + atomlnlike(prob,vecPropCAR1[i],vecPropCAR2[i],vecOutcome[i])
  value
}

hks.ill.sim.2.wave <- function(filename,sim) {

  # Simulate from a similar study design
  # Probably not need for a bit


  # Make basis plot for Supporting Figure
  pdf(file=filename,width=8.7 / cm(1), height=7 / cm(1))
  par(mgp=c(1.5,0.5,0))
  par(mai=(c(0.55,0.55,0.1,0.1)))
  plot(sim$vecinc,type="l",xlab="Time",ylab="Daily incidence",xlim=c(0,400))
  dev.off()

}

hks.plot.long.chp <- function(file,x) {
  pdf(file)
  offset <- 0
  plot(x$month_starting+15,x$TotalA+offset,type="l",col="grey",lwd=4,axes=FALSE,log="")
  axis(1)
  axis(2)
  points(x$month_starting+15,x$sH3N2+offset,type="l",col="green",lwd=2)
  points(x$month_starting+15,x$pH1N1+offset,type="l",col="red",lwd=2)
  points(x$month_starting+15,x$sH1N1+offset,type="l",col="blue",lwd=2)
  dev.off()
}

hks.plot.long.gp <- function(file,x) {
  pdf(file)
  offset <- 0
  plot(x$date.wk.start+3.5,x$gp.A+offset,type="l",col="grey",lwd=6,axes=FALSE,log="")
  axis(1)
  axis(2)
  points(x$date.wk.start+3.5,x$gp.h3n2+offset,type="l",col="green",lwd=2)
  points(x$date.wk.start+3.5,x$gp.h1n1+offset,type="l",col="red",lwd=2)
  dev.off()
}

# Make some confidence bounds for the 4th round
hks.make.cis <- function(vecPos) {
  rtn_pe <- vecPos / sum(vecPos)
  rtn_lb <- srg.ci.binom(rep(sum(vecPos),length(vecPos)),vecPos,P=0.975)
  rtn_ub <- srg.ci.binom(rep(sum(vecPos),length(vecPos)),vecPos,P=0.025)
  return(data.frame(pe=rtn_pe,lb=rtn_lb,ub=rtn_ub))
}

# Print some details needed for the munich meeting
hks.details.ab.munich <- function(x) {
  print(paste("Number of quads ",dim(x)[1]))
  print(paste("Round 1",range(x[,"blood_1_date"])))
  print(paste("Round 2",range(x[,"blood_2_date"])))
  print(paste("Round 3",range(x[,"blood_3_date"],na.rm=TRUE)))
  print(paste("Round 4",range(x[,"blood_4_date"])))
  print(paste("Age",range(x$Age_1,na.rm=TRUE)))
}

# Develop and test the proportion of wave function
# Assumes t1 < t2 and t1 and t2 not in the same dt
hks.prop.wave <- function(vecDtStart,vecCumInc,t1,t2,lag=0) {
  if (length(vecCumInc) != length(vecDtStart)) stop("Problem in hks.plot.waves.aux")
  maxindex <- length(vecDtStart)
  startIndex <- 1
  while (vecDtStart[startIndex] < (t1 - lag) & startIndex <= maxindex) startIndex <- startIndex + 1
  endIndex <- startIndex
  while (vecDtStart[endIndex] < (t2 - lag) & endIndex <= maxindex) endIndex <- endIndex + 1
  (vecCumInc[endIndex] - vecCumInc[startIndex]) / max(vecCumInc)
}

# Take a continuous time series of incidence proxy, dates for wave starts and finishes and create a
# cumulative timeseries for the

hks.making.waves <- function(
  waveStart,
  waveEnd,
  vecIncidenceStarts,
  vecIncidenceEnds,
  vecIncidence) {

  # Check arguments
  notps <- length(vecIncidence)
  if (notps != length(vecIncidenceStarts)) stop("Problem in hks.making.waves")
  rtn <- vector(mode="numeric",length=length(vecIncidenceStarts))

  # Make the timeseries
  timeMask <- (vecIncidenceStarts) <= waveEnd & (vecIncidenceEnds >= waveStart)

  rtn <- vecIncidence * timeMask
  for (i in 2:notps) rtn[i] <- rtn[i-1] + rtn[i]

  # Return vector of cumulative incidence
  rtn
}

# Bernoilli likelihood for a period covering only a single wave
hks.sing.bern <- function(a,vecExp,vecOutcome) {
  prod((a*vecExp)^vecOutcome*(1-a*vecExp)^(1-vecOutcome))
}

# Define a function to calculate a single likelihood of a single obseration
hks.time.func.1 <- function(t,beta_y0,beta_ti,beta_b,beta_w,t_r=21,t_year = 356.25) {
  rtn <- -1
  t_peak <- beta_ti + t_r
  if (t < beta_ti) {
    rtn <- beta_y0
  } else if (t < t_peak) {
    rtn <- beta_y0 + beta_b - (t_peak - t) / t_r * beta_b
  } else {
    rtn <- beta_y0 + beta_b - (t - t_peak) / t_year * beta_w
  }
  if (rtn < 0) rtn <- 0
  rtn
}

# Define a likelihood for a single known individual
hks.time.lnlike.1.1 <- function(y,x) {
  shp <- 10
  rtn <- srg.titre.dens(y,x)
  # rtn <- srg.mod.gamma.dens(y,x,shp=shp)
  # rtn <- dgamma(y+1,shape=shp,scale=(x+1)/shp,log=TRUE) - log(1-dgamma(0,shape=shp,scale=(x+1)/shp,log=FALSE))
  # rtn <- log(pgamma(y+2,shape=shp,scale=(x+1)/shp,log=FALSE)-pgamma(y+1,shape=shp,scale=(x+1)/shp,log=FALSE)) - log(1-pgamma(1,shape=shp,scale=(x+1)/shp,log=FALSE))
  if (rtn < -1e100) stop("alksdjalksd problem here")
  rtn
}

# Calculate a single likelihood for participants
hks.lnlike.full.1.1 <- function(mT,mS,vY0,vT,vB,vW) {
  dimMat <- dim(mT)
  nQ <- dimMat[1]
  nB <- dimMat[2]
  lnlike <- vector(mode="numeric",length=nQ)
  for (i in 1:nQ) {
    for (j in 1:nB) {
      t <- mT[i,j]
      x <- hks.time.func.1(t,vY0[i],vT[i],vB[i],vW[i])
      y <- mS[i,j]
      lnlike[i] <- lnlike[i] + hks.time.lnlike.1.1(y,x)
    }
  }
  lnlike
}

# Simple timestep update
hks.update.ts <- function(vT,mint,maxt,dtmax=7,timeprior=NULL) {
  nQ <- length(vT)
  rtn <- vector(length=nQ,mode="numeric")
  for (i in 1:nQ) {
    dt <- dtmax*2*runif(1)-dtmax
    propt <- vT[i] + dt
    if (propt >= mint & propt <= maxt) {
      rtn[i] <- propt
    } else if (propt < mint) {
      rtn[i] <- 2*mint - propt
    } else {
      rtn[i] <- 2*maxt - propt
    }
  }
  rtn
}

# Time step update using an informative prior
# Note that time prior has to be at least dt larger than [mint,maxt]
# Args:
#  vt - vector of current times
#  mintt - minimum time
#  maxt - maximum time
#  propt - proposal time
hks.update.ts.inform <- function(vT,mint,maxt,timeprior,dtmax=7) {

  # Setup some auxiliiary functions
  nQ <- length(vT)
  rtn <- vector(length=nQ,mode="numeric")
  propratio <- vector(length=nQ,mode="numeric")
  tpt0 <- timeprior$time[1]
  tpdt <- timeprior$time[2]-tpt0

  # Loop through each individual
  for (i in 1:nQ) {

    # Below needs modifying for edge effects
    # Should read them exactly from the prior distribution
    curT <- vT[i]
    mincurdens <- hks.cuminc.density(curT-dtmax,timeprior$time,timeprior$cuminc)
    maxcurdens <- hks.cuminc.density(curT+dtmax,timeprior$time,timeprior$cuminc)
    maxdenscurchange <- (maxcurdens-mincurdens)
    denscurchange <- runif(1)*maxdenscurchange

    propt <- hks.cuminc.time(curT-dtmax,denscurchange,timeprior$time,timeprior$cuminc)

    minpropdens <- hks.cuminc.density(propt-dtmax,timeprior$time,timeprior$cuminc)
    maxpropdens <- hks.cuminc.density(propt+dtmax,timeprior$time,timeprior$cuminc)
    maxdenspropchange <- (maxpropdens-minpropdens)

    # This line is wrong - should be the density lookup for the proposed time
    denspropchange <- hks.cuminc.density(propt,timeprior$time,timeprior$cuminc)

    # problem is here
    # edge effects have to be wrong
    # needs normalizing properly to reflect accurately the probabilty of the forward move
    # relative to the reverse move

    propratio[i] <- log(timeprior$inc[as.integer((curT-tpt0))/tpdt] /
                          timeprior$inc[as.integer((propt-tpt0))/tpdt] *
                          maxdenscurchange /
                          maxdenspropchange)

    # Below here is the same as the uniform version
    if (propt >= mint & propt <= maxt) {
      rtn[i] <- propt
    } else if (propt < mint) {
      rtn[i] <- 2*mint - propt
    } else {
      rtn[i] <- 2*maxt - propt
    }
  }
  list(vt=rtn,pr=propratio)
}

# Prep the data for the MCMC sampler
hks.mcmc.setup.1.1 <- function(datWide,min_bd,max_bd,fnPs,t_growth=21,maxnobloods=4) {

  maskAnyH1Ffold <- datWide$ff_w12_ph1n1 | datWide$ff_w34_ph1n1 | datWide$ff_w24_ph1n1
  fDat <- datWide[maskAnyH1Ffold,]

  noFfold <- dim(fDat)[1]

  # Load in the parameter values
  matP <- read.csv(fnPs,row.names=1)

  # Define vector of infection times and parameters
  vecTs 		<- vector(length=noFfold,mode="numeric")
  vecY0		<- vector(length=noFfold,mode="numeric")

  # Set aux values for mcmc chain
  min_inf_time <- min_bd - t_growth
  max_inf_time <- max_bd

  # Set up initial values for the parameter vectors
  vecTs[]		<- min_inf_time + (max_inf_time - min_inf_time)/2
  vecY0[]		<- 0

  chain <- list(vt=vecTs,y0=vecY0,ps=matP)

  # Make two data mtrices of times and integer values
  matTime <- matrix(nrow=noFfold,ncol=maxnobloods)
  matH1 <- matrix(nrow=noFfold,ncol=maxnobloods)
  for (i in 1:maxnobloods) matTime[,i] <- fDat[,paste("blood_",i,"_date",sep="")]
  for (i in 1:maxnobloods) matH1[,i] <- log(fDat[,paste("H1N1.T",i,sep="")]/5,2)

  # Return chain and data
  list(mt=matTime,mh=matH1,chain=chain)

}

# XX Rewrite here for any strain
hks.extract.ffold.risers <- function(dat,agerange=c(0,110),strain="ph1n1") {

  maskAnyFfold 	<- dat[,paste("ff_r12_",strain,sep="")] | dat[,paste("ff_r24_",strain,sep="")] | dat[,paste("ff_r34_",strain,sep="")]
  maskAge 		<- dat[,"Age_1"] > agerange[1] & dat[,"Age_1"] < agerange[2]
  fDat <- dat[maskAnyFfold & maskAge,]
  fDat

}

# Prep the data for the MCMC sampler
# x <- runif(10)
# x[order(x)]
hks.mcmc.setup.1.2 <- function(
  datWide,
  min_bd,
  max_bd,
  fnPs,
  strain="ph1n1",
  t_growth=21,
  maxnobloods=4,
  agerange=c(0,110)) {

  fDat <- hks.extract.ffold.risers(datWide,agerange,strain)

  noFfold <- dim(fDat)[1]

  # Load in the parameter values
  matP <- read.csv(fnPs,row.names=1)

  # Define vector of infection times and parameters
  vecTs 		<- vector(length=noFfold,mode="numeric")
  vecY0		<- vector(length=noFfold,mode="numeric")

  # Set aux values for mcmc chain
  min_inf_time <- min_bd - t_growth
  max_inf_time <- max_bd

  # Set up initial values for the parameter vectors
  # Arbitrary assumption about starting value here
  vecTs[]		<- min_inf_time + (max_inf_time - min_inf_time)/2
  vecY0[]		<- log(fDat$H1N1.T1 / 5,2) + 0.5

  chain <- list(vt=vecTs,y0=vecY0,ps=matP)

  # Make two data mtrices of times and integer values
  matTime <- matrix(nrow=noFfold,ncol=maxnobloods)
  matH1 <- matrix(nrow=noFfold,ncol=maxnobloods)
  for (i in 1:maxnobloods) matTime[,i] <- fDat[,paste("blood_",i,"_date",sep="")]
  for (i in 1:maxnobloods) matH1[,i] <- log(fDat[,paste("H1N1.T",i,sep="")]/5,2)

  # Return chain and data
  list(mt=matTime,mh=matH1,chain=chain)

}

# Main engine for the mcmc sampler
hks.mcmc.1.1 <- function(
  dataAndChain,
  min_inf_time,
  max_inf_time,
  runstem = "test_",
  outputtop = "./",
  batchsize = 100,
  thin_factor = 1,
  nosamples = 100,
  updatebalance = c(1,1,0),
  dtmax=20,
  debugacceptall=FALSE,
  filemustbeabsent=TRUE,
  timeprior=NULL,
  seed=1234) {

  # Reset seed if required
  if (seed > 0) set.seed(seed)
  lsChain <- dataAndChain$chain
  matT <- dataAndChain$mt
  matH <- dataAndChain$mh

  # Set auxilliary variables and make some checks
  noInds 	<- dim(matT)[1]
  noPs 	<- dim(lsChain$ps)[1]
  headers <- c("sample_ind_right","like",row.names(lsChain$ps),
               paste("tinf.",1:noInds,sep=""))

  # Setup the files - only if not already present
  fnTimes <- paste(runstem,"time.csv",sep="")
  if (fnTimes %in% dir(outputtop) & filemustbeabsent) {
    stop("at least one file already present")
  } else {
    write.table(t(headers),file=paste(outputtop,fnTimes,sep=""),
                row.names=FALSE,col.names=FALSE,sep=",")
  }

  # Simulate timestep updates and accept all
  nots 			<- length(lsChain$vt)
  samplecount 	<- 1
  batchlinecount 	<- 1
  outframe <- matrix(nrow=batchsize,ncol=length(headers),
                     dimnames=list(1:batchsize,headers))

  # Normalize update cumulative
  for (i in 2:length(updatebalance)) {
    updatebalance[i] <- updatebalance[i-1] + updatebalance[i]
  }
  updatebalance <- updatebalance/updatebalance[length(updatebalance)]

  curChain <- lsChain

  curLike <- sum(hks.lnlike.full.1.1(matT,matH,curChain$y0,curChain$vt,
                                     rep(curChain$ps["boost","val"],noInds),
                                     rep(curChain$ps["wane","val"],noInds)))

  while (samplecount <= nosamples) {

    # Write the state of the chain if required
    if ((samplecount %% thin_factor == 0)) {
      outframe[batchlinecount,1] <- samplecount
      outframe[batchlinecount,2] <- curLike
      outframe[batchlinecount,3:(noPs+2)] <- curChain$ps[,"val"]
      outframe[batchlinecount,(noPs+3):(nots+noPs+2)] <- curChain$vt
      batchlinecount <- batchlinecount + 1
    }
    if (batchlinecount > batchsize) {
      write.table(outframe,file=paste(outputtop,fnTimes,sep=""),
                  append=TRUE,col.names=FALSE,row.names=FALSE,sep=",")
      batchlinecount <- 1
    }

    # Decide which update to choose
    # Normalization of updatebalance should make this OK
    rnChoice <- runif(1)
    indChoice <- 1
    while (updatebalance[indChoice] < rnChoice) indChoice <- indChoice + 1

    # Save the current chain
    newChain <- curChain

    # Short range time jump as update one
    if (indChoice==1) {

      newChain$vt <- hks.update.ts(curChain$vt,min_inf_time,
                                   max_inf_time,dtmax=dtmax)
      prop <- 1

      # Parameter updates
    } else if (indChoice==2) {
      newChain$ps[,"val"] <- srg.param.propose.update(curChain$ps)
      prop <- 1

      # Informative time prior updates
      # Function below needs to work of a tmin and tmax with lossibly a large dt ??
    } else if (indChoice==3) {

      # Maybe doesn't work
      # The debug line below have ot have certain properties
      # - for a wide dt, they should return the prior
      # - for a uniform prior, they should not be biased
      # - it must be boundary conditions?

      # XXXXXXXXXXXX
      browser()
      dbarr <- matrix(nrow=1000,ncol=88)
      for (i in 1:1000) {
        tmp <- hks.update.ts.inform(curChain$vt,min_inf_time,max_inf_time,timeprior,dtmax=dtmax)
        dbarr[i,] <- tmp$vt
      }
      # XXXXXXXXXXXX

      tmp <- hks.update.ts.inform(curChain$vt,min_inf_time,max_inf_time,timeprior,dtmax=dtmax)
      newChain$vt <- tmp$vt
      prop <- sum(tmp$pr)

      # For completeness
    } else stop("Problem with ind choice")

    # Test for parameter update or local update
    if (indChoice==1) {

      # Write the standard acceptance for a parameter update
      newLike <- sum(hks.lnlike.full.1.1(
        matT,
        matH,
        newChain$y0,
        newChain$vt,
        rep(newChain$ps["boost","val"],noInds),
        rep(newChain$ps["wane","val"],noInds)))

      # Accept the new chain if necessary
      # This will need to be copied in below for now
      accept <- FALSE
      if (debugacceptall) {
        logaccept <- prop
      } else {
        logaccept <- newLike - curLike + prop
      }
      if (logaccept > 0) {
        accept <- TRUE
      } else if (runif(1) < exp(newLike - curLike)) {
        accept <- TRUE
      }
      if (accept) {
        curChain <- newChain
        curLike <- newLike
      }

      # Accept on a 1-by-1 basis for a individual parameters
      # At the moment, this is just copy-pasted from above, but can be
    } else {

      # Write the standard acceptance for a parameter update
      newLike <- sum(hks.lnlike.full.1.1(
        matT,
        matH,
        newChain$y0,
        newChain$vt,
        rep(newChain$ps["boost","val"],noInds),
        rep(newChain$ps["wane","val"],noInds)))

      # Accept the new chain if necessary
      # This will need to be copied in below for now
      accept <- FALSE
      if (debugacceptall) {
        logaccept <- prop
      } else {
        logaccept <- newLike - curLike + prop
      }
      if (logaccept > 0) {
        accept <- TRUE
      } else if (runif(1) < exp(newLike - curLike)) {
        accept <- TRUE
      }
      if (accept) {
        curChain <- newChain
        curLike <- newLike
      }

    }

    # Increment the sample count
    samplecount <- samplecount + 1

    # Close main sampler loop
  }

  # Test to see if there are leftover samples
  if (batchlinecount > 1) {
    write.table(outframe[1:(batchlinecount-1),],
                file=paste(outputtop,fnTimes,sep=""),
                append=TRUE,col.names=FALSE,row.names=FALSE,sep=",")
  }

}

hks.lasagne.plot <- function(ssDat,strain,fnstem) {

  # Set up some auxilliary variables
  nrows = dim(ssDat)[1]
  lasagcol <- c("white",rev(heat.colors(11)))

  # Define auxilliary routine for each little pixel
  pixel <- function(indper,wave,titre,palette) {
    polygon(
      c(wave-1,wave,wave,wave-1),
      c(indper-1,indper-1,indper,indper),
      col=palette[1+log(titre/5,2)],
      border=NA)
  }

  # Make the main chart
  pdf(paste(fnstem,strain,".pdf",sep=""),height=40/cm(1),width=10/cm(1))
  plot(1:2,ylim=c(0,nrows),xlim=c(0,4),type="n")
  for (i in 1:nrows) {
    pixel(i,1,ssDat[order(ssDat$Age_1)[i],paste(strain,".T1",sep="")],lasagcol)
    pixel(i,2,ssDat[order(ssDat$Age_1)[i],paste(strain,".T2",sep="")],lasagcol)
    pixel(i,3,ssDat[order(ssDat$Age_1)[i],paste(strain,".T3",sep="")],lasagcol)
    pixel(i,4,ssDat[order(ssDat$Age_1)[i],paste(strain,".T4",sep="")],lasagcol)
  }
  dev.off()

  pdf(paste(fnstem,strain,"_legend.pdf",sep=""))
  srg.simple.image.legend(min=0,max=11,colpalette=lasagcol,logscale=FALSE,axlabs=NULL)
  dev.off()

}

hks.month.breaks <- function() {
  y 	<- 2009:2012
  ny 	<- length(y)
  ydash <- NULL
  for (i in y) ydash <- c(ydash,rep(i,12))
  m 	<- rep(c("Jan-","Feb-","Mar-","Apr-","May-","Jun-","Jul-","Aug-","Sep-","Oct-","Nov-","Dec-"),ny)
  d 	<- rep(c("01-"),12*ny)
  rtn <- paste(d,m,ydash,sep="")
}

hks.fig.cards <- function(filestem,dat,strain="H1N1",change="Change") {

  dat$blood_1_date <- as.numeric(dat$blood_1_date) # V added
  dat$blood_2_date <- as.numeric(dat$blood_2_date) # V added
  dat$blood_3_date <- as.numeric(dat$blood_3_date) # V added
  dat$blood_4_date <- as.numeric(dat$blood_4_date) # V added

  min_blood_date <- min(dat$blood_1_date)
  max_blood_date <- max(dat$blood_4_date)

  agecols <- rainbow(100)
  pdf(paste(filestem,strain,"_",change,".pdf",sep=""),height=10/cm(1),width=10/cm(1))
  plot(1:2,ylim=c(0,11),xlim=c(min(dat$blood_1_date),max(dat$blood_4_date)+0),axes=FALSE)
  axis(1)
  axis(2)
  NoChangeMask <- dat[,paste(strain,".T1",sep="")] == dat[,paste(strain,".T2",sep="")] &
    dat[,paste(strain,".T2",sep="")] == dat[,paste(strain,".T3",sep="")] &
    dat[,paste(strain,".T3",sep="")] == dat[,paste(strain,".T4",sep="")]
  ChangeMask <- !(NoChangeMask)
  NumNoChange <- sum(NoChangeMask)
  NumChange <- sum(ChangeMask)
  if (change == "Change") qDatSub <- dat[ChangeMask,] else qDatSub <- dat[NoChangeMask,]
  for (i in 1:NumChange) {
    if (change== "Change") {
      xvals <- c(qDatSub$blood_1_date[i],qDatSub$blood_2_date[i],qDatSub$blood_3_date[i],qDatSub$blood_4_date[i])
    } else {
      xvals <- c(min_blood_date,min_blood_date,max_blood_date,max_blood_date)
    }
    lines(
      xvals,
      log(	c(	qDatSub[i,paste(strain,".T1",sep="")],
              qDatSub[i,paste(strain,".T2",sep="")],
              qDatSub[i,paste(strain,".T3",sep="")],
              qDatSub[i,paste(strain,".T4",sep="")])/5,2) + 0.3*(runif(1)-0.5),
      col=agecols[round(qDatSub$Age_1[i])],
      lwd=0.75)
  }
  dev.off()

  pdf(paste(filestem,strain,"_legend",".pdf",sep=""))
  srg.simple.image.legend(min=0,max=100,colpalette=agecols,logscale=FALSE,axlabs=NULL)
  dev.off()

}

hks.fig.qmh.wave.survival <- function(figfile,datfile) {
  aDatQMH <- hks.read.clean.QMH(datfile)
  test <- hks.aux.survive(aDatQMH$total_positive,waveBounds=c(52,104))
  pdf(figfile)
  plot(aDatQMH$total_positive,type="l",xlim=c(1,104),ylim=c(-100,120),axes=FALSE,main="Wave data and survival functions")
  axis(1)
  axis(2,at=20*(0:6))
  par(new=TRUE)
  plot(test[,1],axes=FALSE,new=TRUE,ylim=c(0,3),type="l",col="red")
  points(test[,2],type="l",col="green")
  axis(2,at=0:10/10)
  dev.off()
}

hks.postproc.mcmc <- function(
  restop,
  resstem,
  chtstem,
  alldat,
  propburn=0.5,
  noparams=0,
  agvar = "ag1",
  wavebounds = as.date(c("01-Apr-2009","01-Dec-2010","31-Dec-2011")),
  suffix="_hks_ii_mcmc.csv") {

  # Load the data
  # XXXX up to here. Need to lad more recent data
  fnRes 	<- paste(restop,resstem,suffix,sep="")
  allres 	<- read.csv(fnRes,header=TRUE)

  # Subset the data
  posdat <- hks.extract.ffold.risers(alldat)

  # Plot whole likelihood
  pdf(paste(chtstem,"mcmc_full_like_",resstem,".pdf",sep=""))
  plot(allres$like,type="l")
  dev.off()

  dimres <- dim(allres)
  nallsamps <- dimres[1]
  npeeps <- dimres[2]-2-noparams
  res <- allres[(nallsamps*propburn):nallsamps,]
  ngoodsamps <- dim(res)[1]
  allinfs <- as.matrix(res[,(3+noparams):(npeeps+2)])

  # Plot whole likelihood
  pdf(paste(chtstem,"mcmc_infs_",resstem,".pdf",sep=""))
  histout <- hist(allinfs,breaks=as.date(hks.month.breaks()),plot = FALSE)
  plot(histout$mids,histout$density,type="l",axes=FALSE)
  axis(1,at=as.date(hks.month.breaks()))
  axis(2)
  dev.off()

  # Calculate age and wave specific attack rates
  # Make mask matrix for the ag of each member of the group
  # Define a function for a single row
  agegroups <- names(table(alldat[,agvar]))
  noagegroups <- length(agegroups)
  nowaves <- length(wavebounds)-1
  vD <- as.vector(table(alldat[,agvar]))

  # Set up a large array for the different age classes, waves and samples
  aR <- array(dim=c(noagegroups,nowaves,ngoodsamps))

  # Calculate the age specific attack rates
  for (i in 1:noagegroups) {

    mask 	<- c(0,0,rep(0,noparams),posdat[,agvar]==agegroups[i])
    agtmp 	<- res[,mask==1]
    for (j in 1:nowaves) {

      # Define waves
      lb 	<- as.numeric(wavebounds[j])
      ub 	<- as.numeric(wavebounds[j+1])

      # Calc vectors
      tmp <- (agtmp >= lb) & (agtmp < ub)
      aR[i,j,] <- apply(tmp,1,sum)

    }
  }

  # Make age group data table
  vP <- as.vector(c(1119*1000, 2105*1000, 2367*1000, 1260*1000))
  vP <- vP/sum(vP)
  arDT = data.frame(
    w1_all=aR[1,1,]/vD[1]*vP[1]+aR[2,1,]/vD[2]*vP[2]+aR[3,1,]/vD[3]*vP[3]+aR[4,1,]/vD[4]*vP[4],
    w2_all=aR[1,2,]/vD[1]*vP[1]+aR[2,2,]/vD[2]*vP[2]+aR[3,2,]/vD[3]*vP[3]+aR[4,2,]/vD[4]*vP[4],
    w1_u19=aR[1,1,]/vD[1],
    w1_19_44=aR[2,1,]/vD[2],
    w1_45_64=aR[3,1,]/vD[3],
    w1_65p=aR[4,1,]/vD[4],
    w2_u19=aR[1,2,]/vD[1],
    w2_19_44=aR[2,2,]/vD[2],
    w2_45_64=aR[3,2,]/vD[3],
    w2_65p=aR[4,2,]/vD[4])

  # Box plot the age specific attack rates
  pdf(paste(chtstem,"mcmc_w1w2_AR_box_",resstem,".pdf",sep=""),width=20/cm(1),height=15/cm(1))
  boxplot(
    arDT,
    ylim = c(0,0.7),
    col=c("red","green",rep("red",4),rep("green",4)),
    names=c("All","All","<19","19-44","45-64","65+","<19","19-44","45-64","65+"),
    range=0,
    xlab="Age group",
    ylab="Attack rate"
  )
  legend(9,0.7,c("Wave 1","Wave 2"),fill=c("red","green"))
  dev.off()

}
# hks.postproc.mcmc("/Users/sriley/Dropbox/results/","20121025",fDat,qDat)

hks.postproc.mcmc.v1.1 <- function(
  restop,
  resstem,
  chtstem,
  alldat,
  propburn=0.5,
  noparams=0,
  agvar = "ag1",
  wavebounds = as.date(c("01-Apr-2009","01-Dec-2010","31-Dec-2011")),
  suffix="_hks_ii_mcmc.csv") {

  # Load the data
  # XXXX up to here. Need to lad more recent data
  fnRes 	<- paste(restop,resstem,suffix,sep="")
  allres 	<- read.csv(fnRes,header=TRUE)

  # Subset the data
  posdat <- hks.extract.ffold.risers(alldat)

  # Plot whole likelihood
  pdf(paste(chtstem,"mcmc_full_like_",resstem,".pdf",sep=""))
  plot(allres$like,type="l")
  dev.off()

  dimres <- dim(allres)
  nallsamps <- dimres[1]
  npeeps <- dimres[2]-2-noparams
  res <- allres[(nallsamps*propburn):nallsamps,]
  ngoodsamps <- dim(res)[1]
  allinfs <- as.matrix(res[,(3+noparams):(npeeps+2)])

  # Plot whole likelihood
  pdf(paste(chtstem,"mcmc_infs_",resstem,".pdf",sep=""))
  histout <- hist(allinfs,breaks=as.date(hks.month.breaks()),plot = FALSE)
  plot(histout$mids,histout$density,type="l",axes=FALSE)
  axis(1,at=as.date(hks.month.breaks()))
  axis(2)
  dev.off()

  # Calculate age and wave specific attack rates
  # Make mask matrix for the ag of each member of the group
  # Define a function for a single row
  agegroups <- names(table(alldat[,agvar]))
  noagegroups <- length(agegroups)
  nowaves <- length(wavebounds)-1
  vD <- as.vector(table(alldat[,agvar]))

  # Set up a large array for the different age classes, waves and samples
  aR <- array(dim=c(noagegroups,nowaves,ngoodsamps))

  # Calculate the age specific attack rates
  for (i in 1:noagegroups) {

    mask 	<- c(0,0,rep(0,noparams),posdat[,agvar]==agegroups[i])
    agtmp 	<- res[,mask==1]
    for (j in 1:nowaves) {

      # Define waves
      lb 	<- as.numeric(wavebounds[j])
      ub 	<- as.numeric(wavebounds[j+1])

      # Calc vectors
      tmp <- (agtmp >= lb) & (agtmp < ub)
      aR[i,j,] <- apply(tmp,1,sum)

    }
  }

  # Make age group data table
  vP <- as.vector(c(1119*1000, 2105*1000, 2367*1000, 1260*1000))
  vP <- vP/sum(vP)

  arDT = data.frame(
    w1_all=aR[1,1,]/vD[1]*vP[1]+aR[2,1,]/vD[2]*vP[2]+aR[3,1,]/vD[3]*vP[3]+aR[4,1,]/vD[4]*vP[4],
    w2_all=aR[1,2,]/vD[1]*vP[1]+aR[2,2,]/vD[2]*vP[2]+aR[3,2,]/vD[3]*vP[3]+aR[4,2,]/vD[4]*vP[4],
    w1_u19=aR[1,1,]/vD[1],
    w1_19_44=aR[2,1,]/vD[2],
    w1_45_64=aR[3,1,]/vD[3],
    w1_65p=aR[4,1,]/vD[4],
    w2_u19=aR[1,2,]/vD[1],
    w2_19_44=aR[2,2,]/vD[2],
    w2_45_64=aR[3,2,]/vD[3],
    w2_65p=aR[4,2,]/vD[4])

  # Box plot the age specific attack rates
  pdf(paste(chtstem,"mcmc_w1w2_AR_box_",resstem,".pdf",sep=""),width=20/cm(1),height=15/cm(1))
  boxplot(
    arDT,
    ylim = c(0,0.7),
    col=c("red","green",rep("red",4),rep("green",4)),
    names=c("All","All","<19","19-44","45-64","65+","<19","19-44","45-64","65+"),
    range=0,
    xlab="Age group",
    ylab="Attack rate"
  )
  legend(9,0.7,c("Wave 1","Wave 2"),fill=c("red","green"))
  dev.off()

}
# hks.postproc.mcmc("/Users/sriley/Dropbox/results/","20121025",fDat,qDat)

# Fit the chain either to the actual data or the simulated data
# Assumes that the ts data and the study data are present in pwd
# runmode is one of "prior","data","sim_actual","sim_perfect"
hks.mcmc.main <- function(
  studDat,
  fnstem="debug",
  agerange=c(0,110),
  nosamples=100,
  thinfactor=1,
  timeprior=NULL,
  runmode="data",
  updatebalance=c(1,1,0),
  paramfile="params.csv",
  checkfile=TRUE,
  min_blood_date = as.date("10-Apr-2009"),
  setseed=-1,
  dtmax = 30) {

  set.seed(setseed)

  if (runmode=="prior") accall <- TRUE else accall <- FALSE

  # Set some aux parameters
  max_blood_date <- max(studDat$blood_4_date)

  # Prep the actual data
  ChainAndData <- hks.mcmc.setup.1.2(studDat,min_blood_date,max_blood_date,paramfile,agerange=agerange)

  # Simulate times of infection using the experiemtnal protocol from the actual data
  # Range found by inspection of the max min times
  if (runmode=="sim_actual" || runmode=="sim_perfect") {
    range_from_chp_data <- 64:83
    vecInfs <- hks.sim.incidence.months(dim(ChainAndData$mt)[1],aDatCHP$month_starting[range_from_chp_data],aDatCHP$month_ending[range_from_chp_data],aDatCHP$TotalA[range_from_chp_data])
  }

  # Simulate the effective antibody concentrations in the same format
  # Control for the timing of measures
  contmeasure <- seq(as.numeric(min_blood_date),as.numeric(max_blood_date),by=7)
  if (runmode=="sim_perfect") {
    measureTimes <- matrix(contmeasure,nrow=length(vecInfs),ncol=length(contmeasure),byrow=TRUE)
  } else {
    measureTimes <- ChainAndData$mt
  }

  if (runmode=="sim_actual" || runmode=="sim_perfect") {
    ChainAndData <- hks.sim.study(measureTimes,vecInfs,t_min=min_blood_date,t_max=max_blood_date,paramfile)
  }

  # If using informative priors, setup the cumulative incidence
  # Put in a single dummy incidence value if there are zeros
  if (!is.null(timeprior)) {
    timeprior$inc <- ifelse(timeprior$inc < 1e-100,1,timeprior$inc)
    timeprior$cuminc <- timeprior$inc
    notps <- dim(timeprior)[1]
    for (i in 2:notps) timeprior$cuminc[i] <- timeprior$cuminc[i-1] + timeprior$inc[i]
    timeprior$cuminc <- timeprior$cuminc / timeprior$cuminc[notps]
  }

  hks.mcmc.1.1(
    ChainAndData,
    min_blood_date-21,
    max_blood_date,
    updatebalance=updatebalance,
    dtmax=dtmax,
    timeprior=timeprior,
    debugacceptall=accall,
    nosamples=nosamples,
    thin_factor=thinfactor,
    batchsize = 100,
    runstem=fnstem,
    filemustbeabsent=checkfile)

}

# By default, assumes it is beiung run in the directory where the results are generated
# and makes the plots there as well
# setwd("/Volumes/HPOffice/Results/2012-12-23")
# setwd("~/Dropbox/tmp/2012-12-23")
# source("~/Dropbox/svneclipse/idsource/R/stevensRfunctions.R")
# hks.postproc.mcmc.v2()
hks.postproc.mcmc.v2 <- function(
  #  	fns = c("ag1_sim_time.csv","ag2_sim_time.csv","ag3_sim_time.csv","ag4_simtime.csv"),
  fns = c("ag1_time.csv","ag2_time.csv","ag3_time.csv","ag4_time.csv"),
  burnin = 0.5) {

  require(gplots)

  # Number of subsets
  nofiles <- length(fns)

  # Read first file
  f1 <- read.csv(fns[1])
  nof1 <- dim(f1)[1]
  fsamp <- ceiling(burnin*nof1)
  firstinf <- match("tinf.1",names(f1))

  noheads <- dim(f1)[2]

  # Setup data
  psamp <- cbind(ag=rep(1,nof1),pindex=1:nof1,f1[,1:(firstinf-1)])
  infsamp <- f1[,firstinf:noheads]
  names(infsamp) <- paste("ag1.",names(infsamp),sep="")

  # Go through the files and make large data bases of the parameter samples
  # and the infection times
  if (nofiles > 1) {
    for (i in 2:nofiles) {
      tmpdat <- read.csv(fns[i])
      tmpfirstinf <- match("tinf.1",names(tmpdat))
      notmp <- dim(tmpdat)[1]
      headstmp <- dim(tmpdat)[2]
      if (tmpfirstinf != firstinf) stop("second file not right")
      if (notmp != nof1) stop("Files have to be the same length")
      tmppsamp <- cbind(ag=rep(i,notmp),pindex=1:nof1,tmpdat[,1:(firstinf-1)])
      psamp <- rbind(psamp,tmppsamp)
      tmpinf <- tmpdat[,firstinf:headstmp]
      names(tmpinf) <- paste("ag",i,".",names(tmpinf),sep="")
      infsamp <- cbind(infsamp,tmpinf)
    }
  }

  # Set up some graphics defaults
  colseq <- c("red","green","blue","cyan","magenta","yellow","black")

  # Make boosting and waning plots
  for (i in 1:nofiles) {
    pdf(paste("boost_wane_ag",i,".pdf",sep=""))
    hist2d(	psamp[psamp$ag==i & psamp$pindex > fsamp,"boost"],
            psamp[psamp$ag==i & psamp$pindex > fsamp,"wane"],
            nbins=80,
            col = c("white",heat.colors(200)),
            xlab = "Boosting over 3 weeks",
            ylab = "Waning per year",
            xlim=c(0,8),
            ylim=c(0,3)
    )
    dev.off()
  }

  # Make horizonal box plots for times of infection
  # See if any are localized
  # Consider sorting by mode
  # Not sure that this plot is useful
  pdf("infection_box_plot.pdf",height=25/cm(1),width=15/cm(1))
  sortedinfs <- order(apply(infsamp,median,MARGIN=2))
  boxplot(infsamp,horizontal=TRUE,axes=FALSE)
  axis(1)
  axis(2,las=2)
  dev.off()

}

hks.postproc.mcmc.v2.1 <- function(
  alldat,
  agvar,
  incdata,
  #		fns = c("ag1_sim_time.csv","ag2_sim_time.csv","ag3_sim_time.csv","ag4_simtime.csv"),
  fns = c("ag1_time.csv","ag2_time.csv","ag3_time.csv","ag4_time.csv"),
  resdir = "./",
  resstem="debug",
  wavebounds = as.date(c("01-Apr-2009","01-Dec-2010","31-Dec-2011")),
  vP = as.vector(c(1119*1000, 2105*1000, 2367*1000, 1260*1000)),
  burnin = 0.5) {

  require(gplots)

  # Number of subsets
  nofiles <- length(fns)

  # Read first file
  f1 <- read.csv(paste(resdir,paste(fns[1]),sep=""))
  nof1 <- dim(f1)[1]
  fsamp <- ceiling(burnin*nof1)
  firstinf <- match("tinf.1",names(f1))

  noheads <- dim(f1)[2]

  # Setup data
  psamp <- cbind(ag=rep(1,nof1),pindex=1:nof1,f1[,1:(firstinf-1)])
  infsamp <- f1[,firstinf:noheads]
  names(infsamp) <- paste("ag1.",names(infsamp),sep="")

  # Go through the files and make large data bases of the parameter samples
  # and the infection times
  if (nofiles > 1) {
    for (i in 2:nofiles) {
      tmpdat <- read.csv(paste(resdir,paste(fns[i]),sep=""))
      tmpfirstinf <- match("tinf.1",names(tmpdat))
      notmp <- dim(tmpdat)[1]
      headstmp <- dim(tmpdat)[2]
      if (tmpfirstinf != firstinf) stop("second file not right")
      if (notmp != nof1) stop("Files have to be the same length")
      tmppsamp <- cbind(ag=rep(i,notmp),pindex=1:nof1,tmpdat[,1:(firstinf-1)])
      psamp <- rbind(psamp,tmppsamp)
      tmpinf <- tmpdat[,firstinf:headstmp]
      names(tmpinf) <- paste("ag",i,".",names(tmpinf),sep="")
      infsamp <- cbind(infsamp,tmpinf)
    }
  }

  # Set up some graphics defaults
  colseq <- c("red","green","blue","cyan","magenta","yellow","black")

  # Make boosting and waning plots
  for (i in 1:nofiles) {
    pdf(paste("./auxilliary/boost_wane_ag",i,".pdf",sep=""))
    hist2d(	psamp[psamp$ag==i & psamp$pindex > fsamp,"boost"],
            psamp[psamp$ag==i & psamp$pindex > fsamp,"wane"],
            nbins=80,
            col = c("white",heat.colors(200)),
            xlab = "Boosting over 3 weeks",
            ylab = "Waning per year",
            xlim=c(0,8),
            ylim=c(0,3)
    )
    dev.off()
  }

  # Make horizonal box plots for times of infection
  # See if any are localized
  # Consider sorting by mode
  # Not sure that this plot is useful and takes a while
  # pdf("auxilliary/infection_box_plot.pdf",height=25/cm(1),width=15/cm(1))
  # plot(1:2,main="chart not useful akljshda?")
  # sortedinfs <- order(apply(infsamp,median,MARGIN=2))
  # boxplot(infsamp,horizontal=TRUE,axes=FALSE)
  # axis(1)
  # axis(2,las=2)
  # dev.off()

  # Plot histogram of infections with a weekly density
  pdf("auxilliary/infection_histogram.pdf",height=25/cm(1),width=15/cm(1))
  trange <- as.date(c("01-Jan-2009","31-Dec-2011"))

  # Convert observed gp incidence into a density
  incrange <- incdata[,1] >= trange[1] & incdata[,1] <= trange[2]
  incdens <- incdata[,2] / sum(incdata[incrange,2]) / 7

  # Make the plot
  atvals <- seq(as.numeric(trange[1]),as.numeric(trange[2]),60)
  par(mai=(c(0,0,0,0)))
  par(fig=c(0.15,0.99,0.15,0.99))
  infhist <- hist(
    as.matrix(infsamp),
    breaks=seq(as.numeric(trange[1]),as.numeric(trange[2]),7),
    plot=FALSE)
  plot(	infhist$mids,
        infhist$density,
        type="l",
        col="red",
        lwd=2,
        axes=FALSE,
        ylim=c(0,max(c(incdens,infhist$denstiy))))
  points(incdata[,1],incdens,type="l")
  axis(2)
  axis(1,las=2,at=atvals,labels=as.date(atvals),padj=5)
  mtext("Date",1,6)
  mtext("Incidence (density per day)",2,3)
  dev.off()

  # Prep the empirical data for attack rate calculations
  posdat <- hks.extract.ffold.risers(alldat)

  # Calculate age and wave specific attack rates
  # Make mask matrix for the ag of each member of the group
  # Define a function for a single row
  ngoodsamps <- dim(infsamp)[1]
  agegroups <- names(table(alldat[,agvar]))
  noagegroups <- length(agegroups)
  nowaves <- length(wavebounds)-1
  vD <- as.vector(table(alldat[,agvar]))

  # Set up a large array for the different age classes, waves and samples
  aR <- array(dim=c(noagegroups,nowaves,ngoodsamps))

  # Calculate the age specific attack rates
  for (i in 1:noagegroups) {

    mask 	<- posdat[,agvar]==agegroups[i]
    agtmp 	<- infsamp[,mask==1]

    for (j in 1:nowaves) {

      # Define waves
      lb 	<- as.numeric(wavebounds[j])
      ub 	<- as.numeric(wavebounds[j+1])

      # Calc vectors
      tmp <- (agtmp >= lb) & (agtmp < ub)
      aR[i,j,] <- apply(tmp,1,sum)

    }
  }


  # Make age group data table
  vP <- vP/sum(vP)

  arDT = data.frame(
    w12_all=(aR[1,1,]+aR[1,2,])/vD[1]*vP[1]+(aR[2,1,]+aR[2,2,])/vD[2]*vP[2]+(aR[3,1,]+aR[3,2,])/vD[3]*vP[3]+(aR[4,1,]+aR[4,2,])/vD[4]*vP[4],
    w1_all=aR[1,1,]/vD[1]*vP[1]+aR[2,1,]/vD[2]*vP[2]+aR[3,1,]/vD[3]*vP[3]+aR[4,1,]/vD[4]*vP[4],
    w2_all=aR[1,2,]/vD[1]*vP[1]+aR[2,2,]/vD[2]*vP[2]+aR[3,2,]/vD[3]*vP[3]+aR[4,2,]/vD[4]*vP[4],
    w1_u19=aR[1,1,]/vD[1],
    w1_19_44=aR[2,1,]/vD[2],
    w1_45_64=aR[3,1,]/vD[3],
    w1_65p=aR[4,1,]/vD[4],
    w2_u19=aR[1,2,]/vD[1],
    w2_19_44=aR[2,2,]/vD[2],
    w2_45_64=aR[3,2,]/vD[3],
    w2_65p=aR[4,2,]/vD[4])

  # Box plot the age specific attack rates
  pdf(paste("auxilliary/mcmc_w1w2_AR_box_",resstem,".pdf",sep=""),width=10/cm(1),height=15/cm(1))
  par(cex=10/12)
  par(mai=(c(0,0,0,0)))
  par(fig=c(0.20,0.99,0.25,0.99))
  boxplot(
    arDT,
    ylim = c(0,0.7),
    col=c("yellow","red","green",rep("red",4),rep("green",4)),
    axes=FALSE,
    range=0,
    xlab="Age group",
    ylab="Attack rate"
  )
  axis(1,las=2,labels=names(arDT),at=1:length(names(arDT)))
  axis(2)
  legend(7,0.7,c("Wave 1","Wave 2"),fill=c("red","green"))
  dev.off()

  # Setup table output
  norows <- dim(arDT)[2]
  tabSum <- data.frame(
    row.names=names(arDT),
    median=apply(arDT,2,median),
    lower95=apply(arDT,2,quantile,probs=c(0.025)),
    upper95=apply(arDT,2,quantile,probs=c(0.975)))
  write.csv(tabSum,file="auxilliary/summary_mcmc_AR.csv")

  # Make boosting and waning plots
  pdf(paste("./auxilliary/w1_w1_ar_corr.pdf",sep=""))
  hist2d(	arDT$w1_all,
          arDT$w2_all,
          nbins=20,
          col = c("white",heat.colors(200)),
          xlab = "Wave 1",
          ylab = "Wave 2",
          xlim=c(0.1,0.3),
          ylim=c(0.1,0.3)
  )
  dev.off()


  # Plot 2D histogram

}

# Gives the cumulative density of the present time
# Needs to be reasonably efficient
# tvec <- c(0,1,2,3)
# dvec <- c(10,20,30,40)
# hks.cuminc.density(2.9999,tvec,dvec)
hks.cuminc.density <- function(pt,vectimes,vecdense) {
  notps <- length(vectimes)
  dt <- vectimes[2]-vectimes[1] # Because the first time slot is sometimes not correct!
  t0 <- vectimes[1]
  if (pt < t0 | pt >= vectimes[notps]) {
    browser()
    stop("pt out of range in hks.cuminc.density")
  }
  baseindex <- as.integer((pt - t0) / dt)
  basetime <- vectimes[baseindex]
  basedens <- vecdense[baseindex]
  rtn <- basedens + (vecdense[baseindex+1]-basedens) * (as.numeric(pt)-as.numeric(basetime)) / dt
  rtn
}

# Needs to return the new time for a known start time and for the proposed
# change in the cumulative density
# this can then be incorporated into the variable time proposal
# This would be a great candidate for NRcppGSL
# tvec <- c(0,1,2,3)
# dvec <- c(0.0,0.2,0.6,1.0)
# hks.cuminc.time(2.5,0.4,tvec,dvec)
hks.cuminc.time <- function(ts,ddens,vectimes,vecdens) {

  notps <- length(vectimes)
  dt <- vectimes[2]-vectimes[1]
  t0 <- vectimes[1]
  if (ts < t0 | ts >= vectimes[notps]) stop("pt out of range in hks.cuminc.density")
  currentindex <- as.integer((ts - t0) / dt) + 1
  currentindextime <- vectimes[currentindex]
  currenttstime <- ts - currentindextime
  basedense <- vecdens[currentindex]
  currentdense <- basedense+(vecdens[currentindex+1]-basedense) * (ts-currentindextime) / dt
  targetdens <- currentdense + ddens
  if (targetdens < 0 | targetdens > 1) stop("Requested density is out of bounds in hks.cuninc.time")
  if (ddens >= 0) {
    # This falls over the first time. needs sorting here.
    while (currentdense+ddens > vecdens[currentindex+1]) {
      ddens <- ddens - (vecdens[currentindex+1]-currentdense)
      currentindex <- currentindex + 1
      currentdense <- vecdens[currentindex]
      currentindextime <- vectimes[currentindex]
      currenttstime <- 0
    }
    rtn <- as.numeric(currentindextime)+currenttstime+ddens/(vecdens[currentindex+1]-vecdens[currentindex])*dt
  } else {
    stop("ddens must be positive")
  }
  rtn
}

hks.round.CARs <- function(vecN,
                           vecn,
                           vecPop,
                           latex=TRUE,
                           pretex="& ",
                           xlabs=c("2-18","19-44","45-64","65+"),
                           width=5/cm(1),
                           height=8/cm(1),
                           filename="~/Dropbox/tmp/hks_round_CARs.pdf",
                           ymax=0.6) {

  nominal.overall = round(sum (vecn * (sum(vecN)*vecPop)  / vecN))

  vecN <- c(sum(vecN),vecN)
  vecn <- c(nominal.overall,vecn)
  xlabs <- c("Overall",xlabs)

  pe <- vecn/vecN
  lb <- srg.ci.binom(vecN,vecn,P=0.975,min=0)
  ub <- srg.ci.binom(vecN,vecn,P=0.025,min=0)

  nopoints <- length(vecN)
  xvals <- 1:nopoints

  pdf(file=filename,width=width,height=height,useDingbats=FALSE)
  par(cex=10/12)
  delta<-0.2
  plot(1:2,type="n", xlim=c(1-delta,nopoints+delta), ylim=c(0,ymax), log="",axes=FALSE)
    axis(1,at=c(xvals),labels=c(xlabs),las=2)
    axis(2,las=1)
    grid(col="lightgrey",lty=1)
    points(xvals,pe,pch=19)
    arrows(xvals,pe,y1=ub,angle=90,length=0.05)
    arrows(xvals,pe,y1=lb,angle=90,length=0.05)
  dev.off()

  if (latex) {
    pe <- round(pe,2)
    lb <- round(lb,2)
    ub <- round(ub,2)
    cat(pretex);for (i in 1:length(vecn)) cat("& ",vecn[i]," ",sep="");cat(" \\\\ \n")
    cat(pretex);for (i in 1:length(vecn)) cat("& ",pe[i],"\\ (",lb[i],"\\, ",ub[i],"\\) ",sep="");cat(" \\\\ \n")
  }
}

hks.coinf.arrows <- function(dat,fname,height=12/cm(1),width=12/cm(1)) {
  # Investigate those that seroconvert to both strains
  pdf(file=fname,height=height,width=width)
  plot(1:2,xlim=c(0,8),ylim=c(0,8),type="n",xlab="Pandemic Titre",ylab="Seasonal Titre",axes=FALSE)
  axis(1,at=0:8,labels=c("<1:10","1:10","1:20","1:40","1:80","1:160","1:320","1:640","1:1280"),las=2)
  axis(2,at=0:8,labels=c("<1:10","1:10","1:20","1:40","1:80","1:160","1:320","1:640","1:1280"),las=2)
  jam <- 0.1
  grid()
  arrows(	jitter(dat[dat$ff_r12_both,"pT1"],amount=jam),
          jitter(dat[dat$ff_r12_both,"sT1"],amount=jam),
          x1=jitter(dat[dat$ff_r12_both,"pT2"],amount=jam),
          y1=jitter(dat[dat$ff_r12_both,"sT2"],amount=jam),
          col="red")
  arrows(	jitter(dat[dat$ff_r23_both,"pT2"],amount=jam),
          jitter(dat[dat$ff_r23_both,"sT2"],amount=jam),
          x1=jitter(dat[dat$ff_r23_both,"pT3"],amount=jam),
          y1=jitter(dat[dat$ff_r23_both,"sT3"],amount=jam),
          col="blue")
  arrows(	jitter(dat[dat$ff_r34_both,"pT3"],amount=jam),
          jitter(dat[dat$ff_r34_both,"sT3"],amount=jam),
          x1=jitter(dat[dat$ff_r34_both,"pT4"],amount=jam),
          y1=jitter(dat[dat$ff_r34_both,"sT4"],amount=jam),
          col="green")
  legend(0,8,c("Rounds 1 to 2","Rounds 2 to 3","Rounds 3 to 4"),col=c("red","blue","green"),lty=1,lwd=2)
  dev.off()
}

hks.check.single.tinf <- function(vect,range=c(as.date("01-Jan-2009"),as.date("31-Dec-2011"))) {
  plot(vect,ylim=range,type="l")
}

# Function to make a table of esitmated rates of clinical disease
# per infection
# Auxcases overrides the surveillance data if not null
hks.wave.strain.severity <- function(
  fnChart,
  fnTable,
  tabDat,
  vecPop,
  auxcases,
  ylims=c(0.2,20),
  totalpop = 6912020) {

  require(xtable)
  require(PropCIs)

  rownames <- c(
    "ff_r12_ph1n1","ff_r23_ph1n1","ff_r34_ph1n1",
    "ff_r12_sh3n2","ff_r23_sh3n2","ff_r34_sh3n2")
  surCols <- c("gp.h1n1","gp.h1n1","gp.h1n1","gp.h3n2","gp.h3n2","gp.h3n2")
  vecStud <- table(tabDat$ag1)

  waveBounds <- c(
    median(tabDat$blood_1_date),
    median(tabDat$blood_2_date),
    median(tabDat$blood_3_date),
    median(tabDat$blood_4_date))

  startIndex <- c(1,2,3,1,2,3)
  endIndex   <- c(2,3,4,2,3,4)
  colnames <- c("inf","case","ratio","lb","ub")
  rtn <- matrix(
    nrow=length(rownames),
    ncol=length(colnames),
    dimnames=list(rownames,colnames))
  rtn[] <- -1

  N_S <- sum(vecStud)

  # Loop around the require outcomes
  for (i in 1:length(rownames)) {

    # All the necessary calculations are in the following block
	# XXXX Need to check exactly what is going on here
	# Needs the
    startdate <- waveBounds[startIndex[i]]
    enddate <- waveBounds[endIndex[i]]
    cases <- auxcases[i]
    vecn <- table(tabDat[,rownames[i]],tabDat$ag1)["TRUE",]
    infections <- round(sum (vecn * (N_S*vecPop)  / vecStud))
    peinf <- infections / N_S
    lbinf <- srg.ci.binom(N_S,infections,P=0.975,min=0)
    ubinf <- srg.ci.binom(N_S,infections,P=0.025,min=0)

    # Now they just need putting into the correct columns
    rtn[i,"inf"] 	<- peinf * totalpop
    rtn[i,"case"] 	<- cases
    rtn[i,"ratio"] 	<- cases / (peinf * totalpop) * 10000
    rtn[i,"lb"] 	<- cases / (ubinf * totalpop) * 10000
    rtn[i,"ub"] 	<- cases / (lbinf * totalpop) * 10000

  }

  # Calculate the ratio for each wave -> this should be printed to the
  # screen to go into the text

  cat("TEXT: Estimates of the ratio of severity between the two different waves. \n")
  startindex <- c(1,2,3,2,1)
  endindex <- c(1,2,3,3,3)
  nrows <- length(startindex)
  names <- paste("from_",startindex,"_to_",endindex+1)
  ratiotab <- data.frame(pe=rep(NA,nrows),lb=rep(NA,nrows),ub=rep(NA,nrows),row.names=names)

  for (i in 1:length(startindex)) {
    range <- startindex[i]:endindex[i]
    n_p <- sum(rtn[range,"inf"] * N_S / totalpop)
    c_p <- sum(rtn[range,"case"])
    n_s <- sum(rtn[3+range,"inf"] * N_S / totalpop)
    c_s <- sum(rtn[3+range,"case"])
    est <- c(n_p/n_s,riskscoreci(n_p,N_S,n_s,N_S,0.95)$conf.int)*c_s/c_p
    ratiotab[i,] <- est
    cat("From round",startindex[i],"to round",endindex[i]+1,est,"\n")
    flush.console()
  }

#  pdf(file=filename,width=width,height=height,useDingbats=FALSE)
#  delta<-0.2
#  plot(1:2,type="n", xlim=c(1-delta,nopoints+delta), ylim=c(0,ymax), log="",axes=FALSE)
#  axis(1,at=c(xvals),labels=c(xlabs),las=2)
#  axis(2,las=1)
#  points(xvals,pe,pch=19)
#  arrows(xvals,pe,y1=ub,angle=90,length=0.05)
#  arrows(xvals,pe,y1=lb,angle=90,length=0.05)
#  dev.off()

  # Make the chart of the ratios
  nopoints <- dim(ratiotab)[1]
  xlims <- c(0,nopoints+1)
  pdf.figure.proof(findex="",file=fnChart,pw=4,ph=7,xpd=FALSE,textscale=10/12)
    plot(1:2,xlim=xlims,ylim=ylims,axes=FALSE,log="y",type="n")
    lines(c(0,nopoints+1),c(1,1),lty="dashed",col="black",lwd=2)
    axis(2,las=2)
    axis(1,at=1:nopoints,labels=row.names(ratiotab),las=2)
    axis(1,at=0:(nopoints+1),labels=rep("",length(0:(nopoints+1))))
    for (x in 1:nopoints) abline(v=x,col="lightgrey")
    for (y in axTicks(2)) abline(h=y,col="lightgrey")
    for (i in 1:nopoints) {
      points(i,ratiotab[i,"pe"],pch=19)
      arrows(c(i,i),ratiotab[i,"pe"],y1=ratiotab[i,"ub"],angle=90,length=0.05)
      arrows(c(i,i),ratiotab[i,"pe"],y1=ratiotab[i,"lb"],angle=90,length=0.05)
    }
  dev.off()

  # Make the table this can be tuned to be prettier once the results are fixed
  # Formatted for the current latex table
  xrtn <- xtable(srg.print.format(rtn))
  norowsxtab <- dim(xrtn)[1]
  for (i in 1:norowsxtab) {
    cat(
      "& & ",
      row.names(xrtn)[i]," & ",
      as.character(xrtn[i,1])," & ",
      as.character(xrtn[i,2])," & ",
      as.character(xrtn[i,3]), " (",
      as.character(xrtn[i,4]),", ",
      as.character(xrtn[i,5]),") \\\\ \n",sep="")
  }

  # print(xrtn,file=fnTable)

  # Return the table
  list(strainspec=rtn,ratios=ratiotab)

}


hks.plot.sev.ratios <- function(filename) {
  pdf(filename)
  plot(1:2)
  dev.off()
}

# Loads up the contact-only data from the serological data and creates some auxilliary data structures
hkc.load.and.clean.hksero.contact <- function(fn="/Volumes/NO NAME/data/influenza/hk_serosurvey/tbase_0722b.csv") {
  rtn <- read.csv(file=fn, head=TRUE)
  rtn <- rtn[rtn[,"A_1" ] > -1,]
  rtn <- rtn[rtn[,"A_2" ] > -1,]
  rtn <- rtn[rtn[,"A_3" ] > -1,]
  rtn <- rtn[rtn[,"A_4" ] > -1,]
  rtn <- rtn[rtn[,"n_location" ] > -1,]
  rtn <- rtn[rtn[,"n_contactp" ] > -1,]
  rtn <- rtn[rtn[,"conage1" ] > -1,]
  rtn <- rtn[rtn[,"conage2" ] > -1,]
  rtn <- rtn[rtn[,"conage3" ] > -1,]
  rtn <- rtn[rtn[,"conage4" ] > -1,]
  rtn$child_p <- relevel(factor(rtn$child_p),ref=2)
  rtn$A234 <-rtn$A_2+rtn$A_3+rtn$A_4
  rtn$A34 <-rtn$A_3+rtn$A_4
  rtn$conage12<-rtn$conage1+rtn$conage2
  rtn$conage123<-rtn$conage1+rtn$conage2+rtn$conage3
  rtn
}
