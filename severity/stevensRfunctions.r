binCINew <- function(N,n,P,min=0.000001) {
	rtn <- NA
	if (N>0) {
		if (N==n) n <- n-1
		if (n==0 && P > 0.5) rtn <- min
		else rtn <- (uniroot(srg.ci.binom.optim,c(0,1),P=P,N=N,n=n))$root
	}
	rtn
}

srg.bin.ci <- function(N,n,P,min=0.000001) {
	rtn <- NA
	if (N>0) {
		if (N==n) n <- n-1
		if (n==0 && P > 0.5) rtn <- min
		else rtn <- (uniroot(srg.ci.binom.optim,c(0,1),P=P,N=N,n=n))$root
	}
	rtn
}

srg.make.binci.table <- function(numerator,denominator,rounding=3) {
	nosamps <- length(numerator)
	ub <- srg.ci.binom(denominator,numerator,P=0.025)
	lb <- srg.ci.binom(denominator,numerator,P=0.975)
	pe <- numerator/denominator
	round(data.frame(n=numerator,N=denominator,pe=pe,lb=lb,ub=ub),rounding)
}

# Computes Binomial 95% Confidence Interval
srg.ci.binom <- function(vecN,vecn,P,min=0) {
    binom.optim <- function(p,P,N,n) {
    guess <- pbinom(n,N,p,lower.tail=TRUE)
    rtn <- P-guess
    rtn
  }
  
	noObs <- length(vecN)
	rtn <- array(NA,c(noObs))
	for (i in 1:noObs) {
		if (vecN[i] > 0) {
			sol <- uniroot(binom.optim,c(0,1),P=P,N=vecN[i],n=vecn[i])
			if (sol$root > min) rtn[i] <- sol$root
			else rtn[i] <- min
		} else rtn[i] <- min
	}
	rtn
}

plot.symptoms <- function(df, N, filename, ylim_in = c(0,1.0), ms=0.3, ...) {
	
	# Define the number of ...
	fw <- 8.3/cm(1)
	fh <- 8.3/cm(1)
	
	
	# Open the pdf file
	pdf(filename,height=fh,width=fw)
	
	# Set some standard parameter options
	par(mai=(c(0.15*fh,0.175*fw,0,0)))
	par(mgp=c(2,0.4,0))
	par(tcl=-0.25)
	par(cex = 10/12)
	
	# Set up the vectors
	vecLev1 <- names(df)
	noLev1 <- length(vecLev1)
	offvec1 <- (1:noLev1 - 0.5)
	
	vecLev2 <- row.names(df)
	noLev2 <- length(vecLev2)
	offvec2 <- ((0:(noLev2-1))/(noLev2-1)*ms - ms/2)
	
	colvec <- c("red","green","blue","cyan","magenta")
	
	
	# Setup the axes
	plot(1:2,type="n", xlim=c(0,noLev1), ylim=ylim_in, axes=FALSE)
	axis(2,las=1)
	veclabs <- rep(" ",noLev1*2+1)
	for (i in 1:noLev1) veclabs[i*2] <- vecLev1[i]
	axis(1,at=(0:(noLev1*2))/2,labels=veclabs)
	
	# Set up a loop for the main offset
	for (i in 1:noLev1) {
		
		# Set up a loop for the secondary offset
		for (j in 1:noLev2) {
			
			# Plot the points and the confidence intervals
			xoff <- offvec1[i] + offvec2[j]
			n <- df[j,i]
			pest <- n/N
			ub <- srg.ci.binom(c(N),c(n),0.025)
			points(c(xoff),c(pest),col=colvec[j],pch=22,bg=colvec[j])
			lb <- srg.ci.binom(c(N),c(n),0.975)
			points(c(xoff,xoff),c(lb,ub),type="l",col=colvec[j])
			
		}
		
	}
	
	legend(0,max(ylim_in),legend=vecLev2,col=colvec,pt.bg=colvec,pch=22,bty="n")
	
	# Close the pdf device
	dev.off()
	
}



pdf.figure.proof <-
		function(findex=1,file=paste("./pdf_sub_",findex,".pdf",sep=""),pw=7,ph=7,textscale=0.6,xpd=NA) {
	plotwidth <- pw/cm(1)
	plotheight <- ph/cm(1)
	margin <- 5/cm(1)
	pdfwidth <- plotwidth+2*margin
	pdfheight <- plotheight+2*margin
	posplot = c(
			margin/pdfwidth,
			(margin+plotwidth)/pdfwidth,
			margin/pdfheight,
			(margin+plotheight)/pdfheight
	)
	pdf(file=file,height=pdfheight,width=pdfwidth,useDingbats=FALSE)
	par(
			cex=textscale,
			fig=posplot,
			mai=c(0,0,0,0),
			xpd=xpd)
}

srg.print.format <- function(x) {
  format(signif(x,3),scientific=FALSE,big.mark=",",drop0trailing=TRUE)
}


srg.simple.image.legend <- function(min=0,max=1,colpalette=heat.colors(10),logscale=FALSE,axlabs=NULL) {
	number=length(colpalette)
	imagedat <- matrix(nrow=number,ncol=2)
	xvals <- c(0,1,2)
	if (logscale) {
		min <- log10(min)
		max <- log10(max)
	}
	yvals <- (0:(number))/(number)*(max-min)+min
	for (i in 1:number) imagedat[i,] <- (yvals[i+1]+yvals[i])/2
	plot(1:2,type="n",axes=FALSE,xlim=c(0,2),ylim=c(min,max),xlab="",ylab="")
	image(xvals,yvals,t(imagedat),add=TRUE,col=colpalette)
	if (is.null(axlabs)) {
		axis(2,las=1)
	} else if (logscale==TRUE) {
		axis(2,las=1,at=log10(axlabs),las=1,labels=axlabs)
	} else {
		axis(2,las=1,at=axlabs,las=1,labels=axlabs)
	}
	
}
