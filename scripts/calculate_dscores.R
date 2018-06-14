######################################################
# To Run: 
# Rscript calculate_dscores.R predicted.txt pred_min svcounts.txt outputfile
#

######################################################
# Read in command line arguments
args <- commandArgs(trailingOnly=TRUE)
predfile= args[1]
pred_min=args[2]
svfile= args[3]
output = args[4]

######################################################
# Load appropriate libraries for doing analysis 
library(VGAM)
library(fdrtool)
library(fitdistrplus) 
library(actuar)

find_best_dist_sv<- function(dat, myxmin=0, myxlim=NULL){
    realcnt <- dat[dat >= myxmin]
    if (is.null(myxlim)){
        qnt <- quantile(realcnt, probs=c(0.25,0.5,0.75))
        iqr <- qnt[3]-qnt[1]
        myxlim <- qnt[2] + 6 *iqr
    }
    mycnts <- realcnt[realcnt < myxlim]
    fl <- fitdist(mycnts, "lnorm", method="mle")
    fll <- fitdist(mycnts, "llogis", method="mle")
    fg <- fitdist(mycnts, "gamma", method="mle")
    myfits = list(fl, fll, fg)
    plot.legend=c("lognormal", "loglogistic", "gamma")
    stats=gofstat(myfits, fitnames=plot.legend)
    # get the pvalues of the data for the distribution with the best fit (using BIC measure)
    bestd <- which.min(stats$bic)
    if (bestd==1){
        pvals <- plnorm(dat, meanlog=fl$estimate[1], sdlog=fl$estimate[2], log.p=TRUE, lower.tail=FALSE)}
    else if (bestd==2){
        pvals <- pllogis(dat, shape=fll$estimate[1], scale=fll$estimate[2], log.p=TRUE, lower.tail=FALSE)}
    else if (bestd==3){
        pvals <- pgamma(dat, shape=fg$estimate[1], rate=fg$estimate[2], log.p=TRUE, lower.tail=FALSE)}
    return(pvals)
}

find_best_dist_dsb<- function(dat, myxmin=0, myxlim=NULL){
    realcnt <- dat[dat >= myxmin]
    if (is.null(myxlim)){
        qnt <- quantile(realcnt, probs=c(0.25,0.5,0.75))
        iqr <- qnt[3]-qnt[1]
        myxlim <- qnt[2] + 6 *iqr
    }
    mycnts <- realcnt[realcnt < myxlim]
    fl <- fitdist(mycnts, "lnorm", method="mle")
    fll <- fitdist(mycnts, "llogis", method="mle")
    fg <- fitdist(mycnts, "gamma", method="mle")
    myfits = list(fl, fll, fg)
    plot.legend=c("lognormal", "loglogistic", "gamma")
    stats=gofstat(myfits, fitnames=plot.legend)
    # get the pvalues of the data for the distribution with the best fit (using BIC measure)
    bestd <- which.min(stats$bic)
    if (bestd==1){
        pvals <- plnorm(dat, meanlog=fl$estimate[1], sdlog=fl$estimate[2], log.p=TRUE, lower.tail=FALSE)}
    else if (bestd==2){
        pvals <- pllogis(dat, shape=fll$estimate[1], scale=fll$estimate[2], log.p=TRUE, lower.tail=FALSE)}
    else if (bestd==3){
        pvals <- pgamma(dat, shape=fg$estimate[1], rate=fg$estimate[2], log.p=TRUE, lower.tail=FALSE)}
    return(pvals)
}

find_best_dist_dscores <- function(mydat, cntmin=1, dsbmin=0, myxmin=NULL, myxlim=NULL){
    dat=mydat$dscore
    realcnt=dat[mydat$realcnt >= cntmin & mydat$predicted >= dsbmin]
    if (is.null(myxlim)){
        qnt <- quantile(realcnt, probs=c(0.25,0.5,0.75))
        iqr <- qnt[3]-qnt[1]
        myxlim <- qnt[2] + 6 * iqr
        myxmin <- qnt[2] - 6 * iqr
    }
    mycnts <- realcnt[realcnt < myxlim & realcnt > myxmin]
    ft <- fitdist(mycnts, "t", start=list(df=4), method="mle")
    fn <- fitdist(mycnts, "norm", method="mle")
    fc <- fitdist(mycnts, "cauchy", method="mle")
    myfits = list(ft, fn, fc)
    plot.legend=c("student-t", "normal", "Cauchy")
    stats=gofstat(myfits, fitnames=plot.legend)
    # get the pvalues of the data for the distribution with the best fit (using BIC measure)
    bestd <- which.min(stats$bic)
    if (bestd==1){
        pvals <- pt(dat, df=ft$estimate[1], log.p=TRUE, lower.tail=FALSE)}
    else if (bestd==2){
        pvals <- pnorm(dat, mean=fn$estimate[1], sd=fn$estimate[2], log.p=TRUE, lower.tail=FALSE)}
    else if (bestd==3){
        pvals <- pcauchy(dat, location=fc$estimate[1], scale=fc$estimate[2], log.p=TRUE, lower.tail=FALSE)}
    return(pvals)
}

make_dscores <- function(pdat, svdat, cntmin=1, dsbmin=0.5){
    ddat=cbind(svdat[,1:3], pdat[,4:6], svdat[,4:5])
    ddat$dscore = pdat$predlogp - svdat$log.p
    # get dscore pvalues 
    pv <- find_best_dist_dscores(mydat=ddat, cntmin=cntmin, dsbmin=dsbmin)
    ddat$dlogp = pv
    return(ddat)
}


preddat=read.table(predfile, header=TRUE)
pv <- find_best_dist_dsb(dat=preddat$predicted, myxmin=pred_min, myxlim=NULL)
preddat$predlogp = pv
dim(preddat)

svdat= read.table(svfile, header=TRUE)
pv <- find_best_dist_sv(dat=svdat$realcnt, myxmin=1, myxlim=NULL)
svdat$log.p = pv
dim(svdat)

ddat <- make_dscores(preddat, svdat, cntmin=1, dsbmin=pred_min)
write.table(ddat, file=output, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)


