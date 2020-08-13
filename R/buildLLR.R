# Copyright (C) 2020  Jochen Weile, Roth Lab
#
# This file is part of maveLLR
#
# maveLLR is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# maveLLR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with maveLLR  If not, see <https://www.gnu.org/licenses/>.

#' Build LLR function using kernel density estimation
#'
#' @param posScores fitness scores in the positive reference set
#' @param negScores fitness scores in the negative reference set
#' @param bw the bandwith to use for kernel density estimation (see kdensity package for more info).
#'   This can either be a numerical value or the name of the algorithm used to automatically choose one.
#' @param kernel the type of kernel to use (see kdensity package for more info)
#'
#' @return a list containing the LLR function in log10 scale and the positive and negative density functions
#' @export
#'
#' @examples
#' posScores <- rnorm(50,0.2,0.3)
#' negScores <- rnorm(30,0.8,0.2)
#' llr <- buildLLR.kernel(posScores,negScores,bw=0.1,kernel="gaussian")
#' drawDensityLLR(c(posScores,negScores),llr$llr,llr$posDens,llr$negDens,posScores,negScores)
#'
buildLLR.kernel <- function(posScores, negScores, bw=0.1, kernel="gaussian") {

  library(kdensity)

  posDens <- kdensity(posScores,bw=bw,kernel=kernel)
  negDens <- kdensity(negScores,bw=bw,kernel=kernel)

  llrFun <- function(score) log10(posDens(score)/negDens(score))

  return(list(llr=llrFun,posDens=posDens,negDens=negDens))
}


#' Build LLR function using Gaussian densities.
#'
#' @param posScores fitness scores in the positive reference set
#' @param negScores fitness scores in the negative reference set
#' @param spline logical; whether or not to apply spline monotonization. TRUE by default.
#'
#' @return a list containing the LLR function in log10 scale and the positive and negative density functions
#' @export
#'
#' @examples
#'
#' posScores <- rnorm(50,0.2,0.3)
#' negScores <- rnorm(30,0.8,0.2)
#' llr <- buildLLR.gauss(posScores,negScores)
#' drawDensityLLR(c(posScores,negScores),llr$llr,llr$posDens,llr$negDens,posScores,negScores)
#'
buildLLR.gauss <- function(posScores, negScores, spline=TRUE) {

  mpos <- mean(posScores)
  spos <- sd(posScores)
  mneg <- mean(negScores)
  sneg <- sd(negScores)

  posDens <- function(x) dnorm(x,mpos,spos)
  negDens <- function(x) dnorm(x,mneg,sneg)

  llrFun <- function(score) log10(posDens(score)/negDens(score))

  if (spline) {
    minPoint <- optimize(llrFun,interval=c(mpos, qnorm(0.999,mneg,sneg)),maximum=FALSE)
    if (minPoint$minimum < mneg) {
      minPoint$minimum <- Inf
    }
    maxPoint <- optimize(llrFun,interval=c(qnorm(0.001,mpos,spos),mneg),maximum=TRUE)
    if (maxPoint$maximum > mpos) {
      maxPoint$maximum <- -Inf
    }

    llrSpline <- function(scores) {
      sapply(scores,function(score) {
        if (is.na(score)) {
          NA
        } else if (score > minPoint$minimum) {
          minPoint$objective
        } else if (score < maxPoint$maximum) {
          maxPoint$objective
        } else {
          llrFun(score)
        }
      })
    }

    return(list(llr=llrSpline,posDens=posDens,negDens=negDens))

  } else {
    return(list(llr=llrFun,posDens=posDens,negDens=negDens))
  }

}


#' Draw a plot that summarizes the densities and LLR calculated by other functions
#'
#' @param scores the numerical scores in the Mave map
#' @param llrFun the LLR function
#' @param posDens the density function of the positive reference set
#' @param negDens the density function of the negative reference set
#' @param posScores the numerical scores in the postive reference set
#' @param negScores the numerical scores in the negative reference set
#'
#' @return nothing
#' @export
#'
#' @examples
#' posScores <- rnorm(50,0.2,0.3)
#' negScores <- rnorm(30,0.8,0.2)
#' llr <- buildLLR.gauss(posScores,negScores)
#' drawDensityLLR(c(posScores,negScores),llr$llr,llr$posDens,llr$negDens,posScores,negScores)
#'
drawDensityLLR <- function(scores, llrFun, posDens, negDens, posScores, negScores) {

  opar <- par(mfrow=c(2,1))
  xlim <- range(scores,na.rm=TRUE,finite=TRUE)
  ymax <- max(c(
    posDens(seq(xlim[[1]],xlim[[2]],length.out = 100)),
    negDens(seq(xlim[[1]],xlim[[2]],length.out = 100))
  ))
  par(mar=c(.1,4,1,1))
  plot(llrFun,from=xlim[[1]],to=xlim[[2]],xlim=xlim,axes=FALSE,xlab="",ylab="LLR")
  abline(h=0,col="gray",lty="dashed")
  axis(2)
  par(mar=c(5,4,.1,1))
  hist(scores,col="gray90",border=NA,freq=FALSE,main="",xlim=xlim,ylim=c(0,ymax))
  plot.function(posDens,from=xlim[[1]],to=xlim[[2]],add=TRUE,col="firebrick3",lwd=2)
  plot.function(negDens,from=xlim[[1]],to=xlim[[2]],add=TRUE,col="darkolivegreen3",lwd=2)
  abline(v=posScores,col="firebrick3")
  abline(v=negScores,col="darkolivegreen3")
  par(opar)

  return(invisible(NULL))
}

#
# mthfr <- read.csv("~/projects/mthfr/folate_response_model5.csv")
# variants <- mthfr[mthfr$type=="substitution","hgvs"]
# scores <- mthfr[mthfr$type=="substitution","m25.score"]
#
# idb <- 337
# prsnrs <- read.csv("~/projects/mthfr/mthfr_prsNrs.csv")
#
# posref <- prsnrs[prsnrs$reference=="positive" & prsnrs$start < idb,"hgvs"]
# negref <- prsnrs[prsnrs$reference=="negative" & prsnrs$start < idb,"hgvs"]
# posScores <- na.omit(scores[variants %in% posref])
# negScores <- na.omit(scores[variants %in% negref])
# llr <- buildLLR.kernel(posScores, negScores,bw=0.1,kernel="gaussian")
# drawDensityLLR(scores,llr$llr,llr$posDens,llr$negDens,posScores,negScores)
# llr <- buildLLR.gauss(posScores, negScores)
# drawDensityLLR(scores,llr$llr,llr$posDens,llr$negDens,posScores,negScores)
#
# posref <- prsnrs[prsnrs$reference=="positive" & prsnrs$start > idb,"hgvs"]
# negref <- prsnrs[prsnrs$reference=="negative" & prsnrs$start > idb,"hgvs"]
# posScores <- na.omit(scores[variants %in% posref])
# negScores <- na.omit(scores[variants %in% negref])
# llr <- buildLLR.kernel(posScores, negScores)
# drawDensityLLR(scores,llr$llr,llr$posDens,llr$negDens,posScores,negScores)
# llr <- buildLLR.gauss(posScores, negScores)
# drawDensityLLR(scores,llr$llr,llr$posDens,llr$negDens,posScores,negScores)
