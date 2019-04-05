#################
# pareto. Produces a Pareto plot of effects.
# Parameters:
# effects - vector or matrix of effects to plot.
# names - vector of names to label the effects.
# xlab - String to display as the x axis label.
# ylab - String to display as the y axis label.
# perlab - Label for the cumulative percentage label.
# heading - Vector of names for plot heading.
#
pareto <- function(effects, names=NULL, xlab=NULL, ylab="Magnitude of Effect",
                   indicate.percent=TRUE, perlab="Cumulative Percentage", heading=NULL)
{
  # set up graphics parameters, note: set las=2 for perpendicular axis.
  oldpar <- par( mar=c(6, 4, 2, 4) + 0.1 , las=3)
  on.exit(par(oldpar))
  if( ! is.matrix(effects)) effects<-as.matrix( effects )
  for( i in 1:ncol(effects) )
  {
    if( i==2 ) oldpar$ask<-par(ask=TRUE)$ask
    # draw bar plot
    eff.ord <- rev(order(abs(effects[,i])))
    ef <- abs(effects[eff.ord,i])
    # plot barplot
    ylimit<-max(ef) + max(ef)*0.19
    ylimit<-c(0,ylimit)
    par( mar=c(6, 4, 2, 4) + 0.1 , las=3)
    x<-barplot(ef, names.arg=names[eff.ord], ylim=ylimit,
               xlab=xlab, ylab=ylab, main=heading[i])
    if( indicate.percent == TRUE ){
      # get cumulative sum of effects
      sumeff <- cumsum(ef)
      m<-max(ef)
      sm<-sum(ef)
      sumeff <- m * sumeff/sm
      # draws curve.
      lines(x, sumeff, lty="solid", lwd=2, col="purple")
      # draw 80% line
      lines( c(0,max(x)), rep(0.8*m,2) )
      # draw axis labling percentage.
      at <- c(0:5)* m/5
      axis(4, at=at,labels=c("0","20","40","60","80","100"), pos=max(x)+.6)
      # add axis lables
      par(las=0)
      mtext(perlab, 4, line=2)
    }
  } # end for each col
}
