set.seed(12)

x1 <- density(rnorm(100, mean=150, sd=30))

par(mar=c(5,5,1,1))

plot(x1, xlim=c(50,350), xlab="Environment", main=NA, axes=F,
     ylab="Density or Suitability")
axis(1, labels=NA)
axis(2, labels=NA)

polygon(x=c(52, 52, 352, 352), y=c(0, max(x1$y), max(x1$y), 0), 
        col='lightgrey', lty=2)

polygon(x1, col=rgb(0,0.2, 0.5, 0.7))

# segments(250,0,250,max(x1$y), lwd=3, lty=2)

# segments(x1$x[which(x1$y==max(x1$y))], 0, 
#          x1$x[which(x1$y==max(x1$y))], max(x1$y), 
#          lwd=3, lty=3)

# arrows(x1$x[which(x1$y==max(x1$y))]+2, quantile(x1$y, 0.9),
#        250-2, quantile(x1$y, 0.9), code=3, len=0.1, 
#        lwd=2, col='blue')

