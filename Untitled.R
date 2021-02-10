require(KernSmooth)

require(logspline)
set.seed(1)
require(ggplot2)
## 1
n=1000
data1=rbeta(n,0.5,0.5)

data2=0.8*rnorm(n,0,1)+0.2*rnorm(n,2,0.04)

data3=rnorm(n,0,1)

####======================================================================================
par(mfrow=c(1,3))
hist(data3, xlab=" ", freq = FALSE,xlim = c(-3,3), ylim=c(0,0.5), main = "Scenario 1: X~N(0,1)");
lines(density(data3,xlab=" ", ylab="", bw = "nrd0", kernel ="gaussian"),col="blue", lty=2);
lines(locpoly(data3,bandwidth = 0.25),col="red", lty=2)
par(new=T)
plot(logspline(data3),xlab=" ", ylab="",col="green", xlim = c(-3,3), ylim=c(0,0.5), lty=2)
par(new=T)
curve(dnorm(x,0,1),xlab=" ", ylab="",col="black",xlim = c(-3,3), ylim=c(0,0.5), lty=1)
letters <- c("KDE", "Local polynomial", "Logspline","True density")
legend("topright", legend = letters, lty = c(2,2,2,1), col = c("blue", "red", "green","black"), cex = 0.9)

hist(data1, xlab=" ", freq = FALSE,xlim = c(0,1), ylim=c(0,3.5), main = "Scenario 2: X~Beta(0.5,0.5)");
lines(density(data1,xlab=" ", ylab="", bw = "nrd0", kernel ="gaussian"),col="blue", lty=2);
lines(locpoly(data1,bandwidth = 0.25),col="red", lty=2)
par(new=T)
plot(logspline(data1),xlab=" ", ylab="",col="green",xlim = c(0,1), ylim=c(0,3.5), lty=2)
par(new=T)
curve(dbeta(x,0.5,0.5),xlab=" ", ylab="",col="black",xlim = c(0,1), ylim=c(0,3.5), lty=1)
letters <- c("KDE", "Local polynomial", "Logspline","True density")
legend("topright", legend = letters, lty = c(2,2,2,1), col = c("blue", "red", "green","black"), cex = 0.9)

hist(data2, xlab=" ", freq = FALSE,xlim = c(-2.4,3.3), ylim=c(0,0.5), main = "Scenario 3: X~0.8×N(0,1) +0.2×N(2,0.04)");
lines(density(data2,xlab=" ", ylab="", bw = "nrd0", kernel ="gaussian"),col="blue", lty=2);
lines(locpoly(data2,bandwidth = 0.25),col="red", lty=2)
par(new=T)
plot(logspline(data2),xlab=" ", ylab="",col="green", xlim = c(-2.4,3.3), ylim=c(0,0.5), lty=2)
par(new=T)
curve(0.8*dnorm(x,0,1)+0.2*dnorm(x,2,0.04),xlab=" ", ylab="",col="black",xlim = c(-2.4,3.3), ylim=c(0,0.5), lty=1)
letters <- c("KDE", "Local polynomial", "Logspline","True density")
legend("topleft", legend = letters, lty = c(2,2,2,1), col = c("blue", "red", "green","black"), cex = 0.9)

##########====================================================











### MC
ise1=matrix(NA, R,3)
ise2=matrix(NA, R,3)
ise3=matrix(NA, R,3)
ise4=matrix(NA, R,3)

## hist
n=1000
R=1000
for (r in 1:R){
  data1=0.8*rnorm(n,0,1)+0.2*rnorm(n,2,0.04)
  a=hist(data1, breaks=12, freq = FALSE)
  coun=a$counts
  densa=a$density
  k=length(a$counts)
  
  yhat=rep(densa[1],coun[1])
  for (i in 2:k){
    yhat=c(yhat,rep(densa[i],coun[i]))
  }
  data=sort(data1)
  ise1[r,3]=sfsmisc::integrate.xy(data,(yhat-0.8*dnorm(data,0,1)-0.2*dnorm(data,2,0.04))^2)
}

mean(ise1[,1])


## kde 
n=1000
R=1000
for (r in 1:R){
  x=0.8*rnorm(n,0,1)+0.2*rnorm(n,2,0.04)
  est=density(x)
  ise2[r,3]=sfsmisc::integrate.xy(est$x,(est$y-0.8*dnorm(est$x,0,1)-0.2*dnorm(est$x,2,0.04))^2)
}

mean(ise2[,1])

## spline
n=250
R=1000
for (r in 1:R){
  x=0.8*rnorm(n,0,1)+0.2*rnorm(n,2,0.04)
  est=logspline(x)
  y=dlogspline(x,est)
  ise3[r,1]=sfsmisc::integrate.xy(x,(y-0.8*dnorm(x,0,1)-0.2*dnorm(x,2,0.04))^2)
}

mean(ise3)
mean(ise3[,1])


## local poly
n=250
R=1000
for (r in 1:R){
  x=0.8*rnorm(n,0,1)+0.2*rnorm(n,2,0.04)
  est=locpoly(x, bandwidth = 0.25)
  ise4[r,1]=sfsmisc::integrate.xy(est$x,(est$y-0.8*dnorm(est$x,0,1)-0.2*dnorm(est$x,2,0.04))^2)
}

mean(ise4)
mean(ise4[,1])

####===============
par(mfrow=c(1,4))
boxplot(ise1,main="ISE for Histogram",outline=FALSE,
        ylab="", xlab="",names=c("n=250", "n=500","n=1000"))
boxplot(ise2,main="ISE for KDE",outline=FALSE,
        ylab="", xlab="",names=c("n=250", "n=500","n=1000"))
boxplot(ise3,main="ISE for logspline",outline=FALSE,
        ylab="", xlab="",names=c("n=250", "n=500","n=1000"))
boxplot(ise4,main="ISE for Local poly",outline=FALSE,
        ylab="", xlab="",names=c("n=250", "n=500","n=1000"))

##==================
x=runif(n,-4,4)
y=dnorm(x,0,1)
est=loess(x,y)


require(KernSmooth)

require(logspline)






