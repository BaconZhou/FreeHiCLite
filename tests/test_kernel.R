x = cars$speed
y = cars$dist
z = rnorm(NROW(cars))

x <- c(0, 1, 2, 3, 4, 5)
y <- c(5, 4, 3, 2, 1, 0)

ksmooth(x, y)

dnorm(0) / dnorm(1)
dnorm(0) / dnorm(2)
41 / 26
41 / 7

c(0, rep(1, ))
dnorm(0:6) / sum(dnorm(0:6))

x <- c(0, 1, 1, 1, 1)
y <- rnorm(5)

band = 1
ksmooth(x, y, kernel = 'normal', bandwidth = band, x.points = c(1, 1))
w = dnorm(1 - x, sd = band / 4)
sum(w * y) / sum(w)



y = c(2.13195724, 0.        , 2.91833491, 1.45916746, 0.        ,
       1.45916746, 0.        , 0.        , 1.45916746, 0.        ,
       1.45916746, 0.        , 1.45916746, 1.45916746, 0.        ,
       0.        , 0.        , 1.45916746, 2.91833491, 5.83666982,
       1.45916746, 0.        , 2.91833491, 0.        , 0.        ,
       1.45916746, 4.37750237, 0.        , 0.        , 0.        ,
       2.91833491, 1.45916746, 1.45916746, 1.45916746, 0.        ,
       2.91833491, 0.        , 0.        , 0.        , 1.45916746)

ksmooth(x, y, kernel = 'normal', bandwidth = 20000)


data <- data.frame(Area = c(11,22,33,44,50,56,67,70,78,89,90,100),        RiverFlow = c(2337,2750,2301,2500,1700,2100,1100,1750,1000,1642, 2000,1932))                                 

x <- data$Area
y <- data$RiverFlow
#function to calculate Gaussian kernel
gausinKernel <- function(x,b){
  K <- (1/((sqrt(2*pi))))*exp(-0.5 *(x/b)^2)
  return(K)
}
b <- 10 #bandwidth
kdeEstimateyX <- x
ykernel <- NULL
for(xesti in kdeEstimateyX){
  xx <-  xesti - x
  K <-gausinKernel(xx,b)
  Ksum <- sum(K)
  weight <- K/Ksum
  yk <- sum(weight*y)
  xkyk <- c(xesti,yk)
  ykernel <- rbind(ykernel,xkyk)
}
ykernel
do.call('cbind', ksmooth(x, y, bandwidth = 10, x.points = x))
