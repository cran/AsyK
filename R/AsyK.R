#' Plot Density by RIG kernel.
#'
#' Plot Kernel density by using Resiprocal Inverse Gaussian Kernel.
#' @param y a numeric vector of positive values.
#' @param k gird points.
#' @param h the bandwidth
#' @import graphics
#' @import stats
#' @examples
#' y<-rexp(23,1)
#' h<-0.79 * IQR(y) * length(y) ^ (-1/5)
#' graphrig(y,80,h)
#' @export

graphrig<-function(y,k,h){
  n <- length(y)
  x <- seq(min(y) + 0.05, max(y), length = k)

  fhat <- rep(0, k)
  KRIG <- matrix(rep(0, k * n), ncol = k)

  for(j in 1:k) {
    for(i in 1:n) {

      KRIG[i, j] <- 1/(sqrt(2 * pi * h * y[i])) * exp(((-1 *(y[i] - h))/(2 * h) * (x[j]/(y[i] - h) - 2 +(y[i] - h)/x[j])))
    }
    fhat[j] <- 1/n * (sum(KRIG[, j]))
    d1<-density(y,bw=h)
    plot(x, fhat, type = "l",ylab = "Density Function", lty = 1, xlab = "Time",ylim=c(min(fhat),max(fhat)))
    lines(d1,type="p",col="red")
    legend("topright", c("Real Density", "Density by  RIG  Kernel"),
           col=c("red", "black"), lty=c(1,2))
  }}
#' Plot Density by Laplace kernel.
#'
#' Plot Kernel density by using Laplace Kernel.
#' @param y a numeric vector of positive values.
#' @param k gird points.
#' @param h the bandwidth
#' @import graphics
#' @import stats
#' @examples
#' y<-rexp(23,1)
#' h<-0.79 * IQR(y) * length(y) ^ (-1/5)
#' graphLap(y,80,h)
#' @export
#'
 graphLap<-function(y,k,h){
n <- length(y)
x <- seq(min(y) + 0.05, max(y), length = k)

fhat <- rep(0, k)
KLaplace <- matrix(rep(0, k * n), ncol = k)

for(j in 1:k) {
  for(i in 1:n) {

    KLaplace[i, j] <-(1/(2*sqrt(h)))*exp(-(abs(y[i]-(x[j])))/sqrt(h))
  }
  fhat[j] <- 1/n * (sum(KLaplace[, j]))
  d1<-density(y,bw=h)
  plot(x,fhat, type="l", ylab = "Density Function",ylim=c(min(fhat),(max(fhat)+0.2))  ,lty = 2, xlab = "Time")
  lines(d1,type="l",col="green")
  legend("topright", c("Real Density", "Density by Laplace Kernel"),
         col=c("green", "black"), lty=c(1,2))
}}

 #' Calculate Mean Squared Error( MSE) when RIG kernel is used.
 #'
 #' Calculate MSE by using Resiprocal Inverse Gaussian Kernel.
 #' @param y a numeric vector of positive values.
 #' @param k gird points.
 #' @param h the bandwidth
 #' @param type mention distribution of vector.If exponential distribution then use "Exp".
 #'     if use gamma distribution then use "Gamma".If Weibull distribution then use "Weibull".
 #' @import graphics
 #' @import stats
 #' @examples
 #' y<-rexp(100,1)
 #' h<-0.79 * IQR(y) * length(y) ^ (-1/5)
 #' mserig(y,200,h,"exp")
 #' @return MSE
 #' @export
 #'
 mserig<-function(y,k,h,type){
   n <- length(y)
   x <- seq(min(y) + 0.05, max(y), length = k)
   ftrue<-switch(type,

                 Exp = dexp(x,(1/mean(x))),
                 Weibull = dweibull(x, ((sd(x)/mean(x))^-1.806), scale = 1, log = FALSE),
                 Gamma = dgamma(x,(mean(x)/(var(x)/mean(x))),(var(x)/mean(x)))

   )

   fhat <- rep(0, k)
   KRIG <- matrix(rep(0, k * n), ncol = k)

   for(j in 1:k) {
     for(i in 1:n) {

       KRIG[i, j] <- 1/(sqrt(2 * pi * h * y[i])) * exp(((-1 *(y[i] - h))/(2 * h) * (x[j]/(y[i] - h) - 2 +(y[i] - h)/x[j])))
     }# 1st loop end
     fhat[j] <- 1/n * (sum(KRIG[, j]))

   }#2nd loop end
   return(mean((ftrue-fhat)^2))#mse of fhat w.r.t. the true density
 }#function end


 #' Calculate Mean Squared Error( MSE) when Laplace Kernel is used.
 #'
 #' Calculate MSE by using Laplace Kernel.
 #' @param y a numeric vector of positive values.
 #' @param k gird points.
 #' @param h the bandwidth
 #' @param type mention distribution of vector.If exponential distribution then use "Exp".
 #'     if use gamma distribution then use "Gamma".If Weibull distribution then use "Weibull".
 #' @import graphics
 #' @import stats
 #' @examples
 #' y<-rexp(100,1)
 #'  h<-0.79 * IQR(y) * length(y) ^ (-1/5)
 #'  mselap(y,200,h,"exp")
 #' @return MSE
 #' @export
 #'
 mselap<-function(y,k,h,type){
   n <- length(y)
   x <- seq(min(y) + 0.05, max(y), length = k)
   ftrue<-switch(type,

                 Exp = dexp(x,(1/mean(x))),
                 Weibull = dweibull(x, ((sd(x)/mean(x))^-1.806), scale = 1, log = FALSE),
                 Gamma = dgamma(x,(mean(x)/(var(x)/mean(x))),(var(x)/mean(x)))
   )

   fhat <- rep(0, k)
   KLaplace <- matrix(rep(0, k * n), ncol = k)

   for(j in 1:k) {
     for(i in 1:n) {

       KLaplace[i, j] <-(1/(2*sqrt(h)))*exp(-(abs(y[i]-(x[j])))/sqrt(h))
     }# 1st loop end
     fhat[j] <- 1/n * (sum(KLaplace[, j]))

   }#2nd loop end
   return(mean((ftrue-fhat)^2))#mse of fhat w.r.t. the true density
 }#function end

 #' Bandwidth Calculation.
 #'
 #' Calculate Bandwidth proposed by Silverman for nonnormal data.
 #' @param y a numeric vector of positive values.
 #' @import graphics
 #' @import stats
 #' @examples
 #' y<-rexp(10,1)
 #'  NSR(y)
 #' @return h
 #' @export
 NSR<-function(y){
   return(0.79 * IQR(y) * length(y) ^ (-1/5))
 }

 #' Calculate MSE with and ranking of Bandwidth with respect to MSE for RIG kernel.
 #'
 #' Caculate MSE with 19 bandwidths by using Resiprocal Inverse Gaussian Kernel.
 #' @param y a numeric vector of positive values.
 #' @param k gird points.
 #' @param type mention distribution of vector.If exponential distribution then use "Exp".
 #'     if use gamma distribution then use "Gamma".If Weibull distribution then use "Weibull".
 #' @import stats
 #' @import KernSmooth
 #' @import decon
 #' @importFrom  kedd h.bcv
 #' @importFrom kedd h.ccv
 #' @importFrom kedd h.mcv
 #' @importFrom kedd h.mlcv
 #' @importFrom kedd h.tcv
 #' @importFrom kedd h.ucv
 #' @import locfit
 #' @importFrom ks hscv
 #' @importFrom  kerdiest ALbw
 #' @importFrom  kerdiest PBbw
 #' @import sm
 #' @import ICV
 #' @import OSCV
 #' @examples
 #'  \donttest{y<-rexp(100,1)
 #'   righcomp(y, 200, "Exp")}
 #' @return MSE withh 19 bandwidths, Ranks, Minimum MSE, Maximum MSE
 #' @export
 righcomp<-function(y,k,type){
 n <- length(y)
 sd<-sd(y)
 x <- seq(min(y) + 0.05, max(y), length = k)

 PI <- abs(dpik(y) )                        ###plug-in
 AL<-abs(ALbw(vec_data=y))		###bandwidth of Altman and Leger
 bw2<-h.bcv(y)                         ###biased cross-validation.
 BCV<-abs(bw2$h)
 bw3<-h.ccv(y)                         ###complete cross-validation
 CCV<-abs(bw3$h)
 RoT <-abs(1.06*sd(y)*(n^(-1/5)))          ###normal scale rule
 bw8<-h.ucv(y)                         ###unbiased cross-validation
 UCV<-abs(bw8$h)
 bw10<-kdeb(y, h0 = 0.01 * sd, h1 = sd, meth = c("AIC", "LCV", "LSCV", "BCV","SJPI", "GKK"), kern = "gauss", gf = 2.5)
 AIC<-abs(bw10[1])			###Akaike information criterion
 GKK<-abs(bw10[6])			###Gasser, Kneip and K¨ohler (1991)
 ICV<-abs(h_ICV(y))			###Indirect cross- validation
 bw8<-cp(y,  sig2=1)			###MAllow
 mallow<-abs(bw8[5])
 bWr<-bw.dboot2(y,sd,h0='dboot1',error='normal',B=1000,grid=100,ub=2)###bootstap with sampling
 bw=gcv(y)
 gcv=abs(bw[4] )                            ###generalized cross-validation
 bw5<-h.mcv(y)                         ###modified cross-validation
 MCV<-abs(bw5$h)
 bw6<-h.mlcv(y)                        ###maximum likelihood cross-validation
 MLCV<-abs(bw6$h)
 bw7<-h.tcv(y)                         ###trimmed cross-validation
 TCV<-abs(bw7$h)
 SCV<-abs(hscv(y))		###smooth cross validation
 bWOr<-bw.dboot1(y,sd, h0="dnrd", error="normal", grid=100, ub=2)###bootstap without sampling
 pb<-PBbw( vec_data=y)		#PB
 oscv<-h_OSCV_dens(y,0)			#Onesided cross validation

 ftrue<-switch(type,

               Exp = dexp(x,(1/mean(x))),
               Gamma = dgamma(x,(mean(x)/(var(x)/mean(x))),(var(x)/mean(x))),
               Weibull = dweibull(x, ((sd(x)/mean(x))^-1.806), scale = 1, log = FALSE)
 )

 fhat1 <- rep(0, k)
 KRIG1 <- matrix(rep(0, k * n), ncol = k)
 fhat2 <- rep(0, k)
 KRIG2 <- matrix(rep(0, k * n), ncol = k)
 fhat3 <- rep(0, k)
 KRIG3 <- matrix(rep(0, k * n), ncol = k)
 fhat4 <- rep(0, k)
 KRIG4 <- matrix(rep(0, k * n), ncol = k)
 fhat5 <- rep(0, k)
 KRIG5 <- matrix(rep(0, k * n), ncol = k)
 fhat6 <- rep(0, k)
 KRIG6 <- matrix(rep(0, k * n), ncol = k)
 fhat7 <- rep(0, k)
 KRIG7 <- matrix(rep(0, k * n), ncol = k)
 fhat8 <- rep(0, k)
 KRIG8 <- matrix(rep(0, k * n), ncol = k)
 fhat9 <- rep(0, k)
 KRIG9 <- matrix(rep(0, k * n), ncol = k)
 fhat10 <- rep(0, k)
 KRIG10 <- matrix(rep(0, k * n), ncol = k)
 fhat11 <- rep(0, k)
 KRIG11 <- matrix(rep(0, k * n), ncol = k)
 fhat12 <- rep(0, k)
 KRIG12 <- matrix(rep(0, k * n), ncol = k)
 fhat13 <- rep(0, k)
 KRIG13 <- matrix(rep(0, k * n), ncol = k)
 fhat14 <- rep(0, k)
 KRIG14 <- matrix(rep(0, k * n), ncol = k)
 fhat15 <- rep(0, k)
 KRIG15 <- matrix(rep(0, k * n), ncol = k)
 fhat16 <- rep(0, k)
 KRIG16 <- matrix(rep(0, k * n), ncol = k)
 fhat17 <- rep(0, k)
 KRIG17 <- matrix(rep(0, k * n), ncol = k)
 fhat18 <- rep(0, k)
 KRIG18 <- matrix(rep(0, k * n), ncol = k)
 fhat19 <- rep(0, k)
 KRIG19 <- matrix(rep(0, k * n), ncol = k)

 for(j in 1:k) {
   for(i in 1:n) {

     KRIG1[i, j] <- 1/(sqrt(2 * pi * PI * y[i])) * exp(((-1 *(y[i] - PI))/(2 * PI) * (x[j]/(y[i] - PI) - 2 +(y[i] - PI)/x[j])))
     KRIG2[i, j] <- 1/(sqrt(2 * pi * AL * y[i])) * exp(((-1 *(y[i] - AL))/(2 *AL) * (x[j]/(y[i] -AL) - 2 +(y[i] - AL)/x[j])))
     KRIG3[i, j] <- 1/(sqrt(2 * pi * BCV * y[i])) * exp(((-1 *(y[i] - BCV))/(2 * BCV) * (x[j]/(y[i] - BCV) - 2 +(y[i] - BCV)/x[j])))
     KRIG4[i, j] <- 1/(sqrt(2 * pi * CCV * y[i])) * exp(((-1 *(y[i] - CCV))/(2 * CCV) * (x[j]/(y[i] - CCV) - 2 +(y[i] - CCV)/x[j])))
     KRIG5[i, j] <- 1/(sqrt(2 * pi * RoT * y[i])) * exp(((-1 *(y[i] -RoT))/(2 *RoT) * (x[j]/(y[i] - RoT) - 2 +(y[i] -RoT)/x[j])))
     KRIG6[i, j] <- 1/(sqrt(2 * pi * UCV * y[i])) * exp(((-1 *(y[i] - UCV))/(2 * UCV) * (x[j]/(y[i] - UCV) - 2 +(y[i] - UCV)/x[j])))
     KRIG7[i, j] <- 1/(sqrt(2 * pi * AIC * y[i])) * exp(((-1 *(y[i] - AIC))/(2 *AIC) * (x[j]/(y[i] - AIC) - 2 +(y[i] - AIC)/x[j])))
     KRIG8[i, j] <- 1/(sqrt(2 * pi * GKK * y[i])) * exp(((-1 *(y[i] - GKK))/(2 * GKK) * (x[j]/(y[i] - GKK) - 2 +(y[i] - GKK)/x[j])))
     KRIG9[i, j] <- 1/(sqrt(2 * pi * ICV * y[i])) * exp(((-1 *(y[i] - ICV))/(2 * ICV) * (x[j]/(y[i] - ICV) - 2 +(y[i] - ICV)/x[j])))
     KRIG10[i, j] <- 1/(sqrt(2 * pi * mallow * y[i])) * exp(((-1 *(y[i] - mallow))/(2 * mallow) * (x[j]/(y[i] - mallow) - 2 +(y[i] - mallow)/x[j])))
     KRIG11[i, j] <- 1/(sqrt(2 * pi * bWr * y[i])) * exp(((-1 *(y[i] -bWr))/(2 *bWr) * (x[j]/(y[i] - bWr) - 2 +(y[i] - bWr)/x[j])))
     KRIG12[i, j] <- 1/(sqrt(2 * pi * gcv * y[i])) * exp(((-1 *(y[i] - gcv))/(2 * gcv) * (x[j]/(y[i] - gcv) - 2 +(y[i] - gcv)/x[j])))
     KRIG13[i, j] <- 1/(sqrt(2 * pi * MCV * y[i])) * exp(((-1 *(y[i] - MCV))/(2 * MCV) * (x[j]/(y[i] - MCV) - 2 +(y[i] - MCV)/x[j])))
     KRIG14[i, j] <- 1/(sqrt(2 * pi * MLCV * y[i])) * exp(((-1 *(y[i] - MLCV))/(2 * MLCV) * (x[j]/(y[i] - MLCV) - 2 +(y[i] - MLCV)/x[j])))
     KRIG15[i, j] <- 1/(sqrt(2 * pi * TCV * y[i])) * exp(((-1 *(y[i] - TCV))/(2 * TCV) * (x[j]/(y[i] - TCV) - 2 +(y[i] - TCV)/x[j])))
     KRIG16[i, j] <- 1/(sqrt(2 * pi * SCV * y[i])) * exp(((-1 *(y[i] -SCV))/(2 * SCV) * (x[j]/(y[i] - SCV) - 2 +(y[i] - SCV)/x[j])))
     KRIG17[i, j] <- 1/(sqrt(2 * pi * bWOr * y[i])) * exp(((-1 *(y[i] - bWOr))/(2 * bWOr) * (x[j]/(y[i] - bWOr) - 2 +(y[i] - bWOr)/x[j])))
     KRIG18[i, j] <- 1/(sqrt(2 * pi * pb * y[i])) * exp(((-1 *(y[i] - pb))/(2 * pb) * (x[j]/(y[i] - pb) - 2 +(y[i] - pb)/x[j])))
     KRIG19[i, j] <- 1/(sqrt(2 * pi *oscv * y[i])) * exp(((-1 *(y[i] - oscv))/(2 * oscv) * (x[j]/(y[i] - oscv) - 2 +(y[i] - oscv)/x[j])))
   }# 1st loop end

   fhat1[j] <- 1/n * (sum(KRIG1[, j]))
   fhat2[j] <- 1/n * (sum(KRIG2[, j]))
   fhat3[j] <- 1/n * (sum(KRIG3[, j]))
   fhat4[j] <- 1/n * (sum(KRIG4[, j]))
   fhat5[j] <- 1/n * (sum(KRIG5[, j]))
   fhat6[j] <- 1/n * (sum(KRIG6[, j]))
   fhat7[j] <- 1/n * (sum(KRIG7[, j]))
   fhat8[j] <- 1/n * (sum(KRIG8[, j]))
   fhat9[j] <- 1/n * (sum(KRIG9[, j]))
   fhat10[j] <- 1/n * (sum(KRIG10[, j]))
   fhat11[j] <- 1/n * (sum(KRIG11[, j]))
   fhat12[j] <- 1/n * (sum(KRIG12[, j]))
   fhat13[j] <- 1/n * (sum(KRIG13[, j]))
   fhat14[j] <- 1/n * (sum(KRIG14[, j]))
   fhat15[j] <- 1/n * (sum(KRIG15[, j]))
   fhat16[j] <- 1/n * (sum(KRIG16[, j]))
   fhat17[j] <- 1/n * (sum(KRIG17[, j]))
   fhat18[j] <- 1/n * (sum(KRIG18[, j]))
   fhat19[j] <- 1/n * (sum(KRIG19[, j]))

 }#2nd loop end
 dpi<-(mean((ftrue-fhat1)^2))#mse of fhat w.r.t. the true density
 al<-(mean((ftrue-fhat2)^2))#mse of fhat w.r.t. the true density
 bcv<-(mean((ftrue-fhat3)^2))#mse of fhat w.r.t. the true density
 ccv<-(mean((ftrue-fhat4)^2))#mse of fhat w.r.t. the true density
 nsr<-(mean((ftrue-fhat5)^2))#mse of fhat w.r.t. the true density
 ucv<-(mean((ftrue-fhat6)^2))#mse of fhat w.r.t. the true density
 aic<-(mean((ftrue-fhat7)^2))#mse of fhat w.r.t. the true density
 gkk<-(mean((ftrue-fhat8)^2))#mse of fhat w.r.t. the true density
 icv<-(mean((ftrue-fhat9)^2))#mse of fhat w.r.t. the true density
 Mallow<-(mean((ftrue-fhat10)^2))#mse of fhat w.r.t. the true density
 bwr<-(mean((ftrue-fhat11)^2))#mse of fhat w.r.t. the true density
 Gcv<-(mean((ftrue-fhat12)^2))#mse of fhat w.r.t. the true density
 mcv<-(mean((ftrue-fhat13)^2))#mse of fhat w.r.t. the true density
 mlcv<-(mean((ftrue-fhat14)^2))#mse of fhat w.r.t. the true density
 tcv<-(mean((ftrue-fhat15)^2))#mse of fhat w.r.t. the true density
 scv<-(mean((ftrue-fhat16)^2))#mse of fhat w.r.t. the true density
 bwor<-(mean((ftrue-fhat17)^2))#mse of fhat w.r.t. the true density
 pB<-(mean((ftrue-fhat18)^2))#mse of fhat w.r.t. the true density
 osCV<-(mean((ftrue-fhat19)^2))#mse of fhat w.r.t. the true density

 v<-c("DPI"=dpi, "AL"=AL, "BCV"=bcv,"CCV"=ccv, "NSR"=nsr,"UCV"=ucv,"AIC"=aic,"GKK"=gkk,"AIC"=aic,"MallowCp"=Mallow,"bWr"=bwr,"GCV"=Gcv, "MCV"=mcv, "MLCV"=mlcv, "TCV"=tcv, "SCV"=scv,"bWOr"=bwor,"PB"=pB,"OSCV"=osCV)
 rank<-rank(v)
 mini<-v[which.min(v)]
 maxi<-v[which.max(v)]
 list<-list("All MSE" = v,  "Rank"=rank,"Minimum MSE" = mini, "Maximum MSE"=maxi)
 return(list)
 }#function end
 #' Calculate MSE with and ranking of Bandwidth with respect to MSE for Laplace Kernel.
 #'
 #' Caculate MSE with 19 bandwidths by using Laplace Kernel.
 #' @param y a numeric vector of positive values.
 #' @param k gird points.
 #' @param type mention distribution of vector.If exponential distribution then use "Exp".
 #'     if use gamma distribution then use "Gamma".If Weibull distribution then use "Weibull".
 #' @import stats
 #' @import KernSmooth
 #' @import decon
 #' @importFrom  kedd h.bcv
 #' @importFrom kedd h.ccv
 #' @importFrom kedd h.mcv
 #' @importFrom kedd h.mlcv
 #' @importFrom kedd h.tcv
 #' @importFrom kedd h.ucv
 #' @import locfit
 #' @importFrom ks hscv
 #' @importFrom  kerdiest ALbw
 #' @importFrom  kerdiest PBbw
 #' @import sm
 #' @import ICV
 #' @import OSCV
 #' @examples
 #'  \donttest{y<-rexp(100,1)
 #'   laphcomp(y, 200, "Exp")}
 #' @return MSE withh 19 bandwidths, Ranks, Minimum MSE, Maximum MSE
 #' @export
 laphcomp<-function(y,k,type){
   n <- length(y)
   sd<-sd(y)
   x <- seq(min(y) + 0.05, max(y), length = k)

   PI <- abs(dpik(y) )                        ###plug-in
   AL<-abs(ALbw(vec_data=y))		###bandwidth of Altman and Leger
   bw2<-h.bcv(y)                         ###biased cross-validation.
   BCV<-abs(bw2$h)
   bw3<-h.ccv(y)                         ###complete cross-validation
   CCV<-abs(bw3$h)
   RoT <-abs(1.06*sd(y)*(n^(-1/5)))          ###normal scale rule
   bw8<-h.ucv(y)                         ###unbiased cross-validation
   UCV<-abs(bw8$h)
   bw10<-kdeb(y, h0 = 0.01 * sd, h1 = sd, meth = c("AIC", "LCV", "LSCV", "BCV","SJPI", "GKK"), kern = "gauss", gf = 2.5)
   AIC<-abs(bw10[1])			###Akaike information criterion
   GKK<-abs(bw10[6])			###Gasser, Kneip and K¨ohler (1991)
   ICV<-abs(h_ICV(y))			###Indirect cross- validation
   bw8<-cp(y,  sig2=1)			###MAllow
   mallow<-abs(bw8[5])
   bWr<-bw.dboot2(y,sd,h0='dboot1',error='normal',B=1000,grid=100,ub=2)###bootstap with sampling
   bw=gcv(y)
   gcv=abs(bw[4] )                            ###generalized cross-validation
   bw5<-h.mcv(y)                         ###modified cross-validation
   MCV<-abs(bw5$h)
   bw6<-h.mlcv(y)                        ###maximum likelihood cross-validation
   MLCV<-abs(bw6$h)
   bw7<-h.tcv(y)                         ###trimmed cross-validation
   TCV<-abs(bw7$h)
   SCV<-abs(hscv(y))		###smooth cross validation
   bWOr<-bw.dboot1(y,sd, h0="dnrd", error="normal", grid=100, ub=2)###bootstap without sampling
   pb<-PBbw( vec_data=y)		#PB
   oscv<-h_OSCV_dens(y,0)			#Onesided cross validation

   ftrue<-switch(type,

                 Exp = dexp(x,(1/mean(x))),

                 Gamma = dgamma(x,(mean(x)/(var(x)/mean(x))),(var(x)/mean(x))),
                 Weibull = dweibull(x, ((sd(x)/mean(x))^-1.806), scale = 1, log = FALSE)
   )

   fhat1 <- rep(0, k)
   KLaplace1 <- matrix(rep(0, k * n), ncol = k)
   fhat2 <- rep(0, k)
   KLaplace2 <- matrix(rep(0, k * n), ncol = k)
   fhat3 <- rep(0, k)
   KLaplace3 <- matrix(rep(0, k * n), ncol = k)
   fhat4 <- rep(0, k)
   KLaplace4 <- matrix(rep(0, k * n), ncol = k)
   fhat5 <- rep(0, k)
   KLaplace5 <- matrix(rep(0, k * n), ncol = k)
   fhat6 <- rep(0, k)
   KLaplace6 <- matrix(rep(0, k * n), ncol = k)
   fhat7 <- rep(0, k)
   KLaplace7 <- matrix(rep(0, k * n), ncol = k)
   fhat8 <- rep(0, k)
   KLaplace8 <- matrix(rep(0, k * n), ncol = k)
   fhat9 <- rep(0, k)
   KLaplace9 <- matrix(rep(0, k * n), ncol = k)
   fhat10 <- rep(0, k)
   KLaplace10 <- matrix(rep(0, k * n), ncol = k)
   fhat11 <- rep(0, k)
   KLaplace11 <- matrix(rep(0, k * n), ncol = k)
   fhat12 <- rep(0, k)
   KLaplace12 <- matrix(rep(0, k * n), ncol = k)
   fhat13 <- rep(0, k)
   KLaplace13 <- matrix(rep(0, k * n), ncol = k)
   fhat14 <- rep(0, k)
   KLaplace14 <- matrix(rep(0, k * n), ncol = k)
   fhat15 <- rep(0, k)
   KLaplace15 <- matrix(rep(0, k * n), ncol = k)
   fhat16 <- rep(0, k)
   KLaplace16 <- matrix(rep(0, k * n), ncol = k)
   fhat17 <- rep(0, k)
   KLaplace17 <- matrix(rep(0, k * n), ncol = k)
   fhat18 <- rep(0, k)
   KLaplace18 <- matrix(rep(0, k * n), ncol = k)
   fhat19 <- rep(0, k)
   KLaplace19 <- matrix(rep(0, k * n), ncol = k)

   for(j in 1:k) {
     for(i in 1:n) {

       KLaplace1[i, j] <-(1/(2*sqrt(PI)))*exp(-(abs(y[i]-(x[j])))/sqrt(PI))
       KLaplace2[i, j] <-(1/(2*sqrt(AL)))*exp(-(abs(y[i]-(x[j])))/sqrt(AL))
       KLaplace3[i, j] <-(1/(2*sqrt(BCV)))*exp(-(abs(y[i]-(x[j])))/sqrt(BCV))
       KLaplace4[i, j] <-(1/(2*sqrt(CCV)))*exp(-(abs(y[i]-(x[j])))/sqrt(CCV))
       KLaplace5[i, j] <-(1/(2*sqrt(RoT)))*exp(-(abs(y[i]-(x[j])))/sqrt(RoT))
       KLaplace6[i, j] <-(1/(2*sqrt(UCV)))*exp(-(abs(y[i]-(x[j])))/sqrt(UCV))
       KLaplace7[i, j] <-(1/(2*sqrt(AIC)))*exp(-(abs(y[i]-(x[j])))/sqrt(AIC))
       KLaplace8[i, j] <-(1/(2*sqrt(GKK)))*exp(-(abs(y[i]-(x[j])))/sqrt(GKK))
       KLaplace9[i, j] <-(1/(2*sqrt(ICV)))*exp(-(abs(y[i]-(x[j])))/sqrt(ICV))
       KLaplace10[i, j] <-(1/(2*sqrt(mallow)))*exp(-(abs(y[i]-(x[j])))/sqrt(mallow))
       KLaplace11[i, j] <-(1/(2*sqrt(bWr)))*exp(-(abs(y[i]-(x[j])))/sqrt(bWr))
       KLaplace12[i, j] <-(1/(2*sqrt(gcv)))*exp(-(abs(y[i]-(x[j])))/sqrt(gcv))
       KLaplace13[i, j] <-(1/(2*sqrt(MCV)))*exp(-(abs(y[i]-(x[j])))/sqrt(MCV))
       KLaplace14[i, j] <-(1/(2*sqrt(MLCV)))*exp(-(abs(y[i]-(x[j])))/sqrt(MLCV))
       KLaplace15[i, j] <-(1/(2*sqrt(TCV)))*exp(-(abs(y[i]-(x[j])))/sqrt(TCV))
       KLaplace16[i, j] <-(1/(2*sqrt(SCV)))*exp(-(abs(y[i]-(x[j])))/sqrt(SCV))
       KLaplace17[i, j] <-(1/(2*sqrt(bWOr)))*exp(-(abs(y[i]-(x[j])))/sqrt(bWOr))
       KLaplace18[i, j] <-(1/(2*sqrt(pb)))*exp(-(abs(y[i]-(x[j])))/sqrt(pb))
       KLaplace19[i, j] <-(1/(2*sqrt(oscv)))*exp(-(abs(y[i]-(x[j])))/sqrt(oscv))
     }# 1st loop end

     fhat1[j] <- 1/n * (sum(KLaplace1[, j]))
     fhat2[j] <- 1/n * (sum(KLaplace2[, j]))
     fhat3[j] <- 1/n * (sum(KLaplace3[, j]))
     fhat4[j] <- 1/n * (sum(KLaplace4[, j]))
     fhat5[j] <- 1/n * (sum(KLaplace5[, j]))
     fhat6[j] <- 1/n * (sum(KLaplace6[, j]))
     fhat7[j] <- 1/n * (sum(KLaplace7[, j]))
     fhat8[j] <- 1/n * (sum(KLaplace8[, j]))
     fhat9[j] <- 1/n * (sum(KLaplace9[, j]))
     fhat10[j] <- 1/n * (sum(KLaplace10[, j]))
     fhat11[j] <- 1/n * (sum(KLaplace11[, j]))
     fhat12[j] <- 1/n * (sum(KLaplace12[, j]))
     fhat13[j] <- 1/n * (sum(KLaplace13[, j]))
     fhat14[j] <- 1/n * (sum(KLaplace14[, j]))
     fhat15[j] <- 1/n * (sum(KLaplace15[, j]))
     fhat16[j] <- 1/n * (sum(KLaplace16[, j]))
     fhat17[j] <- 1/n * (sum(KLaplace17[, j]))
     fhat18[j] <- 1/n * (sum(KLaplace18[, j]))
     fhat19[j] <- 1/n * (sum(KLaplace19[, j]))

   }#2nd loop end
   dpi<-(mean((ftrue-fhat1)^2))#mse of fhat w.r.t. the true density
   al<-(mean((ftrue-fhat2)^2))#mse of fhat w.r.t. the true density
   bcv<-(mean((ftrue-fhat3)^2))#mse of fhat w.r.t. the true density
   ccv<-(mean((ftrue-fhat4)^2))#mse of fhat w.r.t. the true density
   nsr<-(mean((ftrue-fhat5)^2))#mse of fhat w.r.t. the true density
   ucv<-(mean((ftrue-fhat6)^2))#mse of fhat w.r.t. the true density
   aic<-(mean((ftrue-fhat7)^2))#mse of fhat w.r.t. the true density
   gkk<-(mean((ftrue-fhat8)^2))#mse of fhat w.r.t. the true density
   icv<-(mean((ftrue-fhat9)^2))#mse of fhat w.r.t. the true density
   Mallow<-(mean((ftrue-fhat10)^2))#mse of fhat w.r.t. the true density
   bwr<-(mean((ftrue-fhat11)^2))#mse of fhat w.r.t. the true density
   Gcv<-(mean((ftrue-fhat12)^2))#mse of fhat w.r.t. the true density
   mcv<-(mean((ftrue-fhat13)^2))#mse of fhat w.r.t. the true density
   mlcv<-(mean((ftrue-fhat14)^2))#mse of fhat w.r.t. the true density
   tcv<-(mean((ftrue-fhat15)^2))#mse of fhat w.r.t. the true density
   scv<-(mean((ftrue-fhat16)^2))#mse of fhat w.r.t. the true density
   bwor<-(mean((ftrue-fhat17)^2))#mse of fhat w.r.t. the true density
   pB<-(mean((ftrue-fhat18)^2))#mse of fhat w.r.t. the true density
   osCV<-(mean((ftrue-fhat19)^2))#mse of fhat w.r.t. the true density

   v<-c("DPI"=dpi, "AL"=AL, "BCV"=bcv,"CCV"=ccv, "NSR"=nsr,"UCV"=ucv,"AIC"=aic,"GKK"=gkk,"AIC"=aic,"MallowCp"=Mallow,"bWr"=bwr,"GCV"=Gcv, "MCV"=mcv, "MLCV"=mlcv, "TCV"=tcv, "SCV"=scv,"bWOr"=bwor,"PB"=pB,"OSCV"=osCV)
   rank<-rank(v)
   mini<-v[which.min(v)]
   maxi<-v[which.max(v)]
   list<-list("All MSE" = v, "Rank"=rank,"Minimum MSE" = mini, "Maximum MSE"=maxi)
   return(list)
 }

 #' Plot Density by RIG and Laplace kernel.
 #'
 #' Plot densities by using Resiprocal Inverse Gaussian and Laplace Kernel.
 #' @param y a numeric vector of positive values.
 #' @param k gird points.
 #' @param h the bandwidth
 #' @import graphics
 #' @import stats
 #' @examples
 #' y<-rexp(100,1)
 #' h<-0.79 * IQR(y) * length(y) ^ (-1/5)
 #' compLR(y,80,h)
 #' @export
 compLR<-function(y,k,h){
   n <- length(y)
   x <- seq(min(y) + 0.05, max(y), length = k)

   fhat1 <- rep(0, k)
   fhat2 <- rep(0, k)
   KLaplace <- matrix(rep(0, k * n), ncol = k)
   KRIG <- matrix(rep(0, k * n), ncol = k)

   for(j in 1:k) {
     for(i in 1:n) {

       KLaplace[i, j] <-(1/(2*sqrt(h)))*exp(-(abs(y[i]-(x[j])))/sqrt(h))
       KRIG[i, j] <- 1/(sqrt(2 * pi * h * y[i])) * exp(((-1 *(y[i] - h))/(2 * h) * (x[j]/(y[i] - h) - 2 +(y[i] - h)/x[j])))

     }
     fhat1[j] <- 1/n * (sum(KLaplace[, j]))
     fhat2[j] <- 1/n * (sum(KRIG[, j]))
     plot(x,fhat1, type="l", ylab = "Density Function",  ,lty = 2, xlab = "Time")
     par(new=TRUE)
     plot( x, fhat2, type="l", col="red" ,ylab = "Density Function", lty = 1, xlab = "Time", axes = FALSE)
     legend("topright", c("Density by RIG Kernel", "Density by Laplace Kernel"),
            col=c("red", "black"), lty=c(1,2))
   }}
