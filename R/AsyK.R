#' AsyK
#'
#' Kernel Density Estimation and Selection of Optimum Bandwidth
#'
#' @description  A collection of functions related to density estimation by using Chen's (2000) idea. For observing estimated values see \code{\link{Laplace}} and \code{\link{RIG}}. Plots by using these kernels can be drawn by \code{\link{plot.Laplace}} and \code{\link{plot.RIG}}. Additionally, their combined plot is drawn by using \code{\link{compLR}}. Mean squared errors (MSE) can be calculated by \code{\link{mseLap}} and \code{\link{mseRIG}}. Further \code{\link{laphcomp}} and \code{\link{righcomp}} allows to calculate MSE by using
#' 19 different bandwidths for both kernels. Here we also present a normal scale rule bandwidth which is given by Silverman (1986) for nonnormal data.
#'@author Javaria Ahmad Khan, Atif Akbar.
#'@references \itemize{
#' \item Chen, S. X. 2000. Probability density function estimation using Gamma kernels. \emph{Annals of the Institute of Statistical Mathematics} \strong{52} (3), 471-480.
#' \item Silverman, B. W. 1986. \emph{Density Estimation}. Chapman & Hall/ CRC, London.
#' \item Sheather, S.J.; Jones, M.C. 1991. A Reliable Data-Based Bandwidth Selection Method for Kernel Density Estimation. \emph{Journal of the Royal Statistical Society, B}, \strong{53}, 683-690.
#' \item Bowman, A.W.; Hall, P.; Prvan, T. 1998. Cross-validation for The Smoothing of Distribution Functions. \emph{Biometrika}, \strong{85}, 799-808.
#' \item Rudemo, M. 1982. Empirical Choice of Histograms and Kernel Density Estimators. \emph{Scandinavian Journal of Statistics}, \strong{9}, 65–78.
#' \item Feluch, W.; Koronacki, J. 1992. A Note on Modified Cross-Validation in Density Estimation. \emph{Computational Statistics and Data Analysis}, \strong{13}, 143–151.
#' \item Stute, W. 1992. Modified Cross Validation in Density Estimation. \emph{Journal of Statistical Planning and Inference}, \strong{30}, 293–305.
#' \item Miiller, H.-G. 1985. Empirical Bandwidth Choice for Nonparametric Kernel Regression by Means of Pilot Estimators. \emph{Statistics and Decisions}, Supplement No.\strong{2}, 193-206.
#' \item Hall, P.; Marron, J.S.; Park, B.U. 1992. Smoothed Cross-Validation. \emph{Probability Theory and Related Fields}, 92}, 1-20.
#' \item Habbema, J. D. F.; Hermans, J.; Van den Broek, K. 1974. \emph{A Discrimination Analysis Program Using Density Estimation. COMPSTAT 1974: Proceedings in Computational Statistics}. Physica Verlag, Vienna.
#' \item Wang, X.F.; Wang, B. 2011. Deconvolution Estimation in Measurement Error Models: The R package decon. \emph{Journal of Statistical Software}, \strong{39} (10),1- 24.
#' \item Altman, N.; Leger, C. 1995. Bandwidth Selection for Kernel Distribution Function Estimation. \emph{Journal of Statistical Planning and Inference}, \strong{46}, 195–214.
#' \item Jones, M. C.; Kappenman, R. F. 1991. On A Class of Kernel Density Estimate Bandwidth Selectors. \emph{Scandinavian Journal of Statistics}, \strong{19},337–349.
#' \item Scott, D. W., Terrel, G, R. 1987. Biased and Unbiased Cross-Validation in Density Estimation. \emph{Journal of the American Statistical Association}, \strong{82}, 1131–1146.
#' \item Savchuk Y. O.; Jeffrey D.; Hart & Simon P. Sheather 2013. One-sided Cross-Validation for Non-smooth Regression Functions. \emph{Journal of Nonparametric Statistics}, \strong{25} (4), 889-904.
#' \item Savchuk Y. O.; Jeffrey D.; Hart & Simon J. Sheather 2010. Indirect Cross Validation for Density Estimation. \emph{Journal of the American Statistical Association}, \strong{105} (489), 415-423.
#' \item Polansky, A. M.; Baker, E. R. 2000. Multistage plug-in bandwidth selection for kernel distribution function estimates. \emph{Journal of Statistical Computation and Simulation}, \strong{65}, 63–80.
#' \item Akaike, H. 1970. Statistical Predictor Identification. \emph{Annals of the Institute of Statistical Mathematics}, \strong{22}, 203-217.
#' \item Mallows, C. 1973. Some Comments Cp. \emph{Technometrics}, \strong{15}, 661-675.
#' \item Scaillet, O. 2004. Density estimation using inverse and reciprocal inverse Gaussian kernels. \emph{Nonparametric Statistics}, \strong{16}, 217-226.
#' \item Carven, P.; Wahba, G. 1979. Smoothing Noisy Data with Spline Functions. \emph{Numerische Ifathemat Ik}, \strong{31}, 377-403.
#' \item Staniswalis, J.G. 1989a. Local bandwidth selection for kernel estimates. \emph{Journal of the American Statistical Association}, \strong{84}, 284-8.
#' \item Gasser, T.; Kneip, A.; K¨ohler, W.  1991. A Flexible and Fast Method for Automatic Smoothing. \emph{Journal of the American Statistical Association}, \strong{86}, 643–652.
#' \item Khan, J. A.; Akbar, A. Density Estimation by Laplace Kernel. \emph{Working paper,  Department of Statistics, Bahauddin Zakariya University, Multan, Pakistan.}
#' }
"_PACKAGE"
#' Estimated Density Values by Resiprocal Inverse Gaussian kernel
#'
#' Estimated Kernel density values by using Resiprocal Inverse Gaussian Kernel.
#' @details Scaillet 2003. proposed Resiprocal Inverse Gaussian kerenl. He claimed that his proposed kernel share the same properties as those of gamma kernel estimator.
#' \deqn{K_{RIG \left( \ln{ax}4\ln {(\frac{1}{h})} \right)}(y)=\frac{1}{\sqrt {2\pi y}}  exp\left[-\frac{x-h}{2h} \left(\frac{y}{x-h}-2+\frac{x-h}{y}\right)\right]}
#' @param y a numeric vector of positive values.
#' @param k gird points.
#' @param h the bandwidth
#' @import stats
#' @examples
#' y <- rexp(100,1)
#' h <- 0.79 * IQR(y) * length(y) ^ (-1/5)
#' RIG(y,200,h)
#' @return \item{x}{grid points}
#'         \item{y}{estimated values of density}
#' @author Javaria Ahmad Khan, Atif Akbar.
#' @references Scaillet, O. 2004. Density estimation using inverse and reciprocal inverse Gaussian kernels. 	 \emph{Nonparametric Statistics}, \strong{16}, 217-226.
#' @seealso To examine RIG density plot see \code{\link{plot.RIG}} and for Mean Squared Error \code{\link{mseRIG}}. Similarly, for Laplace kernel \code{\link{Laplace}}.
#' @export
RIG<-function(y,k,h){
   n <- length(y)
   x <- seq(min(y) + 0.05, max(y), length = k)

   fhat <- rep(0, k)
   KRIG <- matrix(rep(0, k * n), ncol = k)

   for(j in 1:k) {
      for(i in 1:n) {

         KRIG[i, j] <- 1/(sqrt(2 * pi * h * y[i])) * exp(((-1 *(y[i] - h))/(2 * h) * (x[j]/(y[i] - h) - 2 +(y[i] - h)/x[j])))
      }
      fhat[j] <- 1/n * (sum(KRIG[, j]))
   }
   results <- list(x=x, y=fhat)
   class ( results ) <-c('list', 'RIG')
   results
}

#' Density Plot by Resiprocal Inverse Gaussian kernel
#'
#' Plot density by using Resiprocal Inverse Gaussian Kernel.
#' @param x an object of class "RIG"
#' @param \dots Not presently used in this implementation
#' @import graphics
#' @import stats
#' @examples
#' y <- rexp(100, 1)
#' h <- 0.79 * IQR(y) * length(y) ^ (-1/5)
#' den<-RIG(y,200,h)
#' plot(den, type = "s", ylab = "Density Function", lty = 1, xlab = "Time")
#' d1 <- density(y, bw=h) #To add true density along with estimated
#' lines(d1,type="p",col="red")
#' legend("topright", c("Real Density", "Density by RIG Kernel"), col=c("red", "black"), lty=c(1,2))
#' @return nothing
#' @author Javaria Ahmad Khan, Atif Akbar.
#' @references Scaillet, O. 2004. Density estimation using inverse and reciprocal inverse Gaussian kernels. 	 \emph{Nonparametric Statistics}, \strong{16}, 217-226.
#' @seealso To examine RIG estimated values for density see \code{\link{RIG}} and for Mean Squared Error \code{\link{mseRIG}}. Similarly, for plot of Laplace kernel \code{\link{plot.Laplace}}.
#' @export
plot.RIG <- function(x,...) {
   plot(x$x, x$y,...)
}
#' Estimated Density Values by Laplace kernel
#'
#' Estimated Kernel density values by using Laplace Kernel.
#' @details Laplace kernel is developed by Khan and Akbar. Kernel is developed by using Chen's idea. Laplace kernel is;
#' \deqn{K_{Laplace\left(x,h^{\frac{1}{2}}\right)} (u)=\frac{1}{2\sqrt h}exp \left(-\frac{|t{u-x}|}{\sqrt h}\right)}

#' @param y a numeric vector of positive values.
#' @param k gird points.
#' @param h the bandwidth
#' @import stats
#' @examples
#' y <- rexp(100,1)
#' h <- 0.79 * IQR(y) * length(y) ^ (-1/5)
#' Laplace(y,200,h)
#' @return \item{x}{grid points}
#'         \item{y}{estimated values of density}
#' @author Javaria Ahmad Khan, Atif Akbar.
#' @references Khan, J. A.; Akbar, A. Density Estimation by Laplace Kernel. \emph{Working paper,  Department of Statistics, Bahauddin Zakariya University, Multan, Pakistan.}
#' @seealso To examine laplace density plot see \code{\link{plot.Laplace}} and for Mean Squared Error \code{\link{mseLap}}. Similarly, for RIG kernel \code{\link{RIG}}.
#' @export
Laplace<-function(y,k,h){
   n <- length(y)
   x <- seq(min(y) + 0.05, max(y), length = k)

   fhat <- rep(0, k)
   KLaplace <- matrix(rep(0, k * n), ncol = k)

   for(j in 1:k) {
      for(i in 1:n) {

         KLaplace[i, j] <-(1/(2*sqrt(h)))*exp(-(abs(y[i]-(x[j])))/sqrt(h))
      }
      fhat[j] <- 1/n * (sum(KLaplace[, j]))
   }
   results <- list(x=x, y=fhat)
   class ( results ) <-c('list', 'Laplace')
   results
}

#' Density Plot by Laplace kernel
#'
#' Plot density by using Laplace Kernel.
#' @param x an object of class "Laplace"
#' @param \dots Not presently used in this implementation
#' @import graphics
#' @import stats
#' @examples
#' y <- rexp(100,1)
#' h <- 0.79 * IQR(y) * length(y) ^ (-1/5)
#' den <- Laplace(y, 200, h)
#' plot(den, type = "s", ylab = "Density Function", lty = 1, xlab = "Time")
#' d1 <- density(y, bw=h)
#' lines(d1,type="p",col="red")
#' legend("topright", c("Real Density", "Density by Laplace Kernel"), col=c("red", "black"))
#' @return nothing
#' @author Javaria Ahmad Khan, Atif Akbar.
#' @references Khan, J. A.; Akbar, A. Density Estimation by Laplace Kernel. \emph{Working paper,  Department of Statistics, Bahauddin Zakariya University, Multan, Pakistan.}
#' @seealso To examine Laplace estimated values for density see \code{\link{Laplace}} and for Mean Squared Error \code{\link{mseLap}}. Similarly, for plot of Laplace kernel \code{\link{plot.RIG}}.

#' @export
plot.Laplace <- function(x,...) {
   plot(x$x, x$y,...)
}

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
 #' @references Scaillet, O. 2004. Density estimation using inverse and reciprocal inverse Gaussian kernels. 	 \emph{Nonparametric Statistics}, \strong{16}, 217-226.
 #' @seealso For further MSE by using Laplace kernel see \code{\link{mseLap}}. For density estimation by using RIG Kernel \code{\link{plot.RIG}} and for estimated values
 #' of density \code{\link{RIG}}.
 #' @examples
 #' y<-rexp(100,1)
 #' h<-0.79 * IQR(y) * length(y) ^ (-1/5)
 #' mseRIG(y,200,h,"Exp")
 #' @return MSE
 #' @export
 #'
 mseRIG<-function(y,k,h,type){
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
     }
     fhat[j] <- 1/n * (sum(KRIG[, j]))

   }
   return(mean((ftrue-fhat)^2))#mse of fhat w.r.t. the true density
 }


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
 #' @references Khan, J. A.; Akbar, A. Density Estimation by Laplace Kernel. \emph{Working paper,  Department of Statistics, Bahauddin Zakariya University, Multan, Pakistan.}
 #' @seealso For further MSE by using RIG kernel see \code{\link{mseRIG}}. For density estimation by using Laplace Kernel \code{\link{plot.Laplace}} and for estimated values
 #' of density \code{\link{Laplace}}.
 #' @examples
 #' y<-rexp(100,1)
 #'  h<-0.79 * IQR(y) * length(y) ^ (-1/5)
 #'  mseLap(y,200,h,"Exp")
 #' @return MSE
 #' @export
 #'
 mseLap<-function(y,k,h,type){
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
     }
     fhat[j] <- 1/n * (sum(KLaplace[, j]))

   }
   return(mean((ftrue-fhat)^2))
 }

 #' Bandwidth Calculation.
 #'
 #' Calculate Bandwidth proposed by Silverman for nonnormal data.
 #' @param y a numeric vector of positive values.
 #' @import graphics
 #' @import stats
 #' @author Javaria Ahmad Khan, Atif Akbar.
 #' @references Silverman, B. W. 1986. \emph{Density Estimation}. Chapman & Hall/ CRC, London.
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
 #' @details This function helps to calculate MSE by using 19 different bandwidths which are Normal Sacale Rule (NSR), Complete Cross Validation (CCV), Biased Cross Validation (BCV), Unbiased Cross Validation (UBCV),
 #' Direct Plug-In (DPI), Modified Cross Validation (MCV), Maximum Likelihood Cross Validation (MLCV), Trimmed Cross Validation (TCV),Smooth Cross Validation (SCV), Bootstrap without Sampling (bWOs), Bootstrap with Sampling (bWs),
 #'  Bandwidth of Altman and Leger (AL),One-sided Cross Validation (OCV), Akaike information criterion (AIC),Indirect Cross Validation (ICV), Mallow’ Cp (MallowCp), Generalized Cross Validation (GCV), Polansky and Baker Plug-In (PB),
 #'  and Gasser, Kniep, and Köhler Cross Validation (GKK). For RIG kernel see \code{\link{laphcomp}}
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
 #' @references Scaillet, O. 2004. Density estimation using inverse and reciprocal inverse Gaussian kernels. 	 \emph{Nonparametric Statistics}, \strong{16}, 217-226.
 #' @author Javaria Ahmad Khan, Atif Akbar.
 #' @export
 righcomp<-function(y,k,type){
 n <- length(y)
 sd<-sd(y)
 x <- seq(min(y) + 0.05, max(y), length = k)

 PI <- abs(dpik(y) )
 AL<-abs(ALbw(vec_data=y))
 bw2<-h.bcv(y)
 BCV<-abs(bw2$h)
 bw3<-h.ccv(y)
 CCV<-abs(bw3$h)
 RoT <-abs(1.06*sd(y)*(n^(-1/5)))
 bw8<-h.ucv(y)
 UCV<-abs(bw8$h)
 bw10<-kdeb(y, h0 = 0.01 * sd, h1 = sd, meth = c("AIC", "LCV", "LSCV", "BCV","SJPI", "GKK"), kern = "gauss", gf = 2.5)
 AIC<-abs(bw10[1])
 GKK<-abs(bw10[6])
 ICV<-abs(h_ICV(y))
 bw8<-cp(y,  sig2=1)
 mallow<-abs(bw8[5])
 bWr<-bw.dboot2(y,sd,h0='dboot1',error='normal',B=1000,grid=100,ub=2)
 bw=gcv(y)
 gcv=abs(bw[4] )
 bw5<-h.mcv(y)
 MCV<-abs(bw5$h)
 bw6<-h.mlcv(y)
 MLCV<-abs(bw6$h)
 bw7<-h.tcv(y)
 TCV<-abs(bw7$h)
 SCV<-abs(hscv(y))
 bWOr<-bw.dboot1(y,sd, h0="dnrd", error="normal", grid=100, ub=2)
 pb<-PBbw( vec_data=y)
 oscv<-h_OSCV_dens(y,0)

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
   }

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

 }
 dpi<-(mean((ftrue-fhat1)^2))
 al<-(mean((ftrue-fhat2)^2))
 bcv<-(mean((ftrue-fhat3)^2))
 ccv<-(mean((ftrue-fhat4)^2))
 nsr<-(mean((ftrue-fhat5)^2))
 ucv<-(mean((ftrue-fhat6)^2))
 aic<-(mean((ftrue-fhat7)^2))
 gkk<-(mean((ftrue-fhat8)^2))
 icv<-(mean((ftrue-fhat9)^2))
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
 }
 #' Calculate MSE with and ranking of Bandwidth with respect to MSE for Laplace Kernel.
 #'
 #' Caculate MSE with 19 bandwidths by using Laplace Kernel.
 #' @details This function helps to calculate MSE by using 19 different bandwidths which are Normal Sacale Rule (NSR), Complete Cross Validation (CCV), Biased Cross Validation (BCV), Unbiased Cross Validation (UBCV),
 #' Direct Plug-In (DPI), Modified Cross Validation (MCV), Maximum Likelihood Cross Validation (MLCV), Trimmed Cross Validation (TCV), Smooth Cross Validation (SCV), Bootstrap without Sampling (bWOs), Bootstrap with Sampling (bWs),
 #'  Bandwidth of Altman and Leger (AL), One-sided Cross Validation (OCV), Akaike information criterion (AIC), Indirect Cross Validation (ICV), Mallow’ Cp (MallowCp), Generalized Cross Validation (GCV), Polansky and Baker Plug-In (PB),
 #'  and Gasser, Kniep, and Köhler Cross Validation (GKK). For RIG kernel see \code{\link{righcomp}}
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
 #' @author Javaria Ahmad Khan, Atif Akbar.
 #' @references Khan, J. A.; Akbar, A. Density Estimation by Laplace Kernel. \emph{Working paper,  Department of Statistics, Bahauddin Zakariya University, Multan, Pakistan.}
 #' @export
 laphcomp<-function(y,k,type){
   n <- length(y)
   sd<-sd(y)
   x <- seq(min(y) + 0.05, max(y), length = k)

   PI <- abs(dpik(y) )
   AL<-abs(ALbw(vec_data=y))
   bw2<-h.bcv(y)
   BCV<-abs(bw2$h)
   bw3<-h.ccv(y)
   CCV<-abs(bw3$h)
   RoT <-abs(1.06*sd(y)*(n^(-1/5)))
   bw8<-h.ucv(y)
   UCV<-abs(bw8$h)
   bw10<-kdeb(y, h0 = 0.01 * sd, h1 = sd, meth = c("AIC", "LCV", "LSCV", "BCV","SJPI", "GKK"), kern = "gauss", gf = 2.5)
   AIC<-abs(bw10[1])
   GKK<-abs(bw10[6])
   ICV<-abs(h_ICV(y))
   bw8<-cp(y,  sig2=1)
   mallow<-abs(bw8[5])
   bWr<-bw.dboot2(y,sd,h0='dboot1',error='normal',B=1000,grid=100,ub=2)
   bw=gcv(y)
   gcv=abs(bw[4] )
   bw5<-h.mcv(y)
   MCV<-abs(bw5$h)
   bw6<-h.mlcv(y)
   MLCV<-abs(bw6$h)
   bw7<-h.tcv(y)
   TCV<-abs(bw7$h)
   SCV<-abs(hscv(y))
   bWOr<-bw.dboot1(y,sd, h0="dnrd", error="normal", grid=100, ub=2)
   pb<-PBbw( vec_data=y)
   oscv<-h_OSCV_dens(y,0)

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
     }

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

   }
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
 #' @details It plot the densities by Laplace, RIG kernel and with real densities at the same time.
 #' @param y a numeric vector of positive values.
 #' @param k gird points.
 #' @param h the bandwidth
 #' @import graphics
 #' @import stats
 #' @examples
 #' y<-rexp(100,1)
 #' h<-0.79 * IQR(y) * length(y) ^ (-1/5)
 #' compLR(y,80,h)
 #' @author Javaria Ahmad Khan, Atif Akbar.
 #' @references Khan, J. A.; Akbar, A. Density Estimation by Laplace Kernel. \emph{Working paper,  Department of Statistics, Bahauddin Zakariya University, Multan, Pakistan.}
 #' Scaillet, O. 2004. Density estimation using inverse and reciprocal inverse Gaussian kernels. 	 \emph{Nonparametric Statistics}, \strong{16}, 217-226.
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
