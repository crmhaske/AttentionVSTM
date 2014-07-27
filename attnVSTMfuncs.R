## =============================================================================
## Functions for attnVSTM MS analysis v.1
## By Christie Haskell, 27-07-2014 (v.1)
## Source:
## =============================================================================

#Cohen's D
cohens_d <- function(v1, v2) {
  l1 <- length(v1)- 1
  l2 <- length(v2)- 1
  m  <- abs(mean(v1) - mean(v2))
  sd <- l1 * var(v1) + l2 * var(v2)
  d <- m/sqrt(sd/(l1 + l2))
  return(d)
}

#Selection criteria
#Proportion of trials on which the subect rotated the grating in the correct direction
prop.cor<-function(data) {
   userid<-unique(data$userid)
   sel<-matrix(NA,length(userid),2)
   colnames(sel)<-c("userid","prop.cor")
   sel<-as.data.frame(sel)
   n=1

   for (i in 1:length(userid)) {
       sel[n,1]<-userid[n]

       #userid
       duser<-subset(data,data$userid==sel[n,1])

       #orient < 0
       duser1<-subset(duser,((duser$orient+ 90) %% 180 - 90)<0)
       n.ori.bel<-length(duser1$orient)

       #judgori < 0
       duser2<-subset(duser1,((duser1$judgori+90) %% 180 - 90)<0)
       n.jori.bel<-length(duser2$judgori)

       #orient > 0
       duser3<-subset(duser,((duser$orient+90) %% 180 - 90)>0)
       n.ori.abo<-length(duser3$orient)

       #judgori > 0
       duser4<-subset(duser3,((duser3$judgori+90) %% 180 - 90)>0)
       n.jori.abo<-length(duser4$judgori)

       #prop.cor
       sel[n,2]<-(n.jori.bel+n.jori.abo)/(n.ori.bel+n.ori.abo)
       n=n+1
   }
   return(sel)
}

#Estimated starting parameters
phi.hat <- function(data) {
    n<-length(data)

    LL <- function(par, x) { 
             -sum(log(exp(par[2]*cos(2*x-par[1]))/(pi*besselI(par[2],0)))) 
          }

    vmmle<-optim(c(0,3), LL, x = data)

    muo_hat<-vmmle$par[1]

    nq_max<-vmmle$par[1]+3*pi/8
    nq_min<-vmmle$par[1]-3*pi/8
    NQ1<-data[data>nq_max]
    NQ2<-data[data<nq_min]
    nq<-length(NQ1)+length(NQ2)
    po_hat<-(n-(4*nq))/n
    if(po_hat<0){
         po_hat<-0.5
    }
    #Initial estimate for kappa
    kappao_hat<-vmmle$par[2]

    phi.hat<-c(po_hat, muo_hat, kappao_hat)
    return(phi.hat)
}

#MLE of pfvm+(1-p)fu
mixedfit<- function(data) {

   e<-data*pi/180

   p0<-phi.hat(e)

   LL <- function(par, x) { 
           -sum(log(par[1]*exp(par[3]*cos(2*x-par[2]))/(pi*besselI(par[3],0))+(1-par[1])*1/pi)) 
      }

   #Minimizing the negative sum (i.e. maximizing the log-likelihood function, MLE)
   result<-optim(p0, LL, x = e, method = "L-BFGS-B", lower=c(0,-pi/2,0.1), upper=c(1,pi/2,Inf))

   mu<-round(result$par[2],2)
   kappa<-round(result$par[3],2)
   sigma<-round(sqrt(-log(besselI(result$par[3],1)/besselI(result$par[3],0))),2)
   g<-round(1-result$par[1],2)

   fn<-function(x) {
          result$par[1]*exp(result$par[3]*cos(2*x-result$par[2]))/(pi*besselI(result$par[3],0))+(1-result$par[1])*rep(1/pi,length(x))
       }
   dorig<-hist(e,plot="FALSE",20)
   kstest2<-ks.test(dorig$density, fn(dorig$mids))
   kstestl<-ks.test(dorig$density, fn(dorig$mids),alternative="l")
   kstestg<-ks.test(dorig$density, fn(dorig$mids),alternative="g")
   kstest<-c(kstest2$p.value,kstestl$p.value,kstestg$p.value)
   cor<-max(kstest)

   params<-list("mu" = mu, "kappa" = kappa, "sigma" = sigma,"g" = g, "cor" = cor)

   return(params)

}

#Calculate confidence interval
ci<-function(x,opt_arg) {
   
   if(missing(opt_arg)) { 
      lvl <- .05 
   } else {
      lvl <- opt_arg 
   }   

   if (lvl==.2) {
      z<-1.28
   } else if (lvl ==.1) {
      z<-1.645
   } else if (lvl == .05) {
      z<-1.96
   } else if (lvl ==.02) {
      z<-2.33
   } else if (lvl == .01) {
      z<-2.58
   }

   sd<-sd(x)
   m<-mean(x)
   n<-length(x)

   lower<-m-z*sqrt((sd^2/n))
   upper<-m+z*sqrt((sd^2/n))

   ci<-list("mean" = m, "lower" = lower, "upper" = upper)
   return(ci)
}