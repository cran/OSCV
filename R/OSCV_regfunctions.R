#' Computing the local linear estimate (LLE).
#'
#' Computing the LLE based on data \eqn{(desx,y)} over the given vector of the argument values \eqn{u}. The Gausssian kernel is used. See expression (3) in Savchuk and Hart (2017).
#'
#' Computing the LLE based on the Gaussian kernel for the specified vector of the argument values \eqn{u} and given vectors of design points \eqn{desx} and the corresponding data values \eqn{y}.
#' @param u numerical vector of argument values,
#' @param desx numerical vecror of design points,
#' @param y numerical vecror of data values (corresponding to the specified design points \eqn{desx}),
#' @param h numerical bandwidth value (scalar).
#' @return Numerical vector of the LLE values computed over the specified vector of \eqn{u} points.
#' @references
#' \itemize{
#'   \item Clevelend, W.S. (1979). Robust locally weighted regression and smoothing scatterplots. \emph{Journal of the American Statistical Association}, 74(368), 829-836.
#'   \item  Savchuk, O.Y., Hart, J.D. (2017). Fully robust one-sided cross-validation for regression functions. \emph{Computational Statistics}, doi:10.1007/s00180-017-0713-7.
#' }
#' @seealso \code{\link{OSCV_reg}}, \code{\link{h_OSCV_reg}}, \code{\link{ASE_reg}}, \code{\link{h_ASE_reg}}, \code{\link{CV_reg}}.
#' @examples
#' \dontrun{
#' # Example (simulated data).
#' n=200
#' dx=(1:n-0.5)/n
#' regf=2*dx^10*(1-dx)^2+dx^2*(1-dx)^10
#' u=seq(0,1,len=1000)
#' ydat=regf+rnorm(n,sd=0.002)
#' dev.new()
#' plot(dx,regf,'l',lty="dashed",lwd=3,xlim=c(0,1),ylim=c(1.1*min(ydat),1.1*max(ydat)),
#' cex.axis=1.7,cex.lab=1.7)
#' title(main="Function, generated data, and LLE",cex.main=1.5)
#' points(dx,ydat,pch=20,cex=1.5)
#' lines(u,loclin(u,dx,ydat,0.05),lwd=3,col="blue")
#' legend(0,1.1*max(ydat),legend=c("LLE based on h=0.05","true regression function"),
#' lwd=c(2,3),lty=c("solid","dashed"),col=c("blue","black"),cex=1.5,bty="n")
#' legend(0.7,0.5*min(ydat),legend="n=200",cex=1.7,bty="n")
#' }
#' @export
#' @importFrom "stats" "dnorm"
loclin=function(u,desx,y,h){
  n=length(desx)
  m=length(u)
  arg=matrix(u,m,1)%*%matrix(1,1,n)-matrix(1,m,1)%*%matrix(desx,1,n)
  tn0=dnorm(arg/h)%*%matrix(1,n,1)
  tn1=(dnorm(arg/h)*arg)%*%matrix(1,n,1)
  tn2=(dnorm(arg/h)*arg^2)%*%matrix(1,n,1)
  denom=tn2*tn0-tn1^2
  S1=dnorm(arg/h)%*%matrix(y,n,1)
  S2=(dnorm(arg/h)*arg)%*%matrix(y,n,1)
  numer=tn2*S1-tn1*S2
  as.numeric(numer/denom)
}


#' The family of two-sided cross-validation kernels \eqn{H_I}.
#'
#' The family of two-sided cross-validation kernels \eqn{H_I} defined by equation (15) of Savchuk and Hart (2017).
#'
#' The family of the two-sided cross-validation kernels \eqn{H_I(u;\alpha,\sigma)=(1+\alpha)\phi(u)-\alpha\phi(u/\sigma)/\sigma}, where \eqn{\phi} denotes the Gaussian kernel, \eqn{-\infty<\alpha<\infty} and \eqn{\sigma>0} are the parameters of the kernel. See expression (15) of Savchuk and Hart (2017). The robust kernel plotted in Figure 1 of Savchuk and Hart (2017) is obtained by setting \eqn{\alpha=16.8954588} and \eqn{\sigma=1.01}. Note that the kernels \eqn{H_I} are also used for the bandwidth selection purposes in the indirect cross-validation (ICV) method (see expression (4) of Savchuk, Hart, and Sheather (2010)). The kernel \eqn{H_I} is a two-sided analog of the one-sided kernel \code{\link{L_I}}. The Gaussian kernel \eqn{\phi} is the special case of \eqn{H_I} obtained by either setting \eqn{\alpha=0} or \eqn{\sigma=1}.
#' @param u numerical vector of argument values,
#' @param alpha first parameter of the cross-validation kernel \eqn{H_I},
#' @param sigma second parameter of the cross-validation kernel \eqn{H_I}.
#' @return The value of \eqn{H_I(u;\alpha,\sigma)}.
#' @references
#' \itemize{
#'   \item Savchuk, O.Y., Hart, J.D. (2017). Fully robust one-sided cross-validation for regression functions. \emph{Computational Statistics}, doi:10.1007/s00180-017-0713-7.
#'   \item Savchuk, O.Y., Hart, J.D., Sheather, S.J. (2010). Indirect cross-validation for density estimation. \emph{Journal of the American Statistical Association}, 105(489), 415-423.
#' }
#' @seealso \code{\link{L_I}}, \code{\link{C_smooth}}, \code{\link{OSCV_reg}}, \code{\link{loclin}}.
#' @examples
#' \dontrun{
#' # Plotting the robust kernel from Savchuk and Hart (2017) with alpha=16.8954588 and sigma=1.01.
#' u=seq(-5,5,len=1000)
#' ker=H_I(u,16.8954588,1.01)
#' dev.new()
#' plot(u,ker,'l',lwd=3,cex.axis=1.7, cex.lab=1.7)
#' title(main="Robust kernel H_I along with the Gaussian kernel (phi)",cex=1.7)
#' lines(u,dnorm(u),lty="dashed",lwd=3)
#' legend(-4.85,0.3,lty=c("solid","dashed"),lwd=c(3,3),legend=c("H_I","phi"),cex=1.5)
#' legend(1,0.4,legend=c("alpha=16.8955","sigma=1.01"),cex=1.5,bty="n")
#' }
#' @export
#' @importFrom "stats" "dnorm"
H_I=function(u,alpha,sigma)
  (1+alpha)*dnorm(u)-alpha/sigma*dnorm(u/sigma)


#' The OSCV function in the regression context.
#'
#' Computing \eqn{OSCV(b)}, the value of the OSCV function in the regression context, defined by expression (9) of Savchuk and Hart (2017).
#'
#' Computation of  \eqn{OSCV(b)} for given \eqn{b} (bandwidth vector) and the data values \eqn{y} corresponding to the design points \eqn{desx}. No preliminary sorting of the data (according to the \eqn{desx} variable) is needed. The value of \eqn{m=4} is used. Two choices of the two-sided cross-validation kernel are available: \itemize{\item (\eqn{ktype=0}) Gaussian kernel; \item (\eqn{ktype=1}) robust kernel \code{\link{H_I}} defined by expression (15) of Savchuk and Hart (2017) with \eqn{(\alpha,\sigma)=(16.8954588,1.01)}.}
#' @param b numerical vector of bandwidth values,
#' @param desx numerical vecror of design points,
#' @param y numerical vecror of data points corresponding to the design points \eqn{desx},
#' @param ktype making choice between two cross-validation kernels:  (\eqn{ktype=0}) corresponds to the Gaussian kernel; (\eqn{ktype=1}) corresponds to the robust kernel \code{\link{H_I}} with \eqn{(\alpha,\sigma)=(16.8954588,1.01)}.
#' @return The vector of values of \eqn{OSCV(b)} for the correponsing vector of \eqn{b} values.
#' @references 
#' \itemize{
#' \item Savchuk, O.Y., Hart, J.D. (2017). Fully robust one-sided cross-validation for regression functions. \emph{Computational Statistics}, doi:10.1007/s00180-017-0713-7.
#' \item Hart, J.D. and Yi, S. (1998) One-sided cross-validation. \emph{Journal of the American Statistical Association}, 93(442), 620-631.
#' }
#' @seealso \code{\link{h_OSCV_reg}}, \code{\link{H_I}}, \code{\link{loclin}}, \code{\link{C_smooth}}.
#' @examples
#' \dontrun{
#' # The Old Faithful geyser data set "faithful" is used. The sample size n=272.
#' # The OSCV curves based on the Gaussian kernel and the robust kernel H_I (with 
#' # alpha=16.8954588 and sigma=1.01) are plotted. The horizontal scales of the curves
#' # are changed such that their global minimizers are to be used in computing the
#' # Gaussian local linear estimates of the regression function.
#' xdat=faithful[[2]] #waiting time
#' ydat=faithful[[1]] #eruption duration
#' barray=seq(0.5,10,len=250)
#' C_gauss=C_smooth(1,1)
#' OSCV_gauss=OSCV_reg(barray/C_gauss,xdat,ydat,0)
#' h_gauss=round(h_OSCV_reg(xdat,ydat,0),digits=4)
#' dev.new()
#' plot(barray,OSCV_gauss,'l',lwd=3,cex.lab=1.7,cex.axis=1.7,xlab="h",ylab="OSCV criterion")
#' title(main="OSCV based on the Gaussian kernel",cex.main=1.7)
#' legend(2.5,0.25,legend=paste("h_min=",h_gauss),cex=2,bty="n")
#' C_H_I=C_smooth(16.8954588,1.01)
#' OSCV_H_I=OSCV_reg(barray/C_H_I,xdat,ydat,1)
#' h_H_I=round(barray[which.min(OSCV_H_I)],digits=4)
#' dev.new()
#' plot(barray,OSCV_H_I,'l',lwd=3,cex.lab=1.7,cex.axis=1.7,xlab="h",ylab="OSCV criterion",
#' ylim=c(0.15,0.5))
#' title(main="OSCV based on the robust kernel H_I",cex.main=1.7)
#' legend(2.5,0.4,legend=paste("h_min=",h_H_I),cex=2,bty="n")
#' }
#' @export
#' @importFrom "stats" "dnorm"
OSCV_reg=function(b,desx,y,ktype){
  ind=order(desx)            #Ordering the data
  desx=desx[ind]
  y=y[ind]
  H=function(u,ktype){                #choosing two-sided cross-validation kernel
    if(ktype==0)        dnorm(u)                                                    #Gaussian kernel
    else if(ktype==1)   17.8954588*dnorm(u)-16.8954588/1.01*dnorm(u/1.01)           #H_I (robust)
  }
  n=length(y)
  m=4                                     #the first four points are skipped
  ngrid=length(b)
  oscv=numeric(ngrid)
  for(j in 1:ngrid){
    for(i in (m+1):n){
      arg=(desx[i]-desx)[1:(i-1)]
      kern=H(arg/b[j],ktype)
      tn1=sum(kern*arg)
      tn2=sum(kern*(arg^2))
      weight=kern*(tn2-arg*tn1)
      weight=weight/sum(weight)
      pred=sum(weight*y[1:(i-1)])           #value of one-sided estimator at the i-th point
      oscv[j]=oscv[j]+(pred-y[i])^2
    }
  }
  oscv/(n-m)
}


#' The OSCV bandwidth in the regression context.
#'
#' Computing the OSCV bandwidth for the Gaussian local linear regression estimator. The Gaussian kernel is used in the bandwidth selection stage. The smoothness of the regression function is to be specified by the user.
#'
#' Computing the OSCV bandwidth for the data vector \eqn{(desx,y)}. The Gaussian kernel is used for the cross-validation purposes and in the stage of computing the resulting local linear regression estimate. No additional rescaling of the computed bandwidth is needed. The smoothness of the regression function \eqn{stype}, essentially, determines the value of the bandwidth rescaling constant that is chosen in the body of the function. Thus, the constant is equal to 0.6168471 in the smooth case, and 0.5730 in the nonsmooth case. See Savchuk, Hart and Sheather (2016). The OSCV bandwidth is the minimizer of the OSCV function \code{\link{OSCV_reg}}.
#' @param desx numerical vecror of design points,
#' @param y numerical vecror of data points corresponding to the design points \eqn{desx},
#' @param stype smoothness of the regression function: (\eqn{stype=0}) smooth function; (\eqn{stype=1}) nonsmooth function.
#' @return The OSCV bandwidth (scalar).
#' @references
#' \itemize{
#'   \item Hart, J.D. and Yi, S. (1998). One-sided cross-validation. \emph{Journal of the American Statistical Association}, 93(442), 620-631.
#'   \item Savchuk, O.Y., Hart, J.D., Sheather, S.J. (2013). One-sided cross-validation for nonsmooth regression functions. \emph{Journal of Nonparametric Statistics}, 25(4), 889-904.
#'   \item Savchuk, O.Y., Hart, J.D., Sheather, S.J. (2016). Corrigendum to "One-sided cross-validation for nonsmooth regression functions". \emph{Journal of Nonparametric Statistics}, 28(4), 875-877.
#'   \item Savchuk, O.Y., Hart, J.D. (2017). Fully robust one-sided cross-validation for regression functions. \emph{Computational Statistics}, doi:10.1007/s00180-017-0713-7.
#' }
#' @seealso \code{\link{OSCV_reg}}, \code{\link{loclin}}, \code{\link{C_smooth}}, \code{\link{h_OSCV_dens}}, \code{\link{h_ASE_reg}}.
#' @examples
#' \dontrun{
#' # Example (Old Faithful geyser)
#' xdat=faithful[[2]]     # waiting time
#' ydat=faithful[[1]]     # eruption duration
#' u=seq(40,100,len=1000)
#' h_oscv=round(h_OSCV_reg(xdat,ydat,0),digits=4)
#' l=loclin(u,xdat,ydat,h_oscv)
#' dev.new()
#' plot(xdat,ydat,pch=20,cex=1.5,cex.axis=1.7,cex.lab=1.7,xlab="waiting time",
#' ylab="eruption duration")
#' lines(u,l,'l',lwd=3)
#' title(main="Data and LLE",cex.main=1.7)
#' legend(35,5,legend=paste("h_OSCV=",h_oscv),cex=2,bty="n")
#' legend(80,3,legend="n=272",cex=2,bty="n")
#'}
#' @export
#' @importFrom "stats" "optimize"
h_OSCV_reg=function(desx,y,stype){
  if(stype==0) Const=0.6168471
  else if(stype==1) Const=0.5730
  ind=order(desx)            #Sorting the data
  desx=desx[ind]
  y=y[ind]
  n=length(desx)
  spacemax=max(desx[2:n]-desx[1:(n-1)])    #maximum spacing
  R=range(desx)[2]-range(desx)[1]
  h_G=Const*optimize(OSCV_reg,c(spacemax/(Const*100),R/(Const*4)),tol=0.0001,ktype=0,desx=desx,y=y)$minimum
  h_G
}


#' Nonsmooth regression function with six cusps.
#'
#' Nonsmooth regression function \eqn{r_3} with six cusps used in the simulation studies in Savchuk et al. (2013) and Savchuk et al. (2017).
#'
#' The nonsmooth function \eqn{r_3} can be used in simulation studies.
#' @param u numerical vecror of argument values in the range [0,1].
#' @return The vector of values of \eqn{r_3} corresponding to the values of the vector \eqn{u}.
#' @references
#' \itemize{
#'   \item Savchuk, O.Y., Hart, J.D., Sheather, S.J. (2013). One-sided cross-validation for nonsmooth regression functions. \emph{Journal of Nonparametric Statistics}, 25(4), 889-904.
#'   \item Savchuk, O.Y., Hart, J.D. (2017). Fully robust one-sided cross-validation for regression functions. \emph{Computational Statistics}, doi:10.1007/s00180-017-0713-7.
#' }
#' @examples
#' \dontrun{
#' # n=250 data points are generated from r3 by adding the Gaussian noise with sigma=1/500.
#' # The fixed evenly spaced design is used.
#' u=seq(0,1,len=1000)
#' n=250
#' xdat=(1:n-0.5)/n
#' ydat=reg3(xdat)+rnorm(n,sd=1/500)
#' h_oscv=round(h_OSCV_reg(xdat,ydat,1),digits=4) # L_G-based OSCV based on nonsmooth constant
#' l=loclin(u,xdat,ydat,h_oscv)
#' dev.new()
#' plot(xdat,ydat,pch=20,cex=1.5,cex.axis=1.5,cex.lab=1.5,xlab="x",ylab="y",
#' ylim=c(min(ydat),1.2*max(ydat)))
#' lines(u,l,'l',lwd=3,col="blue")
#' lines(u,reg3(u),lwd=3,lty="dashed")
#' title(main="Data, true regression function and LLE",cex.main=1.7)
#' legend(-0.05,0.003,legend=paste("h_OSCV=",h_oscv),cex=2,bty="n")
#' legend(0.65,0.025, legend="n=250",cex=2,bty="n")
#' legend(0,1.28*max(ydat),legend=c("LLE based on h_OSCV","true regression function"),lwd=c(3,3),
#' lty=c("solid","dashed"),col=c("blue","black"),bty="n",cex=1.5)
#'}
#' @export
reg3=function(u){
  rfun=array(0,length(u))
  for(i in 1:length(u)){
    if(u[i]<0.1) rfun[i]=0.4761904762*sqrt(u[i])
    else if(u[i]<0.3)  rfun[i]=1/10*0.4761904762*exp(-20*u[i]+2)+0.1029656029
    else if(u[i]<0.35) rfun[i]=1.428571429*u[i]-0.3247336524
    else if(u[i]<0.6)  rfun[i]=1.428571425*(u[i]-0.35)*(u[i]-0.45)+0.1752663476
    else if(u[i]<0.7)  rfun[i]=-2.142857143*u[i]+1.514552062
    else if(u[i]<0.8)  rfun[i]=-2.142857143*(u[i]-0.7)^3*(u[i]-0.4)+0.01455206186
    else if(u[i]<=1)   rfun[i]=0.04761904762*log(10*u[i]-7.9)+0.1233418282
  }
  rfun/10
}

#' The ASE function for the local linear estimator (LLE) in the regression context.
#'
#' Computing \eqn{ASE(h)}, the value of the ASE function for the local linear estimator in the regression context, for the given vector of \eqn{h} values.
#'
#' The average squared error (ASE) is used as a measure of performace of the local linear estimator based on the Gaussian kernel.
#' @param h numerical vector of bandwidth values,
#' @param desx numerical vecror of design points,
#' @param y numerical vecror of data points corresponding to the design points \eqn{desx},
#' @param rx numerical vecror of values of the regression function at \eqn{desx}.
#' @return The vector of values of \eqn{ASE(h)} for the correponsing vector of \eqn{h} values.
#' @references 
#' Hart, J.D. and Yi, S. (1998) One-sided cross-validation. \emph{Journal of the American Statistical Association}, 93(442), 620-631.
#' @seealso \code{\link{loclin}}, \code{\link{h_ASE_reg}}, \code{\link{CV_reg}}, \code{\link{OSCV_reg}}.
#' @examples
#' \dontrun{
#' # Example (ASE function for a random sample of size n=100 generated from the function reg3 that
#' # has six cusps. The function originates from the article of Savchuk et al. (2013).
#' # The level of the added Gaussian noise is sigma=1/1000).
#' n=100
#' dx=(1:n-0.5)/n
#' regx=reg3(dx)
#' ydat=regx+rnorm(n,sd=1/1000)
#' harray=seq(0.003,0.05,len=300)
#' ASEarray=ASE_reg(harray,dx,ydat,regx)
#' hmin=round(h_ASE_reg(dx,ydat,regx),digits=4)
#' dev.new()
#' plot(harray,ASEarray,'l',lwd=3,xlab="h",ylab="ASE",main="ASE function for a random sample
#' from r3",cex.lab=1.7,cex.axis=1.7,cex.main=1.5)
#' legend(0.029,0.0000008,legend=c("n=100","sigma=1/1000"),cex=1.7,bty="n")
#' legend(0.005,0.000002,legend=paste("h_ASE=",hmin),cex=2,bty="n")
#' }
#' @export
ASE_reg=function(h,desx,y,rx){
  ngrid=length(h)
  ASE=numeric(ngrid)
  for(i in 1:ngrid)
    ASE[i]=mean((loclin(desx,desx,y,h[i])-rx)^2)
  ASE
}

#' The ASE-optimal bandwidth in the regression context.
#'
#' Computing the ASE-optimal bandwidth for the Gaussian local linear regression estimator.
#'
#' Computing the ASE-optimal bandwidth for the local linear estimator in the regression context. The ASE-optimal bandwidth is the global minimizer of the ASE function \code{\link{ASE_reg}}. This bandwidth is optimal for the data set at hand.
#' @param desx numerical vecror of design points,
#' @param y numerical vecror of data points corresponding to the design points \eqn{desx},
#' @param rx numerical vecror of the regression function values at \eqn{desx}.
#' @return The ASE-optimal bandwidth (scalar).
#' @seealso \code{\link{ASE_reg}}, \code{\link{loclin}}.
#' @examples
#' \dontrun{
#' # Simulated example.
#' n=300
#' dx=runif(n)            #uniform design
#' regx=5*dx^10*(1-dx)^2+2.5*dx^2*(1-dx)^10
#' ydat=regx+rnorm(n,sd=1/250)
#' hase=round(h_ASE_reg(dx,ydat,regx),digits=4)
#' u=seq(0,1,len=1000)
#' fun=5*u^10*(1-u)^2+2.5*u^2*(1-u)^10
#' dev.new()
#' plot(dx,ydat,pch=20,cex=1.5,xlab="argument",ylab="function",cex.lab=1.7,cex.axis=1.7,
#' main="Function, data, and the ASE-optimal bandwidth",cex.main=1.5)
#' lines(u,fun,'l',lwd=3,col="blue")
#' legend(0,0.03,legend=paste("h_ASE=",hase),cex=1.8,bty="n")
#' legend(0.6,-0.002,legend=paste("n=",n),cex=2,bty="n")
#'}
#' @export
#' @importFrom "stats" "optimize"
h_ASE_reg=function(desx,y,rx){
  ind=order(desx)            #Sorting the data
  desx=desx[ind]
  y=y[ind]
  rx=rx[ind]
  n=length(desx)
  spacemax=max(desx[2:n]-desx[1:(n-1)])    #maximum spacing
  R=range(desx)[2]-range(desx)[1]
  optimize(ASE_reg,c(spacemax/100,R/4),tol=0.0001,desx=desx,y=y,rx=rx)$minimum
}


#' The cross-validation (CV) function in the regression context.
#'
#' Computing \eqn{CV(h)}, the value of the CV function in the regression context.
#'
#' The CV function is a measure of fit of the regression estimate to the data. The local linear estimator based on the Gaussian kernel is used. The cross-validation bandwidth is the minimizer of the CV function.
#' @param h numerical vector of bandwidth values,
#' @param desx numerical vecror of design points,
#' @param y numerical vecror of data values corresponding to the design points \eqn{desx}.
#' @return The vector of values of \eqn{CV(h)} for the correponsing vector of \eqn{h} values.
#' @references 
#' Stone, C.J. (1977) Consistent nonparametric regression. \emph{Annals of Statistics}, 5(4), 595-645.
#' @seealso \code{\link{loclin}}, \code{\link{h_ASE_reg}}, \code{\link{ASE_reg}}, \code{\link{OSCV_reg}}.
#' @examples
#' \dontrun{
#' # Example (Old Faithful geyser). Take x=waiting time; y=eruption duration. The sample size n=272.
#' xdat=faithful[[2]]
#' ydat=faithful[[1]]
#' harray=seq(0.5,10,len=100)
#' cv=CV_reg(harray,xdat,ydat)
#' R=range(xdat)
#' h_cv=round(optimize(CV_reg,c(0.01,(R[2]-R[1]/4)),desx=xdat,y=ydat)$minimum,digits=4)
#' dev.new()
#' plot(harray,cv,'l',lwd=3,xlab="h",ylab="CV(h)",main="CV function for the Old Faithful 
#' geyser data", cex.lab=1.7,cex.axis=1.7,cex.main=1.5)
#' legend(6,0.155,legend="n=272",cex=1.8,bty="n")
#' legend(1,0.18,legend=paste("h_CV=",h_cv),cex=2,bty="n")
#' }
#' @export
CV_reg=function(h,desx,y){
  n=length(desx)
  ngrid=length(h)
  CVarray=numeric(ngrid)
  for(j in 1:ngrid){
    cv=0
    for(i in 1:n)
      cv=cv+(loclin(desx[i],desx[-i],y[-i],h[j])-y[i])^2
    CVarray[j]=cv
  }
  CVarray/n
}




  


