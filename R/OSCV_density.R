#' Nonsmooth density function with seven cusps.
#'
#' Nonsmooth density \eqn{f^*} with seven cusps introduced in the article of Savchuk (2017).
#'
#' The function \eqn{f^*} consists of straight lines with different slopes connected together. The support of the density is [-3,3].
#' @param u numerical vecror of argument values in the range [-3,3].
#' @return The vector of values of \eqn{f^*} corresponding to the values of the vector \eqn{u}.
#' @references
#' Savchuk, O.Y. (2017). One-sided cross-validation for nonsmooth density functions, arXiv:1703.05157.
#' @seealso \code{\link{sample_fstar}}, \code{\link{ISE_fstar}}.
#' @examples
#' \dontrun{
#' dev.new()
#' plot(seq(-3.5,3.5,len=1000),fstar(seq(-3.5,3.5,len=1000)),'l',lwd=3,
#' main="Nonsmooth density fstar with seven cusps", xlab="argument", ylab="density",cex.main=1.5,
#' cex.axis=1.7,cex.lab=1.7)
#'}
#' @export
fstar=function(u){
  (1/8*u+3/8)*sapply(u,function(t){if (t>-3 && t<(-3/2)) 1 else 0})+
  (-1/4*u-3/16)*sapply(u,function(t){if (t>-3/2 && t<(-5/4)) 1 else 0})+
  (4/15*u+11/24)*sapply(u,function(t){if (t>-5/4 && t<(-1/2)) 1 else 0})+
  (-2/5*u+1/8)*sapply(u,function(t){if (t>-1/2 && t<0) 1 else 0})+
  (17/48*u+1/8)*sapply(u,function(t){if (t>0 && t<1/2) 1 else 0})+
  (-17/96*u+25/64)*sapply(u,function(t){if (t>1/2 && t<3/2) 1 else 0})+
  (1/4*u-1/4)*sapply(u,function(t){if (t>3/2 && t<2) 1 else 0})+
  (-1/4*u+3/4)*sapply(u,function(t){if (t>2 && t<3) 1 else 0})
}

#' Taking a random sample from \code{\link{fstar}}.
#'
#' Taking a random sample of size \eqn{n} from the density \eqn{f^*} with seven cusps introduced in the article of Savchuk (2017).
#'
#' The density \eqn{f^*} can be used in simulation studies.
#' @param n sample size.
#' @return The numerical vector of size \eqn{n} of the data values.
#' @references
#' Savchuk, O.Y. (2017). One-sided cross-validation for nonsmooth density functions, arXiv:1703.05157.
#' @seealso \code{\link{fstar}}, \code{\link{ISE_fstar}}.
#' @examples
#' \dontrun{
#' dev.new()
#' plot(density(sample_fstar(5000),bw=0.1),lwd=2,ylim=c(0,0.32),xlab="argument",ylab="density",
#' main="KDE and the true density fstar",cex.lab=1.7, cex.axis=1.7,cex.main=1.7)
#' lines(seq(-3.5,3.5,len=1000),fstar(seq(-3.5,3.5,len=1000)),lwd=3,lty="dashed")
#' legend(-3,0.3,legend=c("KDE","True density","h=0.1","n=5000"),lwd=c(2,3),
#' lty=c("solid","dashed"),col=c("black","black","white","white"))
#'}
#' @export
#' @importFrom "mc2d" "rtriang"
#' @importFrom "stats" "runif"
sample_fstar=function(n){
  u=runif(n)
  dat=array(0,n)
  f1=(1:n)[u<=9/16]
  dat[f1]=runif(length(f1),-2,5/2)
  f2=(1:n)[(u>9/16)&(u<=5/8)]
  dat[f2]=rtriang(length(f2),min=-3,mode=-2,max=-2)
  f3=(1:n)[(u>5/8)&(u<=83/128)]
  dat[f3]=rtriang(length(f3),min=-2,mode=-3/2,max=-5/4)
  f4=(1:n)[(u>83/128)&(u<=99/128)]
  dat[f4]=rtriang(length(f4),min=-5/4,mode=-1/2,max=0)
  f5=(1:n)[(u>99/128)&(u<=29/32)]
  dat[f5]=rtriang(length(f5),min=0,mode=1/2,max=3/2)
  f6=(1:n)[(u>29/32)&(u<=31/32)]
  dat[f6]=rtriang(length(f6),min=3/2,mode=2,max=5/2)
  f7=(1:n)[u>31/32]
  dat[f7]=rtriang(length(f7),min=5/2,mode=5/2,max=3)
  dat
}


#' The ISE function in the kernel density estimation (KDE) context in the case when the underlying density is \code{\link{fstar}}.
#'
#' Computing the ISE function for the Gaussian density estimator obtained from a random sample of size \eqn{n} generated from \code{\link{fstar}}.
#'
#' The integrated squared error (ISE) is a measure of closeness of the Gaussian density estimate computed from a data set generated from \code{\link{fstar}} to the true density.
#' @param h numerical vector of bandwidth values,
#' @param n sample size (number of data points generated from \code{\link{fstar}}).
#' @return The vector of values of the ISE function for the correponsing vector of \eqn{h} values.
#' @references 
#' Savchuk, O.Y. (2017). One-sided cross-validation for nonsmooth density functions, arXiv:1703.05157.
#' @seealso \code{\link{fstar}}, \code{\link{sample_fstar}}.
#' @examples
#' \dontrun{
#' dev.new()
#' harray=seq(0.05,1.5,len=1000)
#' ISEarray=ISE_fstar(harray,100)
#' h_ISE=round(harray[which.min(ISEarray)],digits=4)
#' dev.new()
#' plot(harray,ISEarray,lwd=3,'l',xlab="h",ylab="ISE",main="ISE(h)",cex.main=2,cex.lab=1.7,
#' cex.axis=1.7)
#' legend(0.35,ISEarray[5],legend=c("n=100",paste("h_ISE=",h_ISE)),cex=1.8,bty="n")
#' }
#' @export
#' @importFrom "stats" "dnorm" "pnorm"

ISE_fstar=function(h,n){
  dat=sample_fstar(n)
  ngrid=length(h)
  R.f=18671/92160
  arg=matrix(dat,n,1)%*%matrix(1,1,n)-matrix(1,n,1)%*%matrix(dat,1,n)
  a=c(-3,-3/2,-5/4,-1/2,0,1/2,3/2,2) #left endpoint of intervals
  b=c(-3/2,-5/4,-1/2,0,1/2,3/2,2,3) #right endpoint of intervals
  k=c(1/8,-1/4,4/15,-2/5,17/48,-17/96,1/4,-1/4)    #slope coefficients
  C=c(3/8,-3/16,11/24,1/8,1/8,25/64,-1/4,3/4)      #intercepts of lines
  Amatrix=matrix(a,8,1)%*%matrix(1,1,n)
  Bmatrix=matrix(b,8,1)%*%matrix(1,1,n)
  Kmatrix=matrix(k,8,1)%*%matrix(1,1,n)
  Cmatrix=matrix(C,8,1)%*%matrix(1,1,n)
  DATmatrix=matrix(1,8,1)%*%matrix(dat,1,n)
  ISEfstar_vec=1:ngrid
  for(i in 1:ngrid){
    R.fhat=1/n^2*sum(dnorm(arg,sd=h[i]*sqrt(2)))
    M1=(((dnorm(Amatrix-DATmatrix,sd=h[i])-dnorm(Bmatrix-DATmatrix,sd=h[i]))%*%matrix(1,n,1))*matrix(k,8,1))*(h[i]^2)/n
    M2=(((pnorm((Bmatrix-DATmatrix)/h[i])-pnorm((Amatrix-DATmatrix)/h[i]))*(Kmatrix*DATmatrix+Cmatrix))%*%matrix(1,n,1))/n
    mixterm=sum(M1+M2)
    ISEfstar_vec[i]=R.f+R.fhat-2*mixterm
  }
  ISEfstar_vec
}


#' The OSCV function based on the kernel \code{\link{L_I}} in the density estimation (KDE) context.
#'
#' Computing the values of the \eqn{L_I}-based OSCV function in the density estimation context. See Savchuk (2017).
#'
#' Computing the OSCV function for the given vector of bandwidth values \eqn{h} and the data vector \eqn{dat}. The function is based on the one-sided kernel \code{\link{L_I}} that depends on the parameters \eqn{\alpha} and \eqn{\sigma}. The kernel \eqn{L_I} is robust in the special case of \eqn{\alpha=16.8954588} and \eqn{\sigma=1.01}. The other special case is obtained when either of the following holds: \itemize{\item \eqn{\alpha=0} for any \eqn{\sigma>0}; \item \eqn{\sigma=1} for any \eqn{-\infty<\alpha<\infty}.} In the above cases the kernel \eqn{L_I} reduces to the one-sided Gaussian kernel \eqn{L_G}.  The function's minimizer is to be used without additional rescaling to compute the ultimate Gaussian density estimate under the assumption that the underlying density is smooth.
#' @param h numerical vector of bandwidth values,
#' @param dat numerical vecror of data values,
#' @param alpha first parameter of the kernel \eqn{L_I},
#' @param sigma second parameter of the kernel \eqn{L_I}.
#' @return The vector of values of the OSCV function for the correponsing vector of \eqn{h} values.
#' @references
#' Savchuk, O.Y. (2017). One-sided cross-validation for nonsmooth density functions, arXiv:1703.05157.
#' @seealso \code{\link{OSCV_Gauss_dens}}, \code{\link{OSCV_Epan_dens}}, \code{\link{C_smooth}}, \code{\link{L_I}}, \code{\link{H_I}}.
#' @examples
#' \dontrun{
#' # Example 1 (Old Faithful geyser data)
#' dev.new()
#' data=faithful[,1]         # Data on n=272 eruption duration of the Old Faithful geyser.
#' harray=seq(0.025,0.6,len=50)
#' alp=16.8954588
#' sig=1.01
#' plot(harray,OSCV_LI_dens(harray,data,alpha=alp,sigma=sig),lwd=3,'l',xlab="h",
#' ylab="L_I-based OSCV",main="OSCV_LI(h) for eruption duration",cex.main=1.5,cex.lab=1.7,
#' cex.axis=1.7)
#' h_OSCV_LI=round(optimize(OSCV_LI_dens,c(0.001,0.5),tol=0.001,dat=data,alpha=16.8954588,
#' sigma=1.01)$minimum,digits=4)
#' legend(0.01,-0.2,legend=c("n=272",paste("h_OSCV_LI=",h_OSCV_LI)),cex=1.8,bty="n")
#' legend(0.25,-0.33,legend=c("Parameters of L_I:", paste("alpha=",alp),
#' paste("sigma=",sig)),cex=1.7,bty="n")
#' 
#' # Example 2 (Simulated example)
#' dat_norm=rnorm(100)   #generating a random sample of size n=100 from the N(0,1) density
#' harray=seq(0.05,1.5,len=100)
#' OSCVarray=OSCV_LI_dens(harray,dat=dat_norm,16.8954588,1.01)
#' dev.new()
#' plot(harray,OSCVarray,lwd=3,'l',xlab="h",ylab="L_I-based OSCV",
#' main="OSCV_LI(h) for data generated from N(0,1)",cex.main=1.5,cex.lab=1.7,cex.axis=1.7)
#' h_OSCV_LI_norm=round(optimize(OSCV_LI_dens,c(0.001,1),tol=0.001,
#' dat=dat_norm,16.8954588,1.01)$minimum,digits=4)
#' legend(0,OSCVarray[1],legend=c("n=100",paste("h_OSCV_LI=",h_OSCV_LI_norm),
#' "Parameters of the robust kernel L_I:","alpha=16.8954588", "sigma=1.01"),cex=1.5,bty="n")
#' }
#' @export
#' @importFrom "stats" "dnorm" "pnorm" "pgamma"

OSCV_LI_dens=function(h,dat,alpha,sigma){
  Csmooth=C_smooth(alpha,sigma)
  h=h/Csmooth
  n=length(dat)
  ngrid=length(h)
  X1=matrix(dat,n,1) %*% matrix(1,1,n)
  X2=matrix(1,n,1) %*% matrix(dat,1,n)
  X=X1-X2
  denom=pi*alpha*sigma^2+2*alpha^2*sigma^2-4*alpha^2*sigma-pi*alpha+2*alpha^2-4*alpha*sigma-pi+4*alpha+2
  c1=2*pi*(alpha*sigma^2-alpha-1)/denom
  c2=-2*sqrt(2*pi)*(alpha*sigma-alpha-1)/denom
  a=c1*(1+alpha)
  b=c2*(1+alpha)
  c=-c1*alpha
  d=-c2*alpha
  OSCV=1:ngrid
  for(i in 1:ngrid){
    D=abs(X/h[i])
    I1=dnorm(D/sqrt(2))/sqrt(2)*(b^2/4*(1-pgamma(D*D/4,1.5))+sqrt(2)*a*b*dnorm(D/sqrt(2))+(a^2-b^2/2*D*D/2)*(1-pnorm(D/sqrt(2))))
    I2=dnorm(D/sqrt(sigma^2+1))*(b*d*sigma^2/(2*(sigma^2+1)^1.5)*(1-pgamma(D*D/(2*sigma^2*(sigma^2+1)),1.5))+
    dnorm(D/(sigma*sqrt(sigma^2+1)))*sigma/((sigma^2+1)^2)*((a*d+b*c)*(sigma^2+1)+b*d*D*(sigma^2-1))+
    (1-pnorm(D/(sigma*sqrt(sigma^2+1))))/((sigma^2+1)^2.5)*(a*(sigma^2+1)-b*D)*(c*(sigma^2+1)+d*D*sigma^2))
    I3=dnorm(D/sqrt(sigma^2+1))*(b*d*sigma^2/(2*(sigma^2+1)^1.5)*(1-pgamma(D*D*sigma^2/(2*(sigma^2+1)),1.5))+
    (c*(sigma^2+1)-d*D*sigma^2)*(a*(sigma^2+1)+b*D)/((sigma^2+1)^2.5)*(1-pnorm(D*sigma/sqrt(sigma^2+1)))+
    ((b*c+a*d)*sigma*(sigma^2+1)-b*d*D*sigma*(sigma^2-1))/((sigma^2+1)^2)*dnorm(D*sigma/sqrt(sigma^2+1)))
    I4=dnorm(D/(sigma*sqrt(2)))/(sigma*sqrt(2))*(sigma^2*d^2/4*(1-pgamma(D*D/(4*sigma^2),1.5))+
    (c^2-d^2*D^2/4)*(1-pnorm(D/(sigma*sqrt(2))))+sqrt(2)*c*d*sigma*dnorm(D/(sigma*sqrt(2))))
    Rf=sum(I1+I2+I3+I4)/(n*n*h[i])
    term1=(sum((a+b*D)*dnorm(D))-n*a*dnorm(0))/(n*(n-1)*h[i])
    term2=(sum((c+d*D)/sigma*dnorm(D/sigma))-n*c*dnorm(0)/sigma)/(n*(n-1)*h[i])
    OSCV[i]=Rf-(term1+term2)
  }
  h=h*Csmooth
  OSCV
}

#' The OSCV function based on \eqn{L_G}, the one-sided Gaussian kernel, in the kernel density estimation (KDE) context.
#'
#' Computing the values of the \eqn{L_G}-based OSCV function in the density estimation context. See Savchuk (2017).
#'
#' Computing the values of the OSCV function for the given bandwidth vector \eqn{h} and data vector \eqn{dat}. The function is based on the one-sided Gaussian kernel \eqn{L_G}. The (anticipated) smoothness of the underlying density function is to be specified. Thus, \itemize{\item \eqn{stype=0} corresponds to the smooth density; \item \eqn{stype=1} corresponds to the nonsmooth density.} It is usually assumed that the density is smooth if no preliminary information about its nonsmoothness is available. The function's minimizer \code{\link{h_OSCV_dens}} is to be used without additional rescaling to compute the ultimate Gaussian density estimate.
#' @param h numerical vector of bandwidth values,
#' @param dat numerical vecror of data values,
#' @param stype specifies (anticipated) smoothness of the density function. Thus, \eqn{stype=0} corresponds to the \emph{smooth} density, whereas \eqn{stype=1} corresponds to the \emph{nonsmooth} density.
#' @return The vector of values of the OSCV function for the correponsing vector of \eqn{h} values.
#' @references
#' Savchuk, O.Y. (2017). One-sided cross-validation for nonsmooth densty functions, arXiv:1703.05157.
#' @seealso  \code{\link{h_OSCV_dens}}, \code{\link{OSCV_Epan_dens}}, \code{\link{OSCV_LI_dens}}, \code{\link{C_smooth}}.
#' @examples
#' \dontrun{
#' dat_norm=rnorm(300)   #generating random sample of size n=300 from the standard normal density.
#' h_oscv=round(h_OSCV_dens(dat_norm,0),digits=4)
#' y=density(dat_norm,bw=h_oscv)
#' dev.new()
#' plot(y,lwd=3,cex.lab=1.7,cex.axis=1.7,cex.main=1.7,xlab=paste("n=100, h_OSCV=",h_oscv),
#' main="Standard normal density estimate by OSCV",ylim=c(0,0.45),xlim=c(-4.5,4.5))
#' u=seq(-5,5,len=1000)
#' lines(u,dnorm(u),lwd=3,lty="dashed",col="blue")
#' legend(0.75,0.4,legend=c("OSCV estimate","N(0,1) density"),lwd=c(3,3),lty=c("solid","dashed"),
#' col=c("black","blue"),bty="n",cex=1.25)
#' }
#' @export
#' @importFrom "stats" "dnorm" "pnorm" "pgamma"
OSCV_Gauss_dens=function(h,dat,stype){
  a=2*pi/(pi-2)
  b=-2*sqrt(2*pi)/(pi-2)
  n=length(dat)
  ngrid=length(h)
  X1=matrix(dat,n,1) %*% matrix(1,1,n)
  X2=matrix(1,n,1) %*% matrix(dat,1,n)
  X=X1-X2
  if (stype==0) Const=0.6168471
  else if (stype==1) Const=0.5730
  else print("Set stype=0 (smooth case) or stype=1 (nonsmooth case)")
  if((stype==0)|(stype==1)){
    h=h/Const
    OSCV=1:ngrid
    for(i in 1:ngrid){
      B=abs(X/(sqrt(2)*h[i]))
      M1=(b*b/4)*(1-pgamma(B*B/2,1.5))
      OSCV[i]=sum(dnorm(B)*(M1+sqrt(2)*a*b*dnorm(B)+(1-pnorm(B))*(a*a-(b*b/2)*B*B)))
      OSCV[i]=OSCV[i]/(sqrt(2)*n*n*h[i])
      term=(sum((a+b*abs(X/h[i]))*dnorm(X/h[i]))-n*a*dnorm(0))/2
      OSCV[i]=OSCV[i]-2*term/(n*(n-1)*h[i])
    }
    OSCV
  }
}

#' The OSCV bandwidth in the density estimation context.
#'
#' Computing the OSCV bandwidth for the Gaussian density estimator. The one-sided Gaussian kernel \eqn{L_G} is used in the bandwidth selection stage. The (anticipated) smoothness of the density function is to be specified by the user.
#'
#' Computing the OSCV bandwidth for the data vector \eqn{dat}. The one-sided Gaussian kernel \eqn{L_G} is used for the cross-validation purposes and the Gaussian kernel is used for computing the ultimate density estimate. The (anticipated) smoothness of the underlying density function is to be specified. Thus, \itemize{\item \eqn{stype=0} corresponds to the smooth density; \item \eqn{stype=1} corresponds to the nonsmooth density.} It is usually assumed that the density is smooth if no preliminary information about its nonsmoothness is available. No additional rescaling of the computed bandwidth is needed. The smoothness of the density function \eqn{stype}, essentially, determines the value of the bandwidth rescaling constant that is used in the body of the function. Thus, the constant is equal to 0.6168471 in the smooth case, whereas it is equal to 0.5730 in the nonsmooth case. See Savchuk (2017) for details. The OSCV bandwidth is the minimizer of the OSCV function \code{\link{OSCV_Gauss_dens}}.
#' @param dat numerical vecror of data values,
#' @param stype specifies (anticipated) smoothness of the density function. Thus, \eqn{stype=0} corresponds to the \emph{smooth} density, whereas \eqn{stype=1} corresponds to the \emph{nonsmooth} density.
#' @return The OSCV bandwidth (scalar).
#' @references
#' Savchuk, O.Y. (2017). One-sided cross-validation for nonsmooth densty functions, arXiv:1703.05157.
#' @seealso \code{\link{OSCV_Gauss_dens}},  \code{\link{C_smooth}}, \code{\link{h_OSCV_reg}}.
#' @examples
#' \dontrun{
#' data=faithful[,1]         # Data on n=272 eruption duration of the Old Faithful geyser.
#' harray=seq(0.025,0.6,len=100)
#' OSCV_array=OSCV_Gauss_dens(harray,data,0)
#' dev.new()
#' plot(harray,OSCV_array,lwd=3,'l',xlab="h",ylab="L_G-based OSCV",
#' main="OSCV_G(h) for the data on eruption duration",cex.main=1.5,cex.lab=1.7,cex.axis=1.7)
#' h_oscv=round(h_OSCV_dens(data,0),digits=4) #smoothness of the underlying density is assumed
#' legend(0.04,-0.25,legend=c("n=272",paste("h_OSCV=",h_oscv)),cex=2,bty="n")
#'}
#' @export
#' @importFrom "stats" "optimize"
h_OSCV_dens=function(dat,stype){
  R=range(dat)
  up=(R[2]-R[1])/4
  optimize(OSCV_Gauss_dens,c(0.001,up),tol=0.001,dat=dat,stype=stype)$minimum
}


#' The OSCV function based on \eqn{L_E}, the one-sided Epanechnikov kernel, in the kernel density estimation (KDE) context.
#'
#' Computing the values of the \eqn{L_E}-based OSCV function in the density estimation context. See Martinez-Miranda et al. (2009) and Savchuk (2017).
#'
#' Computing the values of the OSCV function for the given bandwidth vector \eqn{h} and data vector \eqn{dat}. The function is based on the one-sided Epanechnikov kernel \eqn{L_E}. The function's minimizer is to be multiplied by the appropriate rescaling constant before it can be used to compute the ultimate kernel density estimate. The formula for the rescaling constant depends on \emph{smothness} of the density and on the \emph{kernel} used in computing the ultimate density estimate.
#' @param h numerical vector of bandwidth values,
#' @param dat numerical vecror of data values.
#' @return The vector of values of the OSCV function for the correponsing vector of \eqn{h} values.
#' @references
#' \itemize{
#'   \item Martinez-Miranda, M.D., Nielsen, J. P., and Sperlich, S. (2009). One sided cross validation for density estimation. In \emph{Operational Risk Towards Basel III: Best Practices and Issues in Modeling, Management and Regulation}, 177-196.
#'   \item Savchuk, O.Y. (2017). One-sided cross-validation for nonsmooth densty functions, arXiv:1703.05157.
#' }
#' @seealso  \code{\link{OSCV_Gauss_dens}}, \code{\link{OSCV_LI_dens}}.
#' @examples
#' \dontrun{
#' # Example 1 (Data on n=272 eruption duration of the Old Faithful geyser).
#' data=faithful[,1]
#' har=seq(0.05,1,len=1000)
#' dev.new()
#' plot(har,OSCV_Epan_dens(har,data),lwd=3,'l',xlab="h",ylab="L_E-based OSCV",
#' main="L_E_based OSCV for the data on eruption duration",cex.main=1.5,cex.lab=1.7,cex.axis=1.7)
#' h_min=round(optimize(OSCV_Epan_dens,c(0.001,1),tol=0.001,dat=data)$minimum, digits=4)
#' legend(0.1,-0.1,legend=c("n=272",paste("h_min=",h_min)),cex=2)
#' # The above graph appears in Savchuk (2017).
#'
#' # Example 2 (Data set of size n=100 is generated from the standard normal density).
#' dat_norm=rnorm(100)
#' harray=seq(0.25,4.25,len=1000)
#' OSCVarray=OSCV_Epan_dens(harray,dat_norm)
#' dev.new()
#' plot(harray,OSCVarray,lwd=3,'l',xlab="h",ylab="L_E-based OSCV",
#' main="L_E-based OSCV for data generated from N(0,1)", cex.main=1.5,cex.lab=1.7,cex.axis=1.7)
#' h_min_norm=round(optimize(OSCV_Epan_dens,c(0.1,4),tol=0.001,dat=dat_norm)$minimum, digits=4)
#' legend(0.5,OSCVarray[1],legend=c("n=100",paste("h_min=",h_min_norm)),cex=2,bty="n")
#' }
#' @export
#' @importFrom "stats" "rnorm"
OSCV_Epan_dens=function(h,dat){
  ngrid=length(h)
  n=length(dat)
  X1=matrix(dat,n,1) %*% matrix(1,1,n)
  X2=matrix(1,n,1) %*% matrix(dat,1,n)
  OSCV=1:ngrid
  for(i in 1:ngrid){
    D=abs((X1-X2)/h[i])
    I=(D<1)
    DI=D[I]
    Rf=(sum(-38/3*DI^2-32*DI+45/28*DI^7-739/30*DI^5+677/12*DI^3+1184/105))*144/(361*n^2*h[i])
    term=(sum((8-15*DI)*(1-DI^2)))*12/(19*h[i]*n*(n-1))-96/(19*h[i]*(n-1))
    OSCV[i]=Rf-term
  }
  OSCV
}