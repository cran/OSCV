#' The OSCV smooth rescaling constant.
#'
#' Computing the OSCV smooth rescaling constant that corresponds to using the two-sided kernel \code{\link{H_I}} for the cross-validation purposes and the Gaussian kernel in the estimation stage. The constant is applicable for the OSCV versions in the regression and kernel density estimation contexts.
#'
#' Computation of the OSCV rescaling constant \eqn{C} (see (10) in Savchuk and Hart (2017) or (3) in Savchuk (2017)). The constant is a function of the parameters \eqn{(\alpha,\sigma)} of the two-sided cross-validation kernel \code{\link{H_I}} defined by expression (15) in Savchuk and Hart (2017). The Gaussian kernel is used for computing the ultimate (regression or density) estimate. The constant is used in the OSCV versions for kernel regression and density estimation. Notice that in the cases \eqn{\alpha=0}, \eqn{\sigma>0} and \eqn{\sigma=1}, \eqn{-\infty<\alpha<\infty} the kernel \code{\link{H_I}} reduces to the Gaussian kernel.
#' @param alpha first parameter of the two-sided cross-validation kernel \code{\link{H_I}},
#' @param sigma second parameter of the two-sided cross-validation kernel \code{\link{H_I}}.
#' @return The OSCV smooth rescaling constant \eqn{C} for the given values of the parameters \eqn{\alpha} and \eqn{\sigma}.
#' @references
#' \itemize{
#'   \item Savchuk, O.Y., Hart, J.D. (2017). Fully robust one-sided cross-validation for regression functions. \emph{Computational Statistics}, doi:10.1007/s00180-017-0713-7.
#'   \item Savchuk, O.Y. (2017). One-sided cross-validation for nonsmooth density functions, arXiv:1703.05157.
#' }
#' @seealso \code{\link{L_I}}, \code{\link{H_I}}, \code{\link{OSCV_reg}}, \code{\link{h_OSCV_reg}}, \code{\link{OSCV_LI_dens}}, \code{\link{OSCV_Gauss_dens}}, \code{\link{h_OSCV_dens}}, \code{\link{loclin}}.
#' @examples
#' # OSCV rescaling constant for the robust cross-validation kernel with 
#' # (alpha,sigma)=(16.8954588,1.01).
#' C_smooth(16.8954588,1.01)
#' # OSCV smooth rescaling constant in the case when the kernel H_I is Gaussian.
#' C_smooth(1,1)
#' @export

C_smooth=function(alpha,sigma){
  N=(pi-4)*(alpha^2*sigma^4+alpha^2+2*alpha+1)-2*alpha*sigma*(alpha+1)*(pi*sigma-2*sigma^2-2)
  D=pi*alpha*(sigma^2-1)+2*alpha^2*(sigma^2+1)-4*alpha*sigma*(alpha+1)-pi+4*alpha+2
  mu2_L_I=-N/D
  sr=sqrt(sigma^2+1)
  N_R=12*sqrt(2)*alpha^2*sigma^3+4*sqrt(2)*alpha*sigma^3-8*sqrt(2)*alpha^2*sigma^4+
  sr*(-7*alpha^2*sigma^3-4*alpha^3*sigma+5*alpha^2*sigma^2-6*alpha^2*sigma+
  2*alpha*sigma^2-4*alpha*sigma-4*alpha*sigma^3-alpha^2*sigma^5+5*alpha^2*sigma^4+2*alpha*sigma^4-alpha^4*sigma^6+
  4*alpha^3*sigma^4-4*alpha^3*sigma^3+alpha^4*sigma^5+alpha^4*sigma^2-alpha^4*sigma+4*alpha^3*sigma^2-pi*alpha^4+
  2*sqrt(2)*sigma^3-2*pi*alpha^3-pi*sigma^3-pi*alpha^2+2*sqrt(2)*sigma-pi*sigma-sigma^3-
  sigma)+2*sqrt(2)*pi*(alpha^4*sigma^7+alpha^3*sigma^7-alpha^4*sigma^5-3*alpha^3*sigma^5-alpha^4*sigma^3-alpha^3*sigma^3+alpha^4*sigma+
  3*alpha^3*sigma)+4*sqrt(2)*alpha^3*sigma^3*(alpha*sigma^2+sigma^2+alpha+3-2*alpha*sigma-4*sigma)+
  sr*(pi*alpha^4*sigma^5-pi*alpha^2*sigma^7+8*sqrt(2)*alpha^4*sigma^4-16*sqrt(2)*alpha^3*sigma^5+2*sqrt(2)*alpha^2*sigma^6+pi*alpha^4*sigma^4+4*pi*alpha^3*sigma^5+
  8*sqrt(2)*alpha^4*sigma^3+12*sqrt(2)*alpha^3*sigma^4-6*sqrt(2)*alpha^2*sigma^5+pi*alpha^4*sigma^3+2*pi*alpha^3*sigma^4+5*pi*alpha^2*sigma^5-
  12*sqrt(2)*alpha^4*sigma^2+20*sqrt(2)*alpha^3*sigma^3+4*sqrt(2)*alpha^2*sigma^4-2*sqrt(2)*alpha*sigma^5+pi*alpha^4*sigma^2+
  2*pi*alpha^3*sigma^3+2*pi*alpha*sigma^5+4*sqrt(2)*alpha^4*sigma-32*sqrt(2)*alpha^3*sigma^2+16*sqrt(2)*alpha^2*sigma^3-
  2*sqrt(2)*alpha*sigma^4-pi*alpha^4*sigma+12*sqrt(2)*alpha^3*sigma-30*sqrt(2)*alpha^2*sigma^2+6*sqrt(2)*alpha*sigma^3-4*pi*alpha^3*sigma-
  pi*alpha^2*sigma^2-2*pi*alpha*sigma^3+14*sqrt(2)*alpha^2*sigma-10*sqrt(2)*alpha*sigma^2-6*pi*alpha^2*sigma+8*sqrt(2)*alpha*sigma-
  4*pi*alpha*sigma-pi*alpha^4*sigma^7+4*sqrt(2)*alpha^4*sigma^6-pi*alpha^4*sigma^6-2*pi*alpha^3*sigma^7-12*sqrt(2)*alpha^4*sigma^5+
  4*sqrt(2)*alpha^3*sigma^6)-4*sqrt(2)*pi*alpha^2*sigma^5+6*sqrt(2)*pi*alpha^2*sigma+2*sqrt(2)*pi*alpha*sigma+2*sqrt(2)*pi*alpha*sigma^3+2*sqrt(2)*pi*alpha^2*sigma^3
  D_R=sigma*(sigma^2+1)^(3/2)*(pi^2*alpha^2*sigma^4+4*pi*alpha^3*sigma^4+4*alpha^4*sigma^4-8*pi*alpha^3*sigma^3-16*alpha^4*sigma^3-2*pi^2*alpha^2*sigma^2-
  8*pi*alpha^2*sigma^3+24*alpha^4*sigma^2-16*alpha^3*sigma^3-2*pi^2*alpha*sigma^2+8*pi*alpha^3*sigma+4*pi*alpha^2*sigma^2-16*alpha^4*sigma+48*alpha^3*sigma^2+
  pi^2*alpha^2-4*pi*alpha^3+16*pi*alpha^2*sigma+4*pi*alpha*sigma^2+4*alpha^4-48*alpha^3*sigma+24*alpha^2*sigma^2+2*pi^2*alpha-12*pi*alpha^2+8*pi*alpha*sigma+
  16*alpha^3-48*alpha^2*sigma+pi^2-12*pi*alpha+24*alpha^2-16*alpha*sigma-4*pi+16*alpha+4)
  RL_I=-N_R*sqrt(pi)/D_R
  (1/2/sqrt(pi)*mu2_L_I^2/RL_I)^(1/5)
}


#' The family of one-sided cross-validation kernels \eqn{L_I}.
#'
#' The one-sided counterpart of the kernel \code{\link{H_I}}. See expressions (15) and (8) of Savchuk and Hart (2017).
#'
#' The family of the one-sided cross-validation kernels \eqn{L_I} indexed by the parameters \eqn{-\infty<\alpha<\infty} and \eqn{\sigma>0}. This family is used in the OSCV implementations in both regression context (see Savchuk and Hart (2017)) and density estimation context (see Savchuk (2017)). The special members of the family: \itemize{\item The \emph{ robust} kernel used in Savchuk and Hart (2017) and Savchuk (2017) is obtained by setting  \eqn{\alpha=16.8954588} and \eqn{\sigma=1.01}; \item  The one-sided Gaussian kernel \eqn{L_G} is obtained by either setting \eqn{\alpha=0} for any \eqn{\sigma>0} or by setting \eqn{\sigma=1} for any \eqn{-\infty<\alpha<\infty}.} The bandwidth selected by \eqn{L_I} should be multiplied by a reascaling constant before it is used in computing the ultimate Gaussian (regression or density) estimate. In the case of a smooth (regression or density) function the rescaling constant is \code{\link{C_smooth}}.
#' @param u numerical vector of argument values,
#' @param alpha first parameter of the cross-validation kernel \eqn{L_I},
#' @param sigma second parameter of the cross-validation kernel \eqn{L_I}.
#' @return The value of \eqn{L_I(u;\alpha,\sigma)}.
#' @references
#' \itemize{
#'   \item Savchuk, O.Y., Hart, J.D. (2017). Fully robust one-sided cross-validation for regression functions. \emph{Computational Statistics}, doi:10.1007/s00180-017-0713-7.
#'   \item Savchuk, O.Y. (2017). One-sided cross-validation for nonsmooth density functions, arXiv:1703.05157.
#' }
#' @seealso \code{\link{H_I}}, \code{\link{C_smooth}}, \code{\link{OSCV_LI_dens}}.
#' @examples
#' \dontrun{
#' # Plotting the robust one-sided kernel from Savchuk and Hart (2017) with 
#' # alpha=16.8954588 and sigma=1.01.
#' u=seq(-1,5,len=1000)
#' rker=L_I(u,16.8954588,1.01)
#' Gker=L_I(u,0,1)
#' dev.new()
#' plot(u,rker,'l',lwd=3,cex.axis=1.7, cex.lab=1.7)
#' title(main="One-sided kernels: L_I (robust) and L_G",cex=1.7)
#' lines(u,Gker,lty="dashed",lwd=3)
#' legend(0.5,2.5,lty=c("solid","dashed"),lwd=c(3,3),legend=c("L_I","L_G"),cex=1.7)
#' legend(2,1.5,legend=c("alpha=16.8955","sigma=1.01"),cex=1.5)
#' }
#' @export
L_I=function(u,alpha,sigma){
  S1=-1/sqrt(2*pi)*(alpha*sigma-alpha-1)
  S2=1/2*(alpha-alpha*sigma^2+1)
  (S2-S1*u)/(S2/2-S1^2)*H_I(u,alpha,sigma)*sapply(u,function(t){if (t>=0) 1 else 0})
}

