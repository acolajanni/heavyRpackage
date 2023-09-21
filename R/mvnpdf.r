################################################################################
# > 20-22 Septembre 2023 - HeavyR
# > Script :
# > Function :
# @ COLAJANNI Antonin
################################################################################



# determinant = det(varcovM)
# solve(varcovM) #solve : inverse de matrice
# transposed_mat = t(varcovM) # inverse row/column
# ident_mat = varcovM %*% inverse_mat | produit matriciel : %*%
#t(x) %*% solve(varcovM) %*% x



#' density of a multivariate normal distribution for on column vector
#'
#' Description
#'
#' Details (opt)
#'
#' @param x matrix : one column
#' @param mu mean (vector of mean)
#' @param varcovM matrix - square matrix (nrow == ncol) strictly positive
#' @param compute_mat Boolean - True to compute inverse and det
#' @param determ if !compute_mat determinant of varcovM (output of det(varcovM))
#' @param inv if !compute_mat inverse of varcovM (output of solve(varcovM))
#'
#' @return vector
gaussian_singleValue=function(x,mu,varcovM , compute_mat = TRUE, determ=NULL, inv=NULL) {
  p = nrow(varcovM)
  if(!compute_mat){
    return( (2*pi)^(-p/2) * det(varcovM)^(-0.5) * exp( -0.5 * t(x - mu) %*% solve(varcovM) %*% (x-mu) ) )  }
  else{
    return( (2*pi)^(-p/2) * determ^(-0.5) * exp( -0.5 * t(x - mu) %*% inv %*% (x-mu) ) )  }
}

#' density of a multivariate normal distribution on Rp at n points
#'
#' @param x matrice, observations: point pour lequel on
#' veut calculer la densité de loi normale
#' @param mu mean (mean vector of a gaussian distribution)
#' @param varcovM matrix - square matrix (nrow == ncol) strictly positive variance covariance matrix
#' @param Log Bool - return log of y or not
#'
#' @return vector of size nrow(x) and x
#' @export
#'
#' @examples
#' mat = matrix(c(1,0,0,1), nrow=2, byrow = TRUE)
#' x=matrix(c(-0.5,1.5,0,1,-1, 1),nrow = 2)
#' mu = matrix(c(0,0),nrow = 2)
#' mvnpdf(x=x, varcovM = mat, mu=mu, Log=FALSE)
mvnpdf = function(x ,  varcovM = diag(nrow(x)), mu = rep(0,nrow(varcovM)), Log=FALSE){

  determ = det(varcovM)
  inv = solve(varcovM)
  n = ncol(x)
  v <- rep(NA, n)
  for (i in 1:n) {
    v[i] = gaussian_singleValue(x[,i,drop=FALSE], mu, varcovM ,compute_mat=FALSE, #[x,y, autres arguments] [] ==> functions
                                  determ=determ, inv=inv)
  }
  if (Log){v=log(v)}
  res=list("x" = x, "y"= v )
  class(res) = "mvnpdf"
  return(res) }


#' Plot of the mvnpdf function
#'
#' @param x an object of class \code{mvnpdf} resulting from a call of
#' \code{mnvpdf()} function.
#' @param ... graphical parameters passed to \code{plot()} function.
#'
#' @return Nothing is returned, only a plot is given.
#' @export
#'
#' @examples
#' pdfvalues <- mvnpdf(x=matrix(seq(-3, 3, by = 0.1), nrow = 1), Log=FALSE)
#' plot(pdfvalues)
plot.mvnpdf <- function(x, ...) {
  plot(x$x, x$y, type = "l", ...)
}

#' @export
mvnpdfoptim <- function(x, mean =  rep(0, nrow(x)),
                        varcovM = diag(nrow(x)), Log=TRUE){

  if(!is.matrix(x)){
    x <- matrix(x, ncol=1)
  }

  n <- ncol(x)
  p <- nrow(x)
  x0 <- x-mean

  Rinv <- backsolve(chol(varcovM), x=diag(p))
  xRinv <- apply(X=x0, MARGIN=2, FUN=crossprod, y=Rinv)
  logSqrtDetvarcovM <- sum(log(diag(Rinv)))

  quadform <- apply(X=xRinv, MARGIN=2, FUN=crossprod)
  y <- (-0.5*quadform + logSqrtDetvarcovM -p*log(2*pi)/2)

  if(!Log){
    y <- exp(y)
  }

  res <- list(x = x, y = y)
  class(res) <- "mvnpdf"
  return(res)
}
