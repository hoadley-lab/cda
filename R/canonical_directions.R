#' Calculate canonical parallel and orthogonal directions
#'
#' Implements the canonical orthogonal and parallel directions
#' from Chapter 2 of Xuxin Liu's 2007 dissertation
#'
#' @param X a p-feature by n-sample \code{matrix}
#' @param Y a p-feature by n-sample \code{matrix},
#' features and samples must align with X
#' @param COD Specifies whether to calculate the COD or not
#' @param method 'eigen' uses slow R matrix multiplication to get the
#' covariance matrix. 'comp' uses the considerably faster prcomp method
#' from irlba (most users will want this option)
#' @returns a list containing:
#' \item{directions}{p by 1 vector of canonical directions}
#' \item{projections}{n by 1 vectors of projections of the data
#'   onto the canonical directions. Example: proj_vcpd_X is the
#'   projection of X onto the CPD direction}
#' @importFrom irlba partial_eigen prcomp_irlba
#' @export

getCanonicalDirections <- function(X, Y, COD = TRUE, method = 'comp') {
  if (method == 'eigen') {
    cov_diff <- (X - Y) %*% t(X - Y)
    v_cpd <- irlba::partial_eigen(cov_diff, n = 1, symmetric = T)$vectors
  } else if (method == "comp") {
    v_cpd <- irlba::prcomp_irlba(t(X - Y), n = 1, center = F)$rotation
    v_cpd <- unname(v_cpd)
  } else {
    stop("Not a valid method ", method)
  }
  rownames(v_cpd) <- rownames(X)
  proj_vcpd_X <- t(X) %*% v_cpd
  proj_vcpd_Y <- t(Y) %*% v_cpd

  if (COD) {
    ## Calculate v_cod
    ## QR decomp of [X-Y, Y]
    N <- ncol(Y)
    H <- cbind(X - Y, Y)
    QR <- qr(H)

    Q2 <- qr.Q(QR)[, seq(N + 1, 2 * N)]

    X_bar <- rowMeans(X)
    c1 <- t(Q2) %*% (X - X_bar) %*% t(X - X_bar) %*% Q2

    C <- irlba::partial_eigen(c1, n = 1, symmetric = T)$vectors

    v_cod <- Q2 %*% C
    rownames(v_cod) <- rownames(X)
    proj_vcod_X <- t(X) %*% v_cod
    proj_vcod_Y <- t(Y) %*% v_cod
  } else {
    v_cod = NA
    proj_vcod_X <- NA
    proj_vcod_Y <- NA
  }

  return(list(
    directions = list(v_cpd = v_cpd, v_cod = v_cod),
    projections = list(
      proj_vcpd_X = proj_vcpd_X,
      proj_vcod_X = proj_vcod_X,
      proj_vcpd_Y = proj_vcpd_Y,
      proj_vcod_Y = proj_vcod_Y
    )
  ))

}


