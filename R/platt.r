##' Simulated Photosynthesis-Irradiance data
##'
##' A simulated dataset consisting of six Photosynthesis-Irradiance
##' relationships.
##'
##' @format A data frame with 124 rows and 3 variables:
##' \describe{
##'   \item{Depth}{Sample identifier}
##'   \item{P}{Photosynthetic rate}
##'   \item{I}{Irradiance}
##' }
"PI"


##' A \code{selfStart} model for Platt's Photosynthesis-Irradiance
##' curve with Intercept.
##'
##' This model differs from Platt (1980) in that it allows for a
##' non-zero intercept.
##' \deqn{P = Pmax (1-\exp(-\alpha*I/Pmax)) \exp(-\beta*I/Pmax) + R}
##'
##' @title Self Starting Modified Platt Model
##' @param I The measured irradiance
##' @param alpha initial slope of the light limited of the P-I curve
##' @param beta the rate of photoinhibition
##' @param Pmax light saturated photosynthetic rate
##' @param R intercept of the P-I curve
##' @seealso \code{\link{SSPlatt0}}
##' @references
##' Platt, T.G.C.L., Gallegos, C.L. and Harrison,
##'   W.G. (1980). Photoinhibition of photosynthesis in natural
##'   assemblages of marine phytoplankton. Journal of Marine Research
##'   (USA).
##'
##' Westwood, K.J., Griffiths, F.B., Meiners, K.M. and Williams,
##'   G.D. (2010). Primary productivity off the Antarctic coast from
##'   30-80 E; BROKE-West survey, 2006. Deep Sea Research Part II:
##'   Topical Studies in Oceanography, 57(9), pp.794-814.
##' @importFrom stats asOneSidedFormula coef setNames lsfit nls
##' @export
SSPlatt <- stats::selfStart(
  ~ Pmax*(1-exp(-alpha*I/Pmax))*exp(-beta*I/Pmax)+R,
  function(mCall,data,LHS,...) {
    ## Extract x and y but do not average replicated x values
    x <- mCall[["I"]]
    y <- LHS
    if (is.language(x) || ((length(x) == 1L) && is.character(x)))
      x <- eval(asOneSidedFormula(x)[[2L]], data)
    x <- as.numeric(x)
    if (is.language(y) || ((length(y) == 1L) && is.character(y)))
      y <- eval(asOneSidedFormula(y)[[2L]], data)
    y <- as.numeric(y)
    keep <- !is.na(x) & !is.na(y)
    x <- x[keep]
    y <- y[keep]
    if(length(unique(x)) < 6) {
      stop("Too few distinct x values to fit a self starting Platt model")
    }
    ord <- order(x)
    x <- x[ord]
    y <- y[ord]
    ## Get initial estimate of smaller rate parameter
    ks <- (length(x)+1)-seq_len(max(4,which(rev(y) >= 0.9*max(y))[1]))
    beta <- -lsfit(x[ks],log(pmax(y[ks]-min(y),1.0E-8)))$coef[2]
    beta <- max(beta,1.0E-8)
    ## Fit partial linear model to estimate second rate parameter
    fit <- nls(y ~ cbind(exp(-x*exp(b))-exp(-x*exp(a)),1),
               data = list(y=y,x=x,b=log(beta)),
               start = list(a=log(beta)+3),
               control=nls.control(maxiter=100,minFactor=1/4096),
               algorithm = "plinear")
    cf <- coef(fit)
    ## Ensure the difference in exponetials is postive
    cf <- if(cf[2]>0) c(cf[1],log(beta),cf[2:3]) else c(log(beta),cf[1],-cf[2:3])
    cf0 <- cf
    ## Refit partial linear model to estimate both parameters
    fit <- tryCatch(nls(y ~ cbind(exp(-x*exp(b))-exp(-x*exp(a)),1),
                        data = list(y=y,x=x),
                        start = list(a=cf[1],b=cf[2]),
                        control=nls.control(maxiter=100,minFactor=1/4096),
                        algorithm = "plinear"),
                    error=function(e) NULL)
    if(!is.null(fit)) {
      cf1 <- coef(fit)
      ## Ensure the difference in exponetials is postive
      cf1 <- if(cf1[3]>0) cf1 else c(cf1[2:1],-cf1[3:4])
      cf <- if(cf1[1] < cf1[2]) cf else cf1
    }
    setNames(c((exp(cf[1])-exp(cf[2]))*cf[3],exp(cf[2])*cf[3],cf[3:4]),
             c("alpha","beta","Pmax","R"))
  },
  c("alpha","beta","Pmax","R"))


##' A \code{selfStart} model for Platt's Photosynthesis-Irradiance
##' curve without intercept.
##'
##' This model is the model decribed by Platt (1980)
##' \deqn{P = Pmax (1-\exp(-\alpha*I/Pmax)) \exp(-\beta*I/Pmax)}
##'
##' @title Self Starting Modified Platt Model
##' @param I The measured irradiance
##' @param alpha initial slope of the light limited of the P-I curve
##' @param beta the rate of photoinhibition
##' @param Pmax light saturated photosynthetic rate
##' @seealso \code{\link{SSPlatt}}
##' @references
##' Platt, T.G.C.L., Gallegos, C.L. and Harrison,
##'   W.G. (1980). Photoinhibition of photosynthesis in natural
##'   assemblages of marine phytoplankton. Journal of Marine Research
##'   (USA).
##'
##' Westwood, K.J., Griffiths, F.B., Meiners, K.M. and Williams,
##'   G.D. (2010). Primary productivity off the Antarctic coast from
##'   30-80 E; BROKE-West survey, 2006. Deep Sea Research Part II:
##'   Topical Studies in Oceanography, 57(9), pp.794-814.
##' @importFrom stats asOneSidedFormula coef setNames lsfit nls
##' @export
SSPlatt0 <- stats::selfStart(
  ~ Pmax*(1-exp(-alpha*I/Pmax))*exp(-beta*I/Pmax),
  function(mCall,data,LHS,...) {
    ## Extract x and y but do not average replicated x values
    x <- mCall[["I"]]
    y <- LHS
    if (is.language(x) || ((length(x) == 1L) && is.character(x)))
      x <- eval(asOneSidedFormula(x)[[2L]], data)
    x <- as.numeric(x)
    if (is.language(y) || ((length(y) == 1L) && is.character(y)))
      y <- eval(asOneSidedFormula(y)[[2L]], data)
    y <- as.numeric(y)
    keep <- !is.na(x) & !is.na(y)
    x <- x[keep]
    y <- y[keep]
    if(length(unique(x)) < 6) {
      stop("Too few distinct x values to fit a self starting Platt model")
    }
    ord <- order(x)
    x <- x[ord]
    y <- y[ord]
    ## Get initial estimate of smaller rate parameter
    ks <- (length(x)+1)-seq_len(max(4,which(rev(y) >= 0.9*max(y))[1]))
    beta <- -lsfit(x[ks],log(pmax(y[ks]-min(y),1.0E-8)))$coef[2]
    beta <- max(beta,1.0E-8)
    ## Fit partial linear model to estimate second rate parameter
    fit <- nls(y ~ exp(-x*exp(b))-exp(-x*exp(a)),
               data = list(y=y,x=x,b=log(beta)),
               start = list(a=log(beta)+3),
               control=nls.control(maxiter=100,minFactor=1/4096),
               algorithm = "plinear")
    cf <- coef(fit)
    ## Ensure the difference in exponetials is postive
    cf <- if(cf[2]>0) c(cf[1],log(beta),cf[2:3]) else c(log(beta),cf[1],-cf[2:3])
    cf0 <- cf
    ## Refit partial linear model to estimate both parameters
    fit <- tryCatch(nls(y ~ exp(-x*exp(b))-exp(-x*exp(a)),
                        data = list(y=y,x=x),
                        start = list(a=cf[1],b=cf[2]),
                        control=nls.control(maxiter=100,minFactor=1/4096),
                        algorithm = "plinear"),
                    error=function(e) NULL)
    if(!is.null(fit)) {
      cf1 <- coef(fit)
      ## Ensure the difference in exponetials is postive
      cf1 <- if(cf1[3]>0) cf1 else c(cf1[2:1],-cf1[3:4])
      cf <- if(cf1[1] < cf1[2]) cf else cf1
    }
    setNames(c((exp(cf[1])-exp(cf[2]))*cf[3],exp(cf[2])*cf[3],cf[3:4]),
             c("alpha","beta","Pmax","R"))
  },
  c("alpha","beta","Pmax"))



##' Estimate respiration corrected Pmax from a fit of a Platt model.
##'
##' Given coefficients from a fitted \code{\link{SSPlatt}} model
##' fitted with \code{nls}, compute Pmax corrected for respiration.
##' @title Respiration Corrected Pmax
##' @param alpha fitted coefficient
##' @param beta fitted coefficient
##' @param Pmax fitted coefficient
##' @return The corrected Pmax
##' @export
correctedPmax <- function(alpha,beta,Pmax) {
  Pmax*(alpha/(alpha+beta))*((beta/alpha+beta))^(beta/alpha)
}
