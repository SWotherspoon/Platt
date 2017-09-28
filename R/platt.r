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
##' @export
SSPlatt <- stats::selfStart(
  ~ Pmax*(1-exp(-alpha*I/Pmax))*exp(-beta*I/Pmax)+R,
  function(mCall,data,LHS) {
    xy <- sortedXyData(mCall[["I"]], LHS, data)
    if(nrow(xy) < 6) {
      stop("Too few distinct x values to fit a self starting Platt model")
    }
    y <- xy[["y"]]
    ## Estimate slope (alpha) and intercept (R) of the initial section
    ## of the curve by fitting a quadratic
    k1 <- max(4,which(y >= 0.8*max(y))[1])
    cs1 <- coef(lm(y~x+I(x^2),data=xy[1:k1,]))
    ## Estimate the rate of decay (beta) from the slope of the log
    ## transformed data after the peak.
    k2 <- nrow(xy)+1-max(4,which(rev(y) >= 0.8*max(y))[1])
    cs2 <- coef(lm(log(pmax(y,1.0E-4))~x,data=xy[k2:nrow(xy),]))
    ## Estimate Pmax from the maximum
    Pmax <- max(xy$y)-cs1[1]
    setNames(c(cs1[2],-cs2[2],Pmax,cs1[1]),c("alpha","beta","Pmax","R"))

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
##' @export
SSPlatt0 <- stats::selfStart(
  ~ Pmax*(1-exp(-alpha*I/Pmax))*exp(-beta*I/Pmax),
  function(mCall,data,LHS) {
    xy <- sortedXyData(mCall[["I"]], LHS, data)
    if(nrow(xy) < 6) {
      stop("Too few distinct x values to fit a self starting Platt model")
    }
    y <- xy[["y"]]
    ## Estimate slope (alpha) of the initial section of the curve by
    ## fitting a quadratic
    k1 <- max(4,which(y >= 0.8*max(y))[1])
    cs1 <- coef(lm(y~x+I(x^2)-1,data=xy[1:k1,]))
    ## Estimate the rate of decay (beta) from the slope of the log
    ## transformed data after the peak.
    k2 <- nrow(xy)+1-max(4,which(rev(y) >= 0.8*max(y))[1])
    cs2 <- coef(lm(log(pmax(y,1.0E-4))~x,data=xy[k2:nrow(xy),]))
    ## Estimate Pmax from the maximum
    Pmax <- max(xy$y)
    setNames(c(cs1[1],-cs2[2],Pmax),c("alpha","beta","Pmax"))

  },
  c("alpha","beta","Pmax"))




##' Estiamte respiration corrected Pmax from a fit of a Platt model.
##'
##' Given a fitted \code{\link{SSPlatt}} model fitted with \code{nls},
##' compute Pmax corrected for respiration.
##' @title Respiration Corrected Pmax
##' @param fit an SSPlatt model fitted with \code{nls}
##' @return The corrected Pmax
##' @export
correctedPmax <- function(fit) {
  with(as.list(setNames(coef(fit),c("alpha","beta","Pmax","R"))),
       Pmax*(alpha/(alpha+beta))*((beta/alpha+beta))^(beta/alpha))
}
