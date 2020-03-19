## convert between sin/cos and amplitude/shift formulation
###################################################
# y = gamma*sin(omega*t)+delta*cos(omega*t)
#   =  A*sin(omega*t + phi)
# with Amplitude A= sqrt(gamma^2+delta^2)
# and shift phi= arctan(delta/gamma)
#################################################
sinCos2amplitudeShift <- function(params){
  # number of sin+cos terms
  lengthParams <- length(params)
  if(lengthParams %% 2 != 0)
    stop("wrong number of params")
  index.sin <- seq(1,lengthParams,by=2)

  one <- function(i=1){
    coef.sin <- params[i]
    coef.cos <- params[i+1]

    amplitude <- sqrt(coef.cos^2+coef.sin^2)
    shift <- atan2(coef.cos, coef.sin)
    return(c(amplitude,shift))
  }
  return(c(sapply(index.sin,one)))
}

amplitudeShift2sinCos <- function(params){
    lengthParams <- length(params)
    if (lengthParams %% 2 != 0)
        stop("wrong number of params")
    index.A <- seq(1, lengthParams, by = 2)
    one <- function(i = 1) {
        coef.A <- params[i]
        coef.shift <- params[i + 1]
        coef.cos <- -coef.A*tan(coef.shift)/sqrt(1+tan(coef.shift)^2)
        coef.sin <- -coef.A/sqrt(1+tan(coef.shift)^2)
        return(c(coef.sin,coef.cos))
    }
    return(c(sapply(index.A, one)))

}

##############################################
# y = gamma*sin(omega*t)+delta*cos(omega*t)
# g(gamma,delta) = [sqrt(gamma^2+delta^2), arctan(delta/gamma) ]'
# compute jacobian (dg_i(x)/dx_j)_ij
#############################################
jacobianAmplitudeShift <- function(params){
  # number of sin+cos terms
  lengthParams <- length(params)
  if(lengthParams %% 2 != 0)
    stop("wrong number of params")
  index.sin <- seq(1,lengthParams,by=2)
  # function to compute jacobian of the transformation sinCos2AmplitudeShift()
  one <- function(i=1){
    coef.sin <- params[i]
    coef.cos <- params[i+1]

    dAmplitude.dcoef.sin <- coef.sin/sqrt(coef.cos^2+coef.sin^2)
    dAmplitude.dcoef.cos <- coef.cos/sqrt(coef.cos^2+coef.sin^2)

    dShift.dcoef.sin <- - coef.cos/(coef.cos^2+coef.sin^2)
    dShift.dcoef.cos <- coef.sin/(coef.cos^2+coef.sin^2)
    return(c(dAmplitude.dcoef.sin,dShift.dcoef.sin,dAmplitude.dcoef.cos,dShift.dcoef.cos))
  }
  jacobi<-sapply(index.sin,one)
  res <- matrix(0,nrow=lengthParams,ncol=lengthParams)
  j<-0
  for (i in index.sin){
    j<-j+1
    res[i:(i+1),i:(i+1)] <- jacobi[,j]
  }
  return(res)
}
