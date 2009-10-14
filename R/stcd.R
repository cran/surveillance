######################################################################
# Shiryaev-Roberts based spatio-temporal cluster detection based
# on the work in Assuncao & Correa (2009). The implementation
# is based on C++ code was originally written by Marcos Oliveira Prates, UMFG,
# Brazil and provided by Thais Correa, UMFG, Brazil during her research
# stay in Munich. This stay was financially supported by the Munich
# Center of Health Sciences.
#
#
# Parameters:
#   x - vector containing spatial x coordinate of the events
#   y - vector containing spatial y coordinate of the events
#   t - vector containing the time points of the events
#   radius - is the radius of the cluster 
#   epsilon - is the relative change of event-intensity within the cluster
#       to detect
#   A - threshold limit for the alarm and should be equal to the desired ARL
######################################################################


stcd <- function(x, y,t,radius,epsilon,threshold) {
  #check that x,y,t are of the same length.
  n <- length(x)
  if ((length(y) != n) | (length(t) != n)) {
    stop("Vectors x,y,t not of same size.")
  }

  res <- .C("SRspacetime", x=as.double(x), y=as.double(y), t=as.double(t), n=as.integer(n), radius=as.double(radius), epsilon=as.double(epsilon), threshold=as.double(threshold),R=as.double(numeric(n)),idxFA=as.integer(-1),idxCC=as.integer(-1),PACKAGE="surveillance")

  #Missing: compute which indices are part of the cluster.
  #Thais
  
  return(list(R=res$R,idxFA=res$idxFA+1,idxCC=res$idxCC+1))
}
