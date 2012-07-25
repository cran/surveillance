### --- Test setup ---
 
if(FALSE) {
  ## Not really needed, but can be handy when writing tests
  library("RUnit")
  library("surveillance")
}

##Simulate data and apply the algorithm
S <- 1 ; t <- 1:120 ; m <- length(t)
beta <- c(1.5,0.6,0.6)
omega <- 2*pi/52
#log mu_{0,t}
alpha <- 0.2
base <- beta[1] + beta[2] * cos(omega*t) + beta[3] * sin(omega*t) 
#Generate example data with changepoint and tau=tau
tau <- 100
kappa <- 0.4
mu0 <- exp(base)
mu1 <- exp(base  + kappa) 

#Generate data (set.seed has problems with RUnit)
#set.seed(42)
#x <- rnbinom(length(t),mu=mu0*(exp(kappa)^(t>=tau)),size=1/alpha)
x <- c(11,5,10,16,19,25,7,11,11,2,13,8,13,12,6,11,10,4,2,0,3,3,1,0,2,2,1,2,0,3,0,1,2,1,3,4,0,1,1,7,1,5,1,2,4,1,2,7,8,9,13,8,1,10,11,14,9,6,10,15,5,1,9,3,26,3,1,8,7,2,5,8,2,2,7,0,3,0,1,0,1,0,1,1,2,2,5,5,3,1,0,2,5,1,4,2,9,4,5,2,5,10,21,6,13,14,19,22,17,7,11,13,17,22,12,13,12,10,4,3)

s.ts <- create.disProg(week=1:length(t),observed=x,state=(t>=tau))

#Define two control objects (should be equivalent)
cntrl1 = list(range=t,c.ARL=5, mu0=mu0, alpha=alpha, change="intercept",ret="value",dir="inc")
cntrl2 = list(range=t,c.ARL=5, M=-1, mu0=mu0, alpha=alpha, change="intercept",ret="value",dir="inc")

#Run
glr.ts1 <- algo.glrnb(s.ts,control=c(cntrl1))
glr.ts2 <- algo.glrnb(s.ts,control=c(cntrl2))

correctUpperbound <- c(0.0933664,0,0.001387989,0.4392282,1.239898,
2.983766,1.954988,1.722341,1.586777,0.7331938,
0.9337575,0.7903225,1.104522,1.425098,1.241290,
1.633672,2.033343,1.788079,1.397671,0.9081794,
0.797097,0.7270934,0.5248943,0.3093548,0.2622768,
0.2301054,0.1595651,0.1484989,0.06889605,0.1504776,
0.04138495,0.02219845,0.0231524,0.00957569,0.1504776,
0.5827537,0.0357062,0.005011513,0,1.390972,
0.3167743,0.5717088,0.1053871,0.003442552,0.0005934715,
0,0,0.05509335,0.1375619,0.2449853,
0.6840703,0.5427538,0.05675776,0.06656547,0.09036596,
0.209314,0.1392091,0.03494786,0.02621600,0.277202,
0.01762547,0,0,0,3.564077,
1.410190,0.290548,0.3740241,0.4269062,0.1296794,
0.1298662,0.6322042,0.2115204,0.1074570,0.9366399,
0.1379007,0.1509654,0.03392803,0.005775552,0,
0,0,0,0,0.001143512,
0.001637927,1.021689,1.965804,1.830440,1.017412,
0.3033473,0.1689957,0.4051742,0.1247774,0.1460143,
0.03590031,0.9459381,0.4189531,0.2637725,0.03925406,
0.01374443,0.2283519,2.535301,1.406133,1.692899,
2.021258,2.951635,4.25683,4.77543,3.90064,
3.646361,3.680106,4.236502,5.522696,0.1221651,
0.4054735,0.6761779,0.803913,0.3913383,0.1261521)

 
### --- Test functions ---
 
test.simple <- function()
{
  upperboundDiff <- as.vector(glr.ts1$upperbound) - correctUpperbound
  checkTrue(sum(abs(upperboundDiff) > 1e-5) == 0)

  #M and Mtilde arguments
  checkTrue(digest(glr.ts1$upperbound) == digest(glr.ts2$upperbound))

}
