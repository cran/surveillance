######################################################################
# Compute a suitable layout for plotting
######################################################################

magic.dim <- function(k){
  if(k==1)
    return(c(1,1))

  #factorize k
  factors <- primeFactors(k)

  #find the best factorization of k into two factors
  res <- bestCombination(factors)

  #if k is a prime or the difference between the two factors of k is too large
  #rather use the roots of the next square number greater than k

  #up is root of the smallest square number >= k
  up <- ceiling(sqrt(k))
  #low is root of the biggest square number < k
  low <- up -1

  if(diff(res) >5){
    # e.g. k=11 is a prime, the next square number is 16 so up=4 and low=3
    # low^2 = 9 < 11 is naturally too small, up^2=16 > 11 so c(4,4) is a solution
    # but low*up = 3*4 = 12 > 11 is also adequate and a better solution
    if((k - low^2) < up)
      res <- c(low,up)
    else
      res <- c(up,up)
  }

  return(sort(res))
}


######################################################################
# Compute the prime number factorization of an integer
######################################################################

primeFactors <- function(x){
  if(x==1)
    return(1)

  factors<- numeric(0)
  i<-1

  #start with i=2 and divide x by i (as often as possible) then try division by i+1
  #until all factors are found, i.e. x=1
  while(i < x){
    i <- i+1

    while((x %% i)==0){
      # each time a new factor i is found, save it and proceed with x = x/i
      # e.g. k=20: 2 is a factor of x=20, continue with x = 10 = 20/2
      #            2 is a factor of x=10, continue with x = 5 = 10/2
      #            3 and 4 are no factors of x = 5
      #            5 is a factor of x = 5, continue with x = 1
      # result: 20 = c(2, 2, 5)
      factors <- c(factors, i)
      x <- x/i
    }
  }
  return(factors)
}


######################################################################
# Given a prime number factorization of a number, e.g. 36
# yields x=c(2,2,3,3)
# and parition x into two groups, such that the product of the numbers
# in group one is as similar as possible to the product
# of the numbers of group two. This is useful in magic.dim
#
# Params:
#  x - the prime number factorization
#
# Returns:
#  c(prod(set1),prod(set2))
######################################################################

bestCombination <- function(x) {
  #Compute the power set of 0:1^length(x), i.e. a binary indicator for
  #variable stating whether to include it in set 1 or not.
  combos <- as.matrix(expand.grid(rep(list(0:1),length(x))))
  mode(combos) <- "logical"

  #Small helper function, given a vector of length(x) stating whether
  #to include an element in set1 or not, compute the product
  #of set1 and set2=x\backslash set1
  #set1: all those for which include is TRUE, set2: all those for which
  #include is FALSE
  setsize <- function(include) { c(prod(x[include]),prod(x[!include])) }

  #Compute the product of set1 and set2 for each possible combination
  sizes <- apply(combos,MARGIN=1,FUN=setsize)
  #Calculate the combination, where x is as close to y as possible
  bestConfig <- combos[which.min(abs(diff(sizes))),]
  #Return this setsize of this configuration
  return(setsize(bestConfig))
}
