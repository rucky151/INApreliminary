    


#' Evaluate nodes for value in invasion detection (based on a specified number of simulations)
#'
#' For an analysis of each potential starting node in turn, outputs a table of number of nodes infected by time detected at each potential sampling node, with rows being starting node and columns being sampling node - for a specified number of simulation
#' @param adjmat adjacency matrix for evaluation
#' @param stoch logical var indicating whether adjacency matrix entries are fixed or probabilities
#' @param nsim number of simulations to be analyzed
#' @keywords prioritization sampling
#' @export
#' @examples
#' multisimtest(adjmat=Amat, wtvec=__, nsim=1)  # multiple sims not needed because not stochastic - BUT, yes, now should be stochastic
#' multisimtest(adjmat=sAmat, wtvec=__, stoch=T, nsim=10)


# to do - GT
# to do - improve examples
# to do - consider skewedness?  kurtosis?
# to do - consider removing option without stochasticity


multisimtest <- function(adjmat, wtvec, stoch=T, nsim=1) {

  dimL <- dim(adjmat)[1]

  outarr <- array(-99, c(dimL,dimL,nsim)) # all results
  meanarr <- matrix(-99, ncol=dimL, nrow=dimL) # mean of results
  vararr <- meanarr # variance of results

  for (i3 in 1:nsim){
    tempout <- multistart(adjmat=adjmat, wtvec=wtvec, stoch=stoch)
    outarr[,,i3] <- tempout
  }

  for (i3 in 1:dimL) {
    for (i4 in 1:dimL) {
      meanarr[i3,i4] <- mean(outarr[i3,i4,])
      vararr[i3,i4] <- var(outarr[i3,i4,])
  } }
# consider potential use of apply to replace the second loop

  list(outarr=outarr, meanarr=meanarr, vararr=vararr)

}


