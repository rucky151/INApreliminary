
#' Evaluate individual nodes for their value in invasion detection (invasion starting from a single node) in one simulation
#'
#' For a given starting node (only one at this point), yields a vector of the number of nodes invaded for each sampling node by the time invasion reaches (is detected at) that sampling node, and the time step of first invasion at the sampling node - for one simulation.  Includes the option for differentially weighted probabilities of the existence of links.
#' @param adjmat adjacency matrix for evaluation
#' @param wtvec vector of weights corresponding to nodes associated with the adjacency matrix adjmat
#' @param start.choice number of node where invasion starts
#' @param stoch logical var indicating whether adjacency matrix entries are fixed or probabilities
#' @keywords prioritization sampling
#' @export
#' @examples
#' onestart(adjmat=Amat, wtvec= ___, start.choice=2) # assumes Amat exists - more complete example to be added
#' onestart(adjmat=sAmat, wtvec= ___, start.choice=2,stoch=T) # assumes Amat exists - more complete example to be added

# to do - GT
# to do - improve examples
# to do - adjust option for deterministic 
# to do - consider adjusting assumption that the matrix would start out as ones and zeroes, and then weights modify

onestart <- function(adjmat, wtvec, start.choice, stoch=T){
   
  dimL <- dim(adjmat)[1] # number of rows of adjacency matrix
  t1 <- matrix(0 * 1:dimL, nrow=1)
  t1[,start.choice] <- 1 # starting node given status 1
  
  # if stochastic , generate the adjacency matrix for this simulation
    # if stochastic, assumes the adjacency matrix is made up of probabilities [might modify]

  if(stoch){ # consider temporal resolution of stochasticity - this version sets to one mat indefinitely within sim
    adjmat <- t(wtvec * t(adjmat)) # weighting
    adjmat <- (matrix(runif(dimL^2), ncol=dimL) < adjmat)
  }

  # generate the infection status for each node at each time point until spread stops

  outmat <- t1
  infcount.pre <- 1
  infcount.post <- -9 # perhaps add improvement 

  while(sum(t1) < dimL & infcount.pre != infcount.post){
    infcount.pre <- sum(t1) # perhaps change to infcount.pre <- infcount.post
    t1 <- as.numeric(t1 %*% adjmat > 0)
    infcount.post <- sum(t1)
    outmat <- rbind(outmat,t1) # row i is the ith time point
  }

  # find the num infected for each single sampling node, 
     # assuming infected node stays infected
	# consider whether to set diag

  # find the number of nodes infected at each time
  inft <- rowSums(outmat)

  # find the first time each node is invaded, and number of nodes infected at that time
  firstt <- -99 + 0*1:dimL # first time each node is invaded
  numt <- firstt # number of nodes infected at that time

  for (i in 1:dimL){ # col i is the ith sample node

    if (sum(outmat[,i] > 0)){ # if node i ever gets infected
      firstt[i] <- min(which(outmat[,i] > 0)) # first time
      numt[i] <- inft[firstt[i]] # how many nodes at that time
    }
    else {
      firstt[i] <- Inf
      numt[i] <- max(inft)
    }

  }

  # for each node, time first invaded and number of nodes infected at that time
  sampnode <- cbind(firstt, numt)

  list(outmat=outmat, sampnode=sampnode)
}

