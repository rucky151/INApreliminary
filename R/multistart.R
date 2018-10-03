


#' Evaluate nodes for value in invasion detection (invasion starting from each node in turn)
#'
#' For an analysis of each potential starting node in turn, generates a table of number of nodes infected by time detected at each potential sampling node, with rows being starting node and columns being sampling node - for one simulation.  Each node can be weighted with a relative likelihood of spread to that node.
#' @param adjmat adjacency matrix for evaluation
#' @param wtvec vector of weights corresponding to nodes associated with adjacency matrix adjmat2
#' @param stoch logical var indicating whether adjacency matrix entries are fixed or probabilities
#' @keywords prioritization sampling
#' @export
#' @examples
#' allstartf.wt(adjmat=Amat,wtvec= ___)
#' allstartf.wt(adjmat=sAmat,wtvec=____, stoch=T)

# to do - GT, especially for weighted components
# to do - consider only stochastic 
# to do - improve examples

multistart <- function(adjmat, wtvec, stoch=T){

  dimL <- dim(adjmat)[1]

  alloutmat <- matrix(-99, ncol=dimL, nrow=dimL)

  for (i in 1:dimL) {
    alloutmat[i,] <- onestart(adjmat=adjmat, wtvec=wtvec, start.choice=i, stoch=stoch)$sampnode[,2]
  }

alloutmat
}


