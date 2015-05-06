#' Beta-adjust Function
#'
#' Given a set of minimum P-values under the different numbers of tests, 
#' returns p-values adjusted using the Beta correction that
#' allows you to control the family wise error rates.
#' @param p a vector of the minimum P-values.
#' @param ntests the numbers of tests for the set of P-values.
#' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
#' @param fixParam an indicator vector (length = 3) of 0 (estimate) or 1 (fix) for alpha, beta and gamma parameters.
#' @param trace printing convergence status
#' @examples
#' p=NULL
#' for(i in 1:100){p = c(p, min(runif(i)))}
#' b.adjust(p, 1:100)

b.adjust <-
function(p, ntests, sub=rep(T,length(p)), fixParam=c(0,0,0), trace=F){
	if((sum(p<0)+sum(p>1))>0){stop("Incorrect P-values!")}
	rmflag = is.na(p) | ntests==0 | p==1 | p==0 | is.na(ntests)
	rmflag[is.na(rmflag)]=T
	rmflag[!sub]=T
	#p=p[!rmflag]
	#ntests=ntests[!rmflag]
	badj=.C("nlmForR", as.double(p[!rmflag]), as.integer(ntests[!rmflag]), as.integer(sum(!rmflag)), estimate=as.double(c(1,1,1)), lkhd=double(1), gradient=double(3), hessian=double(6), double(300), as.integer(fixParam), as.integer(trace))[4:6]
	res = badp(p, ntests, badj)
	attr(res, "b.adjust")=badj
	res
}
