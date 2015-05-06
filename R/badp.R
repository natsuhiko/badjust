#' Beta-adjusted P-values
#'
#' Given the parameter estimate by b.adjust, returns p-values adjusted.
#' @param p a vector of the minimum P-values
#' @param ntests the numbers of tests for the set of P-values
#' @badj the parameter estimate by b.adjust function
#' @examples
#' # not run
#' badp(p, 1:100, badj)


badp <-
function(p, ntests, badj){
	-expm1( pbeta(p, badj[[1]][1], badj[[1]][2], lower=F, log=T) * ntests^badj[[1]][3] )
}
