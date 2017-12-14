#' PeerRank score construction per Walsh (2014).
#'
#' @param A a square matrix of grades having been normalized in interval [0,1].
#' @param alpha an optional parameter. Must be in interval (0,1). Default = .10.
#' @param beta an optional parameter.  Must be in interval (0, 1], as must the
#'        sum of alpha and beta.  The Generalized PeerRank rule degenerates to
#'        the (non-generalized) PeerRank rule when beta is set to zero. Default
#'        = .10.
#' @references Walsh, T. 2014. The PeerRank method for peer assessment. In
#'             Proceedings of the 21st European Conference on Artificial
#'             Intelligence (ECAI), 909-914.
#'             doi: 10.3233/978-1-61499-419-0-909
#' @export

pr <- function(A, alpha=.10, beta=.10, precision=1e-10) {
    if (!is.matrix(A) || nrow(A) != ncol(A) || TRUE %in% is.na(A)
        || min(A) < 0 || max(A) > 1 ) {
        stop("Grades must be given as a square matrix and in interval [0,1].")
    }
    else if (alpha <= 0 || beta < 0 || alpha >= 1 || beta >= 1
             || alpha + beta > 1) {
        stop("Alpha and beta must each be in interval (0,1) with sum <= 1.")
    }
    m <- nrow(A)
    X <- list(rowSums(A)/m)
    while (TRUE) {
        X.last <- X[[length(X)]]
        term1 <- (1 - alpha - beta) * X.last
        term2 <- (alpha/sum(X.last)) * apply(A, 1, function(z) sum(X.last * z))
        term3 <- (beta/m) * apply(A, 2, function(z) sum(1 - abs(z - X.last))) 
        X.new <- term1 + term2 + term3
        X <- c(X,list(X.new))
        if (max(abs(X.new-X.last)) <= precision ) { break; }
    }
    ranks <- list(estimates = X[[length(X)]], iterations = length(X), 
                  alpha = alpha, beta = beta, precision = precision, data = A,
                  agents = if(!is.null(rownames(A))) rownames(A) else c(1:m))
    class(ranks) <- "peerRank"
    return(ranks)
}


#' @export
print.peerRank <- function(X, ...){
    names(X[["estimates"]]) <- X[["agents"]]
    print(X[["estimates"]])
    
}
