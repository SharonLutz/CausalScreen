CSErrorCheck <-
function(X,P1,P2,l,k,y,z=NULL){

  # check X is a matrix or vector
  if(!is.null(X)){

    # Check k is a vector the same length as nrow(X) or length(X)
    if(!(length(k)==nrow(X))){stop("Error: length(k) must equal nrow(X).")}
    
    # Check l is a vector matching k, y is a vector matching l and k
    if(!(length(k)==length(l) & length(l)==length(y))){stop("Error: y, l and k must be vectors of the same length.")}
    
    # Check that neither P1 nor P2 is null
    if(is.null(P1) | is.null(P2)){stop("Error: neither p1 or p2 can be NULL.")}
    # P1 <- as.matrix(P1)
    # P2 <- as.matrix(P2)
    if(!(nrow(P1)==nrow(X) & nrow(P2)==nrow(X) & ncol(P1)==ncol(X) & ncol(P2)==ncol(X))){stop("Error: P1 and P2 must have the same number of dimensions as X.")}
    
    # Check if X has no column names
    if(is.null(colnames(X))){stop("Error: X has no column names.")}
    
    # Check dimensions of z, must be same length as X
    if(!is.null(z)){
      if(is.vector(z)){
        if(length(z)!=nrow(X)){stop("Error: length(z) must equal nrow(X).")}
      }else if(!is.null(nrow(z)) & !is.null(ncol(z))){
        if(nrow(z)!=nrow(X)){stop("Error: nrow(z) must equal nrow(X).")}
      }else{stop("Error: z must either be a vector or matrix of covariates.")}
    }
    
  }else{
    stop("Error: X cannot be NULL")
  }
}
