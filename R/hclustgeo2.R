hclustgeo2<-function (D0, D1 = NULL, alpha = 0, scale = TRUE, wt = NULL, method="ward.D") 
{
res<-NULL    
if (class(D0) != "dist") 
        stop("DO must be of class dist (use as.dist)", 
            call. = FALSE)
    if (!is.null(D1) && (class(D1) != "dist")) 
        stop("D1 must be of class dist (use as.dist)", 
            call. = FALSE)
    n <- as.integer(attr(D0, "Size"))
    if (is.null(n)) 
        stop("invalid dissimilarities", call. = FALSE)
    if (is.na(n) || n > 65536L) 
        stop("size cannot be NA nor exceed 65536", call. = FALSE)
    if (n < 2) 
        stop("must have n >= 2 objects to cluster", call. = FALSE)
    if (!is.null(D1) && length(D0) != length(D1)) 
        stop("the two dissimilarity structures must have the same size", 
            call. = FALSE)
    if ((max(alpha) > 1) || (max(alpha) < 0)) 
        stop("Values alpha must be in [0,1]", call. = FALSE)
    if ((scale == TRUE) && (!is.null(D1))) {
        D0 <- D0/max(D0)
        D1 <- D1/max(D1)
    }
    delta0 <- wardinit(D0, wt)
    if (!is.null(D1)) 
        delta1 <- wardinit(D1, wt)
    else delta1 <- 0
    delta <- (1 - alpha) * delta0 + alpha * delta1
    res$solution <- hclust(delta, method = method, members = wt)
	res$weighteddiss<-delta
    return(res)
}




