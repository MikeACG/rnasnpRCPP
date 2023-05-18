bpdt2bpmat <- function(bpdt, n) {

    bpmat <- matrix(0, nrow = n, ncol = n)
    I <- cbind(bpdt[[1]], bpdt[[2]])
    bpmat[I] <- bpdt[[3]]

    return(bpmat)

}

