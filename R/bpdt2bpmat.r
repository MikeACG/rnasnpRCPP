bpdt2bpmat <- function(bpdt, n) {

    bpmat <- matrix(0, nrow = n, ncol = n)
    I <- rbind(cbind(bpdt[[1]], bpdt[[2]]), cbind(bpdt[[2]], bpdt[[1]]))
    bpmat[I] <- rep(bpdt[[3]], 2)

    return(bpmat)

}

