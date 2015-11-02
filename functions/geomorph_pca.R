geomorph_pca <- function (A, axis1 = 1, axis2 = 2, warpgrids = TRUE, mesh = NULL, 
          label = NULL, groups = NULL, verbose = FALSE) 
{
  if (length(dim(A)) != 3) {
    stop("Data matrix not a 3D array (see 'arrayspecs').")
  }
  if (any(is.na(A)) == T) {
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').")
  }
  k <- dim(A)[2]
  p <- dim(A)[1]
  n <- dim(A)[3]
  ref <- mshape(A)
  x <- two.d.array(A)
  pc.res <- prcomp(x)
  pcdata <- pc.res$x
  
  pcaxis.min.1 <- min(pcdata[, axis1])
  pcaxis.max.1 <- max(pcdata[, axis1])
  pc.min.1 <- pc.max.1 <- rep(0, dim(pcdata)[2])
  pc.min.1[axis1] <- pcaxis.min.1
  pc.max.1[axis1] <- pcaxis.max.1
  shape.min.1 <- arrayspecs(as.matrix(pc.min.1 %*% (t(pc.res$rotation))), 
                            p, k)[, , 1] + ref
  shape.max.1 <- arrayspecs(as.matrix(pc.max.1 %*% (t(pc.res$rotation))), 
                            p, k)[, , 1] + ref
  pcaxis.min.2 <- min(pcdata[, axis2])
  pcaxis.max.2 <- max(pcdata[, axis2])
  pc.min.2 <- pc.max.2 <- rep(0, dim(pcdata)[2])
  pc.min.2[axis2] <- pcaxis.min.2
  pc.max.2[axis2] <- pcaxis.max.2
  shape.min.2 <- arrayspecs(as.matrix(pc.min.2 %*% (t(pc.res$rotation))), 
                            p, k)[, , 1] + ref
  shape.max.2 <- arrayspecs(as.matrix(pc.max.2 %*% (t(pc.res$rotation))), 
                            p, k)[, , 1] + ref
  shapes <- list(shape.min.1, shape.max.1, shape.min.2, shape.max.2)
  names(shapes) <- c(paste("PC", axis1, "min", sep = ""), 
                     paste("PC", axis1, "max", sep = ""), paste("PC", axis2, 
                                                                "min", sep = ""), paste("PC", axis2, "max", sep = ""))
  if (verbose == TRUE) {
    return(list(pc.summary = summary(pc.res), pc.scores = pcdata, 
                pc.shapes = shapes))
  }
  if (verbose == FALSE) {
    return(pc.summary = summary(pc.res))
  }
}
