voomWithQualityWeightsMOD <- function (counts, design = NULL, lib.size = NULL, normalize.method = "none",
                                       plot = FALSE, span = 0.5, var.design = NULL, method = "genebygene", 
                                       maxiter = 50, tol = 1e-10, trace = FALSE, replace.weights = TRUE, 
                                       col = NULL, ...) {
  if (plot) {
    oldpar <- par(mfrow = c(1, 2))
    on.exit(par(oldpar))
  }
  v <- voomMod(counts, design = design, lib.size = lib.size, normalize.method = normalize.method, 
               plot = FALSE, span = span, ...)
  aw <- arrayWeights(v, design = design, method = method, maxiter = maxiter, 
                     tol = tol, var.design = var.design)
  v <- voomMod(counts, design = design, weights = aw, lib.size = lib.size, 
               normalize.method = normalize.method, plot = plot, span = span, 
               ...)
  aw <- arrayWeights(v, design = design, method = method, maxiter = maxiter, 
                     tol = tol, trace = trace, var.design = var.design)
  wts <- asMatrixWeights(aw, dim(v)) * v$weights
  attr(wts, "arrayweights") <- NULL
  if (plot) {
    barplot(aw, names = 1:length(aw), main = "Sample-specific weights", 
            ylab = "Weight", xlab = "Sample", col = col)
    abline(h = 1, col = 2, lty = 2)
  }
  if (replace.weights) {
    v$weights <- wts
    v$sample.weights <- aw
    return(v)
  }
  else {
    return(wts)
  }
}

# Author Julien Roux

voomMod <- function(counts,design=NULL,lib.size=NULL,normalize.method="none",plot=FALSE,span=0.5,...) {
  ## Linear modelling of count data mean-variance modelling at the observational level.
  ## Directly modified from voom code, to use the function cpm() of edgeR, which scales the prior counts used to calculate the log2(cpm)
  out <- list()
  ##	Check counts
  if(is(counts,"DGEList")) {
    out$genes <- counts$genes
    out$targets <- counts$samples
    if(is.null(design) && diff(range(as.numeric(counts$sample$group)))>0) design <- model.matrix(~group,data=counts$samples)
    if(is.null(lib.size)) lib.size <- with(counts$samples,lib.size*norm.factors)
    counts <- counts$counts
  } else {
    isExpressionSet <- suppressPackageStartupMessages(is(counts,"ExpressionSet"))
    if(isExpressionSet) {
      if(length(fData(counts))) out$genes <- fData(counts)
      if(length(pData(counts))) out$targets <- pData(counts)
      counts <- exprs(counts)
    } else {
      counts <- as.matrix(counts)
    }
  }
  ##	Check design
  if(is.null(design)) {
    design <- matrix(1,ncol(counts),1)
    rownames(design) <- colnames(counts)
    colnames(design) <- "GrandMean"
  }
  ##	Check lib.size
  if(is.null(lib.size)) lib.size <- colSums(counts)
  ##	Fit linear model to log2-counts-per-million
  ## y <- t(log2(t(counts+0.5)/(lib.size+1)*1e6)) ## Julien: this is replaced by line below
  y <- cpm(counts, lib.size=lib.size, log=TRUE, prior.count=0.25)
  y <- normalizeBetweenArrays(y,method=normalize.method)
  fit <- lmFit(y,design,...)
  if(is.null(fit$Amean)) fit$Amean <- rowMeans(y,na.rm=TRUE)

  ##	Fit lowess trend to sqrt-standard-deviations by log-count-size
  sx <- fit$Amean+mean(log2(lib.size+1))-log2(1e6)
  sy <- sqrt(fit$sigma)
  allzero <- rowSums(counts)==0
  if(any(allzero)) {
    sx <- sx[!allzero]
    sy <- sy[!allzero]
  }
  l <- lowess(sx,sy,f=span)
  if(plot) {
    plot(sx,sy,xlab="log2( count size + scaled prior )",ylab="Sqrt( standard deviation )",pch=16,cex=0.25)
    title("voom: Mean-variance trend")
    lines(l,col="red")
  }
  ##	Make interpolating rule
  f <- approxfun(l, rule=2)
  ##	Find individual quarter-root fitted counts
  if(fit$rank < ncol(design)) {
    j <- fit$pivot[1:fit$rank]
    fitted.values <- fit$coef[,j,drop=FALSE] %*% t(fit$design[,j,drop=FALSE])
  } else {
    fitted.values <- fit$coef %*% t(fit$design)
  }
  fitted.cpm <- 2^fitted.values
  fitted.count <- 1e-6 * t(t(fitted.cpm)*(lib.size+1))
  fitted.logcount <- log2(fitted.count)
  ##	Apply trend to individual observations
  w <- 1/f(fitted.logcount)^4
  dim(w) <- dim(fitted.logcount)
  ##	Output
  out$E <- y
  out$weights <- w
  out$design <- design
  if(is.null(out$targets))
    out$targets <- data.frame(lib.size=lib.size)
  else
    out$targets$lib.size <- lib.size
  new("EList",out)
}