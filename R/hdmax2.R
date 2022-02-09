##' Epigenome Wide Association Study with both exposure and outcome
##'
##' This function uses lfmm (latent factor mixed models) to estimate
##' the effects of exposures and outcomes on a response matrix.
##'
##'
##' @param M a response variable matrix with n rows and p columns.
##' Each column corresponds to a beta-normalized methylation profile.
##' Response variables must be encoded as numeric. No NAs allowed.
##' @param X an explanatory variable matrix with n rows and d columns.
##' Each column corresponds to a distinct explanatory variable (Exposure).
##' Explanatory variables must be encoded as numeric variables.
##' @param Y an explanatory variable matrix with n rows and d columns.
##' Each column corresponds to a distinct explanatory variable (Outcome).
##' Explanatory variables must be encoded as numeric variables.
##' @param K an integer for the number of latent factors in the regression model.
##' @param conf set of covariable, must be numeric.
##' @return an object with the following attributes:
##'
##'  - U the latent variable score matrix with dimensions n x K.
##'
##'  - B the effect size matrix for the exposure X and the outcome Y.
##'
##'  - score matrix for the exposure X and the outcome Y.
##'
##'  - pValue matrix for the exposure X and the outcome Y.
##'
##'  - calibrated.score2, the calibrated score matrix for the exposure X and the outcome Y.
##'
##'  - calibrated.pvalue, the calibrated pValue matrix for the exposure X and the outcome Y.
##'
##'  - GIF : Genomic Inflation Factor for exposure and outcome
##'
##'  - lfmm : the result of the 2 regressions of lfmm, mod1 for the regression of X on M and mod2 for the regression of Y on M given X.
##'
##' @details
##' The response variable matrix Y and the explanatory variable are centered.
##' Missing values must be imputed. The number of latent factors can be estimated
##' by looking at the screeplot of eigenvalues of a PCA.
##' Possibility of calibrating the scores and pValues by the GIF (Genomic Inflation Factor).
##' See lfmm package for more information.
##' @export
##' @author Basile Jumentier
##' @examples
##'
##' library(hdmax2)
##'
##' # Run mEWAS
##'
##' res <- mEWAS(X = example$X, Y = example$Y, M = example$M, K = 5)
##'
##'
mEWAS <- function(X, Y, M, K, conf = NULL) {

  res <- list()

  # First regression
  dat <- lfmm::lfmm_ridge(Y = M, X = cbind(X, conf), K = K)
  res[[1]] <- dat
  dat <- lfmm::lfmm_test(Y = M, X = cbind(X, conf), lfmm = dat)

  pv1 <- dat$pvalue[, 1]
  sc1 <- dat$score[, 1]

  sc1.cal <- dat$calibrated.score2[, 1]
  pv1.cal <- dat$calibrated.pvalue[, 1]

  gif1 <- dat$gif[1]

  # Second regression
  dat <- lfmm::lfmm_ridge(Y = M, X = cbind(X, Y, conf), K = K)
  res[[2]] <- dat
  # ajout
  U <- dat$U
  dat <- lfmm::lfmm_test(Y = M, X = cbind(X, Y, conf), lfmm = dat)

  pv2 <- dat$pvalue[, 2]
  sc2 <- dat$score[, 2]

  sc2.cal <- dat$calibrated.score2[, 2]
  pv2.cal <- dat$calibrated.pvalue[, 2]

  gif2 <- dat$gif[2]

  names(res) <- c("mod1", "mod2")

  return(list(score = cbind(sc1, sc2),
              pValue = cbind(pv1, pv2),
              calibrated.score2 = cbind(sc1.cal, sc2.cal),
              calibrated.pvalue = cbind(pv1.cal, pv2.cal),
              gif = c(gif1, gif2),
              U = U,
              lfmm = res))
}



##' Compute the squared maximum of two series of pValues
##'
##' This function compute the squared maximum of two series of pValues from the mEWAS() function.
##' The objective of this function is to test all the markers and to determine which could be
##' potential mediators in the exposure-outcome association.
##'
##' @param pval1 vector of pValues (p*1) of exposure.
##' @param pval2 vector of pValues (p*1) of ouctome.
##' @param diagnostic.plot if TRUE the histogram of the p-values together
##' with the estimate of eta0 null line is plotted.
##' This is useful to visually check the fit of the estimated proportion of null p-values.
##' @param ... argument of the fdrtool function from the fdrtool package
##'
##' @return an object with the following attributes:
##'
##'  - a pValue for each markers
##'
##'  - a qValue for each markers
##'
##'  - the eta0 of the set of pValues
##'
##' @details
##' The pValue is computed for each markers following this formula
##'
##' \deqn{pV = max(pVal1, pVal2)^2}
##'
##' This quantity eta0, i.e. the proportion eta0 of null p-values in a given vector of p-values,
##' is an important parameter when controlling the false discovery rate (FDR).
##' A conservative choice is eta0 = 1 but a choice closer to the true value will
##' increase efficiency and power - see Benjamini and Hochberg (1995, 2000) and Storey (2002) for details.
##' We use the fdrtool package to transform pValues into qValues,
##' which allows us to control the FDR.
##' @export
##' @author Basile Jumentier
##' @examples
##'
##' library(hdmax2)
##'
##' # Run mEWAS
##'
##' res <- mEWAS(X = example$X, Y = example$Y, M = example$M, K = 5)
##'
##' # Run max2
##'
##' res <- max2(pval1 = res$pValue[, 1], pval2 = res$pValue[, 2])
##'
##' # Manhattan plot
##'
##' plot(-log10(res$pval), main = paste0("Eta0 = ", round(res$eta0, 3)))
##' abline(h = -log10(0.05 / ncol(example$M)))
##'
max2 <- function(pval1, pval2, diagnostic.plot = F, ...) {

  pval <- apply(cbind(pval1, pval2), 1, max)^2
  eta0 <- fdrtool::pval.estimate.eta0(pval, diagnostic.plot = diagnostic.plot)
  qval <- fdrtool::fdrtool(pval,statistic = "pvalue", plot = F, verbose = F, ...)$qval

  return(list(pval = pval,
              eta0 = eta0,
              qval = qval))
}

##' Run mediation analysis for a set of markers
##'
##' Estimate various quantities for causal mediation analysis for each
##' significant markers, including average causal mediation effects
##' (indirect effect), average direct effects, proportions mediated,
##' and total effect.
##'
##' @param qval set of qValues from max2() function
##' @param M a response variable matrix with n rows and p columns.
##' Each column corresponds to a beta-normalized methylation profile.
##' Response variables must be encoded as numeric. No NAs allowed.
##' @param X an explanatory variable matrix with n rows and d columns.
##' Each column corresponds to a distinct explanatory variable (Exposure).
##' Explanatory variables must be encoded as numeric variables.
##' @param Y an explanatory variable matrix with n rows and d columns.
##' Each column corresponds to a distinct explanatory variable (Outcome).
##' Explanatory variables must be encoded as numeric variables.
##' @param U set of latent factors from mEWAS() function (need include covariable)
##' @param FDR FDR threshold to pass markers in mediation analysis
##' @param sims number of Monte Carlo draws for nonparametric bootstrap or quasi-Bayesian approximation.
##' 10000 is recommended.
##' @param ... argument of the mediate function from the mediation package
##'
##' @return
##' Tables of results of mediation analyzes for markers with a qValue below the FDR threshold.
##' Indirect effect (ACME - average causal mediation effect), ADE (average direct effect),
##' PM (proportion mediated) and TE (total effect). Composition of tables: estimated effect,
##' confidence interval and mediation pValue.
##' We also return, We also return the results of the linear regressions.
##' The xm table corresponds to the regressions of X on Mi and
##' the my table to the regressions of Y on Mi knowing X.
##' With Mi corresponding to the different CpGs tested.
##'
##' @details
##'
##' We use the mediate() function of the mediation package on the set of markers having a qValue lower
##' than the FDR threshold. This function makes it possible to estimate their indirect effects and to
##' test their significance.
##'
##' @export
##' @author Basile Jumentier
##' @examples
##'
##' library(hdmax2)
##'
##' data(example)
##'
##' # Run mEWAS
##'
##' res <- mEWAS(X = example$X, Y = example$Y, M = example$M, K = 5)
##'
##' # Keep latent factors for mediation
##'
##' U <- res$U
##'
##' # Run max2
##'
##' res <- max2(pval1 = res$pValue[, 1], pval2 = res$pValue[, 2])
##'
##' # Run mediation (only 3 simulations for estimate and test indirect effect)
##'
##' res <- wrap_mediation(qval = res$qval,
##'                             X = example$X,
##'                             Y = example$Y,
##'                             M = example$M,
##'                             U = U, sims = 3,
##'                             FDR = 0.5)
##'
##'
wrap_mediation <- function(qval, X, Y, M, U = NULL, FDR = 0.1, sims = 3, ...) {

  if (is.null(colnames(M))) {
    colnames(M) <- 1:ncol(M)
  }

  M <- M[, qval <= FDR]


  # from package mediation
  ACME <- matrix(ncol = 4, nrow = ncol(M))
  ADE <- matrix(ncol = 4, nrow = ncol(M))
  PM <- matrix(ncol = 4, nrow = ncol(M))
  TE <- matrix(ncol = 4, nrow = ncol(M))

  # from linear models
  xm <- matrix(ncol = 4, nrow = ncol(M))
  my <- matrix(ncol = 4, nrow = ncol(M))

  for (i in 1:ncol(M)) {

    dat.x <- data.frame(X = X, Mi = M[, i], covar = U)
    dat.y <- data.frame(X = X, Mi = M[, i], covar = U, Y = Y)

    mod1 <- stats::lm(Mi ~ X + ., data = dat.x)
    mod2 <- stats::lm(Y ~ X + Mi + ., data = dat.y)

    # for linear models
    xm[i, ] <- summary(mod1)$coeff[2, ] # effect of X
    my[i, ] <- summary(mod2)$coeff[3, ] # effect of M

    med <- mediation::mediate(mod1, mod2, sims = sims, treat = "X", mediator = "Mi", ...)

    ACME[i, ] <- c(med$d0, med$d0.ci[1], med$d0.ci[2], med$d0.p)
    ADE[i, ] <- c(med$z0, med$z0.ci[1], med$z0.ci[2], med$z0.p)
    PM[i, ] <- c(med$n0, med$n0.ci[1], med$n0.ci[2], med$n0.p)
    TE[i, ] <- c(med$tau.coef, med$tau0.ci[1], med$tau0.ci[2], med$tau.p)
  }

  ACME <- as.data.frame(ACME)
  ADE <- as.data.frame(ADE)
  PM <- as.data.frame(PM)
  TE <- as.data.frame(TE)
  xm <- as.data.frame(xm)
  my <- as.data.frame(my)

  colnames(ACME) <- c("est", "CI_2.5", "CI_97.5", "pval")
  colnames(ADE) <- c("est", "CI_2.5", "CI_97.5", "pval")
  colnames(PM) <- c("est", "CI_2.5", "CI_97.5", "pval")
  colnames(TE) <- c("est", "CI_2.5", "CI_97.5", "pval")
  colnames(xm) <- c("Estimate", "Std.Error", "t.Value", "pValue")
  colnames(my) <- c("Estimate", "Std.Error", "t.Value", "pValue")

  ACME$CpG <- colnames(M)
  ADE$CpG <- colnames(M)
  PM$CpG <- colnames(M)
  TE$CpG <- colnames(M)
  xm$CpG <- colnames(M)
  my$CpG <- colnames(M)

  return(list(ACME = ACME,
              ADE = ADE,
              PM = PM,
              TE = TE,
              xm = xm,
              my = my))

}


##' Run mediation analysis for a set of markers
##'
##' Function adapt from the combp function() of the ENmix package
##'
##' @param data A data frame from bed format file with colname name
##' "V1","V2", "V3","V4","V5",V1 indicate chromosome (1,2,3,...,X,Y),
##' V2 is chromosome position, V4 is for P value and V5 for name of CpGs.
##' @param dist.cutoff Maximum distance in base pair to combine adjacent DMRs.
##' @param bin.size bin size for autocorrelation calculation.
##' @param seed FDR significance threshold for initial selection of DMR region.
##' @param nCores Number of computer cores used in calculation
##'
##' @return
##' Results of the DMRs analysis.
##'
##' @details
##'
##' The input should be a data frame with column name V1-V5, indicating chromosome, start position,end position,
##' pValues and probe names. The function will use a modified comb-p method to identify
##' differentially methylated regions.
##'
##' @author Basile Jumentier
##'
combp2 <- function (data, dist.cutoff = 1000, bin.size = 310, seed = 0.01, nCores = 10) {

  ##### a function to get a table of p-values for estimating acf
  #####loc should be increasing;
  acf.table<-function(x,loc,dist.cutoff){
    flag=TRUE; lag=1; result=NULL
    while(flag){
      x1=utils::head(x,-lag); x2=utils::tail(x,-lag); dist=diff(loc,lag=lag)
      index=(dist<dist.cutoff)
      if(all(!index)){flag=FALSE}else{
        result=rbind(result,data.frame(x1=x1[index],x2=x2[index],dist=dist[index]))
        lag=lag+1
      }
    }
    return(result)
  }

  ##### a function to estimate acf
  get.acf<-function(data,dist.cutoff,bin.size){
    temp<-NULL
    for (chr in unique(data$V1)){
      y<-data[data$V1==chr,]; y<-y[order(y$V3),]
      temp<-rbind(temp,acf.table(y$V4,y$V3,dist.cutoff))
    }
    bin.label<-findInterval(temp$dist,seq(bin.size,dist.cutoff,bin.size))
    temp.stouffer<-by(temp,bin.label,FUN=function(x){stats::cor.test(stats::qnorm(x$x1),
                                                                     stats::qnorm(x$x2),alternative="greater")},simplify=FALSE)

    cor.stouffer<-sapply(temp.stouffer,function(x){x$estimate})
    p.stouffer<-sapply(temp.stouffer,function(x){x$p.value})

    if (any(p.stouffer>0.05)){
      index=min(which(p.stouffer>0.05))
      cor.stouffer[index:length(cor.stouffer)]=0
    }
    return(cor.stouffer)
  }

  if (nCores > parallel::detectCores()) {
    nCores = parallel::detectCores()
  }
  data = as.data.frame(data)
  acf <- get.acf(data, dist.cutoff, bin.size)
  result <- parallel::mclapply(unique(data$V1), function(chr) {
    y = data[data$V1 == chr, ]
    y = y[order(y$V3), ]
    pos = y$V3
    p = stats::qnorm(y$V4)
    temp = sapply(pos, function(i) {
      index.i = (abs(pos - i) < bin.size)
      if (sum(index.i) > 1) {
        int <- findInterval(c(stats::dist(pos[index.i])), c(bin.size,
                                                            2 * bin.size))
        sd <- sqrt(sum(acf[int + 1]) * 2 + sum(index.i))
        return(stats::pnorm(sum(p[index.i]), mean = 0, sd = sd))
      }
      else {
        return(y$V4[index.i])
      }
    })
    return(data.frame(chr, start = pos, end = pos, s.p = temp))
  }, mc.cores = nCores)
  result <- do.call("rbind", result)
  names(result) = c("chr", "start", "end", "s.p")
  result = result[stats::p.adjust(result$s.p, method = "fdr") < seed,]

  result.fdr = NULL
  if (nrow(result) > 0) {
    for (chr in unique(result$chr)) {
      y = data[data$V1 == chr, ]
      y = y[order(y$V3), ]
      pos = y$V3
      p = stats::qnorm(y$V4)
      result.chr = result[result$chr == chr, ]
      a = IRanges::IRanges(start = result.chr$start, end = result.chr$end)
      b = IRanges::reduce(a, min.gapwidth = dist.cutoff)
      start = IRanges::start(b)
      end = IRanges::end(b)
      region.max <- max(Biostrings::width(b))
      temp = sapply(1:length(b), function(i) {
        index.i = (pos >= start[i] & pos <= end[i])

        # print(sum(index.i))

        if (sum(index.i) > 1) {
          int <- findInterval(c(stats::dist(pos[index.i])),
                              seq(bin.size, region.max + bin.size, bin.size))
          sd <- sqrt(sum(ifelse(int < length(acf), acf[int +
                                                         1], 0)) * 2 + sum(index.i))
          return(stats::pnorm(sum(p[index.i]), mean = 0, sd = sd))
        }
        else {
          return(y$V4[index.i])
        }
      })
      result.fdr = rbind(result.fdr, data.frame(chr, start,
                                                end, p = temp))
    }
    result.fdr$fdr = stats::p.adjust(result.fdr$p, method = "fdr")
    result.fdr <- result.fdr[order(result.fdr$p), ]
    result.fdr$start = (result.fdr$start - 1)
  }

  return(result.fdr)
}


##' Find AMR
##'
##' To identify aggregated methylated regions (AMR) using a modified comb-p method.
##'
##' @param chr chromosomes
##' @param start chromosomal position of markers (start)
##' @param end chromosomal position of markers (end)
##' @param pval pValues for each markers, from the max2 function
##' @param cpg name of each markers
##' @param ... see help of combp of ENmix package
##'
##' @return
##' A set of selected AMRs.
##'
##' @details
##'
##' The function will use a modified comb-p method to identify
##' aggregated methylated regions (AMRs).
##'
##' @export
##' @author Basile Jumentier
##' @examples
##'
##' # see vignette
AMR_search <- function(chr, start, end, pval, cpg, ...) {

  tmp <- data.frame(chr, start, end, pval, cpg)
  colnames(tmp) <- paste0("V", 1:5)

  tmp <- combp2(tmp, ...)

  return(list(res = tmp,
              data = data.frame(chr, start, end, pval, cpg)))
}


##' Build AMR vector
##'
##' To build a vector for each AMR find with AMR_search
##'
##' @param res result object of AMR_search function
##' @param methylation chromosomal position of markers (start)
##' @param nb_cpg chromosomal position of markers (end)
##'
##' @return
##' A set of build AMRs.
##'
##' AMR_mean contains the first components of each PCA for each DMR.
##' CpG_for_each_AMR contains the list of markers (CpGs) present on each AMR.
##'
##'
##' @details
##'
##' We use the series of pValues (one pValue per CpGs) obtained with the mEWAS
##' regression method and the combination of pValue max2.
##' To determine the potential AMRs used the combp method present in the ENmix package (Xu et al. 2016).
##' This method uses the Fisher method to combine the pValues and also the base pair distance (bP)
##' between CpGs (1000 bP maximum between nb_cpg CpGs on the same AMR).
##' The information for each AMR is summarized by doing the mean (by row) of each CpG.
##'
##' @export
##' @author Basile Jumentier
##' @examples
##'
##' # see vignette
AMR_build <- function(res, methylation, nb_cpg = 2) {

  data <- res$data
  res <- res$res

  # Number of CpG per DMR

  nb <- NULL

  for (i in 1:nrow(res)) {

    chri <- as.character(res$chr[i])

    tmp <- dplyr::filter(data, chr == chri)

    nb <- c(nb, sum((res$start[i]:res$end[i]) %in% tmp$start))
  }

  # Select DMRs with nb_cpg CpGs at minimum

  res <- cbind(res, nb)

  res <- dplyr::filter(res, nb >= nb_cpg)

  DMR.select <- list()

  for (i in 1:nrow(res)) {

    chri <- as.character(res$chr[i])

    tmp <- dplyr::filter(data, chr == chri)

    # DMR.select[[i]] <- tmp$cpg[(tmp$start %in% (res$start[i]:res$end[i]))]
    # THE CHANGE
    DMR.select[[i]] <- as.character(tmp$cpg[(tmp$start %in% (res$start[i]:res$end[i]))])
  }

  # Select CpGs values in the methylation matrix

  DMR.meth <- list()

  for (i in 1:length(DMR.select)) {
    DMR.meth[[i]] <- methylation[, DMR.select[[i]]]
  }

  # Built a vector for each DMR with the first component of PCA or with the rowmeans

  DMR.acp <- as.data.frame(matrix(ncol = length(DMR.meth), nrow = nrow(methylation)))
  colnames(DMR.acp) <- paste0("DMR", 1:length(DMR.meth))

  DMR.mean <- as.data.frame(matrix(ncol = length(DMR.meth), nrow = nrow(methylation)))
  colnames(DMR.mean) <- paste0("DMR", 1:length(DMR.meth))

  for (i in 1:length(DMR.meth)) {
    DMR.acp[, i] <- prcomp(DMR.meth[[i]])$x[, 1]
    DMR.mean[, i] <- rowMeans(DMR.meth[[i]])
  }

  # data

  res <- cbind(DMR = colnames(DMR.acp), res)
  names(DMR.select) <- colnames(DMR.acp)

  return(list(AMR_acp = DMR.acp,
              AMR_mean = DMR.mean,
              res = res,
              CpG_for_each_AMR = DMR.select))
}

##' Run mediation analysis on AMR
##'
##' Estimate various quantities for causal mediation analysis for each
##' AMRs, including average causal mediation effects
##' (indirect effect), average direct effects, proportions mediated,
##' and total effect.
##'
##' @param AMR a matrix of DMRs from the AMR_build() function (AMR_mean).
##' @param X an explanatory variable matrix with n rows and d columns.
##' Each column corresponds to a distinct explanatory variable (Exposure).
##' Explanatory variables must be encoded as numeric variables.
##' @param Y an explanatory variable matrix with n rows and d columns.
##' Each column corresponds to a distinct explanatory variable (Outcome).
##' Explanatory variables must be encoded as numeric variables.
##' @param U set of latent factors from mEWAS() function (need include covariable)
##' @param sims number of Monte Carlo draws for nonparametric bootstrap or quasi-Bayesian approximation.
##' 10000 is recommended.
##'
##' @return
##' Tables of results of mediation analyzes for markers with a qValue below the FDR threshold.
##' Indirect effect (ACME - average causal mediation effect), ADE (average direct effect),
##' PM (proportion mediated) and TE (total effect). Composition of tables: estimated effect,
##' confidence interval and mediation pValue.
##' We also return, We also return the results of the linear regressions.
##' The xm table corresponds to the regressions of X on DMRi and
##' the my table to the regressions of Y on AMRi knowing X.
##' With DMRi corresponding to the different AMRs tested.
##'
##' @details
##'
##' We use the mediate() function of the mediation package on the set of selected AMRs.
##' This function makes it possible to estimate their indirect effects and to
##' test their significance.
##'
##' @export
##' @author Basile Jumentier
##' @examples
##'
##' # see vignette
wrap_mediation_AMR <- function(X, Y, AMR, U = NULL, sims = 3) {

  DMR <- AMR

  ACME <- matrix(ncol = 4, nrow = ncol(DMR))
  ADE <- matrix(ncol = 4, nrow = ncol(DMR))
  PM <- matrix(ncol = 4, nrow = ncol(DMR))
  TE <- matrix(ncol = 4, nrow = ncol(DMR))

  # from linear models
  xm <- matrix(ncol = 4, nrow = ncol(DMR))
  my <- matrix(ncol = 4, nrow = ncol(DMR))

  for (i in 1:ncol(DMR)) {

    dat.x <- data.frame(X = X, Mi = DMR[, i], covar = U)
    dat.y <- data.frame(X = X, Mi = DMR[, i], covar = U, Y = Y)

    mod1 <- stats::lm(Mi ~ X + ., data = dat.x)
    mod2 <- stats::lm(Y ~ X + Mi + ., data = dat.y)

    # for linear models
    xm[i, ] <- summary(mod1)$coeff[2, ] # effect of X
    my[i, ] <- summary(mod2)$coeff[3, ] # effect of M

    med <- mediation::mediate(mod1, mod2, sims = sims, treat = "X", mediator = "Mi")

    ACME[i, ] <- c(med$d0, med$d0.ci[1], med$d0.ci[2], med$d0.p)
    ADE[i, ] <- c(med$z0, med$z0.ci[1], med$z0.ci[2], med$z0.p)
    PM[i, ] <- c(med$n0, med$n0.ci[1], med$n0.ci[2], med$n0.p)
    TE[i, ] <- c(med$tau.coef, med$tau0.ci[1], med$tau0.ci[2], med$tau.p)
  }

  ACME <- as.data.frame(ACME)
  ADE <- as.data.frame(ADE)
  PM <- as.data.frame(PM)
  TE <- as.data.frame(TE)
  xm <- as.data.frame(xm)
  my <- as.data.frame(my)

  colnames(ACME) <- c("est", "CI_2.5", "CI_97.5", "pval")
  colnames(ADE) <- c("est", "CI_2.5", "CI_97.5", "pval")
  colnames(PM) <- c("est", "CI_2.5", "CI_97.5", "pval")
  colnames(TE) <- c("est", "CI_2.5", "CI_97.5", "pval")
  colnames(xm) <- c("Estimate", "Std.Error", "t.Value", "pValue")
  colnames(my) <- c("Estimate", "Std.Error", "t.Value", "pValue")

  ACME$DMR <- colnames(DMR)
  ADE$DMR <- colnames(DMR)
  PM$DMR <- colnames(DMR)
  TE$DMR <- colnames(DMR)
  xm$CpG <- colnames(DMR)
  my$CpG <- colnames(DMR)

  return(list(ACME = ACME,
              ADE = ADE,
              PM = PM,
              TE = TE,
              xm = xm,
              my = my))

}


##' Overall Indirect Effect for a set of mediator
##'
##' This function estimate the Overall Indirect Effect
##'
##' @param X exposure
##' @param m a set of mediator
##' @param Y outcome
##' @param C set of covariable (need to include latent factor)
##' @param boots number of bootstrap
##' @return
##' The estimate of Overall Indirect Effect
##'
##' @details
##' A novelty of HDMAX2 is to evaluate a cumulated indirect effect for all CpGs or AMRs identified in by max2() and .
##'
##' @export
##' @author Basile Jumentier
##' @examples
##'
##' # see vignette
est_oie <- function(X, m, Y, C, boots = 100) {

  # bootstrap
  acme_sum <- matrix(nrow = 1, ncol = boots)

  for (i in 1:ncol(acme_sum)) {
    samp <- sample(length(X), replace = T)

    # effet B m -> Y
    dat.1 <- data.frame(X, m, C)
    mod1 <- lm(Y[samp] ~ ., data = dat.1[samp, ])
    B <- as.data.frame(summary(mod1)$coeff[3:(ncol(m) + 2), ])

    # effet A X -> M
    mod2 <- lm(m[samp, ] ~ X[samp] + C[samp])
    A <- t(sapply(summary(mod2), function(x) x$coeff[2, ]))
    A <- data.frame(CpG = rownames(A), A)
    # A <- separate(A, CpG, c("0", "CpG"), " ")[, -1]

    colnames(B) <- c("B", "B_sd", "B_tv", "B_pv")
    colnames(A)[2:5] <- c("A", "A_sd", "A_tv", "A_pv")

    ab <- cbind(A, B)
    rownames(ab) <- NULL

    # effet A*B
    ab$AB <- ab$A * ab$B

    acme_sum[i] <- sum(ab$AB)
  }

  return(oie = as.vector(acme_sum))
}


