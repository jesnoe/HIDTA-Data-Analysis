library(spdep)

localmoran_perm2 <- function (x, listw, nsim = 499L, zero.policy = NULL, na.action = na.fail, 
                              alternative = "two.sided", mlvar = TRUE, spChk = NULL, xx=NULL,
                              adjust.x = FALSE, sample_Ei = TRUE, iseed = NULL) { # returns a list of res and df
  alternative <- match.arg(alternative, c("greater", 
                                          "less", "two.sided"))
  stopifnot(is.vector(x))
  if (!inherits(listw, "listw")) 
    stop(paste(deparse(substitute(listw)), "is not a listw object"))
  if (is.null(zero.policy)) 
    zero.policy <- get("zeroPolicy", envir = .spdepOptions)
  stopifnot(is.logical(zero.policy))
  if (!is.null(attr(listw$neighbours, "self.included")) && 
      attr(listw$neighbours, "self.included")) 
    stop("Self included among neighbours")
  if (is.null(spChk)) 
    spChk <- get.spChkOption()
  if (spChk && !chkIDs(x, listw)) 
    stop("Check of data and weights ID integrity failed")
  if (!is.numeric(x)) 
    stop(paste(deparse(substitute(x)), "is not a numeric vector"))
  NAOK <- deparse(substitute(na.action)) == "na.pass"
  x <- na.action(x)
  na.act <- attr(x, "na.action")
  rn <- attr(listw, "region.id")
  if (!is.null(na.act)) {
    subset <- !(1:length(listw$neighbours) %in% na.act)
    listw <- subset(listw, subset, zero.policy = zero.policy)
    excl <- class(na.act) == "exclude"
  }
  n <- length(listw$neighbours)
  if (n != length(x)) 
    stop("Different numbers of observations")
  res <- matrix(nrow = n, ncol = 9)
  gr <- punif((1:(nsim + 1))/(nsim + 1), 0, 1)
  ls <- rev(gr)
  ts <- (ifelse(gr > ls, ls, gr)) * 2
  if (alternative == "two.sided") {
    probs <- ts
    Prname <- "Pr(z != E(Ii))"
  }else if (alternative == "greater") {
    Prname <- "Pr(z > E(Ii))"
    probs <- gr
  }else {
    Prname <- "Pr(z < E(Ii))"
    probs <- ls
  }
  Prname_rank <- paste0(Prname, " Sim")
  Prname_sim <- "Pr(folded) Sim"
  colnames(res) <- c("Ii", "E.Ii", "Var.Ii", 
                     "Z.Ii", Prname, Prname_rank, Prname_sim, "Skewness", 
                     "Kurtosis")
  if (is.null(xx)) {
    if (adjust.x) {
      nc <- card(listw$neighbours) > 0L
      xx <- mean(x[nc], na.rm = NAOK)
    }
    else {
      xx <- mean(x, na.rm = NAOK)
    }
  }
  
  if (xx == "binary") {
    z <- x
    EIc <- -n*(z^2 * sapply(listw$weights, sum))
  } else {
    z <- x - xx
    EIc <- -(z^2 * sapply(listw$weights, sum))/((n - 1) * (sum(z * 
                                                                 z)/n))
  }
  
  
  lz <- lag.listw(listw, z, zero.policy = zero.policy, NAOK = NAOK)
  lbs <- c("Low", "High")
  quadr_ps <- interaction(cut(z, c(-Inf, 0, Inf), labels = lbs), 
                          cut(lz, c(-Inf, 0, Inf), labels = lbs), sep = "-")
  lx <- lag.listw(listw, x, zero.policy = zero.policy, NAOK = NAOK)
  lxx <- mean(lx, na.rm = NAOK)
  quadr <- interaction(cut(x, c(-Inf, xx, Inf), labels = lbs), 
                       cut(lx, c(-Inf, lxx, Inf), labels = lbs), sep = "-")
  xmed <- median(x, na.rm = NAOK)
  lxmed <- median(lx, na.rm = NAOK)
  quadr_med <- interaction(cut(x, c(-Inf, xmed, Inf), labels = lbs), 
                           cut(lx, c(-Inf, lxmed, Inf), labels = lbs), sep = "-")
  if (mlvar) {
    if (adjust.x) {
      s2 <- sum(z[nc]^2, na.rm = NAOK)/sum(nc)
    }
    else {
      s2 <- sum(z^2, na.rm = NAOK)/n
    }
  }else {
    if (adjust.x) {
      s2 <- sum(z[nc]^2, na.rm = NAOK)/(sum(nc) - 1)
    }
    else {
      s2 <- sum(z^2, na.rm = NAOK)/(n - 1)
    }
  }
  
  if (xx == "binary") {
    res[, 1] <- z * lz
  } else {
    res[, 1] <- (z/s2) * lz
  }
  cores <- get.coresOption()
  if (is.null(cores)) {
    parallel <- "no"
  }else {
    parallel <- ifelse(get.mcOption(), "multicore", 
                       "snow")
  }
  ncpus <- ifelse(is.null(cores), 1L, cores)
  cl <- NULL
  if (parallel == "snow") {
    cl <- get.ClusterOption()
    if (is.null(cl)) {
      parallel <- "no"
      warning("no cluster in ClusterOption, parallel set to no")
    }
  }
  if (!is.null(iseed)) {
    stopifnot(is.numeric(iseed))
    stopifnot(length(iseed) == 1L)
  }
  crd <- card(listw$neighbours)
  permI_int <- function(i, zi, z_i, crdi, wtsi, nsim, Ii, binary) {
    res_i <- rep(as.numeric(NA), 8)
    if (crdi > 0) {
      sz_i <- matrix(sample(z_i, size = crdi * nsim, replace = TRUE), 
                     ncol = crdi, nrow = nsim)
      lz_i <- sz_i %*% wtsi
      
      if (binary == "binary") {
        res_p <- zi * lz_i
      } else {
        res_p <- (zi/s2) * lz_i
      }
      res_i[1] <- mean(res_p)
      res_i[2] <- var(res_p)
      xrank <- rank(c(res_p, Ii))[(nsim + 1L)]
      res_i[5] <- xrank
      rnk0 <- as.integer(sum(res_p >= Ii))
      drnk0 <- nsim - rnk0
      rnk <- ifelse(drnk0 < rnk0, drnk0, rnk0)
      res_i[6] <- rnk0
      res_i[7] <- e1071::skewness(res_p)
      res_i[8] <- e1071::kurtosis(res_p)
    }
    res_i
  }
  lww <- listw$weights
  Iis <- res[, 1]
  if (parallel == "snow") {
    if (requireNamespace("parallel", quietly = TRUE)) {
      sI <- parallel::splitIndices(n, length(cl))
      env <- new.env()
      assign("z", z, envir = env)
      assign("crd", crd, envir = env)
      assign("lww", lww, envir = env)
      assign("nsim", nsim, envir = env)
      assign("Iis", Iis, envir = env)
      parallel::clusterExport(cl, varlist = c("z", 
                                              "crd", "lww", "nsim", "Iis"), 
                              envir = env)
      if (!is.null(iseed)) 
        parallel::clusterSetRNGStream(cl, iseed = iseed)
      oo <- parallel::clusterApply(cl, x = sI, fun = lapply, 
                                   function(i) {
                                     permI_int(i, z[i], z[-i], crd[i], lww[[i]], 
                                               nsim, Iis[i], xx)
                                   })
      out <- do.call("rbind", do.call("c", 
                                      oo))
      rm(env)
    }
    else {
      stop("parallel not available")
    }
  }else if (parallel == "multicore") {
    if (requireNamespace("parallel", quietly = TRUE)) {
      sI <- parallel::splitIndices(n, ncpus)
      oldRNG <- RNGkind()
      RNGkind("L'Ecuyer-CMRG")
      if (!is.null(iseed)) 
        set.seed(iseed)
      oo <- parallel::mclapply(sI, FUN = lapply, function(i) {
        permI_int(i, z[i], z[-i], crd[i], lww[[i]], nsim, 
                  Iis[i], xx)
      }, mc.cores = ncpus)
      RNGkind(oldRNG[1])
      out <- do.call("rbind", do.call("c", 
                                      oo))
    }
    else {
      stop("parallel not available")
    }
  }else {
    if (!is.null(iseed)) 
      set.seed(iseed)
    oo <- lapply(1:n, function(i) permI_int(i, z[i], z[-i], 
                                            crd[i], lww[[i]], nsim, Iis[i], xx))
    out <- do.call("rbind", oo)
  }
  if (sample_Ei) {
    res[, 2] <- out[, 1]
  }else res[, 2] <- EIc
  res[, 3] <- out[, 2]
  res[, 4] <- (res[, 1] - res[, 2])/sqrt(res[, 3])
  if (alternative == "two.sided") 
    res[, 5] <- 2 * pnorm(abs(res[, 4]), lower.tail = FALSE)
  else if (alternative == "greater") 
    res[, 5] <- pnorm(res[, 4], lower.tail = FALSE)
  else res[, 5] <- pnorm(res[, 4])
  res[, 6] <- probs[as.integer(out[, 5])]
  rnk0 <- as.integer(out[, 6]) # out[,6] is sum(res_p >= Ii) or R+1
  drnk0 <- nsim - rnk0
  if (alternative == "two.sided") {
    rnk <- ifelse(drnk0 < rnk0, drnk0, rnk0)
  }else if (alternative == "greater") {
    rnk <- rnk0
  }else {
    rnk <- drnk0
  }
  res[, 7] <- (rnk + 1)/(nsim + 1)
  res[, 8] <- out[, 7]
  res[, 9] <- out[, 8]
  if (!is.null(na.act) && excl) {
    res <- naresid(na.act, res)
  }
  if (!is.null(rn)) 
    rownames(res) <- rn
  attr(res, "call") <- match.call()
  if (!is.null(na.act)) 
    attr(res, "na.action") <- na.act
  df <- data.frame(quadr = quadr, quadr_med = quadr_med, 
                   quadr_ps = quadr_ps)
  class(res) <- c("localmoran", class(res))
  return(out)
}

localmoran_abs <- function (x, listw, nsim = 499L, zero.policy = NULL, na.action = na.fail, 
                              alternative = "two.sided", mlvar = TRUE, spChk = NULL, xx=NULL, # Put a specific xx (sample mean of comparing time period)
                              adjust.x = FALSE, sample_Ei = TRUE, iseed = NULL) {
  alternative <- match.arg(alternative, c("greater", 
                                          "less", "two.sided"))
  stopifnot(is.vector(x))
  if (!inherits(listw, "listw")) 
    stop(paste(deparse(substitute(listw)), "is not a listw object"))
  if (is.null(zero.policy)) 
    zero.policy <- get("zeroPolicy", envir = .spdepOptions)
  stopifnot(is.logical(zero.policy))
  if (!is.null(attr(listw$neighbours, "self.included")) && 
      attr(listw$neighbours, "self.included")) 
    stop("Self included among neighbours")
  if (is.null(spChk)) 
    spChk <- get.spChkOption()
  if (spChk && !chkIDs(x, listw)) 
    stop("Check of data and weights ID integrity failed")
  if (!is.numeric(x)) 
    stop(paste(deparse(substitute(x)), "is not a numeric vector"))
  NAOK <- deparse(substitute(na.action)) == "na.pass"
  x <- na.action(x)
  na.act <- attr(x, "na.action")
  rn <- attr(listw, "region.id")
  if (!is.null(na.act)) {
    subset <- !(1:length(listw$neighbours) %in% na.act)
    listw <- subset(listw, subset, zero.policy = zero.policy)
    excl <- class(na.act) == "exclude"
  }
  n <- length(listw$neighbours)
  if (n != length(x)) 
    stop("Different numbers of observations")
  res <- matrix(nrow = n, ncol = 9)
  gr <- punif((1:(nsim + 1))/(nsim + 1), 0, 1)
  ls <- rev(gr)
  ts <- (ifelse(gr > ls, ls, gr)) * 2
  if (alternative == "two.sided") {
    probs <- ts
    Prname <- "Pr(z != E(Ii))"
  }else if (alternative == "greater") {
    Prname <- "Pr(z > E(Ii))"
    probs <- gr
  }else {
    Prname <- "Pr(z < E(Ii))"
    probs <- ls
  }
  Prname_rank <- paste0(Prname, " Sim")
  Prname_sim <- "Pr(folded) Sim"
  colnames(res) <- c("Ii", "E.Ii", "Var.Ii", 
                     "Z.Ii", Prname, Prname_rank, Prname_sim, "Skewness", 
                     "Kurtosis")
  if (is.null(xx)) {
    if (adjust.x) {
      nc <- card(listw$neighbours) > 0L
      xx <- mean(x[nc], na.rm = NAOK)
    }
    else {
      xx <- mean(x, na.rm = NAOK)
    }
  }
  
  if (xx == "binary") {
    z <- x
    EIc <- -n*(z^2 * sapply(listw$weights, sum))
  } else {
    z <- x - xx
    EIc <- -(z^2 * sapply(listw$weights, sum))/((n - 1) * (sum(z * 
                                                                 z)/n))
  }
  
  
  lz <- lag.listw(listw, z, zero.policy = zero.policy, NAOK = NAOK)
  lbs <- c("Low", "High")
  quadr_ps <- interaction(cut(z, c(-Inf, 0, Inf), labels = lbs), 
                          cut(lz, c(-Inf, 0, Inf), labels = lbs), sep = "-")
  lx <- lag.listw(listw, x, zero.policy = zero.policy, NAOK = NAOK)
  lxx <- mean(lx, na.rm = NAOK)
  quadr <- interaction(cut(x, c(-Inf, xx, Inf), labels = lbs), 
                       cut(lx, c(-Inf, lxx, Inf), labels = lbs), sep = "-")
  xmed <- median(x, na.rm = NAOK)
  lxmed <- median(lx, na.rm = NAOK)
  quadr_med <- interaction(cut(x, c(-Inf, xmed, Inf), labels = lbs), 
                           cut(lx, c(-Inf, lxmed, Inf), labels = lbs), sep = "-")
  if (mlvar) {
    if (adjust.x) {
      s2 <- sum(z[nc]^2, na.rm = NAOK)/sum(nc)
    }
    else {
      s2 <- sum(z^2, na.rm = NAOK)/n
    }
  }else {
    if (adjust.x) {
      s2 <- sum(z[nc]^2, na.rm = NAOK)/(sum(nc) - 1)
    }
    else {
      s2 <- sum(z^2, na.rm = NAOK)/(n - 1)
    }
  }
  
  if (xx == "binary") {
    res[, 1] <- z * lz
  } else {
    res[, 1] <- (z/s2) * lz
  }
  cores <- get.coresOption()
  if (is.null(cores)) {
    parallel <- "no"
  }else {
    parallel <- ifelse(get.mcOption(), "multicore", 
                       "snow")
  }
  ncpus <- ifelse(is.null(cores), 1L, cores)
  cl <- NULL
  if (parallel == "snow") {
    cl <- get.ClusterOption()
    if (is.null(cl)) {
      parallel <- "no"
      warning("no cluster in ClusterOption, parallel set to no")
    }
  }
  if (!is.null(iseed)) {
    stopifnot(is.numeric(iseed))
    stopifnot(length(iseed) == 1L)
  }
  crd <- card(listw$neighbours)
  permI_int <- function(i, zi, z_i, crdi, wtsi, nsim, Ii, binary) {
    res_i <- rep(as.numeric(NA), 8)
    if (crdi > 0) {
      sz_i <- matrix(sample(z_i, size = crdi * nsim, replace = TRUE), 
                     ncol = crdi, nrow = nsim)
      lz_i <- sz_i %*% wtsi
      
      if (binary == "binary") {
        res_p <- zi * lz_i
      } else {
        res_p <- (zi/s2) * lz_i
      }
      res_i[1] <- mean(res_p)
      res_i[2] <- var(res_p)
      xrank <- rank(c(res_p, Ii))[(nsim + 1L)]
      res_i[5] <- xrank
      rnk0 <- as.integer(sum(res_p >= Ii))
      drnk0 <- nsim - rnk0
      rnk <- ifelse(drnk0 < rnk0, drnk0, rnk0)
      res_i[6] <- rnk0
      res_i[7] <- e1071::skewness(res_p)
      res_i[8] <- e1071::kurtosis(res_p)
    }
    res_i
  }
  lww <- listw$weights
  Iis <- res[, 1]
  if (parallel == "snow") {
    if (requireNamespace("parallel", quietly = TRUE)) {
      sI <- parallel::splitIndices(n, length(cl))
      env <- new.env()
      assign("z", z, envir = env)
      assign("crd", crd, envir = env)
      assign("lww", lww, envir = env)
      assign("nsim", nsim, envir = env)
      assign("Iis", Iis, envir = env)
      parallel::clusterExport(cl, varlist = c("z", 
                                              "crd", "lww", "nsim", "Iis"), 
                              envir = env)
      if (!is.null(iseed)) 
        parallel::clusterSetRNGStream(cl, iseed = iseed)
      oo <- parallel::clusterApply(cl, x = sI, fun = lapply, 
                                   function(i) {
                                     permI_int(i, z[i], z[-i], crd[i], lww[[i]], 
                                               nsim, Iis[i], xx)
                                   })
      out <- do.call("rbind", do.call("c", 
                                      oo))
      rm(env)
    }
    else {
      stop("parallel not available")
    }
  }else if (parallel == "multicore") {
    if (requireNamespace("parallel", quietly = TRUE)) {
      sI <- parallel::splitIndices(n, ncpus)
      oldRNG <- RNGkind()
      RNGkind("L'Ecuyer-CMRG")
      if (!is.null(iseed)) 
        set.seed(iseed)
      oo <- parallel::mclapply(sI, FUN = lapply, function(i) {
        permI_int(i, z[i], z[-i], crd[i], lww[[i]], nsim, 
                  Iis[i], xx)
      }, mc.cores = ncpus)
      RNGkind(oldRNG[1])
      out <- do.call("rbind", do.call("c", 
                                      oo))
    }
    else {
      stop("parallel not available")
    }
  }else {
    if (!is.null(iseed)) 
      set.seed(iseed)
    oo <- lapply(1:n, function(i) permI_int(i, z[i], z[-i], 
                                            crd[i], lww[[i]], nsim, Iis[i], xx))
    out <- do.call("rbind", oo)
  }
  if (sample_Ei) {
    res[, 2] <- out[, 1]
  }else res[, 2] <- EIc
  res[, 3] <- out[, 2]
  res[, 4] <- (res[, 1] - res[, 2])/sqrt(res[, 3])
  if (alternative == "two.sided") 
    res[, 5] <- 2 * pnorm(abs(res[, 4]), lower.tail = FALSE)
  else if (alternative == "greater") 
    res[, 5] <- pnorm(res[, 4], lower.tail = FALSE)
  else res[, 5] <- pnorm(res[, 4])
  res[, 6] <- probs[as.integer(out[, 5])]
  rnk0 <- as.integer(out[, 6]) # out[,6] is sum(res_p >= Ii) or R+1
  drnk0 <- nsim - rnk0
  if (alternative == "two.sided") {
    rnk <- ifelse(drnk0 < rnk0, drnk0, rnk0)
  }else if (alternative == "greater") {
    rnk <- rnk0
  }else {
    rnk <- drnk0
  }
  res[, 7] <- (rnk + 1)/(nsim + 1)
  res[, 8] <- out[, 7]
  res[, 9] <- out[, 8]
  if (!is.null(na.act) && excl) {
    res <- naresid(na.act, res)
  }
  if (!is.null(rn)) 
    rownames(res) <- rn
  attr(res, "call") <- match.call()
  if (!is.null(na.act)) 
    attr(res, "na.action") <- na.act
  df <- data.frame(quadr = quadr, quadr_med = quadr_med, 
                   quadr_ps = quadr_ps)
  class(res) <- c("localmoran", class(res))
  return(cbind(res, df))
}

perm_index <- function(z_, crd_, j) {
  n_counties <- length(z_)
  return(sample((1:n_counties)[-j], size = crd_, replace = F))
}
permI_dist <- function(zi, z_i, crdi, wtsi, nsim, Ii, replacement) {
  if (replacement==T) {
    sz_i_w_rep <- matrix(sample(z_i, size = crdi * nsim, replace = T), 
                         ncol = crdi, nrow = nsim)
    lz_i_w_rep <- sz_i_w_rep %*% wtsi
    I_perm <- (zi/s2) * lz_i_w_rep
  } else {
    sz_i_wo_rep <- matrix(rep(0,crdi), ncol=crdi)
    for (i in 1:nsim){
      sz_i_wo_rep <- rbind(sz_i_wo_rep, matrix(sample(z_i, size=crdi, replace=F),
                                               ncol = crdi))
    }
    sz_i_wo_rep <- sz_i_wo_rep[-1,]
    lz_i_w_rep <- sz_i_wo_rep %*% wtsi
    I_perm <- (zi/s2) * lz_i_w_rep
  }
  I_perm <- c(I_perm, Ii)
  result <- data.frame(I_perm=I_perm,
                       observation=c(rep(0,nsim), 1))
  return(result)
}