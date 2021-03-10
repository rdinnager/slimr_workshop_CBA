basic.stats <- function (df_basic.stats , diploid = TRUE, digits = 4) {
  loc.names <- names(df_basic.stats )[-1]
  if (length(table(df_basic.stats [, 1])) < 2)
    df_basic.stats [dim(df_basic.stats )[1] + 1, 1] <- df_basic.stats [dim(df_basic.stats )[1], 1] + 1
  if (dim(df_basic.stats )[2] == 2)
    df_basic.stats  <- data.frame(df_basic.stats , dummy.loc = df_basic.stats [, 2])
  p <- pop.freq(df_basic.stats , diploid)
  n <- t(ind.count(df_basic.stats ))
  if (diploid) {
    dum <- getal.b(df_basic.stats [,-1])
    Ho <- dum[, , 1] == dum[, , 2]
    sHo <-
      (1 - t(apply(
        Ho, 2, fun <- function(x)
          tapply(x, df_basic.stats [, 1], mean, na.rm = TRUE)
      )))
    mHo <- apply(sHo, 1, mean, na.rm = TRUE)
  }
  else {
    sHo <- NA
    mHo <- NA
  }
  sp2 <-
    lapply(p, fun <-
             function(x)
               apply(x, 2, fun2 <- function(x)
                 sum(x ^ 2)))
  sp2 <-
    matrix(unlist(sp2), nrow = dim(df_basic.stats [,-1])[2], byrow = TRUE)
  if (diploid) {
    Hs <- (1 - sp2 - sHo / 2 / n)
    Hs <- n / (n - 1) * Hs
    Fis = 1 - sHo / Hs
  }
  else {
    Hs <- n / (n - 1) * (1 - sp2)
    Fis <- NA
  }
  np <- apply(n, 1, fun <- function(x)
    sum(!is.na(x)))
  mn <- apply(n, 1, fun <- function(x) {
    np <- sum(!is.na(x))
    np / sum(1 / x[!is.na(x)])
  })
  msp2 <- apply(sp2, 1, mean, na.rm = TRUE)
  mp <-
    lapply(p, fun <- function(x)
      apply(x, 1, mean, na.rm = TRUE))
  mp2 <- unlist(lapply(mp, fun1 <- function(x)
    sum(x ^ 2)))
  if (diploid) {
    mHs <- mn / (mn - 1) * (1 - msp2 - mHo / 2 / mn)
    Ht <- 1 - mp2 + mHs / mn / np - mHo / 2 / mn / np
    mFis = 1 - mHo / mHs
  }
  else {
    mHs <- mn / (mn - 1) * (1 - msp2)
    Ht <- 1 - mp2 + mHs / mn / np
    mFis <- NA
  }
  Dst <- Ht - mHs
  Dstp <- np / (np - 1) * Dst
  Htp = mHs + Dstp
  Fst = Dst / Ht
  Fstp = Dstp / Htp
  Dest <- Dstp / (1 - mHs)
  res <- data.frame(cbind(mHo, mHs, Ht, Dst, Htp, Dstp, Fst,
                          Fstp, mFis, Dest))
  names(res) <- c("Ho",
                  "Hs",
                  "Ht",
                  "Dst",
                  "Htp",
                  "Dstp",
                  "Fst",
                  "Fstp",
                  "Fis",
                  "Dest")
  if (diploid) {
    rownames(sHo) <- loc.names
    rownames(Fis) <- loc.names
  }
  is.na(res) <- do.call(cbind, lapply(res, is.infinite))
  overall <- apply(res, 2, mean, na.rm = TRUE)
  overall[7] <- overall[4] / overall[3]
  overall[8] <- overall[6] / overall[5]
  overall[9] <- 1 - overall[1] / overall[2]
  overall[10] <- overall[6] / (1 - overall[2])
  names(overall) <- names(res)
  if (!diploid) {
    overall[-2] <- NA
  }
  all.res <-
    list(
      n.ind.samp = n,
      pop.freq = lapply(p, round, digits),
      Ho = round(sHo, digits),
      Hs = round(Hs, digits),
      Fis = round(Fis, digits),
      perloc = round(res, digits),
      overall = round(overall, digits)
    )
  class(all.res) <- "basic.stats"
  all.res
}

pop.freq <- function (df_pop.freq, diploid = TRUE) {
  nbpop <- length(unique(df_pop.freq[, 1]))
  nbloc <- dim(df_pop.freq)[2] - 1
  if (diploid) {
    df_pop.freq <-
      data.frame(df_pop.freq, dummy.loc = (sample(
        9, replace = TRUE, size = dim(df_pop.freq)[1]
      ) + 100) * 1001)
    df_pop.freq <- getal(df_pop.freq)[,-2]
  }
  else {
    df_pop.freq <-
      data.frame(df_pop.freq, dummy.loc = sample(9, replace = TRUE, size = dim(df_pop.freq)[1]) + 100)
  }
  freq <- function(x) {
    if (is.character(x))
      dum <- table(factor(x), df_pop.freq[, 1])
    else
      dum <- (table(x, df_pop.freq[, 1]))
    eff <- colSums(dum, na.rm = TRUE)
    freq <- sweep(dum, 2, eff, FUN = "/")
  }
  ndat <- df_pop.freq[,-1]
  all.freq <- apply(ndat, 2, freq)
  all.freq <- all.freq[-(nbloc + 1)]
  return(all.freq)
}

getal <- function (df_getal) {
  x <- dim(df_getal)
  if (max(df_getal[, 2], na.rm = TRUE) > 1e+06)
    stop("allele encoding with 3 digits maximum")
  if (max(df_getal[, 2], na.rm = TRUE) < 1e+06)
    modulo <- 1000
  if (max(df_getal[, 2], na.rm = TRUE) < 10000) {
    if (min(df_getal[, 2] %/% 100, na.rm = TRUE) >= 10 &
        max(df_getal[, 2] %% 100, na.rm = TRUE) < 10)
      modulo <- 1000
    else
      modulo <- 100
  }
  if (max(df_getal[, 2], na.rm = TRUE) < 100)
    modulo <- 10
  firstal <- df_getal[,-1] %/% modulo
  secal <- df_getal[,-1] %% modulo
  ind <- vector(length = 0)
  nbpop <- length(unique(df_getal[, 1]))
  for (i in sort(unique(df_getal[, 1]))) {
    dum <- 1:sum(df_getal[, 1] == i)
    if (i == 1)
      ind <- dum
    else
      ind <- c(ind, dum)
  }
  ind <- rep(ind, 2)
  if (x[2] == 2)
    data.al <-
    data.frame(
      pop = rep(df_getal[, 1], 2),
      ind = ind,
      al = c(firstal, secal)
    )
  else
    data.al <-
    data.frame(
      pop = rep(df_getal[, 1], 2),
      ind = ind,
      al = rbind(firstal, secal)
    )
  names(data.al)[-c(1:2)] <- names(df_getal)[-1]
  return(data.al)
}

ind.count <- function (df_ind.count) {
  dum <- function(x) {
    a <- which(!is.na(x))
    tapply(x[a], df_ind.count[a, 1], length)
  }
  df_ind.count[, 1] <- factor(df_ind.count[, 1])
  apply(df_ind.count[,-1], 2, dum)
}

getal.b <- function (df_getal.b) {
  x <- dim(df_getal.b)
  if (max(df_getal.b[, 2], na.rm = TRUE) > 1e+06)
    stop("allele encoding with 3 digits maximum")
  if (max(df_getal.b[, 2], na.rm = TRUE) < 1e+06)
    modulo <- 1000
  if (max(df_getal.b[, 2], na.rm = TRUE) < 10000) {
    if (min(df_getal.b[, 2] %/% 100, na.rm = TRUE) >= 10 &
        max(df_getal.b[, 2] %% 100, na.rm = TRUE) < 10)
      modulo <- 1000
    else
      modulo <- 100
  }
  if (max(df_getal.b[, 2], na.rm = TRUE) < 100)
    modulo <- 10
  firstal <- df_getal.b %/% modulo
  secal <- df_getal.b %% modulo
  y <- array(dim = c(x, 2))
  y[, , 1] <- as.matrix(firstal)
  y[, , 2] <- as.matrix(secal)
  return(y)
}

# Counts the number of copies of the different alleles at each locus and population
allele.count <- function (df_allele.count, diploid = TRUE){
    if (diploid) 
        df_allele.count <- getal(df_getal=df_allele.count)[, -2]
    nbloc <- dim(df_allele.count)[2] - 1
    fun <- function(x) {
        dum <- table(x, df_allele.count[, 1])
    }
    all.count <- apply(as.data.frame(df_allele.count[, -1]), 2, fun)
    if (!is.list(all.count)) {
        my.dim <- dim(table(df_allele.count[, 2], df_allele.count[, 1]))
        my.names <- dimnames(all.count)[[2]]
        dum <- list(nbloc)
        for (i in seq_along(my.names)) {
            dum[[i]] <- matrix(all.count[, i], nrow = my.dim[1])
        }
        names(dum) <- my.names
        all.count <- dum
    }
    return(all.count)
}


nb.alleles <- function (df_nb.alleles, diploid = TRUE){
    dum <- allele.count(df_nb.alleles, diploid)
    if (is.list(dum)) 
        res <- lapply(dum, fun <- function(x) apply(x, 2, fun2 <- function(x) sum(x > 
            0)))
    else res <- apply(dum, 2, fun2 <- function(x) sum(x > 0))
    res <- matrix(unlist(res), nrow = dim(df_nb.alleles)[2] - 1, byrow = TRUE)
    rownames(res) <- names(df_nb.alleles)[-1]
    return(res)
}