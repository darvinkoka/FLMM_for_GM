#utils_fmm
Indic <- function(x,y){
  return(y<=x)
}

add_regressors <- function(df, h_val, Mh_val, Mref_val, period) {
  R <- sqrt(df$djb^2 + h_val^2)
  df[[sprintf("x1_%.2f", period)]] <- (df$mw - Mh_val) * as.integer(df$mw <= Mh_val)
  df[[sprintf("x2_%.2f", period)]] <- (df$mw - Mh_val) * as.integer(df$mw >= Mh_val)
  df[[sprintf("x5_%.2f", period)]] <- (df$mw - Mref_val) * log10(R)
  df[[sprintf("x6_%.2f", period)]] <- log10(R)
  df[[sprintf("x7_%.2f", period)]] <- R
  df$x3 <- ifelse(df$sof == "SS", 1, 0)
  df$x4 <- ifelse(df$sof == "TF", 1, 0)
  df$x8 <- log10(pmin(df$vs30, 1500) / 800)
  return(df)
}

smooth_fd <- function(X, cv, grid, tp) {
  basis <- create.bspline.basis(rangeval = range(tp), breaks = tp, norder = if (cv) 4 else 3)
  if (!cv) {
    return(smooth.basis(tp, X, fdPar(basis, 1, 0.1))$fd)
  }
  gcv <- sapply(grid, function(l) sum(smooth.basis(tp, X, fdPar(basis, 2, l))$gcv))
  smooth.basis(tp, X, fdPar(basis, 2, grid[which.min(gcv)]))$fd
}

build_xlist_fda <- function(data, T.emp, Mh, Mref, h, fixed.idx = c(1, 2),
                            lambda.grid = 10^-seq(-1, 12, by = 1),
                            lambda.grid.mref = 10^-seq(-1, 10, by = 1)) {
  
  ord <- order(T.emp); T.emp <- T.emp[ord]; Mh <- Mh[ord]; Mref <- Mref[ord]; h <- h[ord]
  tp <- log10(T.emp); tp[!is.finite(tp)] <- -2.5
  cols <- names(data)
  
  fcol <- cols[grepl("^x[0-9]+_", cols)]              # es. "x1_0.3" -> funzione del periodo
  fidx <- as.integer(sub("^x([0-9]+)_.*$", "\\1", fcol))
  fper <- as.numeric(sub("^x[0-9]+_", "", fcol))
  
  scol <- cols[grepl("^x[0-9]+$", cols)]               # es. "x3"    -> scalare (nessun suffisso)
  sidx <- as.integer(sub("^x", "", scol))
  
  if (length(fcol) == 0 && length(scol) == 0) stop("Nessuna colonna x<k>[_<T>] trovata.")
  
  xlist <- list(rep(1, nrow(data)))
  for (k in sort(unique(c(fidx, sidx)))) {
    if (k %in% sidx) {
      xlist[[length(xlist) + 1]] <- data[[scol[sidx == k]]]
    } else {
      sel <- sapply(T.emp, function(t) which.min(abs(fper[fidx == k] - t)))
      Xk  <- t(as.matrix(data[, fcol[fidx == k][sel], drop = FALSE]))  # p x n, ordine T.emp
      
      # Passiamo sia lambda.grid che tp
      xlist[[length(xlist) + 1]] <- smooth_fd(Xk, cv = !(k %in% fixed.idx), grid = lambda.grid, tp = tp)
    }
  }
  
  attr(xlist, "params.fd") <- list(
    Mh.fd   = smooth_fd(Mh, cv = FALSE, grid = lambda.grid, tp = tp),
    Mref.fd = smooth_fd(Mref, cv = TRUE, grid = lambda.grid.mref, tp = tp),
    h.fd    = smooth_fd(h, cv = TRUE, grid = lambda.grid, tp = tp)
  )
  
  xlist
}
