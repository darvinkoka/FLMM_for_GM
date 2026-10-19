rm(list=ls())
graphics.off()
cat("\014")
gc()

setwd('~/Documents/PhD/PROJECTS/FLMM_for_GM')
source("Code/utils_fmm.R")
source("Code/weighted-analysis.R")
df_emp<- read.csv("Data/ITA18_SA_flatfile.csv", sep = ";")
library(readxl)
load("Data/ITA18_parameters.RData")

#DATA PREPROCESSING
SA.names.emp<-names(df_emp)[grep("rotD50_T",names(df_emp))][-1]

T.emp <- sub(".*_T", "", SA.names.emp)
T.emp <- as.numeric(sub("_", ".", T.emp))
Mh <- ITA18.parameters$Mh.vec[-1]
h <- ITA18.parameters$h.vec[-1]
Mref <- ITA18.parameters$Mref.vec[-1]
ITA18 <- data.frame(
  event = as.factor(df_emp$event_id),
  mw = df_emp$Mw,
  djb = df_emp$JB_dist,
  vs30 = df_emp$vs30_m_sec,
  sof = df_emp$fm_type_code,
  Station = as.factor(df_emp$station_code)
)
ITA18$djb[is.na(ITA18$djb)]   <- df_emp$epi_dist[is.na(ITA18$djb)]
ITA18$vs30[is.na(ITA18$vs30)] <- df_emp$vs30_m_sec_WA[is.na(ITA18$vs30)]

for(i in 1:length(T.emp)){
  ITA18 <- add_regressors(ITA18, h[i], Mh[i], Mref[i], T.emp[i])
}
SA_data <- df_emp[, SA.names.emp]
SA_data <- cbind(df_emp$rotD50_pga,SA_data)
names(SA_data) <- paste0("sa_", T.emp)
ITA18 <- cbind(ITA18, SA_data)
ITA18 <- ITA18[ITA18$sof %in% c("SS", "TF", "NF"), ]
for(t in T.emp) {
  ITA18[[paste0("logsa_", t)]] <- log10(ITA18[[paste0("sa_", t)]])
}

logsa.names <- names(ITA18)[grep("^log",names(ITA18))]
logsa <- ITA18[,logsa.names]
na_rows <- sapply(1:nrow(logsa), function(i){
  any(is.na(logsa[i,]))
})

recurrence_station <- sapply(unique(ITA18$Station), function(s){
  sum(ITA18$Station == s)
})
too_few <- sapply(ITA18$Station, function(s){
  recurrence_station[which(s==unique(ITA18$Station))]<10
})

ITA18 <- ITA18[!na_rows & !too_few,]
logsa <- logsa[!na_rows & !too_few,]
#BUILD FUNCTIONAL DATASET

#da cambiare
ITA18 <- read.csv("Data/ITA18_source-path.csv")
stopifnot(nrow(ITA18) == nrow(logsa))
library(fda)
xlist <- build_xlist_fda(ITA18, T.emp, Mh, Mref, h)
xlist[[4]] <- as.numeric(xlist[[4]])
xlist[[5]] <- as.numeric(xlist[[5]])
bsplines <- smooth_fd(t(as.matrix(logsa)), TRUE, grid = 10^-seq(-1, 12, by = 1), tp = log10(T.emp))
logsa.fd <- bsplines[[2]]

#ESTIMATE FIXED PARAMETERS - Multistart EAASS algorithm
fPar  <- fdPar(fdobj=logsa.fd$basis, Lfdobj=2, lambda=bsplines[[1]])
blist.default <- list(fPar, fPar, fPar, fPar, fPar, fPar, fPar, fPar, fPar)


calibrate <- FALSE
if(calibrate){
  library(progress)
  # Algorithm parameters
  n.iter <- 15
  B <- 5
  P.size <- 20
  n.par <- length(xlist)
  r=0.5
  perturbation.vec <- c(0.5,2)
  
  n.multistart <- 5
  seed.in <- 22
  
  for(b in 1:n.multistart){ #b <- 1
    seed <- seed.in + b
    print(paste0("Iteration ",b, " of ", n.multistart))
    
    set.seed(seed)
    # Initial population of combinations of tuning parameters
    exp.range <- seq(-5,2,by=1)
    lambda.range <- 10^exp.range
    P.initial <- matrix(data=NA, nrow=P.size, ncol=n.par)
    
    for(i in 1:P.size){
      P.initial[i,] <- sample(lambda.range, size=n.par, replace=TRUE)
    }
    
    # Compute the errors for the initial P via a B-fold CV
    # Errors corresponding to each row of P.initial are saved in vector V.initial
    (Start.Time <- Sys.time())
    V.initial <- c()
    pb <- progress_bar$new(total=(P.size), format = "  computing [:bar] :percent eta: :eta")
    for(i in 1:P.size){# i = 1
      # define blist
      blist <- blist.default
      for(reg in 1:n.par){
        blist[[reg]]$lambda <- P.initial[i,reg]
      }
      
      # evaluate the MSE via B-fold CV
      MSE_pw <- pwMSE(curves    = t(as.matrix(logsa)),
                      curves.fd = logsa.fd,
                      xlist     = xlist,
                      t.points  = log10(T.emp),
                      events    = ITA18$event,
                      blist     = blist,
                      B         = 10,
                      wgts.fd   = NULL,
                      wgts.flag = FALSE,
                      verbose   = TRUE)$MSE_t
      V.initial[i] <- sum(MSE_pw^2)/length(MSE_pw)
      
      pb$tick()
    }
    
    End.Time <- Sys.time()
    round(End.Time - Start.Time, 2)
    #beep()
    
    save(P.initial, V.initial, file=paste0("Results/Calibration/Initial_blist_select_seed",seed,".RData"))
    #load(paste0("data/calibration/Initial_blist_select_seed",seed,".RData"))
    
    V <- V.initial
    P <- P.initial
    
    #set.seed(b)
    (Start.Time <- Sys.time())
    for(k in 1:n.iter){
      print(paste0("Iteration: ",k))
      
      V.sort <- sort(V, decreasing=FALSE)                 # keep the best half (lowest CV error)
      idxs.in <- which(V %in% V.sort[1:(r*P.size)])
      if(length(idxs.in)>(r*P.size)){ #solve ties by randomly selecting only r*P.size indexes
        idxs.in <- sample(idxs.in, size=(r*P.size), replace=FALSE)
      }
      Q <- P[idxs.in,]
      Z <- V[idxs.in]
      
      perturbation <- sample(perturbation.vec, size=(r*P.size*n.par), replace=TRUE)
      pert_mat <- matrix(perturbation, nrow=(r*P.size), ncol=n.par)
      Q_prime <- Q * pert_mat
      
      # ---
      temp <- Q_prime
      idxs.ext <- which(Q_prime < 1e-5 | Q_prime > 1e2)
      temp[idxs.ext] <- sample(lambda.range, size=length(idxs.ext), replace=TRUE)
      Q_prime <- temp
      # ---
      
      pb <- progress_bar$new(total=(r*P.size), format = "  computing [:bar] :percent eta: :eta")
      Z_prime <- c()
      for(i in 1:(r*P.size)){
        # define blist
        blist <- blist.default
        for(reg in 1:n.par){
          blist[[reg]]$lambda <- Q_prime[i,reg]
        }
        
        # evaluate the MSE via B-fold CV
        MSE_pw <- pwMSE(curves    = t(as.matrix(logsa)),
                        curves.fd = logsa.fd,
                        xlist     = xlist,
                        t.points  = log10(T.emp),
                        events    = ITA18$event,
                        blist     = blist,
                        B         = 10,
                        wgts.fd   = NULL,
                        wgts.flag = FALSE,
                        verbose   = TRUE)$MSE_t
        Z_prime[i] <- sum(MSE_pw^2)/length(MSE_pw)
        pb$tick()
      }
      
      P.in <- P[idxs.in,]
      P <- rbind(P.in, Q_prime)
      V.in <- V[idxs.in]
      V <- c(V.in, Z_prime)
    }
    End.Time <- Sys.time()
    round(End.Time - Start.Time, 2)
    #beep()
    
    save(P, V, file=paste0("Results/Calibration/P_V_blist_select_seed",seed,".RData"))
    
    best <- which.min(V)
    lambda.opt <- P[best,]
    blist <- blist.default
    for(reg in 1:n.par){
      blist[[reg]]$lambda <- lambda.opt[reg]
    }
    
    save(blist, file=paste0('Data/Calibration/blist_EAASS_seed',seed,'.RData'))
  }
} else {
  ## Load results and define blist
  n.multistart <- 5
  n.par <- length(xlist)
  P.full <- matrix(data=NA, nrow=0, ncol=n.par)
  V.full <- c()
  seed.in <- 22
  for(b in 1:n.multistart){ # b <- 1
    seed <- seed.in + b
    load(paste0("Results/Calibration/P_V_blist_select_seed",seed,".RData"))
    
    P.full <- rbind(P.full, P)
    V.full <- c(V.full, V)
  }
  
  #plot(1:length(V.full), V.full, type='l')
  #dev.off()
  
  best <- which.min(V.full)
  lambda.opt <- P.full[best,]
  
  #' lambda.opt contains the values reported in Table 1 of the Supplementary Material.
  lambda.opt
  blist <- blist.default
  for(reg in 1:n.par){
    blist[[reg]]$lambda <- lambda.opt[reg]
  }
  
  save(blist, file='Results/Calibration/blist_EAASS_multistart.RData')
  
}

load('Results/Calibration/blist_EAASS_multistart.RData')
fit.fregress <- fRegress(logsa.fd, xlist, blist)
beta.fd.old <- fit.fregress$betaestlist

y.hat <- fit.fregress$yhatfdobj
res.fd <- fit.fregress$yfdobj - fit.fregress$yhatfdobj

################################################################################
# FAMM con refund::pffr
################################################################################
library(refund)
library(mgcv)

## --- dati in formato pffr ------------------------------------------------------
n  <- nrow(logsa); tt <- log10(T.emp)
Xl <- lapply(xlist, function(x) if (is.fd(x)) t(eval.fd(tt, x)) else matrix(x, n, length(tt)))
names(Xl) <- paste0("x", seq_along(Xl))

df <- c(list(Y       = as.matrix(logsa),
             Station = as.factor(ITA18$Station),
             Event   = as.factor(ITA18$event),
             Source  = as.factor(ITA18$source_term)),
        Xl[-1])                                            # x1 = intercetta, gestita da pffr

f <- as.formula(paste("Y ~", paste(c(names(Xl)[-1],
                                     's(Station, bs="re")',
                                     's(Event, bs="re")',
                                     's(Source, bs="re")'), collapse = " + ")))

## --- base in t: B-spline cubiche, un nodo per periodo, penalità su f'' (lambda via REML)
yindex <- tt                                               # nome fisso -> variabile "yindex.vec" in pffr
bsy <- list(bs = "bs", k = length(tt) + 2, m = c(3, 2))
kn <- list(yindex.vec = c(tt[1] - (3:1) * diff(tt[1:2]), tt, tail(tt, 1) + (1:3) * diff(tail(tt, 2))))

## --- fit -----------------------------------------------------------------------
fit.pffr <- pffr(f, yind = yindex, data = df, algorithm = "bam", nthreads = 5,
                 bs.yindex = bsy, bs.int = bsy, knots = kn)
stopifnot(isTRUE(all.equal(fit.pffr$smooth[[2]]$knots[3 + seq_along(tt)], tt)))   # nodi = periodi

## --- beta come oggetti fd (stessa stima, base fda su [-2,1]) + covarianza --------
sm <- fit.pffr$smooth; th <- fit.pffr$coefficients; Vp <- fit.pffr$Vp
rng <- range(tt)
beta.fd <- C.fd <- vector("list", 9)

for (j in 1:9) {
  s   <- sm[[j]]
  br  <- c(rng[1], s$knots[s$knots > rng[1] & s$knots < rng[2]], rng[2])
  bas <- create.bspline.basis(rng, norder = s$m[1] + 1, breaks = br)
  x   <- seq(rng[1], rng[2], length.out = 10 * bas$nbasis)
  nd  <- setNames(data.frame(x), s$term[1]); if (s$by != "NA") nd[[s$by]] <- 1
  X   <- PredictMat(s, nd)
  idx <- s$first.para:s$last.para
  if (j == 1) { X <- cbind(1, X); idx <- c(1, idx) }       # beta_1 = (Intercept) + s(t)
  M   <- qr.solve(eval.basis(x, bas), X)                   # theta_hat -> coeff. fda
  beta.fd[[j]] <- fd(M %*% th[idx], bas)
  C.fd[[j]]    <- M %*% Vp[idx, idx] %*% t(M)
}
names(beta.fd) <- names(C.fd) <- paste0("beta", 1:9)

## --- plot beta: pffr (blu, +-2 se) vs fRegress (nero) ------------------------------
tg <- seq(rng[1], rng[2], length.out = 200)
par(mfrow = c(3, 3))
for (j in 1:9) {
  Phi <- eval.basis(tg, beta.fd[[j]]$basis)
  b   <- eval.fd(tg, beta.fd[[j]])
  se  <- sqrt(rowSums((Phi %*% C.fd[[j]]) * Phi))
  lo  <- b - 2 * se; up <- b + 2 * se
  old <- eval.fd(tg, beta.fd.old[[j]]$fd)
  plot(tg, old, type = "l", ylim = range(old, lo, up), xlab = "log10 T", ylab = "", main = names(beta.fd)[j])
  polygon(c(tg, rev(tg)), c(lo, rev(up)), col = adjustcolor("blue", .2), border = NA)
  lines(tg, b, col = "blue", lwd = 2)
}

## --- componenti di variabilità ---------------------------------------------------
tr <- predict(fit.pffr, type = "terms"); nm <- names(tr)
S  <- tr[[grep("Station", nm)]][!duplicated(df$Station), ]
E  <- tr[[grep("Event",   nm)]][!duplicated(df$Event), ]
L  <- tr[[grep("Source",  nm)]][!duplicated(df$Source), ]
R  <- residuals(fit.pffr)

cols <- c("black", "darkblue", "darkgreen", "darkorange")
par(mfrow = c(1, 1))
matplot(T.emp, cbind(apply(S, 2, sd), apply(E, 2, sd), apply(R, 2, sd), apply(L, 2, sd)),
        type = "l", lty = 1, col = cols, log = "x", xlab = "T [s]", ylab = "sd")
legend("topright", c("dS2S", "dB", "dWS", "dL2L"), col = cols, lty = 1, y.intersp = 0.4)

################################################################################

################################################################################
# Modello scalare: un lmer per periodo
################################################################################
library(lme4)

## --- fit per periodo -------------------------------------------------------------
nT <- length(tt)
fs <- as.formula(paste("y ~", paste(names(Xl)[-1], collapse = " + "),
                       "+ (1|Station) + (1|Event) + (1|Source)"))

fit.lmer <- lapply(1:nT, function(k) {
  dk <- data.frame(y = df$Y[, k], sapply(Xl[-1], function(x) x[, k]),
                   Station = df$Station, Event = df$Event, Source = df$Source)
  lmer(fs, data = dk, REML = TRUE)
})

## --- coefficienti e SE (nT x 9) ---------------------------------------------------
B.lmer  <- t(sapply(fit.lmer, fixef))
SE.lmer <- t(sapply(fit.lmer, function(m) sqrt(diag(vcov(m)))))

## --- componenti di varianza stimate (sigma, non sd dei BLUP) -------------------------
VC.lmer <- t(sapply(fit.lmer, function(m) {
  v <- as.data.frame(VarCorr(m)); setNames(v$sdcor, v$grp)[c("Station", "Event", "Residual", "Source")] }))

## --- plot beta: pffr (blu, +-2 se) vs lmer (rosso, +-2 se per periodo) --------------
par(mfrow = c(3, 3))
for (j in 1:9) {
  Phi <- eval.basis(tg, beta.fd[[j]]$basis)
  b   <- eval.fd(tg, beta.fd[[j]])
  se  <- sqrt(rowSums((Phi %*% C.fd[[j]]) * Phi))
  lo  <- B.lmer[, j] - 2 * SE.lmer[, j]; up <- B.lmer[, j] + 2 * SE.lmer[, j]
  plot(tg, b, type = "n", ylim = range(b - 2 * se, b + 2 * se, lo, up),
       xlab = "log10 T", ylab = "", main = names(beta.fd)[j])
  polygon(c(tg, rev(tg)), c(b - 2 * se, rev(b + 2 * se)), col = adjustcolor("blue", .2), border = NA)
  lines(tg, b, col = "blue", lwd = 2)
  arrows(tt, lo, tt, up, length = 0, col = "red")
  points(tt, B.lmer[, j], pch = 16, cex = .6, col = "red")
}

SD.pffr <- cbind(dS2S = apply(S, 2, sd), dB = apply(E, 2, sd), dWS = apply(R, 2, sd), dL2L = apply(L, 2, sd))
SD.lmer <- t(sapply(fit.lmer, function(m) { re <- ranef(m)
c(dS2S = sd(re$Station[, 1]), dB = sd(re$Event[, 1]), dWS = sd(resid(m)), dL2L = sd(re$Source[, 1])) }))

par(mfrow = c(2, 2))
for (c in colnames(SD.pffr)) {
  matplot(T.emp, cbind(SD.pffr[, c], SD.lmer[, c]), type = "l", lty = 1, lwd = 2,
          col = c("blue", "red"), log = "x", xlab = "T [s]", ylab = "sd", main = c,
          ylim = c(0, max(SD.pffr, SD.lmer)))
}

legend("topright", c("Functional LMM", "Scalar LMM"), col = c("blue", "red"), lty = 1, bty = "n", y.intersp = 0.4)


## --- K-fold per evento (parallelo sui fold) ----------------------------------------
library(parallel)
K <- 10; set.seed(123)
ev   <- levels(df$Event)
fold <- setNames(sample(rep(1:K, length.out = length(ev))), ev)[as.character(df$Event)]
grp  <- c("Station", "Event", "Source")
sub.df <- function(d, i) lapply(d, function(v) if (is.matrix(v)) v[i, , drop = FALSE] else v[i])
lm.df  <- function(d, p) data.frame(y = d$Y[, p], sapply(d[names(Xl)[-1]], function(x) x[, p]), d[grp])

f0   <- as.formula(paste("Y ~", paste(names(Xl)[-1], collapse = " + ")))   # funzionale, solo effetti fissi
fs0  <- as.formula(paste("y ~", paste(names(Xl)[-1], collapse = " + ")))   # scalare, solo effetti fissi
mods <- c("scalar", "functional", "scalar.lmm", "functional.lmm")

cv.fold <- function(k) {
  tr <- which(fold != k); te <- which(fold == k)
  dtr <- sub.df(df, tr); dtr[grp] <- lapply(dtr[grp], droplevels)
  dte <- sub.df(df, te)
  for (g in c("Station", "Source")) dte[[g]] <- factor(dte[[g]], levels = levels(dtr[[g]]))
  stopifnot(!anyNA(dte$Station), !anyNA(dte$Source))           # livelli di test visti in training
  dte$Event <- factor(rep(levels(dtr$Event)[1], length(te)), levels = levels(dtr$Event))  # fittizio, termine escluso
  e <- array(NA, c(length(te), nT, length(mods)), dimnames = list(NULL, NULL, mods))
  
  ## funzionale LMM: fixed + Station + Source (Event escluso)
  m.f <- pffr(f, yind = yindex, data = dtr, algorithm = "bam", bs.yindex = bsy, bs.int = bsy, knots = kn)
  ex  <- grep("Event", sapply(m.f$smooth, `[[`, "label"), value = TRUE)
  e[, , "functional.lmm"] <- dte$Y - predict(m.f, newdata = dte[-1], exclude = ex)
  
  ## funzionale, solo effetti fissi
  m.f0 <- pffr(f0, yind = yindex, data = dtr, algorithm = "bam", bs.yindex = bsy, bs.int = bsy, knots = kn)
  e[, , "functional"] <- dte$Y - predict(m.f0, newdata = dte[-1])
  
  ## scalari: un fit per periodo
  for (p in 1:nT) {
    dp <- lm.df(dtr, p); dq <- lm.df(dte, p)
    e[, p, "scalar.lmm"] <- dq$y - predict(lmer(fs, data = dp), dq, re.form = ~ (1 | Station) + (1 | Source))
    e[, p, "scalar"]     <- dq$y - predict(lm(fs0, data = dp), dq)
  }
  list(te = te, e = e)
}

nc  <- min(K, detectCores() - 1)
out <- mclapply(1:K, cv.fold, mc.cores = nc)
err <- array(NA, c(n, nT, length(mods)), dimnames = list(NULL, NULL, mods))
for (o in out) err[o$te, , ] <- o$e

## --- MSE di test per periodo ----------------------------------------------------------
mse <- apply(err, c(2, 3), function(e) mean(e^2))
lab <- c("Scalar", "Functional", "Scalar LMM", "Functional LMM")
cl  <- c("red", "blue", "red", "blue"); lt <- c(2, 2, 1, 1)
graphics.off()
matplot(T.emp, mse, type = "l", lty = lt, lwd = 2, col = cl, log = "x",
        xlab = "T [s]", ylab = "MSE test")
legend("bottomleft", lab, col = cl, lty = lt, lwd = 2, bty = "n")
setNames(colMeans(mse), lab)
