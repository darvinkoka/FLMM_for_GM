rm(list=ls())
graphics.off()
cat("\014")

setwd("~/Documents/PhD/PROJECTS/FLMM_for_GM")
source("Code/utils_flmm.R")

library(readxl)
library(fda)
library(refund)
library(mgcv)
library(lme4)

#Environment preparation
ITA18 <- read.csv("Data/ITA18_source-path.csv")
Params <- readRDS("Data/Params.RData")
T.emp <- Params$Periods
Mh <- Params$Mh
h <- Params$h
Mref <- Params$Mh
xlist <- build_xlist_fda(ITA18, T.emp, Mh, Mref, h)
logsa.names <- names(ITA18)[grep("^log",names(ITA18))]
logsa <- ITA18[,logsa.names]
bsplines <- smooth_fd(t(as.matrix(logsa)), TRUE, grid = 10^-seq(-1, 12, by = 1), tp = log10(T.emp))
logsa.fd <- bsplines[[2]]

## SCALAR LMM

n  <- nrow(logsa); tt <- log10(T.emp)
Xl <- lapply(xlist, function(x) if (is.fd(x)) t(eval.fd(tt, x)) else matrix(x, n, length(tt)))
names(Xl) <- paste0("x", seq_along(Xl))
nT <- length(T.emp)
fs <- as.formula(paste("y ~", paste(names(Xl)[-1], collapse = " + "),
                      "+ (1|Station) + (1|Event) + (1|Source)"))

df <- c(list(Y       = as.matrix(logsa),
             Station = as.factor(ITA18$Station),
             Event   = as.factor(ITA18$event),
             Source  = as.factor(ITA18$source_term)),
        Xl[-1])     

#per k = 7 (T = 0.15) c'è singolarità
fit.lmer <- lapply(1:nT, function(k) {
  cat("k = ",k,"\n")
  dk <- data.frame(y = df$Y[, k], sapply(Xl[-1], function(x) x[, k]),
                   Station = df$Station, Event = df$Event, Source = df$Source)
  lmer(fs, data = dk, REML = TRUE)
})

B.lmer  <- t(sapply(fit.lmer, fixef))
SE.lmer <- t(sapply(fit.lmer, function(m) sqrt(diag(vcov(m)))))

VC.lmer <- t(sapply(fit.lmer, function(m) {
  v <- as.data.frame(VarCorr(m)); setNames(v$sdcor, v$grp)[c("Station", "Event", "Residual", "Source")] }))

par(mfrow = c(3, 3))
for (j in 1:9) {
  lo  <- B.lmer[, j] - 2 * SE.lmer[, j]; up <- B.lmer[, j] + 2 * SE.lmer[, j]
  plot(tt, B.lmer[, j], pch = 16, cex = .6, col = "red", ylim = c(min(lo), max(up)))
  arrows(tt, lo, tt, up, length = 0, col = "red")
  grid()
}


SD.lmer <- t(sapply(fit.lmer, function(m) { re <- ranef(m)
c(dS2S = sd(re$Station[, 1]), dB = sd(re$Event[, 1]), dWS = sd(resid(m)), dL2L = sd(re$Source[, 1])) 
}))

par(mfrow = c(2, 2))
for (c in colnames(SD.lmer)) {
  matplot(T.emp, SD.lmer[, c], pch = 19, lty = 1, lwd = 2,
          col = "red", log = "x", xlab = "T [s]", ylab = "sd",
          cex.lab = 1.5, cex.axis = 1.5)
  grid()
}

## FUNCTIONL LMM - SPLINE BASIS REPRESENTATION

n  <- nrow(logsa)
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



basi <- list(bs = "bs", k = 20, m = c(3, 2))
fit.pffr <- pffr(f, yind = tt, data = df, algorithm = "bam", bs.int = basi)
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
  if (j == 1) { X <- cbind(1, X); idx <- c(1, idx) }
  M   <- qr.solve(eval.basis(x, bas), X)
  beta.fd[[j]] <- fd(M %*% th[idx], bas)
  C.fd[[j]]    <- M %*% Vp[idx, idx] %*% t(M)
}
names(beta.fd) <- names(C.fd) <- paste0("beta", 1:9)


tg <- seq(rng[1], rng[2], length.out = 200)
par(mfrow = c(3, 3))
for (j in 1:9) {
  Phi <- eval.basis(tg, beta.fd[[j]]$basis)
  b   <- eval.fd(tg, beta.fd[[j]])
  se  <- sqrt(rowSums((Phi %*% C.fd[[j]]) * Phi))
  lo  <- b - 2 * se; up <- b + 2 * se
  plot(tg, b, col = "blue", lwd = 2, ylim = c(min(lo), max(up)))
  polygon(c(tg, rev(tg)), c(lo, rev(up)), col = adjustcolor("blue", .2), border = NA)
  grid()
}

tr <- predict(fit.pffr, type = "terms"); nm <- names(tr)
S  <- tr[[grep("Station", nm)]][!duplicated(df$Station), ]
E  <- tr[[grep("Event",   nm)]][!duplicated(df$Event), ]
L  <- tr[[grep("Source",  nm)]][!duplicated(df$Source), ]
R  <- residuals(fit.pffr)

par(mfrow = c(2, 2))
for (i in 1:4) {
  s <- apply(list(S, E, L, R)[[i]], 2, sd)
  plot(T.emp, s, type = "l", lty = 1, log = "x", xlab = "T [s]", ylab = "sd", cex.lab = 1.5, cex.axis = 1.5)
  grid()
}

par(mfrow = c(2, 2))
for (i in 1:4) {
  s <- apply(list(S, E, L, R)[[i]], 2, sd)
  plot(T.emp, s, type = "l", lty = 1, log = "x", xlab = "T [s]", ylab = "sd",
       ylim = range(s, SD.lmer[, i]), cex.lab = 1.5, cex.axis = 1.5)
  points(T.emp, SD.lmer[, i], pch = 19, col = "red")
  grid()
}


## FUNCTIONL LMM - FPC BASIS REPRESENTATION
basi <- list(bs = "bs", k = 20, m = c(3, 2))
f0   <- as.formula(paste("Y ~", paste(names(Xl)[-1], collapse = " + ")))
fit0 <- pffr(f0, yind = tt, data = df, algorithm = "bam", bs.yindex = basi, bs.int = basi)
R0   <- residuals(fit0)

# FPC di Station, Event, Source
grp <- c("Station", "Event", "Source")
fpc <- fpc_mom(R0, df[grp], pve = 0.95)
sapply(fpc, function(x) length(x$evalues))

# pcre vuole efunctions/evalues come variabili semplici
for (g in grp) {
  assign(paste0("Phi.", g), fpc[[g]]$efunctions)
  assign(paste0("nu.", g),  fpc[[g]]$evalues)
}

re.pc <- paste0("pcre(id = ", grp, ", efunctions = Phi.", grp,
                ", evalues = nu.", grp, ", yind = tt)")
f.pc     <- as.formula(paste("Y ~", paste(c(names(Xl)[-1], re.pc), collapse = " + ")))
fit.pcre <- pffr(f.pc, yind = tt, data = df, algorithm = "bam",  bs.yindex = basi, bs.int = basi)



sm <- fit.pcre$smooth; th <- fit.pcre$coefficients; Vp <- fit.pcre$Vp
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
  if (j == 1) { X <- cbind(1, X); idx <- c(1, idx) }
  M   <- qr.solve(eval.basis(x, bas), X)
  beta.fd[[j]] <- fd(M %*% th[idx], bas)
  C.fd[[j]]    <- M %*% Vp[idx, idx] %*% t(M)
}
names(beta.fd) <- names(C.fd) <- paste0("beta", 1:9)


tg <- seq(rng[1], rng[2], length.out = 200)
par(mfrow = c(3, 3))
for (j in 1:9) {
  Phi <- eval.basis(tg, beta.fd[[j]]$basis)
  b   <- eval.fd(tg, beta.fd[[j]])
  se  <- sqrt(rowSums((Phi %*% C.fd[[j]]) * Phi))
  lo  <- b - 2 * se; up <- b + 2 * se
  plot(tg, b, col = "blue", type = "l",lwd = 1, ylim = c(min(lo), max(up)))
  polygon(c(tg, rev(tg)), c(lo, rev(up)), col = adjustcolor("blue", .2), border = NA)
  points(tt, B.lmer[, j], pch = 16, cex = .6, col = "red", ylim = c(min(lo), max(up)))
  lo  <- B.lmer[, j] - 2 * SE.lmer[, j]; up <- B.lmer[, j] + 2 * SE.lmer[, j]
  arrows(tt, lo, tt, up, length = 0, col = "red")
  grid()
}
