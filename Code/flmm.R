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
library(denseFLMM)
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
  plot(T.emp, B.lmer[, j], pch = 16, cex = .6, col = "red", ylim = c(min(lo), max(up)), log = "x")
  arrows(T.emp, lo, T.emp, up, length = 0, col = "red")
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
# 
# n  <- nrow(logsa)
# Xl <- lapply(xlist, function(x) if (is.fd(x)) t(eval.fd(tt, x)) else matrix(x, n, length(tt)))
# names(Xl) <- paste0("x", seq_along(Xl))
# 
# df <- c(list(Y       = as.matrix(logsa),
#              Station = as.factor(ITA18$Station),
#              Event   = as.factor(ITA18$event),
#              Source  = as.factor(ITA18$source_term)),
#         Xl[-1])                                            # x1 = intercetta, gestita da pffr
# 
# f <- as.formula(paste("Y ~", paste(c(names(Xl)[-1],
#                                      's(Station, bs="re")',
#                                      's(Event, bs="re")',
#                                      's(Source, bs="re")'), collapse = " + ")))
# 
# 
# 
# basi <- list(bs = "bs", k = 20, m = c(3, 2))
# fit.pffr <- pffr(f, yind = tt, data = df, algorithm = "bam", bs.int = basi)
# sm <- fit.pffr$smooth; th <- fit.pffr$coefficients; Vp <- fit.pffr$Vp
# rng <- range(tt)
# beta.fd <- C.fd <- vector("list", 9)
# 
# for (j in 1:9) {
#   s   <- sm[[j]]
#   br  <- c(rng[1], s$knots[s$knots > rng[1] & s$knots < rng[2]], rng[2])
#   bas <- create.bspline.basis(rng, norder = s$m[1] + 1, breaks = br)
#   x   <- seq(rng[1], rng[2], length.out = 10 * bas$nbasis)
#   nd  <- setNames(data.frame(x), s$term[1]); if (s$by != "NA") nd[[s$by]] <- 1
#   X   <- PredictMat(s, nd)
#   idx <- s$first.para:s$last.para
#   if (j == 1) { X <- cbind(1, X); idx <- c(1, idx) }
#   M   <- qr.solve(eval.basis(x, bas), X)
#   beta.fd[[j]] <- fd(M %*% th[idx], bas)
#   C.fd[[j]]    <- M %*% Vp[idx, idx] %*% t(M)
# }
# names(beta.fd) <- names(C.fd) <- paste0("beta", 1:9)
# 
# 
# tg <- seq(rng[1], rng[2], length.out = 200)
# par(mfrow = c(3, 3))
# for (j in 1:9) {
#   Phi <- eval.basis(tg, beta.fd[[j]]$basis)
#   b   <- eval.fd(tg, beta.fd[[j]])
#   se  <- sqrt(rowSums((Phi %*% C.fd[[j]]) * Phi))
#   lo  <- b - 2 * se; up <- b + 2 * se
#   plot(tg, b, col = "blue", lwd = 2, ylim = c(min(lo), max(up)))
#   polygon(c(tg, rev(tg)), c(lo, rev(up)), col = adjustcolor("blue", .2), border = NA)
#   grid()
# }
# 
# tr <- predict(fit.pffr, type = "terms"); nm <- names(tr)
# S  <- tr[[grep("Station", nm)]][!duplicated(df$Station), ]
# E  <- tr[[grep("Event",   nm)]][!duplicated(df$Event), ]
# L  <- tr[[grep("Source",  nm)]][!duplicated(df$Source), ]
# R  <- residuals(fit.pffr)
# 
# par(mfrow = c(2, 2))
# for (i in 1:4) {
#   s <- apply(list(S, E, L, R)[[i]], 2, sd)
#   plot(T.emp, s, type = "l", lty = 1, log = "x", xlab = "T [s]", ylab = "sd", cex.lab = 1.5, cex.axis = 1.5)
#   grid()
# }
# 
# par(mfrow = c(2, 2))
# for (i in 1:4) {
#   s <- apply(list(S, E, L, R)[[i]], 2, sd)
#   plot(T.emp, s, type = "l", lty = 1, log = "x", xlab = "T [s]", ylab = "sd",
#        ylim = range(s, SD.lmer[, i]), cex.lab = 1.5, cex.axis = 1.5)
#   points(T.emp, SD.lmer[, i], pch = 19, col = "red")
#   grid()
# }


## FUNCTIONL LMM - FPC BASIS REPRESENTATION
library(denseFLMM)
basi <- list(bs = "bs", k = length(T.emp)+2, m = c(3, 2))
tt   <- log10(T.emp)
tg   <- seq(min(tt), max(tt), length.out = 100)

# STEP 1: residui in working independence
y.centered.fd <- readRDS("Data/Responce_centered.RData")
Yc <- t(eval.fd(tg, y.centered.fd))

# STEP 2-3: covarianze + FPCA
Effects <- c("Station", "Event", "Source")
groups  <- sapply(df[Effects], as.integer)
Zv      <- rep(list(matrix(1, n, 1)), length(Effects))

fitA <- denseFLMM(Y = Yc, gridpoints = tg, Zvars = Zv, groups = groups, NPC = c(10,10,2,10), smooth = TRUE, smoothalg = "bamREML")
par(mfrow = c(1, 4))
for (g in seq_along(fitA$nu)) {
  plot(cumsum(fitA$nu[[g]]) / sum(fitA$nu[[g]]), type = "b", pch = 19, ylim = c(0, 1), xlab = "PC", ylab = "PVE", main = c(Effects, "E")[g])
  abline(h = .95, lty = 2); grid()
}

fit0    <- denseFLMM(Y = Yc, gridpoints = tg, Zvars = Zv, groups = groups, NPC = c(2,2,2,3), smooth = TRUE, smoothalg = "bamREML")
fit0$exvar
# STEP 4: refit FAMM (beta ri-stimati)

df$res <- factor(1:n)
re  <- c(Effects, "res")
PHI <- lapply(fit0$phi, function(P) apply(P, 2, function(p) approx(tg, p, xout = tt)$y))
df[paste0("PHI.", re)] <- PHI
pc  <- sprintf("pcre(id = %s, efunctions = PHI.%s, evalues = fit0$nu[[%d]], yind = tt)", re, re, 1:4)
f.fpc    <- as.formula(paste("Y ~", paste(c(names(Xl)[-1], pc), collapse = " + ")))
fit.famm <- pffr(f.fpc, yind = tt, data = df, algorithm = "bam",method = "fREML", nthreads = 4, bs.int = basi, bs.yindex = basi, control = gam.control(trace = TRUE))
saveRDS(fit.famm, "Results/famm.RData")
#PLOTS
cf <- coef(fit.famm, n1 = 100)$smterms[1:9]; b0 <- fit.famm$coefficients["(Intercept)"]
par(mfrow = c(3, 3))
for (j in 1:9) {
  Tg <- 10^cf[[j]]$x; b <- cf[[j]]$coef$value + (j == 1) * b0; se <- cf[[j]]$coef$se
  lo <- b - 2 * se; up <- b + 2 * se
  plot(Tg, b, type = "l", lwd = 2, col = "blue", log = "x", ylim = range(lo, up, B.lmer[, j]), xlab = "T [s]", ylab = names(cf)[j])
  polygon(c(Tg, rev(Tg)), c(lo, rev(up)), col = adjustcolor("blue", .2), border = NA)
  points(T.emp, B.lmer[, j], pch = 16, cex = .6, col = "red")
  grid()
}




tr <- predict(fit.famm, type = "terms"); nm <- names(tr)
RE <- c(dS2S = "Station", dB = "Event", dL2L = "Source")
par(mfrow = c(1, 3))
for (r in names(RE)) {
  s <- apply(tr[[grep(RE[r], nm)]][!duplicated(df[[RE[r]]]), ], 2, sd)
  plot(T.emp, s, type = "l", lwd = 2, col = "blue", log = "x", ylim = range(0, s, SD.lmer[, r]), xlab = "T [s]", ylab = r)
  points(T.emp, SD.lmer[, r], pch = 19, col = "red")
  grid()
}
graphics.off()

sdre <- function(x) apply(tr[[grep(x, nm)]][!duplicated(df[[x]]), ], 2, sd)
plot(T.emp, sqrt(sdre("Event")^2 + sdre("Source")^2 + sdre("Station")^2), type = "l", lwd = 2, col = "blue", log = "x", xlab = "T [s]", ylab = "dB + dL2L", ylim = c(0.2,0.6))
points(T.emp, sqrt(SD.lmer[, "dB"]^2 + SD.lmer[, "dL2L"]^2 + SD.lmer[,"dS2S"]^2), pch = 19, col = "red"); grid()


vc <- gam.vcomp(fit.famm); if (is.list(vc)) vc <- vc$vc; if (is.matrix(vc)) vc <- vc[, 1]
SD.famm <- sapply(seq_along(k), function(i) vc[grep(Effects[k[i]], names(vc))] * sqrt(PHI[[i]]^2 %*% fit0$nu[[k[i]]]))
colnames(SD.famm) <- Effects[k]
SD.famm <- cbind(SD.famm, Residual = sqrt(fit.famm$sig2))

lab <- c(Station = expression(phi[S2S]), Event = expression(tau), Source = expression(tau[L2L]), Residual = expression(phi[0]))
tot <- function(V) sqrt(rowSums(V^2))
par(mfrow = c(2, 3))
for (g in colnames(SD.famm)) {
  plot(T.emp, SD.famm[, g], type = "l", lwd = 2, col = "blue", log = "x", ylim = range(0, SD.famm[, g], VC.lmer[, g]), xlab = "T [s]", ylab = lab[g])
  points(T.emp, VC.lmer[, g], pch = 19, col = "red"); grid()
}
plot(T.emp, tot(SD.famm), type = "l", lwd = 2, col = "blue", log = "x", ylim = range(0, tot(SD.famm), tot(VC.lmer)), xlab = "T [s]", ylab = expression(sigma))
points(T.emp, tot(VC.lmer), pch = 19, col = "red"); grid()

