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
