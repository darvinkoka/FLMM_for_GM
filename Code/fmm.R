rm(list=ls())
graphics.off()
cat("\014")
gc()

setwd('~/Documents/PhD/PROJECTS/PRIN')
source("codice/utils_data_integration.R", echo = TRUE)
df_emp<- read.csv("data/ITA18_SA_flatfile.csv", sep = ";")
library(readxl)
event_param<-read_excel("data//BB_SPEEDset_v2.3_flatfile.xlsx", sheet = "Event_source_parameters")
event_proc <- read_excel("data//BB_SPEEDset_v2.3_flatfile.xlsx", sheet = "Event_processing_parameters")
load("data/ITA18_parameters.RData")

#DATA PREPROCESSING
SA.names.emp<-names(df_emp)[grep("rotD50_T",names(df_emp))][-1]
T.emp <- sub(".*_T", "", SA.names.emp)
T.emp <- as.numeric(sub("_",".",T.emp))
Mh <- ITA18.parameters$Mh.vec[-1]
h <- ITA18.parameters$h.vec[-1]
Mref <- ITA18.parameters$Mref.vec[-1]
ITA18 <- data.frame(
  event = df_emp$event_id,
  mw = df_emp$Mw,
  djb = df_emp$JB_dist,
  vs30 = df_emp$vs30_m_sec,
  sof = df_emp$fm_type_code
) 
ITA18$djb[is.na(ITA18$djb)] <- df_emp$epi_dist[is.na(ITA18$djb)]
ITA18$vs30[is.na(ITA18$vs30)]<-df_emp$vs30_m_sec_WA[is.na(ITA18$vs30)]
ITA18$djb[is.na(ITA18$djb)]   <- df_emp$epi_dist[is.na(ITA18$djb)]
ITA18$vs30[is.na(ITA18$vs30)] <- df_emp$vs30_m_sec_WA[is.na(ITA18$vs30)]

for(i in 1:length(T.emp)){
  ITA18 <- add_regressors(ITA18, h[i], Mh[i], Mref[i], T.emp[i])
}
SA_data <- df_emp[, SA.names.emp]
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
ITA18 <- ITA18[!na_rows,]
logsa <- logsa[!na_rows,]
#BUILD FUNCTIONAL DATASET 
xlist <- build_xlist_fda(ITA18, T.emp, Mh, Mref, h)

logsa <- smooth_fd(t(as.matrix(logsa)), TRUE, grid = 10^-seq(-1, 12, by = 1), tp = log10(T.emp))

#ESTIMATE FIXED PARAMETERS
logsa.names
