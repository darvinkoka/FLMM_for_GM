################################################################################
# Reasenberg Clustering Algorithm
################################################################################


# Initial Setup ----------------------------------------------------------------

library(dplyr)
library(sf)
library(oce)

setwd("C:/Users/stefy/OneDrive - Politecnico di Milano/Desktop/Dottorato/Corsi PhD/ASLCD/FLMM_for_GM")
source("Code/fmm.R")


# Preprocessing ----------------------------------------------------------------

View(ITA18)
View(df_emp)

events <- !duplicated(df_emp$event_id)
df_cluster <- data.frame(
  event_id  = df_emp$event_id[events],
  time      = df_emp$event_time[events],
  lat       = df_emp$ev_latitude[events],
  lon       = df_emp$ev_longitude[events],
  depth     = df_emp$ev_depth_km[events],
  Mw        = df_emp$Mw[events],
  nation    = df_emp$ev_nation_code[events],
  row.names = NULL
)
View(df_cluster)
nrow(df_cluster)
anyDuplicated(df_cluster$event_id)
df_cluster$pref <- sub("-.*", "", df_cluster$event_id)
table(df_cluster$pref, df_cluster$nation)

# keep only Italy
df_ita <- df_cluster[grepl("^(IT|EMSC|SI)", df_cluster$event_id),]
View(df_ita)

# fix time-type and sort
str(df_ita)
df_ita$time <- as.POSIXct(df_ita$time, tz = "UTC")
df_ita <- df_ita[order(df_ita$time), ]


# Useful functions -------------------------------------------------------------

# 2D dist func (Reasenberg paper)
dist_2D <- function(j, i, df) {
  geosphere::distHaversine(c(df$lon[i], df$lat[i]), c(df$lon[j], df$lat[j])) / 1000
}

# 3D dist func (extension)
dist_3D <- function(j, i, df) {
  d_epi <- geosphere::distHaversine(c(df$lon[i], df$lat[i]), c(df$lon[j], df$lat[j])) / 1000
  d_depth <- df$depth[j] - df$depth[i]
  sqrt(d_epi^2 + d_depth^2)
}

# dist threshold ri
dist_ri <- function(M) {
  return (10^(0.4 * M - 1.943))
}

# compute delta magnitude
compute_dM <- function(M, Mfix, xfix) {
  return ((1-xfix) * M - Mfix)
}

# compute D parameter
compute_D <- function(deltaM) {
  10^(((max(deltaM,0) - 1) * 2) / 3)
}

# compute tau
compute_tau <- function(prob, deltat, denom) {
  return(- log(1 - prob) * deltat / denom)
}


# Algorithm --------------------------------------------------------------------

# Initial parameters (default version)
tau_min = 1
tau_max = 10
p1 = 0.95
xk = 0.5
Meff = 1.5
rfact = 10

# Reasenberg Algo
reasenberg_cluster <- function(df, rfact = 10, Meff = 1.5, xk = 0.5, tau_min = 1, tau_max = 10, p1 = 0.95){
  
  # Initialization
  N = nrow(df)
  clusters <- rep(NA, N)
  mainshocks <- rep(NA, N)
  background <- rep(FALSE, N)
  reas_background <- rep(FALSE, N)
  n_clusters <- 0
  
  # Loop
  for (i in 1:(N-1)) {
    
    # Compute time window tau_i
    # If event i is not yet in a cluster
    if (is.na(clusters[i])) {
      tau_i <- tau_min
    }
    
    # If event i is already in a cluster
    else {
      cid <- clusters[i]
      idx_main <- mainshocks[cid]
      Mm <- df$Mw[idx_main]
      
      # If event i is the largest in the cluster
      if(df$Mw[i] > Mm) {
        idx_main <- i
        mainshocks[cid] <- i
        Mm <- df$Mw[i]
        tm <- df$time[i]
        tau_i <- tau_min
      }
      
      # If event i is not the largest
      else {
        tm <- df$time[idx_main]
        dt <- as.numeric(difftime(df$time[i], tm, units = "days"))
        dM <- compute_dM(Mm, Meff, xk)
        D <- compute_D(dM)
        tau_i <- compute_tau(p1, dt, D)
        tau_i <- min(max(tau_i, tau_min), tau_max)
      }
    }
    
    # Now consider subsequent events
    j <- i+ 1
    
    while (j <= N) {
      
      dt_ji <- as.numeric(difftime(df$time[j], df$time[i], units = "days")) 
      
      # 0) Event j is too far in time wrt event i --> event j is no longer a clustering candidate
      if (dt_ji >= tau_i)
        break # exit the for loop
      
      # 1) If event i and event j are already in the same cluster --> C = Ci = Cj
      if(!is.na(clusters[i]) && !is.na(clusters[j]) && clusters[i] == clusters[j]){
        j <- j + 1
        next # go to next event j
      }
      
      # Check dist from the last event (event i)
      r_i <- dist_ri(df$Mw[i])
      # d_last <- dist_2D(j, i, df)
      d_last <- dist_3D(j, i, df)
      cond1 <- d_last < rfact * r_i
      
      # Check dist from the main event (event main)
      cond2 <- FALSE
      if (!cond1 && tau_i > tau_min && !is.na(clusters[i])) {
        idx_main <- mainshocks[clusters[i]]
        Mm <- df$Mw[idx_main]
        r_m <- dist_ri(Mm)
        # d_main <- dist_2D(j, idx_main, df)
        d_main <- dist_3D(j, idx_main, df)
        cond2 <- d_main < r_m 
      }
      
      # If the dist-condition is satisfied
      if (cond1 || cond2) {
        
        # 2) If event i belongs to cluster i and j belongs to cluster j --> C = Ci U Cj
        if(!is.na(clusters[i]) && !is.na(clusters[j]) && clusters[i] != clusters[j]){
          
          cid_1 = clusters[i]
          cid_2 = clusters[j]
          
          # # Keep the earliest cluster
          # if(df$time[i] > df$time[j]) {
          #   tmp <- cid_1
          #   cid_1 <- cid_2
          #   cid_2 <- tmp
          # }
          
          clusters[clusters == cid_2] <- cid_1
          idx1 <- mainshocks[cid_1]
          idx2 <- mainshocks[cid_2]
          
          if(df$Mw[idx2] > df$Mw[idx1]) {
            mainshocks[cid_1] <- idx2
          }
          
          mainshocks[cid_2] <- NA 
        }
        
        # 3) If event i belongs to a cluster, while j doesn't --> C = Ci U {j}
        else if(!is.na(clusters[i]) && is.na(clusters[j])) {
          
          cid <- clusters[i]
          clusters[j] <- cid
          idx_main <- mainshocks[cid]
          
          if(df$Mw[j] > df$Mw[idx_main]) {
            mainshocks[cid] <- j
          }
        }
        
        # 4) If event j belongs to a cluster, while i doesn't --> C = Cj U {i}
        else if(!is.na(clusters[j]) && is.na(clusters[i])) {
          
          cid <- clusters[j]
          clusters[i] <- cid
          idx_main <- mainshocks[cid]
          
          if(df$Mw[i] > df$Mw[idx_main]) {
            mainshocks[cid] <- i
          }
        }
        
        # 5) If neither event i nor event j belong to a cluster --> C = {i,j}
        else if(is.na(clusters[i]) && is.na(clusters[j])) {
          
          n_clusters <- n_clusters + 1
          cid <- n_clusters
          clusters[i] <- cid
          clusters[j] <- cid
          
          if(df$Mw[i] > df$Mw[j]) {
            mainshocks[cid] <- i
          }
          else {
            mainshocks[cid] <- j
          }
        }
      }
      
      j <- j + 1
    }
  }
  
  # Final clusters
  df$cluster <- clusters
  
  # Final mainshocks
  mainshocks_final <- rep(NA, max(clusters, na.rm = TRUE))
  for (cid in na.omit(unique(clusters))) {
    idx <- which(clusters == cid)
    mainshocks_final[cid] <- idx[which.max(df$Mw[idx])]
  }
  df$main_final <- seq_len(nrow(df)) %in% mainshocks_final[!is.na(mainshocks_final)]
  
  # Events not in a cluster
  df$background <- is.na(clusters)
  
  # Reasenberg background events: mainshocks for each cluster + events not in a cluster
  df$reas_background <- df$background | df$main_final
  
  # Event !background or !mainshock
  df$is_aftershock <- !df$background & !df$main_final
  
  return(df)
}


# Application ------------------------------------------------------------------

reas_ita <- reasenberg_cluster(df_ita)

# Ita
View(reas_ita)
max(reas_ita$cluster, na.rm = T)

# Save Reasenberg clusters
save(ITA18, df_emp, reas_ita, file="Results/Clustering/clustering_data.RData")

