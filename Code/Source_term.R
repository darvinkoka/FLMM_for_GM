################################################################################
# SOURCE - Corrective Term
################################################################################

# Initial Setup ----------------------------------------------------------------

library(sf)
library(rnaturalearth)
library(ggplot2)
library(stringr)
library(dplyr)
library(lubridate)

setwd("C:/Users/stefy/OneDrive - Politecnico di Milano/Desktop/Dottorato/Corsi PhD/ASLCD/FLMM_for_GM")
load("C:/Users/stefy/OneDrive - Politecnico di Milano/Desktop/Dottorato/Corsi PhD/ASLCD/FLMM_for_GM/Results/Clustering/clustering_data.RData")
source("Code/DISS_SourceAreas.R")

# Preprocessing ----------------------------------------------------------------

# Italy
reas_ita$cluster_final <- reas_ita$cluster
reas_ita_sf <- st_as_sf(reas_ita, coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# Plot (Reasenberg clusters) - Zoom Italy
# Geographic window
world <- ne_countries(scale = "medium", returnclass = "sf")
italy_events <- reas_ita_sf[
  st_coordinates(reas_ita_sf)[, 1] >= 6 &
  st_coordinates(reas_ita_sf)[, 1] <= 19 &
  st_coordinates(reas_ita_sf)[, 2] >= 35 &
  st_coordinates(reas_ita_sf)[, 2] <= 47.5,
]
ggplot() +
  geom_sf(
    data = world,
    fill = "white",
    color = "grey80"
  ) +
  # Background events
  geom_sf(
    data = subset(italy_events, is.na(cluster_final)),
    color = "grey70",
    size = 0.8,
    alpha = 0.6
  ) +
  # Clustered events
  geom_sf(
    data = subset(italy_events, !is.na(cluster_final)),
    aes(color = factor(cluster_final)),
    size = 1.8,
    alpha = 0.9
  ) +
  scale_color_viridis_d(option = "turbo") +
  coord_sf(
    xlim = c(6, 19),
    ylim = c(35, 47.5),
    expand = FALSE
  ) +
  labs(
    color = "Cluster",
    title = "Reasenberg clusters (Italy)"
  ) +
  theme(
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.3, "cm")
  )
# ggsave("Plots/Reasenberg_Italy.png", dpi = 300, width = 10, height = 7)


# Merge with diss dataset - Italy
st_crs(reas_ita_sf) == st_crs(diss) # True
diss <- st_make_valid(diss)
diss <- st_transform(diss, st_crs(reas_ita_sf))
reas_ita_diss <- st_join(reas_ita_sf, diss)
nrow(reas_ita_sf) # 117
nrow(reas_ita_diss) # 122 --> The DISS dataset has 5 additional events --> duplicated events --> some events lies in the intersection of multiple DISS areas


# DISS clusters ----------------------------------------------------------------

# Function to merge Reasenberg clusters that share at least a DISS source

merge_shared_diss_clusters <- function(x) {
  x_df <- if (inherits(x, "sf")) {st_drop_geometry(x)} else {as.data.frame(x)}
  
  required_cols <- c("event_id", "cluster_final", "IDSource")
  missing_cols <- setdiff(required_cols, names(x_df))
  if (length(missing_cols) > 0) {
    stop("Missing columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Reasenberg clusters
  cluster_ids <- unique(x_df$cluster_final[!is.na(x_df$cluster_final)])
  if (length(cluster_ids) == 0) {
    x_out <- x
    x_out$new_cluster <- NA
    event_clusters <- x_df %>%
      distinct(event_id) %>%
      mutate(cluster_final = NA,
             new_cluster = NA
      )
    return(list(joined = x_out, 
                cluster_map = data.frame(cluster_final = cluster_ids, new_cluster = integer(0)),
                event_clusters = event_clusters))
  }
  
  # DISS sources for each R. cluster
  cluster_sources <- lapply(cluster_ids, function(cid) {
    idx <- !is.na(x_df$cluster_final) & x_df$cluster_final == cid
    unique(as.character(x_df$IDSource[idx][!is.na(x_df$IDSource[idx])]))
  })
  
  # Merge of R. clusters with at least 1 common DISS source
  n_clusters <- length(cluster_ids)
  labels <- seq_len(n_clusters)
  changed <- TRUE
  while (changed) {
    old_labels <- labels
    if (n_clusters >= 2) {
      for (i in seq_len(n_clusters - 1)) {
        for (j in seq.int(i + 1, n_clusters)) {
          common <- intersect(cluster_sources[[i]], cluster_sources[[j]])
          if (length(common) > 0) {
            old_ids <- unique(c(labels[i], labels[j]))
            new_id <- min(old_ids)
            labels[labels %in% old_ids] <- new_id
          }
        }
      }
    }
    changed <- !identical(labels, old_labels)
  }
  
  # New clusters
  cluster_map <- data.frame(cluster_final = cluster_ids,
                            new_cluster = as.integer(match(labels, unique(labels))))
  x_df$new_cluster <- NA
  idx_cluster <- !is.na(x_df$cluster_final)
  x_df$new_cluster[idx_cluster] <- cluster_map$new_cluster[match(x_df$cluster_final[idx_cluster], cluster_map$cluster_final)]
  if (inherits(x, "sf")) {
    x_out <- x
    x_out$new_cluster <- x_df$new_cluster
  } else {
    x_out <- x_df
  }
  
  event_clusters <- x_df %>% distinct(event_id, cluster_final, new_cluster)
  conflicting_events <- event_clusters %>%
    count(event_id, name = "n_cluster") %>%
    filter(n_cluster > 1)
  if (nrow(conflicting_events) > 0) {print(conflicting_events)
    stop("Some events belong to different cluster_final.")
  }
  stopifnot(!anyDuplicated(event_clusters$event_id))
  list(
    joined = x_out,
    cluster_map = cluster_map,
    event_clusters = event_clusters
  )
}


# Italy
res_ita <- merge_shared_diss_clusters(reas_ita_diss)
ita_diss <- res_ita$joined
clusters_ita <- res_ita$cluster_map
ita_event_clusters <- res_ita$event_clusters


# Stress Drop clusters ---------------------------------------------------------

# Stress Drop construction (Bindi et al. (2020))
# dsigma_arias = 1/N * sum_j(Y_ij * P(R_ij))

# Fixed parameters
g <- 9.81
beta <- 3.2   
rho <- 2800 
Rtheta <- 0.55
fmax <- 25


# Italy
ita_emp <- df_emp[grepl("^(IT|EMSC|SI)", df_emp$event_id),]

# Compute Yij
ita_emp$M0 <- 10^(1.5*ita_emp$Mw + 16.05) # Seismic moment
C1 <- (2*g*beta/(pi*fmax))^(3/5)
C2 <- ((1e6*rho)/(2*Rtheta*(2*pi)^2))^(6/5)
ita_emp$Yij <- C1 * C2 * (1/(8.47*ita_emp$M0))^(1/5) * (ita_emp$rotD50_ia/ita_emp$rotD50_T90)^(3/5)


### Approach 1: Bindi et al. (2020) ----
sum(is.na(ita_emp$rup_dist)) #3674
ita_emp$Rij <- ita_emp$rup_dist
ita_emp$Rij[is.na(ita_emp$Rij)] <- ita_emp$JB_dist[is.na(ita_emp$Rij)] # Follow Darvin's procedure
ita_emp$Rij[is.na(ita_emp$Rij)] <- ita_emp$epi_dist[is.na(ita_emp$Rij)]
sum(is.na(ita_emp$Rij)) # Ok
ita_emp$distance_type <- ifelse(!is.na(ita_emp$rup_dist), "rupture", ifelse(!is.na(ita_emp$JB_dist), "Joyner-Boore", ifelse(!is.na(ita_emp$epi_dist), "epicentral", NA)))
ita_emp$Rmax <- ifelse(!is.na(ita_emp$Mw) & ita_emp$Mw <= 4.5, 20, !is.na(ita_emp$Mw) & ifelse(ita_emp$Mw >= 7, 140, 20 + (ita_emp$Mw - 4.5) * (140 - 20) / (7 - 4.5)))
ita_emp$geom_corr <- ita_emp$Yij * ita_emp$Rij^(6/5)
ita_emp$keep <- (!is.na(ita_emp$Mw) & ita_emp$Mw >= 4 & !is.na(ita_emp$Rij) & !is.na(ita_emp$Rmax) & ita_emp$Rij <= ita_emp$Rmax & is.finite(ita_emp$geom_corr) & ita_emp$geom_corr > 0)
ita_emp_selected <- ita_emp[ita_emp$keep, ]

ita_stressdrop <- aggregate(geom_corr ~ event_id, data = ita_emp_selected, FUN = mean)
names(ita_stressdrop)[names(ita_stressdrop) == "geom_corr"] <-"stress_drop"
event_info <- aggregate(cbind(Mw) ~ event_id, data = ita_emp_selected, FUN = function(z) z[1])
names(event_info)[names(event_info) == "Mw"] <- "Mw"
n_stations <- aggregate(event_id ~ event_id,data = ita_emp_selected,FUN = length)
names(n_stations)[names(n_stations) == "event_id.1"] <-"n_stations"
ita_stressdrop <- merge(ita_stressdrop, event_info, by = "event_id", all.x = TRUE)
counts <- as.data.frame(table(ita_emp_selected$event_id))
names(counts) <- c("event_id", "n_stations")
counts$event_id <- as.character(counts$event_id)
ita_stressdrop$event_id <- as.character(ita_stressdrop$event_id)
ita_stressdrop <- merge(ita_stressdrop,counts,by = "event_id",all.x = TRUE)
ita_stressdrop$log_stress_drop <- log10(ita_stressdrop$stress_drop)
stopifnot(!anyDuplicated(ita_stressdrop$event_id))
summary(ita_stressdrop$log_stress_drop)
head(ita_stressdrop)
table(ita_stressdrop$n_stations)

ita_sd_cluster <- merge(
  ita_event_clusters,
  ita_stressdrop[ , c("event_id", "stress_drop")],
  by = "event_id",
  all.x = TRUE,
  sort = FALSE
)
stopifnot(!anyDuplicated(ita_sd_cluster$event_id))
ita_sd_cluster$log_stress_drop <- log10(ita_sd_cluster$stress_drop)

## COMMENT: Sgobba et al. don't specify numerically the "similarity in the stress drop".
## As a consequence, we didn't apply any further merging based on the stress drop.


# Final clustering -------------------------------------------------------------

## Italy

# Plot (Reasenberg + DISS clusters)
ggplot() +
  geom_sf(
    data = world,
    fill = "white",
    color = "grey80"
  ) +
  # Background events
  geom_sf(
    data = subset(ita_diss, is.na(new_cluster)),
    color = "grey70",
    size = 0.8,
    alpha = 0.6
  ) +
  # Clustered events
  geom_sf(
    data = subset(ita_diss, !is.na(new_cluster)),
    aes(color = factor(new_cluster)),
    size = 1.8,
    alpha = 0.9
  ) +
  scale_color_viridis_d(option = "turbo") +
  coord_sf(
    xlim = c(6, 19),
    ylim = c(35, 47.5),
    expand = FALSE
  ) +
  labs(
    color = "Cluster",
    title = "Reasenberg + DISS clusters"
  ) +
  theme(
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.3, "cm")
  )
# ggsave("Plots/Reasenberg+DISS_Italy.png", dpi = 300, width = 10, height = 7)


## Remaining non-clustered events ----

# Sgobba et al. group all Reasenberg background events into a unique cluster

## Ita
ita_sd_cluster$source_term <- ifelse(is.na(ita_sd_cluster$new_cluster), max(ita_sd_cluster$new_cluster, na.rm = TRUE)+1, ita_sd_cluster$new_cluster)
stopifnot(!anyDuplicated(ita_sd_cluster$event_id))

# Final merge
events_missing <- setdiff(unique(ITA18$event), ita_sd_cluster$event_id)
length(events_missing)
ITA18 <- merge(
  ITA18,
  ita_sd_cluster[, c("event_id", "source_term")],
  by.x = "event",
  by.y = "event_id",
  all.x = TRUE,
  sort = FALSE
)
ita_emp <- merge(
  ita_emp,
  ita_sd_cluster[, c("event_id", "source_term")],
  by = "event_id",
  all.x = TRUE,
  sort = FALSE
)

# Save Source term
save(ITA18, df_emp, ita_emp, reas_ita, ita_diss, ita_event_clusters, ita_sd_cluster, file="Results/Clustering/clustering_data.RData")
