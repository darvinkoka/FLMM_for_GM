################################################################################
# PATH - Corrective Term
################################################################################

# Initial Setup ----------------------------------------------------------------

library(sf)
library(dplyr)
library(rnaturalearth)
library(ggplot2)
library(viridis)
library(grid)

setwd("C:/Users/stefy/OneDrive - Politecnico di Milano/Desktop/Dottorato/Corsi PhD/ASLCD/FLMM_for_GM")
load("C:/Users/stefy/OneDrive - Politecnico di Milano/Desktop/Dottorato/Corsi PhD/ASLCD/FLMM_for_GM/Results/Clustering/clustering_data.RData")
# source("Code/Source_term.R")


# Preprocessing ----------------------------------------------------------------

path_data <- unique(ita_emp[, c("event_id", "source_term", "station_code", "st_latitude", "st_longitude")])
step <- 0.2

lat_min <- floor(min(path_data$st_latitude) / step) * step # to have "clean" 0.2-cells
lat_max <- ceiling(max(path_data$st_latitude) / step) * step
lon_min <- floor(min(path_data$st_longitude) / step) * step
lon_max <- ceiling(max(path_data$st_longitude) / step) * step

n_lat_cells <- round((lat_max - lat_min) / step)
n_lon_cells <- round((lon_max - lon_min) / step)


# Bounding box (Sgobba et al.) -------------------------------------------------

path_data$lat_cell <- pmin(floor((path_data$st_latitude - lat_min) / step) + 1L, n_lat_cells)
path_data$lon_cell <- pmin(floor((path_data$st_longitude - lon_min) / step) + 1L, n_lon_cells)
path_data$cell_id <- paste(path_data$lat_cell, path_data$lon_cell, sep = "_")
length(unique(path_data$cell_id))
path_data$path_term <- paste(path_data$source_term, path_data$cell_id, sep = "_")
length(unique(path_data$path_term))
stopifnot(!anyDuplicated(path_data[, c("event_id", "station_code")]))
ita_emp <- merge(ita_emp, path_data[,c("event_id", "station_code", "path_term")], by = c("event_id", "station_code"), all.x = TRUE)

path_lookup <- unique(path_data[, c("event_id", "station_code", "path_term")])
anyDuplicated(path_lookup[, c("event_id", "station_code")])
ITA18$event <- as.character(ITA18$event)
ITA18$Station <- as.character(ITA18$Station)
path_lookup$event_id <- as.character(path_lookup$event_id)
path_lookup$station_code <- as.character(path_lookup$station_code)
ITA18 <- merge(
  ITA18,
  path_lookup,
  by.x = c("event", "Station"),
  by.y = c("event_id", "station_code"),
  all.x = TRUE,
  sort = FALSE
)

# Save final data
save(ITA18, df_emp, ita_emp, reas_ita, ita_diss, ita_event_clusters, path_data, ita_sd_cluster, file="Results/Clustering/source-path_data.RData")
# save(ITA18, df_emp, ita_emp, reas_ita, ita_diss, ita_event_clusters, path_data, ita_sd_cluster, file="Risultati/source-path_data.RData")


# Plot ----

# 1) Map of paths:
source_coord <- unique(ita_emp[, c("event_id", "source_term", "ev_latitude", "ev_longitude")])
path_plot <- merge(path_data, source_coord, by = c("event_id", "source_term"))
path_unique <- unique(path_plot[, c("source_term", "ev_latitude", "ev_longitude", "lat_cell", "lon_cell")])

world <- ne_countries(scale = "medium", returnclass = "sf")

ggplot() +
  geom_sf(
    data = world,
    fill = "white",
    color = "grey80"
  ) +
  # path
  geom_segment(
    data = path_plot,
    aes(
      x = ev_longitude,
      y = ev_latitude,
      xend = st_longitude,
      yend = st_latitude,
      color = factor(source_term)
    ),
    alpha = 0.2
  ) +
  # source
  geom_point(
    data = source_coord,
    aes(
      x = ev_longitude,
      y = ev_latitude,
      color = factor(source_term)
    ),
    shape = 24,
    fill = "white",
    size = 2
  ) +
  # station
  geom_point(
    data = path_data,
    aes(
      x = st_longitude,
      y = st_latitude
    ),
    size = 1
  ) +
  coord_sf(
    xlim = c(6, 19),
    ylim = c(35, 48)
  ) +
  theme_bw() +
  labs(
    color = "Map of paths",
    x = "Longitude",
    y = "Latitude"
  )
ggsave("Plots/Path_Italy.png", dpi = 300, width = 10, height = 7)


# 2) Earthquakes in the study area
ggplot() +
  geom_sf(data = world, fill = "grey97", colour = "grey75",
          linewidth = 0.3) +
  geom_sf(data = reas_ita_sf, colour = "firebrick",
          size = 1.5, alpha = 0.7) +
  coord_sf(xlim = c(6, 19), ylim = c(35, 47.5), expand = FALSE) +
  labs(
    title = "Earthquakes in the study area",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal()
ggsave("Plots/Earthquakes_Italy.png", dpi = 300, width = 10, height = 7)


# Save Path term
save(ITA18, df_emp, ita_emp, reas_ita, ita_diss, ita_event_clusters, ita_sd_cluster, path_data, file="Results/Clustering/clustering_data.RData")

