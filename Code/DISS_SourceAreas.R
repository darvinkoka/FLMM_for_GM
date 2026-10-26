################################################################################
# DISS 3.2.1. - Dataset
################################################################################

# Initial Setup ----------------------------------------------------------------
library(sf)
library(ggplot2)
library(rnaturalearth)

setwd("C:/Users/stefy/OneDrive - Politecnico di Milano/Desktop/Dottorato/Corsi PhD/ASLCD/FLMM_for_GM")
st_layers("Data/DISS_3.2.1_shp/SHP") 
#   layer_name geometry_type features fields crs_name
# 1  CSSPLN321       Polygon      188     23   WGS 84
# 2  CSSTOP321   Line String      188     23   WGS 84
# 3     DSS321       Polygon       35      5   WGS 84
# 4     ISS321       Polygon      127     39   WGS 84 # --> this!
# 5 SUBDCNT321   Line String      156      3   WGS 84
# 6 SUBDZON321       Polygon        3     19   WGS 84

diss <- read_sf("Data/DISS_3.2.1_shp/SHP", layer = "ISS321")
italy <- ne_countries(scale = "medium", country = "Italy", returnclass = "sf")

# Plot source zones
ggplot() +
  geom_sf(data = italy, fill = "gray90") +
  geom_sf(data = diss, color = "red", fill = NA, size = 0.3) +
  theme_minimal()
ggsave("Plots/DISS_SourceAreas.png", dpi = 300, width = 10, height = 7)

