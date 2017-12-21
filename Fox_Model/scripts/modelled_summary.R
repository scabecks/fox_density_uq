library(raster)
library(spatial.tools)
library(sf)

setwd("~/Dropbox/UQ_Work/Projects/fox_model_working/Fox Model/")

#Modelled fox densities
foxes_density <- raster("mainland_density_avg_5.06.2017.tif")

#Load in other spatial data used in model
map.r <- raster("BrettMurphy_data_/Spatial data/map_albers") ###Mean annual rainfall
mat.r <- raster("BrettMurphy_data_/Spatial data/mat_albers") ###Mean annual temperature
cover.r <- raster("BrettMurphy_data_/Spatial data/cover_5kmrad") ###MODIS tree cover
rugged.r<-raster("BrettMurphy_data_/Spatial data/rugged_5kmra2") 

# Create synced versions
map.sync <- spatial_sync_raster(map.r, foxes_density, method="bilinear")
mat.sync <- spatial_sync_raster(mat.r, foxes_density, method="bilinear")
cover.sync <- spatial_sync_raster(cover.r, foxes_density, method="bilinear")
rugged.sync <- spatial_sync_raster(rugged.r, foxes_density, method="bilinear")
# May need to add in a landuse layer, e.g. urban/rural type thing

# Load in points data (yay! for Geopackages - down with Shapefile!)
foxes_points <- read_sf("foxes_points.gpkg", layer = "foxes")

#Extract values from rasters at points
foxes_points$modelled_density <- extract(foxes_density,foxes_points, na.rm=T)
foxes_points$rain <- extract(map.sync, foxes_points, na.rm=T)
foxes_points$temp <- extract(mat.sync, foxes_points, na.rm=T)
foxes_points$tree_cover <- extract(cover.sync, foxes_points, na.rm=T)
foxes_points$rugged <- extract(rugged.sync, foxes_points, na.rm=T)
foxes_points$x  <- coordinates(foxes_points)[,1]
foxes_points$y <- coordinates(foxes_points)[,2]

(foxes_points %>% 
  group_by(x, y) %>% 
  summarise(n_readings=n(),mean_mod_den=mean(modelled_density),max_mod_den=max(modelled_density),mean_rain=mean(rain),mean_temp=mean(temp), mean_tree=mean(tree_cover), mean_rugg=mean(rugged)) %>% 
    write_csv("modelled_sumamry.csv")) 

