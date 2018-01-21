# Aggregates data points with 'n' distance as centroids
# To be used to aggregate data values within aggregation distnce
pacman::p_load(sf, rgeos, geosphere, sp, dplyr, rpostgis, ggplot2,raster)

foxes <- read_sf("foxes_points.gpkg")
foxes <- foxes %>%
  mutate(density=as.numeric(ifelse(density=='TBA',NA,density))) %>%
  st_set_crs(3577)

foxes_sp <- as(foxes,"Spatial") # convert sf to sp
mdist <- distm(spTransform(foxes_sp,CRS("+init=epsg:4326"))) # requires input in 4326 / lat-long
hc<-hclust(as.dist(mdist),method = 'complete') # cluster by distance
foxes_sp$clust <- cutree(hc, h = 1000) # 1km disrance

#create matrix
cent <- matrix(ncol=2, nrow=max(foxes_sp$clust))

for (i in 1:max(foxes_sp$clust)){
  # gCentroid from the rgeos package
  cent[i,] <- gCentroid(subset(foxes_sp, clust == i))@coords
}

cent <- as.data.frame(cent)
names(cent) <- c("x","y")

# convert to sf object
cent_sf <- st_as_sf(cent,coords=c("x","y"), crs = 3577)

# save to geopackage
st_write(cent_sf,'foxes_cluster_centroids.gpkg', layer_options = "OVERWRITE=YES")
cent_sf <- read_sf("foxes_cluster_centroids.gpkg") %>% 
  mutate(gid = row_number()) %>% 
  st_set_crs(3577)# if need to read back in; adds a unique ID to aggregate by

# Aggregate values within 1km of aggregated centroids
st_join(foxes, st_buffer(cent_sf, dist = 1000)) %>% # 1km buffer
  group_by(gid) %>%
  summarise(cnt=n(), min_den = min(density), max_den = max(density), mean_density=mean(density)) %>% 
  st_set_crs(3577) -> foxes_agg

write_sf(foxes_agg, "foxes_points_aggregated.gpkg")              

#### some plots; note you need the development version of ggplot to use the sf geoms....
aus_boundary <- st_transform(st_intersection(st_as_sf(raster::getData("GADM", country = "AU", level = 0), crs=4326), st_sfc(st_polygon(list(rbind(c(110,-45), c(110,-8), c(154,-8), c(154,-45), c(110,-45)))), crs=4326)), crs = 3577) # this clips the area down to not include Macquerie and Christmas Islands in the AUS extent

# basic plot
ggplot(aus_boundary) +
  geom_sf() +
  geom_sf(data = foxes_agg, aes(colour=mean_density)) +
  scale_color_gradient(low = 'blue', high = 'red')

