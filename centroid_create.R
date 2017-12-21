# Aggregates data points with 'n' distance as centroids
# To be used to aggregate data values within aggregation distnce

foxes_so <- as(foxes,"Sparial") # convert sf to sp
mist<-distm(spTransform(foxes_sp,4326)) # requires in 4326 / lat-long
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

cent_sf<- st_as_sf(cent,coords=c("x","y"), crs = 3577)

st_write(cent_sf,'/home/scottca/Desktop/test.gpkg')
