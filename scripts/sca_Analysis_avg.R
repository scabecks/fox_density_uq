library(nlme)
library(maptools)
library(sp)
library(rgdal)
library(raster)
library(spatial.tools)
library(readr)
library(dplyr)

aus.r<-raster("BrettMurphy_data_/Spatial data/aus_1_1km")
map.r<-raster("BrettMurphy_data_/Spatial data/map_albers") ###Mean annual rainfall
mat.r<-raster("BrettMurphy_data_/Spatial data/mat_albers") ###Mean annual temperature
cover.r<-raster("BrettMurphy_data_/Spatial data/cover_5kmrad") ###MODIS tree cover
rugged.r<-raster("BrettMurphy_data_/Spatial data/rugged_5kmra2") ###Ruggedness, standard deviation of elevation
fox.r<-raster("BrettMurphy_data_/Spatial data/fox_1km")  #Fox range
#cat.r<-raster("BrettMurphy_data_/Spatial data/cat_dist_1km") ###Islands with cats present
#cat.r<-reclassify(cat.r, c(0.9, 1.1, 0.5, 1.9, 2.1, 1, 2.9, 3.1, 0))
#island.area.r<-raster("BrettMurphy_data_/Spatial data/isl_area_1km") ###Island area
#island.area.r<-reclassify(island.area.r, c(-Inf,0.01, 0.01))
#island.size.prop<-island.area.r/10000
#island.size.prop<-reclassify(island.size.prop, c(1,Inf, 1))
#islands.r<-reclassify(island.area.r, c(-Inf,10000, 1, 10000, Inf, 0))
mainland.r<-reclassify(island.area.r, c(-Inf,10000, 0, 10000, Inf, 1))

map.sync<-spatial_sync_raster(map.r, aus.r, method="bilinear")
mat.sync<-spatial_sync_raster(mat.r, aus.r, method="bilinear")
cover.sync<-spatial_sync_raster(cover.r, aus.r, method="bilinear")
rugged.sync<-spatial_sync_raster(rugged.r, aus.r, method="bilinear")

# Add in Fox Data ####
foxes<- read_csv("Ayesha_data/Fox_densities_7Dec2016 (copy).csv",col_types = 'ccnncccddcccnccc')
foxes <- foxes %>% filter(!is.na(GDALat) & !is.na(`Density (per km2)`)) # removes shitty NAs that were causing setting coordinates to fail.
foxes <- foxes %>% rename(year_pub=`Year of Publication`, year_col=`Year of data collection`,lat=GDALat,long=GDALong,density=`Density (per km2)`,confidence=`CONFIDENCE IN DENSITY ESTIMATE`) # rename because I hate typing; note the GDA geography though
coordinates(foxes)<-~long+lat
projection(foxes)<-CRS("+proj=longlat +datum=WGS84")
foxes<-spTransform(foxes, projection(map.r))
# Yep ####
foxes$map<-extract(map.r, foxes, na.rm=T, buffer=5000, fun=mean)
foxes$mat<-extract(mat.r, foxes, na.rm=T, buffer=5000, fun=mean)
foxes$cover<-extract(cover.r, foxes, na.rm=T, buffer=5000, fun=mean)
foxes$rugged<-extract(rugged.r, foxes, na.rm=T, buffer=5000, fun=mean)

#Correct missing values for remote islands IS THIS APPLICABLE?
#cats$map[cats$location=="Macquarie_Island"]<-977
#cats$mat[cats$location=="Macquarie_Island"]<-4.85
#cats$cover[cats$location=="Macquarie_Island"]<-0		#Assume no trees
#cats$mat[cats$location=="Althorpe_Island"]<-16.64		#Stenhouse Bay, 12.3 km away
#cats$mat[cats$location=="Northwest_Island"]<-23.49		#Heron Island Research Station, 25.8 km away
#cats$cover[cats$location=="Northwest_Island"]<-10		#Closest mainland
#cats$rugged[cats$location=="Macquarie_Island"]<-90

#Change the single density of '0' to '0.01'
#cats$density[cats$density==0]<-0.01


#Island occupancy model
#data.is.cats<-read.table("BrettMurphy_data_/Data/island cats.csv", header=T, sep=",") ###Data on which islands cats are definitely present on or absent from
#data.is.cats$area_km2<-ifelse(data.is.cats$area_km2==0,0.01, data.is.cats$area_km2)
#occ.model<-glm(cats~log(area_km2), binomial, data=data.is.cats)

#pred.occ<-1/(1+1/(exp(coef(occ.model)["(Intercept)"]+coef(occ.model)["log(area_km2)"]*log(island.area.r))))

#occupancy<-cat.r
#occupancy[occupancy==0.5]<-pred.occ[occupancy==0.5]
#writeRaster(occupancy, "temp_outputs/cat_occupancy.tif", format="GTiff", overwrite=TRUE)


#Mainland model ####
data<-foxes@data
data$x_albers<-coordinates(foxes)[,1]
data$y_albers<-coordinates(foxes)[,1]
data$x_albers2<-coordinates(foxes)[,1]+runif(length(data[,1]),0,1000)
data$y_albers2<-coordinates(foxes)[,1]+runif(length(data[,1]),0,1000)

#data.ml<-data[data$island.size>60000,]
data.ml <- data %>% filter(density > 0) # just saving on retyping by renaming to data.ml

AICc<-function(model){
K<-length(coef(model))
(AIC(model)+2*K*(K+1)/(length(resid(model))-K-1))}

m1<-lm(log(density)~1, data=data.ml)
m2<-lm(log(density)~log(map), data=data.ml)
m3<-lm(log(density)~mat, data=data.ml)
m4<-lm(log(density)~log(map)+mat, data=data.ml)
m5<-lm(log(density)~cover, data=data.ml)
m6<-lm(log(density)~log(map)+cover, data=data.ml)
m7<-lm(log(density)~mat+cover, data=data.ml)
m8<-lm(log(density)~log(map)+mat+cover, data=data.ml)
m9<-lm(log(density)~log(map)+mat+log(map):mat, data=data.ml)
m10<-lm(log(density)~log(map)+mat+cover+log(map):mat, data=data.ml)
m1r<-lm(log(density)~rugged, data=data.ml)
m2r<-lm(log(density)~log(map)+rugged, data=data.ml)
m3r<-lm(log(density)~mat+rugged, data=data.ml)
m4r<-lm(log(density)~log(map)+mat+rugged, data=data.ml)
m5r<-lm(log(density)~cover+rugged, data=data.ml)
m6r<-lm(log(density)~log(map)+cover+rugged, data=data.ml)
m7r<-lm(log(density)~mat+cover+rugged, data=data.ml)
m8r<-lm(log(density)~log(map)+mat+cover+rugged, data=data.ml)
m9r<-lm(log(density)~log(map)+mat+log(map):mat+rugged, data=data.ml)
m10r<-lm(log(density)~log(map)+mat+cover+log(map):mat+rugged, data=data.ml)

AICcs<-c(AICc(m1), AICc(m2), AICc(m3), AICc(m4), AICc(m5), AICc(m6), AICc(m7), AICc(m8), AICc(m9), AICc(m10), AICc(m11), AICc(m12), AICc(m1r), AICc(m2r), AICc(m3r), AICc(m4r), AICc(m5r), AICc(m6r), AICc(m7r), AICc(m8r), AICc(m9r), AICc(m10r))
di<-AICcs-min(AICcs)
for.wi<-exp(-0.5*di)
wi<-for.wi/sum(for.wi)

R2<-c(summary(m1)[8], summary(m2)[8], summary(m3)[8], summary(m4)[8], summary(m5)[8], summary(m6)[8], summary(m7)[8], summary(m8)[8], summary(m9)[8], summary(m10)[8], summary(m1r)[8], summary(m2r)[8], summary(m3r)[8], summary(m4r)[8], summary(m5r)[8], summary(m6r)[8], summary(m7r)[8], summary(m8r)[8], summary(m9r)[8], summary(m10r)[8])

write.table(di, "temp_outputs/di_avg.csv", sep=",", row.names=F)
write.table(R2, "temp_outputs/R2_avg.csv", sep=",", row.names=F)

coef.int<-wi[1]*coef(m1)["(Intercept)"]+ wi[2]*coef(m2)["(Intercept)"]+ wi[3]*coef(m3)["(Intercept)"]+ wi[4]*coef(m4)["(Intercept)"]+ wi[5]*coef(m5)["(Intercept)"]+ wi[6]*coef(m6)["(Intercept)"]+ wi[7]*coef(m7)["(Intercept)"]+ wi[8]*coef(m8)["(Intercept)"]+ wi[9]*coef(m9)["(Intercept)"]+ wi[10]*coef(m10)["(Intercept)"]+ wi[11]*coef(m1r)["(Intercept)"]+ wi[12]*coef(m2r)["(Intercept)"]+ wi[13]*coef(m3r)["(Intercept)"]+ wi[14]*coef(m4r)["(Intercept)"]+ wi[15]*coef(m5r)["(Intercept)"]+ wi[16]*coef(m6r)["(Intercept)"]+ wi[17]*coef(m7r)["(Intercept)"]+ wi[18]*coef(m8r)["(Intercept)"]+ wi[19]*coef(m9r)["(Intercept)"]+ wi[20]*coef(m10r)["(Intercept)"]

coef.map<-wi[2]*coef(m2)["log(map)"]+ wi[4]*coef(m4)["log(map)"]+ wi[6]*coef(m6)["log(map)"]+ wi[8]*coef(m8)["log(map)"]+ wi[9]*coef(m9)["log(map)"]+ wi[10]*coef(m10)["log(map)"]+ wi[12]*coef(m2r)["log(map)"]+ wi[14]*coef(m4r)["log(map)"]+ wi[16]*coef(m6r)["log(map)"]+ wi[18]*coef(m8r)["log(map)"]+ wi[20]*coef(m10r)["log(map)"]
#
coef.mat<-wi[3]*coef(m3)["mat"]+ wi[4]*coef(m4)["mat"]+ wi[7]*coef(m7)["mat"]+ wi[8]*coef(m8)["mat"]+ wi[9]*coef(m9)["mat"]+ wi[10]*coef(m10)["mat"]+ wi[13]*coef(m3r)["mat"]+ wi[14]*coef(m4r)["mat"]+ wi[17]*coef(m7r)["mat"]+ wi[18]*coef(m8r)["mat"]+ wi[19]*coef(m9r)["mat"]+ wi[20]*coef(m10r)["mat"]
#
coef.cover<-wi[5]*coef(m5)["cover"]+ wi[6]*coef(m6)["cover"]+ wi[7]*coef(m7)["cover"]+ wi[8]*coef(m8)["cover"]+ wi[10]*coef(m10)["cover"]+ wi[15]*coef(m5r)["cover"]+ wi[16]*coef(m6r)["cover"]+ wi[17]*coef(m7r)["cover"]+ wi[18]*coef(m8r)["cover"]+ wi[20]*coef(m10r)["cover"]

coef.mapXmat<-wi[9]*coef(m9)["log(map):mat"]+ wi[10]*coef(m10)["log(map):mat"]+ wi[19]*coef(m9r)["log(map):mat"]+ wi[20]*coef(m10r)["log(map):mat"]
#
coef.rugged<-wi[11]*coef(m1r)["rugged"]+ wi[12]*coef(m2r)["rugged"]+ wi[13]*coef(m3r)["rugged"]+ wi[14]*coef(m4r)["rugged"]+ wi[15]*coef(m5r)["rugged"]+ wi[16]*coef(m6r)["rugged"]+ wi[17]*coef(m7r)["rugged"]+ wi[18]*coef(m8r)["rugged"]+ wi[19]*coef(m9r)["rugged"]+ wi[20]*coef(m10r)["rugged"]


vary.map<-data.frame(
map=seq(min(data.ml$map), max(data.ml$map), (max(data.ml$map)-min(data.ml$map))/500),
mat=rep(mean(data.ml$mat),501),
cover=rep(mean(data.ml$cover),501),
rugged=rep(mean(data.ml$rugged),501))

fit1<-predict(m1, newdata=vary.map)
fit2<-predict(m2, newdata=vary.map)
fit3<-predict(m3, newdata=vary.map)
fit4<-predict(m4, newdata=vary.map)
fit5<-predict(m5, newdata=vary.map)
fit6<-predict(m6, newdata=vary.map)
fit7<-predict(m7, newdata=vary.map)
fit8<-predict(m8, newdata=vary.map)
fit9<-predict(m9, newdata=vary.map)
fit10<-predict(m10, newdata=vary.map)
fit1r<-predict(m1r, newdata=vary.map)
fit2r<-predict(m2r, newdata=vary.map)
fit3r<-predict(m3r, newdata=vary.map)
fit4r<-predict(m4r, newdata=vary.map)
fit5r<-predict(m5r, newdata=vary.map)
fit6r<-predict(m6r, newdata=vary.map)
fit7r<-predict(m7r, newdata=vary.map)
fit8r<-predict(m8r, newdata=vary.map)
fit9r<-predict(m9r, newdata=vary.map)
fit10r<-predict(m10r, newdata=vary.map)


se1<-predict(m1, newdata=vary.map, se.fit=T)$se
se2<-predict(m2, newdata=vary.map, se.fit=T)$se
se3<-predict(m3, newdata=vary.map, se.fit=T)$se
se4<-predict(m4, newdata=vary.map, se.fit=T)$se
se5<-predict(m5, newdata=vary.map, se.fit=T)$se
se6<-predict(m6, newdata=vary.map, se.fit=T)$se
se7<-predict(m7, newdata=vary.map, se.fit=T)$se
se8<-predict(m8, newdata=vary.map, se.fit=T)$se
se9<-predict(m9, newdata=vary.map, se.fit=T)$se
se10<-predict(m10, newdata=vary.map, se.fit=T)$se
se1r<-predict(m1r, newdata=vary.map, se.fit=T)$se
se2r<-predict(m2r, newdata=vary.map, se.fit=T)$se
se3r<-predict(m3r, newdata=vary.map, se.fit=T)$se
se4r<-predict(m4r, newdata=vary.map, se.fit=T)$se
se5r<-predict(m5r, newdata=vary.map, se.fit=T)$se
se6r<-predict(m6r, newdata=vary.map, se.fit=T)$se
se7r<-predict(m7r, newdata=vary.map, se.fit=T)$se
se8r<-predict(m8r, newdata=vary.map, se.fit=T)$se
se9r<-predict(m9r, newdata=vary.map, se.fit=T)$se
se10r<-predict(m10r, newdata=vary.map, se.fit=T)$se


fit.CI.function<-function(fit,se){
cbind(fit, (fit-1.96*se), (fit+1.96*se))}
p1<-fit.CI.function(fit1, se1)
p2<-fit.CI.function(fit2, se2)
p3<-fit.CI.function(fit3, se3)
p4<-fit.CI.function(fit4, se4)
p5<-fit.CI.function(fit5, se5)
p6<-fit.CI.function(fit6, se6)
p7<-fit.CI.function(fit7, se7)
p8<-fit.CI.function(fit8, se8)
p9<-fit.CI.function(fit9, se9)
p10<-fit.CI.function(fit10, se10)
p1r<-fit.CI.function(fit1r, se1r)
p2r<-fit.CI.function(fit2r, se2r)
p3r<-fit.CI.function(fit3r, se3r)
p4r<-fit.CI.function(fit4r, se4r)
p5r<-fit.CI.function(fit5r, se5r)
p6r<-fit.CI.function(fit6r, se6r)
p7r<-fit.CI.function(fit7r, se7r)
p8r<-fit.CI.function(fit8r, se8r)
p9r<-fit.CI.function(fit9r, se9r)
p10r<-fit.CI.function(fit10r, se10r)

max.predicted.mainland<-exp(max(wi[1]*predict(m1)+ wi[2]*predict(m2)+ wi[3]*predict(m3)+ wi[4]*predict(m4)+ wi[5]*predict(m5)+ wi[6]*predict(m6)+ wi[7]*predict(m7)+ wi[8]*predict(m8)+ wi[9]*predict(m9)+ wi[10]*predict(m10)+ wi[11]*predict(m1r)+ wi[12]*predict(m2r)+ wi[13]*predict(m3r)+ wi[14]*predict(m4r)+ wi[15]*predict(m5r)+ wi[16]*predict(m6r)+ wi[17]*predict(m7r)+ wi[18]*predict(m8r)+ wi[19]*predict(m9r)+ wi[20]*predict(m10r)))

pred.map<-p2
write.table(cbind(vary.map$map, pred.map), "temp_outputs/predictions.map_avg.csv", sep=",", row.names=F)

# Again, will this be relevent?
#Islands model
#data.is<-cats@data
#data.is$x_albers<-coordinates(cats)[,1]
#data.is$y_albers<-coordinates(cats)[,1]
#data.is$x_albers2<-coordinates(cats)[,1]+runif(length(data.is[,1]),0,1000)
#data.is$y_albers2<-coordinates(cats)[,1]+runif(length(data.is[,1]),0,1000)
#data.is$island.numeric<-ifelse(data.is$island=="I",1,0)
#data.is$island.size.prop<-ifelse(data.is$island.size>10000,1,data.is$island.size/10000)

#Work out predictions of fox density using mainland model
#mainland.pred<-coef.int+coef.map*log(data.is$map)+coef.mat*data.is$mat+coef.mapXmat*log(data.is$map)*data.is$mat+coef.cover*data.is$cover+coef.fox*data.is$fox+ coef.rugged*data.is$rugged

#AICc<-function(model){
#K<-length(coef(model))
#(AIC(model)+2*K*(K+1)/(length(resid(model))-K-1))}

#m1<-lm(log(density)~mainland.pred, data=data.is)
#m2<-lm(log(density)~mainland.pred+log(island.size.prop), data=data.is)

#AICcs<-c(AICc(m1), AICc(m2))
#(di<-AICcs-min(AICcs))
#for.wi<-exp(-0.5*di)
#wi<-for.wi/sum(for.wi)
#coef.int.is<-wi[1]*coef(m1)["(Intercept)"]+ wi[2]*coef(m2)["(Intercept)"]
#coef.mainland.pred<-wi[1]*coef(m1)["mainland.pred"]+ wi[2]*coef(m2)["mainland.pred"]
#coef.island.size.prop<-wi[2]*coef(m2)["log(island.size.prop)"]

#R2<-c(summary(m1)[8], summary(m2)[8])

#write.table(di, "temp_outputs/di_island_size_avg.csv", sep=",", row.names=F)
#write.table(R2, "temp_outputs/R2_islands_size_avg.csv", sep=",", row.names=F)

#vary.island.size.prop<-data.frame(
#island.size.prop=exp(seq(log(min(data.is$island.size.prop)), log(max(data.is$island.size.prop)), (log(max(data.is$island.size.prop))-log(min(data.is$island.size.prop)))/500)),
#mainland.pred=rep(mean(mainland.pred),501))

#max.predicted.islands<-exp(max(wi[1]*predict(m1)+ wi[2]*predict(m2)))

#fit2<-predict(m2, newdata=vary.island.size.prop)
#se2<-predict(m2, newdata=vary.island.size.prop, se.fit=T)$se
#fit.CI.function<-function(fit,se){
#cbind(fit, (fit-1.96*se), (fit+1.96*se))}
#p2<-fit.CI.function(fit2, se2)
#pred.island.size.prop<-p2
#write.table(cbind(vary.island.size.prop$island.size.prop, pred.island.size.prop), "temp_outputs/predictions.island.size.prop_avg.csv", sep=",", row.names=F)


#Predict total Australian population and mean density
mainland.density<-overlay(aus.r, map.sync, mat.sync, cover.sync, fox.r, rugged.sync, fun=function(x,y,z,z2,z3,z4){
  x*exp(coef.int+coef.map*log(y)+coef.mat*z+coef.mapXmat*log(y)*z+coef.cover*z2+coef.fox*z3+coef.rugged*z4)  })

island.density<-overlay(aus.r, mainland.density, island.size.prop, fun=function(x,y,z){
x*exp(coef.int.is+coef.mainland.pred*log(y)+coef.island.size.prop*log(z))  })

#Ensure we're not extrapolating beyond our data
mainland.density<-reclassify(mainland.density, c(max.predicted.mainland, Inf, max.predicted.mainland))
island.density<-reclassify(island.density, c(max.predicted.islands, Inf, max.predicted.islands))

#Correct for island occupancy
#island.density.occupancy<-overlay(island.density, occupancy, fun=function(x,y){x*y})

#Combine mainland and island estimates
#density<-overlay(mainland.density, mainland.r, island.density.occupancy, islands.r, fun=function(w,x,y,z){w*x+y*z})

writeRaster(mainland.density, "temp_outputs/mainland_density_avg.tif", format="GTiff", overwrite=TRUE)
#writeRaster(island.density, "temp_outputs/island_density_avg.tif", format="GTiff", overwrite=TRUE)
#writeRaster(island.density.occupancy, "temp_outputs/island_density_occupancy_avg.tif", format="GTiff", overwrite=TRUE)
#writeRaster(density, "temp_outputs/density_avg.tif", format="GTiff", overwrite=TRUE)

plot(mainland.density)
plot(foxes, bg="transparent", add=TRUE)

(mean.density<-cellStats(density, mean))
(total.population<-mean.density*cellStats(aus.r, sum))


#Confidence intervals of continent-wide density (based on boot-strapping dataset)
data<-foxes@data
#data$island.size.prop<-ifelse(data$island.size/10000>1,1,data$island.size/10000)

reduced.res<-aus.r
res(reduced.res)<-10000
aus.coarse<-resample(aus.r, reduced.res, method="bilinear")

map.coarse<-spatial_sync_raster(map.r, aus.coarse, method="bilinear")
mat.coarse<-spatial_sync_raster(mat.r, aus.coarse, method="bilinear")
cover.coarse<-spatial_sync_raster(cover.r, aus.coarse, method="bilinear")
#fox.coarse<-spatial_sync_raster(fox.r, aus.coarse, method="bilinear")
rugged.coarse<-spatial_sync_raster(rugged.r, aus.coarse, method="bilinear")
#island.size.prop.coarse<-spatial_sync_raster(island.size.prop, aus.coarse, method="bilinear")
#occupancy.coarse<-spatial_sync_raster(occupancy, aus.coarse, method="bilinear")

mainland.coarse<-spatial_sync_raster(mainland.r, aus.coarse, method="bilinear")
#islands.coarse<-spatial_sync_raster(islands.r, aus.coarse, method="bilinear")


AICc<-function(model){
K<-length(coef(model))
(AIC(model)+2*K*(K+1)/(length(resid(model))-K-1))}

density.func<-function(x) {
#x[,1]=density, x[,2]=map, x[,3]=mat, x[,4]=cover, x[,5]=fox, x[,6]= island.size.prop

mainland.data<-x[(x[,6])==1,]

#Mainland model
m.density<-mainland.data [,1]
m.map<-mainland.data [,2]
m.mat<-mainland.data [,3]
m.cover<-mainland.data [,4]
m1<-lm(log(m.density)~1)
m2<-lm(log(m.density)~log(m.map))
m3<-lm(log(m.density)~m.mat)
m4<-lm(log(m.density)~log(m.map)+m.mat)
m5<-lm(log(m.density)~m.cover)
m6<-lm(log(m.density)~log(m.map)+m.cover)
m7<-lm(log(m.density)~m.mat+m.cover)
m8<-lm(log(m.density)~log(m.map)+m.mat+m.cover)
m9<-lm(log(m.density)~log(m.map)+m.mat+log(m.map):m.mat)
m10<-lm(log(m.density)~log(m.map)+m.mat+m.cover+log(m.map):m.mat)
AICcs<-c(AICc(m1), AICc(m2), AICc(m3), AICc(m4), AICc(m5), AICc(m6), AICc(m7), AICc(m8), AICc(m9), AICc(m10))
di<-AICcs-min(AICcs)
for.wi<-exp(-0.5*di)
wi<-for.wi/sum(for.wi)

coef.int<-wi[1]*coef(m1)["(Intercept)"]+ wi[2]*coef(m2)["(Intercept)"]+ wi[3]*coef(m3)["(Intercept)"]+ wi[4]*coef(m4)["(Intercept)"]+ wi[5]*coef(m5)["(Intercept)"]+ wi[6]*coef(m6)["(Intercept)"]+ wi[7]*coef(m7)["(Intercept)"]+ wi[8]*coef(m8)["(Intercept)"]+ wi[9]*coef(m9)["(Intercept)"]+ wi[10]*coef(m10)["(Intercept)"]

coef.map<-wi[2]*coef(m2)["log(m.map)"]+ wi[4]*coef(m4)["log(m.map)"]+ wi[6]*coef(m6)["log(m.map)"]+ wi[8]*coef(m8)["log(m.map)"]+ wi[9]*coef(m9)["log(m.map)"]+ wi[10]*coef(m10)["log(m.map)"]

coef.mat<-wi[3]*coef(m3)["m.mat"]+ wi[4]*coef(m4)["m.mat"]+ wi[7]*coef(m7)["m.mat"]+ wi[8]*coef(m8)["m.mat"]+ wi[9]*coef(m9)["m.mat"]+ wi[10]*coef(m102)["m.mat"]

coef.cover<-wi[5]*coef(m5)["m.cover"]+ wi[6]*coef(m6)["m.cover"]+ wi[7]*coef(m7)["m.cover"]+ wi[8]*coef(m8)["m.cover"]+ wi[10]*coef(m10)["m.cover"]

coef.mapXmat<-wi[9]*coef(m9)["log(m.map):m.mat"]+ wi[10]*coef(m10)["log(m.map):m.mat"]

max.predicted.mainland<-exp(max(wi[1]*predict(m1)+ wi[2]*predict(m2)+ wi[3]*predict(m3)+ wi[4]*predict(m4)+ wi[5]*predict(m5)+ wi[6]*predict(m6)+ wi[7]*predict(m7)+ wi[8]*predict(m8)+ wi[9]*predict(m9)+ wi[10]*predict(m10)))

#Island model
#i.density<-x[,1]
#i.map<-x[,2]
#i.mat<-x[,3]
#i.cover<-x[,4]
#i.fox<-x[,5]
#i.island.size.prop<-x[,6]
#mainland.pred<-coef.int+coef.map*log(i.map)+coef.mat*i.mat+coef.mapXmat*log(i.map)*i.mat+coef.cover*i.cover+coef.fox*i.fox
#m1<-lm(log(i.density)~mainland.pred)
#m2<-lm(log(i.density)~mainland.pred+log(i.island.size.prop))
#di<- c(AICc(m1), AICc(m2))-min(c(AICc(m1), AICc(m2)))
#for.wi<-exp(-0.5*di)
#wi<-for.wi/sum(for.wi)
#coef.int.is<-wi[1]*coef(m1)["(Intercept)"]+ wi[2]*coef(m2)["(Intercept)"]
#coef.mainland.pred<-wi[1]*coef(m1)["mainland.pred"]+ wi[2]*coef(m2)["mainland.pred"]
#coef.island.size.prop<-wi[2]*coef(m2)["log(i.island.size.prop)"]
#max.predicted.islands<-exp(max(wi[1]*predict(m1)+ wi[2]*predict(m2)))

#Predictions for the whole continent
mainland.density<-overlay(map.coarse, mat.coarse, cover.coarse, fun=function(y,z,z2,z3){
exp(coef.int+coef.map*log(y)+coef.mat*z+coef.mapXmat*log(y)*z+coef.cover*z2)  })

#island.density<-overlay(mainland.density, island.size.prop.coarse, fun=function(y,z){
#exp(coef.int.is+coef.mainland.pred*log(y)+coef.island.size.prop*log(z))  })

mainland.density[mainland.density>max.predicted.mainland]<-max.predicted.mainland
island.density[island.density>max.predicted.islands]<-max.predicted.islands

density<-mainland.density*mainland.coarse#+island.density*occupancy.coarse*islands.coarse
cellStats(density, mean)}

#density.func(cbind(data$density, data$map, data$mat, data$cover, data$fox, data$island.size.prop))


bootstrap.func<-function(x) {
#x[,1]=density, x[,2]=map, x[,3]=mat, x[,4]=cover, x[,5]=fox, x[,6]=island.size.prop
data.ml<-x[(x[,6])==1,]
n.ml<-length(data.ml[,1])
random.nos<-floor(runif(n.ml,0,n.ml))+1
randomised.ml<-data.ml[random.nos,]
data.is<-x[(x[,6])<1,]
n.is<-length(data.is[,1])
random.nos<-floor(runif(n.is,0,n.is))+1
randomised.is<-data.is[random.nos,]
rbind(randomised.ml, randomised.is)}

#bootstrap.func(cbind(data$density, data$map, data$mat, data$cover, data$fox, data$island.size.prop))

#density.func(bootstrap.func(cbind(data$density, data$map, data$mat, data$cover, data$fox, data$island.size.prop)))


sim.length<-10000
sim.densities<-numeric(sim.length)
for (i in 1:sim.length) sim.densities[i]<-density.func(bootstrap.func(cbind(data$density, data$map, data$mat, data$cover)))#, data$fox, data$island.size.prop)))
quantile(sim.densities,c(0.025, 0.975))
quantile(sim.densities,c(0.025, 0.975))*cellStats(aus.r, sum)
hist(sim.densities*cellStats(aus.r, sum), xlim=c(0,5), breaks=10000)


write.table(c(cellStats(aus.r, sum), sim.densities), "temp_outputs/simulation_avg.csv", sep=",", row.names=F)


