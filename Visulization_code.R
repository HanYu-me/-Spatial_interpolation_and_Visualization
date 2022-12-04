library(data.table)
library(dismo)
library(raster)
library(sp)
library(rgdal)
library(gstat)
library(raster)
library(maptools)
library(RColorBrewer)
library(grid)
library(lattice)

climate=read.csv(file.choose(),header=T)
names(climate)[1:3]=c("id","long","lat")
#set the range of longitude and latitude
climate=climate[climate$long<40&climate$long>-15,]
climate=climate[climate$lat<63&climate$lat>30,]
climate=na.omit(climate)
#covert dataframe to sp form
dsp <- SpatialPoints(climate[,c(2,3)], proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
climate_sp <- SpatialPointsDataFrame(dsp,climate)
climate_sp=spTransform(climate_sp, CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
#load SHP world map and convert to sp form
worldMap=readOGR(file.choose())
worldSP=spTransform(worldMap, CRS=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
#plot point projection
pols <- list("sp.polygons",worldSP, fill = "lightgray")
colors=rev(colorRampPalette(brewer.pal(5,"RdBu"))(20))
spplot(climate_sp,"bio_1", col.regions=colors, sp.layout=pols, pch=20, cex=1)

##interpolation
#convert point data to Voronor geometry
v <- voronoi(climate_sp)
plot(v)
worldMap3<-aggregate(worldSP)#aggregation to reduce resolution
v1<-intersect(v,worldSP)#merge geometry to world map
spplot(v1, "bio_1", col.regions=colors())#draw raw data
#set plot margin
worldSP@bbox[1]=-15
worldSP@bbox[2]=30
worldSP@bbox[3]=40
worldSP@bbox[4]=63
#set raster number as 500*500
blank_raster<-raster(nrow=500,ncol=500,extent(worldSP))
values(blank_raster)<-1
bound_raster<-rasterize(worldSP,blank_raster)
bound_raster[!(is.na(bound_raster))] <- 1
plot(bound_raster)
#LDW
i="bio_1"
#project V1 raster to 500*500 raster
vr <- rasterize(v1,bound_raster,i)
#do interpolation
gs <- gstat(formula=as.formula(paste0("climate_sp$",i,"~1")), locations=climate_sp)
idw <- interpolate(bound_raster, gs)
idwmask<-mask(idw,vr)
colors=rev(c(colorRampPalette(brewer.pal(5,"RdBu"))(23)))
SP <- spplot(idwmask,col.regions=colors, colorkey = list(space = "left", height = 0.4))
SP
## put legent into figure
args <- SP$legend$left$args$key
## Prepare list of arguments needed by `legend=` argument (as described in ?xyplot)
legendArgs <- list(fun = draw.colorkey,args = list(key = args),corner = c(0,.45))
## Call spplot() again, this time passing in to legend the arguments
spplot(idwmask, colorkey = FALSE,legend = list(inside = legendArgs),col.regions=colors,main=i)



