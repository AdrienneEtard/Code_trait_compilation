## Extracting range sizes from IUCN distribution maps for reptiles

# Load packages -----------------------------------------------------------
X <- c("ggmap", "rgdal", "maptools", "rgeos", "dplyr", "raster", "maps", "geosphere", "sp", "sf", "lwgeom", "mapproj")
lapply(X, library, character.only = TRUE)
rm(X)

# Read shapefile object as sf ---------------------------------------------
Reptiles <- st_as_sf(readOGR(dsn="../Data/Range_sizes/Reptiles_IUCN_range_size", layer="REPTILES"))
rownames(Reptiles) <- c(1:nrow(Reptiles))
st_is_longlat(Reptiles)


# Project to conserve surface areas (Albers equal area)
Reptiles_aea <- st_transform(Reptiles, "+proj=aea")

# Areas? (in m2)
Areas <- as.data.frame(Reptiles_aea$binomial)
colnames(Areas) <- "Binomial_name"
Areas$Range_size_m2 <- as.numeric(st_area(Reptiles_aea))

# Group by species name, sum surface areas over species binomial names
Areas  %<>% group_by(Binomial_name) %>% summarise(Range_size_m2=sum(Range_size_m2))

# Export file .csv --------------------------------------------------------
write.csv(Areas, "../../Data/Range_sizes/reptile_range_areas.csv", row.names=FALSE)




