library(leastcostpath)
library(terra)
library(sf)

source("./R/Functions.R")

# aggregation factor for the Digital Elevation Model (DEM), e.g. a value of 10 means that the 10m DEM is scaled to 100m
agg_val <- 10

if(!file.exists(paste0("./Data/SARDINIA_", agg_val*10,"m.tif"))) {
  
  if(agg_val <= 1) { 
    stop("agg_val must have a value of 2 or greater")
    }
  
  library(tinitalyR, lib = "/local/filespace/workspace/alexis/Joe/phd/packages")
  library(raster, lib = "/local/filespace/workspace/alexis/Joe/phd/packages")
  
  tiles <- c("w45540_s10", "w45040_s10", "w44540_s10",
           "w44040_s10", "w43040_s10", "w45045_s10",
           "w44545_s10", "w44045_s10", "w43545_s10",
           "w43045_s10", "w45550_s10", "w45050_s10",
           "w44550_s10", "w44050_s10", "w43550_s10",
           "w43050_s10", "w45055_s10", "w44555_s10",
           "w44055_s10")

  dem <- tinitalyR::download_dem(x = tiles)
  dem1 <- do.call(raster::merge, dem)

  dem2 <- raster::aggregate(x = dem1, agg_val)

  raster::writeRaster(dem2, filename= paste0("./Data/SARDINIA_", agg_val*10,"m.tif"), format="GTiff", overwrite=TRUE)

  sard_dem <- terra::rast(paste0("./Data/SARDINIA_", agg_val*10,"m.tif"))  
  

} else {
  
  sard_dem <- terra::rast(paste0("./Data/SARDINIA_", agg_val*10,"m.tif"))
  
}

sard_roads_republic <- sf::st_read("./Data/Roman_roads_MASTINO_Republic.shp")
sard_roads_empire_1 <- sf::st_read("./Data/Roman_roads_MASTINO_Empire_1.shp")
sard_roads_empire_1_sulci <- sf::st_read("./Data/Roman_roads_MASTINO_Empire_1_sulci.shp")
sard_roads_empire_1_west <- sf::st_read("./Data/Roman_roads_MASTINO_Empire_1_west.shp")
sard_roads_empire_1_central <- sf::st_read("./Data/Roman_roads_MASTINO_Empire_1_central.shp")
sard_roads_empire_2 <- sf::st_read("./Data/Roman_roads_MASTINO_Empire_2.shp")

sard_roads <- sf::st_read("./Data/Roman_roads_MASTINO.shp")
sard_coast <- sf::st_read("./Data/sardegna.shp")
sard_stations <- sf::st_read("./Data/stations.shp")
carales <- sard_stations[sard_stations$Name == "Carales",]
olbia <- sard_stations[sard_stations$Name == "Olbia",]

other_locs <- sf::st_read("./Data/other_locations.shp")
turris_libisonis <- other_locs[other_locs$Name %in% c("Turris Libisonis"),]

sard_dem2 <- sard_dem
sard_dem2[sard_dem2 <= 0 ] <- NA

terra::writeRaster(x = sard_dem2, filename = "./Data/SARDINIA2_100m.tif", overwrite = TRUE)

####################################
#### BASE CONDUCTANCE SURFACE ######
####################################

slope_cs0 <- leastcostpath::create_slope_cs(x = sard_dem2, cost_function = "tobler", neighbours = 16)

slope_cs0_rast <- leastcostpath::rasterise(x = slope_cs0)
terra::writeRaster(x = slope_cs0_rast, filename = paste0("./Output/Rasters/", agg_val * 10, "m_conductance_surface0.tif"), overwrite = TRUE)

accum_cost0a <- create_accum_cost(x = slope_cs0, loc = carales, filename = paste0("./Output/Rasters/", agg_val * 10, "m_accum_cost0_carales.tif"), mask = sard_coast)
accum_cost0b <- create_accum_cost(x = slope_cs0, loc = olbia, filename = paste0("./Output/Rasters/", agg_val * 10, "m_accum_cost0_olbia.tif"), mask = sard_coast)
accum_cost0c<- create_accum_cost(x = slope_cs0, loc = turris_libisonis, filename = paste0("./Output/Rasters/", agg_val * 10, "m_accum_cost0_turris_libisonis.tif"), mask = sard_coast)

accum_cost0_lyrs <- c(accum_cost0a, accum_cost0b, accum_cost0c)
accum_cost0 <- terra::app(accum_cost0_lyrs, fun = function(x) min(x))

terra::writeRaster(x = accum_cost0, filename = paste0("./Output/Rasters/", agg_val * 10, "m_accum_cost0.tif"), overwrite = TRUE)

####################################
######## REPUBLIC ROADS ############
####################################

slope_cs1 <- slope_cs0

# convert 3.75km per hour to metres per second by dividing by 3.6. Divided by maximum resolution of DEM (1000m) to make units comparable with speed calculated using leastcostpath (leastcostpath accounts for distance between cells when calculating speed)
slope_cs1 <- leastcostpath::update_values(x = slope_cs1, sf = sard_roads_republic, 
                                          FUN = function(x) { replace(x = x, values = (3.75/3.6)/max(res(sard_dem)))})

slope_cs1_rast <- leastcostpath::rasterise(x = slope_cs1)
terra::writeRaster(x = slope_cs1_rast, filename = paste0("./Output/Rasters/", agg_val * 10, "m_conductance_surface1.tif"), overwrite = TRUE)

accum_cost1a <- create_accum_cost(x = slope_cs1, loc = carales, filename = paste0("./Output/Rasters/", agg_val * 10, "m_accum_cost1_carales.tif"), mask = sard_coast)
accum_cost1b <- create_accum_cost(x = slope_cs1, loc = olbia, filename = paste0("./Output/Rasters/", agg_val * 10, "m_accum_cost1_olbia.tif"), mask = sard_coast)
accum_cost1c <- create_accum_cost(x = slope_cs1, loc = turris_libisonis, filename = paste0("./Output/Rasters/", agg_val * 10, "m_accum_cost1_turris_libisonis.tif"), mask = sard_coast)

accum_cost1_lyrs <- c(accum_cost1a, accum_cost1b, accum_cost1c)
accum_cost1 <- terra::app(accum_cost1_lyrs, fun = function(x) min(x))

terra::writeRaster(x = accum_cost1, filename = paste0("./Output/Rasters/", agg_val * 10, "m_accum_cost1.tif"), overwrite = TRUE)

accum_cost0_Olbia <- terra::rast(paste0("./Output/Rasters/", agg_val * 10, "m_accum_cost0_olbia.tif"))
accum_cost1_Olbia <- terra::rast(paste0("./Output/Rasters/", agg_val * 10, "m_accum_cost1_olbia.tif"))

diff_accum_cost_01_olbia <- ((accum_cost0_Olbia-accum_cost1_Olbia) / accum_cost1_Olbia) * 100

terra::writeRaster(x = diff_accum_cost_01_olbia, filename = paste0("./Output/Rasters/", agg_val * 10, "m_diff_accum_cost_01_olbia.tif"), overwrite = TRUE)

####################################
######## EMPIRE ROADS ONE ##########
####################################

slope_cs2 <- slope_cs0

slope_cs2 <- leastcostpath::update_values(x = slope_cs2, sf = rbind(sard_roads_republic, sard_roads_empire_1, sard_roads_empire_1_sulci, sard_roads_empire_1_west, sard_roads_empire_1_central), 
                                          FUN = function(x) { replace(x = x, values = (3.75/3.6)/max(res(sard_dem)))})

slope_cs2_rast <- leastcostpath::rasterise(x = slope_cs2)
terra::writeRaster(x = slope_cs2_rast, filename = paste0("./Output/Rasters/", agg_val * 10, "m_conductance_surface2.tif"), overwrite = TRUE)

accum_cost2a <- create_accum_cost(x = slope_cs2, loc = carales, filename = paste0("./Output/Rasters/", agg_val * 10, "m_accum_cost2_carales.tif"), mask = sard_coast)
accum_cost2b <- create_accum_cost(x = slope_cs2, loc = olbia, filename = paste0("./Output/Rasters/", agg_val * 10, "m_accum_cost2_olbia.tif"), mask = sard_coast)
accum_cost2c <- create_accum_cost(x = slope_cs2, loc = turris_libisonis, filename = paste0("./Output/Rasters/", agg_val * 10, "m_accum_cost2_turris_libisonis.tif"), mask = sard_coast)

accum_cost2_lyrs <- c(accum_cost2a, accum_cost2b, accum_cost2c)
accum_cost2 <- terra::app(accum_cost2_lyrs, fun = function(x) min(x))

terra::writeRaster(x = accum_cost2, filename = paste0("./Output/Rasters/", agg_val * 10, "m_accum_cost2.tif"), overwrite = TRUE)

accum_cost_1_Carales <- terra::rast(paste0("./Output/Rasters/", agg_val * 10, "m_accum_cost1_carales.tif"))
accum_cost_2_Carales <- terra::rast(paste0("./Output/Rasters/", agg_val * 10, "m_accum_cost2_carales.tif"))

diff_accum_cost_12_Carales <- ((accum_cost_1_Carales-accum_cost_2_Carales) / accum_cost_2_Carales) * 100

terra::writeRaster(x = diff_accum_cost_12_Carales, filename = paste0("./Output/Rasters/", agg_val * 10, "m_diff_accum_cost_12_Carales.tif"), overwrite = TRUE)

####################################
######## EMPIRE ROADS TWO ##########
####################################

slope_cs3 <- slope_cs0

slope_cs3 <- leastcostpath::update_values(x = slope_cs3, sf = rbind(sard_roads_republic, sard_roads_empire_1, sard_roads_empire_1_sulci, sard_roads_empire_1_west, sard_roads_empire_1_central, sard_roads_empire_2), 
                                          FUN = function(x) { replace(x = x, values = (3.75/3.6)/max(res(sard_dem)))})

slope_cs3_rast <- leastcostpath::rasterise(x = slope_cs3)
terra::writeRaster(x = slope_cs3_rast, filename = paste0("./Output/Rasters/", agg_val * 10, "m_conductance_surface3.tif"), overwrite = TRUE)

accum_cost3a <- create_accum_cost(x = slope_cs3, loc = carales, filename = paste0("./Output/Rasters/", agg_val * 10, "m_accum_cost3_carales.tif"), mask = sard_coast)
accum_cost3b <- create_accum_cost(x = slope_cs3, loc = olbia, filename = paste0("./Output/Rasters/", agg_val * 10, "m_accum_cost3_olbia.tif"), mask = sard_coast)
accum_cost3c <- create_accum_cost(x = slope_cs3, loc = turris_libisonis, filename = paste0("./Output/Rasters/", agg_val * 10, "m_accum_cost3_turris_libisonis.tif"), mask = sard_coast)

accum_cost3_lyrs <- c(accum_cost3a, accum_cost3b, accum_cost3c)
accum_cost3 <- terra::app(accum_cost3_lyrs, fun = function(x) min(x))

terra::writeRaster(x = accum_cost3, filename = paste0("./Output/Rasters/", agg_val * 10, "m_accum_cost3.tif"), overwrite = TRUE)

accum_cost_2_Carales <- terra::rast(paste0("./Output/Rasters/", agg_val * 10, "m_accum_cost2_carales.tif"))
accum_cost_3_Carales <- terra::rast(paste0("./Output/Rasters/", agg_val * 10, "m_accum_cost3_carales.tif"))

diff_accum_cost_23_Carales <- ((accum_cost_2_Carales-accum_cost_3_Carales) / accum_cost_3_Carales) * 100

terra::writeRaster(x = diff_accum_cost_23_Carales, filename = paste0("./Output/Rasters/", agg_val * 10, "m_diff_accum_cost_23_Carales.tif"), overwrite = TRUE)

####################################
######## EMPIRE ROADS THREE ########
####################################

slope_cs4 <- slope_cs0

slope_cs4 <- leastcostpath::update_values(x = slope_cs4, sf = sard_roads, 
                                          FUN = function(x) { replace(x = x, values = (3.75/3.6)/max(res(sard_dem)))})

slope_cs4_rast <- leastcostpath::rasterise(x = slope_cs4)
terra::writeRaster(x = slope_cs4_rast, filename = paste0("./Output/Rasters/", agg_val * 10, "m_conductance_surface4.tif"), overwrite = TRUE)

accum_cost4a <- create_accum_cost(x = slope_cs4, loc = carales, filename = paste0("./Output/Rasters/", agg_val * 10, "m_accum_cost4_carales.tif"), mask = sard_coast)
accum_cost4b <- create_accum_cost(x = slope_cs4, loc = olbia, filename = paste0("./Output/Rasters/", agg_val * 10, "m_accum_cost4_olbia.tif"), mask = sard_coast)
accum_cost4c <- create_accum_cost(x = slope_cs4, loc = turris_libisonis, filename = paste0("./Output/Rasters/", agg_val * 10, "m_accum_cost4_turris_libisonis.tif"), mask = sard_coast)

accum_cost4_lyrs <- c(accum_cost4a, accum_cost4b, accum_cost4c)
accum_cost4 <- terra::app(accum_cost4_lyrs, fun = function(x) min(x))

terra::writeRaster(x = accum_cost4, filename = paste0("./Output/Rasters/", agg_val * 10, "m_accum_cost4.tif"), overwrite = TRUE)

accum_cost_3_Carales <- terra::rast(paste0("./Output/Rasters/", agg_val * 10, "m_accum_cost3_carales.tif"))
accum_cost_4_Carales <- terra::rast(paste0("./Output/Rasters/", agg_val * 10, "m_accum_cost4_carales.tif"))

diff_accum_cost_34_Carales <- ((accum_cost_3_Carales-accum_cost_4_Carales) / accum_cost_4_Carales) * 100

terra::writeRaster(x = diff_accum_cost_34_Carales, filename = paste0("./Output/Rasters/", agg_val * 10, "m_diff_accum_cost_34_Carales.tif"), overwrite = TRUE)