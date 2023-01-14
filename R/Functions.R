create_accum_cost <- function(x, loc, mask, filename) {
  
  cs_rast <- terra::rast(nrow = x$nrow, ncol = x$ncol, extent = x$extent, crs = x$crs)
  
  from_coords <- sf::st_coordinates(loc)[1, 1:2, drop = FALSE]
  from_cell <- terra::cellFromXY(cs_rast, from_coords)
  cm_graph <- igraph::graph_from_adjacency_matrix(x$conductanceMatrix, mode = "directed", weighted = TRUE)
  igraph::E(cm_graph)$weight <- (1/igraph::E(cm_graph)$weight)
  from_distances <- igraph::distances(cm_graph, v = from_cell,  mode="out")
  accum_rast <- terra::setValues(cs_rast, as.numeric(from_distances))
  
  accum_rast[is.infinite(accum_rast)] <- NA
  
  accum_rast <- terra::mask(x = accum_rast, mask = mask)
  
  # convert from seconds to hours
  accum_rast <- accum_rast / 3600

  terra::writeRaster(x = accum_rast, filename = filename, overwrite = TRUE)  
  
  return(accum_rast)
  
}
