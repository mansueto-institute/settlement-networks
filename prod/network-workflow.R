# Author: Nicholas Marchio
# Date: February 2022
# R.Version(): R version 4.1.2 (2021-11-01)

print(R.Version() )
# Data packages
library(tidyverse)
library(dplyr)
library(GGally)
library(units)
library(purrr)

# Vector packages
library(sf)
library(lwgeom)
library(stars)
library(starsExtra)
library(nngeo)
library(osmdata)
library(smoothr)

# Network packages
library(tidygraph)
library(sfnetworks)
library(netrankr)
library(dodgr)

# Raster packages
library(raster)
library(rnaturalearth)
library(elevatr)

# Viz packages
library(ggplot2)
library(viridis)
library(patchwork)
library(ggnewscale)

# Load and clean up archaeological data -----------------------------------

load(paste0(projpath,"/graphs/graph_output.RData"))

sites <- read_csv('/Users/nm/Desktop/Projects/work/turk-networks/ancient-networks/ca_sites_table_022322.csv')

sites_clean <- sites %>%
  select_all(~gsub("\\s+|\\.|\\/", "_", .)) %>%
  rename_all(list(tolower)) %>%
  filter(lat != 0 | long != 0) %>% # filter out 0s
  filter(!is.na(lat) | !is.na(long)) %>%
  st_as_sf(coords = c("long", "lat"), 
           crs = 4326, agr = "constant") 

# Parameters --------------------------------------------------------------
projpath = '/Users/nm/Desktop/Projects/work/turk-networks/ancient-networks/outputs'
phase_list = c('_all', '_aneo', '_pneo', '_ech', '_m_lch', '_cha', '_eba', '_ebi_ii', '_ebiii', '_mba', '_lba', '_iimill', '_iron', '_eia_mia', '_lia')

phase_meta = list('_all' = c('8500-300 BCE', 'Neolithic to Iron Age'), '_aneo' = c('8500-7000 BCE', 'Aceramic Neolithic'), '_pneo' = c('7000-6000 BCE', 'Pottery Neolithic'), '_ech' = c('6000-5500 BCE', 'Early Chalcolithic'), '_m_lch' = c('5500-3300 BCE', 'Mid-Late Chalcolithic'),
                  '_ebi_ii' = c('3300-2500 BCE', 'Early Bronze Age I-II'), '_ebiii' = c('2500-2000 BCE', 'Early Bronze Age III'), '_mba' = c('2000-1650 BCE', 'Middle Bronze Age'), '_lba' = c('1650-1200 BCE', 'Late Bronze Age'),
                  '_eia_mia' = c('1200-700 BCE', 'Early/Mid Iron Age'), '_lia' = c('700-300 BCE', 'Late Iron Age'), '_cha' = c('6000-3300 BCE', 'Chalcolithic'), '_eba' = c('3300-2000 BCE', 'Early Bronze Age'), '_iimill' = c('2000-1200 BCE', '2nd millennium'),  '_iron' = c('1200-300 BCE', 'Iron Age'))

sites_data = sites_clean %>% st_transform(3857)

# Functions ---------------------------------------------------------------

generate_network <- function(sites_input, phase_suffix) {
  
  message(paste0(phase_suffix))
  # Preprocess data
  message('Preprocess data')
  sites_input <- sites_data %>% 
    dplyr::select(c("site_id","name","type", "geometry") | ends_with(phase_suffix)) %>%
    rename_with(~str_remove(.,  phase_suffix)) %>%
    filter(period == 1) %>%
    mutate(phase = paste0(phase_suffix))
  
  # Sites of interest
  message('Sites of interest')
  sites_nodes <- sites_input %>% st_transform(3857) #%>% mutate(id = row_number())
  ggplot() + geom_sf(data = sites_nodes %>% st_transform(4326), alpha = .7, fill ='black', color = 'black', aes(size = area)) +
    geom_raster(data = relief, aes(x = x, y = y, alpha = value)) + scale_alpha(range = c(-1.5, 1)) +
    geom_sf(data = water_poly, fill = '#008CCB', alpha = .6, color = alpha('#008CCB',.5)) #+
  #geom_sf(data = ocean_area, fill = '#008CCB', alpha = .6, color = alpha('#008CCB',.5)) + theme_void()
  
  # Build Delaunay triangulation 
  message('Build Delaunay triangulation')
  sites_triangulated <- st_triangulate(sites_nodes %>% sf::st_union(), dTolerance = 0, bOnlyEdges = TRUE) %>%
    st_as_sf() %>% st_cast("LINESTRING") %>% rename(geometry = x) 
  plot(sites_triangulated) 
  
  # Add elevation attributes 
  message('Add elevation attributes')
  options(scipen = 100)
  start_time <- Sys.time()
  sites_triangulated_weights <- sites_triangulated %>% st_transform(3857) %>%
    mutate(st_extract(x = elevation_terrain_slope, at = sites_triangulated %>% st_transform(3857)) %>% st_as_sf() %>%  rename(slope = names(.)[1]) %>% drop_units()) %>% 
    mutate(st_extract(x = elevation_terrain, at = ., FUN = min) %>% st_as_sf() %>% st_drop_geometry() %>% rename(min = names(.)[1] )) %>%
    mutate(st_extract(x = elevation_terrain, at = ., FUN = max) %>% st_as_sf() %>% st_drop_geometry() %>% rename(max = names(.)[1] )) %>%
    #mutate(st_extract(x = elevation_terrain, at = ., FUN = mean) %>% st_as_sf() %>% st_drop_geometry() %>% rename(mean = names(.)[1] )) %>%
    mutate(vertical_distance = max-min, horizontal_distance = st_length(.)) %>% drop_units() %>%
    mutate(equivalent_distance = horizontal_distance + vertical_distance *  7.92 , # Scarf's equivalence of Naismith's rule (meters)
           slope_dydx = (vertical_distance / horizontal_distance),
           slope_radians = slope*(pi/180),
           walking_speed = round((6 * exp(-3.5*abs(tan(slope_radians)+0.05))),2), # Tobler's hiking function (km / h)
           travel_hours = round(((1/walking_speed) * (equivalent_distance/1000)),2))
  ggplot(data = sites_triangulated_weights %>% st_as_sf() %>% st_drop_geometry() %>% 
           mutate(vertical_climb = case_when(vertical_distance <= 200 ~ '1 - 0-200', vertical_distance > 200 & vertical_distance <= 600 ~ '2 - 200-600', vertical_distance > 600 ~ '3 - 600+')), 
         aes(x=travel_hours)) + geom_histogram(bins = 20, aes(fill = vertical_climb))
  plot(sites_triangulated_weights %>% dplyr::select(travel_hours))
  stop_time <- Sys.time()
  print(stop_time - start_time)
  
  # Convert geospatial triangulation data to weighted undirected graph
  message('Convert geospatial triangulation data to weighted undirected graph')
  sites_triangulated_weights_graph <- sites_triangulated_weights %>%
    as_sfnetwork(directed = FALSE) %>%
    activate("edges") %>%
    convert(to_spatial_explicit) %>%
    mutate(weight = travel_hours) %>%
    # Build clusters for later analysis
    activate("nodes") %>%
    # Ward's minimum variance method aims at finding compact, spherical clusters
    mutate(rank = node_rank_hclust(weights = weight, dist = "shortest", method = "ward.D2")) %>%
    # Map equation to find community structure that minimizes expected length of a random walker trajectory
    mutate(clusters = as.factor(group_infomap(weights = weight, node_weights = rank))) 
  
  # Subgraph of sites that are within local clusters 
  message('Subgraph of sites that are within local clusters')
  sites_triangulated_weights_subgraph <- sites_triangulated_weights_graph  %>%
    activate("edges") %>%
    #convert(to_spatial_explicit) %>%
    filter(travel_hours <= 8) %>% # assumption! max 12 hour trip 
    filter(!edge_crosses(water_barriers)) %>%
    activate("nodes") %>%
    filter(!node_is_isolated()) # %>% 
  #activate("edges")
  plot(sites_triangulated_weights_subgraph)
  
  # Subgraph of detours around water bodies
  message('Subgraph of detours around water bodies')
  water_sites <- sites_triangulated_weights %>% 
    st_filter(x = ., y= water_barriers) %>%
    as_sfnetwork(. , directed = FALSE) %>%
    activate('nodes') %>% 
    st_as_sf() 
  plot(water_sites)
  
  # Connect water boundaries to nearest sites
  message('Connect water boundaries to nearest sites')
  water_sites_edges = st_connect(water_sites, water_boundaries, k = 1, 
                                 maxdist = units::set_units(60000,m)) %>% # assumption! only connect sites within 60km of water bodies
    st_union(., water_boundaries) %>%
    st_collection_extract(. , type = c( "LINESTRING")) %>% 
    st_cast("LINESTRING") %>% st_as_sf()
  
  # Add weights to water boundaries / connecting edges
  message('Add weights to water boundaries / connecting edges')
  start_time <- Sys.time()
  water_sites_edges_weights <- water_sites_edges %>% 
    mutate(st_extract(x = elevation_terrain_slope, at = water_sites_edges %>% st_transform(3857)) %>% st_as_sf() %>%  rename(slope = names(.)[1]) %>% drop_units()) %>% 
    mutate(st_extract(x = elevation_terrain, at = ., FUN = min) %>% st_as_sf() %>% st_drop_geometry() %>% rename(min = names(.)[1] )) %>%
    mutate(st_extract(x = elevation_terrain, at = ., FUN = max) %>% st_as_sf() %>% st_drop_geometry() %>% rename(max = names(.)[1] )) %>%
    #mutate(st_extract(x = elevation_terrain, at = ., FUN = mean) %>% st_as_sf() %>% st_drop_geometry() %>% rename(mean = names(.)[1] )) %>%
    mutate(vertical_distance = max-min, horizontal_distance = st_length(.)) %>% drop_units() %>%
    mutate(equivalent_distance = horizontal_distance + vertical_distance *  7.92 , # Scarf's equivalence of Naismith's rule (meters)
           slope_dydx = (vertical_distance / horizontal_distance),
           slope_radians = slope*(pi/180),
           walking_speed = round((6 * exp(-3.5*abs(tan(slope_radians)+0.05))),2), # Tobler's hiking function (km / h)
           travel_hours = round(((1/walking_speed) * (equivalent_distance/1000)),2))
  stop_time <- Sys.time()
  print(stop_time - start_time)
  
  # Convert geospatial water data to weighted undirected graph
  message('Convert geospatial water data to weighted undirected graph')
  water_sites_graph <- water_sites_edges_weights %>% 
    as_sfnetwork(., directed = FALSE) %>%
    activate('edges') %>%
    mutate(weight = travel_hours) %>% 
    activate("nodes") %>%
    filter(!node_is_isolated())
  plot(water_sites_graph)
  
  # Build minimum spanning tree 
  message('Build minimum spanning tree')
  sites_min_span_tree <- sites_triangulated_weights_graph %>%
    convert(to_spatial_explicit) %>%
    activate("edges") %>%
    filter(!edge_intersects(water_barriers %>% st_transform(3857))) %>%
    convert(to_minimum_spanning_tree, weights = weight) %>%
    st_network_join(., water_sites_graph) %>%
    activate("nodes") %>%
    filter(!node_is_isolated())
  plot(sites_min_span_tree)
  
  # Generate centrality measures 
  message('Generate centrality measures')
  end <- sites_min_span_tree %>%
    activate("edges") %>%
    convert(to_spatial_explicit) %>%
    st_network_join(., sites_triangulated_weights_subgraph) %>%
    mutate(weight = case_when(weight == 0 ~ 0.001, TRUE ~ as.numeric(weight))) %>%
    activate("nodes") %>%
    filter(!node_is_isolated()) %>% 
    mutate(eigenvector_centrality_wt = centrality_eigen(weights = weight, directed = FALSE),
           betweeness_centrality_wt = centrality_betweenness(weights = weight, directed = FALSE),
           betweeness_centrality = centrality_betweenness(directed = FALSE),
           eigen_centrality_wt = centrality_eigen(weights = weight, directed = FALSE),
           #pagerank_wt = centrality_pagerank(weights = weight, directed = FALSE),
           degree_centrality_wt = centrality_degree(weights = weight),
           hub_centrality_wt = centrality_hub(weights= weight),
           authority_centrality_wt = centrality_authority(weights = weight)) %>%
    activate('edges') %>%
    mutate(edge_betweeness_centrality_wt = centrality_edge_betweenness(weights = weight, 
                                                                       directed = FALSE))
  plot(end)
  
  # Convert graph to sf dataframes ------------------------------------------
  message('Convert graph to sf dataframes')
  sites_points <- end %>% activate("nodes") %>% st_as_sf() %>% st_transform(4326) 
  sites_points <- sites_points %>% st_join(x = .,y = sites_nodes %>% st_transform(4326), 
                                           join = st_nn, k = 1, maxdist = units::set_units(1,m)) %>%
    filter(!is.na(name))
  sites_routes <- end %>% activate("edges") %>% st_as_sf() %>% st_make_valid() %>% 
    st_transform(4326) %>%
    mutate(phase = phase_suffix)
  
  result <- list("graph" = end, "nodes" = sites_points, "edges" = sites_routes)
  message('Completed')
  return(result)
}


visualize_network <- function(sites_points, sites_routes, 
                              phase_suffix, phase_label, dir_path) {
  
  sites_points <- sites_points %>%
    filter(phase %in% c(phase_suffix)) %>%
    mutate(betweeness_centrality_wt_scale = (betweeness_centrality_wt - min(betweeness_centrality_wt)) / (max(betweeness_centrality_wt) - min(betweeness_centrality_wt)),
           area_scale = (area - min(area)) / (max(area) - min(area)))
  
  sites_routes <- sites_routes %>%
    filter(phase %in% c(phase_suffix)) %>%
    mutate(edge_betweeness_centrality_wt_scale = (edge_betweeness_centrality_wt - min(edge_betweeness_centrality_wt)) / (max(edge_betweeness_centrality_wt) - min(edge_betweeness_centrality_wt)))
  
  theme_map <- theme_minimal() +
    theme(text = element_text(color = "#22211d"), axis.line = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),
          panel.grid.major = element_line(color = "#ebebe5", size = 0.2), panel.grid.minor = element_blank(), plot.background = element_rect(fill = "#f5f5f2", color = NA), legend.background = element_rect(fill = "#f5f5f2", color = NA), panel.border = element_blank())
  
  correlation <- sites_points %>%
    filter(area > 0) %>%
    dplyr::summarise(correl = paste0('Correlation: ',round(cor(area, betweeness_centrality_wt),4), 
                                     ' | P-Value: ',round(cor.test(sites_points$area, sites_points$betweeness_centrality_wt, method="pearson")[[3]],4)),
                     correl2 = paste0('Correlation: ',round(cor.test(sites_points$area, sites_points$betweeness_centrality_wt, method="pearson")[[4]],4))) %>%
    pull(correl)
  
  (p_betweeness<- ggplot() +
      geom_raster(data = relief, aes(x = x, y = y, alpha = value)) +
      scale_alpha(range = c(-1.5, 1)) +
      geom_sf(data = water_line, fill = '#008CCB', alpha = .6, color = alpha('#008CCB',.5)) +
      geom_sf(data = water_poly, fill = '#008CCB', alpha = .6, color = alpha('#008CCB',.5)) +
      geom_sf(data = ocean_area, fill = '#008CCB', alpha = .6, color = alpha('#008CCB',.5)) +
      geom_sf(data = sites_routes, 
              aes(fill = edge_betweeness_centrality_wt_scale,color = edge_betweeness_centrality_wt_scale), 
              size = .6, alpha = .6) +
      geom_sf(data = sites_points %>% filter(!is.na(name)), 
              aes(color =  betweeness_centrality_wt_scale, fill = betweeness_centrality_wt_scale, 
                  size =  area_scale), alpha = .8) +
      scale_color_viridis(option = 'viridis') +
      scale_fill_viridis(option = 'viridis') + 
      theme_map + theme(legend.position = 'none',
                        plot.background = element_blank(),
                        plot.subtitle = element_text(hjust = .5, vjust = -2, size = 13, face = 'bold')) +
      labs(subtitle = 'Betweeness centrality'))  
  
  (p_size <- ggplot() +
      geom_raster(data = relief, aes(x = x, y = y, alpha = value)) +
      scale_alpha(range = c(-1.5, 1)) +
      geom_sf(data = water_line, fill = '#008CCB', alpha = .6, color = alpha('#008CCB',.5)) +
      geom_sf(data = water_poly, fill = '#008CCB', alpha = .6, color = alpha('#008CCB',.5)) +
      geom_sf(data = ocean_area, fill = '#008CCB', alpha = .6, color = alpha('#008CCB',.5)) +
      geom_sf(data = sites_routes, fill = '#FF9671', color = '#FF9671', 
              size = .6, alpha = .6) +
      geom_sf(data = sites_points, aes(size = area_scale, 
                                       color = area_scale, 
                                       fill =  area_scale), alpha = .8) +
      scale_color_viridis(option = 'viridis') +
      scale_fill_viridis(option = 'viridis') + 
      theme_map + theme(legend.position = 'none',
                        plot.background = element_blank(),
                        plot.subtitle = element_text(hjust = .5, vjust = -2, size = 13, face = 'bold')) +
      labs(subtitle = 'Site size'))
  phase_title =  unlist(phase_label, use.names = FALSE)
  title_label = paste0(phase_title[2],': ',phase_title[1])
  (p_size_betweeness <- p_size + p_betweeness + plot_annotation(
    title = sprintf(title_label),
    caption = correlation) & 
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16, )) ) 
  
  ggsave(plot = p_size_betweeness, filename = paste0(dir_path,'/',phase_suffix,'.png'),
         height = 6, width = 10.66, units = "in", dpi = 320, device = 'png')
  
  return(p_size_betweeness)
}

# Download natural data layers ---------------------------------------------

aoi_box = st_bbox(sites_input %>% st_transform(4326)) %>% st_as_sfc()
aoi_box = (aoi_box - st_centroid(aoi_box)) * 1.1 + st_centroid(aoi_box)
aoi_box <- aoi_box %>% st_set_crs(4326)
Sys.sleep(2)

# Download elevation layer for visualization (higher res)
elevation <- get_elev_raster(locations = aoi_box, z = 7, clip = "bbox") # expand = 4000, 
Sys.sleep(2)

# Convert elevation to Spatial Pixel Dataframe
relief <- as.data.frame(as(elevation, "SpatialPixelsDataFrame")) %>%
  rename(value = names(.)[1] )
# Convert elevation to Spatiotemporal Array format
elevation_terrain <- elevation %>%
  projectRaster(from = ., crs=crs(aoi_box %>% st_transform(3857)) ) %>% st_as_stars()
plot(elevation_terrain, breaks = "equal", col = hcl.colors(11, "Spectral"), main = "input (elevation)")
ggplot() +
  geom_stars(data = elevation_terrain) + 
  scale_fill_gradient(low =  '#36454F', high = 'white' , na.value = 'white') 

# Calculate slope in decimal degrees (0-360 clockwise from north)
elevation_terrain_slope = slope(elevation_terrain)
plot(elevation_terrain_slope, breaks = "equal",  col = hcl.colors(11, "Spectral"), main = "output (slope)")

Sys.sleep(2)
# Download OSM water features 
water <- opq(bbox = st_bbox(aoi_box), memsize = 1e+9) %>%
  add_osm_feature(key = 'water') %>%
  osmdata_sf() 
water_mulitpolygons <- water$osm_multipolygons %>% dplyr::select(osm_id)
water_polygons <- water$osm_polygons %>% dplyr::select(osm_id)
water_lines <- water$osm_lines %>% dplyr::select(osm_id)
Sys.sleep(30)
# Download OSM waterway features
waterway <- opq(bbox = st_bbox(aoi_box), memsize = 1e+9) %>%
  add_osm_feature(key = 'waterway') %>%
  osmdata_sf() 
waterway_mulitpolygons <- waterway$osm_multipolygons %>% dplyr::select(osm_id)
waterway_polygons <- waterway$osm_polygons %>% dplyr::select(osm_id)
waterway_lines <- waterway$osm_lines %>% dplyr::select(osm_id)
waterway_multilines <- waterway$osm_multilines %>% dplyr::select(osm_id)
Sys.sleep(30)
# Download OSM coastline features
coastline <- opq(bbox = st_bbox(aoi_box), memsize = 1e+9) %>%
  add_osm_feature(key = 'natural', value = 'coastline') %>%
  osmdata_sf() %>%
  pluck("osm_lines")%>% dplyr::select(osm_id)
# Parse and combine water linestrings and polygons
water_poly <- rbind(water_mulitpolygons,water_polygons,waterway_mulitpolygons,waterway_polygons) %>%
  st_simplify(., dTolerance = 500) %>% 
  st_intersection(.,st_as_sfc(st_bbox(aoi_box))) %>%
  st_transform(4326) %>% dplyr::select(geometry) 
plot(water_poly)
water_line <- rbind(coastline, water_lines,waterway_lines,waterway_multilines) %>%
  #st_simplify(., dTolerance = 1000)  %>% 
  st_intersection(.,st_as_sfc(st_bbox(aoi_box))) %>%
  st_transform(4326) %>% dplyr::select(geometry)
plot(water_line)
Sys.sleep(5)
# Download ocean feature from naturalearthdata.com
# filedir <- paste0(tempdir())
# unlink(filedir, recursive = TRUE)
# dir.create(filedir)
# ocean_shp <- paste0('https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/110m/physical/ne_110m_ocean.zip')
# download.file(url = ocean_shp, destfile = paste0(filedir, basename(ocean_shp)))
# unzip(paste0(filedir,basename(ocean_shp)), exdir= filedir)
# list.files(path = filedir)

#ocean <- ne_download(type = 'ocean', scale = 'large', category = 'physical', returnclass='sf')
ocean <- st_read(paste0(path_dir,'/ne_10m_ocean/ne_10m_ocean.shp'))
ocean_resid <- ocean %>% st_transform(4326) %>% st_make_valid() %>%
  st_intersection(., aoi_box) %>%
  st_cast("POLYGON") %>% st_as_sf()
ocean_area <- st_difference(aoi_box %>% st_union(), ocean_resid %>% st_union()) %>%
  st_as_sf() %>%
  mutate(geometry_type = st_geometry_type(x)) %>%
  filter(geometry_type == 'POLYGON') %>%
  dplyr::select(x) %>% st_as_sf() %>% rename(geometry = x)
plot(ocean_area)

# Create natural barrier polygons -----------------------------------------

# Build bounding polygon
print('Build bounding polygon')
crop_box = st_bbox(aoi_box) %>% st_as_sfc()
crop_box = (crop_box - st_centroid(crop_box)) * 0.999 + st_centroid(crop_box)
crop_box <- crop_box %>% st_set_crs(4326)

# Clean up OSM water features
print('Clean up OSM water features')
water_barriers = rbind(water_poly %>% st_as_sf() ) %>% 
  st_union() %>% st_as_sf() %>% 
  sf::st_crop(x = . , y = crop_box) %>%
  st_transform(3857) %>%
  st_collection_extract(. , type = c( "POLYGON")) %>%
  st_cast("POLYGON") %>%
  mutate(area_w = st_area(x)) %>%
  filter(area_w >= units::set_units(1e+07,m^2)) %>%
  smooth(., method = "chaikin") %>%
  st_make_valid()
plot(water_barriers)

# Polygons of water features (no oceans)
print('Polygons of water features (no oceans)')
water_boundaries_poly <- water_barriers %>% st_transform(3857) %>%
  st_buffer(x = ., dist = units::set_units(300,m)) %>%
  smooth(., method = "chaikin") %>%
  st_buffer(x = ., dist = units::set_units(1000,m)) %>%
  st_simplify(x = ., preserveTopology = TRUE, dTolerance = units::set_units(500,m)) %>%
  st_make_valid()
plot(water_boundaries_poly)

# Linestrings of water features
print('Linestrings of water features')
water_boundaries <- water_boundaries_poly %>%
  st_boundary() %>%
  st_collection_extract(. , type = c( "LINESTRING")) %>%
  st_cast(., "LINESTRING") %>%
  rename(geometry = x) %>%
  dplyr::select(geometry)
plot(water_boundaries)

ggplot() +
  geom_sf(data = water_boundaries_poly, color = 'blue') +
  geom_sf(data = water_barriers, color = 'green') 


# Deploy functions --------------------------------------------------------

out_full <- map(phase_list, function(x) {
  generate_network(sites_input = sites_data, 
                   phase_suffix = x)
})

# Parse and format files --------------------------------------------------

nodes_full <- map_dfr(.x = 1:length(phase_list), ~ out_full[[.x]]$nodes)
edges_full <- map_dfr(.x = 1:length(phase_list), ~ out_full[[.x]]$edges)

nodes_full_df <- nodes_full %>%
  mutate(lon = map_dbl(geometry, ~st_point_on_surface(.x)[[1]]),
         lat = map_dbl(geometry, ~st_point_on_surface(.x)[[2]])) %>%
  st_drop_geometry() 

edges_full_df <- edges_full %>%
  mutate(lon = map_dbl(geometry, ~st_point_on_surface(.x)[[1]]),
         lat = map_dbl(geometry, ~st_point_on_surface(.x)[[2]])) %>%
  st_drop_geometry() 

# Write files
st_write(obj = nodes_full, dsn = paste0(projpath,'/nodes/site_nodes.shp'))
st_write(obj = nodes_full, dsn = paste0(projpath,'/nodes/site_nodes.geojson'))
write_csv(x = nodes_full_df, file = paste0(projpath,'/nodes/site_nodes.csv'))

st_write(obj = edges_full, dsn = paste0(projpath,'/edges/site_edges.shp'))
st_write(obj = edges_full, dsn = paste0(projpath,'/edges/site_edges.geojson'))
write_csv(x = edges_full_df, file = paste0(projpath,'/edges/site_edges.csv'))

save(out_full, file = paste0(projpath,'/graphs/graph_output.RData'))

# Write maps
for (i in seq_along(phase_list)) {
  visualize_network(sites_points = out_full[[i]]$nodes,
                    sites_routes = out_full[[i]]$edges,
                    phase_suffix = phase_list[i],
                    phase_label = phase_meta[i],
                    dir_path = paste0(projpath,'/maps'))
}

