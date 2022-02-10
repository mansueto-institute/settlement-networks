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

# Load and clean up archaeological data -----------------------------------

# Read in file
sites <- st_read('/Users/nm/Desktop/Projects/work/turk-networks/sites/Ca_sites_271021.shp')

# Clean up labels and coding
sites <- sites %>%
  mutate(absolute_date = case_when(LBA %in% c('*','*?','**') ~ '1650-1200 BCE',
                                   TRUE ~ as.character('')),
         archaeological_phase = case_when(LBA == '*' ~ 'Late Bronze Age (LBA)',
                                          LBA == '*?' ~ 'Late Bronze Age (LBA)',
                                          LBA == '**' ~ 'Late Bronze Age (LBA)',
                                          TRUE ~ as.character('')),
         archaeological_horizon = case_when(LBA == '*' ~ '2nd millennium (II_millBCE)',
                                            TRUE ~ as.character('')))
sites_clean <- sites %>%
  dplyr::select(-one_of(c( 'ANeo', 'PNeo', 'ECh', 'MCh_LCh', 'EBI_II', 'EBIII', 'MBA', 'LBA', 'EIA_MIA', 'LIA', 
                    'ANeo_ha', 'PNeo_ha', 'ECh_ha', 'MCh_LCh_ha', 'EBI_II_ha', 'EBIII_ha', 'MBA_ha', 'LBA_ha', 'EIA_MIA_ha', 'LIA_ha'))) 

# Select time window
sites_clean <- sites_clean %>%
  filter(archaeological_phase == 'Late Bronze Age (LBA)')

# Extract lat lons
sites_clean <- sites_clean %>%
  st_transform(3395)  %>%
  mutate(lon = map_dbl(geometry, ~st_point_on_surface(.x)[[1]]),
         lat = map_dbl(geometry, ~st_point_on_surface(.x)[[2]])) %>%
  st_drop_geometry() %>%
  st_as_sf(coords = c("lon", "lat"), 
           crs = 3395, agr = "constant") %>% 
  st_transform(4326)  %>%
  mutate(lon = map_dbl(geometry, ~st_point_on_surface(.x)[[1]]),
         lat = map_dbl(geometry, ~st_point_on_surface(.x)[[2]]))

# Rename columns to generic labels 
sites_clean <- sites_clean %>%
  rename(site_name = Name,
         site_size = Tot_area_h)

# Data assumptions and issues ---------------------------------------------

# Input data format:
# Simpled features dataframe is called "sites_clean"
# Site name column is labelled "site_name"
# Site size column is labelled "site_size"
# All sites are individual Point geometries
# Not tested on areas outside of Turkey so may not perfectly generalize

# Rule based thresholds (see comments marked with "# assumption!"):
# Local networks are routes that take no more than 12 hours to travel by foot
# Routes connectings sites to water bodies must be within 60 km
# Travel time calculations based on Scarf's equivalence of Naismith's rule and Tobler's hiking function

# Features not in current version:
# Add weights for crossing rivers according to Strahler number
# Develop method to handle oceanic networks and sea port nodes

# Known issues:
# Site data is mapped to the centrality measure data via a st_join using the st_nn function
# The sfnetwork adds pseudonodes when routing around circular linestrings (i.e. water bodies)
# See https://github.com/luukvdmeer/sfnetworks/issues/59

# Download natural data layers ---------------------------------------------

aoi_box = st_bbox(sites_clean %>% st_transform(4326)) %>% st_as_sfc()
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
water <- opq(bbox = st_bbox(aoi_box)) %>%
  add_osm_feature(key = 'water') %>%
  osmdata_sf() 
Sys.sleep(5)
water_mulitpolygons <- water$osm_multipolygons %>% dplyr::select(osm_id)
water_polygons <- water$osm_polygons %>% dplyr::select(osm_id)
water_lines <- water$osm_lines %>% dplyr::select(osm_id)
Sys.sleep(5)
# Download OSM waterway features
waterway <- opq(bbox = st_bbox(aoi_box)) %>%
  add_osm_feature(key = 'waterway') %>%
  osmdata_sf() 
Sys.sleep(5)
waterway_mulitpolygons <- waterway$osm_multipolygons %>% dplyr::select(osm_id)
waterway_polygons <- waterway$osm_polygons %>% dplyr::select(osm_id)
waterway_lines <- waterway$osm_lines %>% dplyr::select(osm_id)
waterway_multilines <- waterway$osm_multilines %>% dplyr::select(osm_id)
Sys.sleep(5)
# Download OSM coastline features
coastline <- opq(bbox = st_bbox(aoi_box)) %>%
  add_osm_feature(key = 'natural', value = 'coastline') %>%
  osmdata_sf() %>%
  pluck("osm_lines")%>% dplyr::select(osm_id)
Sys.sleep(5)
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
ocean <- ne_download(scale = 'large', type = 'ocean', category = 'physical', returnclass='sf')
Sys.sleep(5)
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
crop_box = st_bbox(aoi_box) %>% st_as_sfc()
crop_box = (crop_box - st_centroid(crop_box)) * 0.999 + st_centroid(crop_box)
crop_box <- crop_box %>% st_set_crs(4326)

# Clean up OSM water features
water_barriers = rbind(water_poly %>% st_as_sf() ) %>% 
  st_union() %>% st_as_sf() %>% 
  sf::st_crop(x = . , y = crop_box) %>%
  st_transform(3857) %>%
  st_collection_extract(. , type = c( "POLYGON")) %>%
  st_cast("POLYGON") %>%
  mutate(area = st_area(x)) %>%
  filter(area >= units::set_units(1e+07,m^2)) %>%
  smooth(., method = "chaikin") %>%
  st_make_valid()
plot(water_barriers)

# Polygons of water features (no oceans)
water_boundaries_poly <- water_barriers %>% st_transform(3857) %>%
  st_buffer(x = ., dist = units::set_units(300,m)) %>%
  smooth(., method = "chaikin") %>%
  st_buffer(x = ., dist = units::set_units(1000,m)) %>%
  st_simplify(x = ., preserveTopology = TRUE, dTolerance = units::set_units(500,m)) %>%
  st_make_valid()
plot(water_boundaries_poly)

# Linestrings of water features
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

# Generate network and centrality metrics ---------------------------------

# Sites of interest
sites_nodes <- sites_clean %>% st_transform(3857) %>% mutate(id = row_number())
ggplot() + geom_sf(data = sites_nodes %>% st_transform(4326), alpha = .7, fill ='black', color = 'black', aes(size = site_size)) +
  geom_raster(data = relief, aes(x = x, y = y, alpha = value)) + scale_alpha(range = c(-1.5, 1)) +
  geom_sf(data = water_poly, fill = '#008CCB', alpha = .6, color = alpha('#008CCB',.5)) +
  geom_sf(data = ocean_area, fill = '#008CCB', alpha = .6, color = alpha('#008CCB',.5)) + theme_void()

# Build Delaunay triangulation 
sites_triangulated <- st_triangulate(sites_nodes %>% sf::st_union(), dTolerance = 0, bOnlyEdges = TRUE) %>%
  st_as_sf() %>% st_cast("LINESTRING") %>% rename(geometry = x) 
plot(sites_triangulated) 

# Add elevation attributes 
options(scipen = 100)
sites_triangulated_weights <- sites_triangulated %>% st_transform(3857) %>%
  mutate(st_extract(x = elevation_terrain_slope, at = sites_triangulated %>% st_transform(3857)) %>% st_as_sf() %>%  rename(slope = names(.)[1]) %>% drop_units()) %>% 
  mutate(st_extract(x = elevation_terrain, at = ., FUN = var) %>% st_as_sf() %>% st_drop_geometry() %>% rename(sd = names(.)[1] )) %>%
  mutate(st_extract(x = elevation_terrain, at = ., FUN = min) %>% st_as_sf() %>% st_drop_geometry() %>% rename(min = names(.)[1] )) %>%
  mutate(st_extract(x = elevation_terrain, at = ., FUN = max) %>% st_as_sf() %>% st_drop_geometry() %>% rename(max = names(.)[1] )) %>%
  mutate(st_extract(x = elevation_terrain, at = ., FUN = mean) %>% st_as_sf() %>% st_drop_geometry() %>% rename(mean = names(.)[1] )) %>%
  mutate(vertical_distance = max-min, horizontal_distance = st_length(.)) %>% drop_units() %>%
  mutate(equivalent_distance = horizontal_distance + vertical_distance *  7.92 , # Scarf's equivalence of Naismith's rule (meters)
         walking_speed = round((6 * exp(-3.5*abs(tan(slope*(pi/180))+0.05))),2), # Tobler's hiking function (km / h)
         travel_hours = round(((1/walking_speed) * (equivalent_distance/1000)),2))
ggplot(data = sites_triangulated_weights %>% st_as_sf() %>% st_drop_geometry() %>% 
         mutate(vertical_climb = case_when(vertical_distance <= 200 ~ '1 - 0-200', vertical_distance > 200 & vertical_distance <= 600 ~ '2 - 200-600', vertical_distance > 600 ~ '3 - 600+')), 
       aes(x=travel_hours)) + geom_histogram(bins = 20, aes(fill = vertical_climb))

# Convert geospatial triangulation data to weighted undirected graph
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
sites_triangulated_weights_subgraph <- sites_triangulated_weights_graph  %>%
  activate("edges") %>%
  #convert(to_spatial_explicit) %>%
  filter(travel_hours <= 12) %>% # assumption! max 12 hour trip 
  filter(!edge_crosses(water_barriers)) %>%
  activate("nodes") %>%
  filter(!node_is_isolated()) # %>% 
  #activate("edges")
plot(sites_triangulated_weights_subgraph)

# Subgraph of detours around water bodies
water_sites <- sites_triangulated_weights %>% 
  st_filter(x = ., y= water_barriers) %>%
  as_sfnetwork(. , directed = FALSE) %>%
  activate('nodes') %>% 
  st_as_sf() 
plot(water_sites)

# Connect water boundaries to nearest sites
water_sites_edges = st_connect(water_sites, water_boundaries, k = 1, 
                                maxdist = units::set_units(60000,m)) %>% # assumption! only connect sites within 60km of water bodies
  st_union(., water_boundaries) %>%
  st_collection_extract(. , type = c( "LINESTRING")) %>% 
  st_cast("LINESTRING") %>% st_as_sf()

# Add weights to water boundaries / connecting edges
water_sites_edges_weights <- water_sites_edges %>% 
  mutate(st_extract(x = elevation_terrain_slope, at = water_sites_edges %>% st_transform(3857)) %>% st_as_sf() %>%  rename(slope = names(.)[1]) %>% drop_units()) %>% 
  mutate(st_extract(x = elevation_terrain, at = ., FUN = var) %>% st_as_sf() %>% st_drop_geometry() %>% rename(sd = names(.)[1] )) %>%
  mutate(st_extract(x = elevation_terrain, at = ., FUN = min) %>% st_as_sf() %>% st_drop_geometry() %>% rename(min = names(.)[1] )) %>%
  mutate(st_extract(x = elevation_terrain, at = ., FUN = max) %>% st_as_sf() %>% st_drop_geometry() %>% rename(max = names(.)[1] )) %>%
  mutate(st_extract(x = elevation_terrain, at = ., FUN = mean) %>% st_as_sf() %>% st_drop_geometry() %>% rename(mean = names(.)[1] )) %>%
  mutate(vertical_distance = max-min, horizontal_distance = st_length(.)) %>% drop_units() %>%
  mutate(equivalent_distance = horizontal_distance + vertical_distance *  7.92 , # Scarf's equivalence of Naismith's rule (meters)
         walking_speed = round((6 * exp(-3.5*abs(tan(slope*(pi/180))+0.05))),2), # Tobler's hiking function (km / h)
         travel_hours = round(((1/walking_speed) * (equivalent_distance/1000)),2))

# Convert geospatial water data to weighted undirected graph
water_sites_graph <- water_sites_edges_weights %>% 
  as_sfnetwork(., directed = FALSE) %>%
  activate('edges') %>%
  mutate(weight = travel_hours) %>% 
  activate("nodes") %>%
  filter(!node_is_isolated())
plot(water_sites_graph)

# Build minimum spanning tree 
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
end <- sites_min_span_tree %>%
  activate("edges") %>%
  convert(to_spatial_explicit) %>%
  st_network_join(., sites_triangulated_weights_subgraph) %>%
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
  mutate(edge_betweeness_centrality_wt = centrality_edge_betweenness(weights = travel_hours, 
                                                                     directed = FALSE))
plot(end)

# Convert graph to sf dataframes ------------------------------------------

sites_points <- end %>% activate("nodes") %>% st_as_sf() %>% st_transform(4326) 
sites_points <- sites_points %>% st_join(x = .,y = sites_nodes %>% st_transform(4326), 
                                         join = st_nn, k = 1, maxdist = units::set_units(1,m)) %>%
  filter(!is.na(site_name))
sites_routes <- end %>% activate("edges") %>% st_as_sf() %>% st_make_valid() %>% 
  st_transform(4326) 

# Create maps -------------------------------------------------------------

theme_map <- theme_minimal() +
  theme(text = element_text(color = "#22211d"), axis.line = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),
    panel.grid.major = element_line(color = "#ebebe5", size = 0.2), panel.grid.minor = element_blank(), plot.background = element_rect(fill = "#f5f5f2", color = NA), legend.background = element_rect(fill = "#f5f5f2", color = NA), panel.border = element_blank()
    )

(p_betweeness<- ggplot() +
    geom_raster(data = relief, aes(x = x, y = y, alpha = value)) +
    scale_alpha(range = c(-1.5, 1)) +
    geom_sf(data = water_line, fill = '#008CCB', alpha = .6, color = alpha('#008CCB',.5)) +
    geom_sf(data = water_poly, fill = '#008CCB', alpha = .6, color = alpha('#008CCB',.5)) +
    geom_sf(data = ocean_area, fill = '#008CCB', alpha = .6, color = alpha('#008CCB',.5)) +
    geom_sf(data = sites_routes, aes(color = edge_betweeness_centrality_wt, size = edge_betweeness_centrality_wt*.1), alpha = .6) +
    geom_sf(data = sites_points %>% filter(!is.na(site_name)), 
            aes(color =  betweeness_centrality_wt, 
                fill = betweeness_centrality_wt, 
                size =  betweeness_centrality_wt), alpha = .8) +
    scale_color_viridis(option = 'viridis') +
    scale_fill_viridis(option = 'viridis') + 
    theme_map + theme(legend.position = 'none',
                      plot.subtitle = element_text(hjust = .5, vjust = -2, size = 13, face = 'bold')) +
    labs(subtitle = 'Betweeness centrality'))

(p_size <- ggplot() +
    geom_raster(data = relief, aes(x = x, y = y, alpha = value)) +
    scale_alpha(range = c(-1.5, 1)) +
    geom_sf(data = water_line, fill = '#008CCB', alpha = .6, color = alpha('#008CCB',.5)) +
    geom_sf(data = water_poly, fill = '#008CCB', alpha = .6, color = alpha('#008CCB',.5)) +
    geom_sf(data = ocean_area, fill = '#008CCB', alpha = .6, color = alpha('#008CCB',.5)) +
    geom_sf(data = sites_routes, color = '#FF9671', size = .9, alpha = .5) +
    geom_sf(data = sites_points, aes(size = site_size, color = site_size, fill =  site_size), alpha = .8) +
    scale_color_viridis(option = 'viridis') +
    scale_fill_viridis(option = 'viridis') + 
    theme_map + theme(legend.position = 'none',
                      plot.subtitle = element_text(hjust = .5, vjust = -2, size = 13, face = 'bold')) +
    labs(subtitle = 'Site size'))

(p_size_betweeness <- p_size + p_betweeness)

# Correlations ------------------------------------------------------------

sites_corr <- sites_points %>% st_drop_geometry() %>% 
  #filter(site_size > 0, betweeness_centrality_wt > 0) %>%
  dplyr::select(site_size, betweeness_centrality_wt) %>%
  #dplyr::select(site_size, eigenvector_centrality_wt, betweeness_centrality_wt, betweeness_centrality, eigen_centrality_wt, pagerank_wt, degree_centrality_wt, hub_centrality_wt, authority_centrality_wt)
  rename(`Site size` = site_size,
         `Betweeness centrality` = betweeness_centrality_wt)
(correlo <- ggpairs(sites_corr) + 
    theme(text = element_text(size = 20) ))


routes_corr <- sites_routes %>% st_drop_geometry() %>%
  #filter(site_size > 0, betweeness_centrality_wt > 0) %>%
  dplyr::select(edge_betweeness_centrality_wt,
                travel_hours) %>%
  rename(`Travel hours` =  travel_hours,
         `Edge betweeness centrality` = edge_betweeness_centrality_wt)
(correlo2 <- ggpairs(routes_corr) + 
    theme(text = element_text(size = 20) ))


# Export files ------------------------------------------------------------

# Viz
ggsave(plot = p_size_betweeness, '/Users/nm/Desktop/sites.png')
ggsave(plot = correlo, '/Users/nm/Desktop/sites_corr.png')
ggsave(plot = correlo2, '/Users/nm/Desktop/sites_corr2.png')

# Data
write_csv(sites_points, '/Users/nm/Desktop/sites.csv')
st_write(sites_points, '/Users/nm/Desktop/sites.geojson')
st_write(sites_routes, '/Users/nm/Desktop/sites_routes.geojson')

# Appendix ----------------------------------------------------------------

# Cluster correlations
sites_clusters <- sites_points %>% 
  dplyr::group_by(clusters.x) %>% 
  dplyr::summarise_at(vars(site_size, betweeness_centrality_wt), list(mean), na.rm =TRUE) %>%
  st_cast("POLYGON") %>% st_convex_hull() %>%
  st_buffer(dist = units::set_units(10,km))

clusters_corr <- sites_clusters %>% st_drop_geometry() %>%
  dplyr::select(site_size, betweeness_centrality_wt) %>%
  rename(`Site size` = site_size,
         `Betweeness centrality` = betweeness_centrality_wt)
(correlo3 <- ggpairs(clusters_corr) + 
    theme(text = element_text(size = 20) ))

# Build local clusters from graph
sites_clusters_graph <- sites_triangulated_weights_graph %>%
  activate("edges") %>%
  convert(to_spatial_explicit) %>%
  activate("nodes") %>%
  # Ward's minimum variance method aims at finding compact, spherical clusters
  mutate(rank = node_rank_hclust(weights = weight, dist = "shortest", method = "ward.D2")) %>%
  # Map equation to find community structure that minimizes expected length of a random walker trajectory
  mutate(clusters = as.factor(group_infomap(weights = weight, node_weights = rank))) 

sites_clusters <- sites_clusters_graph %>%
  sf::st_as_sf() %>% 
  dplyr::group_by(clusters) %>% dplyr::summarise() %>%
  st_cast("POLYGON") %>% st_convex_hull() %>%
  st_buffer(dist = units::set_units(10,km))
ggplot() + 
  geom_sf(data = sites_clusters, aes(fill = clusters, color = clusters), alpha = .5) + 
  geom_sf(data = sites_clusters_graph %>% activate('nodes') %>% sf::st_as_sf() , color = 'black', alpha = .5) + theme_void()

