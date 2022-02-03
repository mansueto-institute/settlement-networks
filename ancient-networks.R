
# Data packages
library(tidyverse)
library(dplyr)
library(GGally)

# Vector packages
library(sf)
library(stars)
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

sites <- st_read('/Users/nm/Desktop/Projects/work/turk-networks/sites/Ca_sites_271021.shp')

sites <- sites %>%
  mutate(absolute_date = case_when(#ANeo %in% c('*','*?','**') ~ '8500-7000 BCE', 
                                   #PNeo %in% c('*','*?','**') ~ '7000-6000 BCE',
                                   #ECh %in% c('*','*?','**') ~ '6000-5500 BCE',
                                   #MCh_LCh %in% c('*','*?','**') ~ '5500-3300 BCE',
                                   #EBI_II %in% c('*','*?','**') ~ '3300-2500 BCE',
                                   #EBIII %in% c('*','*?','**') ~ '2500-2000 BCE',
                                   #MBA %in% c('*','*?','**') ~ '2000-1650 BCE',
                                   LBA %in% c('*','*?','**') ~ '1650-1200 BCE',
                                   EIA_MIA %in% c('*','*?','**') ~ '1200-700 BCE',
                                   LIA %in% c('*','*?','**') ~ '700-300 BCE',
                                   TRUE ~ as.character('')),
         archaeological_phase = case_when(ANeo == '*' ~ 'Aceramic Neolithic (ANeo)',
                                          PNeo == '*' ~ 'Pottery Neolithic (PNeo)',
                                          ECh == '*' ~ 'Early Chalcolithic (ECh)',
                                          MCh_LCh == '*' ~ 'Mid-Late Chalcolithic (MCh_LCh)',
                                          EBI_II == '*' ~ 'Early Bronze Age I-II (EBI_II)',
                                          EBIII == '*' ~ 'Early Bronze Age III (EBIII)',
                                          MBA == '*' ~ 'Middle Bronze Age (MBA)',
                                          LBA == '*' ~ 'Late Bronze Age (LBA)',
                                          LBA == '*?' ~ 'Late Bronze Age (LBA)',
                                          LBA == '**' ~ 'Late Bronze Age (LBA)',
                                          EIA_MIA == '*' ~ 'Early/Mid Iron Age (EIA_MIA)',
                                          LIA == '*' ~ 'Late Iron Age (LIA)',
                                          TRUE ~ as.character('')),
         archaeological_horizon = case_when(ANeo == '*' ~ 'Neolithic (Neo)',
                                            PNeo == '*' ~ 'Neolithic (Neo)',
                                            ECh == '*' ~ 'Chalcolithic (Cha)',
                                            MCh_LCh == '*' ~ 'Chalcolithic (Cha)',
                                            EBI_II == '*' ~ 'Early Bronze Age (EBA)',
                                            EBIII == '*' ~ 'Early Bronze Age (EBA)',
                                            MBA == '*' ~ '2nd millennium (II_millBCE)',
                                            LBA == '*' ~ '2nd millennium (II_millBCE)',
                                            EIA_MIA == '*' ~ 'Iron Age (Iron)',
                                            LIA == '*' ~ 'Iron Age (Iron)',
                                            TRUE ~ as.character('')))

sites_clean <- sites %>%
  dplyr::select(-one_of(c( 'ANeo', 'PNeo', 'ECh', 'MCh_LCh', 'EBI_II', 'EBIII', 'MBA', 'LBA', 'EIA_MIA', 'LIA', 
                    'ANeo_ha', 'PNeo_ha', 'ECh_ha', 'MCh_LCh_ha', 'EBI_II_ha', 'EBIII_ha', 'MBA_ha', 'LBA_ha', 'EIA_MIA_ha', 'LIA_ha'))) 

sites_clean <- sites_clean %>%
  filter(archaeological_phase == 'Late Bronze Age (LBA)')

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

# Load additional data layers ---------------------------------------------

aoi_box = st_bbox(sites_clean %>% st_transform(4326)) %>% st_as_sfc()
aoi_box = (aoi_box - st_centroid(aoi_box)) * 1.1 + st_centroid(aoi_box)
st_crs(aoi_box) = 4326
# Download elevation layer for visualization (higher res)
elevation <- get_elev_raster(locations = aoi_box, z = 7, clip = "bbox") # expand = 4000, 
Sys.sleep(2)
relief_spdf <- as(elevation, "SpatialPixelsDataFrame")
relief <- as.data.frame(relief_spdf) %>% 
  rename_at(vars(starts_with("file")), ~paste0("value"))
Sys.sleep(2)
# Download elevation layer for vectorization (lower res)
elevation_raster <- get_elev_raster(locations = aoi_box, z = 2, clip = "bbox", expand = .15) # 
elevation_area =  st_contour(x = st_as_stars(elevation_raster), contour_lines = FALSE, breaks = 1300:4000) %>% 
  filter_at(vars(Min, Max), all_vars(!is.infinite(.))) %>%
  st_union() %>% st_as_sf()
plot(elevation_area)
Sys.sleep(2)
# Download impassable elevation layer for vectorization (lower res)
elevation_raster_impassable <- get_elev_raster(locations = aoi_box, z = 5, clip = "bbox", expand = .15) # 
elevation_area_impassable =  st_contour(x = st_as_stars(elevation_raster_impassable), contour_lines = FALSE, breaks = 1900:4000) %>% 
  filter_at(vars(Min, Max), all_vars(!is.infinite(.))) %>%
  st_union() %>% st_as_sf() %>%
  rename(geometry = x)
plot(elevation_area_impassable)
Sys.sleep(2)
# Download water
water <- opq(bbox = st_bbox(aoi_box)) %>%
  add_osm_feature(key = 'water') %>%
  osmdata_sf() 
water_mulitpolygons <- water$osm_multipolygons %>% dplyr::select(osm_id)
water_polygons <- water$osm_polygons %>% dplyr::select(osm_id)
water_lines <- water$osm_lines %>% dplyr::select(osm_id)
Sys.sleep(2)
# Download waterways
waterway <- opq(bbox = st_bbox(aoi_box)) %>%
  add_osm_feature(key = 'waterway') %>%
  osmdata_sf() 
waterway_mulitpolygons <- waterway$osm_multipolygons %>% dplyr::select(osm_id)
waterway_polygons <- waterway$osm_polygons %>% dplyr::select(osm_id)
waterway_lines <- waterway$osm_lines %>% dplyr::select(osm_id)
waterway_multilines <- waterway$osm_multilines %>% dplyr::select(osm_id)
Sys.sleep(5)
# Download coastline
coastline <- opq(bbox = st_bbox(aoi_box)) %>%
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
  st_simplify(., dTolerance = 1000)  %>% 
  st_intersection(.,st_as_sfc(st_bbox(aoi_box))) %>%
  st_transform(4326) %>% dplyr::select(geometry)
plot(water_line)
# Download ocean
ocean <- ne_download(scale = 'large', type = 'ocean', category = 'physical', returnclass='sf')
Sys.sleep(2)
ocean_resid <- ocean %>% st_transform(4326) %>% st_make_valid() %>%
  st_intersection(., aoi_box) %>%
  st_cast("POLYGON") %>%
  st_as_sf() 
ocean_area <- st_difference(aoi_box %>% st_union(), ocean_resid %>% st_union()) %>%
  st_as_sf() %>%
  mutate(geometry_type = st_geometry_type(x)) %>%
  filter(geometry_type == 'POLYGON') %>%
  dplyr::select(x) %>%
  st_as_sf() %>%
  rename(geometry = x)
plot(ocean_area)

# Create natural barrier polygons -----------------------------------------

# Build bounding polygon
crop_box = st_bbox(aoi_box) %>% st_as_sfc()
crop_box = (crop_box - st_centroid(crop_box)) * 0.999 + st_centroid(crop_box)
st_crs(crop_box) = 4326

impassable_barriers = rbind(water_poly %>% st_as_sf(), ocean_area %>% st_as_sf(), elevation_area_impassable  %>% st_as_sf() ) %>%
  st_union() %>% st_as_sf() %>% 
  sf::st_crop(x = . , y = crop_box) %>%
  st_transform(3857) %>%
  st_collection_extract(. , type = c( "POLYGON")) %>%
  st_cast("POLYGON") %>%
  mutate(area = st_area(x)) %>%
  filter(area >= units::set_units(1e+08,m^2))
plot(impassable_barriers)

# Build network and generate centrality metrics ---------------------------

# Build Delaunay triangulation across all sites
sites_graph <- sites_clean %>% st_transform(3857)
sites_triangulated <- st_triangulate(sites_graph %>% sf::st_union(), dTolerance = 0, bOnlyEdges = TRUE) %>%
  st_as_sf() %>% st_cast("LINESTRING") %>% rename(geometry = x) 
plot(sites_triangulated) 

# Build local site clusters
sites_community_graph <- as_sfnetwork(sites_triangulated, directed = FALSE) %>%
  activate("edges") %>%
  convert(to_spatial_explicit) %>%
  mutate(weight = edge_length()) %>% 
  activate("nodes") %>%
  mutate(rank = node_rank_hclust(weights = weight)) %>%
  mutate(community = as.factor(group_infomap(weights = weight, node_weights = rank))) 

# Convert clusters to polygons
sites_community <- sites_community_graph %>%
  sf::st_as_sf() %>% 
  dplyr::group_by(community) %>% 
  dplyr::summarise() %>%
  st_cast("POLYGON") %>%
  st_convex_hull() %>%
  st_buffer(dist = 20)
plot(sites_community)

# Connect sites that are within local clusters and no further than 50km apart
sites_tri_graph <- as_sfnetwork(sites_triangulated, directed = FALSE) %>%
  activate("edges") %>%
  convert(to_spatial_explicit) %>%
  mutate(weight = edge_length()) %>% 
  filter(weight <= units::set_units(60,km)) %>%
  filter(edge_is_within(sites_community)) %>%
  filter(!edge_crosses(impassable_barriers %>% st_transform(3857)))
plot(sites_tri_graph)
sites_tri_graph 

# Find isolated nodes
sites_cluster_graph <- sites_tri_graph %>%
  activate("nodes") %>%
  filter(!node_is_isolated()) %>% 
  activate("edges") #%>% 
#st_as_sf()
plot(sites_cluster_graph )
sites_cluster_graph 

# Build minimum spanning tree that doesn't cross impassable barriers
sites_min_span_tree <- as_sfnetwork(sites_triangulated, directed = FALSE) %>%
  convert(to_spatial_explicit) %>%
  activate("edges") %>%
  filter(!edge_intersects(impassable_barriers %>% st_transform(3857))) %>%
  mutate(weight = edge_length()) %>%
  convert(to_minimum_spanning_tree, weights = weight) 
plot(sites_min_span_tree)
sites_min_span_tree

# Manually join isolated nodes to nearest node
sites_iso_nodes <- sites_min_span_tree %>%
  activate("nodes") %>%
  filter(node_is_isolated()) 
if (nrow(sites_iso_nodes %>% activate("nodes") %>% st_as_sf()) > 0) {
  sites_iso_edges <- st_connect(x = sites_iso_nodes %>% st_as_sf(), 
                                y = sites_min_span_tree %>% st_as_sf(), k = 1) %>% as_sfnetwork()
  sites_min_span_tree <- sites_min_span_tree %>% 
    activate("nodes") %>% 
    filter(!node_is_isolated()) %>% 
    activate("edges") %>%
    st_network_join(., sites_iso_edges)
}

# Generate centrality measures for nodes
end <- sites_min_span_tree %>%
  st_network_join(., sites_cluster_graph) %>%
  as_sfnetwork(directed = FALSE) %>%
  activate("edges") %>%
  convert(to_spatial_explicit) %>%
  mutate(weight = edge_length()) %>% 
  activate("nodes") %>%
  filter(!node_is_isolated()) %>%
  mutate(eigenvector_centrality_wt = centrality_eigen(weights = weight),
         betweeness_centrality_wt = centrality_betweenness(weights = weight),
         betweeness_centrality = centrality_betweenness(),
         pagerank_wt = centrality_pagerank(weights = weight),
         degree_centrality_wt = centrality_degree(weights = weight),
         hub_centrality_wt = centrality_hub(weights= weight),
         authority_centrality_wt = centrality_authority(weights = weight),
         alpha_centrality = centrality_alpha(),
         integration_centrality = centrality_integration(),
         communicability_centrality = centrality_communicability()) 
plot(end)
end

#random_walk = centrality_random_walk(),
#closeness = centrality_closeness(weights = weight),
#information = centrality_information()

# Visualize the results
theme_map <- theme_minimal() +
  theme(
    text = element_text(color = "#22211d"),
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_line(color = "#ebebe5", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "#f5f5f2", color = NA), 
    panel.background = element_rect(fill = "#f5f5f2", color = NA), 
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.border = element_blank()
  )

(p_betweeness<- ggplot() +
    geom_raster(data = relief, aes(x = x, y = y, alpha = value)) +
    scale_alpha(range = c(-1.5, 1)) +
    geom_sf(data = water_line, fill = '#008CCB', alpha = .6, color = alpha('#008CCB',.5)) +
    geom_sf(data = water_poly, fill = '#008CCB', alpha = .6, color = alpha('#008CCB',.5)) +
    geom_sf(data = ocean_area, fill = '#008CCB', alpha = .6, color = alpha('#008CCB',.5)) +
    geom_sf(data = activate(end, "edges") %>% st_as_sf(), color = '#FF9671', size = .9, alpha = .5) +
    geom_sf(data = activate(end, "nodes") %>% st_as_sf(), 
            aes(color =  betweeness_centrality, fill = betweeness_centrality, size =  betweeness_centrality), 
            alpha = .8) +
    scale_color_viridis(option = 'plasma') +
    scale_fill_viridis(option = 'plasma') + 
    theme_map + theme(legend.position = 'none',
                      plot.subtitle = element_text(hjust = .5, vjust = -2, size = 13, face = 'bold')) +
    labs(subtitle = 'Betweeness centrality'))

(p_size <- ggplot() +
    geom_raster(data = relief, aes(x = x, y = y, alpha = value)) +
    scale_alpha(range = c(-1.5, 1)) +
    geom_sf(data = water_int, fill = '#008CCB', alpha = .6, color = alpha('#008CCB',.5)) +
    geom_sf(data = ocean_area, fill = '#008CCB', alpha = .6, color = alpha('#008CCB',.5)) +
    geom_sf(data = activate(end, "edges") %>% st_as_sf(), color = '#FF9671', size = .9, alpha = .5) +
    geom_sf(data = sites_clean, aes(size = Tot_area_h, color = Tot_area_h, fill =  Tot_area_h), alpha = .8) +
    scale_color_viridis(option = 'plasma') +
    scale_fill_viridis(option = 'plasma') + 
    theme_map + theme(legend.position = 'none',
                      plot.subtitle = element_text(hjust = .5, vjust = -2, size = 13, face = 'bold')) +
    labs(subtitle = 'Site size'))

(p_size_betweeness <- p_size + p_betweeness)

ggsave(plot = p_size_betweeness, '/Users/nm/Desktop/sites.png')

# Export files ------------------------------------------------------------

sites_points <- end %>% activate("nodes") %>% st_as_sf() %>% st_transform(4326) 
sites_points <- sites_points %>% st_join(x = .,y = sites_clean %>% st_transform(4326), join = st_nn)

sites_routes <- end %>%  activate("edges") %>% st_as_sf() %>% st_make_valid() %>% 
  st_transform(4326) 

write_csv(sites_points, '/Users/nm/Desktop/sites.csv')
st_write(sites_points, '/Users/nm/Desktop/sites.geojson')
st_write(sites_routes, '/Users/nm/Desktop/sites_routes.geojson')

# Correlation visualization
sites_corr <- sites_points %>% st_drop_geometry() %>% 
  filter(Tot_area_h > 0, betweeness_centrality > 0) %>%
  mutate(Tot_area_h_ln = log(Tot_area_h),
         betweeness_centrality_ln = log(betweeness_centrality)) %>%
  dplyr::select(Tot_area_h, betweeness_centrality) %>%
  rename(`Site size` = Tot_area_h,
         `Betweeness centrality` = betweeness_centrality)

(correlo <- ggpairs(sites_corr) )
ggsave(plot = correlo, '/Users/nm/Desktop/sites_corr.png')

