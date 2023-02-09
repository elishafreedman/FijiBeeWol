# This function was written by James B Dorey on the 7th of February 2023
# Its purpose is to visualise Fijian bee occurrences
# Please contact jbdorey@me.com for help


FjRasterPointMapR <- function(
    mapData = NULL,
    wolSp = NULL,
    #spColours = NULL,
    filename = NULL,
    outpath = NULL,
    width = 15, height = 9, units = "in",
    dpi = 300,
    # Map extent
    yLimInput = c(-21.5, -16),
    xLimInput = c(180 - 5, 180 + 6),
    yLimBreaks = NULL,
    bg = "white", device = "pdf",
    size = 1,
    # Inset extent
    insetYLim = c(-45, -10), 
    insetXLim = c(145,190),
    #inset size
    insetX = NULL, 
    insetY = NULL,
    insetWidth = NULL, 
    insetHeight = NULL,
    
    # OPTIONAL:
    pointCol = NULL,
    rasterGradient = terrain.colors(10),
    colourDirection = 1,
    wolColour = NULL,
    speciesName = NULL,
    nameColumn = NULL,
    mapAlpha = 0.5,
    mapTitle = "Fijian Homalictus occurrences",
    # Jitter map? enter jitter amount
    jitterValue = NULL,
    colourPalette = NULL,
    legendTitle = NULL,
    ...
){
  require(ggplot2)
  require(dplyr)
  require(stringr)
  require(tidyr)
  require(tidyselect)
  require(forcats)
  require(sf)
  require(geodata)
  require(tidyterra)
  
  
  #### 0.0 Prep ####
  ##### 0.1 errors ####
  ###### a. FATAL errors ####
  if(is.null(mapData)){
    stop(" — Please provide an argument for mapData I'm a program not a magician.")
  }
  if(is.null(outpath)){
    stop(" — Please provide an argument for outpath it seems reckless to let me just guess.")
  }
  ###### b. warnings ####
  if(is.null(filename)){
    writeLines(" — No argument provided for filename. Using default of 'rastFlagsPlot_DATE.pdf'")
    filename = paste0("rastFlagsPlot_", Sys.Date(),".pdf")
  }
  if(!is.null(speciesName) & is.null(nameColumn)){
    writeLines(" — nameColumn is not provided. Defaulting to scientificName.\n")
    nameColumn = "scientificName"
  }
  
  
  #### 1.0 Data ####
    ##### 1.1 Map data ####
  # Download world map using rnaturalearth packages
  # shift coordinates to recenter worldmap
  worldmap <- ggplot2::map_data("world", wrap = c(0, 360))
    # Download the Fijian DEM rasters — one for each side of the dateline and conver to EPSG:3460
      # — Fiji 1986
  FijiRastWest <- geodata::elevation_3s(country='FJI', 
                                        path = RootPath,
                                        mask = TRUE,
                                        lat = c( -16),
                                        lon = c(180),
                                        res = "0.5") %>%
    terra::project(., terra::crs("EPSG:3460"),
                   gdal = TRUE,   method = "near",
                   threads = TRUE, res = 90)
  FijiRastEast <- geodata::elevation_3s(country='FJI', 
                                        path = RootPath,
                                        mask = TRUE,
                                        lat = c( -16),
                                        lon = c(-180),
                                        res = "0.5") %>%
    terra::project(., terra::crs("EPSG:3460"),
                   gdal = TRUE,   method = "near",
                   threads = TRUE, res = 90)
  # Merge the fiji map halves
  FijiMap <- terra::merge(FijiRastWest, FijiRastEast) 
  
  ##### 1.2 Point data ####
  
  # Create the points file
  fijiPoints <- mapData %>% 
    tidyr::drop_na(c(decimalLongitude, decimalLatitude)) %>%
    dplyr::mutate(decimalLongitude = as.numeric(decimalLongitude),
                  decimalLatitude = as.numeric(decimalLatitude)) %>%
    sf::st_as_sf(coords = c("decimalLongitude", "decimalLatitude")) %>%
    # Set the right CRS
    sf::st_set_crs(terra::crs("EPSG:4326")) %>%
    #terra::vect(geom = c("decimalLongitude", "decimalLatitude")) %>%
    sf::st_transform(., crs = terra::crs("EPSG:3460"))
  
  wolPoints <- fijiPoints %>%
    dplyr::filter(Specimen_code %in% wolSp)
  
  # Reproject the Fiji 1986 extent to WGS84 for the inset
  WGS84extent <- sf::st_bbox(c(xmin = xLimInput[1], xmax = xLimInput[2], 
                               ymin = yLimInput[1], ymax = yLimInput[2]),
                             crs = terra::crs("EPSG:3460")) %>%
    sf::st_as_sfc() %>%
    sf::st_transform(., crs = terra::crs("EPSG:4326")) %>%
    st_bbox()
    
    # Extract the limits and format
  WGS84extent2 <- tibble::tibble(
    point = c("xmin", "xmax", "ymin", "ymax"),
    coords = c(WGS84extent[1], WGS84extent[3],WGS84extent[2], WGS84extent[4]) %>%
      as.numeric()) %>%
    dplyr::mutate(coords = dplyr::if_else(coords < 0 & 
                                            stringr::str_detect(point, "^x"),
                                          coords + 360,
                                          coords)  )
    
  
  #### 2.0 Map ####
  
    ##### 2.1 Rast map ####
  (PointMap <- ggplot() +
    geom_spatraster(data = FijiMap, na.rm = TRUE,
                    aes(fill = srtm_72_16)) +
     # Change map colour scheme
     ggplot2::scale_fill_gradientn(colors = c(paletteer::paletteer_c(rasterGradient, n = 256,
                                                           direction = colourDirection)),
                            na.value = NA, aesthetics ="fill") + 
     labs(fill = "m asl")+ 
    xlim(xLimInput) + ylim(yLimInput) +
    geom_sf(data = fijiPoints,
            #mapping = aes(x = "decimalLongitude", y = "decimalLatitude"),
            alpha = mapAlpha, size = size,
            col = pointCol) +
    geom_sf(data = wolPoints,
            #mapping = aes(x = "decimalLongitude", y = "decimalLatitude"),
            alpha = mapAlpha, size = size,
            col = wolColour) +
    # Map formatting
    # Add in the map's north arrow
    theme(panel.grid.major = element_line(color = gray(.1, alpha = 0.1), 
                                          linetype = "dashed", linewidth = 0.5), # Add grid lines
          panel.border = element_rect(color = gray(.1, alpha = 1), 
                                      linetype = "solid", linewidth = 0.5,
                                      fill = NA), # add panel border
          panel.background = element_rect(fill = "aliceblue") ,
          plot.title = element_text(face = "italic"))+ # Add background — colour in the ocean
    # Add in X and Y labels; 1 or -1
    xlab("Longitude") + ylab("Latitude") + 
    # Add in the title
    ggtitle( mapTitle )
   )

  
  ##### 2.2 Inset map ####
  
  (plotInset <- ggplot(data = worldmap, aes(x = long, y = lat, group = group)) +
     # Plot base polygons
     geom_polygon( fill= NA, colour = "black") +
     coord_sf(ylim = insetYLim, xlim = insetXLim) +
     # Set X and Y limits
     # scale_y_continuous(limits = c(-45, -10)) +
     # scale_x_continuous(limits = c(145,190)
     #                    #, labels = c(170, 175, "180", -175)
     # ) +
     geom_rect(
       aes(xmin = WGS84extent2$coords[1], xmax = WGS84extent2$coords[2], 
           ymin = WGS84extent2$coords[3], ymax = WGS84extent2$coords[4]), fill = 'transparent', 
       col = "red", linewidth = 1
     ) +
     theme(
       panel.background = element_rect(fill='white', color=NA),
       plot.background = element_rect(fill = "transparent",colour = NA),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       panel.border = element_rect(color = gray(.1, alpha = 1), 
                                   linetype = "solid", linewidth = 0.5,
                                   fill = 'transparent'),
       axis.text.x = element_blank(),
       axis.text.y = element_blank()
     ) 
  )
  
  
  ##### 2.3 combine ####
  (outMap <- PointMap %>% 
    cowplot::ggdraw() +
    theme(plot.background = element_rect(fill=NA, color = NA)) +
    cowplot::draw_plot(plotInset,
                       x = insetX, 
                       y = insetY,
                       width = insetWidth, 
                       height = insetHeight)
  )
  
  # save as the map as 10*6"
  ggsave(filename, plot = outMap, device = "pdf", 
         width = width, height = height, dpi = dpi, path = outpath,
         bg = bg)
  return(outMap)
  
  
} # END function