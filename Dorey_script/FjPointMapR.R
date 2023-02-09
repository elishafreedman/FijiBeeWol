# This function was written by James B Dorey on the 9th of January 2023
# Its purpose is to visualise Fijian bee occurrences
# Please contact jbdorey@me.com for help


FjPointMapR <- function(
    mapData = NULL,
    #spColours = NULL,
    filename = NULL,
    outpath = NULL,
    width = 15, height = 9, units = "in",
    dpi = 300,
    # Map extent
    yLimInput = c(-21.5, -16),
    xLimInput = c(180 - 5, 180 + 6),
    yLimBreaks = NULL,
    xLimBreaks = seq(175, 190, 5),
    xLimLabels = c("175", "180", "-175", "-170"),
    bg = "white", device = "pdf",
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
    writeLines(" — No argument provided for filename. Using default of 'FlagsPlot_DATE.pdf'")
    filename = paste0("FlagsPlot_", Sys.Date(),".pdf")
  }
  if(!is.null(speciesName) & is.null(nameColumn)){
    writeLines(" — nameColumn is not provided. Defaulting to scientificName.\n")
    nameColumn = "scientificName"
  }
  

#### 1.0 Data ####
  # Download world map using rnaturalearth packages
  # shift coordinates to recenter worldmap
  worldmap <- ggplot2::map_data("world", wrap = c(0, 360))

    # Recenter points
  center <- 180
  # shift coordinates to recenter x
  mapData$decimalLongitude_reCenter <- dplyr::if_else(
    mapData$decimalLongitude %>% as.numeric() < center - 180, 
                                  mapData$decimalLongitude %>% as.numeric() + 360, 
                                  mapData$decimalLongitude %>% as.numeric())
  mapData$decimalLatitude <- mapData$decimalLatitude  %>% as.numeric()

##### 1.1 map ####
  # Create the checklist map
  (PointMap <- ggplot(data = worldmap, aes(x = long, y = lat)) +
      # Plot base polygons
     geom_polygon(aes(group = group), fill="white", colour = "black") +
      # CORE plotting of map and data
      # Set X and Y limits
     scale_y_continuous(limits = yLimInput) +
     scale_x_continuous(limits = xLimInput,
                        breaks = xLimBreaks,
                        labels = xLimLabels
                        #, labels = c(170, 175, "180", -175)
                        ) +
      # Set 1:1 ratio
     coord_equal() +
      # plot point data
      # POINTS IF IS NULL; i.e. DON'T jitter
     ###### a. noJitter + colour ####
      {if(is.null(jitterValue) & !is.null(nameColumn))
        geom_point(data = mapData,
                   mapping = aes(x = decimalLongitude_reCenter, y = decimalLatitude,
                                 col = .data[[nameColumn]]),
                   alpha = mapAlpha)} +
     ###### b. noJitter + setColour ####
   {if(is.null(jitterValue) & !is.null(pointCol))
     geom_point(data = mapData,
                mapping = aes(x = decimalLongitude_reCenter, y = decimalLatitude),
                col = pointCol,
                alpha = mapAlpha)} +
     ###### c. Jitter + colour ####
      # POINTS IF IS NOT NULL; i.e. jitter
      {if(!is.null(jitterValue) & !is.null(nameColumn))
        geom_jitter(mapData, 
                 mapping = aes(x = decimalLongitude_reCenter, y = decimalLatitude,
                               colour = .data[[nameColumn]]),
                 alpha = mapAlpha, width = jitterValue, height = jitterValue)}+ 
     ###### d. Jitter + setColour ####
   # POINTS IF IS NOT NULL; i.e. jitter
   {if(!is.null(jitterValue) & !is.null(pointCol))
     geom_jitter(mapData, 
                 mapping = aes(x = decimalLongitude_reCenter, y = decimalLatitude,
                               col = pointCol),
                 alpha = mapAlpha, width = jitterValue, height = jitterValue)}+ 
     
     
          # Colours
     {if(!is.null(colourPalette))
       scale_color_manual(values = ,
                            name = legendTitle)} +
      # Map formatting
      # Add in the map's north arrow
      theme(panel.grid.major = element_line(color = gray(.1, alpha = 0.1), 
                                            linetype = "dashed", linewidth = 0.5), # Add grid lines
            panel.border = element_rect(color = gray(.1, alpha = 1), 
                                        linetype = "solid", linewidth = 0.5,
                                        fill = NA), # add panel border
            panel.background = element_rect(fill = "aliceblue") ,
            plot.title = element_text(face = "italic"))+ # Add background — colour in the ocean
      # Change map colour scheme
      scale_fill_viridis_d(option = "magma") + # options = "magma", "inferno", "plasma", "cividis"
      # Add in X and Y labels
      xlab("Longitude") + ylab("Latitude") + 
      # Add in the title
      ggtitle( mapTitle) )
  
    #### 1.2 inset ####
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
       aes(xmin = xLimInput[1], xmax = xLimInput[2], 
           ymin = yLimInput[1], ymax = yLimInput[2]), fill = 'transparent', 
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
  
  ##### 1.3 combine ####
  outMap <- PointMap %>% 
    cowplot::ggdraw() +
    theme(plot.background = element_rect(fill=NA, color = NA)) +
    cowplot::draw_plot(plotInset,
                       x = insetX, 
                       y = insetY,
                       width = insetWidth, 
                       height = insetHeight)
  
  # save as the map as 10*6"
  ggsave(filename, plot = outMap, device = "pdf", 
         width = width, height = height, dpi = dpi, path = outpath,
         bg = bg)
return(outMap)


}# END function




