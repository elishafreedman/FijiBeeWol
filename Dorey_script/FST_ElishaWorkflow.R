# This script was started on the 8th of January 2023 to work out FST values for the Fijian Homalictus
  # species for a paper led by Elisha Freedman. This will use only the sequences from 2018 and prior.


#### 0.0 Prep ####
  ##### 0.1 Paths ####
# Choose the path to the root folder in which all other folders can be found (or made by dirMaker)
RootPath <- "/Users/jamesdorey/Desktop/Uni/My_papers/Elisha_Hons_paper/Honours/Dorey_script"
# Set the working directory
setwd(RootPath)
# Install reenv, IF NEEDED
#install.packages("renv")
renv::init() 


##### 0.2 Install packages (if needed) #####
# Install only those packages that are not already present in your system
# Choose packages that need to be installed 
# You may need to install gdal on your computer. This can be done on mac by using
# Homebrew in the terminal and the command "brew install gdal"
list.of.packages <- c("pegas",
                      "apex",
                      "geodata",
                      "dplyr",             #  Part of the tidyverse
                      "adegenet",
                      "mmod",
                      "poppr",
                      "hierfstat",         # For genetics statistics
                      "readr",
                      "rnaturalearth",
                      "rnaturalearthdata",
                      "maps",
                      "terra",
                      "cowplot",
                      "magrittr",          # to use pipes
                      #"ggVennDiagram",     # Extends ggplot2 to make venn diagrams
                      "tibble",            # To use tibbles
                      "forcats",           # tidyverse for working with factors
                      "tidyr",             #  Part of the tidyverse)
                      "tidyselect",        #  Part of the tidyverse
                      "geodata",
                      "tidyterra",
                      "ggspatial")         #  Makes ggplot2 create north arrows or scale bars

# List the new (not installed) packages and then if there are any, install them.
renv::install(packages = c(list.of.packages), 
              rebuild = FALSE) # try changing to TRUE if you're having package troubles
##### 0.3 Load packages ####
# Load all packages from the list specified above
lapply(c(list.of.packages), 
       library, character.only = TRUE)

# Save a snapshot of the environment
renv::snapshot()


#### 1. COI ####
  ##### 1.1 Prep. DNA ####
  # Read in the nexus file
FjHoma <- apex::read.multiFASTA("2018_Fiji_Homalictus.fasta")
  # Data can be visualised here, but what a mess!
# plot(FjHoma)
  # Set locus name because why not
(apex::setLocusNames(FjHoma) <- "COI")
  # Create genind object
FjHoma_genInd <- apex::multidna2genind(FjHoma, mlst = TRUE)


  ##### 1.2 Prep. occurrences + DNA ####
  # Read in the occurrence data
OccData <- readr::read_csv(paste0(RootPath, "/HomalictusCollectionData_2018.csv"),
                           col_types = readr::cols(.default = "c"))
  # Match the species to the collection data
matched <- OccData %>%
  dplyr::right_join(tibble::tibble(labels = FjHoma@labels,
                                    # Set rownumbers to sort by
                                   rownum = row_number(FjHoma@labels)),
                    by = c("Sequence_name" = "labels")) %>%
  # Sort the matched occs
  dplyr::arrange(rownum)

  # Add the species names and codes to the genind object
strata(FjHoma_genInd) <- matched %>% 
    # select wanted columns
  dplyr::select(Specimen_code, Species_name) %>%
    # Make new tibble with this info
  dplyr::rename(
    sequence = Specimen_code,
    populations = Species_name) %>%
  dplyr::mutate(species = populations) %>%
  data.frame()
setPop(FjHoma_genInd) <- ~populations

  ###### a. Filter sample size ####
  # Filter for sample size
Species2remove <- table(FjHoma_genInd@pop) %>% data.frame() %>%
  tibble::tibble() %>%
    # sample size greater than 2
  dplyr::filter(Freq < 3)
  # Choose individuals to remove
individuals2remove <- FjHoma_genInd@strata %>%
  dplyr::right_join(Species2remove, by = c("species" = "Var1"))
  # Remove those individuals
FjHoma_genInd <- FjHoma_genInd[!FjHoma_genInd@strata$sequence %in% individuals2remove$sequence]

##### 1.3 genetic analyses ####
  ###### a. Gst ####
  # Undertake a pairwise Gst Nei analysis
gst <- mmod::pairwise_Gst_Nei(FjHoma_genInd) %>% 
  as.matrix() %>%
  round(3) %>%
  write.csv(paste0(RootPath, "/mmod_Gst_Nei.csv"))

  ###### b. FST ####
  # Undertake pairwise Fst
hierfstat::genet.dist(FjHoma_genInd, method = "Nei87")
Fst_Homa <- hierfstat::pairwise.neifst(dat = FjHoma_genInd, diploid = FALSE) %>%
    # Convert to dataframe
  data.frame() %>%
  round(3) 

colnames(Fst_Homa) <- rownames(Fst_Homa)
  # Save output
write.csv(Fst_Homa, paste0(RootPath, "/hierfstat_FST.csv"))
  

#### 2. Maps ####
  ##### 2.1 prepare data ####
# Read in the occurrence data
OccData <- readr::read_csv(paste0(RootPath, "/HomalictusCollectionData_2018.csv"),
                           col_types = readr::cols(.default = "c")) %>%
  dplyr::rename(
    decimalLongitude = Longitude,
    decimalLatitude = Latitude)
wolbachiaSpecies = readr::read_csv("wolbachiaSpecies.csv")


  ##### 2.2 Raster maps ####
source(paste0(RootPath, "/FjRasterPointMapR.R"))
FjRasterPointMapR(
  mapData = OccData,
  wolSp = wolbachiaSpecies$homaSp,
  #spColours = NULL,
  filename = "rast_FijiHoma_noSpecies.pdf",
  outpath = RootPath,
  width = 7, height = 7, units = "in",
  dpi = 300,
  
  # Map extent — for the raster this is in Fiji 1986 projection (hence the large numbers)
  yLimInput = c(3591385-80000, 4062964+30000),
  xLimInput = c(1871432- 70000, 2298672 + 1.0),
  yLimBreaks = NULL,
  bg = "white", device = "pdf",
  # Inset extent — This is in WGS 84 projection (hence the 360º numbers)
  insetYLim = c(-45, -10), 
  insetXLim = c(145,190),
  #inset size
  insetX = 0.2, 
  insetY = 0.08,
  insetWidth = 0.35, 
  insetHeight = 0.35,
  
  # OPTIONAL:
  pointCol = "#252525",
  wolColour = "#54278f",
  # Raster colours
  naMapCol = "aliceblue",
  rasterGradient = "ggthemes::Blue",
  # 1 or -1 to change directions
  colourDirection = -1,
  # point alpha (opacity)
  mapAlpha = 0.7,
  mapTitle = "Fijian Homalictus occurrences")


##### 2.3 OLD vector maps ####
  ###### a. species-level ####
  # Load in the mapping function
source(paste0(RootPath, "/FjPointMapR.R"))
  # With species coloured
FjPointMapR(
  mapData = OccData %>%
      # Remove NA species
    dplyr::filter(!is.na(Species_name)),
  #spColours = rep("black", length(unique(OccData$Species_name))),
  filename = "FijiHoma_Species.pdf",
  outpath = RootPath,
  width = 16, height = 7, units = "in",
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
  insetX = -0.02, 
  insetY = 0.15,
  insetWidth = 0.3, 
  insetHeight = 0.3,
  
  # OPTIONAL:
  speciesName = NULL,
  nameColumn = "Species_name",
  mapAlpha = 0.5,
  mapTitle = "Fijian Homalictus occurrences",
  # Jitter map? enter jitter amount
  jitterValue = NULL)


###### b. all individuals ####
# WithOUT species coloured
# Load in the mapping function
source(paste0(RootPath, "/FjPointMapR.R"))
FjPointMapR(
  mapData = OccData %>%
    # Remove NA species
    dplyr::filter(!is.na(Species_name)),
  #spColours = rep("black", length(unique(OccData$Species_name))),
  filename = "FijiHoma_noSpecies.pdf",
  outpath = RootPath,
  width = 7, height = 7, units = "in",
  dpi = 300,
  
    # Map extent
  yLimInput = c(-21.5, -16),
  xLimInput = c(180 - 1.0, 180 + 1.0),
  yLimBreaks = NULL,
  xLimBreaks = seq(175, 190, 5),
  xLimLabels = c("175", "180", "-175", "-170"),
  bg = "white", device = "pdf",
    # Inset extent
  insetYLim = c(-45, -10), 
  insetXLim = c(145,190),
  #inset size
  insetX = 0.2, 
  insetY = 0.08,
  insetWidth = 0.35, 
  insetHeight = 0.35,
  
  # OPTIONAL:
  pointCol = "darkgrey",
  speciesName = NULL,
  nameColumn = NULL,
  mapAlpha = 0.5,
  mapTitle = "Fijian Homalictus occurrences",
  # Jitter map? enter jitter amount
  jitterValue = NULL)


  
  


