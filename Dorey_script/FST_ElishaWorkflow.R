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
#

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
FjHoma_genInd <- apex::multidna2genind(FjHoma, mlst = FALSE,
                                       genes = "COI")


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

#### 1.4 Diversity Indexces ####
  # Read in the sheet with the infected/not infected individuals (all tested)
wolbachiaInfected = readr::read_csv("Wolbachia_PositiveNagative.csv")


  ###### a. calculate ####
# Use the genetic information in the matched dataframe to get haplotype statistics
FJHoma_haplotypes <- matched %>%
    # !!! OPTIONAL filter to ONLY Wolbachia individuals
  #dplyr::filter(Specimen_code %in% wolbachiaInfected$seqCode) %>%
    # Group by species name
   # !!! OPTIONAL remove fijiensis
  #dplyr::filter(!Species_name == "Lasioglossum (Homalictus) fijiensis") %>%
  dplyr::group_by(Species_name) %>%
    # Get counts of each haplotype
  dplyr::count(Sequence) %>% 
    # remove the Sequence column as it's not needed for the stats
  dplyr::select(!Sequence) %>%
    # Add row numbers to make them unique
  mutate(row = row_number()) %>%
    # Pivot the tible wider so that each species has it's own column with haplotype counts
  tidyr::pivot_wider(names_from = Species_name,
                     values_from = n,
                     values_fill = 0) 
  
  # Set up formulae from ShannonGen as functions
    # Zahl_1977
Z = function(X) {
  X = X[X > 0]
  Y = X[X > 1]
  n = sum(X)
  -n * sum(X/n * log(X/n)) - (n - 1)/n * sum((n - X) * 
        (-X/(n - 1) * log(X/(n - 1)))) - (n - 1)/n * sum(-Y * 
        (Y - 1)/(n - 1) * log((Y - 1)/(n - 1)))
}

  # Shannon_1949
MLE = function(X) {
  X = X[X > 0]
  n = sum(X)
  -sum(X/n * log(X/n))
}

  # Calculate the statistics 
out_Zahl_1977 <- apply(FJHoma_haplotypes, MARGIN = 2, FUN = Z)
out_Shannon_1949 <- apply(FJHoma_haplotypes, MARGIN = 2, FUN = MLE)

source("ChaoWrapper.R")
# Remove the row column
FJHoma_haplotypes <- FJHoma_haplotypes %>%
  dplyr::select(!row)
# Calculate diversity indices from ChaoSpecies
ChaoResults <- ChaoWrapper(data = FJHoma_haplotypes)
ChaoResults$basicTable %>% readr::write_csv("basicChaoOutputs.csv")
ChaoResults$diversityTable %>% readr::write_csv("diversityChaoOutputs.csv")
  # Get the estimates of Chao1
Chao1 <- ChaoResults$diversityTable %>% dplyr::filter(Name == "Chao1 (Chao, 1984)")




  # Combine the statistics
outCombined <- dplyr::bind_cols(names(out_Zahl_1977), out_Zahl_1977, out_Shannon_1949) %>%
    # Set the column names
  setNames(c("Species_name", "Zahl", "Shannon")) %>%
    # Remove the "row" statistic 
  dplyr::filter(!Species_name == "row") %>%
    # Add haplotype counts
  dplyr::left_join(matched %>%
                     # Group by species name
                     dplyr::group_by(Species_name) %>%
                     dplyr::distinct(Sequence, .keep_all = TRUE) %>%
                     # Get counts of each haplotype
                     dplyr::count(., name = "haplotypeCount"),
                   by = "Species_name") %>%
   # Add sequence counts
  dplyr::left_join(matched %>%
                     # Group by species name
                     dplyr::group_by(Species_name) %>%
                     # Get counts of each haplotype
                     dplyr::count(name = "sequenceCount"),
                   by = "Species_name") %>%
    # Add Wolbachia infection status
  dplyr::mutate(WolbachiaDetected = dplyr::if_else(
    Species_name %in% c("Lasioglossum (Homalictus) ostridorsum", "Lasioglossum (Homalictus) kaicolo",  "Lasioglossum (Homalictus) hadrander", 
                        "Lasioglossum (Homalictus) groomi",  "Lasioglossum (Homalictus) fijiensis",   "Lasioglossum (Homalictus) concavus",  
                        "Lasioglossum (Homalictus) atritergus",  "Lasioglossum (Homalictus) sp. S",   "Lasioglossum (Homalictus) sp. F",  
                        "Lasioglossum (Homalictus) sp. M",  "Lasioglossum (Homalictus) sp. R" ),
    "Infected", "Unknown")) %>%
  dplyr::left_join(Chao1 %>% 
                     dplyr::select(species, Estimate) %>%
                     dplyr::rename(ChaoEstimate = Estimate) %>%
                     dplyr::mutate(log_ChaoEstimate = log(ChaoEstimate)),
                   by = c("Species_name" = "species"))

outCombined_longer <- outCombined %>%
  dplyr::select(!c("haplotypeCount","sequenceCount")) %>%
  tidyr::pivot_longer(cols = c("Zahl","Shannon", "log_ChaoEstimate")) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::arrange(WolbachiaDetected) %>% 
  dplyr::group_by(WolbachiaDetected) %>%
    # Change 0 to NA
  dplyr::mutate(value = dplyr::if_else(value == 0,
                                       NA_integer_, value)) %>%
    # Drop na values for Shannon and Zahl
  tidyr::drop_na(value)
  
  ###### b. diversity plots ####
(diversityPlot <- ggplot2::ggplot(outCombined_longer, aes(fill=name, y=value, 
                                        x= reorder(Species_name, value, decreasing = TRUE))) + 
  ggplot2::geom_bar(position="dodge", stat="identity") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggplot2::xlab("Species") + ggplot2::ylab("Statistic value") +
  ggplot2::facet_wrap(~WolbachiaDetected, drop = FALSE, scale="free_x"))

  ggplot2::ggsave("statistics_perSpecies.pdf", plot = diversityPlot, 
                  width = 10, height = 8, units = "in", dpi = 300)

  ###### c. shapiro.test ####
  # Remove NA from outCombined
outCombined_complete <- outCombined %>%
  dplyr::filter(!c(is.na(Zahl) | is.na(Shannon)| Zahl == 0 | Shannon == 0 ) )
  readr::write_csv(outCombined_complete, "outCombined_complete.csv")

# Test normality
shapiro.test(outCombined_complete$Zahl)
shapiro.test(outCombined_complete$Shannon)
shapiro.test(outCombined_complete$log_ChaoEstimate)
shapiro.test(outCombined_complete$haplotypeCount)
shapiro.test(outCombined_complete$sequenceCount)

  # Get subsets of Unknown and Infected data
UnknownDF <- outCombined_complete %>%
  dplyr::filter(WolbachiaDetected == "Unknown")
InfectedDF <- outCombined_complete %>%
  dplyr::filter(WolbachiaDetected == "Infected")

  # Calculate p-values
    # Wilcoxon rank sum test (equivalent to the Mann-Whitney test: see the Note) is carried out
ZahlP <- wilcox.test(UnknownDF$Zahl, InfectedDF$Zahl, alternative = "two.sided", paired = FALSE,
                     exact = FALSE, conf.int = TRUE)
ShannonP <- wilcox.test(UnknownDF$Shannon, InfectedDF$Shannon, alternative = "two.sided", paired = FALSE,
                   exact = FALSE)
ChaoEstimateP <- wilcox.test(UnknownDF$log_ChaoEstimate, InfectedDF$log_ChaoEstimate,
                                 alternative = "two.sided", paired = FALSE,
                                 exact = FALSE)
sequenceCountP <- wilcox.test(UnknownDF$sequenceCount, InfectedDF$sequenceCount, alternative = "two.sided", paired = FALSE,
                         exact = FALSE)
haplotypeCountP <- wilcox.test(UnknownDF$haplotypeCount, InfectedDF$haplotypeCount, alternative = "two.sided", paired = FALSE,
                          exact = FALSE)

  # Combine t-test outputs into a table
W_testOutput <- tibble::tibble(
  name = c("Zahl", "Shannon", "ChaoEstimate", "sequenceCount", "haplotypeCount"),
  #diffInLocation = c(ZahlP$estimate[1], ShannonP$estimate[1], sequenceCountP$estimate[1], haplotypeCountP$estimate[1]),
  W_statistic = c(ZahlP$statistic, ShannonP$statistic, ChaoEstimateP$statistic,
                  sequenceCountP$statistic, haplotypeCountP$statistic),
  #df = c(ZahlP$parameter, ShannonP$parameter, sequenceCountP$parameter, haplotypeCountP$parameter),
  p_value = c(ZahlP$p.value, ShannonP$p.value, ChaoEstimateP$p.value,
              sequenceCountP$p.value, haplotypeCountP$p.value)
  )

readr::write_csv(W_testOutput, "W_testOutput.csv")

  ###### d. mean plots ####
  # make the data sets for the statistics and for the sampling
outCombined_plot_Stats <- outCombined_complete %>% 
  dplyr::mutate(rowNum = row_number()) %>%
  tidyr::pivot_longer(
    cols =  c("Zahl", "Shannon", "log_ChaoEstimate", "haplotypeCount", "sequenceCount")) %>%
  dplyr::filter(name %in% c("Zahl", "Shannon", "log_ChaoEstimate"))
  # sampling
outCombined_plot_Sampling <- outCombined_complete %>% 
  dplyr::mutate(rowNum = row_number()) %>%
  tidyr::pivot_longer(
    cols =  c("Zahl", "Shannon", "log_ChaoEstimate", "haplotypeCount", "sequenceCount")) %>%
  dplyr::filter(name %in% c("haplotypeCount", "sequenceCount"))

  # Statistic plot
statPlot <- ggplot2::ggplot(outCombined_plot_Stats, 
                            aes(x= name, y=value, fill=WolbachiaDetected)) + 
  ggplot2::geom_boxplot() +
  ggplot2::xlab("") + ggplot2::ylab("Diversity/richness value") +
  ggplot2::theme(legend.position = "none",
                 panel.background = ggplot2::element_rect(fill = "transparent",
                                                          colour = "black",
                                                          linetype = NULL)) +
  ggplot2::scale_fill_manual(values = c("#E97777", "#82AAE3")) +
  ggplot2::scale_x_discrete(limits = c("Zahl", "Shannon", "log_ChaoEstimate"),
                            labels = c("Zahl", "Shannon", "log(Chao)"))
# Sampling plot
samplePlot <- ggplot2::ggplot(outCombined_plot_Sampling, aes(x= name, y= log(value), fill=WolbachiaDetected)) + 
  ggplot2::geom_boxplot() +
  ggplot2::xlab("") + ggplot2::ylab("Log of count") +
  ggplot2::theme(legend.position = "none",
                 panel.background = ggplot2::element_rect(fill = "transparent",
                                                          colour = "black",
                                                          linetype = NULL))+
  ggplot2::scale_fill_manual(values = c("#E97777", "#82AAE3"))
  # Legend
legendPlot <- ggplot2::ggplot(outCombined_complete %>% 
                  dplyr::mutate(rowNum = row_number()) %>%
                  tidyr::pivot_longer(
                    cols =  c("log_ChaoEstimate", "Zahl", "Shannon", "haplotypeCount", "sequenceCount")), 
                aes(x= name, y= value, fill=WolbachiaDetected)) + 
  ggplot2::geom_boxplot() +
  ggplot2::xlab("") + ggplot2::ylab("Log of count") +
  ggplot2::theme(legend.position = "right",
                 panel.background = ggplot2::element_rect(fill = "transparent",
                                                          colour = "black",
                                                          linetype = NULL)) +
  ggplot2::scale_fill_manual(name = "Infection status", 
                             values = c("#E97777", "#82AAE3"),
                             labels = c("Infected", "Unknown"))
  # Combine and save
cowplot::plot_grid(statPlot, samplePlot, cowplot::get_legend(legendPlot), ncol = 3,
                   rel_widths = c(2, 2, 1), labels = c("a", "b", "")) %>%
  cowplot::save_plot(filename = "DiversityPlot.pdf", plot = ., base_width = 8, base_height = 3.5)



###### e. linear models ####
install.packages("Rfit")
# Linear models
# Zahl
Zahl_lm <- Rfit::rfit(formula = Zahl ~ haplotypeCount + sequenceCount + WolbachiaDetected, 
                      data = outCombined)
summary(Zahl_lm)
# Shannon
Shannon_lm <- Rfit::rfit(formula = Shannon ~ haplotypeCount + sequenceCount + WolbachiaDetected, 
                         data = outCombined)
summary(Shannon_lm)
# Chao
Chao_lm <- Rfit::rfit(formula = ChaoEstimate ~ haplotypeCount + sequenceCount + WolbachiaDetected, 
                      data = outCombined)
summary(Chao_lm)


# PLOTS
# Zahl plots
(Zahl_haplo <- ggplot2::ggplot(outCombined, aes(x = haplotypeCount, y = Zahl)) +
    ggplot2::geom_point(aes(colour = WolbachiaDetected)) +                                # scatter plot, coloured by sex
    ggplot2::labs(x = "Haplotype count", y = "Zahl") +
    ggplot2::stat_smooth(method = "lm", aes(fill = WolbachiaDetected, colour = WolbachiaDetected)) +    # adding regression lines for each sex
    ggplot2::scale_colour_manual(values = c("#FFC125", "#36648B")) +
    ggplot2::scale_fill_manual(values = c("#FFC125", "#36648B")) +
    ggplot2::theme(legend.position = c(0.2,0.8),
                   panel.background = ggplot2::element_rect(fill = "transparent",
                                                            colour = "black",
                                                            linetype = NULL)) )
(Zahl_sequence <- ggplot2::ggplot(outCombined, aes(x = sequenceCount, y = Zahl)) +
    ggplot2::geom_point(aes(colour = WolbachiaDetected)) +                                # scatter plot, coloured by sex
    ggplot2::labs(x = "Sequence count", y = "Zahl") +
    ggplot2::stat_smooth(method = "lm", aes(fill = WolbachiaDetected, colour = WolbachiaDetected)) +    # adding regression lines for each sex
    ggplot2::scale_colour_manual(values = c("#FFC125", "#36648B")) +
    ggplot2::scale_fill_manual(values = c("#FFC125", "#36648B")) +
    ggplot2::theme(legend.position = "none",
                   panel.background = ggplot2::element_rect(fill = "transparent",
                                                            colour = "black",
                                                            linetype = NULL)) )


# Shannon plots
(Shannon_haplo <- ggplot2::ggplot(outCombined, aes(x = haplotypeCount, y = Shannon)) +
    ggplot2::geom_point(aes(colour = WolbachiaDetected)) +                                # scatter plot, coloured by sex
    ggplot2::labs(x = "Haplotype count", y = "Shannon") +
    ggplot2::stat_smooth(method = "lm", aes(fill = WolbachiaDetected, colour = WolbachiaDetected)) +    # adding regression lines for each sex
    ggplot2::scale_colour_manual(values = c("#FFC125", "#36648B")) +
    ggplot2::scale_fill_manual(values = c("#FFC125", "#36648B")) +
    ggplot2::theme(legend.position = "none",
                   panel.background = ggplot2::element_rect(fill = "transparent",
                                                            colour = "black",
                                                            linetype = NULL)) )
(Shannon_sequence <- ggplot2::ggplot(outCombined, aes(x = sequenceCount, y = Shannon)) +
    ggplot2::geom_point(aes(colour = WolbachiaDetected)) +                                # scatter plot, coloured by sex
    ggplot2::labs(x = "Sequence count", y = "Shannon") +
    ggplot2::stat_smooth(method = "lm", aes(fill = WolbachiaDetected, colour = WolbachiaDetected)) +    # adding regression lines for each sex
    ggplot2::scale_colour_manual(values = c("#FFC125", "#36648B")) +
    ggplot2::scale_fill_manual(values = c("#FFC125", "#36648B")) +
    ggplot2::theme(legend.position = "none",
                   panel.background = ggplot2::element_rect(fill = "transparent",
                                                            colour = "black",
                                                            linetype = NULL)) ) 

# Chao plots
(ChaoEstimate_haplo <- ggplot2::ggplot(outCombined, aes(x = haplotypeCount, y = ChaoEstimate)) +
    ggplot2::geom_point(aes(colour = WolbachiaDetected)) +                                # scatter plot, coloured by sex
    ggplot2::labs(x = "Haplotype count", y = "ChaoEstimate") +
    ggplot2::stat_smooth(method = "lm", aes(fill = WolbachiaDetected, colour = WolbachiaDetected)) +    # adding regression lines for each sex
    ggplot2::scale_colour_manual(values = c("#FFC125", "#36648B")) +
    ggplot2::scale_fill_manual(values = c("#FFC125", "#36648B")) +
    ggplot2::theme(legend.position = "none",
                   panel.background = ggplot2::element_rect(fill = "transparent",
                                                            colour = "black",
                                                            linetype = NULL)) )
(ChaoEstimate_sequence <- ggplot2::ggplot(outCombined, aes(x = sequenceCount, y = ChaoEstimate)) +
    ggplot2::geom_point(aes(colour = WolbachiaDetected)) +                                # scatter plot, coloured by sex
    ggplot2::labs(x = "Sequence count", y = "ChaoEstimate") +
    ggplot2::stat_smooth(method = "lm", aes(fill = WolbachiaDetected, colour = WolbachiaDetected)) +    # adding regression lines for each sex
    ggplot2::scale_colour_manual(values = c("#FFC125", "#36648B")) +
    ggplot2::scale_fill_manual(values = c("#FFC125", "#36648B")) +
    ggplot2::theme(legend.position = "none",
                   panel.background = ggplot2::element_rect(fill = "transparent",
                                                            colour = "black",
                                                            linetype = NULL)) ) 

# Legend
#legendPlot <- ggplot2::ggplot(outCombined, 
#                             aes(x= sequenceCount, y= WolbachiaDetected, fill=WolbachiaDetected)) + 
# ggplot2::stat_smooth(aes(fill = WolbachiaDetected, colour = WolbachiaDetected)) +
# ggplot2::xlab("") + ggplot2::ylab("Log of count") +
# ggplot2::theme(legend.position = "right",
#                panel.background = ggplot2::element_rect(fill = "transparent",
#                                                         colour = "black",
#                                                         linetype = NULL)) +
#   ggplot2::scale_colour_manual(name = "Infection status",
#                                values = c("#FFC125", "#36648B"),
#                                labels = c("Infected", "Unknown")) +
#   ggplot2::scale_fill_manual(name = "Infection status", 
#                              values = c("#FFC125", "#36648B"),
#                              labels = c("Infected", "Unknown")))
#   

cowplot::plot_grid(Zahl_haplo, Zahl_sequence,
                   Shannon_haplo, Shannon_sequence, 
                   ChaoEstimate_haplo, ChaoEstimate_sequence,
                   #cowplot::get_legend(legendPlot),
                   ncol = 2, labels = c("a", "b", "c", "d", "e", "f", "")) %>%
  cowplot::save_plot(filename = "diversity_sampling.pdf", plot = ., 
                     base_width = 10, base_height = 10)


###### f. explore Wolbachia infections ####
# Find wolbachia-tested species
wolTestedSp <- matched %>% 
  dplyr::filter(Specimen_code %in% wolbachiaInfected$seqCode) %>%
  dplyr::distinct(Species_name)


WolTested <- matched %>%
  # !!! OPTIONAL filter to ONLY Wolbachia individuals
  #dplyr::filter(Specimen_code %in% wolbachiaInfected$seqCode) %>%
  dplyr::left_join(wolbachiaInfected %>%
                     dplyr::distinct(seqCode, .keep_all = TRUE),
                   by = c("Specimen_code" = "seqCode")) %>%
  dplyr::group_by(Species_name) %>%
  dplyr::mutate(percentInfected = round((sum(WolbachiaPositive, na.rm = TRUE)/dplyr::n())*100, 0) ,
                sampleSize = dplyr::n()) %>%
  dplyr::distinct(Species_name, .keep_all = TRUE) %>% 
  dplyr::select(Species_name, percentInfected, sampleSize) %>%
  # Merge with statistics
  dplyr::left_join(outCombined_complete, by = "Species_name")

# Get species with n >= 3 and ONLY for wolbachia-tested species
WolTested_3 <- WolTested %>%
  dplyr::filter(sampleSize > 2) %>%
  dplyr::filter(Species_name %in% wolTestedSp$Species_name)

# infection and sample size + richness
infection_lm <- Rfit::rfit(formula = percentInfected ~ sampleSize + ChaoEstimate, 
                      data = WolTested)
summary(infection_lm)

# infection and diversity + richness
infectionDiv_lm <- Rfit::rfit(formula = percentInfected ~ Zahl + ChaoEstimate, 
                           data = WolTested)
summary(infectionDiv_lm)

shapiro.test( WolTested_3$percentInfected)
shapiro.test( WolTested_3$sampleSize)
shapiro.test( WolTested$ChaoEstimate)
shapiro.test( WolTested_3$Zahl )


cor.test(x= WolTested_3$percentInfected, y = WolTested_3$sampleSize, method = "spearm",
         data = WolTested_3, alternative = "two.sided")
cor.test(x= WolTested_3$percentInfected, y = WolTested_3$Zahl, method = "spearm",
         data = WolTested_3, alternative = "two.sided")
cor.test(x= WolTested_3$percentInfected, y = WolTested_3$Shannon, method = "spearm",
         data = WolTested_3, alternative = "two.sided")
cor.test(x= WolTested_3$percentInfected, y = WolTested_3$ChaoEstimate, method = "spearm",
         data = WolTested_3, alternative = "two.sided")

par(mfrow = c(2, 2))
plot(WolTested_3$percentInfected, WolTested_3$sampleSize)
plot(WolTested_3$percentInfected, WolTested_3$Zahl)
plot(WolTested_3$percentInfected, WolTested_3$Shannon)
plot(WolTested_3$percentInfected, WolTested_3$ChaoEstimate)

  # Diversity vs richness
DivRich <- WolTested %>%
  dplyr::filter(!is.na(WolbachiaDetected)) 

cor.test(x= DivRich$Zahl, y = DivRich$ChaoEstimate, method = "spearm",
         data = DivRich, alternative = "two.sided")
# Diversity and richness
DivRich_lm <- Rfit::rfit(formula = Zahl ~ ChaoEstimate + WolbachiaDetected, 
                         data = DivRich)
summary(DivRich_lm)

  # Plot pattern
(DivRich_fijiensis <- ggplot2::ggplot(DivRich, aes(x = ChaoEstimate, y = Zahl)) +
    ggplot2::geom_point(aes(colour = WolbachiaDetected)) +                                # scatter plot, coloured by sex
    ggplot2::labs(x = "Chao Estimate", y = "Zahl") +
    ggplot2::stat_smooth(method = "lm", formula = y ~ log(x),
                         aes(fill = WolbachiaDetected, colour = WolbachiaDetected)) +    # adding regression lines for each sex
    ggplot2::scale_colour_manual(values = c("#FFC125", "#36648B")) +
    ggplot2::scale_fill_manual(values = c("#FFC125", "#36648B")) +
    ggplot2::theme(legend.position = c(0.2,0.9),
                   panel.background = ggplot2::element_rect(fill = "transparent",
                                                            colour = "black",
                                                            linetype = NULL)) ) 
(DivRich_noFijiensis <- ggplot2::ggplot(DivRich %>%
                                          dplyr::filter(!stringr::str_detect(Species_name, "fijiensis")),
                                        aes(x = ChaoEstimate, y = Zahl)) +
    ggplot2::geom_point(aes(colour = WolbachiaDetected)) +                                # scatter plot, coloured by sex
    ggplot2::labs(x = "Chao Estimate", y = "Zahl") +
    ggplot2::stat_smooth(method = "lm", formula = y ~ log(x), 
                         aes(fill = WolbachiaDetected, colour = WolbachiaDetected)) +    # adding regression lines for each sex
    ggplot2::scale_colour_manual(values = c("#FFC125", "#36648B")) +
    ggplot2::scale_fill_manual(values = c("#FFC125", "#36648B")) +
    ggplot2::theme(legend.position = "none",
                   panel.background = ggplot2::element_rect(fill = "transparent",
                                                            colour = "black",
                                                            linetype = NULL)) ) 
cowplot::plot_grid(DivRich_fijiensis, DivRich_noFijiensis,
                   #cowplot::get_legend(legendPlot_DivRich),
                   ncol = 2, labels = c("a", "b")) %>%
  cowplot::save_plot(filename = "DivRich.pdf", plot = ., 
                     base_width = 10, base_height = 5)






#### 2. Maps ####
  ##### 2.1 prepare data ####
# Read in the occurrence data
OccData <- readr::read_csv(paste0(RootPath, "/HomalictusCollectionData_2018.csv"),
                           col_types = readr::cols(.default = "c")) %>%
  dplyr::rename(
    decimalLongitude = Longitude,
    decimalLatitude = Latitude)
wolbachiaSpecies = readr::read_csv("wolbachiaSpecies.csv")


# Create the points file
fijiPoints <- OccData %>% 
  dplyr::mutate(wolStatus = dplyr::if_else(Specimen_code %in% wolbachiaSpecies$homaSp,
                                           "Infected", "Unknown")) %>%
  tidyr::drop_na(c(decimalLongitude, decimalLatitude)) %>%
  dplyr::mutate(decimalLongitude = as.numeric(decimalLongitude),
                decimalLatitude = as.numeric(decimalLatitude)) %>%
  sf::st_as_sf(coords = c("decimalLongitude", "decimalLatitude")) %>%
  # Set the right CRS
  sf::st_set_crs(terra::crs("EPSG:4326")) %>%
  #terra::vect(geom = c("decimalLongitude", "decimalLatitude")) %>%
  sf::st_transform(., crs = terra::crs("EPSG:3460"))

  ##### 2.2 Raster maps ####
source(paste0(RootPath, "/FjRasterPointMapR.R"))
FjRasterPointMapR(
  mapData = fijiPoints,
  #spColours = NULL,
  filename = "rast_FijiHoma_noSpecies.pdf",
  outpath = RootPath,
  width = 6, height = 6, units = "in",
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
  insetX = 0.1, 
  insetY = 0.1,
  insetWidth = 0.35, 
  insetHeight = 0.35,
  
  # OPTIONAL:
  ptSize = 1.5,
  pchAll = 19,
     # background point and wolbachia point oclours
    colPts = "black",
    colWol = "white",
  # Raster colours
  naMapCol = "aliceblue",
  rasterGradient = "ggthemes::Blue",
  # 1 or -1 to change directions
  colourDirection = -1,
  # point alpha (opacity)
  mapAlpha = 0.8,
  mapTitle = "Fijian Homalictus occurrences")

##### 2.3 Convex hulls ####
  ###### a. get areas ####
  # Make a vector of each species name
fijiSpecies <- unique(fijiPoints$Species_name)
  # Create an empty tibble
loopTibble <- tibble::tibble(Species_name = NA_character_,
                             polygonArea_m2 = NA_integer_)

for(i in 1:length(unique(fijiPoints$Species_name))){
    # select the points for the ith species
  sp_i <- fijiPoints %>%
    dplyr::filter(Species_name == fijiSpecies[[i]])
    # Get the polygon area for that species
  polygonArea_i <- sf::st_convex_hull(sp_i %>% sf::st_union()) %>% sf::st_area() %>%
    as.double()
    # add to loop tibble
  loopTibble <- loopTibble %>%
    dplyr::bind_rows( tibble::tibble(Species_name = fijiSpecies[[i]],
                                     polygonArea_m2 = polygonArea_i/1000000))
}# Finish loop


  ###### b. combine data ####
areaPlus <- WolTested %>%
  dplyr::left_join(loopTibble, by = "Species_name") %>%
    # select only wolbachia-tested species
  dplyr::filter(Species_name %in% wolTestedSp$Species_name)

  ###### c. run tests ####
  # Test area normality _ NOT normal
shapiro.test( areaPlus$polygonArea_m2 )
hist(areaPlus$polygonArea_m2)

  # Simple correlation between area and percentInfected
  # We expect a positive relatinoship
cor.test(x = areaPlus$percentInfected, y = areaPlus$polygonArea_m2, method = "spearm",
         data = areaPlus, alternative = "greater")
  # test area and diversity
cor.test(x = areaPlus$Zahl, y = areaPlus$polygonArea_m2, method = "spearm",
         data = areaPlus, alternative = "greater")
# test area and richness
cor.test(x = areaPlus$ChaoEstimate, y = areaPlus$polygonArea_m2, method = "spearm",
         data = areaPlus, alternative = "greater")




# infection and sample richness, diversity, area
areaInfRate_lm <- Rfit::rfit(formula = percentInfected ~ ChaoEstimate + Zahl + polygonArea_m2, 
                           data = areaPlus)
summary(areaInfRate_lm)
# area and sample richness, diversity
areaDivRich_lm <- Rfit::rfit(formula = polygonArea_m2 ~ ChaoEstimate + Zahl + WolbachiaDetected, 
                             data = areaPlus)
summary(areaDivRich_lm)


##### 2.X OLD vector maps ####
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


  
  


