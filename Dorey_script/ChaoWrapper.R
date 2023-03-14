# This function was written by JAmes B. Dorey on the 14th of March 2023 to wrap ChaoSpecies
# and output table for multiple species at once

ChaoWrapper <- function(
    FJHoma_haplotypes = NULL){

# Duplicate a function within ChaoSpecies to not run data-poor species
f <- function(i, data) {
  length(data[which(data == i)])
}
# make an empty list for the loop
loopList <- tibble::lst()
diversityTable <- tibble::tibble(species = character(),
                                 Name = character(),
                                 Estimate = numeric(),
                                 's.e.' = numeric(),
                                 '95%Lower' = numeric(),
                                 '95%Upper' = numeric())
basicTable <- tibble::tibble()
for(i in 1:ncol(FJHoma_haplotypes)){
  # choose the data for the first species
  inputData_i <- FJHoma_haplotypes[,i]
  # Don't run species where there are no counts that are less than k
  # i.e. ONLY run species that won't throw an error due to weird/poor cases of data
  if(any(inputData_i[[1]] > 0 & inputData_i[[1]] < 5) & 
     !(f(1, inputData_i[[1]]) == sum(inputData_i[[1]]))
  ){
    spOutput <- SpadeR::ChaoSpecies(
      data = inputData_i[[1]],
      datatype = c("abundance"),
      k = 5, conf = 0.95)
    ## 1
    # Get the diversity measures and re-format
    diversityOutput <- spOutput$Species_table %>% as.data.frame() %>% 
      tibble::rownames_to_column() %>% 
      dplyr::mutate(rowname = stringr::str_squish(rowname))
    # Add these data to one table
    diversityTable <- diversityTable %>%
      dplyr::bind_rows(
        tibble::tibble(
          species = names(inputData_i),
          Name = c(names(inputData_i), diversityOutput$rowname),
          Estimate = c(NA_integer_, diversityOutput$Estimate),
          's.e.' = c(NA_integer_, diversityOutput$s.e.),
          '95%Lower' = c(NA_integer_, diversityOutput$`95%Lower`),
          '95%Upper' = c(NA_integer_, diversityOutput$`95%Upper`))
      )
    
    ## 2
    # Get Basic information into one table
    basicOutput <- spOutput$Basic_data_information %>% as.data.frame() %>% 
      tibble::rownames_to_column() %>% 
      dplyr::mutate(rowname = stringr::str_squish(rowname))
    
    # Set up the data with i == 1
    if(i == 1){
      basicTable <- basicOutput
      # rename Value to the species name
      names(basicTable) <- c("rowname", "Variable", names(inputData_i))
    }else{# For i > 1...
      # rename Value to the species name
      names(basicOutput) <- c("rowname", "Variable", names(inputData_i))
      # join the new species' data
      basicTable <- basicTable %>%
        dplyr::left_join(basicOutput %>% dplyr::select(!Variable),
                         by = "rowname"
        )}
  }# END if
}# END loop

  return(list( basicTable, diversityTable) %>%
           setNames(c( "basicTable", "diversityTable")))
} # END function