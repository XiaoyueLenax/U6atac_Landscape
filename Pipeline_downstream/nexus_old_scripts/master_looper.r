# ============================ Library prep ====================================
if (!require("openxlsx")) install.packages("openxlsx")
if (!require("readxl")) install.packages("readxl")
if (!require("dplyr")) install.packages("dplyr")
library(openxlsx)
library(readxl)
library(dplyr)
# ==============================================================================
# =============================== WIP stuff ====================================
# List of worries:
#    1. Why R has no wildcard like python???
#    2. files named differently
#    3. control names hard to get?


# Input
# Define the directory where your files are located, assume to be in the folder.
# Maybe let them name the files manually.
files <- c("Tecan_m1000/1_Fluo_mCherry/20240304-131356_240229MR_2.xlsx",
           "Tecan_m1000/3_Lum_Minor/20240304-140944_240229MR_2.xlsx",
           "Tecan_m1000/2_Lum_Major/20240304-132818_240229MR_2.xlsx")

# Control F to replace:
# place_holder.xlsx

#------------------------------ End of user input ------------------------------
# Function to merge - create a file that combines the mcherry, minor and major.
merge_files <- function(files,out_name){
  data_list <- list()

  # Loop through the files and read each one, adding the data frame to the list
  for (file in files) {
    data_list[[file]] <- read_excel(file)
  }
  # Find the maximum number of rows among all data frames
  max_rows <- max(sapply(data_list, nrow))
  # Adjust each data frame to have the same number of rows by adding 
  # NA rows if necessary
  adjusted_data_list <- lapply(data_list, function(df) {
    nr <- nrow(df)
    if(nr < max_rows) {
      # Create a data frame of NAs with the correct number of missing rows and the 
      # same column names as df
      na_df <- as.data.frame(matrix(ncol = ncol(df), nrow = max_rows - nr))
      names(na_df) <- names(df)
      # Combine df with na_df
      df <- rbind(df, na_df)
    }
    return(df)
  })
  # Combine all adjusted data frames side by side
  input_data <- do.call(cbind, adjusted_data_list)
  # Write the combined data to a new Excel file
  #data <- as.data.frame(data)
  write_xlsx(input_data, out_name)
}

# ------------------------------------------------------------------------------
# Define a function that performs the mCherry division, for each combined file
perform_division <- function(combined_file_path) {
  data <- read_excel(combined_file_path)

  # Remove negative values by replacing with NA
  data[data[[11]] < 0, 11] <- NA  
  data[data[[3]] < 0, 3] <- NA    

  # Calculate ratios
  ratio_minor <- data[[7]] / data[[3]]  
  ratio_major <- data[[11]] / data[[3]]  

  data$ratio_minor <- data$ratio_minor[[1]]
  data$ratio_major <- data$ratio_major[[1]]
  # After correcting the structure, proceed to write the data frame to an Excel file
  write_xlsx(data,combined_file_path)

  # Return both ratios as a list
  return(list(ratio_minor = ratio_minor, ratio_major = ratio_major))
}
#-------------------------------------------------------------------------------
# Define a function that calculates the control averages
ctrl_average <- function(file_path){
  data <- read_excel(file_path)
  data[287,] # This should return L24
  # in this case we didnt have A1?
  #A1 <- data$calc_column_minor[1,]
  #a1 <- data$calc_column_minor[1,]
  A2 <- data$calc_column_minor[1,]
  a2 <- data$calc_column_major[1,]
  # A <- mean(A1, A2)
  # a <- mean(a1, a2)

  D23 <- data$calc_column_minor[94,]
  D24 <- data$calc_column_minor[95,]
  D <- mean(D23,D24)
  d23 <- data$calc_column_major[94,]
  d24 <- data$calc_column_major[95,]
  d <- mean(d23,d24)

  H23 <- data$calc_column_minor[190,]
  H24 <- data$calc_column_minor[191,]
  H <- mean(H23,H24)
  h23 <- data$calc_column_major[190,]
  h24 <- data$calc_column_major[191,]
  h <- mean(h23,h24)

  I1 <- data$calc_column_minor[192,]
  I2 <- data$calc_column_minor[193,]
  I <- mean(I1, I2)
  i1 <- data$calc_column_major[192,]
  i2 <- data$calc_column_major[193,]
  i <- mean(i1, i2)

  L23 <- data$calc_column_minor[286,]
  L24 <- data$calc_column_minor[287,]
  L <- mean(L23, L24)
  l23 <- data$calc_column_major[286,]
  l24 <- data$calc_column_major[287,]
  l <- mean(l23, l24)

  M1 <- data$calc_column_minor[288,]
  M2 <- data$calc_column_minor[289,]
  M <- mean(M1, M2)
  m1 <- data$calc_column_major[288,]
  m2 <- data$calc_column_major[289,]
  m <- mean(m1, m2)

  P24 <- data$calc_column_minor[383,]
  P23 <- data$calc_column_minor[382,]
  P <- mean(P23, P24)
  p24 <- data$calc_column_minor[383,]
  p23 <- data$calc_column_minor[382,]
  p <- mean(p23, p24)

  control_minor <- mean(A2,D,H,I,L,M,P)
  control_major <- mean(a2,d,h,l,m,p)
  return(control_minor, control_major)
}
# ------------------------------------------------------------------------------
# Define a function that calculates the x-fold ratio
perform_x_fold <- function(combined_file_path,control_minor,control_major){
  data <- read_excel(file_path)

  data$x_fold_minor <- 1 / control_minor * data$calc_column_minor 
  data$x_fold_major <- 1/ control_major * data$calc_column_major
  head(data$x_fold_major)
  head(data$x_fold_minor)

  data$x_fold_minor <- data$x_fold_minor[[1]]
  data$x_fold_major <- data$x_fold_major[[1]]
  # After correcting the structure,
  # proceed to write the data frame to an Excel file
  write_xlsx(data, combined_file_path)
}
# ------------------------------------------------------------------------------
# Define a function that does the strict control - WIP
strict_filter <- function(combined_file_path){

}