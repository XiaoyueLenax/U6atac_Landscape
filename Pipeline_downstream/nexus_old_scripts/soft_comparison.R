library(readxl)
library(dplyr)

# Majors
major_2 <- read_excel("filtered/filtered_data_major_2.xlsx")
major_4 <- read_excel("filtered/filtered_data_major_4.xlsx")
major_6 <- read_excel("filtered/filtered_data_major_6.xlsx")
major_3 <- read_excel("filtered/filtered_data_major_3.xlsx")
major_5 <- read_excel("filtered/filtered_data_major_5.xlsx")
major_7 <- read_excel("filtered/filtered_data_major_7.xlsx")

# for minors
minor_2 <- read_excel("filtered_data_minor_2_strict.xlsx")
minor_4 <- read_excel("filtered/filtered_data_minor_4.xlsx")
minor_6 <- read_excel("filtered/filtered_data_minor_6.xlsx")
minor_3 <- read_excel("filtered/filtered_data_minor_3.xlsx")
minor_5 <- read_excel("filtered/filtered_data_minor_5.xlsx")
minor_7 <- read_excel("filtered/filtered_data_minor_7.xlsx")

# Comparison Files: read the combined data
combined_2 <- read_excel("combined_data_2.xlsx")
combined_3 <- read_excel("combined_data_3.xlsx")
combined_4 <- read_excel("combined_data_4.xlsx")
combined_5 <- read_excel("combined_data_5.xlsx")
combined_6 <- read_excel("combined_data_6.xlsx")
combined_7 <- read_excel("combined_data_7.xlsx")
#
combined_2[,14] #major ratio
combined_2[,13] #minor ratio

combined_4[,15] #minor x_fold
combined_4[,16] #major x_fold

# start with 2 -----------------------------------------------------------------
candidates <- minor_2[[1]]
# Filter rows based on the candidate list
# soft to check
# Assuming the first column of combined_3 contains the identifiers
matching_candidates <- combined_4 %>%
  filter(combined_4[[1]] %in% candidates & combined_4[[15]] < 1) %>%
  pull(1)  # Extracts the identifiers from the first column of the filtered data
  matching_candidates # Only finds A24.
#--- 
matching_candidates <- combined_4 %>%
  filter(combined_4[[1]] %in% candidates & 
           combined_4[[15]] < 1 & 
           combined_4[[14]] > 1) %>%
  pull(1)  # Extracts the identifiers from the first column of the filtered data
  print(matching_candidates) # 0, so I guess we dont need to proceed to next file 6...

# start with 3 -----------------------------------------------------------------
candidates <- minor_3[[1]]
matching_candidates <- combined_5 %>%
  filter(combined_5[[1]] %in% candidates & 
           combined_5[[15]] < 1 & 
           combined_5[[14]] > 1) %>%
  pull(1)  # Extracts the identifiers from the first column of the filtered data
matching_candidates # Also 0.

# Try with major? --------------------------------------------------------------
candidates <- major_2[[1]]
matching_candidates <- combined_4 %>%
  filter(combined_4[[1]] %in% candidates & 
           combined_4[[14]] < 1 & 
           combined_4[[15]] > 1) %>%
  pull(1)  # Extracts the identifiers from the first column of the filtered data
print(matching_candidates) # so we could find 2 controls....

# Try again with minor , with different reference ------------------------------
candidates <- minor_4[[1]]
matching_candidates <- combined_2 %>%
  filter(combined_2[[1]] %in% candidates & 
           combined_2[[15]] < 1 & 
           combined_2[[14]] > 1) %>%
  pull(1)  # Extracts the identifiers from the first column of the filtered data
matching_candidates # Still 0. 
