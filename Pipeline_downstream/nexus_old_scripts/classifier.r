# For now, load the manually prepared excel
install.packages("readxl")
# Load the readxl package
library(readxl)
library(readxl)
library(writexl)
# Prepare to filter rows with conditions
library(dplyr)

# load input
library_SelleckChem_df <- read.csv("20240105_L3900_SelleckChem_Clinical.tsv",header = TRUE, sep = "\t")

# Reading an Excel file
input_data <- read_excel("20240304-131356_240229MR_2_example.xlsx", sheet = 1)
head(input_data) # You see NAs here bc of the manual input - can be changed later

# get the x-fole of ctr data, here for minor is column 11
head(input_data[,11])


# for major, is column 20 - making sure they look the same.
head(input_data[,20])



#Implementing rule number 1: less than 0.7 in minor x-folds, across all 3
filtered_data_minor_2 <- input_data[input_data[, 11] <= 0.7 & input_data[, 20] >= 1, ]
filtered_data_major_2 <- input_data[input_data[, 11] >= 1 & input_data[, 20] <= 0.5, ]
# 

write_xlsx(input_data, "x_fold_2_corrected.xlsx")
# This is a list of candicates in the 2 folder. now 
#write.table(filtered_data_major_2, file = "filtered_data_major_2.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
#write.table(filtered_data_minor_2, file = "filtered_data_minor_2.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write_xlsx(filtered_data_major_2, "filtered_data_major_2.xlsx")
write_xlsx(filtered_data_minor_2, "filtered_data_minor_2.xlsx")
#----------------------------------------- for the 3rd folder -----------------------------------

## Merge to prepare the files:
# Specify the paths to the Excel files
files <- c("Tecan_m1000/1_Fluo_mCherry/20240304-131356_240229MR_2.xlsx",
           "Tecan_m1000/3_Lum_Minor/20240304-140944_240229MR_2.xlsx",
           "Tecan_m1000/2_Lum_Major/20240304-132818_240229MR_2.xlsx")

# Initialize an empty list to store data frames
data_list <- list()

# Loop through the files and read each one, adding the data frame to the list
for (file in files) {
  data_list[[file]] <- read_excel(file)
}

# Find the maximum number of rows among all data frames
max_rows <- max(sapply(data_list, nrow))

# Adjust each data frame to have the same number of rows by adding NA rows if necessary
adjusted_data_list <- lapply(data_list, function(df) {
  nr <- nrow(df)
  if(nr < max_rows) {
    # Create a data frame of NAs with the correct number of missing rows and the same column names as df
    na_df <- as.data.frame(matrix(ncol = ncol(df), nrow = max_rows - nr))
    names(na_df) <- names(df)
    # Combine df with na_df
    df <- rbind(df, na_df)
  }
  return(df)
})
# Combine all adjusted data frames side by side
input_data_3 <- do.call(cbind, adjusted_data_list)
# Write the combined data to a new Excel file
write_xlsx(input_data_3, "combined_data_3.xlsx")

input_data_3 <- read_excel("combined_data_3.xlsx", sheet = 1)
#Slice to make sure the correct colums, since this can differ from format
head(input_data_3[,3])
input_data_3[,11]

# Remove negative values, replace with NAs:
# First, divide the value with mcherry
#input_data_3$calc_column_minor <- ifelse(input_data_3[,7] >= 0 & input_data_3[,3] >= 0, input_data_3[,7] / input_data_3[,3], NA)
#input_data_3$calc_column_major <- ifelse(input_data_3[,11] >= 0 & input_data_3[,3] >= 0, input_data_3[,11] / input_data_3[,3], NA)
input_data_3$calc_column_minor <- input_data_3[,7]/input_data_3[,3]
input_data_3$calc_column_major <- input_data_3[,11]/input_data_3[,3]
head(input_data_3$calc_column_minor)
head(input_data_3$calc_column_major)
input_data_3$calc_column_minor[input_data_3$calc_column_minor < 0] <- NA
input_data_3$calc_column_major[input_data_3$calc_column_major < 0] <- NA


# Now look at the controls: A2, D23, D24, H23, H24, I1, I2, L23,L24,M1,M2,P23,P24
input_data_3[287,]
A2 <- input_data_3$calc_column_minor[1,]
a2 <- input_data_3$calc_column_major[1,]

D23 <- input_data_3$calc_column_minor[94,]
D24 <- input_data_3$calc_column_minor[95,]
D <- mean(D23,D24)
d23 <- input_data_3$calc_column_major[94,]
d24 <- input_data_3$calc_column_major[95,]
d <- mean(d23,d24)

H23 <- input_data_3$calc_column_minor[190,]
H24 <- input_data_3$calc_column_minor[191,]
H <- mean(H23,H24)
h23 <- input_data_3$calc_column_major[190,]
h24 <- input_data_3$calc_column_major[191,]
h <- mean(h23,h24)

I1 <- input_data_3$calc_column_minor[192,]
I2 <- input_data_3$calc_column_minor[193,]
I <- mean(I1, I2)
i1 <- input_data_3$calc_column_major[192,]
i2 <- input_data_3$calc_column_major[193,]
i <- mean(i1, i2)

L23 <- input_data_3$calc_column_minor[286,]
L24 <- input_data_3$calc_column_minor[287,]
L <- mean(L23, L24)
l23 <- input_data_3$calc_column_major[286,]
l24 <- input_data_3$calc_column_major[287,]
l <- mean(l23, l24)

M1 <- input_data_3$calc_column_minor[288,]
M2 <- input_data_3$calc_column_minor[289,]
M <- mean(M1, M2)
m1 <- input_data_3$calc_column_major[288,]
m2 <- input_data_3$calc_column_major[289,]
m <- mean(m1, m2)

P24 <- input_data_3$calc_column_minor[383,]
P23 <- input_data_3$calc_column_minor[382,]
P <- mean(P23, P24)
p24 <- input_data_3$calc_column_minor[383,]
p23 <- input_data_3$calc_column_minor[382,]
p <- mean(p23, p24)

control_minor <- mean(A2,D,H,I,L,M,P)
control_major <- mean(a2,d,h,l,m,p)
control_major;control_minor

# Calculation for minor and major
input_data_3$x_fold_minor <- 1 / control_minor * input_data_3$calc_column_minor 
input_data_3$x_fold_major <- 1 / control_minor * input_data_3$calc_column_major



filtered_data_minor_3 <- input_data_3[input_data_3$x_fold_minor <= 0.7 & input_data_3$x_fold_major >= 1, ]
filtered_data_major_3 <- input_data_3[input_data_3$x_fold_minor >= 1 & input_data_3$x_fold_major <= 0.5, ]

# Display the first few rows of filtered_data_minor_3 (corrected object name)
head(filtered_data_minor_3)
head(filtered_data_major_3)

write_xlsx(input_data_3, "x_fold_3_corrected.xlsx")


write_xlsx(filtered_data_major_3, "filtered_data_major_3.xlsx")
write_xlsx(filtered_data_minor_3, "filtered_data_minor_3.xlsx")

# -----------------------------------------------------------------------------------------------
