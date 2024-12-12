# Load required libraries
library(readr)  # For reading data
library(ggplot2)  # For data visualization
library(dplyr)
# 46, 47, and 48. These correspond with the 96-hour 
#siU6atac, scrambled, and ntc data respectively.

# Making a combined file
#data <- read_delim("/Users/xd23w713/Documents/xd_ws/IsoDE_raw_LnCap/46-47-DE_geneFPKM", delim = "\t") %>%
#  mutate(Source="46-47")
#data2 <- read_delim("/Users/xd23w713/Documents/xd_ws/IsoDE_raw_LnCap/46-48-DE_geneFPKM", delim = "\t") %>%
#  mutate(Source="46-48")
#data3 <- read_delim("/Users/xd23w713/Documents/xd_ws/IsoDE_raw_LnCap/47-48-DE_geneFPKM", delim = "\t") %>%
#  mutate(Source="47-48")

colnames(data)
data <- read_delim("/Users/xd23w713/Documents/xd_ws/IsoDE_raw_LnCap/46-47-DE_geneFPKM", delim = "\t") %>%
  rename("Confident log2 FC siU6atac - scr" = "Confident log2 FC") %>%
  rename("Single run log2 FC siU6atac - scr" = " Single run log2 FC") %>%
  rename("c1 average FPKM siU6atact - scr" = " c1 average FPKM") %>%
  rename("c2 average FPKM siU6atac - scr" = "c2 average FPKM")


data2 <- read_delim("/Users/xd23w713/Documents/xd_ws/IsoDE_raw_LnCap/46-48-DE_geneFPKM", delim = "\t") %>%
  rename("Confident log2 FC siU6atac - ntc" = "Confident log2 FC") %>%
  rename("Single run log2 FC siU6atac - ntc" = " Single run log2 FC") %>%
  rename("c1 average FPKM siU6atact - ntc" = " c1 average FPKM") %>%
  rename("c2 average FPKM siU6atac - ntc" = "c2 average FPKM") %>%
  select(-matches("Gene ID", ignore.case = TRUE))

data3 <- read_delim("/Users/xd23w713/Documents/xd_ws/IsoDE_raw_LnCap/47-48-DE_geneFPKM", delim = "\t") %>%
  rename("Confident log2 FC scr - ntc" = "Confident log2 FC") %>%
  rename("Single run log2 FC scr - ntc" = " Single run log2 FC") %>%
  rename("c1 average FPKM scr - ntc" = " c1 average FPKM") %>%
  rename("c2 average FPKM scre - ntc" = "c2 average FPKM") %>%
  select(-matches("Gene ID", ignore.case = TRUE))


combined_data <- cbind(data,data2,data3)
colnames(combined_data)

# Export the combined data as a new file
write_delim(combined_data, "/Users/xd23w713/Documents/xd_ws/IsoDE_raw_LnCap/combined_IsoDE_LnCap_data.csv", delim = ",")