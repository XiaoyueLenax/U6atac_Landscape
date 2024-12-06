# ------ LIBRARY -------
library(readxl)
# -------- QC ----------
data <- read_excel("prostate_rawdata_100424_annon.xlsx")
head(data)
str(data)
colnames(data) #22 colnames
# ----------------------