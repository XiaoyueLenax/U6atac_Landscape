library(readxl)

# Function to read the first column of Excel files within specified ranges
read_first_col <- function(file_path, range = NULL) {
  if (is.null(range)) {
    return(read_excel(file_path)[, 1])
  } else {
    return(read_excel(file_path, range = range)[[1]])
  }
}

# Function to find common strings
find_common_strings <- function(data_list) {
  common_3_times <- Reduce(intersect, data_list)
  pairwise_intersects <- unique(unlist(lapply(1:length(data_list), function(i) {
    unlist(lapply((i+1):length(data_list), function(j) {
      intersect(data_list[[i]], data_list[[j]])
    }))
  })))
  common_2_times <- pairwise_intersects[!pairwise_intersects %in% common_3_times]
  list(common_3_times = common_3_times, common_2_times = common_2_times)
}

# Read data and find common strings for both major and minor
process_data <- function(prefix, indices) {
  data_list <- lapply(indices, function(i) {
    file_path <- paste0("filtered/filtered_data_", prefix, "_", i, ".xlsx")
    if (i == 2 && prefix == "major") {
      read_first_col(file_path)  # major_2 does not specify a range
    } else {
      read_first_col(file_path, range = "A1:A100")
    }
  })
  
  find_common_strings(data_list)
}

# Major data
major_indices <- 2:7
major_results <- process_data("major", major_indices)

# Minor data
minor_indices <- c(2, 4, 6, 3, 5, 7)  # Assuming you want to process these indices for minor
minor_results <- process_data("minor", minor_indices)

# Printing results
print_results <- function(results, prefix, indices) {
  cat("Strings matched 3 times among", prefix, "files", paste(indices, collapse=", "), ":\n")
  print(results$common_3_times)
  
  cat("Strings matched 2 times among", prefix, "files", paste(indices, collapse=", "), ":\n")
  print(results$common_2_times)
}

# Print results for major and minor data
print_results(major_results, "major", major_indices)
# Assuming you want to split the minor results printing
print_results(minor_results, "minor", minor_indices[1:3])
print_results(minor_results, "minor", minor_indices[4:6])

