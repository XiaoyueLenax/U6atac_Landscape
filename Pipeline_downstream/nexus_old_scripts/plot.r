if (!require("ggplot2")) install.packages("ggplot2")
library(ggplot2)
library(readxl)
df <- read_excel("combined_data_7.xlsx", sheet = 1)
df <- df %>%
  mutate(
    minor = x_fold_minor <= 0.7 & x_fold_major >= 1,
    major = x_fold_minor >= 1 & x_fold_major <= 0.5,
    condition = case_when(
      minor & !major ~ "minor",
      major ~ "major",
      TRUE ~ "other"
    )
  )

# Create the plot with the updated conditions
volcano_plot <- ggplot(df, aes(x = x_fold_minor, y = x_fold_major)) +
  geom_point(aes(color = condition), alpha = 0.5) +
  scale_color_manual(values = c("minor" = "red", "major" = "blue", "other" = "black")) +
  theme_minimal() +
  labs(title = "Volcano Plot: x_fold_minor vs. x_fold_major",
       x = "x_fold_minor",
       y = "x_fold_major",
       color = "Condition") +
  geom_text(aes(label = df[[1]]), vjust = 1, hjust = 1, size = 3, check_overlap = TRUE) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "green") +
  geom_vline(xintercept = 0.7, linetype = "dashed", color = "green") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "purple") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "purple") +
  theme(legend.position = "right")

# Display the plot
print(volcano_plot)

# Save the plot
#ggsave("volcano_plot.png", plot = volcano_plot, width = 10, height = 8, dpi = 300)
pdf("sample_7.pdf")
print(volcano_plot)     # Plot 1 --> in the first page of PDF    # Plot 2 ---> in the second page of the PDF
dev.off()
  