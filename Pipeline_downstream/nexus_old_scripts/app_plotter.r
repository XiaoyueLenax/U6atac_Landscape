library(shiny)
library(ggplot2)
library(readxl)
library(dplyr)
#install.packages("readxl")

ui <- fluidPage(
  titlePanel("Dynamic Threshold Checker"),
  sidebarLayout(
    sidebarPanel(
      # Dropdown menu for selecting the sample
      selectInput("sample_selection", "Choose Sample",
                  choices = c("Sample 2" = "2", "Sample 3" = "3", "Sample 4" = "4", 
                              "Sample 5" = "5", "Sample 6" = "6", "Sample 7" = "7")),
      
      sliderInput("minor_threshold_x", "Minor: Minor Upper Threshold",
                  min = 0, max = 1.5, value = 0.7, step = 0.1),
      sliderInput("minor_threshold_y", "Minor: Major Lower Threshold",
                  min = 0, max = 1.5, value = 1, step = 0.1),
      sliderInput("major_threshold_y", "Major: Major Upper Threshold",
                  min = 0, max = 1.5, value = 0.5, step = 0.1),
      sliderInput("major_threshold_x", "Major: Minor Lower Threshold",
                  min = 0, max = 1.5, value = 1, step = 0.1),
      
      # Checkbox for toggling gene names
      checkboxInput("show_gene_names", "Show Gene Names", FALSE)
    ),
    mainPanel(
      plotOutput("volcanoPlot"),
      verbatimTextOutput("listMinor"),
      verbatimTextOutput("listMajor")
    )
  )
)

server <- function(input, output) {
  # Reactive expression to read and process the data
  df_reactive <- reactive({
    file_name <- paste0("nexus/nexus-drug-screen/combined_data_", input$sample_selection, ".xlsx")
    df <- read_excel(file_name, sheet = 1) %>%
      mutate(
        minor = x_fold_minor <= input$minor_threshold_x & x_fold_major >= input$minor_threshold_y,
        major = x_fold_minor >= input$major_threshold_x & x_fold_major <= input$major_threshold_y,
        condition = case_when(
          minor & !major ~ "minor",
          major ~ "major",
          TRUE ~ "other"
        )
      )
    return(df)
  })
  output$volcanoPlot <- renderPlot({
    df <- df_reactive()
    p <- ggplot(df, aes(x = x_fold_minor, y = x_fold_major)) +
      geom_point(aes(color = condition), alpha = 0.5, size = 6) +
      scale_color_manual(values = c("minor" = "red", "major" = "blue", "other" = NA)) +
      theme_minimal() +
      labs(title = paste("Scatter Plot: Sample", input$sample_selection),
           x = "x_fold_minor", y = "x_fold_major", color = "Candidates") +
      geom_hline(yintercept = input$minor_threshold_y, linetype = "dashed", color = "pink") +
      geom_vline(xintercept = input$minor_threshold_x, linetype = "dashed", color = "pink") +
      geom_hline(yintercept = input$major_threshold_y, linetype = "dashed", color = "blue") +
      geom_vline(xintercept = input$major_threshold_x, linetype = "dashed", color = "blue") +
      theme(legend.position = "right")
    # Conditionally add gene names
    if(input$show_gene_names) {
      p <- p + geom_text(aes(label = if_else(condition %in% c("minor", "major"), as.character(df[[1]]), NA_character_)),
                         vjust = 1, hjust = 1, size = 4, check_overlap = TRUE)
    }
    print(p)
  })
  
  # Displaying the list of minor candidate genes
  output$listMinor <- renderText({
    df <- df_reactive()
    genesRed <- df %>% filter(condition == "minor") %>% pull(1)
    if (length(genesRed) == 0) {
      return("No minor candidates found.")
    }
    paste("Candidate genes (Minor):", paste(genesRed, collapse = ", "))
  })
  # Displaying the list of major candidate genes
  output$listMajor <- renderText({
    df <- df_reactive()
    genesBlue <- df %>% filter(condition == "major") %>% pull(1)
    if (length(genesBlue) == 0) {
      return("No major candidates found.")
    }
    paste("Candidate genes (Major):", paste(genesBlue, collapse = ", "))
  })
}


shinyApp(ui = ui, server = server)

#read_excel("nexus/nexus-drug-screen/combined_data_2.xlsx",sheet=1)
