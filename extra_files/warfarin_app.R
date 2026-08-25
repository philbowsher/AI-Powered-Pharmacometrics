# Load required libraries
library(shiny)
library(shinydashboard)
library(plotly)
library(DT)
library(tidyverse)
library(nlmixr2data)
library(patchwork)

# Prepare data inside the app
data(warfarin)

# Create reproducible ADA status assignment
set.seed(2026)
unique_ids <- unique(warfarin$id)

id_mapping <- tibble(
  id = unique_ids,
  ADA_STATUS = sample(
    c("Positive", "Negative"), 
    size = length(unique_ids), 
    replace = TRUE, 
    prob = c(0.20, 0.80)
  )
)

# Create the final dataset with ADA impact
warfarin_app_data <- warfarin |> 
  left_join(id_mapping, by = "id") |> 
  mutate(
    # Apply realistic pharmacokinetic decay for ADA-positive patients
    DV_original = dv,
    dv = case_when(
      dvid == "cp" & ADA_STATUS == "Positive" ~ 
        dv * exp(-0.1 * 3.0 * time),  # 3x higher elimination rate
      TRUE ~ dv
    )
  )

# Create summary statistics for the app
ada_summary <- id_mapping |> 
  count(ADA_STATUS) |> 
  mutate(percentage = round(n / sum(n) * 100, 1))

# Create forest plot data
forest_data <- tibble(
  Variable = c("Weight (100kg vs 70kg)", "Age (70y vs 40y)", "ADA Positive Status"),
  Impact = c(1.25, 1.10, 3.0),
  Lower = c(1.05, 0.95, 2.1),
  Upper = c(1.48, 1.27, 4.3),
  Category = c("Patient Characteristics", "Patient Characteristics", "Immunogenicity"),
  Clinical_Significance = c("Moderate", "Minimal", "High")
)

# UI
ui <- dashboardPage(
  dashboardHeader(title = "Warfarin PK/PD Analysis"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Overview", tabName = "overview", icon = icon("chart-line")),
      menuItem("PK/PD Plots", tabName = "plots", icon = icon("chart-area")),
      menuItem("Patient Data", tabName = "data", icon = icon("table")),
      menuItem("Forest Plot", tabName = "forest", icon = icon("tree"))
    )
  ),
  
  dashboardBody(
    tabItems(
      # Overview Tab
      tabItem(tabName = "overview",
        fluidRow(
          valueBoxOutput("total_patients"),
          valueBoxOutput("ada_positive"),
          valueBoxOutput("observations")
        ),
        fluidRow(
          box(
            title = "Study Overview", status = "primary", solidHeader = TRUE, width = 12,
            h4("Warfarin PK/PD Analysis with ADA Impact"),
            p("This analysis explores the impact of Anti-Drug Antibodies (ADA) on warfarin pharmacokinetics and pharmacodynamics."),
            p("Key findings:"),
            tags$ul(
              tags$li("ADA-positive patients show 3x higher drug clearance"),
              tags$li("Despite lower concentrations, PD effects are preserved"),
              tags$li("Demonstrates PK/PD disconnect in immunogenicity")
            )
          )
        )
      ),
      
      # PK/PD Plots Tab
      tabItem(tabName = "plots",
        fluidRow(
          box(
            title = "Plot Controls", status = "primary", solidHeader = TRUE, width = 3,
            selectInput("ada_filter", "ADA Status Filter:",
                       choices = c("All", "Positive", "Negative"),
                       selected = "All"),
            sliderInput("time_range", "Time Range (hours):",
                       min = 0, max = 120, value = c(0, 120), step = 12),
            checkboxInput("log_scale", "Log Scale for PK", value = FALSE)
          ),
          box(
            title = "PK/PD Plots", status = "primary", solidHeader = TRUE, width = 9,
            plotlyOutput("pkpd_plot", height = "600px")
          )
        )
      ),
      
      # Patient Data Tab
      tabItem(tabName = "data",
        fluidRow(
          box(
            title = "Patient Data Explorer", status = "primary", solidHeader = TRUE, width = 12,
            DT::dataTableOutput("patient_table")
          )
        )
      ),
      
      # Forest Plot Tab
      tabItem(tabName = "forest",
        fluidRow(
          box(
            title = "Clinical Impact Factors", status = "primary", solidHeader = TRUE, width = 12,
            plotOutput("forest_plot", height = "400px"),
            p("Forest plot showing the relative impact of different factors on drug clearance.")
          )
        )
      )
    )
  )
)

# Server
server <- function(input, output) {
  
  # Reactive data filtering
  filtered_data <- reactive({
    data <- warfarin_app_data
    
    if(input$ada_filter != "All") {
      data <- data |> filter(ADA_STATUS == input$ada_filter)
    }
    
    data |> filter(time >= input$time_range[1] & time <= input$time_range[2])
  })
  
  # Value boxes
  output$total_patients <- renderValueBox({
    valueBox(
      value = length(unique(warfarin_app_data$id)),
      subtitle = "Total Patients",
      icon = icon("users"),
      color = "blue"
    )
  })
  
  output$ada_positive <- renderValueBox({
    ada_pos <- ada_summary |> filter(ADA_STATUS == "Positive") |> pull(percentage)
    valueBox(
      value = paste0(ada_pos, "%"),
      subtitle = "ADA Positive",
      icon = icon("exclamation-triangle"),
      color = "yellow"
    )
  })
  
  output$observations <- renderValueBox({
    valueBox(
      value = nrow(warfarin_app_data),
      subtitle = "Total Observations",
      icon = icon("chart-line"),
      color = "green"
    )
  })
  
  # PK/PD Plots
  output$pkpd_plot <- renderPlotly({
    data <- filtered_data()
    
    # PK Plot
    pk_data <- data |> filter(dvid == "cp")
    pk_plot <- ggplot(pk_data, aes(x = time, y = dv, color = ADA_STATUS, group = id)) +
      geom_line(alpha = 0.6) +
      geom_point(alpha = 0.8, size = 1) +
      labs(title = "Pharmacokinetics (Drug Concentration)",
           x = "Time (hours)", y = "Concentration (mg/L)",
           color = "ADA Status") +
      theme_minimal() +
      scale_color_manual(values = c("Positive" = "#e74c3c", "Negative" = "#3498db"))
    
    if(input$log_scale) {
      pk_plot <- pk_plot + scale_y_log10()
    }
    
    # PD Plot
    pd_data <- data |> filter(dvid == "pca")
    pd_plot <- ggplot(pd_data, aes(x = time, y = dv, color = ADA_STATUS, group = id)) +
      geom_line(alpha = 0.6) +
      geom_point(alpha = 0.8, size = 1) +
      labs(title = "Pharmacodynamics (Prothrombin Complex Activity)",
           x = "Time (hours)", y = "PCA (%)",
           color = "ADA Status") +
      theme_minimal() +
      scale_color_manual(values = c("Positive" = "#e74c3c", "Negative" = "#3498db"))
    
    # Combine plots
    combined <- pk_plot / pd_plot
    ggplotly(combined, height = 600)
  })
  
  # Patient Data Table
  output$patient_table <- DT::renderDataTable({
    warfarin_app_data |>
      select(id, time, dvid, dv, amt, wt, age, sex, ADA_STATUS) |>
      arrange(id, time)
  }, options = list(pageLength = 15, scrollX = TRUE))
  
  # Forest Plot
  output$forest_plot <- renderPlot({
    ggplot(forest_data, aes(x = Impact, y = reorder(Variable, Impact))) +
      geom_point(aes(color = Clinical_Significance), size = 4) +
      geom_errorbar(aes(xmin = Lower, xmax = Upper, color = Clinical_Significance), 
                    width = 0.2, orientation = "y") +
      geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.7) +
      scale_color_manual(values = c("High" = "#e74c3c", "Moderate" = "#f39c12", "Minimal" = "#27ae60")) +
      labs(title = "Clinical Impact on Drug Clearance",
           x = "Fold Change in Clearance", y = "Factor",
           color = "Clinical Significance") +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 12))
  })
}

# Run the app
shinyApp(ui = ui, server = server)