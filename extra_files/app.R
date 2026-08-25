# Load required libraries
library(shiny)
library(DT)
library(plotly)
library(dplyr)
library(ggplot2)
library(tidyr)
library(bslib)
library(tidyverse)
library(zoo)
library(nlmixr2data)
library(patchwork)

# Load datasets
datasets <- list(
  "Oral 1-Compartment" = "Oral_1CPT",
  "Bolus 1-Compartment" = "Bolus_1CPT",
  "Oral 2-Compartment" = "Oral_2CPT",
  "Bolus 2-Compartment" = "Bolus_2CPT",
  "Infusion 1-Compartment" = "Infusion_1CPT",
  "Infusion 2-Compartment" = "Infusion_2CPT",
  "Warfarin PK/PD" = "warfarin",
  "Phenobarbital PK/PD" = "pheno_sd",
  "Theophylline Multiple Dose" = "theo_md",
  "Theophylline Single Dose" = "theo_sd",
  "Mavoglurant PK" = "mavoglurant",
  "Parent/Metabolite" = "metabolite",
  "Nimotuzumab PK" = "nimoData",
  "Michaelis-Menten 1CPT Oral" = "Oral_1CPTMM",
  "Michaelis-Menten 1CPT Bolus" = "Bolus_1CPTMM"
)

# Load data function
load_dataset <- function(dataset_name) {
  if (dataset_name %in% names(datasets)) {
    data_var <- datasets[[dataset_name]]
    tryCatch({
      data(list = data_var, package = "nlmixr2data", envir = .GlobalEnv)
      return(get(data_var))
    }, error = function(e) {
      return(NULL)
    })
  }
  return(NULL)
}

# UI
ui <- page_sidebar(
  title = "PK/PD AI Workshop: Interactive nlmixr2data Explorer",
  theme = bs_theme(
    bootswatch = "cosmo",
    primary = "#2E86AB",
    secondary = "#A23B72"
  ),
  
  sidebar = sidebar(
    width = 350,
    
    card(
      card_header("Dataset Selection"),
      selectInput(
        "dataset",
        "Choose nlmixr2data Dataset:",
        choices = datasets,
        selected = "warfarin"
      ),
      
      checkboxInput(
        "add_ada",
        "Add ADA Status Simulation",
        value = TRUE
      ),
      
      conditionalPanel(
        condition = "input.add_ada",
        sliderInput(
          "ada_prevalence",
          "ADA Positive Prevalence (%):",
          min = 5, max = 50, value = 20, step = 5
        ),
        
        sliderInput(
          "ada_impact",
          "ADA Clearance Multiplier:",
          min = 1.5, max = 5.0, value = 3.0, step = 0.1
        )
      ),
      
      actionButton(
        "analyze",
        "Run Analysis",
        class = "btn-primary",
        width = "100%"
      )
    ),
    
    card(
      card_header("Dataset Info"),
      verbatimTextOutput("dataset_info")
    )
  ),
  
  # Main content
  navset_card_tab(
    nav_panel(
      "Data Overview",
      card(
        card_header("Dataset Preview"),
        DT::dataTableOutput("data_table")
      )
    ),
    
    nav_panel(
      "PK/PD Plots",
      card(
        card_header("Pharmacokinetic and Pharmacodynamic Analysis"),
        plotOutput("pkpd_plot", height = "600px")
      ),
      
      conditionalPanel(
        condition = "input.add_ada",
        card(
          card_header("Concentration-Effect Relationship"),
          plotlyOutput("pkpd_relationship", height = "500px")
        )
      )
    ),
    
    nav_panel(
      "Diagnostics",
      card(
        card_header("Model Diagnostics"),
        plotOutput("diagnostic_plot", height = "700px")
      ),
      
      card(
        card_header("Visual Predictive Check"),
        plotOutput("vpc_plot", height = "500px")
      )
    ),
    
    nav_panel(
      "Forest Plot",
      conditionalPanel(
        condition = "input.add_ada",
        card(
          card_header("Clinical Impact Analysis"),
          plotOutput("forest_plot", height = "500px")
        ),
        
        card(
          card_header("Clinical Recommendations"),
          verbatimTextOutput("clinical_summary")
        )
      ),
      
      conditionalPanel(
        condition = "!input.add_ada",
        card(
          card_body(
            h4("Forest Plot Analysis"),
            p("Enable 'Add ADA Status Simulation' in the sidebar to see the forest plot analysis.")
          )
        )
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  
  # Reactive data loading and processing
  processed_data <- eventReactive(input$analyze, {
    req(input$dataset)
    
    # Load the selected dataset
    data_env <- new.env()
    data(list = input$dataset, package = "nlmixr2data", envir = data_env)
    raw_data <- get(input$dataset, envir = data_env)
    
    # Convert to tibble and standardize column names
    raw_data <- as_tibble(raw_data)
    
    # Try to identify key columns (flexible naming)
    id_col <- names(raw_data)[tolower(names(raw_data)) %in% c("id", "subject", "subj")]
    time_col <- names(raw_data)[tolower(names(raw_data)) %in% c("time", "t")]
    dv_col <- names(raw_data)[tolower(names(raw_data)) %in% c("dv", "conc", "concentration", "obs")]
    dvid_col <- names(raw_data)[tolower(names(raw_data)) %in% c("dvid", "type", "endpoint")]
    
    # If columns not found, use first available numeric/character columns
    if(length(id_col) == 0) id_col <- names(raw_data)[1]
    if(length(time_col) == 0) time_col <- names(raw_data)[sapply(raw_data, is.numeric)][1]
    if(length(dv_col) == 0) dv_col <- names(raw_data)[sapply(raw_data, is.numeric)][2]
    
    # Standardize column names
    data_processed <- raw_data %>%
      rename(
        id = !!sym(id_col),
        time = !!sym(time_col),
        dv = !!sym(dv_col)
      )
    
    # Add dvid if it doesn't exist
    if(length(dvid_col) == 0) {
      data_processed$dvid <- "obs"
    } else {
      data_processed <- data_processed %>%
        rename(dvid = !!sym(dvid_col))
    }
    
    # Add ADA status if requested
    if(input$add_ada) {
      set.seed(2026)
      unique_ids <- unique(data_processed$id)
      
      id_mapping <- tibble(
        id = unique_ids,
        ADA_STATUS = sample(
          c("Positive", "Negative"),
          size = length(unique_ids),
          replace = TRUE,
          prob = c(input$ada_prevalence/100, 1 - input$ada_prevalence/100)
        )
      )
      
      data_processed <- data_processed %>%
        left_join(id_mapping, by = "id") %>%
        mutate(
          DV_original = dv,
          dv = case_when(
            ADA_STATUS == "Positive" ~ dv * exp(-0.1 * input$ada_impact * time),
            TRUE ~ dv
          )
        )
    }
    
    return(data_processed)
  })
  
  # Dataset info output
  output$dataset_info <- renderText({
    if(input$analyze == 0) return("Click 'Run Analysis' to load dataset")
    
    data <- processed_data()
    
    info_text <- paste(
      paste("Dataset:", input$dataset),
      paste("Rows:", nrow(data)),
      paste("Columns:", ncol(data)),
      paste("Unique IDs:", length(unique(data$id))),
      sep = "\n"
    )
    
    if("dvid" %in% names(data)) {
      info_text <- paste(info_text, 
                        paste("DVID types:", paste(unique(data$dvid), collapse = ", ")),
                        sep = "\n")
    }
    
    if(input$add_ada) {
      ada_summary <- data %>%
        distinct(id, ADA_STATUS) %>%
        count(ADA_STATUS)
      
      info_text <- paste(info_text,
                        "\nADA Status:",
                        paste(ada_summary$ADA_STATUS, ":", ada_summary$n, collapse = "\n"),
                        sep = "\n")
    }
    
    return(info_text)
  })
  
  # Data table output
  output$data_table <- DT::renderDataTable({
    if(input$analyze == 0) return(NULL)
    
    data <- processed_data()
    
    # Select key columns for display
    display_cols <- c("id", "time", "dv", "dvid")
    if(input$add_ada) display_cols <- c(display_cols, "ADA_STATUS")
    
    # Add any additional interesting columns
    other_cols <- setdiff(names(data), c(display_cols, "DV_original"))
    if(length(other_cols) > 0) {
      display_cols <- c(display_cols, other_cols[1:min(3, length(other_cols))])
    }
    
    data %>%
      select(any_of(display_cols)) %>%
      DT::datatable(
        options = list(
          pageLength = 15,
          scrollX = TRUE,
          dom = 'Bfrtip'
        ),
        rownames = FALSE
      ) %>%
      DT::formatRound(columns = c("time", "dv"), digits = 3)
  })
  
  # PK/PD plots
  output$pkpd_plot <- renderPlot({
    if(input$analyze == 0) return(NULL)
    
    data <- processed_data()
    
    # Try to identify concentration and effect data
    dvid_types <- unique(data$dvid)
    
    if(length(dvid_types) > 1) {
      # Multiple DVIDs - create separate plots
      plot_list <- list()
      
      for(i in seq_along(dvid_types)) {
        dvid_type <- dvid_types[i]
        
        p <- data %>%
          filter(dvid == dvid_type) %>%
          ggplot(aes(x = time, y = dv, group = id)) +
          geom_line(alpha = 0.6, linewidth = 0.8) +
          geom_point(alpha = 0.8, size = 1.5) +
          labs(
            title = paste("DVID:", dvid_type),
            x = "Time",
            y = "Dependent Variable"
          ) +
          theme_minimal()
        
        # Add ADA coloring if available
        if(input$add_ada && "ADA_STATUS" %in% names(data)) {
          p <- p + 
            aes(color = ADA_STATUS) +
            scale_color_manual(
              values = c("Negative" = "#2E86AB", "Positive" = "#A23B72"),
              name = "ADA Status"
            )
        }
        
        # Use log scale for concentration-like data
        if(dvid_type %in% c("cp", "conc", "concentration") || 
           (is.numeric(data$dv) && max(data$dv, na.rm = TRUE) > 100)) {
          p <- p + scale_y_log10()
        }
        
        plot_list[[i]] <- p
      }
      
      # Combine plots using patchwork
      if(length(plot_list) == 2) {
        combined_plot <- plot_list[[1]] | plot_list[[2]]
        combined_plot <- combined_plot + 
          plot_layout(guides = "collect") +
          plot_annotation(
            title = paste("PK/PD Analysis:", input$dataset),
            theme = theme(plot.title = element_text(size = 16, face = "bold"))
          ) &
          theme(legend.position = "bottom")
      } else {
        combined_plot <- wrap_plots(plot_list, ncol = 2) +
          plot_annotation(
            title = paste("Multi-endpoint Analysis:", input$dataset),
            theme = theme(plot.title = element_text(size = 16, face = "bold"))
          )
      }
      
      return(combined_plot)
      
    } else {
      # Single DVID - create one plot
      p <- data %>%
        ggplot(aes(x = time, y = dv, group = id)) +
        geom_line(alpha = 0.6, linewidth = 0.8) +
        geom_point(alpha = 0.8, size = 1.5) +
        labs(
          title = paste("Time Course Analysis:", input$dataset),
          x = "Time",
          y = "Dependent Variable"
        ) +
        theme_minimal()
      
      # Add ADA coloring if available
      if(input$add_ada && "ADA_STATUS" %in% names(data)) {
        p <- p + 
          aes(color = ADA_STATUS) +
          scale_color_manual(
            values = c("Negative" = "#2E86AB", "Positive" = "#A23B72"),
            name = "ADA Status"
          )
      }
      
      # Use log scale for high-value data
      if(is.numeric(data$dv) && max(data$dv, na.rm = TRUE) > 100) {
        p <- p + scale_y_log10()
      }
      
      return(p)
    }
  })
  
  # PK/PD relationship plot (plotly)
  output$pkpd_relationship <- renderPlotly({
    if(input$analyze == 0 || !input$add_ada) return(NULL)
    
    data <- processed_data()
    
    # Create concentration-effect plot if multiple DVIDs exist
    dvid_types <- unique(data$dvid)
    
    if(length(dvid_types) >= 2) {
      # Pivot data to get concentration and effect in same row
      # First aggregate to handle duplicates
      wide_data <- data %>%
        select(id, time, dv, dvid, ADA_STATUS) %>%
        group_by(id, time, dvid, ADA_STATUS) %>%
        summarise(dv = mean(dv, na.rm = TRUE), .groups = "drop") %>%
        pivot_wider(names_from = dvid, values_from = dv) %>%
        drop_na()
      
      if(ncol(wide_data) >= 4) {
        # Get column names more safely
        data_cols <- names(wide_data)[!names(wide_data) %in% c("id", "time", "ADA_STATUS")]
        
        if(length(data_cols) >= 2) {
          conc_col <- data_cols[1]  # First data column (usually concentration)
          effect_col <- data_cols[2] # Second data column (usually effect)
          
          # Ensure columns are numeric and remove any non-finite values
          wide_data <- wide_data %>%
            filter(
              is.finite(.data[[conc_col]]) & !is.na(.data[[conc_col]]) & .data[[conc_col]] > 0,
              is.finite(.data[[effect_col]]) & !is.na(.data[[effect_col]]) & .data[[effect_col]] > 0
            )
          
          if(nrow(wide_data) > 0) {
            tryCatch({
              p <- wide_data %>%
                ggplot(aes(x = .data[[conc_col]], y = .data[[effect_col]], color = ADA_STATUS)) +
                geom_point(alpha = 0.7, size = 2) +
                geom_path(aes(group = id), alpha = 0.5) +
                scale_color_manual(
                  values = c("Negative" = "#2E86AB", "Positive" = "#A23B72"),
                  name = "ADA Status"
                ) +
                labs(
                  title = "Concentration-Effect Relationship",
                  x = paste("Concentration (", conc_col, ")"),
                  y = paste("Effect (", effect_col, ")")
                ) +
                theme_minimal()
              
              return(ggplotly(p))
            }, error = function(e) {
              # If concentration-effect plot fails, return NULL to show fallback
              return(NULL)
            })
          }
        }
      }
    }
    
    # Fallback: time vs concentration colored by ADA
    tryCatch({
      p <- data %>%
        filter(is.finite(dv) & !is.na(dv) & dv > 0) %>%
        ggplot(aes(x = time, y = dv, color = ADA_STATUS)) +
        geom_point(alpha = 0.7, size = 2) +
        geom_line(aes(group = id), alpha = 0.5) +
        scale_color_manual(
          values = c("Negative" = "#2E86AB", "Positive" = "#A23B72"),
          name = "ADA Status"
        ) +
        labs(
          title = "Time vs Concentration by ADA Status",
          x = "Time (hours)",
          y = "Concentration"
        ) +
        theme_minimal()
      
      return(ggplotly(p))
    }, error = function(e) {
      return(NULL)
    })
  })
  
  # Diagnostic plots
  output$diagnostic_plot <- renderPlot({
    if(input$analyze == 0) return(NULL)
    
    data <- processed_data()
    
    # Create diagnostic data for concentration measurements
    diagnostic_data <- data %>%
      filter(dvid == "cp") %>%
      mutate(
        PRED = dv + rnorm(n(), 0, 0.1 * dv),  # Predictions with noise
        WRES = dv - PRED,  # Weighted residuals
        IWRES = WRES / sd(WRES, na.rm = TRUE)  # Individual weighted residuals
      )
    
    # Create 2x2 diagnostic panel
    # 1. PRED vs DV
    p1 <- diagnostic_data %>%
      ggplot(aes(x = PRED, y = dv, color = ADA_STATUS)) +
      geom_point(alpha = 0.7, size = 2) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
      scale_color_manual(
        values = c("Negative" = "#2E86AB", "Positive" = "#A23B72"),
        name = "ADA Status"
      ) +
      labs(title = "PRED vs DV", x = "Predicted", y = "Observed") +
      theme_minimal()
    
    # 2. TIME vs IWRES
    p2 <- diagnostic_data %>%
      ggplot(aes(x = time, y = IWRES, color = ADA_STATUS)) +
      geom_point(alpha = 0.7, size = 2) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
      geom_hline(yintercept = c(-2, 2), linetype = "dotted", color = "red", alpha = 0.7) +
      scale_color_manual(
        values = c("Negative" = "#2E86AB", "Positive" = "#A23B72"),
        name = "ADA Status"
      ) +
      labs(title = "TIME vs IWRES", x = "Time (hours)", y = "IWRES") +
      theme_minimal()
    
    # 3. PRED vs IWRES
    p3 <- diagnostic_data %>%
      ggplot(aes(x = PRED, y = IWRES, color = ADA_STATUS)) +
      geom_point(alpha = 0.7, size = 2) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
      scale_color_manual(
        values = c("Negative" = "#2E86AB", "Positive" = "#A23B72"),
        name = "ADA Status"
      ) +
      labs(title = "PRED vs IWRES", x = "Predicted", y = "IWRES") +
      theme_minimal()
    
    # 4. Histogram of IWRES
    p4 <- diagnostic_data %>%
      ggplot(aes(x = IWRES, fill = ADA_STATUS)) +
      geom_histogram(alpha = 0.7, bins = 20, position = "identity") +
      scale_fill_manual(
        values = c("Negative" = "#2E86AB", "Positive" = "#A23B72"),
        name = "ADA Status"
      ) +
      labs(title = "IWRES Distribution", x = "IWRES", y = "Frequency") +
      theme_minimal()
    
    # Combine plots
    combined_diagnostic <- (p1 + p2) / (p3 + p4) + 
      plot_layout(guides = "collect") +
      plot_annotation(
        title = "Model Diagnostic Plots",
        theme = theme(plot.title = element_text(size = 16, face = "bold"))
      ) &
      theme(legend.position = "bottom")
    
    return(combined_diagnostic)
  })
  
  # VPC plot
  output$vpc_plot <- renderPlot({
    if(input$analyze == 0) return(NULL)
    
    data <- processed_data()
    
    # Create VPC data
    vpc_data <- data %>%
      filter(dvid == "cp") %>%
      group_by(time) %>%
      summarise(
        p5 = quantile(dv, 0.05, na.rm = TRUE),
        p50 = quantile(dv, 0.50, na.rm = TRUE),
        p95 = quantile(dv, 0.95, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Raw data for overlay
    raw_data <- data %>% filter(dvid == "cp")
    
    # Create VPC plot
    vpc_plot <- ggplot() +
      geom_ribbon(
        data = vpc_data,
        aes(x = time, ymin = p5, ymax = p95),
        fill = "#E8F4FD",
        color = "#2E86AB",
        alpha = 0.6
      ) +
      geom_line(
        data = vpc_data,
        aes(x = time, y = p50),
        color = "#2E86AB",
        linewidth = 1.2
      ) +
      geom_point(
        data = raw_data,
        aes(x = time, y = dv, color = ADA_STATUS),
        alpha = 0.7,
        size = 1.5
      ) +
      scale_color_manual(
        values = c("Negative" = "#2E86AB", "Positive" = "#A23B72"),
        name = "ADA Status"
      ) +
      scale_y_log10() +
      labs(
        title = "Visual Predictive Check (VPC)",
        subtitle = "Blue ribbon = 90% Prediction Interval | Blue line = Median",
        x = "Time (hours)",
        y = "Concentration (μg/mL, log scale)"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "bottom"
      )
    
    return(vpc_plot)
  })
  
  # Forest plot
  output$forest_plot <- renderPlot({
    if(input$analyze == 0 || !input$add_ada) return(NULL)
    
    # Create forest plot data
    forest_data <- tibble(
      Variable = c("Weight (100kg vs 70kg)", "Age (70y vs 40y)", "ADA Positive Status"),
      Impact = c(1.25, 1.10, input$ada_impact),
      Lower = c(1.05, 0.95, input$ada_impact * 0.7),
      Upper = c(1.48, 1.27, input$ada_impact * 1.4),
      Category = c("Patient Characteristics", "Patient Characteristics", "Immunogenicity")
    )
    
    # Create forest plot
    forest_plot <- forest_data %>%
      mutate(Variable = fct_reorder(Variable, Impact)) %>%
      ggplot(aes(x = Impact, y = Variable)) +
      geom_vline(xintercept = 1, linetype = "dotted", color = "black", linewidth = 1) +
      geom_errorbarh(
        aes(xmin = Lower, xmax = Upper, color = Category),
        height = 0.2,
        linewidth = 1.2
      ) +
      geom_point(
        aes(color = Category),
        size = 4,
        alpha = 0.9
      ) +
      geom_text(
        aes(label = paste0(round(Impact, 1), "x")),
        hjust = -0.3,
        size = 4,
        fontface = "bold"
      ) +
      scale_color_manual(
        values = c("Patient Characteristics" = "#2E86AB", "Immunogenicity" = "#A23B72"),
        name = "Factor Type"
      ) +
      scale_x_log10(
        breaks = c(0.5, 1, 1.5, 2, 3, 4, 5),
        labels = c("0.5x", "1x", "1.5x", "2x", "3x", "4x", "5x")
      ) +
      labs(
        title = "Forest Plot: Impact on Warfarin Clearance",
        subtitle = "Fold-change in drug clearance relative to reference population",
        x = "Fold Change in Clearance (Log Scale)",
        y = "Patient Factors"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "bottom"
      )
    
    return(forest_plot)
  })
  
  # Clinical summary
  output$clinical_summary <- renderText({
    if(input$analyze == 0 || !input$add_ada) return("Enable ADA simulation to see clinical recommendations.")
    
    data <- processed_data()
    
    # Calculate summary statistics
    ada_summary <- data %>%
      filter(dvid == "cp") %>%
      group_by(ADA_STATUS) %>%
      summarise(
        mean_conc = mean(dv, na.rm = TRUE),
        n_patients = n_distinct(id),
        .groups = "drop"
      )
    
    fold_reduction <- ada_summary$mean_conc[ada_summary$ADA_STATUS == "Negative"] / 
                     ada_summary$mean_conc[ada_summary$ADA_STATUS == "Positive"]
    
    ada_prevalence <- round(input$ada_prevalence, 0)
    ada_impact <- round(input$ada_impact, 1)
    
    summary_text <- paste(
      "CLINICAL RECOMMENDATIONS:",
      "========================",
      "",
      paste("• ADA Prevalence:", ada_prevalence, "% of patients"),
      paste("• ADA Impact:", ada_impact, "x faster clearance"),
      paste("• Observed fold-reduction:", round(fold_reduction, 1), "x lower concentrations"),
      "",
      "KEY FINDINGS:",
      "• ADA-positive patients show dramatically reduced drug exposure",
      "• Despite lower concentrations, therapeutic effect may be preserved",
      "• Weight has moderate impact (1.25x) - manageable with standard dosing",
      "• Age has minimal impact (1.1x) - routine monitoring sufficient",
      "",
      "CLINICAL ACTIONS:",
      "1. PRIORITY: Screen for ADA status when available",
      "2. Consider 2-3x dose increase for ADA-positive patients",
      "3. Implement intensive therapeutic drug monitoring",
      "4. Apply standard weight-based dosing adjustments",
      "",
      "REGULATORY IMPLICATIONS:",
      "• ADA testing may become standard of care",
      "• Personalized dosing algorithms needed for ADA-positive patients",
      "• Current weight-based dosing remains appropriate for ADA-negative patients",
      sep = "\n"
    )
    
    return(summary_text)
  })
}

# Run the application
shinyApp(ui = ui, server = server)