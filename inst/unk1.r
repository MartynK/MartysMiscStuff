library(shiny)
library(ggplot2)
library(Hotelling)  # Ensure this package is installed

ui <- fluidPage(
  titlePanel("Comparing Test A and Test B"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("diff_mean",
                  "Difference in Means (B - A)",
                  min = -2, max = 2, value = 0, step = 0.1),
      sliderInput("sd_rat",
                  "Multiplicative Ratio for SD (B vs A)",
                  min = 0.1, max = 10, value = 1, step = 0.1),
      helpText("Note: Test A has mean = 1 and SD = 1.",
               "Test B has mean = 1 + (diff_mean) and SD = 1 * (sd_rat).",
               "Red triangle: true parameter difference.",
               "Simulated replicate points will be blue if the two-sample Hotelling T-test is significant (p <= 0.05) and grey otherwise.")
    ),
    mainPanel(
      plotOutput("scatterPlot")
    )
  )
)

server <- function(input, output) {
  output$scatterPlot <- renderPlot({
    # Parameters for Test A (fixed)
    mean_A <- 1
    sd_A <- 1

    # For Test B: Only multiplicative scaling on the SD is applied.
    mean_B <- 1 + input$diff_mean
    sd_B <- sd_A * input$sd_rat

    # Compute the "true" differences (Test B - Test A)
    true_diff_mean <- mean_B - mean_A   # This equals input$diff_mean
    true_diff_sd   <- sd_B - sd_A         # This equals (sd_rat - 1)

    # Simulation settings
    nRep <- 40  # number of replicate experiments
    nObs <- 40  # number of observations per test in each experiment

    # Data frames to store simulation results and summary statistics
    sim_results <- data.frame(diff_mean = numeric(nRep),
                              diff_sd = numeric(nRep))
    statsA <- data.frame(mean = numeric(nRep), sd = numeric(nRep))
    statsB <- data.frame(mean = numeric(nRep), sd = numeric(nRep))

    # Simulate replicate experiments
    for(i in 1:nRep) {
      # Generate observations for each test
      obs_A <- rnorm(nObs, mean = mean_A, sd = sd_A)
      obs_B <- rnorm(nObs, mean = mean_B, sd = sd_B)

      # Compute summary statistics for each test
      statsA$mean[i] <- mean(obs_A)
      statsA$sd[i]   <- sd(obs_A)
      statsB$mean[i] <- mean(obs_B)
      statsB$sd[i]   <- sd(obs_B)

      # Compute differences (Test B - Test A) for plotting
      sim_results$diff_mean[i] <- statsB$mean[i] - statsA$mean[i]
      sim_results$diff_sd[i]   <- statsB$sd[i]   - statsA$sd[i]
    }

    # Perform a two-sample Hotelling T² test on the bivariate (mean, SD) outcomes
    test_result <- hotelling.test(statsA, statsB)
    p_value <- test_result$pval

    # Set color: blue if significant (p <= 0.05), grey otherwise.
    pointColor <- ifelse(p_value > 0.05, "grey", "blue")

    # Create the scatter plot
    ggplot(sim_results, aes(x = diff_sd, y = diff_mean)) +
      geom_point(color = pointColor, size = 3) +
      geom_point(aes(x = true_diff_sd, y = true_diff_mean),
                 color = "red", size = 5, shape = 17) +
      labs(x = "Difference in SD (Test B - Test A)",
           y = "Difference in Mean (Test B - Test A)",
           title = "Simulated Differences from 40 Experiments",
           subtitle = paste("Hotelling T-test p-value:", round(p_value, 3),
                            "| Red triangle: true parameter difference")) +
      theme_minimal() +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_vline(xintercept = 0, linetype = "dashed") +
      geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
      xlim(min(sim_results$diff_sd, true_diff_sd, -0.5) - 0.5,
           max(sim_results$diff_sd, true_diff_sd, 0) + 0.5) +
      ylim(min(sim_results$diff_mean, true_diff_mean, -0.5) - 0.5,
           max(sim_results$diff_mean, true_diff_mean, 0) + 0.5)
  })
}

shinyApp(ui = ui, server = server)
