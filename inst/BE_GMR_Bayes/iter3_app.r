library(shiny)
library(rstan)
library(ggplot2)
library(dplyr)

ui <- fluidPage(
  titlePanel("GMR Bayesian Dashboard"),

  sidebarLayout(
    sidebarPanel(
      numericInput("sample_size", "Sample Size:", value = 24, min = 3, max = 90),
      numericInput("true_gmr", "True GMR:", value = 1.05, min = 0.7, max = 1/0.7, step=0.01),
      numericInput("cv", "CV%:", value = 36, min = 1, max = 100, step=1),
      numericInput("seed", "Random Seed:", value = 42, min = 1, max = 100, step=1),
      selectInput("prior_choice", "Prior Distribution:",
                  choices = c("Uniform (0.7, 1.43)" = "uniform",
                              "Normal (centered normal)" = "extreme",
                              "Wide Beta (0,2)" = "wide_beta"
                  )),
      selectInput("credence_level", "Credence Level:",
                  choices = c("95%" = 0.95, "90%" = 0.90, "80%" = 0.80),
                  selected = 0.90),
      actionButton("execute", "Execute")
    ),

    mainPanel(
      plotOutput("ci_distribution_plot"),
      plotOutput("prior_posterior_plot")
    )
  )
)

server <- function(input, output) {
  stan_model_reactive <- reactiveVal(NULL)

  observeEvent(input$prior_choice, {
    os_type <- Sys.info()[["sysname"]]
    model_file <- switch(input$prior_choice,
                         "uniform" = ifelse(os_type == "Windows",
                                            "stan_uniform.rdata",
                                            "uniform_ubuntu.rdata"),
                         "extreme" = ifelse(os_type == "Windows",
                                            "stan_extreme.rdata",
                                            "extreme_ubuntu.rdata"),
                         "wide_beta" = ifelse(os_type == "Windows",
                                              "stan_wide_beta.rdata",
                                              "wide_beta_ubuntu.rdata"))

    load(file.path("comp_models", model_file))  # Load precompiled model from file
    stan_model_reactive(stan_model)
  }, ignoreNULL = FALSE)

  results <- eventReactive(input$execute, {
    N <- input$sample_size
    true_mu <- log(input$true_gmr)
    true_sigma <- sqrt(log(1 + (input$cv / 100)^2))



    set.seed(input$seed)
    simulated_data <- rnorm(N, mean = true_mu, sd = true_sigma)

    fit <- sampling(stan_model_reactive(),
                    data = list(N = N, y = simulated_data),
                    iter = 2000, chains = 4, refresh=0, init_r=0.5)

    posterior_samples <- extract(fit)$GMR

    credence <- as.numeric(input$credence_level)
    alpha <- (1 - credence) / 2
    ci <- quantile(posterior_samples, c(alpha, 1 - alpha))

    # Iterating over seeds
    ci_lower <- numeric(100)
    ci_upper <- numeric(100)
    medians <- numeric(100)

    withProgress(message = 'Running iterations', value = 0, {
      for (i in 1:100) {
        set.seed(i)
        sim_data_iter <- rnorm(N, mean = true_mu, sd = true_sigma)
        fit_iter <- sampling(stan_model_reactive(),
                             data = list(N = N, y = sim_data_iter),
                             iter = 2000, chains = 4, refresh=0)
        posterior_iter <- extract(fit_iter)$GMR
        ci_iter <- quantile(posterior_iter, c(0.025, 0.975))
        ci_lower[i] <- ci_iter[1]
        ci_upper[i] <- ci_iter[2]
        medians[i] <- median(posterior_iter)
        incProgress(1/100)
      }
    })

    list(posterior_samples = posterior_samples,
         ci = ci,
         prior = input$prior_choice,
         ci_lower = ci_lower,
         ci_upper = ci_upper,
         median = medians,
         simulated_data = simulated_data)
  })

  output$prior_posterior_plot <- renderPlot({
    req(results())

    posterior_samples <- results()$posterior_samples
    ci <- results()$ci
    prior_type <- results()$prior

    df_prior <- data.frame(GMR=seq(0.5, 1.5, length.out=200))
    df_prior$density <- if(prior_type == "uniform") {
        dunif(df_prior$GMR, 0.7, 1.43)
      } else if(prior_type == "extreme") {
        dnorm(log(df_prior$GMR), log(1), 0.5)#/df_prior$GMR
      } else {
        dbeta(df_prior$GMR/2, 5, 5)/2
      }

    ggplot() +
      geom_area(data=df_prior, aes(x=GMR, y=density), fill="blue", alpha=0.2) +
      geom_density(data=data.frame(GMR=posterior_samples), aes(x=GMR), fill="red", alpha=0.4) +
      geom_vline(xintercept=ci, linetype="dashed") +
      #geom_vline(xintercept = mean(exp(results()$simulated_data))) +
      scale_x_continuous(limits = c(0.5, 1.5), breaks = c(seq(0.5,1.5,length.out = 5),0.9,1.11)) +
      labs(title="Example Prior (blue) and Posterior (red) Distributions for one simulated study\nwith 90% CI", x="GMR", y="Density") + theme_minimal()
  })

  output$ci_distribution_plot <- renderPlot({
    req(results())

    df_ci <- data.frame(lower=results()$ci_lower,
                        upper=results()$ci_upper,
                        median=results()$median)

    credence <- as.numeric(input$credence_level)
    alpha <- (1 - credence) / 2

    ggplot() +
      geom_density(data=df_ci, aes(x=lower), fill="blue", alpha=0.4) +
      geom_density(data=df_ci, aes(x=upper), fill="green", alpha=0.4) +
      geom_density(data=df_ci, aes(x=median), fill="red", alpha=0.4) +
      geom_vline(xintercept = quantile(df_ci$median, probs = c(alpha, 1 - alpha)),
                 linetype = "dashed") +
      geom_vline(xintercept = mean(results()$simulated_data)) +
      scale_x_continuous(limits = c(0.5, 1.5), breaks = c(seq(0.5,1.5,length.out = 5),0.9,1.11)) +
      labs(title="Distribution of 95% CI Limits and Median GMR over 100 Seeds",
           x="GMR", y="Density") +
      theme_minimal() +
      annotate("text", x=c(0.75,1.25,1), y=Inf, label=c("Lower CI","Upper CI","Median"),
               vjust=2, color=c("blue","green","red"))+
      annotate("text", x= quantile(df_ci$median, probs=c(alpha, 1 - alpha)),
               vjust = 3, y=Inf,
               label= round(quantile(df_ci$median, probs=c(alpha, 1 - alpha)), digits=3))
  })
}

shinyApp(ui, server)
