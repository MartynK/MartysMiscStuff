library(shiny)
library(rstan)
library(ggplot2)

ui <- fluidPage(
  titlePanel("Bioequivalence GMR Bayesian Dashboard"),

  sidebarLayout(
    sidebarPanel(
      numericInput("sample_size", "Sample Size:", value = 33, min = 3, max = 90),
      numericInput("true_gmr", "True GMR:", value = 1.05, min = 0.7, max = 1/0.7, step=0.01),
      numericInput("cv", "CV%:", value = 36, min = 1, max = 100, step=1),
      numericInput("seed", "Random Seed:", value = 42, min = 1, max = 100, step=1),
      selectInput("prior_choice", "Prior Distribution:",
                  choices = c("Uniform (0.7, 1.3)" = "uniform",
                              "Beta" = "extreme",
                              "Wide Beta (0,2)" = "wide_beta")),
      actionButton("execute", "Execute")
    ),

    mainPanel(
      plotOutput("prior_posterior_plot")
    )
  )
)

server <- function(input, output) {
  stan_model_reactive <- reactiveVal(NULL)

  observeEvent(input$prior_choice, {
    prior_code <- switch(input$prior_choice,
                         "uniform" = "mu ~ uniform(log(0.7), log(1.3));",
                         "extreme" = "mu ~ beta(2, 2);",
                         "wide_beta" = "mu ~ beta(5, 5);")

    stan_model_code <- paste0("
data {
  int<lower=1> N;
  vector[N] y;
}
parameters {
  real<lower=", ifelse(input$prior_choice == "wide_beta", "log(0)", "log(0.7)"), ", upper=", ifelse(input$prior_choice == "wide_beta", "log(2)", "log(1.3)"), "> mu;
  real<lower=0> sigma;
}
model {
  ", prior_code, "
  sigma ~ cauchy(0, 2.5);
  y ~ normal(mu, sigma);
}
generated quantities {
  real GMR = exp(mu);
}")

    stan_model_reactive(rstan::stan_model(model_code = stan_model_code))
  }, ignoreNULL = FALSE)

  results <- eventReactive(input$execute, {
    N <- input$sample_size
    true_mu <- log(input$true_gmr)
    true_sigma <- sqrt(log(1 + (input$cv / 100)^2))

    set.seed(input$seed)
    simulated_data <- rnorm(N, mean = true_mu, sd = true_sigma)

    fit <- sampling(stan_model_reactive(),
                    data = list(N = N, y = simulated_data),
                    iter = 1000, chains = 2, refresh=0, init_r=0.5)

    posterior_samples <- extract(fit)$GMR
    ci <- quantile(posterior_samples, c(0.025, 0.975))

    list(posterior_samples = posterior_samples, ci = ci, prior = input$prior_choice)
  })

  output$prior_posterior_plot <- renderPlot({
    req(results())

    posterior_samples <- results()$posterior_samples
    ci <- results()$ci
    prior_type <- results()$prior

    df_prior <- data.frame(GMR=seq(if(prior_type=="wide_beta") 0 else 0.7, if(prior_type=="wide_beta") 2 else 1.3, length.out=100))
    df_prior$density <- if(prior_type == "uniform") {
      dunif(df_prior$GMR, min=0.7, max=1.3)
    } else if(prior_type == "extreme") {
      dbeta((df_prior$GMR - 0.7)/0.6, 2, 2) / 0.6
    } else {
      dbeta(df_prior$GMR/2, 5, 5) / 2
    }

    df_posterior <- data.frame(GMR=posterior_samples)

    ggplot() +
      geom_area(data=df_prior, aes(x=GMR, y=density), fill="blue", alpha=0.2) +
      geom_density(data=df_posterior, aes(x=GMR, y=..density..), fill="red", alpha=0.4) +
      geom_vline(xintercept=ci, linetype="dashed", color="black") +
      annotate("text", x=ci, y=0, label=round(ci, 2), vjust=-0.5, hjust=c(1.1,-0.1)) +
      labs(title="Prior (blue) and Posterior (red) Distributions",
           x="GMR", y="Density") +
      scale_x_continuous(limits = c(.5,1.5), breaks = c(seq(0.5,1.5,length.out = 5),1.11,0.9)) +
      theme_minimal()
  })
}

shinyApp(ui, server)
