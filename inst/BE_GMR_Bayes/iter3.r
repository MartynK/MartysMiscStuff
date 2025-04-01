library(shiny)
library(rstan)
library(ggplot2)

ui <- fluidPage(
  titlePanel("Bioequivalence GMR Bayesian Dashboard"),

  sidebarLayout(
    sidebarPanel(
      numericInput("sample_size", "Sample Size:", value = 24, min = 3, max = 90),
      numericInput("true_gmr", "True GMR:", value = 1.05, min = 0.7, max = 1/0.7, step=0.01),
      numericInput("cv", "CV%:", value = 36, min = 1, max = 100, step=1),
      selectInput("prior_choice", "Prior Distribution:",
                  choices = c("Normal (centered normal)" = "extreme",
                              "Uniform (0.7, 1.3)" = "uniform",
                              "Wide Beta (0,2)" = "wide_beta")),
      actionButton("execute", "Execute")
    ),

    mainPanel(
      plotOutput("prior_posterior_plot"),
      plotOutput("ci_distribution_plot")
    )
  )
)

server <- function(input, output) {
  stan_model_reactive <- reactiveVal(NULL)

  observeEvent(input$prior_choice, {
    stan_model_code <- switch(input$prior_choice,
                              "uniform" = "
data {
  int<lower=1> N;
  vector[N] y;
}
parameters {
  real<lower=log(0.7), upper=log(1.3)> mu;
  real<lower=0> sigma;
}
model {
  mu ~ uniform(log(0.7), log(1.3));
  sigma ~ cauchy(0, 2.5);
  y ~ normal(mu, sigma);
}
generated quantities {
  real GMR = exp(mu);
}",
"extreme" = "
data {
  int<lower=1> N;
  vector[N] y;
}
parameters {
  real<lower=log(0.7), upper=log(1.3)> mu;
  real<lower=0> sigma;
}
model {
  mu ~ normal(log(1), 0.5);
  sigma ~ cauchy(0, 2.5);
  y ~ normal(mu, sigma);
}
generated quantities {
  real GMR = exp(mu);
}",
"wide_beta" = "
data {
  int<lower=1> N;
  vector[N] y;
}
parameters {
  real<lower=0, upper=2> gmr_raw;
  real<lower=0> sigma;
}
transformed parameters {
  real mu = log(gmr_raw);
}
model {
  gmr_raw ~ beta(5, 5);
  sigma ~ cauchy(0, 2.5);
  y ~ normal(mu, sigma);
}
generated quantities {
  real GMR = gmr_raw;
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

    # Iterating over seeds
    ci_lower <- numeric(100)
    ci_upper <- numeric(100)
    withProgress(message = 'Running iterations', value = 0, {
      for (i in 1:100) {
        set.seed(i)
        sim_data_iter <- rnorm(N, mean = true_mu, sd = true_sigma)
        fit_iter <- sampling(stan_model_reactive(),
                             data = list(N = N, y = sim_data_iter),
                             iter = 1000, chains = 2, refresh=0)
        posterior_iter <- extract(fit_iter)$GMR
        ci_iter <- quantile(posterior_iter, c(0.025, 0.975))
        ci_lower[i] <- ci_iter[1]
        ci_upper[i] <- ci_iter[2]
        incProgress(1/100)
      }
    })

    list(posterior_samples = posterior_samples, ci = ci, prior = input$prior_choice,
         ci_lower = ci_lower, ci_upper = ci_upper)
  })

  output$prior_posterior_plot <- renderPlot({
    req(results())

    posterior_samples <- results()$posterior_samples
    ci <- results()$ci

    df_posterior <- data.frame(GMR=posterior_samples)

    ggplot(df_posterior, aes(x=GMR)) +
      geom_density(fill="red", alpha=0.4) +
      geom_vline(xintercept=ci, linetype="dashed", color="black") +
      annotate("text", x=ci, y=0, label=round(ci, 2), vjust=-0.5, hjust=c(1.1,-0.1)) +
      labs(title="Posterior Distribution with Credible Interval",
           x="GMR", y="Density") +
      theme_minimal()
  })

  output$ci_distribution_plot <- renderPlot({
    req(results())

    df_ci <- data.frame(lower=results()$ci_lower, upper=results()$ci_upper)

    ggplot(df_ci) +
      geom_histogram(aes(x=lower), fill="blue", alpha=0.5, bins=30) +
      geom_histogram(aes(x=upper), fill="green", alpha=0.5, bins=30) +
      labs(title="Distribution of Credible Interval Limits over 100 seeds",
           x="GMR", y="Frequency") +
      theme_minimal()
  })
}

shinyApp(ui, server)
