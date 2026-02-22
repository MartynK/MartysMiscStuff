########################  ui.R  ########################
# -----------------------------------------------------
# Shiny UI for pooling ISCVs from 2×2 crossover trials
# -----------------------------------------------------
library(shiny)

shinyUI(
  fluidPage(
    titlePanel("Bayesian pooling of within-subject CV (2×2 crossover)"),
    sidebarLayout(
      sidebarPanel(
        # -------------------------------------------------
        # 1. How many studies? (build a grid of inputs)
        # -------------------------------------------------
        numericInput(
          "k", "Number of studies (1–20):",
          value = 2, min = 1, max = 20, step = 1
        ),
        uiOutput("dynInputs"),   # << rendered in server

        # -------------------------------------------------
        # 2. MCMC settings
        # -------------------------------------------------
        numericInput(
          "ndraw",
          "Posterior draws per chain (>= 1000):",
          value = 10000, min = 1000, step = 1000
        ),

        actionButton("update", "Update analysis", class = "btn-primary")
      ),

      mainPanel(
        tabsetPanel(
          tabPanel("Prior & posterior",
                   plotOutput("priorPlot", height = 280),
                   plotOutput("postPlot",  height = 280)),
          tabPanel("Posterior summary",
                   DT::DTOutput("postSum"))
        )
      )
    )
  )
)
