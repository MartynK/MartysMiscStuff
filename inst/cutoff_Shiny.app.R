library(shiny)
library(ggplot2)

ui <- fluidPage(
  titlePanel("Acceptance/Rejection Regions"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("a", "Coefficient a", min = -5, max = 5, value = 1, step = 0.1),
      sliderInput("b", "Coefficient b", min = -5, max = 5, value = 1, step = 0.1),
      sliderInput("c", "Exponent c", min = 0.1, max = 5, value = 1.1, step = 0.02),
      sliderInput("cutoff", "Cutoff", min = 0, max = 10, value = 0, step = 0.1)
    ),
    mainPanel(
      plotOutput("plot")
    )
  )
)

server <- function(input, output) {
  output$plot <- renderPlot({
    set.seed(123)  # For reproducibility
    n <- 20000
    x <- runif(n, -2, 2)
    y <- runif(n, -2, 2)

    df <- data.frame(x = x, y = y)

    # Calculate the acceptance condition: a * x^c + b * y^c > cutoff
    df$val <- input$a * (as.complex(df$x) ^ input$c) +
                input$b * (as.complex(df$y) ^ input$c)

    # Color points based on the condition: red if accepted, blue if rejected
    df$color <- ifelse(
      abs(df$val ^ input$c)  > input$cutoff, "red", "blue")

    ggplot(df, aes(x = x, y = y)) +
      geom_point(aes(color = color), alpha = .3, size = 1.5) +
      scale_color_identity() +
      labs(
        title = "Acceptance/Rejection Regions",
        subtitle = "Red: Accepted (val > cutoff) | Blue: Rejected",
        x = "x",
        y = "y"
      ) +
      theme_minimal()
  })
}

shinyApp(ui = ui, server = server)
