library(shiny)

ui <- fluidPage(
  selectInput("cov", "Covariance Function", choices = list("Gauss", "Exp", "Unif")),
  sliderInput("noise", "noise", min = 0, max = 10, value = 0.1),
  checkboxInput("drawrand", "Draw training points at random"),
  sliderInput("xlim", "xlim", min = -10, max = 10, value = c(-5,5)),
  plotOutput("plot1")
)
server <- function(input, output){
  output$plot1 <- renderPlot({
    if (input$drawrand){
      X <- reactive(matrix(runif(10,input$xlim[1],input$xlim[2]), nrow = 1))
    }
    else{
      X <- reactive(matrix(seq(input$xlim[1],input$xlim[2],by = (input$xlim[2] - input$xlim[1])/10), nrow = 1))
    }
    y <- reactive(c(X()^2 + rnorm(length(X()),0,sqrt(input$noise))))
    #hyperparameters
    sigma_1 <- 1
    sigma_2 <- 0.1
    l <- 1
    kappa <- function(x,y) sigma_1 * exp(-(1/(2*l^2))*(x - y)^2) + sigma_2 * ifelse(isTRUE(all.equal(x,y)),1,0)
    Gaussian <- reactive(GPR$new(X(), y(), kappa, input$noise))
    Gaussian()$plot(seq(input$xlim[1],input$xlim[2], by = 0.5))
  })
}

shinyApp(ui, server)
