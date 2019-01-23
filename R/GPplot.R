library(shiny)

ui <- fluidPage(
  selectInput("cov", "Covariance Function", choices = list("Gauss", "Exp", "Unif")),
  sliderInput("noise", "noise", min = -0.1, max = 10, value = 0.1),
  sliderInput("sigma_1", "sigma_1", min = 0.01, max = 3, value = 0.1),
  sliderInput("sigma_2", "sigma_2", min = 0.00001, max = 1, value = 0.1),
  sliderInput("l", "l", min = 0.1, max = 3, value = 0.1),
  checkboxInput("drawrand", "Draw training points at random"),
  sliderInput("xlim", "xlim", min = -10, max = 10, value = c(-5,5)),
  plotOutput("plot1")
)
server <- function(input, output){
  output$plot1 <- renderPlot({
    n <- 5 #number of training points in the intervall
    if (input$drawrand){
      X <- reactive(matrix(runif(n,input$xlim[1],input$xlim[2]), nrow = 1))
    }
    else{
      X <- reactive(matrix(seq(input$xlim[1],input$xlim[2],by = (input$xlim[2] - input$xlim[1])/n), nrow = 1))
    }
    y <- reactive(c(sin(X()) + rnorm(length(X()),0,sqrt(input$noise))))
    #hyperparameters
    sigma_1 <- 1
    sigma_2 <- 0.1
    l <- 1
    kappa <- function(x,y){
      input$sigma_1 * exp(-(1/(2*input$l^2))*(x - y)^2) + input$sigma_2 * ifelse(isTRUE(all.equal(x,y)),1,0)
    }
    Gaussian <- reactive(GPR$new(X(), y(), kappa, input$noise))
    Gaussian()$plot(seq(input$xlim[1],input$xlim[2], by = 0.1))
  })
}

shinyApp(ui, server)
