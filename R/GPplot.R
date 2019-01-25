library(shiny)

ui <- fluidPage(
  fluidRow(
    column(width = 6,
      sliderInput("noise", "noise", min = 0, max = 1, value = 0.1),
      sliderInput("gennoise", "generating noise", min = 0, max = 1, value = 0.1),
      numericInput("n", "number of training points", min = 1, value = 10),
      checkboxInput("drawrand", "Choose training points at random"),
      sliderInput("xlim", "xlim", min = -10, max = 10, value = c(-5,5))
      ),
    
      
    column(width = 6,
      selectInput("cov", "Covariance Function", choices = list("Squared Exponential", "Gamma Exponential", 
                                                                    "Constant", "Linear","Polynomial", 
                                                                    "Rational Quadratic")),
      uiOutput("selectors")
      )
  ),
  fluidRow(
    column(width = 12,  
      plotOutput("plot1")))
)

#Helper function to remove sliders
switchUI <- function(i){
  n <- 6 #number of options
  for (j in (1:n)[-i]) removeUI(selector = paste0("div#selectdiv",j))
}

server <- function(input, output){
  observeEvent(input$cov, {
    if (input$cov == "Squared Exponential"){
      switchUI(1)
      output$selectors <- renderUI({
        tags$div(id = "selctdiv1",
                 withMathJax("$$\\huge{\\sigma_1 \\cdot \\text{exp} \\left( \\frac{|x_p-x_q|^2}{2 \\ell^2} \\right) + \\sigma_2 \\delta_{pq}}$$"),
                 sliderInput("sigma_1", withMathJax("$$\\huge{\\sigma_1}$$"), min = 0.01, max = 3, value = 1),
                 sliderInput("sigma_2", withMathJax("$$\\huge{\\sigma_2}$$"), min = 0.00001, max = 1, value = 0.00001),
                 sliderInput("l", withMathJax("$$\\huge{\\ell}$$"), min = 0.1, max = 3, value = 1)
                  
        )
      })
    }
    if (input$cov == "Constant"){
      switchUI(2)
      output$selectors <- renderUI({
        tags$div(id = "selctdiv2", 
                 withMathJax("c"),
                 sliderInput("c", "c", min = 0, max = 10, value = 1))
      })
    }
    if (input$cov == "Linear"){
      switchUI(3)
      output$selectors <- renderUI({
        tags$div(id = "selctdiv3", 
                 withMathJax("$$\\huge{c \\cdot \\sum_{d=1}^D x_d x_d^\\top}$$"),
                 sliderInput("c", "c", min = 0, max = 10, value = 1))
      })
    }
    if (input$cov == "Polynomial"){
      switchUI(4)
      output$selectors <- renderUI({
        tags$div(id = "selctdiv4", 
                 withMathJax("$$\\huge{(x \\cdot x'^\\top + \\sigma)^p}$$"),
                 sliderInput("sigma_3", withMathJax("$$\\huge{\\sigma}$$"), min = 0, max = 10, value = 1),
                 sliderInput("p", "p", min = 1, max = 10, value = 1))
      })
    }
    if (input$cov == "Gamma Exponential"){
      switchUI(5)
      output$selectors <- renderUI({
        tags$div(id = "selctdiv5", 
                 withMathJax("$$\\huge{\\text{exp} (- \\frac{|x-x'|}{\\ell})^\\gamma}$$"),
                 sliderInput("sigma_4", withMathJax("$$\\huge{\\gamma}$$"), min = 0, max = 10, value = 1),
                 sliderInput("l2", withMathJax("$$\\huge{\\ell}$$"), min = 1, max = 10, value = 1))
      })
    }
    if (input$cov == "Rational Quadratic"){
      switchUI(6)
      output$selectors <- renderUI({
        tags$div(id = "selctdiv6", 
                 withMathJax("$$\\huge{\\left( 1 + \\frac{|x-x'|^2}{2 \\alpha \\ell^2}\\right)^{-\\alpha}}$$"),
                 sliderInput("alpha", withMathJax("$$\\huge{\\alpha}$$"), min = 0, max = 10, value = 1),
                 sliderInput("l3", withMathJax("$$\\huge{\\ell}$$"), min = 1, max = 10, value = 1))
      })
    }
  })
  output$plot1 <- renderPlot({
    if (input$drawrand){
      X <- reactive(matrix(runif(input$n,input$xlim[1],input$xlim[2]), nrow = 1))
    }
    else{
      X <- reactive(matrix(seq(input$xlim[1],input$xlim[2],
                               by = (input$xlim[2] - input$xlim[1])/input$n), nrow = 1))
    }
    y <- reactive(c(sin(X()) + rnorm(length(X()),0,sqrt(input$gennoise))))
    Gaussian <- reactive(GPR.sqrexp$new(X(), y(), l = 1, noise = input$noise))
    if (input$cov == "Squared Exponential"){
      if ((!is.null(input$sigma_1)) & (!is.null(input$sigma_2)) & (!is.null(input$l))){
        kappa <- reactive(function(x,y){
        input$sigma_1 * exp(-(1/(2*input$l^2))*(x - y)^2) + input$sigma_2 * ifelse(isTRUE(all.equal(x,y)),1,0)
        })
        Gaussian <- reactive(GPR$new(X(), y(), kappa(), input$noise))
      }
    }
    if (input$cov == "Constant"){
      if (!is.null(input$c)) Gaussian <- reactive(GPR.constant$new(X(), y(), input$c, input$noise))
    }
    if (input$cov == "Linear"){
      if (!is.null(input$c)) Gaussian <- reactive(GPR.linear$new(X(), y(), input$c, input$noise))
    }
    if (input$cov == "Polynomial"){
      if (!is.null(input$sigma_3) & !is.null(input$p)) Gaussian <- reactive(GPR.polynomial$new(X(), y(), input$sigma_3, input$p, input$noise))
    }
    if (input$cov == "Gamma Exponential"){
      if (!is.null(input$sigma_4) & !is.null(input$l2)){
      Gaussian <- reactive(GPR.gammaexp$new(X(), y(), input$sigma_4, input$l2, input$noise))
      }
    }
    if (input$cov == "Rational Quadratic"){
      if (!is.null(input$alpha) & !is.null(input$l3)){
      Gaussian <- reactive(GPR.rationalquadratic$new(X(), y(), input$alpha, input$l3, input$noise))
    }
    }
    Gaussian()$plot(seq(input$xlim[1] - 1, input$xlim[2] + 1, by = 0.1))
  })
}

shinyApp(ui, server)

