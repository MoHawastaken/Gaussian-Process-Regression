library(shiny)

seed <- runif(1,0,100)
f <- sin #function used to create datapoints
cov_names <- list("Squared Exponential" = "sqrexp", "Gamma Exponential" = "gammaexp", 
                  "Constant" = "constant", "Linear" = "linear", "Polynomial" = "polynomial", 
                  "Rational Quadratic" = "rationalquadratic")

ui <- fluidPage(
  fluidRow(
    column(width = 6,
      sliderInput("noise", "noise", min = 0, max = 1, value = 0),
      sliderInput("gennoise", "generating noise", min = 0, max = 1, value = 0),
      numericInput("n", "number of training points", min = 1, value = 10),
      checkboxInput("drawrand", "Choose training points at random"),
      sliderInput("xlim", "xlim", min = -10, max = 10, value = c(-5,5))
      ),
      
    column(width = 6,
      selectInput("cov", "Covariance Function", choices = names(cov_names)),
      uiOutput("selectors"),
      actionButton("opthyp","Optimize Hyperparameters")
      )
  ),
  fluidRow(
    column(width = 12,  
      plotOutput("plot1")))
)

#Helper function to remove sliders
switchrenderUI <- function(i, session, min_noise, kdesc, ...){
  n <- 6 #number of options
  for (j in (1:n)[-i]) removeUI(selector = paste0("div#selectdiv", j)) #remove other UIs
  #change minimum value for noise
  updateSliderInput(session, "noise", value = min_noise,
                    min = min_noise)
  #return new UI
  renderUI({
    tags$div(id = paste0("selctdiv",i),withMathJax(paste("$$\\huge{",kdesc,"}$$")),...)
  })
}

server <- function(input, output,session){
  observeEvent(input$cov, {
    #switch formula and sliders for parameters for each covariance function
    switch(input$cov, 
    "Squared Exponential" = {
      output$selectors <- switchrenderUI(1,session,0,
        "\\sigma \\cdot \\text{exp} \\left( \\frac{|x_p-x_q|^2}{2 \\ell^2} \\right)",
        sliderInput("par1", withMathJax("$$\\huge{\\ell}$$"), min = 0.01, max = 3, value = 1),
        sliderInput("par2", withMathJax("$$\\huge{\\sigma}$$"), min = 0.1, max = 3, value = 1))
    },
    "Constant" = {
      output$selectors <- switchrenderUI(2,session,0.1,"c",
                     sliderInput("par1", "c", min = 0.1, max = 10, value = 1))
    },
    "Linear" = {
      output$selectors <- switchrenderUI(3,session,0.1,
                 "c \\cdot \\sum_{d=1}^D x_d x_d^\\top",
                 sliderInput("par1", "c", min = 0, max = 10, value = 1))
    },
    "Polynomial" = {
      output$selectors <- switchrenderUI(4,session,0.1, 
                 "(x \\cdot x'^\\top + \\sigma)^p",
                 sliderInput("par1", withMathJax("$$\\huge{\\sigma}$$"), min = 0, max = 10, value = 1),
                 sliderInput("par2", withMathJax("$$\\huge{p}$$"), min = 1, max = 10, value = 5))
    },
    "Gamma Exponential" = {
      output$selectors <- switchrenderUI(5,session,0, 
                 "\\text{exp} \\left( - \\left(\\frac{|x-x'|}{\\ell}\\right)^\\gamma \\right)",
                 sliderInput("par1", withMathJax("$$\\huge{\\gamma}$$"), min = 0, max = 10, value = 1),
                 sliderInput("par2", withMathJax("$$\\huge{\\ell}$$"), min = 0.01, max = 3, value = 1))
    },
    "Rational Quadratic" = {
      output$selectors <- switchrenderUI(6,session,0, 
                 "\\left( 1 + \\frac{|x-x'|^2}{2 \\alpha \\ell^2}\\right)^{-\\alpha}",
                 sliderInput("par1", withMathJax("$$\\huge{\\alpha}$$"), min = 0.5, max = 10, value = 1),
                 sliderInput("par2", withMathJax("$$\\huge{\\ell}$$"), min = 0.01, max = 3, value = 1))
    }
    )
  })
  observeEvent(input$opthyp,{
    print(cov_names[[input$cov]])
    z <- fit(X,y,input$noise+0.1,list(cov_names[[input$cov]]))
    for (i in seq_along(z$par)) updateSliderInput(session, sprintf("par%s",i), value = z$par[i])
    print(z$par)
  })
  output$plot1 <- renderPlot({
    if (input$drawrand){
      set.seed(seed) #keep the randomly generated datapoints fixed so changing other parameter doesn't change them
      X <- reactive(matrix(runif(input$n,input$xlim[1],input$xlim[2]), nrow = 1))
    }
    else{
      X <- reactive(matrix(seq(input$xlim[1],input$xlim[2],
                               by = (input$xlim[2] - input$xlim[1])/input$n), nrow = 1))
    }
    #generate datapoints using function f
    y <- reactive(c(f(X()) + rnorm(length(X()),0,sqrt(input$gennoise))))
    #standard plot if nothing is selected:
    Gaussian <- reactive(GPR.sqrexp$new(X(), y(), l = 1, noise = input$noise)) 
    switch(input$cov, 
      "Squared Exponential" = {
        if (!is.null(input$par1) & !is.null(input$par2)){
          kappa <- reactive(function(x,y){
          input$par2 * exp(-(1/(2*input$par1^2))*(x - y)^2)
        })
        Gaussian <- reactive(GPR$new(X(), y(), kappa(), input$noise))
        }
      },
     "Constant" = {
        if (!is.null(input$par1)) Gaussian <- reactive(GPR.constant$new(X(), y(), input$par1, input$noise))
      },
      "Linear" = {
        if (!is.null(input$par1)) Gaussian <- reactive(GPR.linear$new(X(), y(), input$par1, input$noise))
      },
      "Polynomial" = {
        if (!is.null(input$par1) & !is.null(input$par2)){
        Gaussian <- reactive(GPR.polynomial$new(X(), y(), input$par1, input$par2, input$noise))}
      },
      "Gamma Exponential" = {
        if (!is.null(input$par1) & !is.null(input$par2)){
        Gaussian <- reactive(GPR.gammaexp$new(X(), y(), input$par1, input$par2, input$noise))}
      },
      "Rational Quadratic" = {
      if (!is.null(input$par1) & !is.null(input$par2)){
        Gaussian <- reactive(GPR.rationalquadratic$new(X(), y(), input$par1, input$par2, input$noise))}
      })
    Gaussian()$plot(seq(input$xlim[1], input$xlim[2], by = 0.1))
  })
}

shinyApp(ui, server)


