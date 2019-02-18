library(shiny)

`%then%` <- shiny:::`%OR%`
seed <- runif(1, 0, 100)

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
                    .shiny-output-error-validation {
                    color: red;
                    }
                    "))
    ),
  fluidRow(
    column(width = 6, wellPanel(
      h4("Datapoints"),
      textInput("func", "Data function of x", value = "sin(x)"),
      sliderInput("noise", "Noise", min = 0, max = 1, value = 0),
      sliderInput("gennoise", "Generating noise", min = 0, max = 1, value = 0),
      numericInput("n", "Number of training points", min = 1, value = 10),
      checkboxInput("drawrand", "Choose training points at random"),
      sliderInput("xlim", "Data Boundaries", min = -10, max = 10, value = c(-5, 5))
      )),
      
    column(width = 6, wellPanel(
      h4("Covariance function"),
      selectInput("cov", "", choices = cov_df$display),
      uiOutput("selectors"),
      actionButton("opthyp", "Optimize Hyperparameters")
      ))
  ),
  fluidRow(
    column(width = 12,  
      plotOutput("plot1")))
)

#Helper function to remove sliders
switchrenderUI <- function(i, session, min_noise, act_noise, kdesc, ...){
  n <- nrow(cov_df) #number of options
  for (j in (1:n)[-i]) removeUI(selector = paste0("div#selectdiv", j)) #remove other UIs
  #change minimum value for noise
  updateSliderInput(session, "noise", min = min_noise)
  if (act_noise == 0) updateSliderInput(session, "noise", value = min_noise,
                                        min = min_noise)
  #return new UI
  renderUI({
    tags$div(id = paste0("selctdiv", i), withMathJax(paste("$$\\huge{", kdesc, "}$$")), ...)
  })
}

server <- function(input, output, session){
  #function used to create datapoints
  f <- reactive({
    validate(
      need(try(eval({x <- 0; parse(text = input$func)})), "Input a valid function definition") %then%
      need(input$func != "", "Input a valid function definition"))
    function(x) eval(parse(text = input$func))
  })
  observeEvent(input$cov, {
    #switch formula and sliders for parameters for each covariance function
    switch(input$cov, 
    "Squared Exponential" = {
      output$selectors <- switchrenderUI(1, session, 0, input$noise,
        "\\text{exp} \\left(- \\frac{||x-x'||^2}{2 \\ell^2} \\right)",
        sliderInput("par1", withMathJax("$$\\huge{\\ell}$$"), min = 0.01, max = 10, value = 1))
    },
    "Constant" = {
      output$selectors <- switchrenderUI(2, session, 0.1, input$noise, "c",
                     sliderInput("par1", "c", min = 0.1, max = 10, value = 1))
    },
    "Linear" = {
      output$selectors <- switchrenderUI(3, session, 0.1, input$noise,
                 "c \\cdot \\sum_{d=1}^D x_d x_d'",
                 sliderInput("par1", "c", min = 0, max = 10, value = 1.2, step = 0.1))
    },
    "Polynomial" = {
      output$selectors <- switchrenderUI(4, session, 0.1, input$noise,
                 "(x \\cdot x'^\\top + \\sigma)^p",
                 sliderInput("par1", withMathJax("$$\\huge{\\sigma}$$"), min = 0, max = 10, value = 1, step = 0.01),
                 sliderInput("par2", withMathJax("$$\\huge{p}$$"), min = 1, max = 10, value = 5))
    },
    "Gamma Exponential" = {
      output$selectors <- switchrenderUI(5, session, 0, input$noise,
                 "\\exp \\left( - \\left(\\frac{||x-x'||}{\\ell}\\right)^\\gamma \\right)",
                 sliderInput("par1", withMathJax("$$\\huge{\\gamma}$$"), min = 1, max = 10, value = 2),
                 sliderInput("par2", withMathJax("$$\\huge{\\ell}$$"), min = 0.01, max = 3, value = 1))
    },
    "Rational Quadratic" = {
      output$selectors <- switchrenderUI(6, session, 0, input$noise,
                 "\\left( 1 + \\frac{||x-x'||^2}{2 \\alpha \\ell^2}\\right)^{-\\alpha}",
                 sliderInput("par1", withMathJax("$$\\huge{\\alpha}$$"), min = 0.1, max = 3, value = 1, step = 0.1),
                 sliderInput("par2", withMathJax("$$\\huge{\\ell}$$"), min = 0.01, max = 7, value = 1))
    }
    )
  })
  X <- reactive({validate(need(input$n, "Number of datapoints can't be empty"))
    if (input$drawrand){
    set.seed(seed) #keep the randomly generated datapoints fixed so changing other parameter doesn't change them
    matrix(runif(input$n, input$xlim[1], input$xlim[2]), nrow = 1)
  }
  else{
    matrix(seq(input$xlim[1],input$xlim[2], by = (input$xlim[2] - input$xlim[1])/input$n), nrow = 1)
  }})
  #generate datapoints using function f
  set.seed(seed)
  y <- reactive(c(f()(X()) + rnorm(length(X()), 0, sqrt(input$gennoise))))
  observeEvent(input$opthyp,{
    z <- fit(X(), y(), input$noise + 0.1, list(cov_df$name[cov_df$display == input$cov]))
    for (i in seq_along(z$par)) updateSliderInput(session, sprintf("par%s", i), value = z$par[i])
    print(paste("Optimal parameter: ", z$par))
  })
  output$plot1 <- renderPlot({
    #standard plot if nothing is selected:
    Gaussian <- reactive(GPR.sqrexp$new(X(), y(), l = 1, noise = input$noise)) 
    switch(input$cov, 
      "Squared Exponential" = {validate(need(input$par1, "Invalid parameters"))
        Gaussian <- reactive(GPR.sqrexp$new(X(), y(), input$par1, input$noise))
      },
     "Constant" = {validate(need(input$par1, "Invalid parameters") %then% 
                              need(input$par1 > 0, "Invalid parameters") )
        Gaussian <- reactive(GPR.constant$new(X(), y(), input$par1, input$noise))
      },
      "Linear" = {validate(need(input$par1, "Invalid parameters") %then% 
                             need(input$noise > 0, "Invalid parameters"))
        Gaussian <- reactive(GPR.linear$new(X(), y(), input$par1, input$noise))
      },
      "Polynomial" = {validate(need(input$par1, "Invalid parameters") %then% 
                                 need(input$noise > 0, "Invalid parameters") %then% 
                                 need(input$par2, "Invalid parameters") )
        Gaussian <- reactive(GPR.polynomial$new(X(), y(), input$par1, input$par2, input$noise))
      },
      "Gamma Exponential" = {validate(need(input$par1, "Invalid parameters") %then% 
                                        need(input$par2, "Invalid parameters"))
        Gaussian <- reactive(GPR.gammaexp$new(X(), y(), input$par1, input$par2, input$noise))
      },
      "Rational Quadratic" = {validate(need(input$par1, "Invalid parameters") %then% 
                                         need(input$par2, "Invalid parameters"))
        Gaussian <- reactive(GPR.rationalquadratic$new(X(), y(), input$par1, input$par2, input$noise))
      })
    Gaussian()$plot(seq(input$xlim[1], input$xlim[2], by = 0.1))
  })
}

shinyApp(ui, server)


