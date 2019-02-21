#'@export
GPplot2 <- function(){
  library(shiny)
  
  `%then%` <- shiny:::`%OR%`
  
  ui <- fluidPage(tabsetPanel(tabPanel("Regression",
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
                                           div(style="display: inline-block;vertical-align:top; width: 300px;", 
                                               checkboxInput("drawrand", "Choose training points at random")),
                                           div(style="display: inline-block;vertical-align:top; width: 50px;", 
                                               actionButton("refrng", "", icon = icon("sync"))),
                                           checkboxInput("drawtrue", "Draw true function"),
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
                                         ), tabPanel("Classification",
                                                     fluidRow(column(width = 6, wellPanel(
                                                       h4("Datapoints"),
                                                       numericInput("n2", "Number of training points", min = 1, value = 10)
                                                     )),
                                                     column(width = 6, wellPanel(
                                                       h4("Datapoints"),
                                                       textOutput("info"),
                                                       radioButtons("label", "Label of added points", choices = c(-1, 1), inline = T),
                                                       actionButton("refresh", "Refresh data boundaries", icon = icon("sync"))
                                                     ))),
                                                     fluidRow(column(width = 12, plotOutput("plot2", click = "plot_click")))
                                         )
                                       ))
  
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
      tags$div(id = paste0("selctdiv", i), kdesc, ...)
    })
  }
  
  server <- function(input, output, session){
    #refresh random generator with button
    seed <- eventReactive(input$refrng, {
      round(as.integer(runif(1, 0, 100)))
    }, ignoreNULL = FALSE)
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
                                                  "exp(-sum((x - y)^2)/(2 * l^2))",
                                                  sliderInput("par1", "l", min = 0.01, max = 10, value = 1))
             },
             "Constant" = {
               output$selectors <- switchrenderUI(2, session, 0.1, input$noise, "c",
                                                  sliderInput("par1", "c", min = 0.1, max = 10, value = 1))
             },
             "Linear" = {
               output$selectors <- switchrenderUI(3, session, 0.1, input$noise,
                                                  "c * sum(x * y)",
                                                  sliderInput("par1", "c", min = 0, max = 10, value = 1.2, step = 0.1))
             },
             "Polynomial" = {
               output$selectors <- switchrenderUI(4, session, 0.1, input$noise,
                                                  "(x %*% y + sigma)^p",
                                                  sliderInput("par1", "sigma", min = 0, max = 10, value = 1, step = 0.01),
                                                  sliderInput("par2", "p", min = 1, max = 10, value = 5))
             },
             "Gamma Exponential" = {
               output$selectors <- switchrenderUI(5, session, 0, input$noise,
                                                  "exp(-(sqrt(sum((x - y)^2)) / l) ^ gamma)",
                                                  sliderInput("par1", "gamma", min = 1, max = 10, value = 2),
                                                  sliderInput("par2", "l", min = 0.01, max = 3, value = 1))
             },
             "Rational Quadratic" = {
               output$selectors <- switchrenderUI(6, session, 0, input$noise,
                                                  "(1 + sum((x - y)^2) / (2 * alpha * l^2))^(-alpha)",
                                                  sliderInput("par1", "alpha", min = 0.1, max = 3, value = 1, step = 0.1),
                                                  sliderInput("par2", "l", min = 0.01, max = 7, value = 1))
             }
      )
    })
    X <- reactive({validate(need(input$n, "Number of datapoints can't be empty"))
      if (input$drawrand){
        set.seed(seed()) #keep the randomly generated datapoints fixed so changing other parameter doesn't change them
        matrix(runif(input$n, input$xlim[1], input$xlim[2]), nrow = 1)
      }
      else{
        matrix(seq(input$xlim[1],input$xlim[2], by = (input$xlim[2] - input$xlim[1])/input$n), nrow = 1)
      }})
    #generate datapoints using function f
    y <- reactive({set.seed(seed())
      c(f()(X()) + rnorm(length(X()), 0, sqrt(input$gennoise)))})
    observeEvent(input$opthyp,{
      z <- fit(X(), y(), input$noise + 0.1, list(cov_df$name[cov_df$display == input$cov]))
      for (i in seq_along(z$par)) updateSliderInput(session, sprintf("par%s", i), value = z$par[i])
    })
    output$plot1 <- renderPlot({
      #standard plot if nothing is selected:
      Gaussian <- reactive(GPR.sqrexp$new(X(), y(), noise = input$noise, l = 1)) 
      switch(input$cov, 
             "Squared Exponential" = {validate(need(input$par1, "Invalid parameters"))
               Gaussian <- reactive(GPR.sqrexp$new(X(), y(), input$noise, input$par1))
             },
             "Constant" = {validate(need(input$par1, "Invalid parameters") %then% 
                                      need(input$par1 > 0, "Invalid parameters") )
               Gaussian <- reactive(GPR.constant$new(X(), y(), input$noise, input$par1))
             },
             "Linear" = {validate(need(input$par1, "Invalid parameters") %then% 
                                    need(input$noise > 0, "Invalid parameters"))
               Gaussian <- reactive(GPR.linear$new(X(), y(), input$noise, input$par1))
             },
             "Polynomial" = {validate(need(input$par1, "Invalid parameters") %then% 
                                        need(input$noise > 0, "Invalid parameters") %then% 
                                        need(input$par2, "Invalid parameters") )
               Gaussian <- reactive(GPR.polynomial$new(X(), y(), input$noise, input$par1, input$par2))
             },
             "Gamma Exponential" = {validate(need(input$par1, "Invalid parameters") %then% 
                                               need(input$par2, "Invalid parameters"))
               Gaussian <- reactive(GPR.gammaexp$new(X(), y(), input$noise, input$par1, input$par2))
             },
             "Rational Quadratic" = {validate(need(input$par1, "Invalid parameters") %then% 
                                                need(input$par2, "Invalid parameters"))
               Gaussian <- reactive(GPR.rationalquadratic$new(X(), y(), input$noise, input$par1, input$par2))
             })
      X_points <- reactive(seq(input$xlim[1], input$xlim[2], by = 0.1))
      updateSliderInput(session, "noise", value = Gaussian()$noise)
      p <- Gaussian()$plot()$plot
      if (input$drawtrue){
        p <- p + ggplot2::geom_line(data = data.frame(x = X_points(), y = sapply(X_points(), f())), 
                                    ggplot2::aes(x = x, y = y), linetype = "dashed")
      }
      p
    })
    #Classification Panel
    click_saved <- reactiveValues(singleclick = NULL)
    rv = reactiveValues(m=data.frame(x = 0, y = 0, label = -1))
    observeEvent(input$plot_click,{
      click_saved$singleclick <- input$plot_click
      rv$m <- rbind(rv$m, c(input$plot_click$x, input$plot_click$y, as.double(input$label)))
    })
    
    X_c <- eventReactive(input$refresh, {
      set.seed(0)
      cbind(cbind(multivariate_normal(input$n2, c(0.5, 0.5), diag(c(0.1, 0.1))), 
                  multivariate_normal(input$n2, c(-0.5, -0.5), diag(c(0.1, 0.1)))), t(as.matrix(rv$m[, 1:2])))
    }, ignoreNULL = FALSE)
    y_c <- eventReactive(input$refresh, {
      c(rep(c(1, -1), each = input$n2), rv$m$label)
    }, ignoreNULL = FALSE)
    output$info <- renderText({
      paste0(input$plot_click$x, input$plot_click$y)
    })
    dat <- eventReactive(c(input$refresh,input$n2),{
      kappa <- function(x, y) sqrexp(x, y, l = 1)
      gaussian_classifier <- GPC$new(X_c(), y_c(), 1e-5, kappa)
      s <- seq(min(X_c()), max(X_c()), by = 0.1)
      testpoints <- matrix(c(rep(s, each = length(s)), rep(s, times = length(s))), nrow = 2, byrow = T)
      predictions <- gaussian_classifier$predict_class(testpoints)
      data.frame(x.1 = testpoints[1,], x.2 = testpoints[2,], y_pred = 2*as.integer(predictions >= 0.5) - 1)
    })
    output$plot2 <- renderPlot({
      ggplot2::ggplot(dat(), inherit.aes = F, ggplot2::aes(x = x.1, y = x.2, fill = factor(y_pred))) +
        ggplot2::theme_classic() +
        ggplot2::scale_y_continuous(expression("x_2")) +
        ggplot2::scale_x_continuous(expression("x_1")) +
        ggplot2::geom_tile() + 
        ggplot2::scale_fill_manual(values = c("red", "blue")) +
        ggplot2::guides(fill = ggplot2::guide_legend(title = "Labels")) +
        ggplot2::geom_point(inherit.aes = F, data = data.frame(xpoints = c(X_c()[1,]), ypoints = c(X_c()[2,]), ylab = y_c()), 
                            mapping = ggplot2::aes(x = xpoints, y = ypoints, shape = factor(ylab))) +
        ggplot2::scale_shape_manual(values = c(4, 2)) +
        ggplot2::guides(shape = ggplot2::guide_legend(title = "Testpoints")) +
        ggplot2::scale_color_manual(values = c("red", "blue")) + 
        ggplot2::geom_point(inherit.aes = F, data = rv$m[-1,], 
                            mapping = ggplot2::aes(x = x, y = y, shape = factor(label)))
    }, width = 600, height = 512)
  }
  
  shinyApp(ui, server)
}
