require(shiny)
require(deSolve)
require(phaseR)
library(shiny)
library(ggplot2)
library(gridExtra)
# Define UI
ui = pageWithSidebar(
    headerPanel("The Covid19-SEIRVT model"),
    sidebarPanel(
        #Sliders:
        sliderInput("beta_s","beta_s", 2.6695711,
                    min = .6,
                    max = 5
        ),
        sliderInput("beta_a","beta_a", 1.98927728,
                    min = .3,
                    max = 3
        ),
        sliderInput("epsilon", "epsilon", 0.12578,
                    min = 0.0001, max = 0.5
        ),
        sliderInput("delta_e", "delta_e", 0.125,
                    min = 0.0001, max = 0.5
        ),
        sliderInput("delta_v", "delta_v", 0.05,
                    min = 0.0001, max = 0.1
        ),
        sliderInput("p", "p", 0.701194082,
                    min = 0.0, max = 1.0
        ),
        sliderInput("q", "q", 0.701194082,
                    min = 0.0, max = 1.0
        ),
####alpha####
        sliderInput("alpha_a", "alpha_a", 0.16750418760,
                    min = 0.0, max = 1.0
        ),
        sliderInput("alpha_t", "alpha_t", 0.05,
            min = 0.0, max = 1.0
        ),
        sliderInput("alpha_s", "alpha_s", 2.99518149,
            min = 0.0, max = 4.0
        ),
####mu####
        sliderInput("mu_inverse", "1/mu", 70,
                    min = 40, max = 80
        ),
        sliderInput("mu_s", "mu_s", 0.11,
                    min = 0.05, max = 0.3
        ),
        sliderInput("mu_a", "mu_a", 0.0001,
                    min = 0.00001, max = 0.001
        ),
        sliderInput("lambda_v", "lambda_v", 0.00258,
                    min = 0.001 , max = 0.005
        ),
        sliderInput("lambda_t", "lambda_t", 0.00047895,
            min = 0.001 , max = 0.005
        ),
        sliderInput("n_whole", "n_whole", 905263,
                    min = 905000 , max = 910000
        ),
        sliderInput("T", "Time range:",
                    min = 0, max = 365, value = c(0,132))
    ),
    #Main panel for figures and equations
    mainPanel(
        #Multiple tabs in main panel
        tabsetPanel(
            #Tab 1: Time plot (plot1 from server)
            tabPanel("Time", plotOutput("plot1")),
            #Tab 2: Phase plot (plot2 from server)
            #Tab 3: MathJax typeset equations
            tabPanel("Equations",
                     withMathJax(
                         helpText("Susceptible
                            $$
                                \\frac{dS}{dt} =
                                \\mu (N - S) - \\frac{\\beta I S}{N}
                            $$
                         "),
                         helpText("Infecitous
                            $$
                                \\frac{dI}{dt} =
                                \\frac{\\beta I S}{N} -
                                (\\mu +  \\sigma) I
                            $$
                         "),
                         helpText("Removed $$\\frac{dR}{dt} =
                                    \\gamma I - \\mu R$$"),
                         helpText("Reproductive ratio $$R_0 =
                                  \\frac{1}{\\gamma+\\mu}
                                  \\frac{\\beta N}{N}$$")
                     )
            )
        )
    )
) #End of ui()

# Define server logic to plot various variables against mpg ----
server = function(input, output){
####Gradient function####
    seirvt_mod = function(t, state, parameters){
            with(as.list(c(state, parameters)), {
                ####
                n_bar <- s + e + i_s + i_a + r + v + treat
                force_infection = (beta_s * i_s + beta_a * i_a) / n_bar
                rhs_s = mu * n_bar - force_infection * s -
                    (mu + lambda_v) * s + delta_v * v
                rhs_e = force_infection * (epsilon * v + s) - (mu + delta_e) * e
                rhs_i_s = p * delta_e * e -
                    (mu + mu_s + alpha_s + lambda_t) * i_s -
                    (1.0 - q) * alpha_t * treat
                rhs_i_a = (1 - p) * delta_e * e -
                    (mu + mu_a + alpha_a) * i_a
                rhs_r = alpha_s * i_s + alpha_a * i_a +
                    q * alpha_t * treat - mu * r
                rhs_d = mu_s * i_s + mu_a * i_a
                rhs_v = lambda_v * s - epsilon * force_infection * v -
                    (mu + delta_v) * v
                rhs_treat = lambda_t * i_s -
                    (mu + alpha_t) * treat
                rhs = c(rhs_s,
                        rhs_e,
                        rhs_i_s,
                        rhs_i_a,
                        rhs_r,
                        rhs_d,
                        rhs_v,
                        rhs_treat)
                return(list(rhs))
            })
    }
    # Plot1: renderPlot to be passed to UI tab 1
    output$plot1 = renderPlot({
        parms = c(beta_s = input$beta_s,
                  beta_a = input$beta_a,
                  epsilon = input$epsilon,
                  delta_e = input$delta_e,
                  delta_v = input$delta_v,
                  p = input$p,
                  q = input$q,
                  alpha_a = input$alpha_a,
                  alpha_t = input$alpha_t,
                  alpha_s = input$alpha_s,
                  mu = 1.0 / input$mu_inverse,
                  mu_s = input$mu_s,
                  mu_a = input$mu_a,
                  lambda_v = input$lambda_v,
                  lambda_t = input$lambda_t,
                  n_whole = input$n_whole,
                  T = input$T
                )
####Initial Conditions####
        n_whole = input$n_whole
        initial_conditions = c(s = (n_whole - 2)/n_whole,
                               e = 0.0,
                               i_s = 1 / n_whole,
                               i_a = 1 / n_whole,
                               r = 0.0,
                               d = 0.0,
                               v = 0.0,
                               treat =  0.0
                               )
        times = seq(0, input$T[2], by = 1/1000)
        R0 = round(with(as.list(parms),
                         beta_s + beta_a /
                             (alpha_a + alpha_s + alpha_t + mu + mu_s + mu_a),
                         ),1
         )
        out = ode(y = initial_conditions,
                  times = times,
                  func = seirvt_mod,
                  parms = parms)
        out = as.data.frame(out)
    #Plot1
    p1 <- ggplot(data = out, aes(x = time, y = s)) +
        geom_line() +
        labs(y = "Suceptible fraction") +
        labs(x = "time") +
        labs(title = paste("R0=", R0))
    #
    p2 <- ggplot(data = out, aes(x = time, y = e)) +
        geom_line() +
        labs(y = "Exposed fraction") +
        labs(x = "time")
    #
    p3 <- ggplot(data = out, aes(x = time, y = i_s)) +
        geom_line() +
        labs(y = "Symptomatic infected") +
        labs(x = "time")
    #
    p4 <- ggplot(data = out, aes(x = time, y = i_a)) +
        geom_line() +
        labs(y = "Asymtomatic infected") +
        labs(x = "time")
    #
    p5 <- ggplot(data = out, aes(x = time, y = r)) +
        geom_line() +
        labs(y = "Recovered") +
        labs(x = "time")
    #
    p6 <- ggplot(data = out, aes(x = time, y = d)) +
        geom_line() +
        labs(y = "Death") +
        labs(x = "time")
    #
    p7 <- ggplot(data = out, aes(x = time, y = v)) +
        geom_line() +
        labs(y = "Vaccinated") +
        labs(x = "time")
    #
    p8 <-ggplot(data = out, aes(x = time, y = treat)) +
        geom_line() +
        labs(y = "Treats") +
        labs(x = "time")
    #
    p <- grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, nrow = 3)



####TODO: Make subpltos of each population####
    })
}
shinyApp(ui, server)
