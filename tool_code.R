# Norovirus SEIRS Modeling Tool
# R Shiny application code lives in this file. app.R sources this file for deployment.

library(shiny)
library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)

seirs_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    N <- S + E + I + R

    dS <- -beta * S * I / N + omega * R
    dE <- beta * S * I / N - sigma * E
    dI <- sigma * E - gamma * I
    dR <- gamma * I - omega * R
    dCumulative_Infections <- sigma * E

    list(c(dS, dE, dI, dR, dCumulative_Infections))
  })
}

run_seirs <- function(population,
                      initial_exposed,
                      initial_infectious,
                      initial_recovered,
                      r0,
                      incubation_period,
                      infectious_period,
                      immunity_duration,
                      simulation_days,
                      time_step = 0.1) {
  initial_susceptible <- population - initial_exposed - initial_infectious - initial_recovered

  if (initial_susceptible < 0) {
    stop("Initial exposed, infectious, and recovered counts must not exceed population size.")
  }

  parameters <- c(
    beta = r0 * (1 / infectious_period),
    sigma = 1 / incubation_period,
    gamma = 1 / infectious_period,
    omega = 1 / immunity_duration
  )

  state <- c(
    S = initial_susceptible,
    E = initial_exposed,
    I = initial_infectious,
    R = initial_recovered,
    Cumulative_Infections = initial_infectious
  )

  times <- seq(0, simulation_days, by = time_step)

  output <- deSolve::ode(
    y = state,
    times = times,
    func = seirs_model,
    parms = parameters
  )

  as.data.frame(output)
}

summarize_seirs <- function(results) {
  peak_row <- results[which.max(results$I), ]
  final_row <- results[nrow(results), ]

  data.frame(
    Metric = c(
      "Peak infectious count",
      "Day of peak infectious count",
      "Cumulative infections",
      "Final susceptible count",
      "Final exposed count",
      "Final infectious count",
      "Final recovered count"
    ),
    Value = round(c(
      peak_row$I,
      peak_row$time,
      final_row$Cumulative_Infections,
      final_row$S,
      final_row$E,
      final_row$I,
      final_row$R
    ), 2),
    check.names = FALSE
  )
}

make_seirs_plot <- function(results) {
  results |>
    tidyr::pivot_longer(
      cols = c(S, E, I, R),
      names_to = "Compartment",
      values_to = "Count"
    ) |>
    ggplot2::ggplot(ggplot2::aes(x = time, y = Count, color = Compartment)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::labs(
      title = "Norovirus SEIRS model",
      x = "Time, days",
      y = "Number of people",
      color = "Compartment"
    ) +
    ggplot2::theme_minimal(base_size = 13)
}

make_infectious_plot <- function(results) {
  ggplot2::ggplot(results, ggplot2::aes(x = time, y = I)) +
    ggplot2::geom_line(color = "#B22222", linewidth = 1.1) +
    ggplot2::labs(
      title = "Currently infectious individuals over time",
      x = "Time, days",
      y = "Infectious individuals"
    ) +
    ggplot2::theme_minimal(base_size = 13)
}

ui <- shiny::fluidPage(
  shiny::titlePanel("Norovirus SEIRS Modeling Tool"),

  shiny::sidebarLayout(
    shiny::sidebarPanel(
      shiny::helpText("Adjust model inputs to explore simplified norovirus transmission dynamics."),
      shiny::numericInput("population", "Population size", value = 1000, min = 1, step = 1),
      shiny::numericInput("initial_exposed", "Initial exposed", value = 5, min = 0, step = 1),
      shiny::numericInput("initial_infectious", "Initial infectious", value = 5, min = 0, step = 1),
      shiny::numericInput("initial_recovered", "Initial recovered", value = 0, min = 0, step = 1),
      shiny::sliderInput("r0", "Basic reproduction number, R0", min = 0.1, max = 6, value = 2, step = 0.1),
      shiny::numericInput("incubation_period", "Incubation period, days", value = 1.2, min = 0.1, step = 0.1),
      shiny::numericInput("infectious_period", "Infectious period, days", value = 2.0, min = 0.1, step = 0.1),
      shiny::numericInput("immunity_duration", "Immunity duration, days", value = 180, min = 1, step = 1),
      shiny::numericInput("simulation_days", "Simulation length, days", value = 120, min = 1, step = 1),
      shiny::downloadButton("download_results", "Download CSV")
    ),

    shiny::mainPanel(
      shiny::h3("Model output"),
      shiny::plotOutput("seirs_plot", height = "420px"),
      shiny::plotOutput("infectious_plot", height = "320px"),
      shiny::h3("Summary"),
      shiny::tableOutput("summary_table"),
      shiny::h3("Model equations"),
      shiny::withMathJax(
        shiny::p("$$\\frac{dS}{dt} = -\\beta SI/N + \\omega R$$"),
        shiny::p("$$\\frac{dE}{dt} = \\beta SI/N - \\sigma E$$"),
        shiny::p("$$\\frac{dI}{dt} = \\sigma E - \\gamma I$$"),
        shiny::p("$$\\frac{dR}{dt} = \\gamma I - \\omega R$$")
      ),
      shiny::h3("Disclaimer"),
      shiny::p("This simplified deterministic SEIRS model is for educational and exploratory use only. It is not intended for clinical decision-making, outbreak response, or public health policy without expert review and validation.")
    )
  )
)

server <- function(input, output, session) {
  results <- shiny::reactive({
    shiny::validate(
      shiny::need(
        input$initial_exposed + input$initial_infectious + input$initial_recovered <= input$population,
        "Initial exposed, infectious, and recovered counts must not exceed population size."
      ),
      shiny::need(input$incubation_period > 0, "Incubation period must be greater than 0."),
      shiny::need(input$infectious_period > 0, "Infectious period must be greater than 0."),
      shiny::need(input$immunity_duration > 0, "Immunity duration must be greater than 0.")
    )

    run_seirs(
      population = input$population,
      initial_exposed = input$initial_exposed,
      initial_infectious = input$initial_infectious,
      initial_recovered = input$initial_recovered,
      r0 = input$r0,
      incubation_period = input$incubation_period,
      infectious_period = input$infectious_period,
      immunity_duration = input$immunity_duration,
      simulation_days = input$simulation_days
    )
  })

  output$seirs_plot <- shiny::renderPlot({
    make_seirs_plot(results())
  })

  output$infectious_plot <- shiny::renderPlot({
    make_infectious_plot(results())
  })

  output$summary_table <- shiny::renderTable({
    summarize_seirs(results())
  })

  output$download_results <- shiny::downloadHandler(
    filename = function() {
      paste0("norovirus_seirs_results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(results(), file, row.names = FALSE)
    }
  )
}

app <- shiny::shinyApp(ui = ui, server = server)
