# Norovirus SEIRS Modeling Tool

Repository: <https://github.com/geithkel/epid799_noro_seirs_tool>

This repository contains an interactive R Shiny application for exploring simplified norovirus transmission dynamics using an SEIRS compartmental model.

## Purpose

The app allows users to adjust epidemiologic inputs such as population size, initial exposed and infectious individuals, incubation period, infectious period, immunity duration, and basic reproduction number. It then displays simulated susceptible, exposed, infectious, and recovered populations over time.

## Model

The deterministic SEIRS model is:

\[
\frac{dS}{dt} = -\beta \frac{SI}{N} + \omega R
\]

\[
\frac{dE}{dt} = \beta \frac{SI}{N} - \sigma E
\]

\[
\frac{dI}{dt} = \sigma E - \gamma I
\]

\[
\frac{dR}{dt} = \gamma I - \omega R
\]

Where:

| Parameter | Meaning |
|---|---|
| `S` | Susceptible individuals |
| `E` | Exposed/incubating individuals |
| `I` | Infectious individuals |
| `R` | Recovered/temporarily immune individuals |
| `beta` | Transmission rate |
| `sigma` | Rate of progression from exposed to infectious |
| `gamma` | Recovery rate |
| `omega` | Waning immunity rate |
| `N` | Total population |

The app derives rates as:

- `sigma = 1 / incubation_period`
- `gamma = 1 / infectious_period`
- `omega = 1 / immunity_duration`
- `beta = R0 * gamma`

## Default values

The starting defaults are intended for educational exploration and are not outbreak-specific estimates.

| Input | Default |
|---|---:|
| Population size | 1,000 |
| Initial exposed | 5 |
| Initial infectious | 5 |
| Initial recovered | 0 |
| Basic reproduction number, R0 | 2.0 |
| Incubation period | 1.2 days |
| Infectious period | 2.0 days |
| Immunity duration | 180 days |
| Simulation length | 120 days |

## Files

- `tool_code.R`: Main Shiny app and SEIRS model code.
- `app.R`: Shiny entry point that sources `tool_code.R`.
- `.gitignore`: Files and folders excluded from version control.

## Run locally

Install required packages in R:

```r
install.packages(c("shiny", "deSolve", "ggplot2", "dplyr", "tidyr"))
```

Then run the app from the repository root:

```r
shiny::runApp()
```

Alternatively, open `app.R` or `tool_code.R` in RStudio/Positron and click **Run App**.

## Deployment

This app can be deployed to shinyapps.io or Posit Connect. For shinyapps.io:

```r
install.packages("rsconnect")
rsconnect::deployApp()
```

After deployment, add the live app URL here.

## Disclaimer

This model is for educational and exploratory use only. It is not intended for clinical decision-making, outbreak response, or public health policy without expert review and validation. Norovirus transmission can include environmental persistence, asymptomatic infection, heterogeneity in contact patterns, and setting-specific control measures that are not represented in this first simplified model.
