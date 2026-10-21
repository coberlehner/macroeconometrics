library(zoo)
library(dplyr)
library(BVAR)
library(vars)
library(ggplot2)
library(gridExtra)
library(gt)

# ------------------------------------ #
# ------- Data preparation ----------- #
# ------------------------------------ #

# Load data
gdp <- read.csv("data/ECB Data Portal_gdp.csv")
cpi <- read.csv("data/ECB Data Portal_cpi.csv")
irt <- read.csv("data/ECB Data Portal_irt.csv")
unemp <- read.csv("data/ECB Data Portal_unemp.csv")
debt <- read.csv("data/ECB Data Portal_debt.csv")

colnames(irt)[3] <- "interest_rate"
colnames(cpi)[3] <- "infl"
colnames(gdp)[3] <- "growth"
colnames(unemp)[3] <- "unemp"
colnames(debt)[3] <- "debt"

# Change frequencies to quarterly by averaging and bring in same format
irt$DATE <- as.Date(irt$DATE)
irt <- aggregate(interest_rate ~ as.yearqtr(DATE), data = irt, FUN = mean) %>%
  rename("DATE" = "as.yearqtr(DATE)")

cpi$DATE <- as.Date(cpi$DATE)
cpi <- aggregate(infl ~ as.yearqtr(DATE), data = cpi, FUN = mean) %>% 
  rename("DATE" = "as.yearqtr(DATE)")

unemp$DATE <- as.Date(unemp$DATE)
unemp <- aggregate(unemp ~ as.yearqtr(DATE), data = unemp, FUN = mean) %>% 
  rename("DATE" = "as.yearqtr(DATE)")

gdp <- gdp %>% 
  mutate("DATE" = gsub("([0-9]{4})Q([1-4])", "\\1 Q\\2", TIME.PERIOD)) %>%
  dplyr::select(-TIME.PERIOD)

debt <- debt %>% 
  mutate("DATE" = gsub("([0-9]{4})Q([1-4])", "\\1 Q\\2", TIME.PERIOD)) %>%
  dplyr::select(-TIME.PERIOD)

# Merge data 
data <- Reduce(function(x, y) merge(x, y, by = "DATE"),
       list(gdp, cpi, irt, debt, unemp))
data <- na.omit(data)
rownames(data) <- NULL


# ---------------------------------------- #
# ------- Forecast evaluation ------------ #
# ---------------------------------------- #

# Split into M1 and M2
data_M1 <- data[,1:4]

# Select lag for Models
VARselect(data_M1[,2:4], lag.max = 15) # Selects lag 2
VARselect(data[,2:6], lag.max = 15) # Selects lag 2


# Transform data into ts
Y1 <- ts(data[, c("growth", "infl", "interest_rate")],
         start = c(2000,1),
         frequency = 4)
Y2 <- ts(data[, c("growth", "infl", "interest_rate", "unemp", "debt")],
         start = c(2000,1),
         frequency = 4)


# Create function that runs forecast evaluation and calculates RMSEs for both models
forecast_rmse <- function(Y, p, holdout = 8){
  
  n <- nrow(Y)
  k <- ncol(Y)
  
  # Create data matrices 
  forecasts <- matrix(NA, holdout, k)
  actuals   <- matrix(NA, holdout, k)
  colnames(forecasts) <- colnames(Y)
  colnames(actuals)   <- colnames(Y)
  
  for(i in 1:holdout){
    # Expanding estimation sample
    train <- Y[1:(n - holdout + i - 1), ]
    # Estimate Bayesian VAR
    fit <- bvar(train, lags = p)
    # One-step-ahead forecast
    fcst <- predict(fit, horizon = 1)
    # Posterior mean forecast
    forecasts[i, ] <- apply(fcst$fcast[, 1, ], 2, mean)
    # Actual realization
    actuals[i, ] <- Y[n - holdout + i, ]
  }
  
  # Compute RMSE
  rmse <- sqrt(colMeans((actuals - forecasts)^2))
  
  list(forecasts = forecasts,
       actuals = actuals,
       RMSE = rmse)
}

# compare both models in terms of RMSE as the chosen measure
results_M1 <- forecast_rmse(Y1, 2)
results_M2 <- forecast_rmse(Y2, 2)

comparison <- data.frame(Variable = c("Growth", "Inflation", "Interest rate"),
                         RMSE_M1 = results_M1$RMSE,
                         RMSE_M2 = results_M2$RMSE[1:3])

comparison$Winner <- ifelse(comparison$RMSE_M1 < comparison$RMSE_M2, "M1", "M2") # chooses M1


# Export table
comparison %>%
  gt() %>%
  gtsave("RMSE_comparison.html")


# Produce quarterly forecasts of GDP and Inflation based on the chosen model M1
forecast_mod <- bvar(Y1, lags = 2)
forecast <- predict(forecast_mod, h=8)

# ---------------------------------- #
# -------- GDP forecast ------------ #
# ---------------------------------- #

gdp_data <- data.frame(Time = 100:(nrow(gdp)-1),
                       Period = gdp[100:(nrow(gdp)-1),]$DATE,
                       Growth = gdp[100:(nrow(gdp)-1),]$growth,
                       Type = "data")
forecast_gdp <- data.frame(Time = nrow(gdp):(nrow(gdp)+7),
                           Period = c("2026 Q1", "2026 Q2","2026 Q3","2026 Q4",
                                      "2027 Q1","2027 Q2","2027 Q3","2027 Q4"),
                           Growth = forecast$quants[2, , 1],
                           Lower = forecast$quants[1, , 1],
                           Upper = forecast$quants[3, , 1],
                           Type = "forecast")
gdp_line <- data.frame(Time = c(gdp_data$Time, forecast_gdp$Time),
                        Growth = c(gdp_data$Growth, forecast_gdp$Growth),
                        Type = c(rep("Historical", nrow(gdp_data)),
                                 rep("Forecast", nrow(forecast_gdp))))
gdp_plot <- ggplot() +
  # confidence intervals of forecast
  geom_ribbon(data = forecast_gdp, aes(x = Time, ymin = Lower, ymax = Upper), 
                  fill = "purple", alpha = 0.2) +
  # actual data and forecast
  geom_line(data = gdp_line, aes(x = Time, y = Growth), linewidth = 0.75) +
  geom_point(data = forecast_gdp, aes(x = Time, y = Growth), size = 2, color= "purple4") +
  # axis labels
  scale_x_continuous(breaks = seq(min(gdp_data$Time), max(forecast_gdp$Time), by = 2),
                     labels = c(gdp_data$Period, forecast_gdp$Period)[seq(1, length(c(gdp_data$Period, forecast_gdp$Period)), by=2)]) +
  labs(title = "Real GDP Growth Forecast", x = "Quarter", y = "Growth (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 20))

# ---------------------------------------- #
# -------- Inflation forecast ------------ #
# ---------------------------------------- #

infl_data <- data.frame(Time = 93:(nrow(cpi)-1),
                        Period = as.character(cpi[93:(nrow(cpi)-1),]$DATE),
                        Inflation = cpi[93:(nrow(cpi)-1),]$infl,
                        Type = "data")
forecast_infl <- data.frame(Time = nrow(cpi):(nrow(cpi)+7),
                            Period = c("2026 Q1", "2026 Q2","2026 Q3","2026 Q4",
                                        "2027 Q1","2027 Q2","2027 Q3","2027 Q4"),
                            Inflation = forecast$quants[2, , 2],
                            Lower = forecast$quants[1, , 2],
                            Upper = forecast$quants[3, , 2],
                            Type = "forecast")
infl_line <- data.frame(Time = c(infl_data$Time, forecast_infl$Time),
                       Inflation = c(infl_data$Inflation, forecast_infl$Inflation),
                       Type = c(rep("Historical", nrow(infl_data)),
                                rep("Forecast", nrow(forecast_infl))))
infl_plot <- ggplot() +
  # confidence intervals of forecast
  geom_ribbon(data = forecast_infl, aes(x = Time, ymin = Lower, ymax = Upper), 
              fill = "purple", alpha = 0.2) +
  # actual data and forecast
  geom_line(data = infl_line, aes(x = Time, y = Inflation), linewidth = 0.75) +
  geom_point(data = forecast_infl, aes(x = Time, y = Inflation), size = 2, color= "purple4") +
  # axis labels
  scale_x_continuous(breaks = seq(min(infl_data$Time), max(forecast_infl$Time), by = 2),
                     labels = c(infl_data$Period, forecast_infl$Period)[seq(1, length(c(infl_data$Period, forecast_infl$Period)), by=2)]) +
  labs(title = "Inflation Forecast", x = "Quarter", y = "Inflation (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 20))


grid.arrange(gdp_plot, infl_plot, ncol = 2)
