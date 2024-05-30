#Leafs_data Residual
#Initialise workspace and data
set.seed(123)
library(knitr)
library(moderndive)
library(caret)
library(ggplot2)
library(lattice)
library(splines)
library(tidyverse)
library(readxl)
leafs_data <- read_excel("Downloads/Matematik KU /Bachelor/leafs.data.xlsx")
B <- 999
k <- 5
M <- 1
colnames(leafs_data)<- c("sc", "bfkg")
folds <- createFolds(1:nrow(leafs_data), k = k)
fold_results <- list()
r_squared_values<-numeric(k)
mse_values<-numeric(k)
suppressWarnings({
for (fold_idx in 1:k) {
  fold <- folds[[fold_idx]]
  train <- leafs_data[-fold, ]
  test <- leafs_data[fold, ]
  
  # Fit model 
  model <- lm(bfkg ~ bs(sc, knots=c(quantile(train$sc, probs=c(0.1,0.5,0.9))),
                        degree=3), data = train)
  
  # Sample residuals 
  residuals_model <- residuals(model)
  fitted_model <- fitted(model)
  
  # Generate bootstrap sets and models
  bootstrap_models <- lapply(1:B, function(b) {
    resampled_residuals <- sample(residuals_model, replace = TRUE)
    bootstrap_data <- data.frame(
      sc = train$sc,
      bfkg = fitted_model + resampled_residuals
    )
    lm(bfkg ~ bs(sc, knots=c(quantile(bootstrap_data$sc, probs=c(0.1, 0.5, 0.9)))
                 ,degree=3), data = bootstrap_data)
  })
  
  # Calculate residuals and r_bar
  x_bar <- mean(train$sc)
  ssx <- sum((train$sc - x_bar)^2)
  h <- (1 / nrow(train)) + ((train$sc - x_bar)^2) / ssx
  r <- residuals_model / sqrt(1 - h)
  r <- na.omit(r)
  r_bar <- mean(r)
  
  alpha <- 0.025
  quantile_high <- (B * M + 1) * alpha
  quantile_low <- (B * M + 1) * (1 - alpha)
  
  lower_bounds_fold <- numeric(nrow(test))
  upper_bounds_fold <- numeric(nrow(test))
  predictions_fold <- numeric(nrow(test))
  
  for (j in 1:nrow(test)) {
    prediction_errors <- sapply(1:B, function(b) {
      predict(bootstrap_models[[b]], newdata = data.frame(sc = test$sc[j])) -
        (predict(model, newdata = data.frame(sc = test$sc[j])) +
           sample(r - r_bar, M)) #sample(residuals(bootstrap_models[[b]]),M)
    })
    sorted_prediction_errors <- sort(prediction_errors)
    y_hat <- predict(model, newdata = data.frame(sc = test$sc[j]))
    lower_bounds_fold[j] <- y_hat - sorted_prediction_errors[quantile_low]
    upper_bounds_fold[j] <- y_hat - sorted_prediction_errors[quantile_high]
    predictions_fold[j] <- y_hat
  }
  
  fold_results[[fold_idx]] <- data.frame(
    sc = test$sc,
    bfkg = test$bfkg,
    lower = lower_bounds_fold,
    upper = upper_bounds_fold,
    prediction = predictions_fold,
    fold = fold_idx
  )
  # Calculate R-squared for this fold 
  ss_total<- sum((test$bfkg-mean(test$bfkg))^2)
  ss_residual<- sum((test$bfkg-predictions_fold)^2)
  r_squared_values[fold_idx]<-1-(ss_residual/ss_total)
  #Calculate MSE for this fold 
  mse_values[fold_idx]<- mean((test$bfkg-predictions_fold)^2)
}
})
# Combine results into a single data frame for easier plotting
results <- do.call(rbind, fold_results)

# Calculate coverage probability
results$within_interval <- with(results, bfkg >= lower & bfkg <= upper)
coverage_probability <- mean(results$within_interval)

# Calculate average interval width
results$interval_width <- results$upper - results$lower
average_interval_width <- mean(results$interval_width)

#Calculate R-squared 
average_r_squared<- mean(r_squared_values)

#Calculate MSE
overall_mse<- mean(mse_values)

# Print the evaluation metrics
cat("Coverage Probability:", coverage_probability, "\n")
cat("Average Interval Width:", average_interval_width, "\n")
cat("R-squared:", average_r_squared,"\n")
cat("MSE:", overall_mse,"\n")

# Create a combined plot with facets for each fold
ggplot(results, aes(x = sc, y = bfkg)) +
  geom_point(color = "red", alpha = 0.5, size=0.1) +  # Data points
  geom_line(aes(y = prediction), color = "blue", size = 0.5) +  # Fitted line
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "green") +  # Prediction intervals
  labs(title = "Prediction Intervals for leafs_data using Residual Bootstrap for Each Fold",
       x = "Predictor (sc)",
       y = "Response (bfkg)") +
  theme_minimal() +
  facet_wrap(~ fold)
