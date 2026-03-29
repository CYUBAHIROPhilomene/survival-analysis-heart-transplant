rm(list = ls())

library(survival)
library(flexsurv)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(scales)
library(viridis)
library(gt)
library(patchwork)
library(rms)
library(survminer)
library(muhaz)

# DATA EXPLORATION

data(stanford2)

head(stanford2)
n_total    <- nrow(stanford2)
n_events   <- sum(stanford2$status == 1)
n_censored <- sum(stanford2$status == 0)

cat("Total observations:", n_total, "\n")
cat("Deaths (events):", n_events, "\n")
cat("Censored:", n_censored, "\n")


cat("\n Summary of Survival Time (days) \n")
cat("Min:",    min(stanford2$time),                "\n")
cat("Max:",    max(stanford2$time),                "\n")
cat("Mean:",   round(mean(stanford2$time), 2),     "\n")
cat("Median:", median(stanford2$time),             "\n")


cat("\n Summary of Age (years) \n")
cat("Min:",  min(stanford2$age),             "\n")
cat("Max:",  max(stanford2$age),             "\n")
cat("Mean:", round(mean(stanford2$age), 2),  "\n")


cat("\n Missing Values (Before Cleaning) \n")
print(colSums(is.na(stanford2)))



## DATA PREPARATION

table(stanford2$status)
stanford2_clean <- stanford2[!is.na(stanford2$time) & !is.na(stanford2$status), ]
cat("Observations after cleaning:", nrow(stanford2_clean), "\n")

# Handle missing values in t5
table(stanford2$status)

stanford2_clean <- stanford2[!is.na(stanford2$time) & !is.na(stanford2$status), ]
cat("Observations after cleaning:", nrow(stanford2_clean), "\n")

# t5 is right-skewed with outliers -> use median for imputation
cat("\nt5 before imputation:\n")
cat("Missing:", sum(is.na(stanford2_clean$t5)), "\n")
cat("Mean:",   round(mean(stanford2_clean$t5,   na.rm = TRUE), 4), "\n")
cat("Median:", round(median(stanford2_clean$t5, na.rm = TRUE), 4), "\n")

t5_median <- median(stanford2_clean$t5, na.rm = TRUE)
stanford2_clean$t5[is.na(stanford2_clean$t5)] <- t5_median

cat("\nMissing values after imputation:\n")
print(colSums(is.na(stanford2_clean)))


# DESCRIPTIVE STATISTICS TABLE

desc_stats <- function(x) {
  c(
    N       = sum(!is.na(x)),
    Missing = sum(is.na(x)),
    Mean    = mean(x, na.rm = TRUE),
    SD      = sd(x, na.rm = TRUE),
    Median  = median(x, na.rm = TRUE),
    Min     = min(x, na.rm = TRUE),
    Max     = max(x, na.rm = TRUE)
  )
}

vars <- c("time", "age", "t5")

desc_table <- lapply(stanford2_clean[vars], desc_stats) %>%
  do.call(rbind, .) %>%
  as.data.frame()

desc_table$Variable <- rownames(desc_table)
rownames(desc_table) <- NULL

desc_table <- desc_table %>%
  select(Variable, N, Missing, Mean, SD, Median, Min, Max) %>%
  mutate(
    across(c(Mean, SD, Median, Min, Max), ~ round(., 2)),
    `Mean (SD)` = paste0(Mean, " (", SD, ")")
  ) %>%
  select(Variable, N, Missing, `Mean (SD)`, Median, Min, Max)

desc_table %>%
  gt() %>%
  tab_header(
    title = md("**Table 1. Descriptive Statistics of Key Variables**")
  ) %>%
  cols_align(
    align = "center",
    columns = -Variable
  )



# DATA VISUALIZATION

# Create grouped status variable
stanford2_plot <- stanford2_clean %>%
  mutate(
    status_group = factor(
      status,
      levels = c(0, 1),
      labels = c("Censored", "Event")
    )
  )

# Survival time distribution
plot_time_facet <- ggplot(stanford2_plot, aes(x = time, fill = status_group)) +
  geom_histogram(
    bins = 15,
    color = "white",
    alpha = 0.8
  ) +
  facet_wrap(~ status_group, ncol = 2) +
  labs(
    title = "Survival Time Distribution for Each Status Group",
    x = "Survival Time (days)",
    y = "Count"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

plot_time_facet

# Age distribution 
plot_age_facet_loess <- ggplot(stanford2_plot, aes(x = age)) +
  geom_histogram(
    bins = 15,
    fill = "skyblue",
    color = "white",
    alpha = 0.8
  ) +
  geom_density(
    aes(y = after_stat(count)),
    color = "red",
    linewidth = 1
  ) +
  facet_wrap(~ status_group, ncol = 2) +
  labs(
    title = "Age Distribution for Each Status Group",
    x = "Age (years)",
    y = "Count"
  ) +
  theme_minimal(base_size = 14)

plot_age_facet_loess

ggsave(
  filename = "survival_time_facet.png",
  plot = plot_time_facet,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  filename = "age_facet_loess.png",
  plot = plot_age_facet_loess,
  width = 8,
  height = 5,
  dpi = 300
)


# FIT MULTIPLE PARAMETRIC MODELS

SurvObj <- Surv(time = stanford2_clean$time, event = stanford2_clean$status)

fits <- list(
  Exponential = flexsurvreg(SurvObj ~ 1, dist = "exp"),
  Weibull     = flexsurvreg(SurvObj ~ 1, dist = "weibull"),
  LogNormal   = flexsurvreg(SurvObj ~ 1, dist = "lnorm"),
  LogLogistic = flexsurvreg(SurvObj ~ 1, dist = "llogis"),
  Gamma       = flexsurvreg(SurvObj ~ 1, dist = "gamma"),
  GenGamma    = flexsurvreg(SurvObj ~ 1, dist = "gengamma")
)

fits

time_grid <- seq(0, max(stanford2_clean$time), length.out = 300)


param_surv_list <- lapply(names(fits), function(model_name) {
  s <- summary(fits[[model_name]], type = "survival", t = time_grid)[[1]]
  data.frame(
    time  = s$time,
    surv  = s$est,
    Model = model_name
  )
})

param_surv_df <- bind_rows(param_surv_list)


#  MODEL SELECTION (AIC)

aic_values    <- sapply(fits, AIC)
loglik_values <- sapply(fits, function(m) m$loglik)

model_comparison <- data.frame(
  Model         = names(fits),
  LogLikelihood = round(loglik_values, 3),
  AIC           = round(aic_values, 3)
)
model_comparison <- model_comparison[order(model_comparison$AIC), ]

cat("\n Model Comparison Table (sorted by AIC) \n")
print(model_comparison)


best_name  <- model_comparison$Model[1]
best_model <- fits[[best_name]]

best_surv <- summary(best_model, type = "survival", t = time_grid)[[1]]
best_df   <- data.frame(
  time = best_surv$time,
  surv = best_surv$est
)


# LIFE FUNCTIONS OF THE LOG-NORMAL MODEL

fit_lnorm <- flexsurvreg(SurvObj ~ 1, dist = "lnorm")
print(fit_lnorm)

mu    <- fit_lnorm$res["meanlog", "est"]
sigma <- fit_lnorm$res["sdlog",   "est"]

cat("\nMLE of meanlog (mu):",  round(mu, 3),    "\n")
cat("MLE of sdlog  (sigma):", round(sigma, 3), "\n")

time_seq <- seq(0, max(stanford2_clean$time), by = 10)

# Life functions
S_t <- 1 - pnorm((log(time_seq + 1e-10) - mu) / sigma)
f_t <- dlnorm(time_seq + 1e-10, meanlog = mu, sdlog = sigma)
h_t <- f_t / (S_t + 1e-10)
F_t <- 1 - S_t

# Mean and variance
mean_T <- exp(mu + sigma^2 / 2)
var_T  <- (exp(sigma^2) - 1) * exp(2 * mu + sigma^2)

cat("\nMean Survival Time:", round(mean_T, 2), "days\n")
cat("Variance:",           round(var_T, 2),  "days^2\n")

# Plots
par(mfrow = c(2, 2))

plot(time_seq, S_t, type = "l", col = "blue", lwd = 2,
     xlab = "Time (days)", ylab = "S(t)",
     main = "Survival Function (Log-Normal)")

plot(time_seq, h_t, type = "l", col = "red", lwd = 2,
     xlab = "Time (days)", ylab = "h(t)",
     main = "Hazard Function (Log-Normal)")

plot(time_seq, f_t, type = "l", col = "darkgreen", lwd = 2,
     xlab = "Time (days)", ylab = "f(t)",
     main = "Density Function (Log-Normal)")

plot(time_seq, F_t, type = "l", col = "purple", lwd = 2,
     xlab = "Time (days)", ylab = "F(t)",
     main = "CDF (Log-Normal)")

par(mfrow = c(1, 1))



# MLE PARAMETER ESTIMATES

cat("\n MLE Estimates for Log-Normal Model \n")
print(fit_lnorm$res)


# NON-PARAMETRIC ANALYSIS — KAPLAN-MEIER ESTIMATOR

# Overall KM Estimate
km_fit <- survfit(SurvObj ~ 1, data = stanford2_clean)
summary(km_fit)

# Key statistics
cat("\n--- KM Overall Summary ---\n")
cat("Median survival time:", km_fit$table["median"], "days\n")
cat("95% CI for median:  [",
    km_fit$table["0.95LCL"], ",",
    km_fit$table["0.95UCL"], "] days\n")


#the percentile time points from the observed data
t25 <- quantile(stanford2_clean$time, 0.25)
t50 <- quantile(stanford2_clean$time, 0.50)
t75 <- quantile(stanford2_clean$time, 0.75)

km_at_times <- summary(km_fit, times = c(t25, t50, t75))

surv_prob_table <- data.frame(
  Percentile       = c("25th", "50th", "75th"),
  Time_days        = c(t25, t50, t75),
  Survival_Prob    = round(km_at_times$surv, 4),
  Lower_95CI       = round(km_at_times$lower, 4),
  Upper_95CI       = round(km_at_times$upper, 4)
)

print(surv_prob_table)

surv_prob_table %>%
  gt() %>%
  tab_header(
    title = md("**Table: KM Survival Probabilities at Key Time Points**")
  ) %>%
  cols_label(
    Percentile    = "Percentile",
    Time_days     = "Time (days)",
    Survival_Prob = "S(t)",
    Lower_95CI    = "95% CI Lower",
    Upper_95CI    = "95% CI Upper"
  ) %>%
  cols_align(align = "center", columns = -Percentile)


# Stratified KM Curves by Age Group

stanford2_clean <- stanford2_clean %>%
  mutate(
    age_group = factor(
      ifelse(age <= 40, "Young (≤40)", "Older (>40)"),
      levels = c("Young (≤40)", "Older (>40)")
    )
  )

# Stratified KM fit
km_fit_grouped <- survfit(
  Surv(time, status) ~ age_group,
  data = stanford2_clean
)

print(km_fit_grouped)
summary(km_fit_grouped)$table


#  Overall KM Curve 
km_df <- data.frame(
  time  = km_fit$time,
  surv  = km_fit$surv,
  lower = km_fit$lower,
  upper = km_fit$upper
)

plot_km_overall <- ggplot(km_df, aes(x = time, y = surv)) +
  geom_step(color = "#1b7837", linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.15, fill = "#1b7837") +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             color = "red", linewidth = 0.7) +
  annotate("text", x = max(km_df$time) * 0.6, y = 0.52,
           label = paste("Median =", km_fit$table["median"], "days"),
           color = "red", size = 4) +
  labs(
    title    = "Kaplan-Meier Survival Curve (Overall)",
    subtitle = "Shaded region = 95% confidence interval",
    x        = "Time (days)",
    y        = "Survival Probability S(t)"
  ) +
  theme_minimal(base_size = 13)

plot_km_overall


# Stratified KM Curves by Age Group 
km_group_df <- data.frame(
  time      = km_fit_grouped$time,
  surv      = km_fit_grouped$surv,
  lower     = km_fit_grouped$lower,
  upper     = km_fit_grouped$upper,
  age_group = rep(names(km_fit_grouped$strata),
                  km_fit_grouped$strata)
)

# Clean up strata labels
km_group_df$age_group <- gsub("age_group=", "", km_group_df$age_group)

plot_km_stratified <- ggplot(km_group_df,
                             aes(x = time, y = surv,
                                 color = age_group,
                                 fill  = age_group)) +
  geom_step(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.12, color = NA) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             color = "grey40", linewidth = 0.6) +
  scale_color_manual(values = c("#2166ac", "#d6604d")) +
  scale_fill_manual(values  = c("#2166ac", "#d6604d")) +
  labs(
    title    = "Stratified Kaplan-Meier Curves by Age Group",
    subtitle = "Shaded region = 95% confidence interval",
    x        = "Time (days)",
    y        = "Survival Probability S(t)",
    color    = "Age Group",
    fill     = "Age Group"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

plot_km_stratified


ggsave("km_overall.png",    plot = plot_km_overall,    width = 8, height = 5, dpi = 300)
ggsave("km_stratified.png", plot = plot_km_stratified, width = 8, height = 5, dpi = 300)