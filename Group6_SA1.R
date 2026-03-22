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

# DATA EXPLORATION
data(stanford2)

head(stanford2)
#View(stanford2)

# Total observations
n_total    <- nrow(stanford2)
n_events   <- sum(stanford2$status == 1)
n_censored <- sum(stanford2$status == 0)

cat("Total observations:", n_total, "\n")
cat("Deaths (events):", n_events, "\n")
cat("Censored:", n_censored, "\n")

# Summary of survival time
cat(" Summary of Survival Time (days) \n")
cat("Min:",    min(stanford2$time),                "\n")
cat("Max:",    max(stanford2$time),                "\n")
cat("Mean:",   round(mean(stanford2$time), 2),     "\n")
cat("Median:", median(stanford2$time),             "\n")

# Summary of age
cat("\n Summary of Age (years) \n")
cat("Min:",  min(stanford2$age),             "\n")
cat("Max:",  max(stanford2$age),             "\n")
cat("Mean:", round(mean(stanford2$age), 2),  "\n")

# Missing values
cat("\n Missing Values \n")
print(colSums(is.na(stanford2)))


# DATA PREPARATION

table(stanford2$status)

# Remove rows with missing time or status only
stanford2_clean <- stanford2[!is.na(stanford2$time) & !is.na(stanford2$status), ]
cat("Observations after cleaning:", nrow(stanford2_clean), "\n")



# Descriptive statistics table

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
    title = md("**Table 1. Descriptive Statistics of Key Variables**")  ) %>%
  cols_align(
    align = "center",
    columns = -Variable
  )

gtsave(desc_table %>%
         gt() %>%
         tab_header(
           title = md("**Table 1. Descriptive Statistics of Key Variables**")         ) %>%
         cols_align(
           align = "center",
           columns = -Variable
         ),
       filename = "descriptive_table.png")




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

# Age distribution - separate panels + loess line
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


# Create survival object
SurvObj <- Surv(time = stanford2_clean$time, event = stanford2_clean$status)
#View(SurvObj)

# FIT MULTIPLE PARAMETRIC MODELS

fits <- list(
  Exponential = flexsurvreg(SurvObj ~ 1, dist = "exp"),
  Weibull     = flexsurvreg(SurvObj ~ 1, dist = "weibull"),
  LogNormal   = flexsurvreg(SurvObj ~ 1, dist = "lnorm"),
  LogLogistic = flexsurvreg(SurvObj ~ 1, dist = "llogis"),
  Gamma       = flexsurvreg(SurvObj ~ 1, dist = "gamma"),
  GenGamma    = flexsurvreg(SurvObj ~ 1, dist = "gengamma")
)

fits

# MODEL SELECTION (AIC)

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

# LIFE FUNCTIONS OF THE LOG-NORMAL MODEL

# Fit Log-Normal model
fit_lnorm <- flexsurvreg(SurvObj ~ 1, dist = "lnorm")
print(fit_lnorm)

# Extract parameters
mu    <- fit_lnorm$res["meanlog", "est"]
sigma <- fit_lnorm$res["sdlog",   "est"]

cat("\nMLE of meanlog (mu):",   round(mu, 3),    "\n")
cat("MLE of sdlog  (sigma):",  round(sigma, 3), "\n")

# Time sequence
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

# 2x2 plot
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
