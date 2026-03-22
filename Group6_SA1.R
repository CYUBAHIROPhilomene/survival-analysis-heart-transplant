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

