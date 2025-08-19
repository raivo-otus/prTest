# R Script to generate plots for poster
# Author: Rasmus Hindström

# Load libraries
library(brms)
library(mia)
library(tidybayes)
library(ggplot2)
library(patchwork)
library(dplyr)

# Set opts and seed
set.seed(66642)
options(mc.cores = 4)

# Fetch dataset
data("GlobalPatterns", package = "mia")
data("peerj13075", package = "mia")

# Prepare data
tse <- peerj13075
tse <- addAlpha(tse, assay.type = "counts", index = "shannon")

# Cast to dataframe
df <- as_tibble(colData(tse))

# Rename vars
df <- df %>%
  rename(
    response = shannon,
    group = Age
  )

# Partial pooling model
ppM <- brm(
  formula = bf(
    response ~ (1 | group),
    sigma ~ (1 | group)
  ),
  data = df,
  family = student(),
  algorithm = "sampling",
  iter = 1e4,
  cores = 4,
  control = list(adapt_delta = 0.96), # Tighter acceptance

  # Cache model
  file = "model/partial_pooling"
)
