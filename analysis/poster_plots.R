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
data("peerj13075", package = "mia")

# Prepare data
tse <- peerj13075
tse <- addAlpha(tse, assay.type = "counts", index = "shannon")

# Cast to dataframe
df <- as_tibble(colData(tse))

# Prepare data
lvl <- c("Adult", "Middle_age", "Elderly")
df <- df %>%
  rename(
    response = shannon,
    group = Age
  ) %>%
  mutate(group = factor(group, level = lvl))

# Partial pooling model
mod <- brm(
  formula = bf(
    response ~ (1 | group),
    sigma ~ (1 | group)
  ),
  data = df,
  family = student(),
  algorithm = "sampling",
  iter = 2e4, # Longer chains
  cores = 4,
  control = list(adapt_delta = 0.98) # Tighter acceptance
)

## All together
raw_means <- df %>%
  group_by(group) %>%
  summarise(raw_mean = mean(response))
grid <- df %>% modelr::data_grid(group)
means <- grid %>% add_epred_draws(mod)
pred <- grid %>% add_predicted_draws(mod)

theme_set(theme_minimal(base_size = 16))
# Plot1 for ppc only
plot1 <- df %>%
  ggplot(aes(y = group, x = response)) +
  stat_interval(aes(x = .prediction), data = pred) +
  geom_point() +
  scale_color_brewer() +
  labs(
    title = "Posterior Predictive check",
    x = "Shannon Index",
    y = "Group",
    caption = NULL
  ) 

# plot2 for shrinkage
plot2 <- df %>%
  ggplot(aes(x = response, y = group)) +
  stat_pointinterval(
    aes(x = .epred),
    data = means,
    .width = c(0.66, 0.95)
  ) +
  geom_point(
    aes(x = raw_mean, y = group),
    data = raw_means,
    shape = 5,
    size = 3,
    stroke = 1.5
  ) +
  geom_vline(xintercept = mean(df$response), linetype = "dashed") +
  labs(
    title = "Mean Shannon Index estimates",
    x = "Shannon Index",
    y = "Group",
    caption = NULL
  )

# Save the plots
ggsave("~/Projects/prTest/analysis/plots/ppc.png", plot1, width = 7, height = 6, dpi = 300)
ggsave("~/Projects/prTest/analysis/plots/group_means.png", plot2, width = 7, height = 6, dpi = 300)
