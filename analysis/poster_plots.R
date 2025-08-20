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

# Prepare data
lvl <- c("Adult", "Middle_age", "Elderly")
df <- df %>%
  rename(
    response = shannon,
    group = Age
  ) %>%
  mutate(group = factor(group, level = lvl))

# Partial pooling model
ppM <- brm(
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
means <- grid %>% add_epred_draws(ppM)
pred <- grid %>% add_predicted_draws(ppM)

# plot caption
cap <- stringr::str_wrap(
  paste(
    "Fig.1: Mean Shannon index estimates for groups. Dashed vertical line is",
    "total population mean. Point interval bars show model estimates, thicker",
    "line represents 66% Credible Interval (CI), the thinner line represents ",
    "95% CI. Original datapoints are displayed as dots, with model predictive",
    "intervals. Vertical ticks are the raw means from data."
  )
)

theme_set(theme_minimal(base_size = 16))
plot1 <- df %>%
  ggplot(aes(y = group, x = response)) +
  stat_interval(aes(x = .prediction), data = pred) +
  stat_pointinterval(
    aes(x = .epred),
    data = means,
    .width = c(0.66, 0.95),
    position = position_nudge(y = -0.1)
  ) +
  geom_point() +
  geom_point(
    aes(x = raw_mean, y = group),
    data = raw_means,
    shape = 3,
    size = 4,
    stroke = 1.5,
    position = position_nudge(y = -0.1)
  ) +
  scale_color_brewer() +
  geom_vline(xintercept = mean(df$response), linetype = "dashed") +
  labs(
    title = "Mean Shannon index estimates",
    subtitle = "Partial pooling model",
    x = "Shannon Index",
    y = "Group",
    caption = cap
  ) +
  theme(plot.caption = element_text(hjust = 0))

plot1
ggsave("plots/group_means.png", plot1, width = 7, height = 6, dpi = 300)
ggsave("plots/group_means.pdf", plot1, width = 7, height = 6)
