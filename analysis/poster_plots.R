# R Script to generate plots for poster
# Author: Rasmus Hindström

# Load libraries
library(brms)
library(mia)
library(tidybayes)
library(ggplot2)
library(dplyr)
library(broom)
library(gt)
library(bayestestR)

# Set opts and seed
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
formula1 <- bf(
  response ~ (1 | group),
  sigma ~ (1 | group)
)

# Switch compile/update
update <- TRUE

if (update == FALSE) {
  mod <- brm(
    formula = formula1,
    data = df,
    family = student(),
    cores = 4,
  )
} else {
  mod <- update(
    mod,
    iter = 2e4,
    cores = 4,
    control = list(adapt_delta = 0.98)
  )
}

## All together
raw_means <- df %>%
  group_by(group) %>%
  summarise(raw_mean = mean(response))
grid <- df %>% modelr::data_grid(group)
means <- grid %>% add_epred_draws(mod)
pred <- grid %>% add_predicted_draws(mod)

theme_set(theme_linedraw(base_size = 16))
# Plot1 for ppc only
plot1 <- df %>%
  ggplot(aes(y = group, x = response)) +
  stat_interval(aes(x = .prediction),
    data = pred,
    size = 6
  ) +
  geom_point(size = 3) +
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
    .width = c(0.66, 0.95),
  ) +
  geom_point(
    aes(x = raw_mean, y = group),
    data = raw_means,
    shape = 5,
    size = 3,
    stroke = 1
  ) +
  geom_vline(xintercept = mean(df$response), linetype = "dashed") +
  labs(
    title = "Mean Shannon Index estimates",
    x = "Shannon Index",
    y = "Group",
    caption = NULL
  )

# Save the plots
ggsave("~/Projects/prTest/analysis/plots/ppc.png", plot1, width = 7, height = 6, dpi = 600)
ggsave("~/Projects/prTest/analysis/plots/group_means.png", plot2, width = 7, height = 6, dpi = 600)

# Annova tables
fit_aov <- aov(response ~ group, data = df)
anova_tidy <- broom::tidy(fit_aov) %>%
  filter(term != "Residuals") %>%
  mutate(
    statistic = round(statistic, 2),
    p.value = scales::pvalue(p.value)
  ) %>%
  select(Term = term, Df = df, `F` = statistic, `p-value` = p.value)

anova_tbl <- gt(anova_tidy) |>
  tab_header(title = md("**One-way ANOVA**: Shannon ~ Age Group"))

gt::gtsave(anova_tbl, "aovtable.html", inline_css = TRUE)
tryCatch(
  webshot2::webshot("aovtable.html", file = "aovtable.png", selector = ".gt_table", zoom = 2),
  error = function(e) webshot2::webshot("aovtable.html", file = "aovtable.png", zoom = 2)
)

# Tukey's HSD
fit_tukey <- TukeyHSD(fit_aov)
tukey_tidy <- broom::tidy(fit_tukey) %>%
  filter(term != "null.value") %>%
  mutate(ci = sprintf("[%.2f, %.2f]", conf.low, conf.high)) %>%
  select(Comparison = contrast, difference = estimate, `95% CI` = ci, `Adjusted p-value` = adj.p.value)

tukey_tbl <- gt(tukey_tidy) %>%
  tab_header(title = md("**Tukey's HSD**")) %>%
  fmt_number(
    columns = c(difference, `Adjusted p-value`),
    decimals = 3
  )

gt::gtsave(tukey_tbl, "tuktable.html", inline_css = TRUE)
tryCatch(
  webshot2::webshot("tuktable.html", file = "tuktable.png", selector = ".gt_table", zoom = 2),
  error = function(e) webshot2::webshot("tuktable.html", file = "tuktable.png", zoom = 2)
)

# Get random intercepts
re <- ranef(mod, summary = FALSE)$group[, , "Intercept"]
# Pairwise differences
diff_ma <- re[, "Middle_age"] - re[, "Adult"]
diff_ea <- re[, "Elderly"] - re[, "Adult"]
diff_em <- re[, "Elderly"] - re[, "Middle_age"]
# Probability of direction
pd_ma <- p_direction(diff_ma)
pd_ea <- p_direction(diff_ea)
pd_em <- p_direction(diff_em)
# Log2fc
post <- as_draws_df(mod)
post_a <- post$b_Intercept + post$"r_group[Adult,Intercept]"
post_m <- post$b_Intercept + post$"r_group[Middle_age,Intercept]"
post_e <- post$b_Intercept + post$"r_group[Elderly,Intercept]"
diff_ma <- post_m - post_a
diff_ea <- post_e - post_a
diff_em <- post_e - post_m
ci_ma <- quantile(diff_ma, probs = c(0.05, 0.95))
ci_ea <- quantile(diff_ea, probs = c(0.05, 0.95))
ci_em <- quantile(diff_em, probs = c(0.05, 0.95))
# Collect into a dataframe
tbl <- tibble(
  comp = c("Middle_age-Adult", "Elderly-Adult", "Elderly-Middle_age"),
  pd = c(pd_ma$pd, pd_ea$pd, pd_em$pd),
  diff = c(mean(diff_ma), mean(diff_ea), mean(diff_em)),
  ci_low = c(ci_ma[1], ci_ea[1], ci_em[1]),
  ci_upper = c(ci_ma[2], ci_ea[2], ci_em[2])
)


pr_table <- tbl %>%
  mutate(ci = sprintf("[%.2f, %.2f]", ci_low, ci_upper)) %>%
  select(comp, diff, ci, pd) %>%
  gt() %>%
  tab_header(
    title = md("**Probabilistic Comparisons of Groups**"),
  ) %>%
  fmt_number(
    columns = c(diff, pd),
    decimals = 3
  ) %>%
  cols_label(
    comp = "Comparison",
    diff = "Est. Difference",
    ci = "95% CrI",
    pd = "Pr-Direction"
  )

# Step 1: HTML (self-contained)
gt::gtsave(pr_table, "pr_table.html", inline_css = TRUE)

# Step 2: PNG via webshot2 (try selector first, then fallback)
tryCatch(
  webshot2::webshot("pr_table.html", file = "pr_table.png", selector = ".gt_table", zoom = 2),
  error = function(e) webshot2::webshot("pr_table.html", file = "pr_table.png", zoom = 2)
)
