# nested_beta_binomial_pipeline.R
# robust nested (county within state) Beta-Binomial pipeline (grid approx per-state)

library(tidyverse)
library(matrixStats)
library(ggplot2)
library(lubridate)

# reproducible random draws for later
set.seed(2025)

# setup

data_file <- "stp505_pivot_clean_dataset.csv"  
output_summary_csv <- "beta_binomial_posteriors_summary.csv"
output_samples_rds <- "beta_binomial_samples.rds"
plot_prefix <- ""          # optional
S_draws <- 2000            # posterior draws per-state and rows in theta_draws_matrix (moderate)
grid_length_default <- 120 # resolution (size) of grid approximation that estimates state level beta hypers
# plot sizing variables
plot_width <- 10; plot_height <- 6; dpi <- 300

# waic and independent draws utility/helper functions

# compute WAIC given a log-likelihood matrix (S x N), rows=draws, cols=obs
# numerically stable waic using log-sum-exp and pointwise variance for p_waic
compute_waic <- function(log_lik_mat) {
  S <- nrow(log_lik_mat); N <- ncol(log_lik_mat)
  col_log_sum_exp <- matrixStats::colLogSumExps(log_lik_mat)
  lppd <- sum(col_log_sum_exp - log(S))
  p_waic <- sum(apply(log_lik_mat, 2, var))
  waic <- -2 * (lppd - p_waic)
  list(waic = waic, elpd = lppd - p_waic, p_waic = p_waic)
}

# independent posterior draws (Beta(1,1) per county)
# generates S draws of county level posterior thetas assuming independent priors
# returns SxJ matrix
# alpha and beta are clamped to smallpositive numbers to avoid invalid beta params
make_independent_draws <- function(y_vec, n_vec, S = S_draws) {
  mat <- matrix(NA_real_, nrow = S, ncol = length(y_vec))
  for (j in seq_along(y_vec)) {
    a <- 1 + y_vec[j]; b <- 1 + n_vec[j] - y_vec[j]
    a <- pmax(a, 1e-8); b <- pmax(b, 1e-8)
    mat[, j] <- rbeta(S, shape1 = a, shape2 = b)
  }
  mat
}

# approximates state-level posteriors for alpha_s, beta_s using grid parameterized in log_ratio and log_sum for numeric stability
# computes marginal loglikelihood of data for each grid cell (vector) and adds a weak prior on alpha+beta then converts to sampling weights
# samples s grid cells by posterior weights and for each sampled (a,b) draws J county theta values from Beta(a+y, beta+n-y)
# returns SxJ matrix of theta draws (columns correspond to counties in that state)
# if we get non-finite values NULL is returned 
refit_state_beta_grid <- function(state_name, df, S = S_draws, grid_length = grid_length_default,
                                  log_ratio_range = c(-3, 3),
                                  log_sum_range = c(log(0.5), log(200))) {
  idx_state <- which(df$state == state_name)
  yv <- df$y[idx_state]
  nv <- df$n[idx_state]
  J <- length(yv)
  if (J == 0) stop("State not found:", state_name)
  
  # single-county fallback: Beta posterior with flat prior
  if (J == 1) {
    mat <- matrix(NA_real_, nrow = S, ncol = 1)
    for (s in seq_len(S)) mat[s, 1] <- rbeta(1, 1 + yv, 1 + nv - yv)
    return(mat)
  }
  
  # Build grid (log_ratio, log_sum) -> alpha,beta
  log_ratio_grid <- seq(log_ratio_range[1], log_ratio_range[2], length.out = grid_length)
  log_sum_grid <- seq(log_sum_range[1], log_sum_range[2], length.out = grid_length)
  grid <- expand.grid(log_ratio = log_ratio_grid, log_sum = log_sum_grid)
  
  # parameterize alpha and beta safely (clamp to small positive)
  alpha <- pmax(exp(grid$log_sum) * exp(grid$log_ratio) / (1 + exp(grid$log_ratio)), 1e-12)
  beta <- pmax(exp(grid$log_sum) / (1 + exp(grid$log_ratio)), 1e-12)
  
  # vectorized marginal log likelihood computation
  # create matrices (grid_cells x J)
  A <- matrix(alpha, nrow = length(alpha), ncol = J)
  B <- matrix(beta, nrow = length(beta), ncol = J)
  Y <- matrix(rep(yv, each = length(alpha)), nrow = length(alpha))
  Nmat <- matrix(rep(nv, each = length(alpha)), nrow = length(alpha))
  
  # compute log terms
  # lgamma(A + Y) + lgamma(B + N - Y) - lgamma(A + B + N)
  log_terms <- lgamma(A + Y) + lgamma(B + Nmat - Y) - lgamma(A + B + Nmat)
  log_marg <- rowSums(log_terms) - J * (lgamma(alpha) + lgamma(beta) - lgamma(alpha + beta))
  
  # weak prior on total strength alpha+beta
  log_prior <- -2.5 * log(alpha + beta + 1e-12)
  #log_prior <- 0
  log_post_unnorm <- log_prior + log_marg
  
  # stabilize with max and compute weights
  maxlog <- max(log_post_unnorm, na.rm = TRUE)
  if (!is.finite(maxlog)) return(NULL)
  w_unnorm <- exp(log_post_unnorm - maxlog)
  if (!is.finite(sum(w_unnorm)) || sum(w_unnorm) == 0) return(NULL)
  w <- w_unnorm / sum(w_unnorm)
  
  # sample grid indices and draw theta per sample
  set.seed(2025)
  samp_idx <- sample(seq_along(w), size = S, replace = TRUE, prob = w)
  alpha_s <- alpha[samp_idx]; beta_s <- beta[samp_idx]
  alpha_s <- pmax(alpha_s, 1e-12); beta_s <- pmax(beta_s, 1e-12)
  
  theta_mat <- matrix(NA_real_, nrow = S, ncol = J)
  for (s in seq_len(S)) {
    theta_mat[s, ] <- rbeta(J, shape1 = alpha_s[s] + yv, shape2 = beta_s[s] + nv - yv)
  }
  
  if (any(!is.finite(theta_mat))) return(NULL)
  theta_mat
}


# reads cleaned data and checks col names creating standardized cols, then orders by state, fips, and adds row_id for tracing
# we remove all commas
# pop is numeric
# pop_num is rounded to an integer 
df <- read_csv(data_file, show_col_types = FALSE)

if (!all(c("FIPS","Recip_County","Recip_State","Series_Complete_Yes","pop") %in% names(df))) {
  stop("Data missing required columns. Required: FIPS, Recip_County, Recip_State, Series_Complete_Yes, pop")
}

df <- df %>%
  mutate(
    FIPS = as.character(FIPS),
    county = as.character(Recip_County),
    state = as.character(Recip_State),
    y = as.integer(str_replace_all(Series_Complete_Yes, ",", "")),
    pop_num = as.numeric(pop),
    n = as.integer(round(pop_num))
  ) %>%
  arrange(state, FIPS) %>%
  mutate(row_id = row_number())

# sanity checks, not needed but i was having issues at first
if (any(is.na(df$y))) stop("Some y are NA after coercion; inspect data.")
if (any(is.na(df$n))) stop("Some n are NA after coercion; inspect data.")
if (any(df$y < 0)) stop("Some y < 0.")
if (any(df$n < 1)) stop("Some n < 1. You may want to inspect tiny denominators.")

# if rounding created y > n for any rows above, raise n to y (preserve count)
over_idx <- which(df$y > df$n)
if (length(over_idx) > 0) {
  message("Adjusting ", length(over_idx), " rows where y > n by setting n = y to preserve counts.")
  df$n[over_idx] <- df$y[over_idx]
}

# records sizing for sanity 
N_counties <- nrow(df)
states <- unique(df$state)
S_states <- length(states)
cat("Loaded data: counties =", N_counties, "states =", S_states, "\n")

# prepares storage (list and matrix (SxN))
state_results <- list()          # store per-state metadata
theta_draws_matrix <- matrix(NA_real_, nrow = S_draws, ncol = N_counties)

# for each state we use grid to put the theta draws into the global matrix
pb <- txtProgressBar(min = 0, max = length(states), style = 3)
col_index_by_state <- split(seq_len(N_counties), df$state)

for (i in seq_along(states)) {
  st <- states[i]
  setTxtProgressBar(pb, i)
  idxs <- col_index_by_state[[st]]
  theta_mat <- tryCatch(
    refit_state_beta_grid(st, df, S = S_draws, grid_length = grid_length_default,
                          log_ratio_range = c(-3,3), log_sum_range = c(log(0.5), log(200))),
    error = function(e) NULL,
    warning = function(w) NULL
  )
  
  state_results[[st]] <- list(method = method_used, theta_samples = theta_mat)
  theta_draws_matrix[, idxs] <- theta_mat
}
close(pb)
cat("\nPer-state fitting complete.\n")


# this computes posterior means and 95% CIs for each county under hierarchial model (theta_draws matrix) and indpendent (beta(1,1)) model (ind_theta_samples)
# results are numeric vectors aligned to rows of df
# just a saved summary of hierarchial and independent post distributions 
theta_mean <- colMeans(theta_draws_matrix)
theta_CI <- t(apply(theta_draws_matrix, 2, quantile, probs = c(0.025, 0.975)))

# independent (no-pooling) posterior draws and summaries
ind_theta_samples <- make_independent_draws(df$y, df$n, S = S_draws)
ind_post_mean <- colMeans(ind_theta_samples)
ind_post_CI <- t(apply(ind_theta_samples, 2, quantile, probs = c(0.025, 0.975)))

# for each post draw we get predicted counts Y_tilde ~ binom(n, theta) for every county (hierarchial and independent)
# computes predictive means and 95% predictive intervals counts)
# posterior predictive draws (counts) for hierarchical
pred_counts <- matrix(NA_integer_, nrow = S_draws, ncol = N_counties)
for (s in seq_len(S_draws)) {
  probs <- theta_draws_matrix[s, ]
  pred_counts[s, ] <- rbinom(N_counties, size = df$n, prob = probs)
}
pred_mean_counts <- colMeans(pred_counts)
pred_CI_counts <- t(apply(pred_counts, 2, quantile, probs = c(0.025, 0.975)))


# WAIC calculations

# hierarchical
safe_theta <- theta_draws_matrix
safe_theta[safe_theta <= 0] <- .Machine$double.xmin
safe_theta[safe_theta >= 1] <- 1 - .Machine$double.xmin
log_lik_hier <- matrix(NA_real_, nrow = S_draws, ncol = N_counties)
for (s in seq_len(S_draws)) {
  log_lik_hier[s, ] <- dbinom(df$y, size = df$n, prob = safe_theta[s, ], log = TRUE)
}
waic_hier <- compute_waic(log_lik_hier)

# independent
ind_safe <- ind_theta_samples
ind_safe[ind_safe <= 0] <- .Machine$double.xmin
ind_safe[ind_safe >= 1] <- 1 - .Machine$double.xmin
log_lik_ind <- matrix(NA_real_, nrow = S_draws, ncol = N_counties)
for (s in seq_len(S_draws)) log_lik_ind[s, ] <- dbinom(df$y, size = df$n, prob = ind_safe[s, ], log = TRUE)
waic_ind <- compute_waic(log_lik_ind)

cat("WAIC independent:", round(waic_ind$waic,2), " pWAIC:", round(waic_ind$p_waic,2), "\n")
cat("WAIC hierarchical:", round(waic_hier$waic,2), " pWAIC:", round(waic_hier$p_waic,2), "\n")

# assembly
results <- df %>%
  mutate(ind_post_mean = ind_post_mean,
         ind_CI_low = ind_post_CI[,1],
         ind_CI_high = ind_post_CI[,2],
         hier_post_mean = theta_mean,
         hier_CI_low = theta_CI[,1],
         hier_CI_high = theta_CI[,2],
         pred_mean = pred_mean_counts,
         pred_PI_low = pred_CI_counts[,1],
         pred_PI_high = pred_CI_counts[,2],
         obs_rate = y / n,
         label = paste0(county, " | ", FIPS)
  )

write_csv(results, output_summary_csv)
saveRDS(list(theta_draws = theta_draws_matrix,
             ind_theta_samples = ind_theta_samples,
             pred_counts = pred_counts,
             pred_counts_ind = pred_counts_ind,
             state_results = state_results,
             waic_ind = waic_ind, waic_hier = waic_hier,
             df = df, results = results),
        output_samples_rds)
# tells me where
cat("Saved outputs:", output_summary_csv, "and", output_samples_rds, "\n")

# plots

# order labels by observed rate (small to large)
results <- results %>%
  mutate(label = factor(label, levels = results$label[order(results$obs_rate, decreasing = FALSE)]))

# WAIC comparison
waic_df <- tibble(
  model = c("Independent", "Hierarchical"),
  WAIC = c(waic_ind$waic, waic_hier$waic),
  pWAIC = c(waic_ind$p_waic, waic_hier$p_waic)
)
p_waic <- ggplot(waic_df, aes(x = model, y = WAIC)) +
  geom_col() +
  geom_text(aes(label = round(WAIC,1)), vjust = -0.5) +
  labs(title = "WAIC comparison", y = "WAIC (lower is better)") +
  theme_minimal()
ggsave(paste0(plot_prefix, "waic_comparison.pdf"), p_waic, width = 6, height = 4, dpi = dpi)


# predictive histogram for first row
sel <- 1
if (sel > ncol(pred_counts)) sel <- 1
df_sel <- data.frame(pred = pred_counts[, sel])
p_hist <- ggplot(df_sel, aes(x = pred)) +
  geom_histogram(bins = 30) +
  geom_vline(xintercept = df$y[sel], linetype = "dashed") +
  labs(title = paste0("Posterior predictive draws for ", df$county[sel], " | ", df$FIPS[sel]),
       x = "Predicted count", y = "Frequency") +
  theme_minimal()
ggsave(paste0(plot_prefix, "predictive_hist_selected.pdf"), p_hist, width = 6, height = 4, dpi = dpi)

cat("Plots saved. Done.\n")

# makes a 6-county shrinkage snapshot that hits interpretability points
library(dplyr); library(ggplot2)
set.seed(2025)

snap <- results %>%
  mutate(obs_rate = y / n,
         abs_diff = abs(obs_rate - hier_post_mean),
         ind_CI_width = ind_CI_high - ind_CI_low,
         hier_CI_width = hier_CI_high - hier_CI_low)

# specifies quantile cuts
q_n <- quantile(snap$n, probs = c(0.2, 0.8), na.rm = TRUE)
n_small_cut <- q_n[1]; n_large_cut <- q_n[2]

# picks one from each bucket + top & bottom
small_high <- snap %>% filter(n <= n_small_cut) %>% arrange(desc(abs_diff)) %>% slice_head(n = 1)
large_low <- snap %>% filter(n >= n_large_cut) %>% arrange(abs_diff) %>% slice_head(n = 1)
moderate <- snap %>% filter(n > n_small_cut & n < n_large_cut) %>% arrange(desc(abs_diff)) %>% slice_head(n = 1)
top_obs <- snap %>% arrange(desc(obs_rate)) %>% slice_head(n = 1)
bottom_obs <- snap %>% arrange(obs_rate) %>% slice_head(n = 1)

selected6 <- bind_rows(small_high, large_low, moderate, top_obs, bottom_obs) %>%
  distinct(FIPS, .keep_all = TRUE)

# safety -- if returned fewer than 6 unique, fill with next most informative by abs_diff
if (nrow(selected6) < 6) {
  needed <- 6 - nrow(selected6)
  filler <- snap %>% filter(!FIPS %in% selected6$FIPS) %>% arrange(desc(abs_diff)) %>% slice_head(n = needed)
  selected6 <- bind_rows(selected6, filler)
}

# final order low to high observed rate
selected6 <- selected6 %>% mutate(label = paste0(county, " (", state, " | ", FIPS, ")")) %>% arrange(obs_rate)
selected6$label <- factor(selected6$label, levels = selected6$label)

p6 <- ggplot(selected6, aes(x = label)) +
  geom_point(aes(y = obs_rate), size = 3) +
  geom_errorbar(aes(ymin = ind_CI_low, ymax = ind_CI_high), width = 0.04, position = position_nudge(x = -0.12), color = "gray40") +
  geom_point(aes(y = ind_post_mean), position = position_nudge(x = -0.12), shape = 1, color = "gray40") +
  geom_errorbar(aes(ymin = hier_CI_low, ymax = hier_CI_high), width = 0.04, position = position_nudge(x = 0.12), color = "steelblue4") +
  geom_point(aes(y = hier_post_mean), position = position_nudge(x = 0.12), shape = 17, color = "steelblue4") +
  coord_flip() +
  labs(title = "Compact shrinkage snapshot (6 counties)",
       subtitle = "Selected to show small/large/moderate, top/bottom & shrinkage behavior",
       x = "", y = "Vaccination proportion") +
  theme_minimal(base_size = 12)

ggsave("shrinkage_snapshot_6_counties.pdf", p6, width = 8.5, height = 5, dpi = 300)
print(p6)
cat("Saved shrinkage_snapshot_6_counties.pdf\n")

# interpretive notes (printed)
cat("- Small counties with large differences between observed and hierarchical estimates indicate noisy observed rates; hierarchical estimates shrink toward the state mean.\n")
cat("- Large counties with small differences show counties where observed rates are reliable and the model performs little shrinkage.\n")
cat("- Top/bottom observed counties reveal persistent high/low performance; check whether hierarchical estimates remain extreme (real signal) or shrink (likely noise).\n")
cat("- Compare CI widths: hierarchicalCI << independentCI indicates the model reduces uncertainty by borrowing strength; similar widths suggest strong data for that county.\n")

# model efficacy indicated by shrinkage focused neighborhood zoom
snap <- results %>%
  mutate(
    obs_rate = y / n,
    abs_shrink = abs(obs_rate - hier_post_mean),
    ind_CI_width = ind_CI_high - ind_CI_low,
    hier_CI_width = hier_CI_high - hier_CI_low,
    CI_improvement = ind_CI_width - hier_CI_width,   # how much the CI shrank
    shrink_ratio = hier_CI_width / ind_CI_width
  )

# pick best example county based on combined shrinkage and CI improvement
snap <- snap %>%
  mutate(model_efficacy_score = scale(abs_shrink) + scale(CI_improvement))

center <- snap %>% arrange(desc(model_efficacy_score)) %>% slice(1)
center_FIPS <- center$FIPS
cat("Center county:", center$county, center$state, "(FIPS:", center_FIPS, ")\n")

# ---- 2) Choose neighbors emphasizing contrast ----
# Priority:
#   - High shrinkage
#   - Low shrinkage (model barely changed them)
#   - High CI reduction
#   - Low CI reduction
#   - Nearest in observed rate (for context)
k <- 8   # number of counties including center

# top shrinkers (besides center)
high_shrink_neighbors <- snap %>%
  filter(FIPS != center_FIPS) %>%
  arrange(desc(abs_shrink)) %>%
  slice_head(n = 3)

# low shrink neighbors (large n, stable counts)
low_shrink_neighbors <- snap %>%
  filter(FIPS != center_FIPS) %>%
  arrange(abs_shrink) %>%
  slice_head(n = 2)

# neighbors close in observed rate
nearest_neighbors <- snap %>%
  filter(FIPS != center_FIPS) %>%
  arrange(abs(obs_rate - center$obs_rate)) %>%
  slice_head(n = 3)

# combine 
selected_eff <- bind_rows(center, high_shrink_neighbors,
                          low_shrink_neighbors, nearest_neighbors) %>%
  distinct(FIPS, .keep_all = TRUE) %>%
  slice_head(n = k)

# order for plot
selected_eff <- selected_eff %>%
  mutate(label = paste0(county, " (", state, " | ", FIPS, ")")) %>%
  arrange(obs_rate)
selected_eff$label <- factor(selected_eff$label, levels = selected_eff$label)

# plot itself
p_eff <- ggplot(selected_eff, aes(x = label)) +
  # observed
  geom_point(aes(y = obs_rate), size = 3, color = "black") +
  
  # independent CI + mean (left)
  geom_errorbar(aes(ymin = ind_CI_low, ymax = ind_CI_high),
                width = 0.05, size = 0.8,
                position = position_nudge(x = -0.2), color = "gray40") +
  geom_point(aes(y = ind_post_mean),
             position = position_nudge(x = -0.2), size = 2,
             shape = 1, color = "gray40") +
  
  # hierarchical CI + mean (right)
  geom_errorbar(aes(ymin = hier_CI_low, ymax = hier_CI_high),
                width = 0.05, size = 0.8,
                position = position_nudge(x = 0.2), color = "steelblue4") +
  geom_point(aes(y = hier_post_mean),
             position = position_nudge(x = 0.2), size = 3, shape = 17,
             color = "steelblue4") +
  
  coord_flip() +
  labs(title = "Model Efficacy Zoom: How Hierarchical Shrinkage Improves Estimates",
       subtitle = paste0("Center county: ", center$county, " (", center$state, ") — largest shrinkage + CI improvement"),
       x = "", y = "Vaccination proportion") +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(size = 10),
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10))

ggsave("shrinkage_model_efficacy_zoom.pdf", p_eff, width = 10, height = 6, dpi = 300)
print(p_eff)

cat("Saved shrinkage_model_efficacy_zoom.pdf\n")

# ensure predictive rates exist for next portion
results <- results %>%
  mutate(pred_mean_rate = as.numeric(pred_mean) / n,
         pred_low_rate  = as.numeric(pred_PI_low) / n,
         pred_high_rate = as.numeric(pred_PI_high) / n,
         obs_rate = y / n)

# compute whether observed y inside predictive 95% interval (counts)
results <- results %>%
  mutate(obs_in_pred = (y >= pred_PI_low & y <= pred_PI_high))
# empirical predictive coverage (proportion inside 95% PI)
empirical_coverage <- mean(results$obs_in_pred)
cat("Empirical predictive 95% coverage:", empirical_coverage, "\n")


snap <- results %>%
  dplyr::mutate(abs_diff_pred = abs(obs_rate - pred_mean_rate),
                pred_CI_width = pred_high_rate - pred_low_rate)

# quantiles to pick small vs large
q_n <- quantile(snap$n, probs = c(0.2, 0.8), na.rm = TRUE)
n_small_cut <- q_n[1]; n_large_cut <- q_n[2]

small <- snap %>% filter(n <= n_small_cut) %>% arrange(desc(abs_diff_pred)) %>% slice_head(n = 1)
large <- snap %>% filter(n >= n_large_cut) %>% arrange(abs_diff_pred) %>% slice_head(n = 1)
moderate <- snap %>% filter(n > n_small_cut & n < n_large_cut) %>% arrange(desc(abs_diff_pred)) %>% slice_head(n = 1)
top_pred <- snap %>% arrange(desc(pred_mean_rate)) %>% slice_head(n = 1)
bottom_pred <- snap %>% arrange(pred_mean_rate) %>% slice_head(n = 1)

selected6_pred <- bind_rows(small, large, moderate, top_pred, bottom_pred) %>%
  distinct(FIPS, .keep_all = TRUE)

# if < 6, fill by largest pred CI width
if (nrow(selected6_pred) < 6) {
  need <- 6 - nrow(selected6_pred)
  filler <- snap %>% filter(!FIPS %in% selected6_pred$FIPS) %>% arrange(desc(pred_CI_width)) %>% slice_head(n = need)
  selected6_pred <- bind_rows(selected6_pred, filler)
}

selected6_pred <- selected6_pred %>%
  dplyr::mutate(label = paste0(county, " (", state, " | ", FIPS, ")")) %>%
  arrange(pred_mean_rate)
selected6_pred$label <- factor(selected6_pred$label, levels = selected6_pred$label)

p6_pred <- ggplot(selected6_pred, aes(x = label)) +
  # predictive interval
  geom_linerange(aes(ymin = pred_low_rate, ymax = pred_high_rate), size = 1.0, color = "steelblue4") +
  geom_point(aes(y = pred_mean_rate), shape = 3, size = 3, color = "steelblue4") +
  # observed
  geom_point(aes(y = obs_rate, color = obs_in_pred), size = 3) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "red"), guide = FALSE) +
  coord_flip() +
  labs(title = "Compact predictive-interval snapshot (6 counties)",
       subtitle = "Cross = predictive mean; line = 95% posterior predictive interval; dot = observed (red = outside interval)",
       x = "", y = "Vaccination proportion") +
  theme_minimal(base_size = 12)

ggsave("predictive_snapshot_6_counties.pdf", p6_pred, width = 9, height = 5, dpi = 300)
print(p6_pred)
cat("Saved predictive_snapshot_6_counties.pdf\n")

# modelefficacy predictive zoom neighborhood

library(dplyr)
library(ggplot2)

set.seed(2025)

snap <- results %>%
  dplyr::mutate(
    obs_rate = y / n,
    pred_mean_rate = pred_mean / n,
    pred_low_rate  = pred_PI_low / n,
    pred_high_rate = pred_PI_high / n,
    abs_pred_diff = abs(obs_rate - pred_mean_rate),
    pred_CI_width = pred_high_rate - pred_low_rate   # predictive interval width
  )

# center county (most model relevant -- has the biggest abs difference btwn observed and predictive mean)
center <- snap %>% arrange(desc(abs_pred_diff)) %>% slice(1)
center_FIPS <- center$FIPS
cat("Center county for predictive zoom:", center$county, center$state, "\n")

# neighbors based on observed-rate proximity 
k <- 8   # total counties including center
neighbors <- snap %>%
  filter(FIPS != center_FIPS) %>%
  arrange(abs(obs_rate - center$obs_rate)) %>%
  slice_head(n = k)

selected_neigh_pred <- bind_rows(center, neighbors) %>%
  distinct(FIPS, .keep_all = TRUE) %>% 
  slice_head(n = k)

# order
selected_neigh_pred <- selected_neigh_pred %>%
  dplyr::mutate(label = paste0(county, " (", state, " | ", FIPS, ")")) %>%
  arrange(pred_mean_rate)
selected_neigh_pred$label <- factor(selected_neigh_pred$label, 
                                    levels = selected_neigh_pred$label)

# plotting pred intervals
p_neigh_pred <- ggplot(selected_neigh_pred, aes(x = label)) +
  geom_linerange(aes(ymin = pred_low_rate, ymax = pred_high_rate), 
                 size = 1.0, color = "steelblue4") +
  geom_point(aes(y = pred_mean_rate), shape = 3, size = 3, color = "steelblue4") +
  geom_point(aes(y = obs_rate, color = (obs_rate < pred_low_rate | obs_rate > pred_high_rate)),
             size = 3) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"), guide = FALSE) +
  coord_flip() +
  labs(
    title = paste0("Predictive intervals — Neighborhood Zoom: ", 
                   center$county, " (", center$state, ")"),
    subtitle = "Line = 95% predictive interval; cross = predictive mean; dot = observed (red = observed outside interval)",
    x = "", y = "Vaccination proportion"
  ) +
  theme_minimal(base_size = 12)

ggsave("predictive_model_efficacy_zoom.pdf", p_neigh_pred, 
       width = 10, height = 6, dpi = 300)
print(p_neigh_pred)

cat("Saved predictive_model_efficacy_zoom.pdf\n")





