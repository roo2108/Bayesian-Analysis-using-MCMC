# Bayesian MCMC for Gamma severity model (Metropolis-Hastings)
library(coda)
library(ggplot2)
library(gridExtra)

# 1) Data generation (use exact code from assignment)
seed <- "Enter Nummber Here"
set.seed(seed)
n_claims <- 200
claims_data_severity <- data.frame(
     claim_id = 1:n_claims,
     severity = rgamma(n_claims, shape = 3, scale = 2) # mean ~ 6, var ~ 12
)
hist(claims_data_severity$severity, xlab = "Claim Severity", main = "Histogram of Claim Severity", col = "white", border = "black", ylim = c(0, 60), xlim = c(0, 25))
head(claims_data_severity$severity)
summary(claims_data_severity$severity)
x_all <- claims_data_severity$severity

# Quick sample mean (to be used for credibility comparison)
sample_mean_all <- mean(x_all)
sample_mean_all

# 2) Log posterior function (we parametrize Gamma as shape = alpha, scale = theta)
log_posterior_transformed <- function(log_alpha, log_theta, data) { # run MH proposing on log(alpha) and log (theta)
     alpha <- exp(log_alpha)
     theta <- exp(log_theta)
     ll <- sum(dgamma(data, shape = alpha, scale = theta, log = TRUE)) # log-likelihood
     lp_alpha <- dgamma(alpha, shape = 2, scale = 1, log = TRUE) # log-prior: Gamma(2,1) on alpha and theta with scale parametrisation
     lp_theta <- dgamma(theta, shape = 2, scale = 1, log = TRUE)
     log_jacobian <- log(alpha) + log(theta) # Jacobian terms for transformation of alpha and theta
     return(ll + lp_alpha + lp_theta + log_jacobian) # total log posterior
}

# 3) Metropolis-Hastings sampler (random-walk on log-scale)
metropolis_mh <- function(data, n_iter = 50000, burn_in = 10000,
                          prop_sd_log_alpha = 0.02, prop_sd_log_theta = 0.02,
                          init_alpha = 3, init_theta = 2) {
     # initialize
     samples <- matrix(NA, nrow = n_iter, ncol = 2)
     colnames(samples) <- c("alpha", "theta")
     log_alpha_cur <- log(init_alpha)
     log_theta_cur <- log(init_theta)
     post_cur <- log_posterior_transformed(log_alpha_cur, log_theta_cur, data)
     accept_count <- 0

     for (i in seq_len(n_iter)) {
          log_alpha_prop <- rnorm(1, mean = log_alpha_cur, sd = prop_sd_log_alpha) # Proposing new log-parameters
          log_theta_prop <- rnorm(1, mean = log_theta_cur, sd = prop_sd_log_theta)
          post_prop <- log_posterior_transformed(log_alpha_prop, log_theta_prop, data)

          # acceptance probability (log scale)
          log_alpha_ratio <- post_prop - post_cur
          if (log(runif(1)) < log_alpha_ratio) {
               # accept
               log_alpha_cur <- log_alpha_prop
               log_theta_cur <- log_theta_prop
               post_cur <- post_prop
               accept_count <- accept_count + 1
          }
          # store current (transformed back)
          samples[i, "alpha"] <- exp(log_alpha_cur)
          samples[i, "theta"] <- exp(log_theta_cur)
     }

     acc_rate <- accept_count / n_iter
     samples_mcmc <- mcmc(samples)
     list(
          samples = samples_mcmc, accept_rate = acc_rate,
          burn_in = burn_in, n_iter = n_iter
     )
}

# 4) Run MCMC
set.seed(1234)
mh_out <- metropolis_mh(
     data = x_all,
     n_iter = 50000,
     burn_in = 10000,
     prop_sd_log_alpha = 0.1,
     prop_sd_log_theta = 0.1,
     init_alpha = 3, init_theta = 2
)

Acceptance_rate <- mh_out$accept_rate
print(Acceptance_rate)

# Discard burn-in and get posterior samples
post_samples <- window(mh_out$samples, start = mh_out$burn_in + 1)
alpha_post <- as.numeric(post_samples[, "alpha"])
theta_post <- as.numeric(post_samples[, "theta"])
Eseverity_post <- alpha_post * theta_post # expected claim severity per posterior sample

# Prior Mean severity
prior_mean_E <- 2 * 2 # E[alpha] * E[theta] = 4
prior_mean_E

# Posterior summaries
# Posterior of alpha
mean(alpha_post)
sd(alpha_post)

# Posterior of theta
mean(theta_post)
sd(theta_post)

# Posterior severity
mean(Eseverity_post)
sd(Eseverity_post)

# 5) Credibility factor Z
# Express posterior mean as Z * sample_mean + (1-Z) * prior_mean => solve for Z
posterior_mean_E <- mean(Eseverity_post)
prior_mean_E <- 2 * 2 # E[alpha] * E[theta] = 4
Z_hat <- (posterior_mean_E - prior_mean_E) / (sample_mean_all - prior_mean_E)
cat("\nCredibility factor Z (posterior interpreted as credibility-weighted):", Z_hat, "\n")


# 6) Convergence Diagnostics: visual and Quantitative
# EFFECTIVE SAMPLE SIZE
effectiveSize(cbind(alpha_post, theta_post))

# Integrated Autocorrelation Time (IACT)
n_samples <- 40000
ESS_alpha <- effectiveSize(alpha_post)
ESS_theta <- effectiveSize(theta_post)

IACT_alpha <- n_samples / ESS_alpha
IACT_theta <- n_samples / ESS_theta

print(IACT_alpha)
print(IACT_theta)

# 95_pct credible intervals
quantile(alpha_post, probs = c(0.025, 0.5, 0.975))
quantile(theta_post, probs = c(0.025, 0.5, 0.975))
quantile(Eseverity_post, probs = c(0.025, 0.5, 0.975))

# --- Traceplots ---
par(mfrow = c(2, 1))
plot(alpha_post,
     type = "l", col = "black",
     main = "Traceplot of alpha (α)", ylab = "alpha (α)", xlab = "Number of Iterations"
)

plot(theta_post,
     type = "l", col = "black",
     main = "Traceplot of theta (θ)", ylab = "theta (θ)", xlab = "Number of Iterations"
)

# --- Kernel density plots ---
par(mfrow = c(2, 1))
plot(density(alpha_post),
     col = "black", lwd = 2,
     main = expression(paste("Posterior Density plot for ", alpha)), xlab = expression(alpha), ylab = "Density", yaxt = "n", cex.lab = 1.3
)
abline(v = mean(alpha_post), col = "black", lty = 2, lwd = 2)

plot(density(theta_post),
     col = "black", lwd = 2,
     main = expression(paste("Posterior Density plot for ", theta)), xlab = expression(theta), ylab = "Density", yaxt = "n", cex.lab = 1.3
)
abline(v = mean(theta_post), col = "black", lty = 2, lwd = 2)

# --- Ergodic plots---
# Ergodic / running mean plots
par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))

# Running mean function
running_mean <- function(x) cumsum(x) / seq_along(x)

plot(running_mean(alpha_post),
     type = "l", col = "black",
     xlab = "Number of Iterations", ylab = "Running mean",
     main = "Ergodic plot for alpha (α)"
)
abline(h = mean(alpha_post), col = "black", lty = 2)

plot(running_mean(theta_post),
     type = "l", col = "black",
     xlab = "Number of Iterations", ylab = "Running mean",
     main = "Ergodic plot for theta (θ)"
)
abline(h = mean(theta_post), col = "black", lty = 2)


# Run 3 chains with different inits
chain1 <- metropolis_mh(x_all, n_iter = 50000, burn_in = 10000, init_alpha = 1, init_theta = 1)$samples
chain2 <- metropolis_mh(x_all, n_iter = 50000, burn_in = 10000, init_alpha = 3, init_theta = 2)$samples
chain3 <- metropolis_mh(x_all, n_iter = 50000, burn_in = 10000, init_alpha = 5, init_theta = 3)$samples

# Combine as mcmc.list
chains <- mcmc.list(mcmc(chain1), mcmc(chain2), mcmc(chain3))

# Gelman-Rubin diagnostic
gelman.diag(chains)

# 7) Credibility vs sample size: subsample the simulated dataset and compute Z(n)
compute_Z_for_n <- function(full_data, n_values, mcmc_iters = 20000, burn = 4000) {
     Zs <- numeric(length(n_values))
     for (i in seq_along(n_values)) {
          n <- n_values[i]
          data_n <- sample(full_data, n) # sample without replacement (observational)
          out_n <- metropolis_mh(
               data = data_n,
               n_iter = mcmc_iters,
               burn_in = burn,
               prop_sd_log_alpha = 0.04,
               prop_sd_log_theta = 0.04,
               init_alpha = 3, init_theta = 2
          )
          post_n <- window(out_n$samples, start = out_n$burn_in + 1)
          E_n <- as.numeric(post_n[, "alpha"]) * as.numeric(post_n[, "theta"])
          post_mean_E_n <- mean(E_n)
          sample_mean_n <- mean(data_n)
          # avoid dividing by zero when sample_mean_n == prior_mean_E
          if (abs(sample_mean_n - prior_mean_E) < 1e-8) {
               Zs[i] <- NA
          } else {
               Zs[i] <- (post_mean_E_n - prior_mean_E) / (sample_mean_n - prior_mean_E)
          }
          cat("n =", n, "Z =", Zs[i], "\n")
     }
     data.frame(n = n_values, Z = Zs)
}

# choose a grid of sample sizes to explore
n_grid <- c(5, 20, 40, 80, 120, 160, 200)
set.seed(2025)
z_vs_n <- compute_Z_for_n(x_all, n_grid, mcmc_iters = 15000, burn = 3000)

# Plot Z vs n
par(mfrow = c(1, 1))
plot(z_vs_n$n, z_vs_n$Z,
     type = "b", pch = 16,
     xlab = "Sample size (n)", ylab = "Credibility factor Z",
     main = "Credibility factor vs Sample Size", ylim = c(0, 1.5)
)
abline(h = 1, lty = 2, col = "black") # reference line for full credibility

# Monte Carlo Standard Error (MCSE)
ess <- effectiveSize(Eseverity_post)
posterior_sd <- sd(Eseverity_post)
mcse <- posterior_sd / sqrt(ess)
print(mcse)
