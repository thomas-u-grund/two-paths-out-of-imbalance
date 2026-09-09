# ------------------------------------------------------------------
# 20: Endogeneity/identification checks. Structural covariates
#     (embeddedness, actor constraint) are
#     computed from the same yearly network that generates the outcome,
#     so they cannot be treated as exogenously manipulated. Two checks:
#
#     (1) Lagged-covariate robustness: re-estimate the primary hazard
#         models using each triad's structural covariates from ONE YEAR
#         BEFORE the (already-predetermined) row year used in the
#         primary specification -- i.e., embeddedness at t-1 predicting
#         the same t -> t+1 transition the primary model uses embeddedness
#         at t to predict. If the effect survives a second year of lag,
#         that is real (if partial) evidence against the estimates being
#         driven by short-run reverse causation/simultaneity.
#
#     (2) Cinelli-Hazlett (2020) omitted-variable sensitivity bounds:
#         computed on a linear-probability-model (LPM) approximation of
#         the primary specification, since the closed-form Cinelli-
#         Hazlett bounds are derived for OLS. This is a standard, widely
#         used practical workaround for applying the method to a
#         binary/logistic outcome, reported here explicitly as an
#         approximation, not an exact result for the logistic hazard
#         model itself. Reports the robustness value (RV) -- the minimum
#         strength (partial R^2 with both treatment and outcome) an
#         unobserved confounder would need to fully explain away the
#         embeddedness effect -- benchmarked against the strongest
#         observed covariate (capability) as an intuitive yardstick.
# ------------------------------------------------------------------
suppressMessages({ library(dplyr); library(sandwich); library(lmtest); library(sensemakr) })

# Run from the replication pack root, e.g.:
#   cd replication_pack && Rscript scripts/20_endogeneity_checks.R
# (relative paths below assume that working directory; no setwd() needed)

py <- readRDS("data/person_year_risk_table_with_duration.rds") %>%
  mutate(dur_bin = cut(duration, breaks = c(0, 1, 2, 3, 5, 10, Inf),
                        labels = c("yr1", "yr2", "yr3", "yr4-5", "yr6-10", "yr11+")))

# ------------------------------------------------------------------
# Check 1: lagged covariates (t-1)
# ------------------------------------------------------------------
tr <- readRDS("data/triads_all_years.rds") %>%
  select(year, triad_id, tie_btw_mean, tie_emb_mean, actor_load_mean, triad_log_cinc_mean) %>%
  rename(lag_tie_btw_mean = tie_btw_mean, lag_tie_emb_mean = tie_emb_mean,
         lag_actor_load_mean = actor_load_mean, lag_triad_log_cinc_mean = triad_log_cinc_mean)

py_lag <- py %>%
  mutate(lag_year = year - 1) %>%
  left_join(tr, by = c("triad_id" = "triad_id", "lag_year" = "year"))

n_missing_lag <- sum(is.na(py_lag$lag_tie_emb_mean))
cat("Rows lacking a t-1 covariate value (triad not a closed triad the prior year, dropped):",
    n_missing_lag, "of", nrow(py_lag), "\n")

py_lag <- py_lag %>% filter(!is.na(lag_tie_emb_mean)) %>%
  mutate(z_tie_btw_lag = as.numeric(scale(lag_tie_btw_mean)),
         z_tie_emb_lag = as.numeric(scale(lag_tie_emb_mean)),
         z_actor_load_lag = as.numeric(scale(lag_actor_load_mean)),
         z_cinc_lag = as.numeric(scale(lag_triad_log_cinc_mean)))

cat("Lagged-covariate sample:", nrow(py_lag), "triad-years\n")

m_dis_lag <- glm(event_dissolve ~ z_tie_btw_lag + z_tie_emb_lag + z_actor_load_lag + z_cinc_lag + era + dur_bin,
                  data = py_lag, family = binomial())
m_flip_lag <- glm(event_flip ~ z_tie_btw_lag + z_tie_emb_lag + z_actor_load_lag + z_cinc_lag + era + dur_bin,
                   data = py_lag, family = binomial())
ct_dis_lag <- coeftest(m_dis_lag, vcov = vcovCL(m_dis_lag, cluster = py_lag$triad_id))
ct_flip_lag <- coeftest(m_flip_lag, vcov = vcovCL(m_flip_lag, cluster = py_lag$triad_id))

# same-sample contemporaneous (t) comparison, so the comparison isn't
# confounded by the change in sample from dropping rows with no t-1 value
m_dis_t <- glm(event_dissolve ~ z_tie_btw + z_tie_emb + z_actor_load + z_cinc + era + dur_bin,
                data = py_lag, family = binomial())
m_flip_t <- glm(event_flip ~ z_tie_btw + z_tie_emb + z_actor_load + z_cinc + era + dur_bin,
                 data = py_lag, family = binomial())
ct_dis_t <- coeftest(m_dis_t, vcov = vcovCL(m_dis_t, cluster = py_lag$triad_id))
ct_flip_t <- coeftest(m_flip_t, vcov = vcovCL(m_flip_t, cluster = py_lag$triad_id))

cat("\n=== DISSOLUTION: contemporaneous (t) vs. lagged (t-1) embeddedness, same sample ===\n")
cat("t:   b=", coef(m_dis_t)["z_tie_emb"], " SE=", ct_dis_t["z_tie_emb","Std. Error"], " p=", ct_dis_t["z_tie_emb","Pr(>|z|)"], "\n")
cat("t-1: b=", coef(m_dis_lag)["z_tie_emb_lag"], " SE=", ct_dis_lag["z_tie_emb_lag","Std. Error"], " p=", ct_dis_lag["z_tie_emb_lag","Pr(>|z|)"], "\n")
cat("\n=== REALIGNMENT: contemporaneous (t) vs. lagged (t-1) embeddedness, same sample ===\n")
cat("t:   b=", coef(m_flip_t)["z_tie_emb"], " SE=", ct_flip_t["z_tie_emb","Std. Error"], " p=", ct_flip_t["z_tie_emb","Pr(>|z|)"], "\n")
cat("t-1: b=", coef(m_flip_lag)["z_tie_emb_lag"], " SE=", ct_flip_lag["z_tie_emb_lag","Std. Error"], " p=", ct_flip_lag["z_tie_emb_lag","Pr(>|z|)"], "\n")

# ------------------------------------------------------------------
# Check 2: Cinelli-Hazlett sensitivity bounds (LPM approximation)
# ------------------------------------------------------------------
m_dis_lpm <- lm(event_dissolve ~ z_tie_btw + z_tie_emb + z_actor_load + z_cinc + era + dur_bin, data = py)
m_flip_lpm <- lm(event_flip ~ z_tie_btw + z_tie_emb + z_actor_load + z_cinc + era + dur_bin, data = py)

sens_dis <- sensemakr(model = m_dis_lpm, treatment = "z_tie_emb",
                       benchmark_covariates = "z_cinc", kd = 1, ky = 1, q = 1, alpha = 0.05)
sens_flip <- sensemakr(model = m_flip_lpm, treatment = "z_tie_emb",
                        benchmark_covariates = "z_cinc", kd = 1, ky = 1, q = 1, alpha = 0.05)

cat("\n\n=== Cinelli-Hazlett sensitivity, DISSOLUTION (LPM approximation) ===\n")
print(summary(sens_dis))
cat("\n=== Cinelli-Hazlett sensitivity, REALIGNMENT (LPM approximation) ===\n")
print(summary(sens_flip))

sink("results/endogeneity_checks_output.txt")
cat("=== Lagged-covariate check ===\n")
cat("Rows lacking a t-1 value:", n_missing_lag, "of", nrow(py) + 0, "\n")
cat("Sample used:", nrow(py_lag), "triad-years\n\n")
cat("--- DISSOLUTION ---\n")
print(ct_dis_t); cat("\n"); print(ct_dis_lag)
cat("\n--- REALIGNMENT ---\n")
print(ct_flip_t); cat("\n"); print(ct_flip_lag)

cat("\n\n=== Cinelli-Hazlett sensitivity bounds (LPM approximation) ===\n")
cat("\n--- DISSOLUTION ---\n"); print(summary(sens_dis))
cat("\n--- REALIGNMENT ---\n"); print(summary(sens_flip))
sink()

saveRDS(list(m_dis_lag = m_dis_lag, m_flip_lag = m_flip_lag, ct_dis_lag = ct_dis_lag, ct_flip_lag = ct_flip_lag,
             ct_dis_t = ct_dis_t, ct_flip_t = ct_flip_t, n_missing_lag = n_missing_lag, n_lag_sample = nrow(py_lag),
             sens_dis = sens_dis, sens_flip = sens_flip),
        "results/endogeneity_checks.rds")
cat("\n20_endogeneity_checks.R complete. See results/endogeneity_checks_output.txt\n")
