# ------------------------------------------------------------------
# 04: Empirical models on the fixed 5-year resolution panel (script 03).
#     Outcome is three-way: an unbalanced triad, by t+5, either dissolved
#     (a tie disappeared), flipped (became balanced), or froze (stayed
#     closed and unbalanced) -- rather than a collapsed binary "resolved"
#     outcome, since dissolution and flipping are conceptually distinct
#     resolutions.
# ------------------------------------------------------------------
suppressMessages({
  library(dplyr); library(readr); library(broom); library(sandwich); library(lmtest)
  library(car); library(nnet)
})

# Run from the replication pack root, e.g.:
#   cd replication_pack && Rscript scripts/04_analysis.R
# (relative paths below assume that working directory; no setwd() needed)
dir.create("results", showWarnings = FALSE)

tr <- readRDS("data/triads_all_years.rds")
status_lookup <- tr %>% select(year, triad_id, balanced)

# Restrict to genuine SPELL ONSETS, not every year a triad happens to still
# be unbalanced. Without this, a triad that stays unbalanced for (say) 20
# consecutive years contributes ~20 highly non-independent "starting point"
# rows -- one per year -- almost all coded the same way, while a triad that
# resolves in 1 year contributes exactly 1. Clustering SEs by triad_id
# corrects inference for within-triad correlation but does NOT fix this:
# the point estimates themselves are still computed by treating each
# repeated row as separately informative, mechanically oversampling
# chronic/persistent cases (in this data, 476 of ~2,800 unique triads,
# each appearing >5 times, contributed 52% of all rows). Fixed by keeping
# only triad-years where the triad was NOT already unbalanced the
# previous year (i.e., it just became unbalanced, or is newly observed).
unbalanced_all <- tr %>% filter(!balanced)
prior_year_status <- status_lookup %>% mutate(year = year + 1) %>% rename(balanced_prior = balanced)
unbalanced <- unbalanced_all %>%
  left_join(prior_year_status, by = c("year", "triad_id")) %>%
  filter(is.na(balanced_prior) | balanced_prior)  # NA = triad not closed last year (new); TRUE = was balanced last year (newly unbalanced)
cat("Spell-onset filter: ", nrow(unbalanced_all), "raw unbalanced triad-years ->",
    nrow(unbalanced), "genuine onsets (", round(100*nrow(unbalanced)/nrow(unbalanced_all),1), "%)\n")
tie_trade <- read_csv("data/tie_trade.csv", show_col_types = FALSE) %>%
  mutate(ccode_low = as.character(ccode_low), ccode_high = as.character(ccode_high)) %>%
  select(year, ccode_low, ccode_high, trade_dependence)
get_trade <- function(df, a, b, newcol) {
  df %>% mutate(lo = pmin(.data[[a]], .data[[b]]), hi = pmax(.data[[a]], .data[[b]])) %>%
    left_join(tie_trade, by = c("year" = "year", "lo" = "ccode_low", "hi" = "ccode_high")) %>%
    rename(!!newcol := trade_dependence) %>% select(-lo, -hi)
}

df <- unbalanced %>%
  mutate(target_year = year + 5) %>%
  left_join(status_lookup, by = c("target_year" = "year", "triad_id" = "triad_id"), suffix = c("", "_t5")) %>%
  mutate(
    outcome = factor(case_when(
      is.na(balanced_t5) ~ "dissolved",
      balanced_t5         ~ "flipped",
      !balanced_t5         ~ "frozen"
    ), levels = c("flipped", "dissolved", "frozen")),  # "flipped" = reference (classical balance-theory outcome)
    frozen = outcome == "frozen",
    era = cut(year, breaks = c(1815, 1945, 1989, 2015), labels = c("pre-1945", "1945-1989", "post-1989"))
  ) %>%
  get_trade("node1", "node2", "trade12") %>% get_trade("node1", "node3", "trade13") %>% get_trade("node2", "node3", "trade23") %>%
  rowwise() %>% mutate(tie_trade_mean = mean(c(trade12, trade13, trade23), na.rm = TRUE)) %>% ungroup()

write_csv(df, "data/outcome_panel.csv")
saveRDS(df, "data/outcome_panel.rds")

cat("N unbalanced triad-years:", nrow(df), "\n")
print(table(df$outcome))
cat(sprintf("dissolved %.1f%% | flipped %.1f%% | frozen %.1f%%\n",
            100 * mean(df$outcome == "dissolved"), 100 * mean(df$outcome == "flipped"), 100 * mean(df$outcome == "frozen")))

zscore <- function(x) as.numeric(scale(x))
df <- df %>% mutate(
  z_actor_load = zscore(actor_load_mean),
  z_tie_btw    = zscore(tie_btw_mean),
  z_tie_emb    = zscore(tie_emb_mean),
  z_cinc       = zscore(triad_log_cinc_mean)
)

# ------------------------------------------------------------------
# PRIMARY MODEL: multinomial logit, outcome in {flipped, dissolved, frozen}
# reference category = flipped (the classical balance-theory outcome)
# ------------------------------------------------------------------
m_multi <- multinom(outcome ~ z_tie_btw + z_tie_emb + I(z_tie_emb^2) + z_actor_load + z_cinc + era,
                     data = df, trace = FALSE)
multi_sum <- summary(m_multi)
z_multi <- multi_sum$coefficients / multi_sum$standard.errors
p_multi <- 2 * (1 - pnorm(abs(z_multi)))

# ------------------------------------------------------------------
# H1: TWO-EQUATION DECOMPOSITION.
# A single collapsed "frozen vs. everything else" model with a quadratic
# tie-load term produces a spurious-looking curvilinear result that is
# actually a compositional artifact of pooling two categories (dissolved,
# flipped) that both move against frozen at different rates. The correct
# primary specification decomposes into two clean pairwise contrasts:
#   (i)  dissolved vs. flipped  -- does load prevent a tie from lapsing?
#   (ii) frozen vs. flipped     -- among triads that stay intact, does
#        load prevent the flip to balance?
# Both are estimated on their respective two-category subsamples (each
# triad-year excluded from the contrast it isn't part of), with clustered
# SEs, and a quadratic term is tested in each but not assumed.
# ------------------------------------------------------------------
df_diss_flip   <- df %>% filter(outcome %in% c("dissolved", "flipped"))
df_frozen_flip <- df %>% filter(outcome %in% c("frozen", "flipped"))

m_dissolved_lin  <- glm(I(outcome == "dissolved") ~ z_tie_btw + z_tie_emb + z_actor_load + z_cinc + era, data = df_diss_flip, family = binomial())
m_dissolved_quad <- glm(I(outcome == "dissolved") ~ z_tie_btw + z_tie_emb + I(z_tie_emb^2) + z_actor_load + z_cinc + era, data = df_diss_flip, family = binomial())
m_frozen_lin  <- glm(I(outcome == "frozen") ~ z_tie_btw + z_tie_emb + z_actor_load + z_cinc + era, data = df_frozen_flip, family = binomial())
m_frozen_quad <- glm(I(outcome == "frozen") ~ z_tie_btw + z_tie_emb + I(z_tie_emb^2) + z_actor_load + z_cinc + era, data = df_frozen_flip, family = binomial())

ct_dissolved_lin  <- coeftest(m_dissolved_lin,  vcov = vcovCL(m_dissolved_lin,  cluster = df_diss_flip$triad_id))
ct_dissolved_quad <- coeftest(m_dissolved_quad, vcov = vcovCL(m_dissolved_quad, cluster = df_diss_flip$triad_id))
ct_frozen_lin  <- coeftest(m_frozen_lin,  vcov = vcovCL(m_frozen_lin,  cluster = df_frozen_flip$triad_id))
ct_frozen_quad <- coeftest(m_frozen_quad, vcov = vcovCL(m_frozen_quad, cluster = df_frozen_flip$triad_id))

lrt_dissolved <- anova(m_dissolved_lin, m_dissolved_quad, test = "Chisq")
lrt_frozen    <- anova(m_frozen_lin, m_frozen_quad, test = "Chisq")

# retained for the descriptive Figure 2 (predicted share of "frozen" among
# ALL unbalanced triads) -- explicitly a compositional/derived quantity,
# not treated as evidence of a direct nonlinear mechanism
m_frozen_pooled_quad <- glm(frozen ~ z_tie_btw + z_tie_emb + I(z_tie_emb^2) + z_actor_load + z_cinc + era, data = df, family = binomial())

# ------------------------------------------------------------------
# H2: within triads that DO flip, does resolution concentrate on the
# lowest-load tie? (proxy: min of tie_btw/tie_emb; same logic as before,
# tested only among triads that flipped)
# ------------------------------------------------------------------
flipped_df <- df %>% mutate(min_tie_load_proxy = pmin(tie_btw_mean, tie_emb_mean, na.rm = TRUE),
                             z_min_tie = zscore(min_tie_load_proxy))
m2 <- glm(I(outcome == "flipped") ~ z_min_tie + z_actor_load + z_cinc + era, data = flipped_df, family = binomial())
ct2 <- coeftest(m2, vcov = vcovCL(m2, cluster = flipped_df$triad_id))

# ------------------------------------------------------------------
# H3: actor load and initiation — tested on the SAME two clean contrasts
# as H1, for consistency (both bivariate and conditional specifications
# reported, for transparency)
# ------------------------------------------------------------------
m3  <- glm(I(outcome == "frozen") ~ z_actor_load + z_cinc + era, data = df_frozen_flip, family = binomial())
ct3 <- coeftest(m3, vcov = vcovCL(m3, cluster = df_frozen_flip$triad_id))
m3b <- glm(I(outcome == "frozen") ~ z_actor_load, data = df_frozen_flip, family = binomial())
ct3b <- coeftest(m3b, vcov = vcovCL(m3b, cluster = df_frozen_flip$triad_id))
m3_diss  <- glm(I(outcome == "dissolved") ~ z_actor_load + z_cinc + era, data = df_diss_flip, family = binomial())
ct3_diss <- coeftest(m3_diss, vcov = vcovCL(m3_diss, cluster = df_diss_flip$triad_id))

# ------------------------------------------------------------------
# Robustness: trade-based tie load (1870-2012 subsample), frozen-vs-flipped
# ------------------------------------------------------------------
df_trade <- df_frozen_flip %>% filter(!is.na(tie_trade_mean), year >= 1870) %>% mutate(z_trade = zscore(tie_trade_mean))
m_robust <- glm(I(outcome == "frozen") ~ z_trade + z_tie_emb + z_actor_load + z_cinc + era, data = df_trade, family = binomial())
ct_robust <- coeftest(m_robust, vcov = vcovCL(m_robust, cluster = df_trade$triad_id))

# ------------------------------------------------------------------
# Save everything
# ------------------------------------------------------------------
sink("results/model_output.txt")
cat("=== Outcome distribution ===\n"); print(table(df$outcome)); cat("\n")

cat("=== Multinomial logit, descriptive only (ref = flipped) ===\n")
cat("(Reported for completeness; the two-equation decomposition below is\n")
cat(" the primary specification.)\n")
cat("\n-- Coefficients --\n"); print(round(multi_sum$coefficients, 3))
cat("\n-- p-values --\n"); print(round(p_multi, 4))
cat("\nN =", nrow(df), " AIC =", AIC(m_multi), "\n")

cat("\n\n=== PRIMARY H1a: dissolved vs flipped — linear ===\n")
print(ct_dissolved_lin)
cat("\n\n=== PRIMARY H1a: dissolved vs flipped — quadratic ===\n")
print(ct_dissolved_quad)
cat("\nLRT linear vs quadratic (dissolved vs flipped):\n"); print(lrt_dissolved)

cat("\n\n=== PRIMARY H1b: frozen vs flipped — linear ===\n")
print(ct_frozen_lin)
cat("\n\n=== PRIMARY H1b: frozen vs flipped — quadratic ===\n")
print(ct_frozen_quad)
cat("\nLRT linear vs quadratic (frozen vs flipped):\n"); print(lrt_frozen)

cat("\n\n=== [Descriptive only] Pooled frozen-vs-rest quadratic model ===\n")
cat("(used only to generate Fig. 2's descriptive P(frozen) curve; NOT\n")
cat(" interpreted as a direct nonlinear mechanism)\n")
print(summary(m_frozen_pooled_quad)$coefficients)

cat("\n\n=== H2: lowest-load tie predicts flip-not-freeze ===\n")
print(ct2)

cat("\n\n=== H3: actor load, frozen vs flipped (conditional) ===\n")
print(ct3)
cat("\n\n=== H3: actor load, frozen vs flipped (unconditional/bivariate) ===\n")
print(ct3b)
cat("\n\n=== H3 (extension): actor load, dissolved vs flipped ===\n")
print(ct3_diss)

cat("\n\n=== Robustness: trade-based tie load, 1870-2012 subsample, frozen vs flipped ===\n")
print(ct_robust)
cat("\nN =", nobs(m_robust), "\n")
sink()

saveRDS(list(m_multi = m_multi, multi_sum = multi_sum, p_multi = p_multi,
             m_dissolved_lin = m_dissolved_lin, ct_dissolved_lin = ct_dissolved_lin,
             m_dissolved_quad = m_dissolved_quad, ct_dissolved_quad = ct_dissolved_quad, lrt_dissolved = lrt_dissolved,
             m_frozen_lin = m_frozen_lin, ct_frozen_lin = ct_frozen_lin,
             m_frozen_quad = m_frozen_quad, ct_frozen_quad = ct_frozen_quad, lrt_frozen = lrt_frozen,
             m_frozen_pooled_quad = m_frozen_pooled_quad,
             m2 = m2, ct2 = ct2, m3 = m3, ct3 = ct3, m3b = m3b, ct3b = ct3b,
             m3_diss = m3_diss, ct3_diss = ct3_diss,
             m_robust = m_robust, ct_robust = ct_robust, df = df,
             df_diss_flip = df_diss_flip, df_frozen_flip = df_frozen_flip),
        "results/models.rds")

cat("\n04_analysis.R complete. See results/model_output.txt\n")
