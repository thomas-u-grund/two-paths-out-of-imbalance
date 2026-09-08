# ------------------------------------------------------------------
# 23: Two follow-up items from a second round of critical review of
#     script 22's diagnostics:
#
#  (A) The same-sign comparison (75.2% of same-sign tie-pairs go to the
#      higher-unbalanced-load tie, n=137) is now the paper's primary
#      evidential claim for H2, but had not itself been checked for the
#      dyad-year dependence documented for the pooled 96.6% figure.
#      Count unique dyad-year transitions underlying the 137 pairs and
#      cluster-bootstrap the same-sign hit rate at that level.
#
#      Also: the 80.1% "dyad-year-level collapse" diagnostic from
#      script 22 (changed dyad had max load in ANY of several triads it
#      resolved) was compared against an implicit 33.3% baseline in the
#      write-up; that comparison is invalid once "any of several" no
#      longer has a clean 1/3 null (the null probability depends on how
#      many triads each transition touches). Recomputed here as a
#      transition-size-weighted null instead.
#
#  (B) A multinomial discrete-time robustness model (persist / dissolve
#      / realign, one joint model) as an appendix-only check against
#      the primary two-separate-cause-specific-logits specification.
# ------------------------------------------------------------------
suppressMessages({ library(dplyr); library(tidyr); library(nnet) })

# Run from the replication pack root, e.g.:
#   cd replication_pack && Rscript scripts/23_h2_dependence_and_multinomial.R
# (relative paths below assume that working directory; no setwd() needed)

# ============================================================
# (A) DEPENDENCE-CORRECT THE SAME-SIGN TEST
# ============================================================
cat("############################################################\n")
cat("(A) DEPENDENCE CORRECTION FOR THE SAME-SIGN TEST\n")
cat("############################################################\n\n")

dep <- readRDS("results/h2_sign_and_dependence.rds")
single <- readRDS("results/h2_unbalanced_load.rds")$single

# Rebuild the exact same-sign-pair long frame as script 22, but keep
# the triad_id/t1 (-> dyad-year transition id) attached to each pair.
long <- single %>%
  mutate(event_id = row_number()) %>%
  select(event_id, triad_id, t1, ul12, ul13, ul23, sign12_t1, sign13_t1, sign23_t1,
         changed12, changed13, changed23, node1, node2, node3) %>%
  pivot_longer(cols = c(ul12, ul13, ul23), names_to = "tie", values_to = "unbal_load") %>%
  mutate(tie = sub("ul", "", tie)) %>%
  left_join(
    single %>% mutate(event_id = row_number()) %>%
      select(event_id, sign12_t1, sign13_t1, sign23_t1) %>%
      pivot_longer(cols = c(sign12_t1, sign13_t1, sign23_t1), names_to = "tie2", values_to = "sign") %>%
      mutate(tie2 = sub("sign", "", tie2), tie2 = sub("_t1", "", tie2)) %>%
      select(event_id, tie2, sign),
    by = c("event_id", "tie" = "tie2")
  ) %>%
  left_join(
    single %>% mutate(event_id = row_number()) %>%
      select(event_id, changed12, changed13, changed23) %>%
      pivot_longer(cols = c(changed12, changed13, changed23), names_to = "tie3", values_to = "changed") %>%
      mutate(tie3 = sub("changed", "", tie3)) %>%
      select(event_id, tie3, changed),
    by = c("event_id", "tie" = "tie3")
  ) %>%
  mutate(
    tie_n1 = case_when(tie == "12" ~ node1, tie == "13" ~ node1, tie == "23" ~ node2),
    tie_n2 = case_when(tie == "12" ~ node2, tie == "13" ~ node3, tie == "23" ~ node3),
    dyad_lo = pmin(tie_n1, tie_n2), dyad_hi = pmax(tie_n1, tie_n2)
  )

same_sign_pairs <- long %>%
  group_by(event_id) %>%
  filter(n() == 3) %>%
  group_modify(~ {
    df <- .x
    pairs <- combn(1:3, 2, simplify = FALSE)
    purrr::map_dfr(pairs, function(idx) {
      a <- df[idx[1], ]; b <- df[idx[2], ]
      if (a$sign == b$sign && a$changed != b$changed && a$unbal_load != b$unbal_load) {
        # transition_id = the CHANGED tie's own dyad-year (the political
        # event this comparison is actually about)
        changed_row <- if (a$changed) a else b
        tibble(
          transition_id = paste(changed_row$dyad_lo, changed_row$dyad_hi, changed_row$t1, sep = "_"),
          higher_load_changed = (a$unbal_load > b$unbal_load) == a$changed
        )
      } else {
        tibble(transition_id = NA_character_, higher_load_changed = NA)
      }
    })
  }) %>%
  ungroup() %>%
  filter(!is.na(higher_load_changed))

n_unique_transitions_137 <- n_distinct(same_sign_pairs$transition_id)
cat("The 137 same-sign comparable pairs correspond to", n_unique_transitions_137,
    "unique (changed dyad, year) transitions.\n\n")

per_transition_137 <- same_sign_pairs %>% count(transition_id, name = "n_pairs")
cat("Distribution of same-sign pairs per unique transition:\n")
print(summary(per_transition_137$n_pairs))

cat("\nObserved hit rate (same as script 22):",
    round(mean(same_sign_pairs$higher_load_changed) * 100, 1), "%\n")

set.seed(43)
B <- 5000
transitions_137 <- unique(same_sign_pairs$transition_id)
boot_hits_137 <- numeric(B)
for (b in seq_len(B)) {
  samp_trans <- sample(transitions_137, length(transitions_137), replace = TRUE)
  idx <- unlist(lapply(samp_trans, function(t) which(same_sign_pairs$transition_id == t)))
  boot_hits_137[b] <- mean(same_sign_pairs$higher_load_changed[idx])
}
cat("\nCluster bootstrap (5,000 resamples, resampling dyad-year transitions with replacement):\n")
ci_lo <- quantile(boot_hits_137, 0.025); ci_hi <- quantile(boot_hits_137, 0.975)
cat("Bootstrap 95% CI:", round(ci_lo * 100, 1), "% to", round(ci_hi * 100, 1), "%\n")
cat("(Naive exact-binomial CI treating all 137 pairs as independent: [67.1%, 82.2%])\n")

# Also fit the clogit-equivalent (paired logistic) with clustered SEs at
# the transition level, as an alternative to the bootstrap.
same_sign_pairs_id <- same_sign_pairs %>% mutate(pair_id = row_number())
m_glm <- glm(higher_load_changed ~ 1, data = same_sign_pairs_id, family = binomial())
library(sandwich); library(lmtest)
ct_clustered <- coeftest(m_glm, vcov = vcovCL(m_glm, cluster = same_sign_pairs_id$transition_id))
cat("\nIntercept-only logit, SEs clustered on dyad-year transition:\n")
print(ct_clustered)
p_hat <- plogis(coef(m_glm)[1])
se_logit <- ct_clustered[1, "Std. Error"]
ci_logit <- plogis(coef(m_glm)[1] + c(-1.96, 1.96) * se_logit)
cat("Implied hit rate:", round(p_hat * 100, 1), "%  cluster-robust 95% CI: [",
    round(ci_logit[1] * 100, 1), "%,", round(ci_logit[2] * 100, 1), "%]\n")

# ---- Fix the invalid 33.3%-baseline framing on the 80.1% diagnostic ----
cat("\n--- Correcting the dyad-year-level 80.1% diagnostic's implied baseline ---\n")
strict_b <- dep$single_b %>% filter(n_distinct_ul == 3)
dyad_level <- strict_b %>% group_by(transition_id) %>%
  summarise(is_max_ul_any = any(is_max_ul), n_triads = n(), .groups = "drop")
# Under a null where each triad independently assigns "is_max_ul" to the
# changed tie with probability 1/3, the chance that AT LEAST ONE of a
# transition's n_triads triads shows this by chance is 1-(2/3)^n_triads,
# which rises with n_triads -- so the correct baseline is transition-size
# weighted, not a flat 33.3%.
dyad_level <- dyad_level %>% mutate(chance_baseline = 1 - (2/3)^n_triads)
weighted_null <- mean(dyad_level$chance_baseline)
cat("Observed: 80.1% of transitions have max UL in at least one resolved triad.\n")
cat("Size-weighted chance baseline (not a flat 33.3%): would be",
    round(weighted_null * 100, 1), "% if unbalanced load carried no information\n")
cat("(most transitions resolve many triads, so EVEN UNDER THE NULL, 'at least one hit'\n")
cat(" is already likely by chance -- this diagnostic should be read descriptively,\n")
cat(" not as a hypothesis test against 33.3%).\n")

# ============================================================
# (B) MULTINOMIAL DISCRETE-TIME ROBUSTNESS MODEL
# ============================================================
cat("\n\n############################################################\n")
cat("(B) MULTINOMIAL DISCRETE-TIME MODEL (persist / dissolve / realign)\n")
cat("############################################################\n\n")

py <- readRDS("data/person_year_risk_table_with_duration.rds") %>%
  mutate(
    dur_bin = cut(duration, breaks = c(0, 1, 2, 3, 5, 10, Inf),
                  labels = c("yr1", "yr2", "yr3", "yr4-5", "yr6-10", "yr11+")),
    outcome3 = factor(case_when(
      event_dissolve == 1 ~ "dissolve",
      event_flip == 1 ~ "realign",
      TRUE ~ "persist"
    ), levels = c("persist", "dissolve", "realign"))
  )
cat("Outcome distribution:\n"); print(table(py$outcome3))

m_multi <- multinom(outcome3 ~ z_tie_btw + z_tie_emb + z_actor_load + z_cinc + era + dur_bin,
                     data = py, trace = FALSE, maxit = 300)
cat("\n--- Multinomial logit (reference category: persist) ---\n")
print(summary(m_multi))

co <- summary(m_multi)$coefficients
se <- summary(m_multi)$standard.errors
z <- co / se
p <- 2 * (1 - pnorm(abs(z)))
cat("\n--- z_tie_emb (triadic embeddedness), primary comparison ---\n")
cat("dissolve: b=", round(co["dissolve", "z_tie_emb"], 3),
    " SE=", round(se["dissolve", "z_tie_emb"], 3),
    " p=", signif(p["dissolve", "z_tie_emb"], 3), "\n")
cat("realign:  b=", round(co["realign", "z_tie_emb"], 3),
    " SE=", round(se["realign", "z_tie_emb"], 3),
    " p=", signif(p["realign", "z_tie_emb"], 3), "\n")
cat("\nCompare to primary two-separate-logit specification (main text Table 1):\n")
cat("dissolution b=-2.508, realignment b=+0.439\n")

sink("results/h2_dependence_and_multinomial_output.txt")
cat("=== (A) Same-sign test, dependence-corrected ===\n")
cat("137 pairs ->", n_unique_transitions_137, "unique dyad-year transitions\n")
cat("Observed hit rate: 75.2%\n")
cat("Cluster bootstrap 95% CI:", round(ci_lo*100,1), "% to", round(ci_hi*100,1), "%\n")
cat("Cluster-robust logit CI: [", round(ci_logit[1]*100,1), "%,", round(ci_logit[2]*100,1), "%]\n")
cat("\n80.1% diagnostic: size-weighted chance baseline =", round(weighted_null*100,1), "%",
    "(not 33.3%)\n\n")
cat("=== (B) Multinomial model ===\n")
cat("z_tie_emb dissolve: b=", round(co["dissolve","z_tie_emb"],3), " p=", signif(p["dissolve","z_tie_emb"],3), "\n")
cat("z_tie_emb realign:  b=", round(co["realign","z_tie_emb"],3), " p=", signif(p["realign","z_tie_emb"],3), "\n")
sink()

saveRDS(list(same_sign_pairs = same_sign_pairs, n_unique_transitions_137 = n_unique_transitions_137,
             boot_hits_137 = boot_hits_137, ci_logit = ci_logit, dyad_level = dyad_level,
             weighted_null = weighted_null, m_multi = m_multi),
        "results/h2_dependence_and_multinomial.rds")
cat("\n23_h2_dependence_and_multinomial.R complete.\n")
