# ------------------------------------------------------------------
# 22: Three diagnostics requested directly by a critical external
#     review of H2 (the unbalanced-load tie-choice result):
#
#  (A) Sign-composition test. Does unbalanced load predict which tie
#      changes beyond simply identifying "the unique negative tie in a
#      ++- triad"? Compare a sign-only model, load-only, and both
#      jointly; also restrict to comparisons among same-sign ties only
#      (where a sign indicator carries zero information) to see whether
#      load still discriminates.
#
#  (B) Dyad-year dependence / pseudoreplication. How many of the 5,414
#      single-tie-change realignment "events" are actually distinct
#      political transitions, versus the same underlying dyad-year
#      sign change appearing once per unbalanced triad it happens to
#      also resolve? Count unique (dyad, year) transitions; report the
#      distribution of triads-per-transition; re-run the primary H2
#      test clustering at the dyad-year-transition level.
#
#  (C) Net balance gain. The reviewer's point that a tie's B (balanced
#      overlapping triangles, which a flip would newly UNbalance) should
#      offset its U (unbalanced overlapping triangles, which a flip
#      would resolve) -- test U-B (equivalently 2U-E, E = total
#      embeddedness = U+B) as an alternative/additional predictor to
#      unbalanced load alone.
#
#     Uses the same single-tie-change realignment sample as scripts 09
#     and 12 (results/h2_unbalanced_load.rds), not a new sample.
# ------------------------------------------------------------------
suppressMessages({ library(dplyr); library(tidyr); library(survival); library(purrr) })

# Run from the replication pack root, e.g.:
#   cd replication_pack && Rscript scripts/22_h2_sign_and_dependence.R
# (relative paths below assume that working directory; no setwd() needed)

h2 <- readRDS("results/h2_unbalanced_load.rds")
single <- h2$single  # 5,414 single-tie-change realignment events

# ============================================================
# (A) SIGN-COMPOSITION TEST
# ============================================================
cat("############################################################\n")
cat("(A) SIGN-COMPOSITION TEST\n")
cat("############################################################\n\n")

# For each event, is the changed tie the negative one at t1?
single_a <- single %>%
  mutate(
    changed_sign = case_when(changed12 ~ sign12_t1, changed13 ~ sign13_t1, changed23 ~ sign23_t1),
    n_negative = (sign12_t1 == -1) + (sign13_t1 == -1) + (sign23_t1 == -1),
    changed_is_negative = changed_sign == -1
  )

cat("Starting configuration of the 5,414 single-tie-change events (by count of negative ties at t1):\n")
print(table(single_a$n_negative))
cat("\n(1 negative = the standard ++- configuration; 3 negative = ---, the all-hostile case)\n\n")

only1neg <- single_a %>% filter(n_negative == 1)
cat("Restricting to the dominant ++- configuration (n=", nrow(only1neg), "):\n", sep = "")
cat("  Q1: In ++- events, how often IS the negative tie the max-unbalanced-load tie?\n")
q1 <- mean(only1neg$is_max_ul)
cat("      ", round(q1*100, 1), "%\n", sep = "")
cat("  Q2: In ++- events, how often DOES the negative tie actually change?\n")
q2 <- mean(only1neg$changed_is_negative)
cat("      ", round(q2*100, 1), "%\n", sep = "")

# Conditional logit: sign-only, load-only, both -- same 3-alternatives-per-event
# structure as scripts 09/12.
long <- single %>%
  mutate(event_id = row_number()) %>%
  select(event_id, ul12, ul13, ul23, sign12_t1, sign13_t1, sign23_t1,
         changed12, changed13, changed23) %>%
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
  mutate(is_negative = sign == -1)

cat("\n--- Conditional logit: sign indicator ONLY ---\n")
clog_sign <- clogit(changed ~ is_negative + strata(event_id), data = long)
print(summary(clog_sign))

cat("\n--- Conditional logit: unbalanced load ONLY (same as script 12) ---\n")
clog_load <- clogit(changed ~ unbal_load + strata(event_id), data = long)
print(summary(clog_load))

cat("\n--- Conditional logit: sign + unbalanced load JOINTLY ---\n")
clog_both <- clogit(changed ~ is_negative + unbal_load + strata(event_id), data = long)
print(summary(clog_both))

cat("\n--- Concordance comparison ---\n")
cat("sign only:      ", round(summary(clog_sign)$concordance[1], 4), "\n")
cat("load only:      ", round(summary(clog_load)$concordance[1], 4), "\n")
cat("sign + load:    ", round(summary(clog_both)$concordance[1], 4), "\n")

# Q3/Q4: restrict comparisons to ties that SHARE a sign (where a sign
# indicator is uninformative by construction) -- does load still
# discriminate among ties of the same sign?
cat("\n--- Restricting to same-sign tie *pairs* only: does load still discriminate ---\n")
cat("    where a sign indicator carries ZERO information within the pair?\n")
same_sign_pairs <- long %>%
  group_by(event_id) %>%
  filter(n() == 3) %>%
  group_modify(~ {
    df <- .x
    pairs <- combn(1:3, 2, simplify = FALSE)
    map_dfr(pairs, function(idx) {
      a <- df[idx[1], ]; b <- df[idx[2], ]
      if (a$sign == b$sign && a$changed != b$changed) {
        # exactly one of the pair changed (the third tie is irrelevant here)
        tibble(higher_load_changed = if (a$unbal_load == b$unbal_load) NA
               else (a$unbal_load > b$unbal_load) == a$changed)
      } else {
        tibble(higher_load_changed = NA_real_)
      }
    })
  }) %>%
  ungroup() %>%
  filter(!is.na(higher_load_changed))

cat("N same-sign, exactly-one-changed comparable pairs (excl. load ties):", nrow(same_sign_pairs), "\n")
cat("Share where the higher-unbalanced-load tie of the pair is the one that changed:",
    round(mean(same_sign_pairs$higher_load_changed) * 100, 1), "% (null = 50%)\n")
print(binom.test(sum(same_sign_pairs$higher_load_changed), nrow(same_sign_pairs), p = 0.5))

# ============================================================
# (B) DYAD-YEAR DEPENDENCE / PSEUDOREPLICATION
# ============================================================
cat("\n\n############################################################\n")
cat("(B) DYAD-YEAR DEPENDENCE / PSEUDOREPLICATION\n")
cat("############################################################\n\n")

single_b <- single %>%
  mutate(
    changed_n1 = case_when(changed12 ~ node1, changed13 ~ node1, changed23 ~ node2),
    changed_n2 = case_when(changed12 ~ node2, changed13 ~ node3, changed23 ~ node3),
    dyad_lo = pmin(changed_n1, changed_n2), dyad_hi = pmax(changed_n1, changed_n2),
    transition_id = paste(dyad_lo, dyad_hi, t1, sep = "_")
  )

n_unique_transitions <- n_distinct(single_b$transition_id)
cat("5,414 single-tie-change realignment events arose from", n_unique_transitions,
    "unique (changed dyad, year) transitions.\n\n")

per_transition <- single_b %>% count(transition_id, name = "n_triads_resolved")
cat("Distribution of triads resolved per unique dyad-year transition:\n")
print(summary(per_transition$n_triads_resolved))
cat("\nMedian:", median(per_transition$n_triads_resolved),
    " | Max:", max(per_transition$n_triads_resolved), "\n")
cat("\nShare of transitions that resolve exactly 1 triad:",
    round(mean(per_transition$n_triads_resolved == 1) * 100, 1), "%\n")
cat("Share of transitions that resolve 10+ triads simultaneously:",
    round(mean(per_transition$n_triads_resolved >= 10) * 100, 1), "%\n")
cat("Top 10 largest transitions (dyad-year, n triads resolved):\n")
print(head(per_transition %>% arrange(desc(n_triads_resolved)), 10))

# Re-run the PRIMARY H2 binomial test (strict: all 3 unbalanced-load
# values distinct), but clustered at the dyad-year-transition level via
# a cluster bootstrap, instead of treating each triad-event as
# independent.
strict_b <- single_b %>% filter(n_distinct_ul == 3)
cat("\nPrimary H2 sample (all 3 UL values distinct), n =", nrow(strict_b),
    ", from", n_distinct(strict_b$transition_id), "unique dyad-year transitions.\n")

set.seed(42)
B <- 2000
transitions <- unique(strict_b$transition_id)
boot_hits <- numeric(B)
for (b in seq_len(B)) {
  samp_trans <- sample(transitions, length(transitions), replace = TRUE)
  # rebuild a resampled event set by pulling each sampled transition's rows
  idx <- unlist(lapply(samp_trans, function(t) which(strict_b$transition_id == t)))
  boot_hits[b] <- mean(strict_b$is_max_ul[idx])
}
cat("\nCluster bootstrap (2,000 resamples, resampling dyad-year transitions with replacement):\n")
cat("Point estimate (observed):", round(mean(strict_b$is_max_ul) * 100, 1), "%\n")
cat("Bootstrap 95% CI:", round(quantile(boot_hits, 0.025) * 100, 1), "% to",
    round(quantile(boot_hits, 0.975) * 100, 1), "%\n")
cat("(Compare to the naive binomial CI reported in the main text/appendix, which treats\n")
cat(" every triad-event as an independent draw.)\n")

# Dyad-year-level reframing: among unique dyad-year transitions, does
# the MAX unbalanced load across the triads that dyad touches predict
# whether it was in fact this dyad (vs its co-triad-members' other
# edges) that changed? -- collapse to one row per transition using the
# transition's own changed edge and its unbalanced load in whichever
# triad gave it the highest reported load, as an upper-bound-power
# dyad-level check.
dyad_level <- strict_b %>%
  group_by(transition_id) %>%
  summarise(is_max_ul_any = any(is_max_ul), n_triads = n(), .groups = "drop")
cat("\nDyad-year-level collapse (one row per unique transition, n=", nrow(dyad_level), "):\n", sep = "")
cat("Share of unique transitions where the changed dyad had the max unbalanced load in AT LEAST ONE\n")
cat("of the triads it resolved:", round(mean(dyad_level$is_max_ul_any) * 100, 1), "%\n")

# ============================================================
# (C) NET BALANCE GAIN (U - B, equivalently 2U - E)
# ============================================================
cat("\n\n############################################################\n")
cat("(C) NET BALANCE GAIN (U-B) AS AN ALTERNATIVE TO UNBALANCED LOAD ALONE\n")
cat("############################################################\n\n")

single_c <- single %>%
  mutate(
    nbg12 = 2 * ul12 - emb12,  # U - B = U - (E - U) = 2U - E
    nbg13 = 2 * ul13 - emb13,
    nbg23 = 2 * ul23 - emb23,
    changed_nbg = case_when(changed12 ~ nbg12, changed13 ~ nbg13, changed23 ~ nbg23),
    max_nbg = pmax(nbg12, nbg13, nbg23),
    n_distinct_nbg = pmap_int(list(nbg12, nbg13, nbg23), function(a, b, c) length(unique(c(a, b, c)))),
    is_max_nbg = changed_nbg == max_nbg
  )

strict_nbg <- single_c %>% filter(n_distinct_nbg == 3)
cat("Net balance gain (U-B) tie structure: all 3 distinct in", nrow(strict_nbg),
    "of", nrow(single_c), sprintf("(%.1f%%)\n", 100 * nrow(strict_nbg) / nrow(single_c)))
cat("\nPRIMARY test (all 3 U-B values distinct, n=", nrow(strict_nbg), "):\n", sep = "")
cat("Changed tie has MAX (U-B) in", round(mean(strict_nbg$is_max_nbg) * 100, 1), "% of events (null=33.3%)\n")
print(binom.test(sum(strict_nbg$is_max_nbg), nrow(strict_nbg), p = 1/3))

long_c <- single_c %>%
  mutate(event_id = row_number()) %>%
  select(event_id, nbg12, nbg13, nbg23, ul12, ul13, ul23, changed12, changed13, changed23) %>%
  pivot_longer(cols = c(nbg12, nbg13, nbg23), names_to = "tie", values_to = "net_balance_gain") %>%
  mutate(tie = sub("nbg", "", tie)) %>%
  left_join(
    single_c %>% mutate(event_id = row_number()) %>%
      select(event_id, ul12, ul13, ul23) %>%
      pivot_longer(cols = c(ul12, ul13, ul23), names_to = "tie2", values_to = "unbal_load") %>%
      mutate(tie2 = sub("ul", "", tie2)) %>% select(event_id, tie2, unbal_load),
    by = c("event_id", "tie" = "tie2")
  ) %>%
  left_join(
    single_c %>% mutate(event_id = row_number()) %>%
      select(event_id, changed12, changed13, changed23) %>%
      pivot_longer(cols = c(changed12, changed13, changed23), names_to = "tie3", values_to = "changed") %>%
      mutate(tie3 = sub("changed", "", tie3)) %>% select(event_id, tie3, changed),
    by = c("event_id", "tie" = "tie3")
  )

cat("\n--- Conditional logit: net balance gain (U-B) ONLY ---\n")
clog_nbg <- clogit(changed ~ net_balance_gain + strata(event_id), data = long_c)
print(summary(clog_nbg))

cat("\n--- Conditional logit: unbalanced load vs. net balance gain, head to head ---\n")
clog_compare <- clogit(changed ~ unbal_load + net_balance_gain + strata(event_id), data = long_c)
print(summary(clog_compare))

cat("\n--- AIC comparison (unbalanced load alone vs. net balance gain alone) ---\n")
cat("Unbalanced load only AIC:", round(AIC(clog_load), 2), "\n")
cat("Net balance gain only AIC:", round(AIC(clog_nbg), 2), "\n")

# ------------------------------------------------------------------
sink("results/h2_sign_and_dependence_output.txt")
cat("=== (A) Sign-composition ===\n")
cat("++- configuration Q1 (neg tie has max UL):", round(q1*100,1), "%\n")
cat("++- configuration Q2 (neg tie actually changes):", round(q2*100,1), "%\n")
cat("Concordance -- sign only:", round(summary(clog_sign)$concordance[1],4),
    " load only:", round(summary(clog_load)$concordance[1],4),
    " sign+load:", round(summary(clog_both)$concordance[1],4), "\n")
cat("Same-sign-pair comparison: n=", nrow(same_sign_pairs),
    ", higher-load-changed=", round(mean(same_sign_pairs$higher_load_changed)*100,1), "%\n\n", sep="")
cat("=== (B) Dyad-year dependence ===\n")
cat("Unique dyad-year transitions:", n_unique_transitions, "of", nrow(single_b), "triad-events\n")
cat("Median triads/transition:", median(per_transition$n_triads_resolved),
    " Max:", max(per_transition$n_triads_resolved), "\n")
cat("Bootstrap 95% CI (dyad-year clustered):", round(quantile(boot_hits, 0.025)*100,1), "% to",
    round(quantile(boot_hits, 0.975)*100,1), "%\n")
cat("Dyad-year-level collapse hit rate:", round(mean(dyad_level$is_max_ul_any)*100,1), "%\n\n")
cat("=== (C) Net balance gain ===\n")
cat("Primary test n=", nrow(strict_nbg), ", hit rate=", round(mean(strict_nbg$is_max_nbg)*100,1), "%\n", sep="")
cat("AIC unbalanced load:", round(AIC(clog_load),2), " AIC net balance gain:", round(AIC(clog_nbg),2), "\n")
sink()

saveRDS(list(single_a = single_a, clog_sign = clog_sign, clog_load = clog_load, clog_both = clog_both,
             same_sign_pairs = same_sign_pairs, single_b = single_b, per_transition = per_transition,
             boot_hits = boot_hits, dyad_level = dyad_level, strict_nbg = strict_nbg,
             clog_nbg = clog_nbg, clog_compare = clog_compare),
        "results/h2_sign_and_dependence.rds")
cat("\n22_h2_sign_and_dependence.R complete. See results/h2_sign_and_dependence_output.txt\n")
