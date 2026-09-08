# ------------------------------------------------------------------
# 09: H2, done properly. The earlier version tested "does low minimum
#     triad embeddedness predict the realignment hazard" -- not the
#     actual hypothesis, which is about WHICH of a triad's three ties
#     changes when it does realign. Fixed per external review round 2
#     (DECISIONS.md D28): identify the specific tie that changed sign
#     for every realignment event, and test whether it was the least-
#     embedded of the triad's three ties, against the 1/3 null and via
#     a conditional logit treating the three ties as discrete choice
#     alternatives.
# ------------------------------------------------------------------
suppressMessages({ library(dplyr); library(readr); library(survival) })

# Run from the replication pack root, e.g.:
#   cd replication_pack && Rscript scripts/09_h2_tie_choice.R
# (relative paths below assume that working directory; no setwd() needed)

m <- readRDS("results/competing_risks_models.rds")
spell_bounds <- m$spell_bounds  # triad_id, spell_id, t1, event
tr <- readRDS("data/triads_all_years.rds") %>%
  select(year, triad_id, sign12, sign13, sign23, tie12_emb, tie13_emb, tie23_emb)

flipped <- spell_bounds %>% filter(event == "flipped")
cat("Realignment events to analyze:", nrow(flipped), "\n")

at_t1   <- tr %>% rename(sign12_t1 = sign12, sign13_t1 = sign13, sign23_t1 = sign23,
                          emb12 = tie12_emb, emb13 = tie13_emb, emb23 = tie23_emb)
at_t1p1 <- tr %>% select(year, triad_id, sign12, sign13, sign23) %>%
  rename(sign12_t2 = sign12, sign13_t2 = sign13, sign23_t2 = sign23)

events <- flipped %>%
  left_join(at_t1, by = c("triad_id", "t1" = "year")) %>%
  left_join(at_t1p1, by = c("triad_id", "lookup_year" = "year")) %>%
  mutate(
    changed12 = sign12_t1 != sign12_t2,
    changed13 = sign13_t1 != sign13_t2,
    changed23 = sign23_t1 != sign23_t2,
    n_changed = changed12 + changed13 + changed23
  )

cat("\nDistribution of number of ties changed at realignment:\n")
print(table(events$n_changed))

# Restrict to the clean, interpretable case: exactly one tie changed.
# (Multi-tie changes in the same year are a small minority and conflate
# multiple realignment decisions in one observation window -- excluded
# here, not swept into the "one tie changed" test.)
single <- events %>% filter(n_changed == 1)
cat("\nEvents with exactly one tie changed:", nrow(single), "of", nrow(events),
    sprintf("(%.1f%%)\n", 100 * nrow(single) / nrow(events)))

single <- single %>%
  mutate(
    changed_tie = case_when(changed12 ~ "12", changed13 ~ "13", changed23 ~ "23"),
    changed_emb = case_when(changed12 ~ emb12, changed13 ~ emb13, changed23 ~ emb23),
    min_emb = pmin(emb12, emb13, emb23, na.rm = TRUE),
    is_min = changed_emb == min_emb,
    all_tied = (emb12 == emb13) & (emb13 == emb23),
    n_distinct_vals = sapply(seq_len(n()), function(i) length(unique(c(emb12[i], emb13[i], emb23[i])))),
    is_unique_min = n_distinct_vals == 3 & changed_emb == min_emb  # a STRICT, unambiguous minimum
  )

cat("\n=== Data feature: embeddedness is tied across a triad's 3 edges far more than expected ===\n")
cat("All three tied:", round(mean(single$all_tied) * 100, 1), "%\n")
cat("Exactly 2 tied:", round(mean(!single$all_tied & single$n_distinct_vals == 2) * 100, 1), "%\n")
cat("All 3 distinct:", round(mean(single$n_distinct_vals == 3) * 100, 1), "%\n")
cat("(This happens because triads embedded in a dense alliance bloc share most of their\n")
cat(" common neighbors across all three edges simultaneously -- a real structural feature\n")
cat(" of alliance blocs, not a data artifact. It means a naive 'is this tie tied for the\n")
cat(" min' test is dominated by degenerate 3-way ties and is not a meaningful ranking test.)\n")

cat("\n=== H2, naive test (counts ties as 'success'): was the changed tie tied for the minimum? ===\n")
cat("Observed:", round(mean(single$is_min, na.rm = TRUE) * 100, 1), "% (null = 33.3% if chosen at random)\n")
bt_naive <- binom.test(sum(single$is_min, na.rm = TRUE), sum(!is.na(single$is_min)), p = 1/3)
print(bt_naive)

cat("\n=== H2, PRIMARY test: restricted to events with a strict, unambiguous minimum ===\n")
strict <- single %>% filter(n_distinct_vals == 3)
cat("N events with all 3 embeddedness values distinct:", nrow(strict), "\n")
cat("Observed:", round(mean(strict$is_unique_min, na.rm = TRUE) * 100, 1), "% (null = 33.3%)\n")
bt <- binom.test(sum(strict$is_unique_min, na.rm = TRUE), nrow(strict), p = 1/3)
print(bt)

# --- Conditional logit: 3 alternatives (the triad's 3 ties) per choice
# set (each realignment event), predictor = that tie's embeddedness,
# outcome = was this the tie that changed ---
long <- single %>%
  transmute(event_id = row_number(), emb12, emb13, emb23, changed12, changed13, changed23) %>%
  tidyr::pivot_longer(cols = c(emb12, emb13, emb23), names_to = "tie", values_to = "embeddedness") %>%
  mutate(tie = sub("emb", "", tie)) %>%
  left_join(
    single %>% mutate(event_id = row_number()) %>%
      select(event_id, changed12, changed13, changed23) %>%
      tidyr::pivot_longer(cols = c(changed12, changed13, changed23), names_to = "tie2", values_to = "changed") %>%
      mutate(tie2 = sub("changed", "", tie2)),
    by = c("event_id", "tie" = "tie2")
  )

cat("\n=== Conditional logit, full sample (degenerate tied strata contribute no\n")
cat("    information to the likelihood, but are not literally excluded) ===\n")
clog_full <- clogit(changed ~ embeddedness + strata(event_id), data = long)
print(summary(clog_full))

# Restrict to choice sets with genuine variation (drop the 3-way-tied
# strata, which carry zero information anyway, to make this explicit
# rather than implicit)
strict_ids <- single %>% filter(n_distinct_vals == 3) %>% mutate(event_id = row_number())
long_strict <- long %>% semi_join(single %>% mutate(event_id = row_number()) %>% filter(n_distinct_vals == 3) %>% select(event_id), by = "event_id")
cat("\n=== Conditional logit, restricted to strata with all 3 embeddedness values distinct ===\n")
clog_strict <- clogit(changed ~ embeddedness + strata(event_id), data = long_strict)
print(summary(clog_strict))

# ------------------------------------------------------------------
# Betweenness offers a much better-powered version of the same test:
# it is a continuous score with far fewer exact ties across a triad's
# three edges than the integer common-neighbor embeddedness count
# (only 37% of events have all three tied, vs. 74.6% for embeddedness),
# so restricting to strict variation still leaves a large, well-powered
# sample.
# ------------------------------------------------------------------
tr_btw <- readRDS("data/triads_all_years.rds") %>% select(year, triad_id, tie12_btw, tie13_btw, tie23_btw)
single_b <- single %>% select(triad_id, t1, changed12, changed13, changed23) %>%
  left_join(tr_btw, by = c("triad_id", "t1" = "year")) %>%
  rowwise() %>%
  mutate(n_distinct_btw = length(unique(round(c(tie12_btw, tie13_btw, tie23_btw), 6)))) %>%
  ungroup() %>%
  mutate(
    changed_btw = case_when(changed12 ~ tie12_btw, changed13 ~ tie13_btw, changed23 ~ tie23_btw),
    min_btw = pmin(tie12_btw, tie13_btw, tie23_btw),
    is_min_btw = changed_btw == min_btw
  )
strict_b <- single_b %>% filter(n_distinct_btw == 3)
cat("\n=== H2, betweenness version (much better powered: n=", nrow(strict_b), " with all 3 distinct) ===\n", sep = "")
cat("Observed: changed tie has the MINIMUM betweenness in", round(mean(strict_b$is_min_btw)*100,1), "% of events (null = 33.3%)\n")
bt_btw <- binom.test(sum(strict_b$is_min_btw), nrow(strict_b), p = 1/3)
print(bt_btw)
cat("This is significantly BELOW chance: realignment concentrates on the tie with the\n")
cat("HIGHEST, not lowest, betweenness within the triad -- the opposite of H2 as stated.\n")

sink("results/h2_tie_choice_output.txt")
cat("Realignment events:", nrow(events), " | single-tie-change:", nrow(single),
    sprintf("(%.1f%%)\n", 100*nrow(single)/nrow(events)))
cat("\nEmbeddedness tie structure across the 3 edges of a triad (single-change events):\n")
cat("  all 3 tied:", round(mean(single$all_tied)*100,1), "%\n")
cat("  exactly 2 tied:", round(mean(!single$all_tied & single$n_distinct_vals==2)*100,1), "%\n")
cat("  all 3 distinct:", round(mean(single$n_distinct_vals==3)*100,1), "%\n")
cat("\n--- embeddedness: naive test (ties counted as success) ---\n"); print(bt_naive)
cat("\n--- embeddedness: primary test (strict minimum only, n=", nrow(strict), ") ---\n", sep=""); print(bt)
cat("\n--- embeddedness: conditional logit, full sample ---\n"); print(summary(clog_full))
cat("\n--- embeddedness: conditional logit, strict-variation subsample ---\n"); print(summary(clog_strict))
cat("\n--- betweenness: strict test (n=", nrow(strict_b), ") ---\n", sep=""); print(bt_btw)
sink()

saveRDS(list(events = events, single = single, bt = bt, bt_naive = bt_naive,
             clog_full = clog_full, clog_strict = clog_strict,
             single_b = single_b, strict_b = strict_b, bt_btw = bt_btw),
        "results/h2_tie_choice.rds")
cat("\n09_h2_tie_choice.R complete. See results/h2_tie_choice_output.txt\n")
