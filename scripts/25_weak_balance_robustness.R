# ------------------------------------------------------------------
# 25: Robustness check -- does H1 survive redefining imbalance under
#     Davis's (1967) weak-balance / clusterability criterion instead of
#     Cartwright and Harary's (1956) strong-balance criterion?
#
#     Strong balance treats both the all-negative triad (---) and the
#     two-positive/one-negative triad (++-) as unbalanced (sign product
#     negative). Weak balance treats only ++- as unbalanced; an
#     all-negative triad is balanced under weak balance (three mutual
#     enemies can each belong to their own faction; enemies of an enemy
#     need not be friends). This script reclassifies balance under the
#     weak criterion, rebuilds the unbalanced-spell structure exactly as
#     in 07_competing_risks.R, and re-estimates the cause-specific hazard
#     models to see whether the paper's H1 result is an artifact of
#     including all-negative triads as "unbalanced."
#
#     FIX (external review, second pass): this script originally omitted
#     the spell-duration fixed effects (dur_bin) that 14_duration_
#     dependence.R added to the PRIMARY specification reported in Table
#     2, and compared against 07_competing_risks.R's no-duration output
#     instead of 14's actual primary (with-duration) models. Both sides
#     of the Appendix A13/Table A14 comparison were therefore silently
#     using the pre-duration-control specification. Fixed here: duration
#     bins are reconstructed the same way as in script 14, added to both
#     hazard formulas, and the comparison now pulls the correct primary
#     estimates from results/duration_dependence_models.rds (m_dis1/
#     m_flip1) instead of results/competing_risks_models.rds.
# ------------------------------------------------------------------
suppressMessages({
  library(dplyr); library(readr); library(sandwich); library(lmtest); library(tidyr)
})


tr_raw <- readRDS("data/triads_all_years.rds")

neg_count <- with(tr_raw, (sign12 == -1) + (sign13 == -1) + (sign23 == -1))
cat("Among all triad-years classified unbalanced under strong balance:\n")
cat("  neg_count==1 (++-):", sum(!tr_raw$balanced & neg_count == 1), "\n")
cat("  neg_count==3 (---):", sum(!tr_raw$balanced & neg_count == 3), "\n")

tr <- tr_raw %>%
  mutate(balanced_strong = balanced,
         balanced = balanced_strong | (neg_count == 3)) %>%   # weak balance: --- reclassified as balanced
  select(year, triad_id, balanced, actor_load_mean, tie_btw_mean, tie_emb_mean,
         triad_log_cinc_mean) %>%
  arrange(triad_id, year)

cat("\nTriad-years reclassified from unbalanced (strong) to balanced (weak):",
    sum(tr_raw$balanced == FALSE & neg_count == 3), "of", sum(!tr_raw$balanced),
    "strong-unbalanced triad-years (",
    round(100 * sum(tr_raw$balanced == FALSE & neg_count == 3) / sum(!tr_raw$balanced), 2), "%)\n")

LAST_YEAR <- max(tr$year)  # 2012

# --- identify spells (vectorized within each triad_id) -- identical logic to 07 ---
tr <- tr %>%
  group_by(triad_id) %>%
  mutate(
    prev_year = lag(year), prev_balanced = lag(balanced),
    is_onset = !balanced & (is.na(prev_year) | prev_year != year - 1 | prev_balanced),
    spell_id_local = cumsum(is_onset)
  ) %>%
  ungroup() %>%
  mutate(spell_id = paste(triad_id, spell_id_local, sep = "_S"))

spells <- tr %>% filter(!balanced)
cat("\n[Weak balance] Person-year candidate rows (all unbalanced triad-years):", nrow(spells), "\n")
cat("[Weak balance] Unique spells:", n_distinct(spells$spell_id), "\n")

spell_bounds <- spells %>% group_by(spell_id, triad_id) %>% summarise(t1 = max(year), .groups = "drop")

next_status <- tr %>% select(year, triad_id, balanced) %>% rename(next_year = year, balanced_next = balanced)

spell_bounds <- spell_bounds %>%
  mutate(lookup_year = t1 + 1) %>%
  left_join(next_status, by = c("triad_id" = "triad_id", "lookup_year" = "next_year")) %>%
  mutate(
    event = case_when(
      t1 >= LAST_YEAR                 ~ "censored_end_of_panel",
      is.na(balanced_next)            ~ "dissolved",
      balanced_next                   ~ "flipped",
      TRUE                            ~ "ERROR_still_unbalanced"
    )
  )
stopifnot(!any(spell_bounds$event == "ERROR_still_unbalanced"))
cat("\n[Weak balance] Spell terminal events:\n"); print(table(spell_bounds$event))

person_years <- spells %>%
  left_join(spell_bounds %>% select(spell_id, t1, event), by = "spell_id") %>%
  mutate(
    is_terminal_row = year == t1,
    row_event = ifelse(is_terminal_row, event, "continues")
  ) %>%
  filter(row_event != "censored_end_of_panel")

cat("\n[Weak balance] Usable person-year rows:", nrow(person_years), "\n")
print(table(person_years$row_event))

# --- spell-duration fixed effects, reconstructed exactly as in 14_duration_dependence.R ---
onset <- spells %>% group_by(spell_id) %>% summarise(onset_year = min(year), .groups = "drop")
person_years <- person_years %>%
  left_join(onset, by = "spell_id") %>%
  mutate(duration = year - onset_year + 1,
         dur_bin = cut(duration, breaks = c(0, 1, 2, 3, 5, 10, Inf),
                        labels = c("yr1", "yr2", "yr3", "yr4-5", "yr6-10", "yr11+")))
cat("\n[Weak balance] Duration bin distribution:\n"); print(table(person_years$dur_bin))

zscore <- function(x) as.numeric(scale(x))
person_years <- person_years %>% mutate(
  z_actor_load = zscore(actor_load_mean),
  z_tie_btw    = zscore(tie_btw_mean),
  z_tie_emb    = zscore(tie_emb_mean),
  z_cinc       = zscore(triad_log_cinc_mean),
  era = cut(year, breaks = c(1815, 1945, 1989, 2015), labels = c("pre-1945", "1945-1989", "post-1989")),
  event_dissolve = as.numeric(row_event == "dissolved"),
  event_flip     = as.numeric(row_event == "flipped")
)

m_haz_dissolve <- glm(event_dissolve ~ z_tie_btw + z_tie_emb + z_actor_load + z_cinc + era + dur_bin,
                       data = person_years, family = binomial())
m_haz_flip <- glm(event_flip ~ z_tie_btw + z_tie_emb + z_actor_load + z_cinc + era + dur_bin,
                   data = person_years, family = binomial())

ct_haz_dissolve <- coeftest(m_haz_dissolve, vcov = vcovCL(m_haz_dissolve, cluster = person_years$triad_id))
ct_haz_flip     <- coeftest(m_haz_flip,     vcov = vcovCL(m_haz_flip,     cluster = person_years$triad_id))

# --- side-by-side comparison against the paper's ACTUAL primary (strong-balance,
#     WITH duration fixed effects) estimates -- the fix: this used to read from
#     competing_risks_models.rds (script 07), which has no duration control and
#     does not match Table 2. The correct primary models (m_dis1/m_flip1) live in
#     duration_dependence_models.rds (script 14). ---
orig_all <- readRDS("results/duration_dependence_models.rds")
orig <- list(ct_haz_dissolve = orig_all$ct_dis1, ct_haz_flip = orig_all$ct_flip1)

compare_coef <- function(term, orig_ct, new_ct) {
  o <- orig_ct[term, c("Estimate", "Pr(>|z|)")]
  n <- new_ct[term, c("Estimate", "Pr(>|z|)")]
  cat(sprintf("  %-14s strong: b=%7.4f p=%6.4f   |   weak: b=%7.4f p=%6.4f\n",
              term, o[1], o[2], n[1], n[2]))
}

sink("results/weak_balance_robustness_output.txt")
cat("=== Weak-balance robustness: DISSOLUTION hazard ===\n")
print(ct_haz_dissolve)
cat("N =", nobs(m_haz_dissolve), "\n\n")
cat("=== Weak-balance robustness: FLIP-TO-BALANCE (realignment) hazard ===\n")
print(ct_haz_flip)
cat("N =", nobs(m_haz_flip), "\n\n")

cat("=== Side-by-side: strong balance (paper) vs. weak balance (this check) ===\n")
cat("\n--- Dissolution hazard ---\n")
for (term in c("z_tie_btw", "z_tie_emb", "z_actor_load", "z_cinc")) {
  compare_coef(term, orig$ct_haz_dissolve, ct_haz_dissolve)
}
cat("\n--- Flip-to-balance (realignment) hazard ---\n")
for (term in c("z_tie_btw", "z_tie_emb", "z_actor_load", "z_cinc")) {
  compare_coef(term, orig$ct_haz_flip, ct_haz_flip)
}
sink()
cat(readLines("results/weak_balance_robustness_output.txt"), sep = "\n")

saveRDS(list(m_haz_dissolve = m_haz_dissolve, ct_haz_dissolve = ct_haz_dissolve,
             m_haz_flip = m_haz_flip, ct_haz_flip = ct_haz_flip,
             person_years = person_years, spell_bounds = spell_bounds),
        "results/weak_balance_robustness.rds")

cat("\n\n25_weak_balance_robustness.R complete. See results/weak_balance_robustness_output.txt\n")
