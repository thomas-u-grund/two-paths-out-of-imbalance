# ------------------------------------------------------------------
# 07: Competing-risks (discrete-time, cause-specific hazard) model.
#     This is the paper's primary specification, in preference to the
#     fixed 5-year-window logistic decomposition (scripts 03-04): it uses
#     the full duration of each spell rather than a single fixed window.
#
#     Each unbalanced spell is expanded into one person-year row per
#     year at risk, with TIME-VARYING load covariates (updated every
#     year the spell continues, not frozen at onset). Two competing
#     events can end a spell: dissolution (a tie disappears) or flipping
#     (the triad becomes balanced). Continued imbalance at the end of
#     the observation window (2012) is right-censoring, not a third
#     outcome category -- "freezing" is degrees of survival in the
#     unbalanced state, visualized via a cumulative incidence curve,
#     not a separately-hazarded event.
# ------------------------------------------------------------------
suppressMessages({
  library(dplyr); library(readr); library(sandwich); library(lmtest); library(tidyr)
})

# Run from the replication pack root, e.g.:
#   cd replication_pack && Rscript scripts/07_competing_risks.R
# (relative paths below assume that working directory; no setwd() needed)

tr <- readRDS("data/triads_all_years.rds") %>%
  select(year, triad_id, balanced, actor_load_mean, tie_btw_mean, tie_emb_mean,
         triad_log_cinc_mean) %>%
  arrange(triad_id, year)

LAST_YEAR <- max(tr$year)  # 2012

# --- identify spells (vectorized within each triad_id) ---
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
cat("Person-year candidate rows (all unbalanced triad-years):", nrow(spells), "\n")
cat("Unique spells:", n_distinct(spells$spell_id), "\n")

# --- terminal year & event per spell ---
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
      TRUE                            ~ "ERROR_still_unbalanced"  # should not occur by construction
    )
  )
stopifnot(!any(spell_bounds$event == "ERROR_still_unbalanced"))
cat("\nSpell terminal events:\n"); print(table(spell_bounds$event))

# --- person-year risk table ---
person_years <- spells %>%
  left_join(spell_bounds %>% select(spell_id, t1, event), by = "spell_id") %>%
  mutate(
    is_terminal_row = year == t1,
    row_event = ifelse(is_terminal_row, event, "continues")
  ) %>%
  filter(row_event != "censored_end_of_panel")  # can't observe the Y->Y+1 transition for these

cat("\nUsable person-year rows (transition to Y+1 observed):", nrow(person_years), "\n")
print(table(person_years$row_event))

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

write_csv(person_years, "data/person_year_risk_table.csv")
saveRDS(person_years, "data/person_year_risk_table.rds")

# ------------------------------------------------------------------
# Cause-specific discrete-time hazard models
# ------------------------------------------------------------------
m_haz_dissolve <- glm(event_dissolve ~ z_tie_btw + z_tie_emb + z_actor_load + z_cinc + era,
                       data = person_years, family = binomial())
m_haz_flip <- glm(event_flip ~ z_tie_btw + z_tie_emb + z_actor_load + z_cinc + era,
                   data = person_years, family = binomial())

ct_haz_dissolve <- coeftest(m_haz_dissolve, vcov = vcovCL(m_haz_dissolve, cluster = person_years$triad_id))
ct_haz_flip     <- coeftest(m_haz_flip,     vcov = vcovCL(m_haz_flip,     cluster = person_years$triad_id))

# H2 on hazard framing: minimum tie load predicts the flip hazard specifically
person_years <- person_years %>% mutate(
  min_tie_load_proxy = pmin(tie_btw_mean, tie_emb_mean, na.rm = TRUE),
  z_min_tie = zscore(min_tie_load_proxy)
)
m_haz_flip_h2 <- glm(event_flip ~ z_min_tie + z_actor_load + z_cinc + era, data = person_years, family = binomial())
ct_haz_flip_h2 <- coeftest(m_haz_flip_h2, vcov = vcovCL(m_haz_flip_h2, cluster = person_years$triad_id))

sink("results/competing_risks_output.txt")
cat("=== Person-year risk table ===\n")
cat("Rows:", nrow(person_years), " | unique spells:", n_distinct(person_years$spell_id), "\n")
print(table(person_years$row_event))

cat("\n\n=== Cause-specific hazard: DISSOLUTION (event_dissolve ~ ...) ===\n")
print(ct_haz_dissolve)
cat("N =", nobs(m_haz_dissolve), "\n")

cat("\n\n=== Cause-specific hazard: FLIP TO BALANCE (event_flip ~ ...) ===\n")
print(ct_haz_flip)
cat("N =", nobs(m_haz_flip), "\n")

cat("\n\n=== H2 (hazard framing): minimum tie load predicts flip hazard ===\n")
print(ct_haz_flip_h2)
sink()

# ------------------------------------------------------------------
# Cumulative incidence function (descriptive): among spells starting in a
# given embeddedness tercile (measured at onset), the probability of
# having dissolved / flipped / still ongoing, by spell duration
# ------------------------------------------------------------------
onset_tercile <- spells %>% filter(year == ave(year, spell_id, FUN = min)) %>%
  distinct(spell_id, .keep_all = TRUE) %>%
  mutate(emb_tercile = ntile(tie_emb_mean, 3)) %>%
  select(spell_id, emb_tercile, onset_year = year)

spell_duration <- spell_bounds %>%
  left_join(onset_tercile, by = "spell_id") %>%
  mutate(duration = t1 - onset_year + 1) %>%
  filter(!is.na(emb_tercile))

MAXDUR <- 20
cif <- expand.grid(emb_tercile = 1:3, dur = 1:MAXDUR) %>%
  rowwise() %>%
  mutate(
    n_total = sum(spell_duration$emb_tercile == emb_tercile),
    n_dissolved = sum(spell_duration$emb_tercile == emb_tercile & spell_duration$event == "dissolved" & spell_duration$duration <= dur),
    n_flipped   = sum(spell_duration$emb_tercile == emb_tercile & spell_duration$event == "flipped"   & spell_duration$duration <= dur),
    p_dissolved = n_dissolved / n_total,
    p_flipped   = n_flipped / n_total
  ) %>% ungroup() %>%
  mutate(p_ongoing = 1 - p_dissolved - p_flipped)

write_csv(cif, "results/cumulative_incidence.csv")

saveRDS(list(m_haz_dissolve = m_haz_dissolve, ct_haz_dissolve = ct_haz_dissolve,
             m_haz_flip = m_haz_flip, ct_haz_flip = ct_haz_flip,
             m_haz_flip_h2 = m_haz_flip_h2, ct_haz_flip_h2 = ct_haz_flip_h2,
             person_years = person_years, spell_bounds = spell_bounds, cif = cif),
        "results/competing_risks_models.rds")

cat("\n07_competing_risks.R complete. See results/competing_risks_output.txt\n")
