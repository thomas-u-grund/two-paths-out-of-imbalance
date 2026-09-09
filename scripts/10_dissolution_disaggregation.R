# ------------------------------------------------------------------
# 10: Disaggregate "dissolution" by which specific tie disappeared and
#     what type it was (alliance vs. militarized dispute) at the time,
#     and whether it coincides with one of the two endpoint states
#     exiting the interstate system entirely. "Dissolution" as used
#     elsewhere pools alliance expiration, MID non-recurrence, and state
#     death, so it's worth checking these aren't driving the result in
#     substantively different ways.
# ------------------------------------------------------------------
suppressMessages({ library(dplyr); library(readr) })

# Run from the replication pack root, e.g.:
#   cd replication_pack && Rscript scripts/10_dissolution_disaggregation.R
# (relative paths below assume that working directory; no setwd() needed)

m <- readRDS("results/competing_risks_models.rds")
spell_bounds <- m$spell_bounds
tr <- readRDS("data/triads_all_years.rds") %>%
  select(year, triad_id, node1, node2, node3, sign12, sign13, sign23)
state_years <- read_csv("data/state_years.csv", show_col_types = FALSE) %>%
  mutate(ccode = as.character(ccode))

dissolved <- spell_bounds %>% filter(event == "dissolved")
cat("Dissolution events to analyze:", nrow(dissolved), "\n")

at_t1 <- tr %>% rename(sign12_t1 = sign12, sign13_t1 = sign13, sign23_t1 = sign23)

# For a dissolved spell, the triad is NOT present at all at lookup_year
# (that's the definition of dissolution -- see script 07). Identify which
# specific edge disappeared by checking, in the raw signed_ties panel
# (not just the triad-closure table), which of the 3 pairs is no longer
# signed at lookup_year.
signed_ties <- read_csv("data/signed_ties.csv", show_col_types = FALSE) %>%
  mutate(ccode_low = as.character(ccode_low), ccode_high = as.character(ccode_high))

check_tie_present <- function(df, a, b, year_col, newcol) {
  df %>% mutate(lo = pmin(.data[[a]], .data[[b]]), hi = pmax(.data[[a]], .data[[b]])) %>%
    left_join(signed_ties %>% select(year, ccode_low, ccode_high, sign) %>% rename(chk_sign = sign),
              by = c(setNames("year", year_col), "lo" = "ccode_low", "hi" = "ccode_high")) %>%
    rename(!!newcol := chk_sign) %>% select(-lo, -hi)
}

events <- dissolved %>%
  left_join(at_t1, by = c("triad_id", "t1" = "year")) %>%
  mutate(check_year = lookup_year) %>%
  check_tie_present("node1", "node2", "check_year", "sign12_t2") %>%
  check_tie_present("node1", "node3", "check_year", "sign13_t2") %>%
  check_tie_present("node2", "node3", "check_year", "sign23_t2") %>%
  mutate(
    gone12 = is.na(sign12_t2), gone13 = is.na(sign13_t2), gone23 = is.na(sign23_t2),
    n_gone = gone12 + gone13 + gone23
  )

cat("\nDistribution of how many of the triad's 3 ties disappeared:\n")
print(table(events$n_gone))

single <- events %>% filter(n_gone == 1) %>%
  mutate(
    dissolved_sign = case_when(gone12 ~ sign12_t1, gone13 ~ sign13_t1, gone23 ~ sign23_t1),
    dissolved_tie_type = ifelse(dissolved_sign == 1, "alliance (positive tie)", "MID (negative tie)"),
    endpoint1 = case_when(gone12 ~ node1, gone13 ~ node1, gone23 ~ node2),
    endpoint2 = case_when(gone12 ~ node2, gone13 ~ node3, gone23 ~ node3)
  )
cat("\nEvents with exactly one tie disappearing:", nrow(single), "of", nrow(events),
    sprintf("(%.1f%%)\n", 100 * nrow(single) / nrow(events)))

cat("\n=== Which type of tie disappears when a spell dissolves? ===\n")
print(table(single$dissolved_tie_type))
print(prop.table(table(single$dissolved_tie_type)))

# --- did either endpoint exit the state system that year? ---
last_year_in_system <- state_years %>% group_by(ccode) %>% summarise(last_year = max(year), .groups = "drop")
single <- single %>%
  left_join(last_year_in_system, by = c("endpoint1" = "ccode")) %>% rename(last1 = last_year) %>%
  left_join(last_year_in_system, by = c("endpoint2" = "ccode")) %>% rename(last2 = last_year) %>%
  mutate(state_exit = (t1 == last1) | (t1 == last2))

cat("\n=== Of tie disappearances, how many coincide with a state exiting the system? ===\n")
cat(sum(single$state_exit, na.rm = TRUE), "of", nrow(single),
    sprintf("(%.1f%%)\n", 100 * mean(single$state_exit, na.rm = TRUE)))

cat("\n=== Cross-tab: tie type x state exit ===\n")
print(table(single$dissolved_tie_type, single$state_exit))

sink("results/dissolution_disaggregation_output.txt")
cat("Dissolution events:", nrow(events), " | single-tie-disappearance:", nrow(single),
    sprintf("(%.1f%%)\n\n", 100*nrow(single)/nrow(events)))
cat("Tie type that disappeared:\n"); print(table(single$dissolved_tie_type))
print(round(prop.table(table(single$dissolved_tie_type)) * 100, 1))
cat("\nCoincides with a state exiting the interstate system:\n")
cat(sum(single$state_exit, na.rm=TRUE), "of", nrow(single),
    sprintf("(%.1f%%)\n", 100*mean(single$state_exit, na.rm=TRUE)))
cat("\nCross-tab (tie type x state exit):\n"); print(table(single$dissolved_tie_type, single$state_exit))
sink()

saveRDS(list(events = events, single = single), "results/dissolution_disaggregation.rds")
cat("\n10_dissolution_disaggregation.R complete.\n")
