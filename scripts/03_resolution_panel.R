# ------------------------------------------------------------------
# 03: Build the triad-year resolution panel used to test H1-H3
#     Unbalanced triad at year t -> did it resolve by t+5?
# ------------------------------------------------------------------
suppressMessages({ library(dplyr); library(readr); library(tidyr) })

# Run from the replication pack root, e.g.:
#   cd replication_pack && Rscript scripts/03_resolution_panel.R
# (relative paths below assume that working directory; no setwd() needed)

triads <- readRDS("data/triads_all_years.rds")
tie_trade <- read_csv("data/tie_trade.csv", show_col_types = FALSE)

WINDOW <- 5

unbalanced <- triads %>% filter(!balanced)
cat("Unbalanced triad-years (starting points):", nrow(unbalanced), "\n")

status_lookup <- triads %>% select(year, triad_id, balanced)

# For each unbalanced triad at year t, look up status at t+WINDOW (if the
# same triad_id — i.e. the same three states — still forms a *closed* triad
# then; if it no longer appears as a closed triad, that itself counts as
# "no longer an active unbalanced configuration" and is coded separately)
res <- unbalanced %>%
  mutate(target_year = year + WINDOW) %>%
  left_join(status_lookup, by = c("target_year" = "year", "triad_id" = "triad_id"),
            suffix = c("", "_t5")) %>%
  mutate(
    still_closed_t5 = !is.na(balanced_t5),
    resolved = case_when(
      still_closed_t5 & balanced_t5  ~ TRUE,   # became balanced
      still_closed_t5 & !balanced_t5 ~ FALSE,  # still unbalanced
      !still_closed_t5               ~ NA      # triad dissolved / one tie disappeared: ambiguous, handled separately
    )
  )

cat("Of", nrow(res), "unbalanced triad-years:\n")
cat("  still a closed triad 5y later:", sum(res$still_closed_t5), "\n")
cat("  resolved to balance (of those still closed):", sum(res$resolved, na.rm = TRUE), "\n")
cat("  dissolved (no longer closed triad) at t+5:", sum(!res$still_closed_t5), "\n")

# Primary analysis sample: triads that remain a closed triad at t+5
# (so "resolved" is a clean binary: flipped to balance vs. stayed unbalanced).
# Dissolution (a tie disappearing entirely) is a distinct outcome and is
# reported separately as a robustness check, not pooled into "resolved",
# since losing a tie entirely is a different kind of change than flipping
# its sign.
main_sample <- res %>% filter(still_closed_t5)

# --- attach trade-based tie load (secondary measure, 1870-2014 only) ---
tie_trade_lookup <- tie_trade %>%
  mutate(ccode_low = as.character(ccode_low), ccode_high = as.character(ccode_high)) %>%
  select(year, ccode_low, ccode_high, trade_dependence)

get_trade <- function(df, a, b, newcol) {
  df %>%
    mutate(lo = pmin(.data[[a]], .data[[b]]), hi = pmax(.data[[a]], .data[[b]])) %>%
    left_join(tie_trade_lookup, by = c("year" = "year", "lo" = "ccode_low", "hi" = "ccode_high")) %>%
    rename(!!newcol := trade_dependence) %>%
    select(-lo, -hi)
}

main_sample <- main_sample %>%
  get_trade("node1", "node2", "trade12") %>%
  get_trade("node1", "node3", "trade13") %>%
  get_trade("node2", "node3", "trade23") %>%
  rowwise() %>%
  mutate(tie_trade_mean = mean(c(trade12, trade13, trade23), na.rm = TRUE)) %>%
  ungroup()

write_csv(main_sample, "data/resolution_panel.csv")
saveRDS(main_sample, "data/resolution_panel.rds")

cat("Main analysis sample (still-closed at t+5):", nrow(main_sample), "\n")
cat("Resolution rate:", round(mean(main_sample$resolved) * 100, 1), "%\n")
cat("03_resolution_panel.R complete.\n")
