# ------------------------------------------------------------------
# 01: Build signed dyad-year panel from Correlates of War data
# ------------------------------------------------------------------
library(dplyr)
library(readr)
library(tidyr)

# Run from the replication pack root, e.g.:
#   cd replication_pack && Rscript scripts/01_build_dyad_year.R
# (relative paths below assume that working directory; no setwd() needed)
raw <- "raw"

# 1. Alliances -> positive ties (defense, neutrality, nonaggression, entente)
alliances <- read_csv(file.path(raw, "alliances/version4.1_csv/alliance_v4.1_by_dyad_yearly.csv"),
                       show_col_types = FALSE) %>%
  mutate(ccode_low = pmin(ccode1, ccode2), ccode_high = pmax(ccode1, ccode2)) %>%
  filter(defense == 1 | neutrality == 1 | nonaggression == 1 | entente == 1) %>%
  mutate(sign = 1) %>%
  select(year, ccode_low, ccode_high, sign) %>%
  distinct()

# 2. MIDs -> negative ties (hostility level >= 3: display of force or higher)
mids <- read_csv(file.path(raw, "mids/dyadic_mid_4.03_update/dyadic_mid_4.03.csv"),
                  show_col_types = FALSE) %>%
  mutate(ccode_low = pmin(statea, stateb), ccode_high = pmax(statea, stateb)) %>%
  filter(hihost >= 3) %>%
  mutate(sign = -1) %>%
  select(year, ccode_low, ccode_high, sign) %>%
  distinct()

# 3. Merge: if a dyad-year has both an alliance and a MID, code as negative
#    (conflict dominates cooperation in the same year — conservative choice)
signed_ties <- bind_rows(alliances, mids) %>%
  group_by(year, ccode_low, ccode_high) %>%
  summarise(sign = ifelse(any(sign == -1), -1, 1),
            multiplex = n_distinct(sign) > 1, .groups = "drop")

cat("Signed ties built:", nrow(signed_ties), "\n")
cat("  positive:", sum(signed_ties$sign == 1), " negative:", sum(signed_ties$sign == -1), "\n")
cat("  multiplex (both alliance+MID same year):", sum(signed_ties$multiplex), "\n")

# 4. State system membership -> restrict universe to actual states in that year
states <- read_csv(file.path(raw, "states/States2024/statelist2024.csv"), show_col_types = FALSE)
state_years <- states %>%
  rowwise() %>%
  mutate(yrs = list(seq(styear, endyear))) %>%
  ungroup() %>%
  select(ccode, yrs) %>%
  unnest(yrs) %>%
  rename(year = yrs) %>%
  distinct()

write_csv(signed_ties, "data/signed_ties.csv")
write_csv(state_years, "data/state_years.csv")

# ------------------------------------------------------------------
# 5. Tie load (secondary/robustness): trade dependence, 1870-2014 only
# ------------------------------------------------------------------
dyadic_trade <- read_csv(file.path(raw, "trade/COW_Trade_4.0/Dyadic_COW_4.0.csv"), show_col_types = FALSE) %>%
  mutate(across(c(flow1, flow2, smoothtotrade), ~ na_if(., -9))) %>%
  mutate(ccode_low = pmin(ccode1, ccode2), ccode_high = pmax(ccode1, ccode2)) %>%
  group_by(year, ccode_low, ccode_high) %>%
  summarise(total_dyadic_trade = sum(smoothtotrade, na.rm = TRUE), .groups = "drop")

national_trade <- read_csv(file.path(raw, "trade/COW_Trade_4.0/National_COW_4.0.csv"), show_col_types = FALSE) %>%
  mutate(total_trade = imports + exports,
         # floor: require at least $50m (in COW's reporting units) of recorded
         # national trade before using it as a denominator, else shares from
         # tiny/misreported totals blow up above 1 (documented data-quality
         # issue found in the original prototype)
         total_trade = ifelse(total_trade < 50, NA, total_trade)) %>%
  select(ccode, year, total_trade)

tie_trade <- dyadic_trade %>%
  left_join(national_trade, by = c("ccode_low" = "ccode", "year")) %>%
  rename(total_trade_low = total_trade) %>%
  left_join(national_trade, by = c("ccode_high" = "ccode", "year")) %>%
  rename(total_trade_high = total_trade) %>%
  mutate(
    overflow = total_dyadic_trade > pmin(total_trade_low, total_trade_high, na.rm = TRUE),
    trade_share_low  = ifelse(overflow, NA, total_dyadic_trade / total_trade_low),
    trade_share_high = ifelse(overflow, NA, total_dyadic_trade / total_trade_high),
    trade_dependence = ifelse(!is.na(trade_share_low) & !is.na(trade_share_high),
                               sqrt(pmax(trade_share_low, 0) * pmax(trade_share_high, 0)), NA),
    log_trade_dependence = ifelse(!is.na(trade_dependence) & trade_dependence > 0,
                                   log(trade_dependence), NA)
  ) %>%
  select(year, ccode_low, ccode_high, trade_dependence, log_trade_dependence, overflow)

n_overflow <- sum(tie_trade$overflow, na.rm = TRUE)
cat("Trade dyad-years:", nrow(tie_trade), " | overflow (excluded, share>1):", n_overflow,
    sprintf("(%.2f%%)\n", 100 * n_overflow / nrow(tie_trade)))

write_csv(tie_trade, "data/tie_trade.csv")

cat("01_build_dyad_year.R complete.\n")
