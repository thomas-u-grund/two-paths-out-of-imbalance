# ------------------------------------------------------------------
# 14: Add spell-duration fixed effects to the primary competing-risks
#     hazard models (dissolution and realignment), and re-run the
#     trade/contiguity/regime robustness checks with the same duration
#     control added, for full internal consistency with the primary
#     specification.
#
#     The discrete-time hazard framework needs a baseline-hazard
#     specification (duration dummies, in this case), or a covariate that
#     happens to correlate with spell duration can absorb what is really
#     just "transitions cluster early in a spell's life" (visible
#     descriptively in Figure 1's cumulative incidence curves). This
#     checks whether H1 survives once duration is modeled explicitly.
#     It does.
# ------------------------------------------------------------------
suppressMessages({
  library(dplyr); library(readr); library(readxl); library(sandwich); library(lmtest); library(tidyr)
})

# Run from the replication pack root, e.g.:
#   cd replication_pack && Rscript scripts/14_duration_dependence.R
# (relative paths below assume that working directory; no setwd() needed)

# --- reconstruct spell onset year per row (person_year_risk_table.rds lacks it) ---
tr <- readRDS("data/triads_all_years.rds") %>%
  select(year, triad_id, balanced) %>% arrange(triad_id, year)
tr <- tr %>%
  group_by(triad_id) %>%
  mutate(prev_year = lag(year), prev_balanced = lag(balanced),
         is_onset = !balanced & (is.na(prev_year) | prev_year != year - 1 | prev_balanced),
         spell_id_local = cumsum(is_onset)) %>%
  ungroup() %>%
  mutate(spell_id = paste(triad_id, spell_id_local, sep = "_S"))
onset <- tr %>% filter(!balanced) %>% group_by(spell_id) %>% summarise(onset_year = min(year), .groups = "drop")

py <- readRDS("data/person_year_risk_table.rds") %>%
  left_join(onset, by = "spell_id") %>%
  mutate(duration = year - onset_year + 1,
         dur_bin = cut(duration, breaks = c(0, 1, 2, 3, 5, 10, Inf),
                        labels = c("yr1", "yr2", "yr3", "yr4-5", "yr6-10", "yr11+")))
cat("Duration bin distribution:\n"); print(table(py$dur_bin))

zscore <- function(x) as.numeric(scale(x))

# ------------------------------------------------------------------
# Primary hazard models, without vs. with duration fixed effects
# ------------------------------------------------------------------
m_dis0  <- glm(event_dissolve ~ z_tie_btw + z_tie_emb + z_actor_load + z_cinc + era, data = py, family = binomial())
m_flip0 <- glm(event_flip     ~ z_tie_btw + z_tie_emb + z_actor_load + z_cinc + era, data = py, family = binomial())
m_dis1  <- glm(event_dissolve ~ z_tie_btw + z_tie_emb + z_actor_load + z_cinc + era + dur_bin, data = py, family = binomial())
m_flip1 <- glm(event_flip     ~ z_tie_btw + z_tie_emb + z_actor_load + z_cinc + era + dur_bin, data = py, family = binomial())

ct_dis0  <- coeftest(m_dis0,  vcov = vcovCL(m_dis0,  cluster = py$triad_id))
ct_flip0 <- coeftest(m_flip0, vcov = vcovCL(m_flip0, cluster = py$triad_id))
ct_dis1  <- coeftest(m_dis1,  vcov = vcovCL(m_dis1,  cluster = py$triad_id))
ct_flip1 <- coeftest(m_flip1, vcov = vcovCL(m_flip1, cluster = py$triad_id))

lrt_dis  <- anova(m_dis0, m_dis1, test = "Chisq")
lrt_flip <- anova(m_flip0, m_flip1, test = "Chisq")

# ------------------------------------------------------------------
# Robustness checks re-run with duration control added
# ------------------------------------------------------------------
tr_lookup <- readRDS("data/triads_all_years.rds") %>% select(year, triad_id, node1, node2, node3)
py2 <- py %>% left_join(tr_lookup, by = c("year", "triad_id"))

# --- trade dependence, 1870-2012 subsample ---
tie_trade <- read_csv("data/tie_trade.csv", show_col_types = FALSE) %>%
  mutate(ccode_low = as.character(ccode_low), ccode_high = as.character(ccode_high)) %>%
  select(year, ccode_low, ccode_high, trade_dependence)
get_trade <- function(df, a, b, newcol) {
  df %>% mutate(lo = pmin(.data[[a]], .data[[b]]), hi = pmax(.data[[a]], .data[[b]])) %>%
    left_join(tie_trade, by = c("year" = "year", "lo" = "ccode_low", "hi" = "ccode_high")) %>%
    rename(!!newcol := trade_dependence) %>% select(-lo, -hi)
}
py2 <- py2 %>%
  get_trade("node1", "node2", "trade12") %>% get_trade("node1", "node3", "trade13") %>% get_trade("node2", "node3", "trade23") %>%
  rowwise() %>% mutate(tie_trade_mean = mean(c(trade12, trade13, trade23), na.rm = TRUE)) %>% ungroup()

py_trade <- py2 %>% filter(!is.na(tie_trade_mean), year >= 1870) %>% mutate(z_trade = zscore(tie_trade_mean))
cat("\nTrade robustness sample:", nrow(py_trade), "person-years,", n_distinct(py_trade$triad_id), "triads\n")

m_d_trade <- glm(event_dissolve ~ z_trade + z_tie_emb + z_actor_load + z_cinc + era + dur_bin, data = py_trade, family = binomial())
m_f_trade <- glm(event_flip     ~ z_trade + z_tie_emb + z_actor_load + z_cinc + era + dur_bin, data = py_trade, family = binomial())
ct_d_trade <- coeftest(m_d_trade, vcov = vcovCL(m_d_trade, cluster = py_trade$triad_id))
ct_f_trade <- coeftest(m_f_trade, vcov = vcovCL(m_f_trade, cluster = py_trade$triad_id))

# --- contiguity + regime type ---
contig <- read_csv("raw/contiguity/DirectContiguity320/contdird.csv", show_col_types = FALSE) %>%
  mutate(ccode_low = as.character(pmin(state1no, state2no)), ccode_high = as.character(pmax(state1no, state2no))) %>%
  distinct(year, ccode_low, ccode_high) %>% mutate(contiguous = 1L)
get_contig <- function(df, a, b, newcol) {
  df %>% mutate(lo = pmin(.data[[a]], .data[[b]]), hi = pmax(.data[[a]], .data[[b]])) %>%
    left_join(contig, by = c("year" = "year", "lo" = "ccode_low", "hi" = "ccode_high")) %>%
    rename(!!newcol := contiguous) %>% select(-lo, -hi)
}
py2 <- py2 %>%
  get_contig("node1", "node2", "c12") %>% get_contig("node1", "node3", "c13") %>% get_contig("node2", "node3", "c23") %>%
  mutate(across(c(c12, c13, c23), ~ replace_na(.x, 0))) %>%
  rowwise() %>% mutate(triad_contig_share = mean(c(c12, c13, c23))) %>% ungroup()

polity <- read_excel("raw/polity/p5v2018.xls") %>%
  transmute(ccode = as.character(ccode), year, polity2) %>%
  filter(!is.na(polity2)) %>% distinct(ccode, year, .keep_all = TRUE)
get_polity <- function(df, node, newcol) {
  df %>% left_join(polity, by = c("year" = "year", setNames("ccode", node))) %>% rename(!!newcol := polity2)
}
py2 <- py2 %>%
  get_polity("node1", "p1") %>% get_polity("node2", "p2") %>% get_polity("node3", "p3") %>%
  rowwise() %>%
  mutate(triad_min_polity2 = min(c(p1, p2, p3), na.rm = ifelse(all(is.na(c(p1,p2,p3))), TRUE, FALSE))) %>%
  ungroup() %>% mutate(triad_min_polity2 = ifelse(is.infinite(triad_min_polity2), NA, triad_min_polity2))

py_ir <- py2 %>% filter(!is.na(triad_min_polity2)) %>%
  mutate(z_contig = zscore(triad_contig_share), z_polity = zscore(triad_min_polity2))
cat("IR-controls robustness sample:", nrow(py_ir), "person-years,", n_distinct(py_ir$triad_id), "triads\n")

m_d_ir <- glm(event_dissolve ~ z_tie_btw + z_tie_emb + z_actor_load + z_cinc + z_contig + z_polity + era + dur_bin,
              data = py_ir, family = binomial())
m_f_ir <- glm(event_flip ~ z_tie_btw + z_tie_emb + z_actor_load + z_cinc + z_contig + z_polity + era + dur_bin,
              data = py_ir, family = binomial())
ct_d_ir <- coeftest(m_d_ir, vcov = vcovCL(m_d_ir, cluster = py_ir$triad_id))
ct_f_ir <- coeftest(m_f_ir, vcov = vcovCL(m_f_ir, cluster = py_ir$triad_id))

# ------------------------------------------------------------------
sink("results/duration_dependence_output.txt")
cat("=== DISSOLUTION, no duration control ===\n"); print(ct_dis0)
cat("\n=== DISSOLUTION, with duration dummies (PRIMARY) ===\n"); print(ct_dis1)
cat("\n=== REALIGNMENT, no duration control ===\n"); print(ct_flip0)
cat("\n=== REALIGNMENT, with duration dummies (PRIMARY) ===\n"); print(ct_flip1)
cat("\n=== LRT, dissolution: duration dummies jointly significant? ===\n"); print(lrt_dis)
cat("\n=== LRT, realignment: duration dummies jointly significant? ===\n"); print(lrt_flip)
cat("\n\n=== Trade robustness, WITH duration control, dissolution ===\n"); cat("N =", nobs(m_d_trade), "\n"); print(ct_d_trade)
cat("\n=== Trade robustness, WITH duration control, realignment ===\n"); cat("N =", nobs(m_f_trade), "\n"); print(ct_f_trade)
cat("\n\n=== Contiguity+regime robustness, WITH duration control, dissolution ===\n"); cat("N =", nobs(m_d_ir), "\n"); print(ct_d_ir)
cat("\n=== Contiguity+regime robustness, WITH duration control, realignment ===\n"); cat("N =", nobs(m_f_ir), "\n"); print(ct_f_ir)
sink()

cat(readLines("results/duration_dependence_output.txt"), sep = "\n")

saveRDS(list(m_dis0 = m_dis0, m_flip0 = m_flip0, m_dis1 = m_dis1, m_flip1 = m_flip1,
             ct_dis1 = ct_dis1, ct_flip1 = ct_flip1, lrt_dis = lrt_dis, lrt_flip = lrt_flip,
             m_d_trade = m_d_trade, m_f_trade = m_f_trade, ct_d_trade = ct_d_trade, ct_f_trade = ct_f_trade,
             m_d_ir = m_d_ir, m_f_ir = m_f_ir, ct_d_ir = ct_d_ir, ct_f_ir = ct_f_ir),
        "results/duration_dependence_models.rds")
cat("\n14_duration_dependence.R complete. See results/duration_dependence_output.txt\n")
