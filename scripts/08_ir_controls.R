# ------------------------------------------------------------------
# 08: Robustness -- add contiguity and regime-type (Polity5) controls
#     to the primary competing-risks hazard models, since both are
#     standard confounds in interstate-conflict/alliance research and
#     their omission would be a natural first objection.
# ------------------------------------------------------------------
suppressMessages({
  library(dplyr); library(readr); library(readxl); library(sandwich); library(lmtest); library(tidyr)
})

# Run from the replication pack root, e.g.:
#   cd replication_pack && Rscript scripts/08_ir_controls.R
# (relative paths below assume that working directory; no setwd() needed)

py <- readRDS("data/person_year_risk_table.rds")
tr_lookup <- readRDS("data/triads_all_years.rds") %>% select(year, triad_id, node1, node2, node3)
py <- py %>% left_join(tr_lookup, by = c("year", "triad_id"))

# --- Contiguity: share of the triad's 3 dyads that are directly
# contiguous that year (COW Direct Contiguity v3.2; presence in the
# file = some form of land/water contiguity, the standard operationalization) ---
contig <- read_csv("raw/contiguity/DirectContiguity320/contdird.csv", show_col_types = FALSE) %>%
  mutate(ccode_low = as.character(pmin(state1no, state2no)), ccode_high = as.character(pmax(state1no, state2no))) %>%
  distinct(year, ccode_low, ccode_high) %>%
  mutate(contiguous = 1L)

get_contig <- function(df, a, b, newcol) {
  df %>% mutate(lo = pmin(.data[[a]], .data[[b]]), hi = pmax(.data[[a]], .data[[b]])) %>%
    left_join(contig, by = c("year" = "year", "lo" = "ccode_low", "hi" = "ccode_high")) %>%
    rename(!!newcol := contiguous) %>% select(-lo, -hi)
}
py <- py %>%
  get_contig("node1", "node2", "c12") %>% get_contig("node1", "node3", "c13") %>% get_contig("node2", "node3", "c23") %>%
  mutate(across(c(c12, c13, c23), ~ replace_na(.x, 0))) %>%
  rowwise() %>% mutate(triad_contig_share = mean(c(c12, c13, c23))) %>% ungroup()

# --- Regime type: Polity5 polity2 score, min across the triad's 3 members
# ("weakest link" operationalization, standard in the joint-democracy
# literature and robust to a single well-governed outlier) ---
polity <- read_excel("raw/polity/p5v2018.xls") %>%
  transmute(ccode = as.character(ccode), year, polity2) %>%
  filter(!is.na(polity2)) %>%
  distinct(ccode, year, .keep_all = TRUE)

get_polity <- function(df, node, newcol) {
  df %>% left_join(polity, by = c("year" = "year", setNames("ccode", node))) %>%
    rename(!!newcol := polity2)
}
py <- py %>%
  get_polity("node1", "p1") %>% get_polity("node2", "p2") %>% get_polity("node3", "p3") %>%
  rowwise() %>%
  mutate(triad_min_polity2 = min(c(p1, p2, p3), na.rm = ifelse(all(is.na(c(p1,p2,p3))), TRUE, FALSE))) %>%
  ungroup() %>%
  mutate(triad_min_polity2 = ifelse(is.infinite(triad_min_polity2), NA, triad_min_polity2))

cat("N with contiguity data:", sum(!is.na(py$triad_contig_share)), "of", nrow(py), "\n")
cat("N with polity data:", sum(!is.na(py$triad_min_polity2)), "of", nrow(py), "\n")

zscore <- function(x) as.numeric(scale(x))
py_ir <- py %>% filter(!is.na(triad_min_polity2)) %>%  # Polity coverage starts later than the full panel
  mutate(z_contig = zscore(triad_contig_share), z_polity = zscore(triad_min_polity2))

cat("\nRobustness sample (non-missing IR controls):", nrow(py_ir), "person-years,",
    n_distinct(py_ir$triad_id), "triads, years", min(py_ir$year), "-", max(py_ir$year), "\n")

m_d_ir <- glm(event_dissolve ~ z_tie_btw + z_tie_emb + z_actor_load + z_cinc + z_contig + z_polity + era,
              data = py_ir, family = binomial())
m_f_ir <- glm(event_flip ~ z_tie_btw + z_tie_emb + z_actor_load + z_cinc + z_contig + z_polity + era,
              data = py_ir, family = binomial())
ct_d_ir <- coeftest(m_d_ir, vcov = vcovCL(m_d_ir, cluster = py_ir$triad_id))
ct_f_ir <- coeftest(m_f_ir, vcov = vcovCL(m_f_ir, cluster = py_ir$triad_id))

sink("results/ir_controls_output.txt")
cat("=== Dissolution hazard, with contiguity + regime-type controls ===\n")
cat("N =", nobs(m_d_ir), "\n")
print(ct_d_ir)
cat("\n\n=== Realignment hazard, with contiguity + regime-type controls ===\n")
cat("N =", nobs(m_f_ir), "\n")
print(ct_f_ir)
sink()

cat(readLines("results/ir_controls_output.txt"), sep = "\n")

saveRDS(list(m_d_ir = m_d_ir, ct_d_ir = ct_d_ir, m_f_ir = m_f_ir, ct_f_ir = ct_f_ir, py_ir = py_ir),
        "results/ir_controls_models.rds")
cat("\n08_ir_controls.R complete.\n")
