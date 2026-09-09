# ------------------------------------------------------------------
# 24: H1 realignment hazard, controlling for the triad's CURRENT
#     unbalanced load (mean, across the triad's 3 edges, of the number
#     of OTHER closed triads sharing that edge that are themselves
#     unbalanced in year t), alongside total triadic embeddedness.
#
#     The paper's own H1 argument reframes total embeddedness as a
#     triad-level EXPOSURE
#     measure (how many surrounding configurations could eventually turn
#     unbalanced), distinct from unbalanced load's INSTANTANEOUS pressure
#     (Section 2.2-2.3). That raises the natural question of whether
#     embeddedness's positive coefficient on the realignment hazard (H1)
#     is partly just picking up triads that already happen to have high
#     CURRENT unbalanced load right now -- i.e., whether H1 and H2 are
#     silently the same mechanism at different levels of aggregation.
#     This script tests that directly: does z_tie_emb survive adding
#     the triad's own current unbalanced load as a covariate in the
#     same discrete-time hazard specification used for Table 2?
# ------------------------------------------------------------------
suppressMessages({
  library(dplyr); library(readr); library(sandwich); library(lmtest); library(tidyr)
})

# Run from the replication pack root, e.g.:
#   cd replication_pack && Rscript scripts/24_h1_unbalanced_load_control.R
# (relative paths below assume that working directory; no setwd() needed)

# --- full (year, edge) -> unbalanced-triangle-count lookup, built from
#     the COMPLETE triad panel (not just the H2 realignment subsample
#     script 12 uses) so it covers every triad-year at risk in H1 ---
tr_full <- readRDS("data/triads_all_years.rds") %>%
  select(year, triad_id, node1, node2, node3, balanced)

edge_long <- bind_rows(
  tr_full %>% transmute(year, triad_id, balanced, e_lo = pmin(node1, node2), e_hi = pmax(node1, node2)),
  tr_full %>% transmute(year, triad_id, balanced, e_lo = pmin(node1, node3), e_hi = pmax(node1, node3)),
  tr_full %>% transmute(year, triad_id, balanced, e_lo = pmin(node2, node3), e_hi = pmax(node2, node3))
)
edge_summary <- edge_long %>%
  group_by(year, e_lo, e_hi) %>%
  summarise(n_unbalanced = sum(!balanced), .groups = "drop")
cat("Unique (year, edge) cells:", nrow(edge_summary), "\n")

get_unbal_load <- function(df, a, b, newcol) {
  df %>% mutate(lo = pmin(.data[[a]], .data[[b]]), hi = pmax(.data[[a]], .data[[b]])) %>%
    left_join(edge_summary, by = c("year" = "year", "lo" = "e_lo", "hi" = "e_hi")) %>%
    rename(!!newcol := n_unbalanced) %>%
    mutate(!!newcol := replace_na(!!sym(newcol), 0)) %>%
    select(-lo, -hi)
}

# --- primary H1 person-year table (same rows/spec as script 14 / Table 2) ---
py <- readRDS("data/person_year_risk_table_with_duration.rds") %>%
  left_join(tr_full %>% select(year, triad_id, node1, node2, node3), by = c("year", "triad_id")) %>%
  mutate(dur_bin = cut(duration, breaks = c(0, 1, 2, 3, 5, 10, Inf),
                        labels = c("yr1", "yr2", "yr3", "yr4-5", "yr6-10", "yr11+")))

py <- py %>%
  get_unbal_load("node1", "node2", "ul12") %>%
  get_unbal_load("node1", "node3", "ul13") %>%
  get_unbal_load("node2", "node3", "ul23") %>%
  mutate(
    # subtract 1 from each edge's count: every row here is itself an
    # unbalanced triad by construction (it's a person-year of a live
    # imbalance spell), so it is counted once in its own edges' totals
    # and must be removed to leave OTHER unbalanced triads only, exactly
    # as in script 12's H2 construction.
    ul12 = pmax(ul12 - 1, 0), ul13 = pmax(ul13 - 1, 0), ul23 = pmax(ul23 - 1, 0),
    triad_unbal_load_mean = (ul12 + ul13 + ul23) / 3
  )

zscore <- function(x) as.numeric(scale(x))
py <- py %>% mutate(z_unbal_load = zscore(triad_unbal_load_mean))

cat("\nCorrelation, triadic embeddedness vs. current triad-level unbalanced load: r =",
    round(cor(py$z_tie_emb, py$z_unbal_load, use = "complete.obs"), 3), "\n")
cat("Summary of triad_unbal_load_mean:\n"); print(summary(py$triad_unbal_load_mean))

# ------------------------------------------------------------------
# Realignment ("flip") hazard: baseline (Table 2 spec) vs. + unbalanced load
# ------------------------------------------------------------------
m_flip_base <- glm(event_flip ~ z_tie_btw + z_tie_emb + z_actor_load + z_cinc + era + dur_bin,
                    data = py, family = binomial())
m_flip_ul   <- glm(event_flip ~ z_tie_btw + z_tie_emb + z_actor_load + z_cinc + era + dur_bin + z_unbal_load,
                    data = py, family = binomial())

ct_flip_base <- coeftest(m_flip_base, vcov = vcovCL(m_flip_base, cluster = py$triad_id))
ct_flip_ul   <- coeftest(m_flip_ul,   vcov = vcovCL(m_flip_ul,   cluster = py$triad_id))

lrt_flip <- anova(m_flip_base, m_flip_ul, test = "Chisq")

# --- same check on the dissolution hazard, for completeness/symmetry ---
m_dis_base <- glm(event_dissolve ~ z_tie_btw + z_tie_emb + z_actor_load + z_cinc + era + dur_bin,
                   data = py, family = binomial())
m_dis_ul   <- glm(event_dissolve ~ z_tie_btw + z_tie_emb + z_actor_load + z_cinc + era + dur_bin + z_unbal_load,
                   data = py, family = binomial())
ct_dis_base <- coeftest(m_dis_base, vcov = vcovCL(m_dis_base, cluster = py$triad_id))
ct_dis_ul   <- coeftest(m_dis_ul,   vcov = vcovCL(m_dis_ul,   cluster = py$triad_id))
lrt_dis <- anova(m_dis_base, m_dis_ul, test = "Chisq")

sink("results/h1_unbalanced_load_control_output.txt")
cat("=== Correlation: z_tie_emb vs. z_unbal_load (triad-level, current year) ===\n")
cat("r =", round(cor(py$z_tie_emb, py$z_unbal_load, use = "complete.obs"), 3), "\n\n")

cat("=== REALIGNMENT hazard, baseline (Table 2 spec) ===\n"); print(ct_flip_base)
cat("N =", nobs(m_flip_base), "\n\n")

cat("=== REALIGNMENT hazard, + triad-level current unbalanced load ===\n"); print(ct_flip_ul)
cat("N =", nobs(m_flip_ul), "\n\n")

cat("=== LRT: does adding current unbalanced load improve fit (realignment)? ===\n")
print(lrt_flip)

cat("\n\n=== DISSOLUTION hazard, baseline (Table 2 spec) ===\n"); print(ct_dis_base)
cat("N =", nobs(m_dis_base), "\n\n")

cat("=== DISSOLUTION hazard, + triad-level current unbalanced load ===\n"); print(ct_dis_ul)
cat("N =", nobs(m_dis_ul), "\n\n")

cat("=== LRT: does adding current unbalanced load improve fit (dissolution)? ===\n")
print(lrt_dis)
sink()

cat(readLines("results/h1_unbalanced_load_control_output.txt"), sep = "\n")

saveRDS(list(py = py, m_flip_base = m_flip_base, m_flip_ul = m_flip_ul,
             ct_flip_base = ct_flip_base, ct_flip_ul = ct_flip_ul, lrt_flip = lrt_flip,
             m_dis_base = m_dis_base, m_dis_ul = m_dis_ul,
             ct_dis_base = ct_dis_base, ct_dis_ul = ct_dis_ul, lrt_dis = lrt_dis),
        "results/h1_unbalanced_load_control.rds")
cat("\n24_h1_unbalanced_load_control.R complete. See results/h1_unbalanced_load_control_output.txt\n")
