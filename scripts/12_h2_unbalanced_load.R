# ------------------------------------------------------------------
# 12: H2, tested with a measure specific to realignment pressure. Total
#     embeddedness (all overlapping triangles, balanced + unbalanced) and
#     betweenness (bridging importance) capture structural position
#     generally, not pressure toward realignment specifically, which
#     should come from OTHER UNBALANCED triads sharing the tie. A tie can
#     be highly embedded while under little consistency pressure if most
#     of what it's embedded in is already balanced.
#
#     Measure: UnbalancedLoad(e) at year t is the number of OTHER closed
#     triads (besides the focal one) at year t that (a) contain edge e
#     and (b) are themselves unbalanced. H2: within a realigning triad,
#     the tie that changes is the one with the HIGHEST UnbalancedLoad
#     (most other unbalanced triads pulling on it).
# ------------------------------------------------------------------
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(survival); library(purrr) })

# Run from the replication pack root, e.g.:
#   cd replication_pack && Rscript scripts/12_h2_unbalanced_load.R
# (relative paths below assume that working directory; no setwd() needed)

tr <- readRDS("data/triads_all_years.rds") %>% select(year, triad_id, node1, node2, node3, balanced)
tr_btw <- readRDS("data/triads_all_years.rds") %>% select(year, triad_id, tie12_btw, tie13_btw, tie23_btw)

# --- per-(year, edge) unbalanced-triangle count ---
edge_long <- bind_rows(
  tr %>% transmute(year, triad_id, balanced, e_lo = pmin(node1, node2), e_hi = pmax(node1, node2)),
  tr %>% transmute(year, triad_id, balanced, e_lo = pmin(node1, node3), e_hi = pmax(node1, node3)),
  tr %>% transmute(year, triad_id, balanced, e_lo = pmin(node2, node3), e_hi = pmax(node2, node3))
)
cat("Edge-triad membership rows:", nrow(edge_long), "\n")

edge_summary <- edge_long %>%
  group_by(year, e_lo, e_hi) %>%
  summarise(n_triads = n(), n_unbalanced = sum(!balanced), .groups = "drop")
cat("Unique (year, edge) cells:", nrow(edge_summary), "\n")

get_unbal_load <- function(df, a, b, year_col, newcol) {
  df %>% mutate(lo = pmin(.data[[a]], .data[[b]]), hi = pmax(.data[[a]], .data[[b]])) %>%
    left_join(edge_summary, by = c(setNames("year", year_col), "lo" = "e_lo", "hi" = "e_hi")) %>%
    rename(!!newcol := n_unbalanced) %>%
    mutate(!!newcol := !!sym(newcol) - 1) %>%  # subtract the focal triad itself
    select(-lo, -hi, -n_triads)
}

# --- reuse the single-tie-change realignment events from script 09 ---
h2 <- readRDS("results/h2_tie_choice.rds")
node_lookup <- tr %>% distinct(triad_id, node1, node2, node3)
single <- h2$single %>% left_join(node_lookup, by = "triad_id") %>%
  left_join(tr_btw, by = c("triad_id", "t1" = "year"))  # has triad_id, t1, changed12/13/23, emb12/13/23, tie12_btw etc.

single <- single %>%
  get_unbal_load("node1", "node2", "t1", "ul12") %>%
  get_unbal_load("node1", "node3", "t1", "ul13") %>%
  get_unbal_load("node2", "node3", "t1", "ul23") %>%
  mutate(across(c(ul12, ul13, ul23), ~ replace_na(.x, 0)))  # edge present in 0 OTHER triads -> 0 unbalanced load

cat("\nUnbalancedLoad summary (all three ties, all events):\n")
print(summary(c(single$ul12, single$ul13, single$ul23)))

single <- single %>%
  mutate(
    changed_ul = case_when(changed12 ~ ul12, changed13 ~ ul13, changed23 ~ ul23),
    max_ul = pmax(ul12, ul13, ul23),
    n_distinct_ul = pmap_int(list(ul12, ul13, ul23), function(a, b, c) length(unique(c(a, b, c)))),
    is_max_ul = changed_ul == max_ul
  )

cat("\n=== H2: naive test (ties counted as success) ===\n")
cat("Changed tie has MAX UnbalancedLoad in", round(mean(single$is_max_ul) * 100, 1), "% of events (null=33.3%)\n")
print(binom.test(sum(single$is_max_ul), nrow(single), p = 1/3))

cat("\n=== UnbalancedLoad tie structure across the 3 edges ===\n")
cat("All three equal:", round(mean(single$n_distinct_ul == 1) * 100, 1), "%\n")
cat("Exactly two equal:", round(mean(single$n_distinct_ul == 2) * 100, 1), "%\n")
cat("All three distinct:", round(mean(single$n_distinct_ul == 3) * 100, 1), "%\n")

strict <- single %>% filter(n_distinct_ul == 3)
cat("\n=== H2, PRIMARY test: restricted to events with all 3 UnbalancedLoad values distinct (n=", nrow(strict), ") ===\n", sep = "")
cat("Changed tie has MAX UnbalancedLoad in", round(mean(strict$is_max_ul) * 100, 1), "% of events (null=33.3%)\n")
bt <- binom.test(sum(strict$is_max_ul), nrow(strict), p = 1/3)
print(bt)

# --- conditional logit: 3 alternatives per event, predictor = UnbalancedLoad,
#     with and without embeddedness/betweenness as controls ---
long <- single %>%
  transmute(event_id = row_number(), ul12, ul13, ul23, emb12, emb13, emb23, tie12_btw, tie13_btw, tie23_btw,
            changed12, changed13, changed23) %>%
  pivot_longer(cols = c(ul12, ul13, ul23), names_to = "tie", values_to = "unbal_load") %>%
  mutate(tie = sub("ul", "", tie)) %>%
  left_join(
    single %>% mutate(event_id = row_number()) %>%
      pivot_longer(cols = c(emb12, emb13, emb23), names_to = "tie2", values_to = "embeddedness") %>%
      mutate(tie2 = sub("emb", "", tie2)) %>% select(event_id, tie2, embeddedness),
    by = c("event_id", "tie" = "tie2")
  ) %>%
  left_join(
    single %>% mutate(event_id = row_number()) %>%
      pivot_longer(cols = c(tie12_btw, tie13_btw, tie23_btw), names_to = "tie3", values_to = "betweenness") %>%
      mutate(tie3 = sub("tie", "", tie3), tie3 = sub("_btw", "", tie3)) %>% select(event_id, tie3, betweenness),
    by = c("event_id", "tie" = "tie3")
  ) %>%
  left_join(
    single %>% mutate(event_id = row_number()) %>%
      pivot_longer(cols = c(changed12, changed13, changed23), names_to = "tie4", values_to = "changed") %>%
      mutate(tie4 = sub("changed", "", tie4)) %>% select(event_id, tie4, changed),
    by = c("event_id", "tie" = "tie4")
  )

cat("\n=== Conditional logit: P(this tie changed) ~ UnbalancedLoad ===\n")
clog_ul <- clogit(changed ~ unbal_load + strata(event_id), data = long)
print(summary(clog_ul))

cat("\n=== Conditional logit: UnbalancedLoad controlling for embeddedness and betweenness ===\n")
clog_full <- clogit(changed ~ unbal_load + embeddedness + betweenness + strata(event_id), data = long)
print(summary(clog_full))

# --- pressure proportion refinement: UnbalancedLoad / total embeddedness ---
single <- single %>% mutate(
  pp12 = ifelse(emb12 > 0, ul12 / emb12, 0),
  pp13 = ifelse(emb13 > 0, ul13 / emb13, 0),
  pp23 = ifelse(emb23 > 0, ul23 / emb23, 0)
)
strict_pp <- single %>% mutate(
  changed_pp = case_when(changed12 ~ pp12, changed13 ~ pp13, changed23 ~ pp23),
  max_pp = pmax(pp12, pp13, pp23),
  n_distinct_pp = pmap_int(list(round(pp12,6), round(pp13,6), round(pp23,6)), function(a,b,c) length(unique(c(a,b,c)))),
  is_max_pp = abs(changed_pp - max_pp) < 1e-9
) %>% filter(n_distinct_pp == 3)
cat("\n=== Refinement: Pressure Proportion (UnbalancedLoad / total embeddedness), n=", nrow(strict_pp), " ===\n", sep = "")
cat("Changed tie has MAX pressure proportion in", round(mean(strict_pp$is_max_pp) * 100, 1), "% of events (null=33.3%)\n")
bt_pp <- binom.test(sum(strict_pp$is_max_pp), nrow(strict_pp), p = 1/3)
print(bt_pp)

sink("results/h2_unbalanced_load_output.txt")
cat("Realignment events (single tie changed):", nrow(single), "\n\n")
cat("=== Naive test (ties counted as success) ===\n")
cat(round(mean(single$is_max_ul)*100,1), "% (null=33.3%)\n\n")
cat("UnbalancedLoad tie structure: all equal", round(mean(single$n_distinct_ul==1)*100,1),
    "% | two equal", round(mean(single$n_distinct_ul==2)*100,1),
    "% | all distinct", round(mean(single$n_distinct_ul==3)*100,1), "%\n\n")
cat("=== PRIMARY test (all 3 distinct, n=", nrow(strict), ") ===\n", sep=""); print(bt)
cat("\n=== Conditional logit, UnbalancedLoad only ===\n"); print(summary(clog_ul))
cat("\n=== Conditional logit, UnbalancedLoad + embeddedness + betweenness ===\n"); print(summary(clog_full))
cat("\n=== Pressure Proportion refinement (n=", nrow(strict_pp), ") ===\n", sep=""); print(bt_pp)
sink()

saveRDS(list(single = single, strict = strict, bt = bt, clog_ul = clog_ul, clog_full = clog_full,
             strict_pp = strict_pp, bt_pp = bt_pp), "results/h2_unbalanced_load.rds")
cat("\n12_h2_unbalanced_load.R complete. See results/h2_unbalanced_load_output.txt\n")
