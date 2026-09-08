# ------------------------------------------------------------------
# 13: H2 heterogeneity extension. Does unbalanced load predict the
#     changed tie equally well when that tie was an alliance
#     (positive->negative, "alliance collapse") vs. a dispute
#     (negative->positive, "reconciliation")? And does unbalanced load
#     predict the direction of change toward balance specifically?
# ------------------------------------------------------------------
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(survival); library(purrr) })

# Run from the replication pack root, e.g.:
#   cd replication_pack && Rscript scripts/13_h2_heterogeneity.R
# (relative paths below assume that working directory; no setwd() needed)

h2ul <- readRDS("results/h2_unbalanced_load.rds")
single <- h2ul$single

# --- restrict to the primary (all-3-distinct) sample used for the main H2 test ---
strict <- single %>% filter(n_distinct_ul == 3)
cat("Strict (all 3 unbalanced-load values distinct) sample:", nrow(strict), "\n")

strict <- strict %>% mutate(
  changed_type_sign = case_when(changed12 ~ sign12_t1, changed13 ~ sign13_t1, changed23 ~ sign23_t1),
  changed_tie_type = ifelse(changed_type_sign == 1, "alliance (pos->neg)", "dispute (neg->pos)"),
  # triad configuration at t1 (before change): count of positive vs negative ties
  n_pos_t1 = (sign12_t1 == 1) + (sign13_t1 == 1) + (sign23_t1 == 1),
  triad_config = case_when(
    n_pos_t1 == 2 ~ "two_pos_one_neg",
    n_pos_t1 == 0 ~ "zero_pos_three_neg",
    TRUE ~ paste0("other(n_pos=", n_pos_t1, ")")
  )
)

cat("\n=== Triad configuration at t1 (before change) ===\n")
print(table(strict$triad_config))

cat("\n=== Changed tie type ===\n")
print(table(strict$changed_tie_type))

# ------------------------------------------------------------------
# 1. Alliance vs. dispute heterogeneity
# ------------------------------------------------------------------
run_group_test <- function(df, label) {
  n <- nrow(df)
  hits <- sum(df$is_max_ul)
  bt <- binom.test(hits, n, p = 1/3)
  tibble(group = label, n = n, hit_rate = hits / n,
         ci_low = bt$conf.int[1], ci_high = bt$conf.int[2], p_value = bt$p.value)
}

het_results <- bind_rows(
  run_group_test(strict, "All events"),
  run_group_test(strict %>% filter(changed_tie_type == "alliance (pos->neg)"), "Alliance events (pos->neg)"),
  run_group_test(strict %>% filter(changed_tie_type == "dispute (neg->pos)"), "Dispute events (neg->pos)")
)
cat("\n=== Heterogeneity by changed-tie type ===\n")
print(het_results)

# formal test of whether the hit rate differs by group (proportions test)
tab <- strict %>% count(changed_tie_type, is_max_ul) %>% pivot_wider(names_from = is_max_ul, values_from = n, values_fill = 0)
cat("\n2x2 table (rows=tie type, cols=hit):\n"); print(tab)
prop_test <- prop.test(x = c(sum(strict$is_max_ul[strict$changed_tie_type == "alliance (pos->neg)"]),
                              sum(strict$is_max_ul[strict$changed_tie_type == "dispute (neg->pos)"])),
                        n = c(sum(strict$changed_tie_type == "alliance (pos->neg)"),
                              sum(strict$changed_tie_type == "dispute (neg->pos)")))
cat("\n=== Formal test: hit rate, alliance vs. dispute events ===\n")
print(prop_test)

# ------------------------------------------------------------------
# 2. Direction-of-change cross-tab
# ------------------------------------------------------------------
cat("\n=== Cross-tab: triad configuration x changed-tie type x hit rate ===\n")
crosstab <- strict %>% group_by(triad_config, changed_tie_type) %>%
  summarise(n = n(), hit_rate = mean(is_max_ul), .groups = "drop")
print(crosstab)

# ------------------------------------------------------------------
# 3. Conditional logit with unbalanced_load x tie_type interaction
#    (fit on the full single-change sample, all events, matching script 12's
#    primary conditional logit sample -- not restricted to the distinct-UL
#    subsample, consistent with Table 2 in the main text)
# ------------------------------------------------------------------
node_lookup <- readRDS("data/triads_all_years.rds") %>% distinct(triad_id, node1, node2, node3)
full_single <- single  # all 5,414 single-tie-change events

long <- full_single %>%
  transmute(event_id = row_number(), triad_id, ul12, ul13, ul23,
            sign12_t1, sign13_t1, sign23_t1, changed12, changed13, changed23) %>%
  pivot_longer(cols = c(ul12, ul13, ul23), names_to = "tie", values_to = "unbal_load") %>%
  mutate(tie = sub("ul", "", tie)) %>%
  left_join(
    full_single %>% mutate(event_id = row_number()) %>%
      pivot_longer(cols = c(sign12_t1, sign13_t1, sign23_t1), names_to = "tie2", values_to = "sign_t1") %>%
      mutate(tie2 = sub("sign", "", tie2), tie2 = sub("_t1", "", tie2)) %>% select(event_id, tie2, sign_t1),
    by = c("event_id", "tie" = "tie2")
  ) %>%
  left_join(
    full_single %>% mutate(event_id = row_number()) %>%
      pivot_longer(cols = c(changed12, changed13, changed23), names_to = "tie3", values_to = "changed") %>%
      mutate(tie3 = sub("changed", "", tie3)) %>% select(event_id, tie3, changed),
    by = c("event_id", "tie" = "tie3")
  ) %>%
  mutate(tie_type_alliance = as.numeric(sign_t1 == 1))  # 1 = this tie is currently an alliance, 0 = dispute

cat("\nPerson-tie rows for interaction model:", nrow(long), " | events:", n_distinct(long$event_id), "\n")

clog_main <- clogit(changed ~ unbal_load + tie_type_alliance + strata(event_id), data = long)
cat("\n=== Conditional logit: main effects (unbalanced load + tie type) ===\n")
print(summary(clog_main))

clog_int <- clogit(changed ~ unbal_load * tie_type_alliance + strata(event_id), data = long)
cat("\n=== Conditional logit: unbalanced load x tie-type interaction ===\n")
print(summary(clog_int))

lrt <- anova(clog_main, clog_int, test = "Chisq")
cat("\nLikelihood ratio test, main-effects vs. interaction model:\n"); print(lrt)

# ------------------------------------------------------------------
# 4. Descriptive table + save everything
# ------------------------------------------------------------------
write_csv(het_results, "results/h2_heterogeneity_table.csv")
write_csv(crosstab, "results/h2_direction_crosstab.csv")

sink("results/h2_heterogeneity_output.txt")
cat("=== Table: hit rate by group ===\n"); print(het_results)
cat("\n=== 2x2 table ===\n"); print(tab)
cat("\n=== Proportions test, alliance vs dispute ===\n"); print(prop_test)
cat("\n=== Direction/configuration cross-tab ===\n"); print(crosstab)
cat("\n=== Conditional logit: main effects ===\n"); print(summary(clog_main))
cat("\n=== Conditional logit: interaction ===\n"); print(summary(clog_int))
cat("\n=== LRT main vs interaction ===\n"); print(lrt)
sink()

saveRDS(list(strict = strict, het_results = het_results, crosstab = crosstab,
             clog_main = clog_main, clog_int = clog_int, lrt = lrt),
        "results/h2_heterogeneity.rds")
cat("\n13_h2_heterogeneity.R complete. See results/h2_heterogeneity_output.txt\n")
