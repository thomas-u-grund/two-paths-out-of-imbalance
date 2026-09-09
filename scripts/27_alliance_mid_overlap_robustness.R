# ------------------------------------------------------------------
# 27: Alliance/MID-overlap robustness check (external review, third
#     pass). The precedence rule fixed in script 01 (dispute wins over
#     an active alliance in the same dyad-year: 870 of 63,941 dyad-
#     years, 1.36%) means "realignment" refers to a change in the
#     dyad's CODED signed state, not necessarily to the underlying
#     treaty or dispute itself starting or ending. Since one dyad
#     participates in many triads, a small share of ambiguous dyad-
#     years could in principle generate a disproportionate share of
#     triad-level events.
#
#     Two checks:
#     (A) H1: drop every dyad-year with an alliance/MID overlap
#         entirely (treat as if the dyad had no signed tie that year,
#         same as any other dyad with neither), rebuild the closed-
#         triad panel and the spell/risk-table structure from scratch
#         on this restricted signed-tie set, and re-estimate the
#         primary (duration-adjusted) hazard models.
#     (B) H2: identify, among the 5,414 single-tie-change realignment
#         events, how many involve a changing dyad that was itself an
#         alliance/MID overlap case at either endpoint of the
#         transition (t1 or t1+1). Report the count, then re-run the
#         primary unbalanced-load test and the same-sign-pairs test
#         excluding those events.
# ------------------------------------------------------------------
suppressMessages({
  library(dplyr); library(readr); library(tidyr); library(igraph)
  library(sandwich); library(lmtest); library(purrr)
})


# ==================================================================
# PART A: H1, alliance/MID-overlap dyad-years dropped entirely
# ==================================================================

signed_ties_full <- read_csv("data/signed_ties.csv", show_col_types = FALSE)
cat("Full signed dyad-years:", nrow(signed_ties_full), "\n")
cat("Multiplex (alliance+MID overlap) dyad-years:", sum(signed_ties_full$multiplex), "\n")

signed_ties <- signed_ties_full %>%
  filter(!multiplex) %>%                      # drop overlap dyad-years entirely
  filter(year >= 1816, year <= 2012) %>%
  select(year, ccode_low, ccode_high, sign)
cat("Signed dyad-years after dropping overlaps:", nrow(signed_ties), "\n\n")

# --- triangulation, identical logic to 02_triads_and_load.R -------
nmc <- read_csv("raw/nmc/NMCv7/abridged/NMC-70-abridged.csv", show_col_types = FALSE) %>%
  select(ccode, year, cinc) %>%
  mutate(log_cinc = log(pmax(cinc, 1e-6)))

years <- sort(unique(signed_ties$year))
cat("Years to process:", min(years), "-", max(years), "(", length(years), "years )\n")

process_year <- function(yr, ties) {
  d <- ties %>% filter(year == yr) %>%
    mutate(ccode1 = as.character(ccode_low), ccode2 = as.character(ccode_high)) %>%
    select(ccode1, ccode2, sign) %>%
    distinct(ccode1, ccode2, .keep_all = TRUE)
  if (nrow(d) < 3) return(NULL)
  g <- graph_from_data_frame(d[, c("ccode1", "ccode2")], directed = FALSE)
  E(g)$sign <- d$sign
  tri <- igraph::triangles(g)
  if (length(tri) == 0) return(NULL)
  tri_mat <- matrix(V(g)$name[tri], ncol = 3, byrow = TRUE)

  deg <- degree(g)
  con <- tryCatch(constraint(g), error = function(e) rep(NA_real_, vcount(g)))
  node_tbl <- tibble(node = V(g)$name, degree = deg, constraint = con)
  get_node <- function(x, col) node_tbl[[col]][match(x, node_tbl$node)]

  el <- as_data_frame(g, what = "edges") %>% mutate(e1 = pmin(from, to), e2 = pmax(from, to))
  edge_btw <- edge_betweenness(g)
  embeddedness <- sapply(seq_len(ecount(g)), function(i) {
    ends_i <- ends(g, i)
    length(intersect(neighbors(g, ends_i[1]), neighbors(g, ends_i[2])))
  })
  edge_tbl <- tibble(e1 = el$e1, e2 = el$e2, sign = E(g)$sign, edge_btw = edge_btw, embeddedness = embeddedness)

  n1 <- tri_mat[, 1]; n2 <- tri_mat[, 2]; n3 <- tri_mat[, 3]
  ids_sorted <- t(apply(tri_mat, 1, sort))
  s1 <- ids_sorted[, 1]; s2 <- ids_sorted[, 2]; s3 <- ids_sorted[, 3]
  idx12 <- match(paste(s1, s2), paste(edge_tbl$e1, edge_tbl$e2))
  idx13 <- match(paste(s1, s3), paste(edge_tbl$e1, edge_tbl$e2))
  idx23 <- match(paste(s2, s3), paste(edge_tbl$e1, edge_tbl$e2))

  sign12 <- edge_tbl$sign[idx12]; sign13 <- edge_tbl$sign[idx13]; sign23 <- edge_tbl$sign[idx23]
  balanced <- (sign12 * sign13 * sign23) > 0

  tibble(
    year = yr,
    node1 = ids_sorted[, 1], node2 = ids_sorted[, 2], node3 = ids_sorted[, 3],
    sign12, sign13, sign23, balanced = balanced,
    actor_load_mean = (get_node(n1, "constraint") + get_node(n2, "constraint") + get_node(n3, "constraint")) / 3,
    tie12_emb = edge_tbl$embeddedness[idx12], tie13_emb = edge_tbl$embeddedness[idx13], tie23_emb = edge_tbl$embeddedness[idx23],
    tie12_btw = edge_tbl$edge_btw[idx12], tie13_btw = edge_tbl$edge_btw[idx13], tie23_btw = edge_tbl$edge_btw[idx23],
    tie_btw_mean = (tie12_btw + tie13_btw + tie23_btw) / 3,
    tie_emb_mean = (tie12_emb + tie13_emb + tie23_emb) / 3,
    triad_id = paste(ids_sorted[, 1], ids_sorted[, 2], ids_sorted[, 3], sep = "-")
  )
}

results <- vector("list", length(years))
t0 <- Sys.time()
for (i in seq_along(years)) {
  results[[i]] <- tryCatch(process_year(years[i], signed_ties), error = function(e) {
    message("YEAR ", years[i], " FAILED: ", conditionMessage(e)); NULL
  })
}
cat("Triangulation elapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")

triads_ov <- bind_rows(results)
cat("Closed triads (overlap-dropped panel):", nrow(triads_ov),
    " (full panel: 378,623)\n")
cat("Unbalanced:", sum(!triads_ov$balanced), " (full panel: 10,350)\n\n")

nmc_lookup <- nmc %>% mutate(ccode = as.character(ccode)) %>% select(ccode, year, log_cinc)
triads_ov <- triads_ov %>%
  left_join(nmc_lookup, by = c("node1" = "ccode", "year")) %>% rename(cinc1 = log_cinc) %>%
  left_join(nmc_lookup, by = c("node2" = "ccode", "year")) %>% rename(cinc2 = log_cinc) %>%
  left_join(nmc_lookup, by = c("node3" = "ccode", "year")) %>% rename(cinc3 = log_cinc) %>%
  rowwise() %>% mutate(triad_log_cinc_mean = mean(c(cinc1, cinc2, cinc3), na.rm = TRUE)) %>%
  ungroup() %>% select(-cinc1, -cinc2, -cinc3)

# --- spell construction, identical logic to 07_competing_risks.R --
tr <- triads_ov %>%
  select(year, triad_id, balanced, actor_load_mean, tie_btw_mean, tie_emb_mean, triad_log_cinc_mean) %>%
  arrange(triad_id, year)
LAST_YEAR <- max(tr$year)

tr <- tr %>%
  group_by(triad_id) %>%
  mutate(prev_year = lag(year), prev_balanced = lag(balanced),
         is_onset = !balanced & (is.na(prev_year) | prev_year != year - 1 | prev_balanced),
         spell_id_local = cumsum(is_onset)) %>%
  ungroup() %>%
  mutate(spell_id = paste(triad_id, spell_id_local, sep = "_S"))

spells <- tr %>% filter(!balanced)
cat("Unbalanced-spell candidate rows (overlap-dropped panel):", nrow(spells), "\n")
cat("Unique spells:", n_distinct(spells$spell_id), "\n")

spell_bounds <- spells %>% group_by(spell_id, triad_id) %>% summarise(t1 = max(year), .groups = "drop")
next_status <- tr %>% select(year, triad_id, balanced) %>% rename(next_year = year, balanced_next = balanced)
spell_bounds <- spell_bounds %>%
  mutate(lookup_year = t1 + 1) %>%
  left_join(next_status, by = c("triad_id" = "triad_id", "lookup_year" = "next_year")) %>%
  mutate(event = case_when(
    t1 >= LAST_YEAR      ~ "censored_end_of_panel",
    is.na(balanced_next) ~ "dissolved",
    balanced_next        ~ "flipped",
    TRUE                 ~ "ERROR_still_unbalanced"
  ))
stopifnot(!any(spell_bounds$event == "ERROR_still_unbalanced"))
cat("Spell terminal events (overlap-dropped panel):\n"); print(table(spell_bounds$event))

person_years <- spells %>%
  left_join(spell_bounds %>% select(spell_id, t1, event), by = "spell_id") %>%
  mutate(is_terminal_row = year == t1, row_event = ifelse(is_terminal_row, event, "continues")) %>%
  filter(row_event != "censored_end_of_panel")
cat("Usable person-year rows (overlap-dropped panel):", nrow(person_years), "\n")
print(table(person_years$row_event))

# --- duration bins, identical logic to 14_duration_dependence.R ---
onset <- spells %>% group_by(spell_id) %>% summarise(onset_year = min(year), .groups = "drop")
zscore <- function(x) as.numeric(scale(x))
person_years <- person_years %>%
  left_join(onset, by = "spell_id") %>%
  mutate(duration = year - onset_year + 1,
         dur_bin = cut(duration, breaks = c(0, 1, 2, 3, 5, 10, Inf),
                        labels = c("yr1", "yr2", "yr3", "yr4-5", "yr6-10", "yr11+")),
         z_actor_load = zscore(actor_load_mean), z_tie_btw = zscore(tie_btw_mean),
         z_tie_emb = zscore(tie_emb_mean), z_cinc = zscore(triad_log_cinc_mean),
         era = cut(year, breaks = c(1815, 1945, 1989, 2015), labels = c("pre-1945", "1945-1989", "post-1989")),
         event_dissolve = as.numeric(row_event == "dissolved"),
         event_flip = as.numeric(row_event == "flipped"))

m_dis_ov  <- glm(event_dissolve ~ z_tie_btw + z_tie_emb + z_actor_load + z_cinc + era + dur_bin,
                  data = person_years, family = binomial())
m_flip_ov <- glm(event_flip     ~ z_tie_btw + z_tie_emb + z_actor_load + z_cinc + era + dur_bin,
                  data = person_years, family = binomial())
ct_dis_ov  <- coeftest(m_dis_ov,  vcov = vcovCL(m_dis_ov,  cluster = person_years$triad_id))
ct_flip_ov <- coeftest(m_flip_ov, vcov = vcovCL(m_flip_ov, cluster = person_years$triad_id))

orig <- readRDS("results/duration_dependence_models.rds")

cat("\n=== H1: primary specification, alliance/MID-overlap dyad-years dropped entirely ===\n")
cat("N (person-years):", nobs(m_dis_ov), "(full panel: 10,170)\n\n")
cat("--- Dissolution hazard ---\n")
print(ct_dis_ov["z_tie_emb", ])
cat("(full panel, Table 2): b=-2.508, SE=0.161\n\n")
cat("--- Realignment hazard ---\n")
print(ct_flip_ov["z_tie_emb", ])
cat("(full panel, Table 2): b=+0.439, SE=0.046\n")

# ==================================================================
# PART B: H2, flag and exclude overlap-affected realignment events
# ==================================================================
cat("\n\n=== H2: how many single-tie-change events involve an overlap-affected changing dyad? ===\n")

h2 <- readRDS("results/h2_unbalanced_load.rds")
single <- h2$single  # 5,414 single-tie-change realignment events (original, full panel)

overlap_lookup <- signed_ties_full %>%
  filter(multiplex) %>%
  transmute(year, lo = pmin(ccode_low, ccode_high), hi = pmax(ccode_low, ccode_high)) %>%
  mutate(key = paste(year, lo, hi))

is_overlap <- function(node_a, node_b, yr) {
  lo <- pmin(as.numeric(node_a), as.numeric(node_b))
  hi <- pmax(as.numeric(node_a), as.numeric(node_b))
  paste(yr, lo, hi) %in% overlap_lookup$key
}

single <- single %>%
  mutate(
    changed_n1 = case_when(changed12 ~ node1, changed13 ~ node1, changed23 ~ node2),
    changed_n2 = case_when(changed12 ~ node2, changed13 ~ node3, changed23 ~ node3),
    overlap_at_t1 = is_overlap(changed_n1, changed_n2, t1),
    overlap_at_t2 = is_overlap(changed_n1, changed_n2, lookup_year),
    overlap_affected = overlap_at_t1 | overlap_at_t2
  )

n_affected <- sum(single$overlap_affected)
cat("Single-tie-change events:", nrow(single), "\n")
cat("Of these, the changing dyad was itself an alliance/MID overlap case at t1 or t1+1 in:",
    n_affected, sprintf("(%.2f%%)\n", 100 * n_affected / nrow(single)))

clean <- single %>% filter(!overlap_affected)
cat("\nRemaining after excluding overlap-affected events:", nrow(clean), "\n")

# --- primary unbalanced-load test, re-run on the clean subsample ---
strict_clean <- clean %>% filter(n_distinct_ul == 3)
cat("\n--- Primary test (all 3 UnbalancedLoad values distinct), overlap-affected events excluded ---\n")
cat("N =", nrow(strict_clean), "(original strict sample: 2,516)\n")
cat("Changed tie has MAX UnbalancedLoad in", round(mean(strict_clean$is_max_ul) * 100, 1),
    "% of events (original: 96.6%; null = 33.3%)\n")
print(binom.test(sum(strict_clean$is_max_ul), nrow(strict_clean), p = 1/3))

# --- same-sign-pairs test, re-run on the clean subsample ---
long_clean <- clean %>%
  mutate(event_id = row_number()) %>%
  select(event_id, ul12, ul13, ul23, sign12_t1, sign13_t1, sign23_t1, changed12, changed13, changed23) %>%
  pivot_longer(cols = c(ul12, ul13, ul23), names_to = "tie", values_to = "unbal_load") %>%
  mutate(tie = sub("ul", "", tie)) %>%
  left_join(
    clean %>% mutate(event_id = row_number()) %>%
      select(event_id, sign12_t1, sign13_t1, sign23_t1) %>%
      pivot_longer(cols = c(sign12_t1, sign13_t1, sign23_t1), names_to = "tie2", values_to = "sign") %>%
      mutate(tie2 = sub("sign", "", tie2), tie2 = sub("_t1", "", tie2)) %>% select(event_id, tie2, sign),
    by = c("event_id", "tie" = "tie2")
  ) %>%
  left_join(
    clean %>% mutate(event_id = row_number()) %>%
      select(event_id, changed12, changed13, changed23) %>%
      pivot_longer(cols = c(changed12, changed13, changed23), names_to = "tie3", values_to = "changed") %>%
      mutate(tie3 = sub("changed", "", tie3)) %>% select(event_id, tie3, changed),
    by = c("event_id", "tie" = "tie3")
  )

same_sign_pairs_clean <- long_clean %>%
  group_by(event_id) %>%
  filter(n() == 3) %>%
  group_modify(~ {
    df <- .x
    pairs <- combn(1:3, 2, simplify = FALSE)
    map_dfr(pairs, function(idx) {
      a <- df[idx[1], ]; b <- df[idx[2], ]
      if (a$sign == b$sign && a$changed != b$changed) {
        tibble(higher_load_changed = if (a$unbal_load == b$unbal_load) NA
               else (a$unbal_load > b$unbal_load) == a$changed)
      } else tibble(higher_load_changed = NA_real_)
    })
  }) %>%
  ungroup() %>% filter(!is.na(higher_load_changed))

cat("\n--- Same-sign-pairs test, overlap-affected events excluded ---\n")
cat("N comparable pairs =", nrow(same_sign_pairs_clean), "(original: 137)\n")
cat("Higher-unbalanced-load tie changed in", round(mean(same_sign_pairs_clean$higher_load_changed) * 100, 1),
    "% of cases (original: 75.2%; null = 50%)\n")
print(binom.test(sum(same_sign_pairs_clean$higher_load_changed), nrow(same_sign_pairs_clean), p = 0.5))

sink("results/alliance_mid_overlap_robustness_output.txt")
cat("=== PART A: H1, alliance/MID-overlap dyad-years dropped ===\n")
cat("N (person-years):", nobs(m_dis_ov), " | full panel: 10,170\n\n")
cat("Dissolution hazard, z_tie_emb:\n"); print(ct_dis_ov["z_tie_emb", ]); cat("full panel: b=-2.508, SE=0.161\n\n")
cat("Realignment hazard, z_tie_emb:\n"); print(ct_flip_ov["z_tie_emb", ]); cat("full panel: b=+0.439, SE=0.046\n\n")
cat("Full dissolution model:\n"); print(ct_dis_ov)
cat("\nFull realignment model:\n"); print(ct_flip_ov)

cat("\n\n=== PART B: H2, overlap-affected realignment events ===\n")
cat("Single-tie-change events:", nrow(single), "\n")
cat("Overlap-affected:", n_affected, sprintf("(%.2f%%)\n", 100 * n_affected / nrow(single)))
cat("\nPrimary test (all 3 UL distinct), excluding overlap-affected: N=", nrow(strict_clean),
    ", hit rate=", round(mean(strict_clean$is_max_ul)*100,1), "% (original 96.6%)\n", sep = "")
print(binom.test(sum(strict_clean$is_max_ul), nrow(strict_clean), p = 1/3))
cat("\nSame-sign-pairs test, excluding overlap-affected: N=", nrow(same_sign_pairs_clean),
    ", hit rate=", round(mean(same_sign_pairs_clean$higher_load_changed)*100,1), "% (original 75.2%)\n", sep = "")
print(binom.test(sum(same_sign_pairs_clean$higher_load_changed), nrow(same_sign_pairs_clean), p = 0.5))
sink()
cat(readLines("results/alliance_mid_overlap_robustness_output.txt"), sep = "\n")

saveRDS(list(triads_ov = triads_ov, person_years = person_years,
             m_dis_ov = m_dis_ov, m_flip_ov = m_flip_ov, ct_dis_ov = ct_dis_ov, ct_flip_ov = ct_flip_ov,
             single = single, n_affected = n_affected, strict_clean = strict_clean,
             same_sign_pairs_clean = same_sign_pairs_clean),
        "results/alliance_mid_overlap_robustness.rds")
cat("\n27_alliance_mid_overlap_robustness.R complete.\n")
