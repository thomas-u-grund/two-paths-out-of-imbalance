# ------------------------------------------------------------------
# 17: Cheap proxy check for a two-mode (state x organization) view of
#     alliance structure, without redefining balance itself.
#
#     The one-mode projection used throughout the paper treats a
#     multilateral alliance with N members as a complete dyadic clique:
#     every pair of members gets a full-strength alliance tie, and a
#     tie's "embeddedness" is the number of triangles it anchors in that
#     projected graph. Appendix A0b showed H1 survives this, but the
#     representation itself still overstates how much independent
#     structural weight a huge, diffuse organization like NATO or the
#     Rio Pact contributes relative to a genuine bilateral treaty.
#
#     A full two-mode treatment (states and organizations as distinct
#     node types, closed triads redefined as state-org-state paths)
#     would fix this at the representational level, but requires
#     re-deriving every measure in the paper on a new unit of analysis
#     -- out of scope here. This script instead asks a narrower,
#     answerable question: if a tie's "embeddedness" is measured as its
#     WEIGHTED count of shared organizational memberships -- where a
#     bilateral treaty (2 members) counts fully and an N-member
#     multilateral pact counts as 1/(N-1) per member, so a huge alliance
#     contributes proportionally less "thin" structural weight per
#     co-membership rather than being dyadized at full strength -- does
#     the H1 pattern (embeddedness reduces dissolution, increases
#     realignment) still hold?
# ------------------------------------------------------------------
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(sandwich); library(lmtest) })

# Run from the replication pack root, e.g.:
#   cd replication_pack && Rscript scripts/17_twomode_proxy.R
# (relative paths below assume that working directory; no setwd() needed)

mem <- read_csv("raw/alliances/version4.1_csv/alliance_v4.1_by_member.csv", show_col_types = FALSE) %>%
  mutate(ccode = as.character(ccode),
         mem_end_year = ifelse(is.na(mem_end_year) | mem_end_year == 0, 9999, mem_end_year))

YEARS <- 1816:2012

# --- for every (version4id, year), get the active member set and its size ---
org_year <- mem %>%
  select(version4id, ccode, mem_st_year, mem_end_year) %>%
  distinct()

cat("Computing organization-weighted co-membership for every year... (this enumerates member pairs per org-year, not all states, so it stays fast)\n")

edge_weights <- lapply(YEARS, function(y) {
  active <- org_year %>% filter(mem_st_year <= y, mem_end_year >= y)
  if (nrow(active) == 0) return(NULL)
  active <- active %>% group_by(version4id) %>% filter(n() >= 2) %>% ungroup()
  if (nrow(active) == 0) return(NULL)
  active %>%
    group_by(version4id) %>%
    group_modify(~ {
      states <- sort(unique(.x$ccode))
      n <- length(states)
      if (n < 2) return(tibble())
      pairs <- as.data.frame(t(combn(states, 2)))
      names(pairs) <- c("ccode_lo", "ccode_hi")
      pairs %>% mutate(w = 1 / (n - 1))
    }) %>% ungroup() %>%
    mutate(year = y) %>%
    group_by(year, ccode_lo, ccode_hi) %>%
    summarise(org_emb = sum(w), .groups = "drop")
}) %>% bind_rows()

cat("Edge-year rows with org-weighted embeddedness:", nrow(edge_weights), "\n")
write_csv(edge_weights, "data/org_weighted_embeddedness.csv")

# --- attach to every triad-year's three edges ---
tr <- readRDS("data/triads_all_years.rds") %>%
  select(year, triad_id, node1, node2, node3, balanced, tie_emb_mean, tie_btw_mean,
         actor_load_mean, triad_log_cinc_mean)

get_org_emb <- function(df, a, b, newcol) {
  df %>% mutate(lo = pmin(.data[[a]], .data[[b]]), hi = pmax(.data[[a]], .data[[b]])) %>%
    left_join(edge_weights, by = c("year" = "year", "lo" = "ccode_lo", "hi" = "ccode_hi")) %>%
    rename(!!newcol := org_emb) %>%
    mutate(!!newcol := replace_na(!!sym(newcol), 0)) %>%
    select(-lo, -hi)
}

tr2 <- tr %>%
  get_org_emb("node1", "node2", "oe12") %>%
  get_org_emb("node1", "node3", "oe13") %>%
  get_org_emb("node2", "node3", "oe23") %>%
  rowwise() %>%
  mutate(tie_orgemb_mean = mean(c(oe12, oe13, oe23))) %>%
  ungroup()

cat("\ncor(tie_emb_mean [triangle count on projected graph], tie_orgemb_mean [org-weighted]):",
    cor(tr2$tie_emb_mean, tr2$tie_orgemb_mean, use = "complete.obs"), "\n")
cat("\ntie_orgemb_mean summary:\n"); print(summary(tr2$tie_orgemb_mean))

saveRDS(tr2, "data/triads_with_org_emb.rds")

# ------------------------------------------------------------------
# Rebuild the person-year risk table with the alternative measure and
# refit the primary hazard models
# ------------------------------------------------------------------
tr3 <- tr2 %>% arrange(triad_id, year) %>%
  group_by(triad_id) %>%
  mutate(prev_year = lag(year), prev_balanced = lag(balanced),
         is_onset = !balanced & (is.na(prev_year) | prev_year != year - 1 | prev_balanced),
         spell_id_local = cumsum(is_onset)) %>% ungroup() %>%
  mutate(spell_id = paste(triad_id, spell_id_local, sep = "_S"))

spells <- tr3 %>% filter(!balanced)
LAST_YEAR <- max(tr3$year)
spell_bounds <- spells %>% group_by(spell_id, triad_id) %>% summarise(t1 = max(year), .groups = "drop")
next_status <- tr3 %>% select(year, triad_id, balanced) %>% rename(next_year = year, balanced_next = balanced)
spell_bounds <- spell_bounds %>%
  mutate(lookup_year = t1 + 1) %>%
  left_join(next_status, by = c("triad_id" = "triad_id", "lookup_year" = "next_year")) %>%
  mutate(event = case_when(t1 >= LAST_YEAR ~ "censored_end_of_panel",
                            is.na(balanced_next) ~ "dissolved",
                            balanced_next ~ "flipped",
                            TRUE ~ "ERROR"))
stopifnot(!any(spell_bounds$event == "ERROR"))

py <- spells %>% left_join(spell_bounds %>% select(spell_id, t1, event), by = "spell_id") %>%
  mutate(is_terminal_row = year == t1, row_event = ifelse(is_terminal_row, event, "continues")) %>%
  filter(row_event != "censored_end_of_panel") %>%
  left_join(spells %>% group_by(spell_id) %>% summarise(onset_year = min(year), .groups = "drop"), by = "spell_id") %>%
  mutate(duration = year - onset_year + 1,
         dur_bin = cut(duration, breaks = c(0, 1, 2, 3, 5, 10, Inf),
                        labels = c("yr1", "yr2", "yr3", "yr4-5", "yr6-10", "yr11+")))

zscore <- function(x) as.numeric(scale(x))
py <- py %>% mutate(
  z_orgemb = zscore(tie_orgemb_mean),
  z_tie_emb = zscore(tie_emb_mean),
  z_tie_btw = zscore(tie_btw_mean),
  z_actor_load = zscore(actor_load_mean),
  z_cinc = zscore(triad_log_cinc_mean),
  era = cut(year, breaks = c(1815, 1945, 1989, 2015), labels = c("pre-1945", "1945-1989", "post-1989")),
  event_dissolve = as.numeric(row_event == "dissolved"),
  event_flip = as.numeric(row_event == "flipped")
)

cat("\nProxy-measure risk table:", nrow(py), "triad-years,", n_distinct(py$spell_id), "spells\n")

m_dis_proxy <- glm(event_dissolve ~ z_tie_btw + z_orgemb + z_actor_load + z_cinc + era + dur_bin, data = py, family = binomial())
m_flip_proxy <- glm(event_flip ~ z_tie_btw + z_orgemb + z_actor_load + z_cinc + era + dur_bin, data = py, family = binomial())
ct_dis_proxy <- coeftest(m_dis_proxy, vcov = vcovCL(m_dis_proxy, cluster = py$triad_id))
ct_flip_proxy <- coeftest(m_flip_proxy, vcov = vcovCL(m_flip_proxy, cluster = py$triad_id))

# horse race: both measures in the same model
m_dis_both <- glm(event_dissolve ~ z_tie_btw + z_tie_emb + z_orgemb + z_actor_load + z_cinc + era + dur_bin, data = py, family = binomial())
m_flip_both <- glm(event_flip ~ z_tie_btw + z_tie_emb + z_orgemb + z_actor_load + z_cinc + era + dur_bin, data = py, family = binomial())
ct_dis_both <- coeftest(m_dis_both, vcov = vcovCL(m_dis_both, cluster = py$triad_id))
ct_flip_both <- coeftest(m_flip_both, vcov = vcovCL(m_flip_both, cluster = py$triad_id))

sink("results/twomode_proxy_output.txt")
cat("cor(tie_emb_mean, tie_orgemb_mean):", cor(tr2$tie_emb_mean, tr2$tie_orgemb_mean, use = "complete.obs"), "\n\n")
cat("=== DISSOLUTION, org-weighted embeddedness ONLY ===\n"); print(ct_dis_proxy); cat("N =", nobs(m_dis_proxy), "\n")
cat("\n=== REALIGNMENT, org-weighted embeddedness ONLY ===\n"); print(ct_flip_proxy); cat("N =", nobs(m_flip_proxy), "\n")
cat("\n\n=== DISSOLUTION, both measures (horse race) ===\n"); print(ct_dis_both)
cat("\n=== REALIGNMENT, both measures (horse race) ===\n"); print(ct_flip_both)
sink()
cat(readLines("results/twomode_proxy_output.txt"), sep = "\n")

saveRDS(list(edge_weights = edge_weights, tr2 = tr2, py = py,
             m_dis_proxy = m_dis_proxy, m_flip_proxy = m_flip_proxy,
             ct_dis_proxy = ct_dis_proxy, ct_flip_proxy = ct_flip_proxy,
             m_dis_both = m_dis_both, m_flip_both = m_flip_both,
             ct_dis_both = ct_dis_both, ct_flip_both = ct_flip_both),
        "results/twomode_proxy.rds")
cat("\n17_twomode_proxy.R complete. See results/twomode_proxy_output.txt\n")
