# ------------------------------------------------------------------
# 15: Diagnose the source of the pre-1945 volatility in Figure A1
#     (share of closed triads unbalanced, by year). The appendix had
#     attributed this to generic "small annual denominators." That
#     explanation turns out to be wrong, or at least badly incomplete:
#     the real cause is that COW's alliance dataset codes multilateral
#     treaties as a COMPLETE dyadic clique among all signatories -- a
#     single N-member pact contributes C(N,2) alliance dyads and C(N,3)
#     closed triads simultaneously, almost all of them balanced (since
#     the pact itself is all-positive). A handful of large multilateral
#     instruments therefore dominate the pre-1945 closed-triad count and
#     their formation/dissolution produces mechanical step-changes in
#     the year-by-year series that have nothing to do with a genuine
#     shift in how much triadic imbalance existed that year.
#
#     This script (a) identifies the specific instruments responsible,
#     (b) quantifies how much of the pre-1945 population they account
#     for, (c) checks whether this contaminates the spell-based hazard
#     analysis the paper's actual results rest on (it does not, to a
#     reassuring degree), and (d) reruns the primary hazard models
#     excluding every clique-touching spell as an explicit robustness
#     check.
# ------------------------------------------------------------------
suppressMessages({ library(dplyr); library(readr); library(sandwich); library(lmtest) })

# Run from the replication pack root, e.g.:
#   cd replication_pack && Rscript scripts/15_clique_diagnostic.R
# (relative paths below assume that working directory; no setwd() needed)

mem <- read_csv("raw/alliances/version4.1_csv/alliance_v4.1_by_member.csv", show_col_types = FALSE)
dyad_yr <- read_csv("raw/alliances/version4.1_csv/alliance_v4.1_by_dyad_yearly.csv", show_col_types = FALSE)

# --- step 1: which multilateral instruments have >=8 members? ---
big_pacts <- mem %>% group_by(version4id, all_st_year) %>%
  summarise(n_members = n_distinct(ccode), .groups = "drop") %>%
  filter(n_members >= 8, all_st_year < 1945) %>% arrange(all_st_year)
cat("=== Multilateral alliance instruments with >=8 members, pre-1945 ===\n")
print(big_pacts)

for (id in big_pacts$version4id) {
  sub <- mem %>% filter(version4id == id)
  cat("\nversion4id", id, "| years", unique(sub$all_st_year), "-", unique(sub$all_end_year),
      "| type: defense=", unique(sub$defense), "neutrality=", unique(sub$neutrality),
      "nonaggression=", unique(sub$nonaggression), "entente=", unique(sub$entente), "\n")
  cat("Members:", paste(sort(unique(sub$state_name)), collapse = ", "), "\n")
}

# --- step 2: decompose every pre-1945 closed triad into "inside one of
#     these cliques" vs "other", year by year ---
tr <- readRDS("data/triads_all_years.rds")
CLIQUE_IDS <- big_pacts$version4id  # 3 (German Confederation), 165 (Buenos Aires), 169 (Nyon)

get_members_at <- function(id, yr) {
  sub <- mem %>% filter(version4id == id, mem_st_year <= yr,
                         (mem_end_year >= yr | is.na(mem_end_year) | mem_end_year == 0))
  unique(sub$ccode)
}

years <- 1816:1944
decomp <- lapply(years, function(y) {
  ty <- tr %>% filter(year == y)
  if (nrow(ty) == 0) return(tibble(year = y, total = 0, unbalanced = 0, clique_total = 0, clique_unbalanced = 0))
  in_clique <- Reduce(`|`, lapply(CLIQUE_IDS, function(id) {
    m <- get_members_at(id, y)
    rowSums(cbind(ty$node1 %in% m, ty$node2 %in% m, ty$node3 %in% m)) == 3
  }))
  tibble(year = y, total = nrow(ty), unbalanced = sum(!ty$balanced),
         clique_total = sum(in_clique), clique_unbalanced = sum(in_clique & !ty$balanced))
}) %>% bind_rows()

cat("\n\n=== Pre-1945 closed-triad population: clique-origin vs. other ===\n")
cat("Total closed triad-years, 1816-1944:", sum(decomp$total), "\n")
cat("Inside one of the", length(CLIQUE_IDS), "big multilateral cliques:", sum(decomp$clique_total),
    sprintf("(%.1f%%)\n", 100 * sum(decomp$clique_total) / sum(decomp$total)))
cat("Total unbalanced triad-years, 1816-1944:", sum(decomp$unbalanced), "\n")
cat("Of which inside a clique:", sum(decomp$clique_unbalanced),
    sprintf("(%.1f%%)\n", 100 * sum(decomp$clique_unbalanced) / sum(decomp$unbalanced)))

write_csv(decomp, "results/clique_decomposition_by_year.csv")

# --- step 3: does this contaminate the SPELL-based hazard analysis? ---
tr2 <- tr %>% arrange(triad_id, year)
years_present <- sort(unique(tr2$year))
tag <- lapply(years_present, function(y) {
  in_clique <- Reduce(`|`, lapply(CLIQUE_IDS, function(id) {
    m <- get_members_at(id, y)
    ty <- tr2 %>% filter(year == y)
    rowSums(cbind(ty$node1 %in% m, ty$node2 %in% m, ty$node3 %in% m)) == 3
  }))
  ty <- tr2 %>% filter(year == y)
  tibble(triad_id = ty$triad_id, year = y, in_clique = in_clique)
}) %>% bind_rows()

tr3 <- tr2 %>% left_join(tag, by = c("triad_id", "year")) %>%
  group_by(triad_id) %>%
  mutate(prev_year = lag(year), prev_balanced = lag(balanced),
         is_onset = !balanced & (is.na(prev_year) | prev_year != year - 1 | prev_balanced),
         spell_id_local = cumsum(is_onset)) %>% ungroup() %>%
  mutate(spell_id = paste(triad_id, spell_id_local, sep = "_S"))

spell_summary <- tr3 %>% filter(!balanced) %>% group_by(spell_id) %>%
  summarise(onset_year = min(year), ever_clique = any(in_clique), .groups = "drop")

cat("\n\n=== Does this reach the spell-based hazard analysis? ===\n")
cat("Total spells in the primary analysis:", nrow(spell_summary), "\n")
cat("Spells that ever touch one of the 3 cliques:", sum(spell_summary$ever_clique),
    sprintf("(%.1f%% of all spells)\n", 100 * mean(spell_summary$ever_clique)))
cat("These are mechanically rare because a multilateral peace pact is almost entirely\n")
cat("all-positive, so it only rarely becomes unbalanced in the first place -- a spell\n")
cat("requires a negative tie, which the pact itself does not supply.\n")

clique_spell_ids <- spell_summary %>% filter(ever_clique) %>% pull(spell_id)

# --- step 4: robustness -- rerun the primary hazard models excluding
#     every clique-touching spell entirely ---
py <- readRDS("data/person_year_risk_table_with_duration.rds") %>%
  mutate(dur_bin = cut(duration, breaks = c(0, 1, 2, 3, 5, 10, Inf),
                        labels = c("yr1", "yr2", "yr3", "yr4-5", "yr6-10", "yr11+")))
py_clean <- py %>% filter(!(spell_id %in% clique_spell_ids))
cat("\nTriad-years dropped by excluding clique-touching spells:", nrow(py) - nrow(py_clean), "of", nrow(py), "\n")

m_dis  <- glm(event_dissolve ~ z_tie_btw + z_tie_emb + z_actor_load + z_cinc + era + dur_bin, data = py_clean, family = binomial())
m_flip <- glm(event_flip     ~ z_tie_btw + z_tie_emb + z_actor_load + z_cinc + era + dur_bin, data = py_clean, family = binomial())
ct_dis  <- coeftest(m_dis,  vcov = vcovCL(m_dis,  cluster = py_clean$triad_id))
ct_flip <- coeftest(m_flip, vcov = vcovCL(m_flip, cluster = py_clean$triad_id))

sink("results/clique_diagnostic_output.txt")
cat("=== Big multilateral instruments ===\n"); print(big_pacts)
cat("\n=== Pre-1945 population decomposition ===\n")
cat("Total closed triad-years 1816-1944:", sum(decomp$total), "\n")
cat("Inside a big clique:", sum(decomp$clique_total), sprintf("(%.1f%%)\n", 100*sum(decomp$clique_total)/sum(decomp$total)))
cat("Total unbalanced triad-years 1816-1944:", sum(decomp$unbalanced), "\n")
cat("Inside a big clique:", sum(decomp$clique_unbalanced), sprintf("(%.1f%%)\n", 100*sum(decomp$clique_unbalanced)/sum(decomp$unbalanced)))
cat("\n=== Spell-level exposure ===\n")
cat("Total spells:", nrow(spell_summary), " | clique-touching:", sum(spell_summary$ever_clique), "\n")
cat("\n=== DISSOLUTION hazard, excluding clique-touching spells ===\n"); print(ct_dis)
cat("N =", nobs(m_dis), "\n")
cat("\n=== REALIGNMENT hazard, excluding clique-touching spells ===\n"); print(ct_flip)
cat("N =", nobs(m_flip), "\n")
sink()

saveRDS(list(decomp = decomp, spell_summary = spell_summary, clique_spell_ids = clique_spell_ids,
             m_dis = m_dis, m_flip = m_flip, ct_dis = ct_dis, ct_flip = ct_flip),
        "results/clique_diagnostic.rds")
cat("\n15_clique_diagnostic.R complete. See results/clique_diagnostic_output.txt\n")
