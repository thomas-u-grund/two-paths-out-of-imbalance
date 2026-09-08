# ------------------------------------------------------------------
# 16: Follow-up to script 15. That script showed 3 pre-1945 multilateral
#     alliances dominate the pre-1945 closed-triad population, but only
#     checked those 3 instruments -- it did not check the much larger
#     set of post-1945 multilateral security blocs (NATO, the Warsaw
#     Pact, the Rio Pact/OAS, the Arab League, SEATO, CIS/CSTO, and
#     several regional African pacts). Checking against ALL 19
#     multilateral alliances with >=8 members across the full panel
#     shows a much bigger footprint: 87.6% of triad-years and 86.5% of
#     the 6,277 spells touch at least one of them.
#
#     Unlike the pre-1945 cases, these are not spurious artifacts --
#     NATO and the Warsaw Pact really were the organizing structure of
#     Cold War international relations. But their scale raises two
#     distinct, legitimate questions this script answers directly
#     rather than asserting an answer to:
#     (1) Is the H1 embeddedness effect actually explained by bloc
#         membership, or does it survive independent of it?
#     (2) Are triad-clustered standard errors (the primary specification)
#         understating uncertainty, given many different triads share
#         membership in the same handful of blocs and are not truly
#         independent draws?
# ------------------------------------------------------------------
suppressMessages({ library(dplyr); library(readr); library(sandwich); library(lmtest) })

# Run from the replication pack root, e.g.:
#   cd replication_pack && Rscript scripts/16_bloc_robustness.R
# (relative paths below assume that working directory; no setwd() needed)

mem <- read_csv("raw/alliances/version4.1_csv/alliance_v4.1_by_member.csv", show_col_types = FALSE)
big <- mem %>% group_by(version4id) %>% summarise(n_members = n_distinct(ccode), .groups = "drop") %>% filter(n_members >= 8)
ALL_IDS <- big$version4id
cat("Multilateral alliances with >=8 members, full panel:", length(ALL_IDS), "\n")

tr <- readRDS("data/triads_all_years.rds") %>% arrange(triad_id, year)
get_members_at <- function(id, yr) {
  sub <- mem %>% filter(version4id == id, mem_st_year <= yr, (mem_end_year >= yr | is.na(mem_end_year) | mem_end_year == 0))
  unique(sub$ccode)
}

# --- tag every triad-year: in ANY big bloc, and its PRIMARY bloc id ---
years_present <- sort(unique(tr$year))
tag <- lapply(years_present, function(y) {
  ty <- tr %>% filter(year == y)
  if (nrow(ty) == 0) return(NULL)
  primary_bloc <- rep(NA_integer_, nrow(ty))
  for (id in ALL_IDS) {
    m <- get_members_at(id, y)
    if (length(m) == 0) next
    hit <- rowSums(cbind(ty$node1 %in% m, ty$node2 %in% m, ty$node3 %in% m)) == 3
    primary_bloc[is.na(primary_bloc) & hit] <- id
  }
  tibble(triad_id = ty$triad_id, year = y,
         in_bloc = !is.na(primary_bloc),
         primary_bloc = ifelse(is.na(primary_bloc), "none", as.character(primary_bloc)))
}) %>% bind_rows()

tr2 <- tr %>% left_join(tag, by = c("triad_id", "year"))
cat("\nTriad-years inside ANY of the", length(ALL_IDS), "big blocs:", sum(tr2$in_bloc),
    sprintf("(%.1f%% of %d)\n", 100 * mean(tr2$in_bloc), nrow(tr2)))
cat("Correlation between in_bloc and tie_emb_mean:", cor(tr2$in_bloc, tr2$tie_emb_mean, use = "complete.obs"), "\n")
cat("SD of tie_emb_mean within bloc triads:", sd(tr2$tie_emb_mean[tr2$in_bloc], na.rm = TRUE),
    " | outside:", sd(tr2$tie_emb_mean[!tr2$in_bloc], na.rm = TRUE), "\n")

# --- spell-level exposure ---
tr3 <- tr2 %>% group_by(triad_id) %>%
  mutate(prev_year = lag(year), prev_balanced = lag(balanced),
         is_onset = !balanced & (is.na(prev_year) | prev_year != year - 1 | prev_balanced),
         spell_id_local = cumsum(is_onset)) %>% ungroup() %>%
  mutate(spell_id = paste(triad_id, spell_id_local, sep = "_S"))
spell_summary <- tr3 %>% filter(!balanced) %>% group_by(spell_id) %>%
  summarise(ever_bloc = any(in_bloc), .groups = "drop")
cat("\nSpells ever touching a big bloc:", sum(spell_summary$ever_bloc), "of", nrow(spell_summary),
    sprintf("(%.1f%%)\n", 100 * mean(spell_summary$ever_bloc)))
bloc_spell_ids <- spell_summary %>% filter(ever_bloc) %>% pull(spell_id)

# --- primary risk table with bloc tags ---
py <- readRDS("data/person_year_risk_table_with_duration.rds") %>%
  mutate(dur_bin = cut(duration, breaks = c(0, 1, 2, 3, 5, 10, Inf),
                        labels = c("yr1", "yr2", "yr3", "yr4-5", "yr6-10", "yr11+"))) %>%
  left_join(tag, by = c("triad_id", "year")) %>%
  mutate(in_bloc = ifelse(is.na(in_bloc), FALSE, in_bloc),
         primary_bloc = ifelse(is.na(primary_bloc), "none", primary_bloc))

# --- check 1: exclude every bloc-touching spell entirely ---
py_nonbloc <- py %>% filter(!(spell_id %in% bloc_spell_ids))
cat("\n=== Check 1: excluding ALL bloc-touching spells ===\n")
cat("Triad-years remaining:", nrow(py_nonbloc), "of", nrow(py), "\n")
m_dis_nb <- glm(event_dissolve ~ z_tie_btw + z_tie_emb + z_actor_load + z_cinc + era + dur_bin, data = py_nonbloc, family = binomial())
m_flip_nb <- glm(event_flip ~ z_tie_btw + z_tie_emb + z_actor_load + z_cinc + era + dur_bin, data = py_nonbloc, family = binomial())
ct_dis_nb <- coeftest(m_dis_nb, vcov = vcovCL(m_dis_nb, cluster = py_nonbloc$triad_id))
ct_flip_nb <- coeftest(m_flip_nb, vcov = vcovCL(m_flip_nb, cluster = py_nonbloc$triad_id))

# --- check 2: bloc membership as an additive control ---
m_dis_ctrl <- glm(event_dissolve ~ z_tie_btw + z_tie_emb + z_actor_load + z_cinc + era + dur_bin + in_bloc, data = py, family = binomial())
m_flip_ctrl <- glm(event_flip ~ z_tie_btw + z_tie_emb + z_actor_load + z_cinc + era + dur_bin + in_bloc, data = py, family = binomial())
ct_dis_ctrl <- coeftest(m_dis_ctrl, vcov = vcovCL(m_dis_ctrl, cluster = py$triad_id))
ct_flip_ctrl <- coeftest(m_flip_ctrl, vcov = vcovCL(m_flip_ctrl, cluster = py$triad_id))

# --- check 3: embeddedness x in_bloc interaction ---
m_dis_int <- glm(event_dissolve ~ z_tie_btw + z_tie_emb * in_bloc + z_actor_load + z_cinc + era + dur_bin, data = py, family = binomial())
m_flip_int <- glm(event_flip ~ z_tie_btw + z_tie_emb * in_bloc + z_actor_load + z_cinc + era + dur_bin, data = py, family = binomial())
ct_dis_int <- coeftest(m_dis_int, vcov = vcovCL(m_dis_int, cluster = py$triad_id))
ct_flip_int <- coeftest(m_flip_int, vcov = vcovCL(m_flip_int, cluster = py$triad_id))

# --- check 4: two-way clustering (triad_id + primary_bloc) on the primary spec ---
m_dis_primary <- glm(event_dissolve ~ z_tie_btw + z_tie_emb + z_actor_load + z_cinc + era + dur_bin, data = py, family = binomial())
m_flip_primary <- glm(event_flip ~ z_tie_btw + z_tie_emb + z_actor_load + z_cinc + era + dur_bin, data = py, family = binomial())
cl2 <- py[, c("triad_id", "primary_bloc")]
ct_dis_1way <- coeftest(m_dis_primary, vcov = vcovCL(m_dis_primary, cluster = py$triad_id))
ct_dis_2way <- coeftest(m_dis_primary, vcov = vcovCL(m_dis_primary, cluster = cl2, multi0 = TRUE))
ct_flip_1way <- coeftest(m_flip_primary, vcov = vcovCL(m_flip_primary, cluster = py$triad_id))
ct_flip_2way <- coeftest(m_flip_primary, vcov = vcovCL(m_flip_primary, cluster = cl2, multi0 = TRUE))

sink("results/bloc_robustness_output.txt")
cat("=== Big multilateral alliances (>=8 members) ===\n"); print(big)
cat("\n=== Population concentration ===\n")
cat("Triad-years in a big bloc:", sum(tr2$in_bloc), sprintf("(%.1f%%)\n", 100 * mean(tr2$in_bloc)))
cat("Spells touching a big bloc:", sum(spell_summary$ever_bloc), sprintf("(%.1f%%)\n", 100 * mean(spell_summary$ever_bloc)))
cat("cor(in_bloc, tie_emb_mean):", cor(tr2$in_bloc, tr2$tie_emb_mean, use = "complete.obs"), "\n")

cat("\n\n=== Check 1: DISSOLUTION, excluding all bloc-touching spells (n=", nrow(py_nonbloc), ") ===\n", sep=""); print(ct_dis_nb)
cat("\n=== Check 1: REALIGNMENT, excluding all bloc-touching spells ===\n"); print(ct_flip_nb)

cat("\n\n=== Check 2: DISSOLUTION, with in_bloc as additive control ===\n"); print(ct_dis_ctrl)
cat("\n=== Check 2: REALIGNMENT, with in_bloc as additive control ===\n"); print(ct_flip_ctrl)

cat("\n\n=== Check 3: DISSOLUTION, embeddedness x in_bloc interaction ===\n"); print(ct_dis_int)
cat("\n=== Check 3: REALIGNMENT, embeddedness x in_bloc interaction ===\n"); print(ct_flip_int)

cat("\n\n=== Check 4: DISSOLUTION, z_tie_emb under 1-way vs 2-way clustering ===\n")
cat("1-way (triad only):    b=", coef(m_dis_primary)["z_tie_emb"], " SE=", ct_dis_1way["z_tie_emb","Std. Error"], " p=", ct_dis_1way["z_tie_emb","Pr(>|z|)"], "\n")
cat("2-way (triad + bloc):  b=", coef(m_dis_primary)["z_tie_emb"], " SE=", ct_dis_2way["z_tie_emb","Std. Error"], " p=", ct_dis_2way["z_tie_emb","Pr(>|z|)"], "\n")
cat("\n=== Check 4: REALIGNMENT, z_tie_emb under 1-way vs 2-way clustering ===\n")
cat("1-way (triad only):    b=", coef(m_flip_primary)["z_tie_emb"], " SE=", ct_flip_1way["z_tie_emb","Std. Error"], " p=", ct_flip_1way["z_tie_emb","Pr(>|z|)"], "\n")
cat("2-way (triad + bloc):  b=", coef(m_flip_primary)["z_tie_emb"], " SE=", ct_flip_2way["z_tie_emb","Std. Error"], " p=", ct_flip_2way["z_tie_emb","Pr(>|z|)"], "\n")
sink()

cat(readLines("results/bloc_robustness_output.txt"), sep = "\n")

saveRDS(list(big = big, tr2 = tr2, spell_summary = spell_summary,
             ct_dis_nb = ct_dis_nb, ct_flip_nb = ct_flip_nb,
             ct_dis_ctrl = ct_dis_ctrl, ct_flip_ctrl = ct_flip_ctrl,
             ct_dis_int = ct_dis_int, ct_flip_int = ct_flip_int,
             ct_dis_1way = ct_dis_1way, ct_dis_2way = ct_dis_2way,
             ct_flip_1way = ct_flip_1way, ct_flip_2way = ct_flip_2way),
        "results/bloc_robustness.rds")
cat("\n16_bloc_robustness.R complete. See results/bloc_robustness_output.txt\n")
