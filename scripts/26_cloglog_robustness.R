# ------------------------------------------------------------------
# 26: Cloglog robustness check for the discrete-time hazard link
#     (online appendix, Section A15).
#
#     Motivated by external review: the main text's original
#     justification for a logistic link ("hazard probabilities well
#     below where the links diverge") does not hold -- the realignment
#     hazard's overall annual event rate is 53.3%, nowhere near rare.
#     This re-estimates the primary hazard specification (Table 2:
#     triadic embeddedness, tie betweenness, actor constraint,
#     capability, era and spell-duration fixed effects) with a
#     complementary log-log link instead of logit, on the identical
#     person-year risk table, to check whether the link choice matters.
#     It does not: both hazards keep the same sign and remain highly
#     significant under either link.
# ------------------------------------------------------------------
suppressMessages({
  library(dplyr); library(sandwich); library(lmtest)
})

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

event_rate_dissolve <- mean(py$event_dissolve)
event_rate_flip <- mean(py$event_flip)
cat("Overall annual dissolution event rate:", round(event_rate_dissolve, 4), "\n")
cat("Overall annual realignment event rate:", round(event_rate_flip, 4), "\n\n")

m_dis_logit  <- glm(event_dissolve ~ z_tie_btw + z_tie_emb + z_actor_load + z_cinc + era + dur_bin,
                     data = py, family = binomial(link = "logit"))
m_flip_logit <- glm(event_flip     ~ z_tie_btw + z_tie_emb + z_actor_load + z_cinc + era + dur_bin,
                     data = py, family = binomial(link = "logit"))
m_dis_clog   <- glm(event_dissolve ~ z_tie_btw + z_tie_emb + z_actor_load + z_cinc + era + dur_bin,
                     data = py, family = binomial(link = "cloglog"))
m_flip_clog  <- glm(event_flip     ~ z_tie_btw + z_tie_emb + z_actor_load + z_cinc + era + dur_bin,
                     data = py, family = binomial(link = "cloglog"))

ct_dis_logit  <- coeftest(m_dis_logit,  vcov = vcovCL(m_dis_logit,  cluster = py$triad_id))
ct_flip_logit <- coeftest(m_flip_logit, vcov = vcovCL(m_flip_logit, cluster = py$triad_id))
ct_dis_clog   <- coeftest(m_dis_clog,   vcov = vcovCL(m_dis_clog,   cluster = py$triad_id))
ct_flip_clog  <- coeftest(m_flip_clog,  vcov = vcovCL(m_flip_clog,  cluster = py$triad_id))

sink("results/cloglog_robustness_output.txt")
cat("Overall annual dissolution event rate:", round(event_rate_dissolve, 4), "\n")
cat("Overall annual realignment event rate:", round(event_rate_flip, 4), "\n\n")
cat("=== DISSOLUTION: logit (primary) ===\n"); print(ct_dis_logit)
cat("\n=== DISSOLUTION: cloglog ===\n"); print(ct_dis_clog)
cat("\n=== REALIGNMENT: logit (primary) ===\n"); print(ct_flip_logit)
cat("\n=== REALIGNMENT: cloglog ===\n"); print(ct_flip_clog)
sink()
cat(readLines("results/cloglog_robustness_output.txt"), sep = "\n")

saveRDS(list(event_rate_dissolve = event_rate_dissolve, event_rate_flip = event_rate_flip,
             m_dis_logit = m_dis_logit, m_flip_logit = m_flip_logit,
             m_dis_clog = m_dis_clog, m_flip_clog = m_flip_clog,
             ct_dis_logit = ct_dis_logit, ct_flip_logit = ct_flip_logit,
             ct_dis_clog = ct_dis_clog, ct_flip_clog = ct_flip_clog),
        "results/cloglog_robustness.rds")

cat("\n26_cloglog_robustness.R complete. See results/cloglog_robustness_output.txt\n")
