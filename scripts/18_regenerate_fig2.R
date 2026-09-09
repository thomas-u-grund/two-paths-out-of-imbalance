# ------------------------------------------------------------------
# 18: Generate Figure 3 (results/fig3_hazard_coefficients.png):
#     coefficient plot for the primary duration-adjusted hazard models
#     (script 14) -- dissolution and realignment side by side, greyscale
#     and shape-coded rather than color-coded, with variable labels and
#     terminology matched to Table 1 and the main text (tie betweenness,
#     triadic embeddedness, actor constraint, triad-years).
# ------------------------------------------------------------------
suppressMessages({ library(dplyr); library(ggplot2); library(sandwich); library(lmtest) })

# Run from the replication pack root, e.g.:
#   cd replication_pack && Rscript scripts/18_regenerate_fig2.R
# (relative paths below assume that working directory; no setwd() needed)

py <- readRDS("data/person_year_risk_table_with_duration.rds") %>%
  mutate(dur_bin = cut(duration, breaks = c(0, 1, 2, 3, 5, 10, Inf),
                        labels = c("yr1", "yr2", "yr3", "yr4-5", "yr6-10", "yr11+")))

m_dis  <- glm(event_dissolve ~ z_tie_btw + z_tie_emb + z_actor_load + z_cinc + era + dur_bin, data = py, family = binomial())
m_flip <- glm(event_flip     ~ z_tie_btw + z_tie_emb + z_actor_load + z_cinc + era + dur_bin, data = py, family = binomial())
ct_dis  <- coeftest(m_dis,  vcov = vcovCL(m_dis,  cluster = py$triad_id))
ct_flip <- coeftest(m_flip, vcov = vcovCL(m_flip, cluster = py$triad_id))

cat("Sanity check against Table 1 (should match to 3 decimals):\n")
cat("Dissolution z_tie_emb:", round(ct_dis["z_tie_emb", "Estimate"], 3), "(Table 1: -2.508)\n")
cat("Realignment z_tie_emb:", round(ct_flip["z_tie_emb", "Estimate"], 3), "(Table 1: 0.439)\n")

build_coef_df <- function(ct, label) {
  df <- as.data.frame(unclass(ct))
  colnames(df) <- c("estimate", "se", "z", "p")
  df %>%
    tibble::rownames_to_column("term") %>%
    filter(term %in% c("z_tie_btw", "z_tie_emb", "z_actor_load", "z_cinc", "era1945-1989", "erapost-1989")) %>%
    mutate(contrast = label,
           label = recode(term,
             z_tie_btw = "Tie betweenness",
             z_tie_emb = "Triadic embeddedness",
             z_actor_load = "Actor constraint",
             z_cinc = "Capability (log CINC)",
             `era1945-1989` = "Era: 1945-1989", `erapost-1989` = "Era: post-1989"),
           lo = estimate - 1.96 * se, hi = estimate + 1.96 * se)
}

coef_df <- bind_rows(
  build_coef_df(ct_dis, "Dissolution hazard"),
  build_coef_df(ct_flip, "Realignment hazard")
)

p <- ggplot(coef_df, aes(x = estimate, y = reorder(label, estimate), shape = contrast)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.15, position = position_dodge(width = 0.5), color = "grey20") +
  geom_point(size = 4, position = position_dodge(width = 0.5), color = "grey10", fill = "white", stroke = 1.1) +
  scale_shape_manual(values = c("Dissolution hazard" = 16, "Realignment hazard" = 24)) +
  labs(x = "Cause-specific hazard model coefficient (95% CI, clustered SE by triad)", y = NULL, shape = NULL) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 15)
  )

ggsave("results/fig3_hazard_coefficients.png", p, width = 9.5, height = 5.5, dpi = 300)
cat("\nFigure 3 regenerated (greyscale): results/fig3_hazard_coefficients.png\n")
