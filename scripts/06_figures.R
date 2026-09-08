suppressMessages({ library(dplyr); library(readr); library(ggplot2); library(broom); library(scales) })
# Run from the replication pack root, e.g.:
#   cd replication_pack && Rscript scripts/06_figures.R
# (relative paths below assume that working directory; no setwd() needed)

# --- Figure 1: balanced vs unbalanced closed triads over time ---
by_year <- read_csv("data/triads_summary_by_year.csv", show_col_types = FALSE) %>%
  mutate(total = balanced + unbalanced, pct_unbalanced = unbalanced / total)

p1 <- ggplot(by_year, aes(x = year, y = pct_unbalanced)) +
  geom_line(color = "#a83232", linewidth = 0.4, alpha = 0.5) +
  geom_point(aes(size = total), color = "#a83232", alpha = 0.6) +
  scale_size_area(max_size = 4, name = "Closed triads\nthat year") +
  scale_y_continuous(labels = percent) +
  labs(x = NULL, y = "Share of closed triads unbalanced",
       title = "Unbalanced closed triads in the international system, 1816-2012",
       subtitle = "Point size = closed triads that year;\nsmall early points reflect thin data, not unstable trends") +
  theme_minimal(base_size = 12)
ggsave("results/fig1_imbalance_over_time.png", p1, width = 8, height = 4.5, dpi = 300)

# --- Figure 2: descriptive P(frozen) across embeddedness (compositional
#     quantity -- see DECISIONS.md D19: NOT a directly-estimated nonlinear
#     effect; both component contrasts below are monotonic) ---
m <- readRDS("results/models.rds")
mq <- m$m_frozen_pooled_quad
df <- m$df

newdat <- expand.grid(
  z_tie_emb = seq(quantile(df$z_tie_emb, .02), quantile(df$z_tie_emb, .98), length.out = 100),
  z_tie_btw = 0, z_actor_load = 0, z_cinc = 0, era = "1945-1989"
)
pred <- predict(mq, newdata = newdat, type = "link", se.fit = TRUE)
newdat$fit <- plogis(pred$fit)
newdat$lo <- plogis(pred$fit - 1.96 * pred$se.fit)
newdat$hi <- plogis(pred$fit + 1.96 * pred$se.fit)

p2 <- ggplot(newdat, aes(x = z_tie_emb, y = fit)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#2c5f8a", alpha = 0.15) +
  geom_line(color = "#2c5f8a", linewidth = 1) +
  scale_y_continuous(labels = percent) +
  labs(x = "Tie/triad embeddedness (standardized)", y = "Predicted share of all unbalanced triads that freeze by t+5",
       title = "Freezing's share peaks at intermediate embeddedness",
       subtitle = "A compositional pattern: dissolution and flipping both rise with embeddedness\n(at different rates), leaving frozen as a shrinking-then-recovering residual") +
  theme_minimal(base_size = 12)
ggsave("results/fig2_curvilinear_frozen.png", p2, width = 7.5, height = 5, dpi = 200)

# --- Figure 3: outcome shares by embeddedness tercile (illustrative;
#     purely descriptive, unaffected by D19) ---
df3 <- df %>% mutate(emb_tercile = ntile(tie_emb_mean, 3),
                      emb_label = factor(emb_tercile, labels = c("Low embeddedness", "Medium embeddedness", "High embeddedness")))
share_df <- df3 %>% count(emb_label, outcome) %>% group_by(emb_label) %>% mutate(share = n / sum(n))
share_df$outcome <- factor(share_df$outcome, levels = c("dissolved", "frozen", "flipped"))

p3 <- ggplot(share_df, aes(x = emb_label, y = share, fill = outcome)) +
  geom_col(position = "stack", width = 0.6) +
  scale_y_continuous(labels = percent) +
  scale_fill_manual(values = c(dissolved = "#c9a227", frozen = "#a83232", flipped = "#3a7d44"),
                     labels = c(dissolved = "Dissolved (tie disappears)", frozen = "Frozen (stays unbalanced)", flipped = "Flipped to balance")) +
  labs(x = NULL, y = "Share of unbalanced triad-years", fill = NULL,
       title = "How tension exits the system, by tie/triad embeddedness") +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom")
ggsave("results/fig3_outcome_shares.png", p3, width = 7, height = 5, dpi = 200)

# --- Figure 4: coefficient plot, the two PRIMARY contrasts side by side
#     (dissolved-vs-flipped, frozen-vs-flipped), linear specification ---
build_coef_df <- function(ct, contrast_label) {
  as.data.frame(ct[-1, , drop = FALSE]) %>%
    tibble::rownames_to_column("term") %>%
    rename(estimate = Estimate, se = `Std. Error`) %>%
    mutate(contrast = contrast_label,
           label = recode(term,
             z_tie_btw = "Tie load (bridging)",
             z_tie_emb = "Tie/triad load (embeddedness)",
             z_actor_load = "Actor load (constraint)",
             z_cinc = "Capability (log CINC)",
             `era1945-1989` = "Era: 1945-1989", `erapost-1989` = "Era: post-1989"),
           lo = estimate - 1.96 * se, hi = estimate + 1.96 * se) %>%
    filter(!is.na(label))
}
coef_df <- bind_rows(
  build_coef_df(m$ct_dissolved_lin, "Dissolved vs. flipped"),
  build_coef_df(m$ct_frozen_lin, "Frozen vs. flipped")
)

p4 <- ggplot(coef_df, aes(x = estimate, y = reorder(label, estimate), color = contrast)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("Dissolved vs. flipped" = "#c9a227", "Frozen vs. flipped" = "#a83232")) +
  labs(x = "Logit coefficient (95% CI, clustered SE by triad)", y = NULL, color = NULL,
       title = "What predicts NOT reaching classical balance-theoretic resolution",
       subtitle = "Two clean, monotonic contrasts against the reference category (flipped to balance)") +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom")
ggsave("results/fig4_coefficients.png", p4, width = 8, height = 4.5, dpi = 200)

cat("Figures 1-4 written to results/\n")
