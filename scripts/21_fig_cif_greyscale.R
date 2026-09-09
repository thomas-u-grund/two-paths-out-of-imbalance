# ------------------------------------------------------------------
# 21: Generate Figure 2 (results/fig_cif.png) -- cumulative incidence
#     of dissolution/realignment by embeddedness tercile at onset.
#     Uses results/cumulative_incidence.csv (tercile x years-since-onset
#     shares), computed by 07_competing_risks.R. Greyscale only (no
#     color), with larger axis/facet/legend text and a simplified,
#     minimal legend.
# ------------------------------------------------------------------
suppressMessages({ library(dplyr); library(ggplot2); library(readr) })

# Run from the replication pack root, e.g.:
#   cd replication_pack && Rscript scripts/21_fig_cif_greyscale.R
# (relative paths below assume that working directory; no setwd() needed)

cif <- read_csv("results/cumulative_incidence.csv", show_col_types = FALSE)

long <- cif %>%
  transmute(
    emb_tercile,
    dur,
    Dissolved = p_dissolved,
    `Flipped to balance` = p_flipped,
    `Still unbalanced (ongoing)` = p_ongoing
  ) %>%
  tidyr::pivot_longer(c(Dissolved, `Flipped to balance`, `Still unbalanced (ongoing)`),
                       names_to = "outcome", values_to = "share") %>%
  mutate(
    outcome = factor(outcome, levels = c("Dissolved", "Still unbalanced (ongoing)", "Flipped to balance")),
    emb_label = factor(emb_tercile, labels = c("Low embeddedness", "Medium embeddedness", "High embeddedness"))
  )

p <- ggplot(long, aes(x = dur, y = share, fill = outcome)) +
  geom_area(position = "stack", color = "grey20", linewidth = 0.15) +
  facet_wrap(~emb_label, nrow = 1) +
  scale_fill_manual(values = c(
    "Dissolved" = "grey75",
    "Still unbalanced (ongoing)" = "grey45",
    "Flipped to balance" = "grey15"
  )) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Years since spell onset", y = "Cumulative share of spells", fill = NULL) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 13),
    legend.text = element_text(size = 14),
    panel.spacing = unit(1.2, "lines"),
    panel.grid.minor = element_blank()
  )

ggsave("results/fig_cif.png", p, width = 11, height = 5, dpi = 300)
cat("Figure 2 written: results/fig_cif.png\n")
