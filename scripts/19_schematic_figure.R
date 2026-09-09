# ------------------------------------------------------------------
# 19: Schematic figure explaining triadic embeddedness and unbalanced
#     load, built up incrementally across several small, simple panels
#     rather than one figure trying to show the whole worked example
#     at once.
#
#     Panel sequence:
#       (a) a bare tie (two nodes, one edge) -- the starting point
#       (b) the same tie with ONE overlapping triangle -- embeddedness 1
#       (c) the same tie with TWO overlapping triangles -- embeddedness 2
#       (d) a tie with a BALANCED overlapping triangle -- unbalanced
#           load = 0
#       (e) the same tie, changing only the sign of one satellite edge,
#           making the SAME overlapping triangle UNBALANCED --
#           unbalanced load = 1. (d) and (e) are a deliberate minimal
#           pair: everything is held fixed except the one sign that
#           flips balance.
#       (f) the same two-triangle topology as (c), but now BOTH
#           satellite triangles are unbalanced -- unbalanced load = 2.
#           (c) and (f) are the parallel minimal pair for embeddedness
#           vs. load at count = 2: identical topology, differing only
#           in how many of the two overlapping triangles are unbalanced.
#
#     Panels (a)-(c) are sign-agnostic (embeddedness counts overlapping
#     triangles regardless of sign) and drawn with plain lines. Panels
#     (d)-(f) introduce the solid=positive/dashed=negative convention
#     used throughout the paper's other figures, only once sign starts
#     to matter for the concept being shown.
#
#     Caption text is drawn INSIDE each panel's own data space (a
#     geom_text below the diagram), not as a ggplot plot.title: plot
#     titles are not clipped to their own panel's width by
#     patchwork/ggplot2 the way in-panel data is, so a long caption in
#     the title would overflow sideways into the neighbouring panel.
#     Titles are just the bare panel letter.
# ------------------------------------------------------------------
suppressMessages({ library(ggplot2); library(dplyr); library(tibble); library(patchwork) })

# Run from the replication pack root, e.g.:
#   cd replication_pack && Rscript scripts/19_schematic_figure.R
# (relative paths below assume that working directory; no setwd() needed)

NODE_SIZE <- 11
NODE_STROKE <- 1.3
LABEL_SIZE <- 5.6
LETTER_SIZE <- 15
CAPTION_SIZE <- 4.9

base_theme <- theme_void(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5, size = LETTER_SIZE, face = "bold",
                                   margin = margin(b = 2)),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(4, 6, 4, 6))

# P, Q are the focal tie's endpoints in every panel, fixed in place so
# the eye has a stable anchor across all five panels. Every panel uses
# the SAME coordinate range (set in xy_range) so the "zoom level" and
# node size are visually consistent across the whole figure.
P <- c(x = -0.6, y = 0); Q <- c(x = 0.6, y = 0)
CAPTION_Y <- -1.35

node_df <- function(...) {
  pts <- list(...)
  tibble(name = names(pts), x = sapply(pts, `[`, "x"), y = sapply(pts, `[`, "y"))
}
seg <- function(a, b, linetype = "solid", weight = 1.6, color = "grey15") {
  tibble(x = a["x"], y = a["y"], xend = b["x"], yend = b["y"], linetype = linetype,
         linewidth = weight, color = color)
}

draw_nodes <- function(p, nodes) {
  p + geom_point(data = nodes, aes(x = x, y = y), shape = 21, fill = "white",
                 color = "grey15", size = NODE_SIZE, stroke = NODE_STROKE) +
    geom_text(data = nodes, aes(x = x, y = y, label = name), fontface = "bold", size = LABEL_SIZE)
}
draw_edges <- function(p, edges) {
  p + geom_segment(data = edges, aes(x = x, y = y, xend = xend, yend = yend,
                                      linetype = linetype, linewidth = linewidth, color = I(color))) +
    scale_linewidth_identity() +
    scale_linetype_identity()
}
add_caption <- function(p, text) {
  p + annotate("text", x = 0, y = CAPTION_Y, label = text, size = CAPTION_SIZE, lineheight = 0.95)
}

xy_range <- coord_fixed(clip = "off", xlim = c(-1.15, 1.15), ylim = c(-1.7, 1.25))

# --- (a) a bare tie ---
nodes_a <- node_df(P = P, Q = Q)
edges_a <- seg(P, Q)
pa <- ggplot() %>% draw_edges(edges_a) %>% draw_nodes(nodes_a)
pa <- pa %>% add_caption("A tie between\ntwo states")
pa <- pa + xy_range + labs(title = "(a)") + base_theme

# --- (b) one overlapping triangle ---
R <- c(x = 0, y = 0.85)
nodes_b <- node_df(P = P, Q = Q, R = R)
edges_b <- bind_rows(seg(P, Q), seg(P, R, weight = 1.1), seg(Q, R, weight = 1.1))
pb <- ggplot() %>% draw_edges(edges_b) %>% draw_nodes(nodes_b)
pb <- pb %>% add_caption("One overlapping triangle\nEmbeddedness = 1")
pb <- pb + xy_range + labs(title = "(b)") + base_theme

# --- (c) two overlapping triangles ---
R1 <- c(x = 0, y = 0.85); R2 <- c(x = 0, y = -0.85)
nodes_c <- node_df(P = P, Q = Q, R1 = R1, R2 = R2)
edges_c <- bind_rows(seg(P, Q), seg(P, R1, weight = 1.1), seg(Q, R1, weight = 1.1),
                      seg(P, R2, weight = 1.1), seg(Q, R2, weight = 1.1))
pc <- ggplot() %>% draw_edges(edges_c) %>% draw_nodes(nodes_c)
pc <- pc %>% add_caption("Two overlapping triangles\nEmbeddedness = 2")
pc <- pc + xy_range + labs(title = "(c)") + base_theme

# --- (d)/(e): minimal pair -- same P-Q, same P-R; only Q-R's sign
#     changes, flipping the triangle from balanced to unbalanced ---
tri_poly <- tibble(x = c(P["x"], Q["x"], R["x"]), y = c(P["y"], Q["y"], R["y"]))
nodes_de <- node_df(P = P, Q = Q, R = R)

make_de_panel <- function(qr_sign, fill, title, caption) {
  edges <- bind_rows(
    seg(P, Q, linetype = "solid"),                                  # P-Q always positive
    seg(P, R, linetype = "solid", weight = 1.3),                    # P-R always positive
    seg(Q, R, linetype = if (qr_sign == 1) "solid" else "22", weight = 1.3)  # Q-R varies
  )
  base_p <- ggplot() + geom_polygon(data = tri_poly, aes(x = x, y = y), fill = fill, color = NA)
  p <- base_p %>% draw_edges(edges) %>% draw_nodes(nodes_de)
  p <- p %>% add_caption(caption)
  p + xy_range + labs(title = title) + base_theme
}

pd <- make_de_panel(qr_sign = 1, fill = "white", title = "(d)",
                     caption = "Balanced overlap\nUnbalanced load = 0")
pe <- make_de_panel(qr_sign = -1, fill = "grey55", title = "(e)",
                     caption = "Unbalanced overlap\nUnbalanced load = 1")

# --- (f): same two-triangle topology as (c), but BOTH satellite
#     triangles are now unbalanced (each has one negative edge off the
#     positive P-Q tie), so unbalanced load = 2. Parallels (c) exactly
#     in topology, differing only in sign. ---
tri_poly_f1 <- tibble(x = c(P["x"], Q["x"], R1["x"]), y = c(P["y"], Q["y"], R1["y"]), grp = "t1")
tri_poly_f2 <- tibble(x = c(P["x"], Q["x"], R2["x"]), y = c(P["y"], Q["y"], R2["y"]), grp = "t2")
tri_poly_f <- bind_rows(tri_poly_f1, tri_poly_f2)
nodes_f <- node_df(P = P, Q = Q, R1 = R1, R2 = R2)
edges_f <- bind_rows(
  seg(P, Q, linetype = "solid"),
  seg(P, R1, linetype = "solid", weight = 1.1), seg(Q, R1, linetype = "22", weight = 1.1),
  seg(P, R2, linetype = "solid", weight = 1.1), seg(Q, R2, linetype = "22", weight = 1.1)
)
pf <- ggplot() + geom_polygon(data = tri_poly_f, aes(x = x, y = y, group = grp), fill = "grey55", color = NA)
pf <- pf %>% draw_edges(edges_f) %>% draw_nodes(nodes_f)
pf <- pf %>% add_caption("Two unbalanced overlaps\nUnbalanced load = 2")
pf <- pf + xy_range + labs(title = "(f)") + base_theme

# --- legend row, spanning the full width beneath both panel rows ---
legend_df <- tibble(
  label = c("Positive (alliance)", "Negative (dispute)"),
  x = c(0.15, 0.15),
  y = c(0.15, -0.15),
  linetype = c("solid", "22")
)
p_legend <- ggplot() +
  geom_segment(data = legend_df, aes(x = x - 0.35, xend = x - 0.05, y = y, yend = y, linetype = linetype),
               color = "grey15", linewidth = 1.1) +
  geom_text(data = legend_df, aes(x = x, y = y, label = label), hjust = 0, size = 4.9) +
  scale_linetype_identity() +
  coord_cartesian(clip = "off", xlim = c(-0.6, 2.4), ylim = c(-0.5, 0.5)) +
  theme_void(base_size = 15) +
  theme(plot.margin = margin(2, 6, 2, 6))

top_row <- pa | pb | pc
mid_row <- pd | pe | pf
combined <- top_row / mid_row / p_legend + plot_layout(heights = c(1, 1, 0.22))

ggsave("results/fig_schematic_measures.png", combined, width = 10.5, height = 7.6, dpi = 300, bg = "white")
cat("Saved results/fig_schematic_measures.png\n")
