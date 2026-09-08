# ------------------------------------------------------------------
# 05: Compact agent-based generative-sufficiency check (replaces full
#     NetLogo model -- see DECISIONS.md D6). This is now a direct
#     simulation of the exit-vs-adaptation theory (DECISIONS.md D27),
#     not the earlier "blocked flip" version: that version only let load
#     block or permit a sign flip, with no dissolution channel at all,
#     so it could not actually speak to the paper's real claim (that low
#     load produces dissolution and high load produces realignment) and
#     produced a confusing, hard-to-interpret pattern on larger networks
#     as a result (see D27 for the full diagnosis).
#
#     Two rules compared on each of four real historical interstate
#     networks (1914, 1949, 1985, 2005):
#       Model A: classical -- every tick, a randomly chosen unbalanced
#                closed triad has a randomly chosen edge FLIPPED. No
#                exit option exists, matching classical balance theory.
#       Model B: exit-vs-adaptation -- same triad/edge selection, but the
#                edge is DISSOLVED (removed from the graph) if its
#                combined load is in the lowest quartile of that
#                network's own load distribution, and FLIPPED otherwise.
#     Tracked over time: the share of the network's ORIGINAL closed
#     triads that are (a) still closed and balanced, (b) still closed
#     and unbalanced, (c) no longer closed (a member edge dissolved) --
#     directly analogous to the empirical flipped/frozen/dissolved
#     partition (Figure 2), which this simulation is meant to
#     generatively reproduce.
# ------------------------------------------------------------------
suppressMessages({ library(dplyr); library(readr); library(igraph); library(ggplot2) })

# Run from the replication pack root, e.g.:
#   cd replication_pack && Rscript scripts/05_abm_robustness.R
# (relative paths below assume that working directory; no setwd() needed)
set.seed(20260823)

SNAPSHOT_YEARS <- c(1914, 1949, 1985, 2005)
N_TICKS <- 6000
N_REPS  <- 40

signed_ties <- read_csv("data/signed_ties.csv", show_col_types = FALSE)

run_snapshot <- function(yr) {
  d0 <- signed_ties %>% filter(year == yr) %>%
    mutate(ccode1 = as.character(ccode_low), ccode2 = as.character(ccode_high)) %>%
    select(ccode1, ccode2, sign) %>% distinct(ccode1, ccode2, .keep_all = TRUE)
  g0 <- graph_from_data_frame(d0[, c("ccode1", "ccode2")], directed = FALSE)
  E(g0)$sign <- d0$sign
  n_edges <- ecount(g0)

  node_constraint <- tryCatch(constraint(g0), error = function(e) rep(0, vcount(g0)))
  node_constraint[is.na(node_constraint)] <- median(node_constraint, na.rm = TRUE)
  names(node_constraint) <- V(g0)$name

  edge_btw <- edge_betweenness(g0)
  edge_emb <- sapply(seq_len(n_edges), function(i) {
    e <- ends(g0, i); length(intersect(neighbors(g0, e[1]), neighbors(g0, e[2])))
  })
  tie_load_raw <- scale(edge_btw)[, 1] + scale(edge_emb)[, 1]
  actor_load_raw <- scale(node_constraint)[, 1]; names(actor_load_raw) <- V(g0)$name

  el <- as_data_frame(g0, what = "edges")
  edge_load_combined <- tie_load_raw + actor_load_raw[el$from] + actor_load_raw[el$to]
  threshold_q25 <- quantile(edge_load_combined, 0.25)
  low_load_edge <- edge_load_combined <= threshold_q25  # dissolution-eligible edges

  tri <- igraph::triangles(g0)
  if (length(tri) == 0) return(NULL)
  tri_mat <- matrix(tri, ncol = 3, byrow = TRUE)
  eid <- function(a, b) igraph::get_edge_ids(g0, c(a, b))
  tri_edges <- t(apply(tri_mat, 1, function(v) c(eid(v[1], v[2]), eid(v[1], v[3]), eid(v[2], v[3]))))
  n_tri <- nrow(tri_edges)
  init_signs <- E(g0)$sign

  run_sim <- function(allow_dissolution) {
    s <- init_signs
    active <- rep(TRUE, n_edges)
    out <- matrix(NA_real_, nrow = N_TICKS, ncol = 3)  # cols: balanced, unbalanced, dissolved (shares of n_tri)
    for (t in seq_len(N_TICKS)) {
      closed <- active[tri_edges[, 1]] & active[tri_edges[, 2]] & active[tri_edges[, 3]]
      prods <- s[tri_edges[, 1]] * s[tri_edges[, 2]] * s[tri_edges[, 3]]
      unbal_idx <- which(closed & prods <= 0)
      out[t, 1] <- sum(closed & prods > 0) / n_tri
      out[t, 2] <- length(unbal_idx) / n_tri
      out[t, 3] <- sum(!closed) / n_tri
      if (length(unbal_idx) == 0) {
        out[t:N_TICKS, 1] <- out[t, 1]; out[t:N_TICKS, 2] <- 0; out[t:N_TICKS, 3] <- out[t, 3]
        break
      }
      pick <- unbal_idx[sample.int(length(unbal_idx), 1)]
      which_edge <- sample.int(3, 1); e <- tri_edges[pick, which_edge]
      if (allow_dissolution && low_load_edge[e]) {
        active[e] <- FALSE
      } else {
        s[e] <- -s[e]
      }
    }
    out
  }

  reps_a <- replicate(N_REPS, run_sim(FALSE), simplify = "array")  # ticks x 3 x reps
  reps_b <- replicate(N_REPS, run_sim(TRUE),  simplify = "array")

  summarize_reps <- function(reps, model_label) {
    tibble(
      tick = rep(1:N_TICKS, 3),
      state = rep(c("Flipped to balance", "Still unbalanced", "Dissolved"), each = N_TICKS),
      mean_share = c(apply(reps[, 1, ], 1, mean), apply(reps[, 2, ], 1, mean), apply(reps[, 3, ], 1, mean)),
      model = model_label
    )
  }

  list(
    year = yr, n_nodes = vcount(g0), n_edges = n_edges, n_tri = n_tri,
    summary_df = bind_rows(summarize_reps(reps_a, "A: Classical (flip only)"),
                            summarize_reps(reps_b, "B: Exit-vs-adaptation")),
    final_a = apply(reps_a[N_TICKS, , ], 1, mean),
    final_b = apply(reps_b[N_TICKS, , ], 1, mean)
  )
}

t0 <- Sys.time()
results <- lapply(SNAPSHOT_YEARS, run_snapshot)
cat("Total sim time:", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "sec\n")

for (r in results) {
  cat(sprintf("Year %d: %d nodes, %d edges, %d triangles\n", r$year, r$n_nodes, r$n_edges, r$n_tri))
  cat(sprintf("  Model A final -- balanced: %.1f%% | unbalanced: %.1f%% | dissolved: %.1f%%\n",
              100*r$final_a[1], 100*r$final_a[2], 100*r$final_a[3]))
  cat(sprintf("  Model B final -- balanced: %.1f%% | unbalanced: %.1f%% | dissolved: %.1f%%\n",
              100*r$final_b[1], 100*r$final_b[2], 100*r$final_b[3]))
}

all_summary <- bind_rows(lapply(results, function(r) r$summary_df %>% mutate(year = r$year)))
all_summary$state <- factor(all_summary$state, levels = c("Dissolved", "Still unbalanced", "Flipped to balance"))
write_csv(all_summary, "results/abm_summary.csv")

p <- ggplot(all_summary, aes(x = tick, y = mean_share, fill = state)) +
  geom_area(position = "stack", alpha = 0.9) +
  facet_grid(model ~ year, labeller = labeller(year = function(y) paste0(y))) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("Dissolved" = "#c9a227", "Still unbalanced" = "#a83232", "Flipped to balance" = "#3a7d44")) +
  labs(x = "Simulation tick", y = "Share of the network's original closed triads", fill = NULL,
       title = "Classical dynamics vs. exit-vs-adaptation dynamics, across four historical networks",
       subtitle = "Model A has no exit option; Model B lets low-load ties dissolve and high-load ties realign") +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
ggsave("results/fig_abm.png", p, width = 11, height = 5.5, dpi = 200)

saveRDS(results, "results/abm_results.rds")
cat("05_abm_robustness.R complete.\n")
