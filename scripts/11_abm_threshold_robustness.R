# ------------------------------------------------------------------
# 11: ABM threshold robustness -- does the qualitative result (Model B
#     reaches higher final balance / lower residual imbalance than
#     Model A on every network) hold across different dissolution-
#     eligibility thresholds, not just the 25th percentile used as the
#     primary specification?
# ------------------------------------------------------------------
suppressMessages({ library(dplyr); library(readr); library(igraph); library(ggplot2) })

# Run from the replication pack root, e.g.:
#   cd replication_pack && Rscript scripts/11_abm_threshold_robustness.R
# (relative paths below assume that working directory; no setwd() needed)
set.seed(20260823)

SNAPSHOT_YEARS <- c(1914, 1949, 1985, 2005)
THRESHOLDS <- c(0.10, 0.25, 0.50)
N_TICKS <- 6000
N_REPS  <- 30  # reduced from 40 given 3x more conditions; still ample for stable means

signed_ties <- read_csv("data/signed_ties.csv", show_col_types = FALSE)

run_snapshot <- function(yr, q) {
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
  threshold <- quantile(edge_load_combined, q)
  low_load_edge <- edge_load_combined <= threshold

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
    final_state <- NULL
    for (t in seq_len(N_TICKS)) {
      closed <- active[tri_edges[, 1]] & active[tri_edges[, 2]] & active[tri_edges[, 3]]
      prods <- s[tri_edges[, 1]] * s[tri_edges[, 2]] * s[tri_edges[, 3]]
      unbal_idx <- which(closed & prods <= 0)
      if (length(unbal_idx) == 0) {
        final_state <- c(balanced = sum(closed & prods > 0) / n_tri, unbalanced = 0, dissolved = sum(!closed) / n_tri)
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
    if (is.null(final_state)) {
      closed <- active[tri_edges[, 1]] & active[tri_edges[, 2]] & active[tri_edges[, 3]]
      prods <- s[tri_edges[, 1]] * s[tri_edges[, 2]] * s[tri_edges[, 3]]
      final_state <- c(balanced = sum(closed & prods > 0) / n_tri, unbalanced = sum(closed & prods <= 0) / n_tri, dissolved = sum(!closed) / n_tri)
    }
    final_state
  }

  final_a <- rowMeans(replicate(N_REPS, run_sim(FALSE)))
  final_b <- rowMeans(replicate(N_REPS, run_sim(TRUE)))
  tibble(year = yr, threshold_pct = q * 100,
         model = c("A: Classical", "B: Exit-vs-adaptation"),
         balanced = c(final_a["balanced"], final_b["balanced"]),
         unbalanced = c(final_a["unbalanced"], final_b["unbalanced"]),
         dissolved = c(final_a["dissolved"], final_b["dissolved"]))
}

t0 <- Sys.time()
results <- bind_rows(lapply(SNAPSHOT_YEARS, function(yr) bind_rows(lapply(THRESHOLDS, function(q) run_snapshot(yr, q)))))
cat("Total time:", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "sec\n")

write_csv(results, "results/abm_threshold_robustness.csv")
print(results %>% arrange(year, threshold_pct, model), n = 50)

cat("\n=== Qualitative check: does Model B always reach higher balanced share than Model A? ===\n")
wide <- results %>% select(year, threshold_pct, model, balanced) %>%
  tidyr::pivot_wider(names_from = model, values_from = balanced)
names(wide) <- c("year", "threshold_pct", "A", "B")
wide <- wide %>% mutate(B_higher = B > A)
print(wide)
cat("\nB > A in", sum(wide$B_higher), "of", nrow(wide), "network x threshold combinations\n")

cat("\n11_abm_threshold_robustness.R complete.\n")
