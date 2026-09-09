# ------------------------------------------------------------------
# 02: Enumerate closed triads per year, classify balance, compute load
# ------------------------------------------------------------------
suppressMessages({
  library(dplyr); library(readr); library(tidyr); library(igraph)
})

# Run from the replication pack root, e.g.:
#   cd replication_pack && Rscript scripts/02_triads_and_load.R
# (relative paths below assume that working directory; no setwd() needed)

signed_ties <- read_csv("data/signed_ties.csv", show_col_types = FALSE) %>%
  # 1816-2012, not 2014: alliance data (v4.1) right-censors at 2012 while
  # MID data runs to 2014, so 2013-2014 have conflict ties but essentially
  # no alliance ties -- a data-boundary artifact, not signal.
  # This filter also drops 2 stray year=0 rows (raw MID data artifact).
  filter(year >= 1816, year <= 2012)
nmc <- read_csv("raw/nmc/NMCv7/abridged/NMC-70-abridged.csv", show_col_types = FALSE) %>%
  select(ccode, year, cinc) %>%
  mutate(log_cinc = log(pmax(cinc, 1e-6)))

years <- sort(unique(signed_ties$year))
cat("Years to process:", min(years), "-", max(years), "(", length(years), "years )\n")

process_year <- function(yr) {
  d <- signed_ties %>% filter(year == yr) %>%
    mutate(ccode1 = as.character(ccode_low), ccode2 = as.character(ccode_high)) %>%
    select(ccode1, ccode2, sign) %>%
    distinct(ccode1, ccode2, .keep_all = TRUE)

  if (nrow(d) < 3) return(NULL)

  # edge order in graph == row order in d (graph_from_data_frame preserves it)
  g <- graph_from_data_frame(d[, c("ccode1", "ccode2")], directed = FALSE)
  E(g)$sign <- d$sign

  tri <- igraph::triangles(g)
  if (length(tri) == 0) return(NULL)
  tri_mat <- matrix(V(g)$name[tri], ncol = 3, byrow = TRUE)

  # --- node-level metrics on the full network at year t ---
  deg <- degree(g)
  con <- tryCatch(constraint(g), error = function(e) rep(NA_real_, vcount(g)))
  node_tbl <- tibble(node = V(g)$name, degree = deg, constraint = con)
  get_node <- function(x, col) node_tbl[[col]][match(x, node_tbl$node)]

  # --- edge-level metrics, keyed the same way for sign / betweenness / embeddedness ---
  el <- as_data_frame(g, what = "edges") %>%
    mutate(e1 = pmin(from, to), e2 = pmax(from, to))
  edge_btw <- edge_betweenness(g)
  embeddedness <- sapply(seq_len(ecount(g)), function(i) {
    ends_i <- ends(g, i)
    length(intersect(neighbors(g, ends_i[1]), neighbors(g, ends_i[2])))
  })
  edge_tbl <- tibble(e1 = el$e1, e2 = el$e2, sign = E(g)$sign,
                      edge_btw = edge_btw, embeddedness = embeddedness)
  get_edge <- function(a, b, col) {
    key <- paste(pmin(a, b), pmax(a, b))
    edge_tbl[[col]][match(key, paste(edge_tbl$e1, edge_tbl$e2))]
  }

  n1 <- tri_mat[, 1]; n2 <- tri_mat[, 2]; n3 <- tri_mat[, 3]

  # IMPORTANT: index/label edges off the SORTED node1<node2<node3 output
  # ordering, not off the arbitrary n1/n2/n3 order igraph::triangles()
  # happens to enumerate in. The balance classification itself doesn't
  # care (the sign product is order-invariant), but tracking a SPECIFIC
  # tie's identity across years (needed for H2 and dissolution
  # disaggregation) does -- an inconsistent labeling would silently
  # compare different physical edges across years for the same triad_id.
  ids_sorted <- t(apply(tri_mat, 1, sort))
  s1 <- ids_sorted[, 1]; s2 <- ids_sorted[, 2]; s3 <- ids_sorted[, 3]
  idx12 <- match(paste(s1, s2), paste(edge_tbl$e1, edge_tbl$e2))  # edge (node1,node2)
  idx13 <- match(paste(s1, s3), paste(edge_tbl$e1, edge_tbl$e2))  # edge (node1,node3)
  idx23 <- match(paste(s2, s3), paste(edge_tbl$e1, edge_tbl$e2))  # edge (node2,node3)

  sign12 <- edge_tbl$sign[idx12]; sign13 <- edge_tbl$sign[idx13]; sign23 <- edge_tbl$sign[idx23]
  balanced <- (sign12 * sign13 * sign23) > 0

  # triad load: overlap = how many OTHER closed triads share each edge
  edge_in_triangle <- table(c(idx12, idx13, idx23))
  tri_overlap_mean <- (
    as.numeric(edge_in_triangle[as.character(idx12)]) +
      as.numeric(edge_in_triangle[as.character(idx13)]) +
      as.numeric(edge_in_triangle[as.character(idx23)])
  ) / 3 - 1

  tibble(
    year = yr,
    node1 = ids_sorted[, 1], node2 = ids_sorted[, 2], node3 = ids_sorted[, 3],
    sign12, sign13, sign23,
    balanced = balanced,
    actor_load_mean = (get_node(n1, "constraint") + get_node(n2, "constraint") + get_node(n3, "constraint")) / 3,
    actor_load_min  = pmin(get_node(n1, "constraint"), get_node(n2, "constraint"), get_node(n3, "constraint")),
    deg_mean = (get_node(n1, "degree") + get_node(n2, "degree") + get_node(n3, "degree")) / 3,
    # per-edge values (not just the triad mean) -- needed to identify which
    # specific tie changed at realignment (H2) or disappeared at dissolution
    tie12_emb = edge_tbl$embeddedness[idx12], tie13_emb = edge_tbl$embeddedness[idx13], tie23_emb = edge_tbl$embeddedness[idx23],
    tie12_btw = edge_tbl$edge_btw[idx12], tie13_btw = edge_tbl$edge_btw[idx13], tie23_btw = edge_tbl$edge_btw[idx23],
    tie_btw_mean = (tie12_btw + tie13_btw + tie23_btw) / 3,
    tie_emb_mean = (tie12_emb + tie13_emb + tie23_emb) / 3,
    triad_overlap = tri_overlap_mean,
    triad_id = paste(ids_sorted[, 1], ids_sorted[, 2], ids_sorted[, 3], sep = "-")
  )
}

results <- vector("list", length(years))
t0 <- Sys.time()
for (i in seq_along(years)) {
  yr <- years[i]
  res <- tryCatch(process_year(yr), error = function(e) {
    message("YEAR ", yr, " FAILED: ", conditionMessage(e)); NULL
  })
  results[[i]] <- res
  if (i %% 20 == 0) cat("...", yr, " (", i, "/", length(years), ") elapsed:",
                         round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
}
cat("Total elapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")

triads_all <- bind_rows(results)
cat("Total closed triads across all years:", nrow(triads_all), "\n")
cat("Balanced:", sum(triads_all$balanced), " Unbalanced:", sum(!triads_all$balanced), "\n")

# --- Attach capability control (mean log CINC of the 3 members) ---
nmc_lookup <- nmc %>% mutate(ccode = as.character(ccode)) %>% select(ccode, year, log_cinc)
triads_all <- triads_all %>%
  left_join(nmc_lookup, by = c("node1" = "ccode", "year")) %>% rename(cinc1 = log_cinc) %>%
  left_join(nmc_lookup, by = c("node2" = "ccode", "year")) %>% rename(cinc2 = log_cinc) %>%
  left_join(nmc_lookup, by = c("node3" = "ccode", "year")) %>% rename(cinc3 = log_cinc) %>%
  rowwise() %>%
  mutate(triad_log_cinc_mean = mean(c(cinc1, cinc2, cinc3), na.rm = TRUE)) %>%
  ungroup() %>%
  select(-cinc1, -cinc2, -cinc3)

write_csv(triads_all, "data/triads_all_years.csv")
saveRDS(triads_all, "data/triads_all_years.rds")

by_year <- triads_all %>% count(year, balanced) %>%
  pivot_wider(names_from = balanced, values_from = n, values_fill = 0) %>%
  rename(unbalanced = `FALSE`, balanced = `TRUE`)
write_csv(by_year, "data/triads_summary_by_year.csv")

cat("02_triads_and_load.R complete.\n")
