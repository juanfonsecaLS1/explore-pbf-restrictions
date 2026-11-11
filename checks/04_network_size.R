# Size and complexity of the two representations of the same road network:
#   primal graph:  junctions are vertices, directed road edges are edges
#                  (no turn restrictions: it cannot hold them)
#   line graph:    directed road edges are vertices, allowed turns are edges
#                  (restrictions, U-turn rules and via-way copies included)
# for the study polygon and for larger areas (buffers around it), so the
# growth with network size can be seen.
# Outputs: outputs/network_size.csv, outputs/network_size.png
source("processing/config.R")
suppressMessages(library(ggplot2))

BUFFERS_M <- c(0, 1000, 3000, 6000)
N_SOURCES <- 100  # one-to-all shortest-path trees timed per graph

graph_stats <- function(g, sources) {
  dout <- degree(g, mode = "out"); din <- degree(g, mode = "in")
  cc <- components(g, mode = "strong")
  t <- system.time(distances(g, v = sources, mode = "out"))[["elapsed"]]
  data.table(
    vertices = vcount(g), edges = ecount(g),
    edges_per_vertex = ecount(g) / vcount(g),
    max_out_degree = max(dout), max_in_degree = max(din),
    density = ecount(g) / (vcount(g) * (vcount(g) - 1)),
    largest_scc_pct = 100 * max(cc$csize) / vcount(g),
    memory_mb = as.numeric(object.size(g)) / 2^20,
    sssp_ms = 1000 * t / length(sources)
  )
}

pbf <- colombia_pbf()
poly <- st_transform(study_polygon(), 3116)
set.seed(1)
res <- list(); deg <- list()
for (b in BUFFERS_M) {
  area <- if (b == 0) poly else st_buffer(poly, b)
  area_km2 <- as.numeric(st_area(area)) / 1e6
  net <- osm_driving_network(pbf, st_transform(area, 4326), osmium = OSMIUM)
  lg <- build_line_graph(net$edges, net$restrictions, nodes = net$nodes)
  e <- as.data.table(st_drop_geometry(net$edges))
  tt <- edge_travel_time(net$edges)

  gp <- graph_from_data_frame(data.frame(from = as.character(e$from),
                                         to = as.character(e$to), weight = tt),
                              directed = TRUE)
  gl <- line_graph_igraph(lg)
  # same origins in both: an edge in the line graph, its start node in the
  # primal graph
  src_e <- sample(lg$vertices$vertex_id[!grepl("#", lg$vertices$vertex_id)], N_SOURCES)
  src_n <- unique(as.character(e$from[match(src_e, e$edge_id)]))

  # Turns: every in-edge x out-edge pair at a node, and what the rules removed
  pairs <- sum(as.numeric(degree(gp, mode = "in")) * degree(gp, mode = "out"))
  uturn_pairs <- nrow(merge(e[, .(node = to, link_id)], e[, .(node = from, link_id)],
                            by = c("node", "link_id")))
  tp <- lg$turns
  res[[length(res) + 1]] <- rbind(
    cbind(buffer_m = b, area_km2 = area_km2, graph = "primal",
          graph_stats(gp, src_n)),
    cbind(buffer_m = b, area_km2 = area_km2, graph = "line graph",
          graph_stats(gl, src_e)),
    fill = TRUE
  )[, `:=`(
    in_x_out_pairs = pairs, uturn_pairs = uturn_pairs,
    uturns_kept = sum(grepl("uturn", tp$turn_type)),
    removed_by_node_via = lg$report$turns_removed_node_via,
    via_way_copies = sum(grepl("#", lg$vertices$vertex_id))
  )]
  if (b == 0) {
    deg <- rbind(
      data.table(graph = "primal", out_degree = degree(gp, mode = "out")),
      data.table(graph = "line graph", out_degree = degree(gl, mode = "out"))
    )
  }
  cat(sprintf("buffer %5d m: %6.1f km2, %6d junctions, %6d edges, %6d turns\n",
              b, area_km2, vcount(gp), ecount(gp), ecount(gl)))
}
res <- rbindlist(res)
fwrite(res, file.path(OUT_DIR, "network_size.csv"))

cat("\n== Size and structure ==\n")
print(res[, .(buffer_m, area_km2 = round(area_km2), graph, vertices, edges,
              e_per_v = round(edges_per_vertex, 2), max_out = max_out_degree,
              scc_pct = round(largest_scc_pct, 1), memory_mb = round(memory_mb, 2),
              sssp_ms = round(sssp_ms, 1))])

cat("\n== Line graph relative to primal ==\n")
w <- dcast(res, buffer_m ~ graph, value.var = c("vertices", "edges", "memory_mb", "sssp_ms"))
print(w[, .(buffer_m,
            vertices_x = round(`vertices_line graph` / vertices_primal, 2),
            edges_x = round(`edges_line graph` / edges_primal, 2),
            memory_x = round(`memory_mb_line graph` / memory_mb_primal, 2),
            sssp_x = round(`sssp_ms_line graph` / sssp_ms_primal, 2))])

cat("\n== Where the turns come from (line-graph edges) ==\n")
print(unique(res[graph == "line graph", .(buffer_m, in_x_out_pairs, uturn_pairs,
                                          uturns_kept, removed_by_node_via,
                                          via_way_copies, turns = edges)]))

cat("\n== Out-degree distribution, study polygon ==\n")
print(dcast(deg[, .N, by = .(graph, out_degree)], out_degree ~ graph, value.var = "N", fill = 0))

# Figure: growth with area, and out-degree distribution
long <- melt(res[, .(area_km2, graph, vertices, edges)], id.vars = c("area_km2", "graph"))
p1 <- ggplot(long, aes(area_km2, value, colour = graph, linetype = variable)) +
  geom_line() + geom_point() +
  scale_colour_manual(values = c(primal = "#0072B2", `line graph` = "#D55E00")) +
  labs(x = "area (km²)", y = "count", colour = NULL, linetype = NULL,
       title = "Size against area") +
  theme_minimal()
p2 <- ggplot(deg[, .N, by = .(graph, out_degree)][, share := N / sum(N), by = graph],
             aes(factor(out_degree), share, fill = graph)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c(primal = "#0072B2", `line graph` = "#D55E00")) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "out-degree", y = "share of vertices", fill = NULL,
       title = "Out-degree, study polygon") +
  theme_minimal()
ggsave(file.path(OUT_DIR, "network_size.png"), patchwork::wrap_plots(p1, p2),
       width = 10, height = 4, dpi = 150, bg = "white")
cat("saved", file.path(OUT_DIR, "network_size.csv"), "and network_size.png\n")
