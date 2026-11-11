# Routability checks on the line graph built by processing/01_build_network.R:
#   1. connectivity (largest strongly connected component)
#   2. a route through every turn restriction must avoid the banned movement
#   3. random O-D pairs: reachability, validity, comparison with a naive graph
#   4. edges outside the main component (sinks, sources, isolated)
# Writes layers routes, isolated_edges and sink_source_edges to
# outputs/routability_examples.gpkg (the file is recreated).
source("processing/config.R")

x <- readRDS(NETWORK_RDS); net <- x$net; lg <- x$lg
poly <- study_polygon()
g <- line_graph_igraph(lg)
e <- as.data.table(st_drop_geometry(net$edges))
tt <- edge_travel_time(net$edges); names(tt) <- e$edge_id

cat("== 1. Connectivity ==\n")
cc <- components(g, mode = "strong")
inbig <- cc$membership == which.max(cc$csize)
orig_in <- tapply(inbig, sub("#.*", "", V(g)$name), any)
cat("vertices in largest strongly connected component:", sum(inbig), "of",
    vcount(g), sprintf("(%.1f%%)", 100 * mean(inbig)), "\n")
out_edges <- names(orig_in)[!orig_in]
edge_inside <- as.vector(st_within(net$edges, st_geometry(poly), sparse = FALSE))
names(edge_inside) <- e$edge_id
cat("road edges outside it:", length(out_edges),
    " of which fully inside the polygon:", sum(edge_inside[out_edges]), "\n")
oi <- e[edge_id %in% out_edges & edge_inside[edge_id]]
print(table(highway = oi$highway, oneway = oi$oneway))

cat("\n== 2. Routing through every restriction ==\n")
rx <- car_restrictions(net$restrictions)
rn <- rx[via_type == "node" & lengths(via) == 1,
         .(relation_id, restriction, kind, from = unlist(from),
           via = unlist(via), to = unlist(to))]
rn <- merge(merge(rn, e[, .(from = way_id, via = to, e1 = edge_id)],
                  by = c("from", "via")),
            e[, .(to = way_id, via = from, e2 = edge_id)], by = c("to", "via"))
# Banned (from, to) road-edge pairs: 'no' turns plus the turns 'only' excludes
all_turns <- merge(e[, .(e1 = edge_id, node = to, l1 = link_id)],
                   e[, .(e2 = edge_id, node = from, l2 = link_id)],
                   by = "node", allow.cartesian = TRUE)[e1 != e2]
only_ex <- merge(rn[kind == "only", .(relation_id, restriction, e1, keep = e2, node = via)],
                 all_turns, by = c("e1", "node"))[e2 != keep & l1 != l2]
tests <- rbind(rn[kind == "no", .(relation_id, restriction, e1, e2, type = "node")],
               only_ex[, .(relation_id, restriction, e1, e2, type = "node")])
way_paths <- lg$report$way_via_paths
tests <- rbind(tests, data.table(
  relation_id = as.numeric(names(way_paths)),
  restriction = rx$restriction[match(as.numeric(names(way_paths)), rx$relation_id)],
  e1 = sapply(way_paths, `[`, 1),
  e2 = sapply(way_paths, function(p) p[length(p)]),
  type = "way"
))
contains_seq <- function(r, s) {
  n <- length(s)
  if (length(r) < n) return(FALSE)
  any(sapply(seq_len(length(r) - n + 1), function(i) all(r[i:(i + n - 1)] == s)))
}
res <- rbindlist(lapply(seq_len(nrow(tests)), function(i) {
  t <- tests[i]
  r <- suppressWarnings(route_edges(g, t$e1, t$e2))
  banned <- if (t$type == "way") way_paths[[as.character(t$relation_id)]] else c(t$e1, t$e2)
  data.table(
    relation_id = t$relation_id, restriction = t$restriction, type = t$type,
    e1 = t$e1, e2 = t$e2, reachable = !is.null(r),
    uses_banned = !is.null(r) && contains_seq(r, banned),
    n_edges = if (is.null(r)) NA else length(r),
    cost = if (is.null(r)) NA else attr(r, "cost"),
    # travel time of the banned direct movement (edge times only)
    direct = if (t$type == "way") sum(tt[banned]) - 0.5 * (tt[t$e1] + tt[t$e2])
             else 0.5 * (tt[t$e1] + tt[t$e2])
  )
}))
res[, detour_s := cost - direct]
print(res[, .(tests = .N, reachable = sum(reachable),
              uses_banned = sum(uses_banned, na.rm = TRUE),
              median_detour_s = round(median(detour_s, na.rm = TRUE)),
              max_detour_s = round(max(detour_s, na.rm = TRUE))),
          by = .(type, kind = sub("_.*", "", restriction))])
cat("unreachable tests:\n")
print(res[reachable == FALSE, .(relation_id, restriction, e1, e2,
                                e1_in_scc = as.vector(orig_in[e1]),
                                e2_in_scc = as.vector(orig_in[e2]))])

cat("\n== 3. Random O-D pairs ==\n")
set.seed(1)
orig_v <- lg$vertices[!grepl("#", vertex_id), vertex_id]
od <- data.table(o = sample(orig_v, 1000, TRUE), d = sample(orig_v, 1000, TRUE))[o != d]
D <- distances(g, v = unique(od$o), mode = "out")
# reaching any copy of the destination edge counts
dest_cols <- split(seq_len(ncol(D)), sub("#.*", "", colnames(D)))
od[, cost := mapply(function(o, d) min(D[o, dest_cols[[d]]]), o, d)]
cat("reachable:", sum(is.finite(od$cost)), "of", nrow(od), "\n")
cat("unreachable explained by an edge outside the largest component:",
    sum(!is.finite(od$cost) & (!orig_in[od$o] | !orig_in[od$d])), "\n")
# Naive graph (all turns, U-turns everywhere, no restrictions, no penalties)
naive <- all_turns[, .(from = e1, to = e2, weight = 0.5 * tt[e1] + 0.5 * tt[e2])]
gn <- graph_from_data_frame(naive, directed = TRUE,
                            vertices = data.frame(name = e$edge_id))
Dn <- distances(gn, v = unique(od$o), mode = "out")
od[, cost_naive := Dn[cbind(o, d)]]
ok <- is.finite(od$cost)
cat("routes longer than naive:", sum(od$cost[ok] > od$cost_naive[ok] + 1e-6),
    "of", sum(ok), " median extra:", round(median((od$cost - od$cost_naive)[ok]), 1),
    "s  max extra:", round(max((od$cost - od$cost_naive)[ok])), "s\n")
cat("never shorter than naive:", all(od$cost[ok] >= od$cost_naive[ok] - 1e-6), "\n")

# Path validity on a sample: consecutive vertices are turns in the graph and
# no banned pair or sequence is used
tk <- paste(lg$turns$from, lg$turns$to)
banned_pairs <- paste(tests[type == "node", e1], tests[type == "node", e2])
bad <- 0
for (i in which(ok)[1:200]) {
  r <- route_edges(g, od$o[i], od$d[i]); v <- attr(r, "vertices")
  if (length(v) > 1 && !all(paste(head(v, -1), v[-1]) %in% tk)) bad <- bad + 1
  if (length(r) > 1 && any(paste(head(r, -1), r[-1]) %in% banned_pairs)) bad <- bad + 1
  if (any(sapply(way_paths, function(p) contains_seq(r, p)))) bad <- bad + 1
}
cat("invalid routes in 200 checked:", bad, "\n")

# Example routes for QGIS / checks/03_plot_routes.R: the longest detour of each
# restriction type plus the next four via-node ones
ex <- res[reachable == TRUE][order(-detour_s)][, .SD[1], by = type]
ex <- rbind(ex, res[reachable == TRUE & type == "node"][order(-detour_s)][2:5])
routes <- rbindlist(lapply(seq_len(nrow(ex)), function(i) {
  r <- route_edges(g, ex$e1[i], ex$e2[i])
  data.table(relation_id = ex$relation_id[i], restriction = ex$restriction[i],
             seq = seq_along(r), edge_id = as.character(r))
}))
rg <- merge(net$edges[, c("edge_id", "highway", "name")], routes, by = "edge_id")
rg <- rg[order(rg$relation_id, rg$seq), ]
st_write(rg, EXAMPLES_GPKG, layer = "routes", quiet = TRUE, delete_dsn = TRUE)
cat("exported", nrow(ex), "example routes; relations:",
    paste(ex$relation_id, collapse = ", "), "\n")

cat("\n== 4. Edges outside the main component ==\n")
gv <- which(inbig)[1]
fwd <- names(subcomponent(g, gv, mode = "out"))
bwd <- names(subcomponent(g, gv, mode = "in"))
oe <- out_edges[edge_inside[out_edges]]
cls <- ifelse(oe %in% fwd, "can enter, cannot leave (sink)",
              ifelse(oe %in% bwd, "can leave, cannot enter (source)", "isolated"))
print(table(cls))
iso <- oe[cls == "isolated"]
cat("isolated component sizes:\n")
print(table(cc$csize[cc$membership[match(iso, V(g)$name)]]))
ways_cut <- e$way_id[!edge_inside]
cat("sink/source edges whose OSM way leaves the polygon:",
    sum(e$way_id[match(oe[cls != "isolated"], e$edge_id)] %in% ways_cut), "of",
    sum(cls != "isolated"), "\n")
st_write(net$edges[net$edges$edge_id %in% iso, ], EXAMPLES_GPKG,
         layer = "isolated_edges", quiet = TRUE, append = FALSE)
sink_e <- net$edges[net$edges$edge_id %in% oe[cls != "isolated"], ]
sink_e$class <- cls[cls != "isolated"][match(sink_e$edge_id, oe[cls != "isolated"])]
st_write(sink_e, EXAMPLES_GPKG, layer = "sink_source_edges", quiet = TRUE,
         append = FALSE)
cat("saved", EXAMPLES_GPKG, "\n")
