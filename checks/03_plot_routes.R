# Static maps of the example routes (layer "routes" of
# outputs/routability_examples.gpkg): route, origin and destination, the banned
# direct movement and the restriction's via member, over an OSM basemap
# (tiles downloaded by tmap / maptiles). One PNG per route in
# outputs/route_plots/. These are sense-check plots, not publication figures.
source("processing/config.R")
suppressMessages(library(tmap))
dir.create(PLOT_DIR, showWarnings = FALSE)

x <- readRDS(NETWORK_RDS); net <- x$net; lg <- x$lg
g <- line_graph_igraph(lg)
routes <- st_read(EXAMPLES_GPKG, "routes", quiet = TRUE)
rx <- car_restrictions(net$restrictions)
edges <- st_transform(net$edges, 3116)            # MAGNA-SIRGAS Bogotá, metres
e <- as.data.table(st_drop_geometry(edges))
tt <- edge_travel_time(net$edges); names(tt) <- e$edge_id
geom_of <- function(ids) st_geometry(edges)[match(ids, edges$edge_id)]

# Point at fraction f along a line, moved `off` metres to the right of the
# direction of travel, so the two directions of a two-way road do not overlap
pt_right <- function(line, f, off = 14) {
  s <- st_coordinates(st_cast(st_line_sample(line, sample = pmin(pmax(c(f - 0.03, f, f + 0.03), 0), 1)), "POINT"))
  d <- s[3, 1:2] - s[1, 1:2]; d <- d / sqrt(sum(d^2))
  st_point(s[2, 1:2] + off * c(d[2], -d[1]))
}
end_xy <- function(line) { xy <- st_coordinates(line); xy[nrow(xy), 1:2] }

for (id in unique(routes$relation_id)) {
  r <- routes[routes$relation_id == id, ]
  path <- r$edge_id[order(r$seq)]
  rel <- rx[relation_id == id]
  re <- edges[match(path, edges$edge_id), ]
  o <- path[1]; d <- path[length(path)]
  cost <- attr(route_edges(g, o, d), "cost")
  len_km <- sum(e$length_m[match(path, e$edge_id)]) / 1000
  banned <- lg$report$way_via_paths[[as.character(id)]]
  if (is.null(banned)) banned <- c(o, d)
  direct_s <- if (rel$via_type == "way") sum(tt[banned]) - 0.5 * (tt[o] + tt[d]) else 0.5 * (tt[o] + tt[d])

  o_pt <- pt_right(geom_of(o), 0.5); d_pt <- pt_right(geom_of(d), 0.5)
  # Very short edges: push the two markers apart so neither hides the other
  v <- as.numeric(d_pt) - as.numeric(o_pt); dist <- sqrt(sum(v^2))
  if (dist < 35) {
    u <- if (dist > 1) v / dist else c(1, 0)
    o_pt <- st_point(as.numeric(o_pt) - u * (35 - dist) / 2)
    d_pt <- st_point(as.numeric(d_pt) + u * (35 - dist) / 2)
  }
  # Banned movement drawn as O -> via -> D
  if (rel$via_type == "node") {
    via_xy <- end_xy(geom_of(o))
    banned_xy <- rbind(st_coordinates(st_sfc(o_pt)), via_xy, st_coordinates(st_sfc(d_pt)))
    via <- st_sf(label = "via node", geometry = st_sfc(st_point(via_xy), crs = 3116))
  } else {
    vg <- geom_of(banned[2:(length(banned) - 1)])
    via_xy <- do.call(rbind, lapply(vg, function(l) st_coordinates(l)[, 1:2]))
    banned_xy <- rbind(st_coordinates(st_sfc(o_pt)), via_xy, st_coordinates(st_sfc(d_pt)))
    via <- st_sf(label = "via way", geometry = vg)
  }
  banned_ln <- st_sf(geometry = st_sfc(st_linestring(banned_xy), crs = 3116))
  od <- st_sf(label = c("O", "D"), role = c("origin", "destination"),
              geometry = st_sfc(o_pt, d_pt, crs = 3116))

  # numbered steps along the route to show the direction of travel
  cum <- cumsum(e$length_m[match(path, e$edge_id)]); cum <- cum / cum[length(cum)]
  at <- c(0.2, 0.4, 0.6, 0.8)
  steps <- do.call(c, lapply(at, function(a) {
    k <- which(cum >= a)[1]; prev <- if (k == 1) 0 else cum[k - 1]
    st_sfc(pt_right(geom_of(path[k]), (a - prev) / (cum[k] - prev), off = 0), crs = 3116)
  }))
  steps <- st_sf(label = as.character(seq_along(at)), geometry = steps)

  # Square extent around the route, at least 400 m across
  bb <- st_bbox(st_buffer(st_as_sfc(st_bbox(re)), 100))
  ctr <- c((bb[1] + bb[3]) / 2, (bb[2] + bb[4]) / 2)
  half <- max(200, (bb[3] - bb[1]) / 2, (bb[4] - bb[2]) / 2)
  bb[c(1, 3)] <- ctr[1] + c(-half, half); bb[c(2, 4)] <- ctr[2] + c(-half, half)
  zoom <- if (half > 900) 16 else 17

  m <- tm_shape(st_as_sfc(bb)) + tm_borders(col = NA) +
    tm_basemap("OpenStreetMap", zoom = zoom) +
    tm_shape(re) + tm_lines(col = "#0072B2", lwd = 5, col_alpha = 0.75) +
    tm_shape(via) +
    (if (rel$via_type == "node") tm_symbols(fill = "#E69F00", col = "black", size = 0.9)
     else tm_lines(col = "#E69F00", lwd = 9, col_alpha = 0.9)) +
    tm_shape(banned_ln) + tm_lines(col = "#D55E00", lwd = 3, lty = "dashed") +
    tm_shape(steps) + tm_symbols(fill = "white", col = "#0072B2", size = 0.75) +
    tm_text("label", size = 0.65, fontface = "bold", col = "#0072B2") +
    tm_shape(od) + tm_symbols(fill = "role", size = 1.3, col = "black",
                              fill.scale = tm_scale(values = c(origin = "#009E73", destination = "#CC79A7")),
                              fill.legend = tm_legend_hide()) +
    tm_text("label", size = 0.85, fontface = "bold", col = "white") +
    tm_title(sprintf("r%s  %s (via %s)", id, rel$restriction, rel$via_type), size = 1.1) +
    tm_add_legend(type = "lines", labels = c("Fastest legal route", "Banned direct movement", paste("Via", rel$via_type)),
                  col = c("#0072B2", "#D55E00", "#E69F00"), lty = c("solid", "dashed", "solid"), lwd = c(4, 3, 6)) +
    tm_add_legend(type = "symbols", labels = c("O  origin edge", "D  destination edge", "1-4  direction of travel"),
                  fill = c("#009E73", "#CC79A7", "white"), col = c("black", "black", "#0072B2"), size = 0.8) +
    tm_legend(position = tm_pos_in("right", "top"), bg.color = "white", bg.alpha = 0.85, text.size = 0.6, frame = TRUE) +
    tm_credits(sprintf(paste0("Route: %d edges, %.2f km, %.1f min including turn penalties\n",
                              "Banned direct movement: %.0f s\n",
                              "Basemap \u00a9 OpenStreetMap contributors"),
                       length(path), len_km, cost / 60, direct_s),
               position = c("left", "bottom"), size = 0.6, bg.color = "white", bg.alpha = 0.85) +
    tm_scalebar(position = c("right", "bottom"), text.size = 0.5) +
    tm_layout(frame = TRUE, inner.margins = 0)
  f <- file.path(PLOT_DIR, sprintf("route_r%s.png", id))
  suppressMessages(tmap_save(m, f, width = 1800, height = 1900, dpi = 220))
  cat("saved", f, "| zoom", zoom, "\n")
}
