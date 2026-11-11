# Toy networks and plotting helpers for the turn-restriction schematics.
#
# Toy networks are planar (no CRS) sfnetwork objects. Directed edges carry the
# same columns as osm_driving_network()$edges, so build_line_graph() from
# line_graph.R runs on them unchanged.

library(sf)
library(sfnetworks)
library(ggplot2)
library(patchwork)
library(data.table)
source("processing/line_graph.R")

# Okabe-Ito colours: readable for colour-blind viewers and on white
pal <- c(
  base = "grey55", from = "#0072B2", via = "#E69F00", to = "#CC79A7",
  banned = "#D55E00", copy = "#009E73", invented = "#D55E00"
)

# Shift a 2-point line to its right (right-hand traffic) and trim its ends so
# arrowheads stop short of the nodes
lane <- function(p1, p2, shift = 0.06, trim = 0.12) {
  d <- p2 - p1
  len <- sqrt(sum(d^2))
  u <- d / len
  n <- c(u[2], -u[1])
  a <- p1 + n * shift + u * min(trim, len / 3)
  b <- p2 + n * shift - u * min(trim, len / 3)
  st_linestring(rbind(a, b))
}

# Quadratic curve from p1 to p2 pulled towards `ctrl` (the junction node), so
# a turn is drawn passing through the intersection where it happens
curve <- function(p1, p2, ctrl, pull = 0.8, n = 24) {
  ctrl <- (p1 + p2) / 2 + pull * (ctrl - (p1 + p2) / 2)
  t <- seq(0, 1, length.out = n)
  st_linestring(cbind(
    (1 - t)^2 * p1[1] + 2 * (1 - t) * t * ctrl[1] + t^2 * p2[1],
    (1 - t)^2 * p1[2] + 2 * (1 - t) * t * ctrl[2] + t^2 * p2[2]
  ))
}

#' Directed toy network
#' @param nodes data.frame(id, x, y, label)
#' @param links data.frame(way, a, b, oneway, fwd, bwd): fwd/bwd are the
#'              labels of the a->b and b->a edges (bwd ignored if oneway)
toy_network <- function(nodes, links) {
  xy <- as.matrix(nodes[, c("x", "y")])
  rownames(xy) <- nodes$id
  rows <- list()
  geoms <- list()
  for (i in seq_len(nrow(links))) {
    l <- links[i, ]
    two <- !l$oneway
    pa <- xy[as.character(l$a), ]
    pb <- xy[as.character(l$b), ]
    rows[[length(rows) + 1]] <- data.frame(
      edge_id = l$fwd, link_id = paste(l$way, i, sep = "_"), way_id = l$way,
      from = l$a, to = l$b, label = l$fwd
    )
    geoms[[length(geoms) + 1]] <- lane(pa, pb, shift = if (two) 0.06 else 0)
    if (two) {
      rows[[length(rows) + 1]] <- data.frame(
        edge_id = l$bwd, link_id = paste(l$way, i, sep = "_"), way_id = l$way,
        from = l$b, to = l$a, label = l$bwd
      )
      geoms[[length(geoms) + 1]] <- lane(pb, pa)
    }
  }
  edges <- st_sf(do.call(rbind, rows), geometry = st_sfc(geoms))
  node_sf <- st_as_sf(nodes, coords = c("x", "y"))
  # edge geometries are offset lanes, so their ends do not touch the nodes;
  # force = TRUE keeps the from/to topology without that check
  edges$from <- match(edges$from, nodes$id)
  edges$to <- match(edges$to, nodes$id)
  net <- sfnetwork(node_sf, edges, directed = TRUE, force = TRUE)
  net
}

# Edges as the line-graph code expects them (OSM node ids in from/to)
toy_edges <- function(net) {
  n <- st_as_sf(net, "nodes")
  e <- st_as_sf(net, "edges")
  e$from <- n$id[e$from]
  e$to <- n$id[e$to]
  e
}

#' Restrictions in the format of oe_restriction_members()
#' @param ... lists(id, restriction, from, via, to, via_type)
toy_restrictions <- function(...) {
  r <- list(...)
  list(
    relations = data.frame(
      osm_id = vapply(r, `[[`, numeric(1), "id"),
      type = "restriction",
      restriction = vapply(r, `[[`, character(1), "restriction"),
      except = NA_character_
    ),
    members = do.call(rbind, lapply(r, function(x) data.frame(
      relation_id = x$id, seq = 1:3,
      member_type = c("way", x$via_type, "way"),
      member_id = c(x$from, x$via, x$to),
      role = c("from", "via", "to"), lon = NA_real_, lat = NA_real_
    )))
  )
}

# Line graph as an sfnetwork: vertices at the middle of their road edge
# (copies drawn further out to the right), turns as curves
line_graph_sfnet <- function(net, lg) {
  e <- st_as_sf(net, "edges")
  mid <- st_coordinates(st_line_sample(st_geometry(e), sample = 0.5))[, 1:2]
  rownames(mid) <- e$edge_id
  v <- lg$vertices
  depth <- lengths(regmatches(v$vertex_id, gregexpr("#", v$vertex_id)))
  base <- mid[v$edge_id, , drop = FALSE]
  # unit normal (to the right) of each road edge, for placing copies
  nrm <- t(vapply(e$edge_id, function(id) {
    xy <- st_coordinates(st_geometry(e)[e$edge_id == id])[, 1:2]
    d <- xy[nrow(xy), ] - xy[1, ]
    d <- d / sqrt(sum(d^2))
    c(d[2], -d[1])
  }, numeric(2)))
  pos <- base + nrm[v$edge_id, , drop = FALSE] * 0.16 * depth
  lab <- e$label[match(v$edge_id, e$edge_id)]
  v_sf <- st_sf(
    vertex_id = v$vertex_id, edge_id = v$edge_id,
    label = paste0(lab, strrep("′", depth)), copy = depth > 0,
    geometry = st_sfc(lapply(seq_len(nrow(pos)), function(i) st_point(pos[i, ])))
  )
  t <- lg$turns
  rownames(pos) <- v$vertex_id
  # junction of each turn = node where the from-edge ends
  nd <- st_as_sf(net, "nodes")
  nxy <- st_coordinates(nd)
  junction <- nxy[e$to[match(orig_id(t$from), e$edge_id)], , drop = FALSE]
  geom <- st_sfc(lapply(seq_len(nrow(t)), function(i) {
    curve(pos[t$from[i], ], pos[t$to[i], ], junction[i, ])
  }))
  t_sf <- st_sf(from = match(t$from, v$vertex_id), to = match(t$to, v$vertex_id),
                from_id = t$from, to_id = t$to, turn_type = t$turn_type,
                geometry = geom)
  sfnetwork(v_sf, t_sf, directed = TRUE, force = TRUE)
}

# All turns of a line graph before any restriction, for showing what was
# removed (same rules as build_line_graph() steps 1-2 when uturns = FALSE)
all_turns <- function(net, uturns = TRUE) {
  e <- as.data.table(st_drop_geometry(toy_edges(net)))
  t <- merge(e[, .(from = edge_id, node = to, l1 = link_id)],
             e[, .(to = edge_id, node = from, l2 = link_id)],
             by = "node", allow.cartesian = TRUE)[from != to]
  t[, turn_type := fifelse(l1 == l2, "uturn", "turn")]
  if (!uturns) t <- t[turn_type == "turn"]
  list(vertices = data.table(vertex_id = e$edge_id, edge_id = e$edge_id),
       turns = t[, .(from, to, turn_type)])
}

# ---- Plotting ---------------------------------------------------------------

arr <- arrow(length = unit(2.2, "mm"), type = "closed")
fill_pal <- replace(pal, "base", "white")

# Label positions: middle of each line, pushed `shift` to its right
mid_points <- function(lines, shift = 0) {
  st_sfc(lapply(st_geometry(lines), function(l) {
    xy <- st_coordinates(l)[, 1:2]
    d <- xy[nrow(xy), ] - xy[1, ]
    d <- d / sqrt(sum(d^2))
    st_point((xy[1, ] + xy[nrow(xy), ]) / 2 + c(d[2], -d[1]) * shift)
  }))
}

theme_schematic <- function() {
  theme_void(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 10),
          plot.subtitle = element_text(size = 8.5, colour = "grey30"),
          legend.position = "none",
          plot.margin = margin(4, 4, 4, 4))
}

#' Primal network
#' @param roles named character: edge label -> "from" / "via" / "to" / ...
#' @param node_roles named character: node label -> role (e.g. via node)
plot_primal <- function(net, title, subtitle = NULL, roles = character(),
                        node_roles = character(), dim_others = FALSE) {
  e <- st_as_sf(net, "edges")
  n <- st_as_sf(net, "nodes")
  e$role <- ifelse(e$label %in% names(roles), roles[e$label], "base")
  n$role <- ifelse(n$label %in% names(node_roles), node_roles[n$label], "base")
  lab <- st_sf(label = e$label, role = e$role, geometry = mid_points(e, 0.13))
  ggplot() +
    geom_sf(data = e, aes(colour = role, linewidth = role != "base"),
            arrow = arr, alpha = if (dim_others) 0.9 else 1) +
    geom_sf(data = n, aes(fill = role), shape = 21, size = 3.2,
            colour = "grey20", stroke = 0.4) +
    geom_sf_text(data = n, aes(label = label), size = 2.8, fontface = "bold",
                 nudge_x = -0.14, nudge_y = 0.14, colour = "grey15") +
    geom_sf_text(data = lab, aes(label = label, colour = role), size = 2.4) +
    scale_colour_manual(values = pal) +
    scale_fill_manual(values = fill_pal) +
    scale_linewidth_manual(values = c(`FALSE` = 0.45, `TRUE` = 1.1)) +
    labs(title = title, subtitle = subtitle) +
    coord_sf(datum = NA, expand = TRUE) +
    theme_schematic()
}

#' Line graph drawn over the faint primal network
#' @param removed turns (from_id, to_id) to draw dashed in the banned colour
#' @param highlight named character: vertex label -> role
plot_line_graph <- function(net, lgnet, title, subtitle = NULL,
                            removed = NULL, highlight = character(),
                            invented = NULL) {
  e <- st_as_sf(net, "edges")
  v <- st_as_sf(lgnet, "nodes")
  t <- st_as_sf(lgnet, "edges")
  t$kind <- ifelse(t$turn_type == "turn", "turn", "uturn")
  v$role <- ifelse(v$copy, "copy",
                   ifelse(v$label %in% names(highlight), highlight[v$label], "base"))
  p <- ggplot() +
    geom_sf(data = e, colour = "grey88", linewidth = 2.2) +
    geom_sf(data = st_as_sf(net, "nodes"), colour = "grey60", size = 1.4) +
    geom_sf(data = t, aes(linetype = kind), colour = "grey35",
            linewidth = 0.35, arrow = arrow(length = unit(1.6, "mm"), type = "closed"))
  if (!is.null(removed) && nrow(removed)) {
    p <- p + geom_sf(data = removed, colour = pal["banned"], linewidth = 0.6,
                     linetype = "22",
                     arrow = arrow(length = unit(1.8, "mm"), type = "closed"))
  }
  if (!is.null(invented) && nrow(invented)) {
    p <- p + geom_sf(data = invented, colour = pal["invented"], linewidth = 0.7,
                     arrow = arrow(length = unit(1.8, "mm"), type = "closed"))
  }
  p +
    geom_sf(data = v, aes(fill = role), shape = 21, size = 5.2,
            colour = "grey20", stroke = 0.35) +
    geom_sf_text(data = v, aes(label = label), size = 2.1, fontface = "bold") +
    scale_fill_manual(values = fill_pal) +
    scale_linetype_manual(values = c(turn = "solid", uturn = "41")) +
    labs(title = title, subtitle = subtitle) +
    coord_sf(datum = NA, expand = TRUE) +
    theme_schematic()
}

# Add a parallel copy of a road edge (GraphHopper's artificial edge)
add_edge_copy <- function(net, lab, new_lab, shift = 0.14) {
  n <- st_as_sf(net, "nodes")
  e <- st_as_sf(net, "edges")
  xy <- st_coordinates(n)
  ce <- e[e$label == lab, ]
  ce$label <- new_lab
  ce$edge_id <- new_lab
  st_geometry(ce) <- st_sfc(lane(xy[ce$from, ], xy[ce$to, ], shift = shift))
  sfnetwork(n, rbind(e, ce), directed = TRUE, force = TRUE)
}

# Turns drawn on the primal network: from the end of one edge, through the
# junction, to the start of the next. turns = data.frame(from, to, label)
turn_marks <- function(net, turns) {
  e <- st_as_sf(net, "edges")
  nxy <- st_coordinates(st_as_sf(net, "nodes"))
  geoms <- lapply(seq_len(nrow(turns)), function(i) {
    fe <- e[e$label == turns$from[i], ]
    te <- e[e$label == turns$to[i], ]
    # start part-way along the from-edge, end part-way along the to-edge
    a <- st_coordinates(st_line_sample(st_geometry(fe), sample = 0.6))[1, 1:2]
    b <- st_coordinates(st_line_sample(st_geometry(te), sample = 0.4))[1, 1:2]
    curve(a, b, nxy[fe$to, ], pull = 1)
  })
  lines <- st_sf(label = turns$label, geometry = st_sfc(geoms))
  # labels near the end of each turn, so turns from one edge do not overlap
  pts <- st_sf(label = turns$label,
               geometry = st_line_sample(st_geometry(lines), sample = 0.85) |>
                 st_cast("POINT"))
  list(lines = lines, labels = pts)
}

# Primal plot with turns marked and labelled with their cost. An optional
# `turns$kind` column ("forbidden" / "allowed") sets the colour; default is
# forbidden (infinite cost)
plot_turn_costs <- function(net, turns, title, subtitle = NULL, roles = character(),
                            node_roles = character()) {
  if (is.null(turns$kind)) turns$kind <- "forbidden"
  m <- turn_marks(net, turns)
  m$lines$kind <- turns$kind
  m$labels$kind <- turns$kind
  kind_col <- c(forbidden = unname(pal["banned"]), allowed = "grey30")
  plot_primal(net, title, subtitle, roles = roles, node_roles = node_roles) +
    geom_sf(data = m$lines, aes(linetype = kind), colour = kind_col[m$lines$kind],
            linewidth = 0.7, arrow = arrow(length = unit(1.8, "mm"), type = "closed")) +
    geom_sf_label(data = m$labels, aes(label = label), colour = kind_col[m$labels$kind],
                  fill = "white", size = 4.2, fontface = "bold", linewidth = 0,
                  label.padding = unit(0.4, "mm")) +
    scale_linetype_manual(values = c(forbidden = "22", allowed = "solid"))
}

# Turns present in `before` but not in `after`, drawn with `before`'s geometry
removed_turns <- function(lg_before_net, lg_after) {
  t <- st_as_sf(lg_before_net, "edges")
  keep <- !paste(t$from_id, t$to_id) %in% paste(lg_after$turns$from, lg_after$turns$to)
  t[keep, ]
}
