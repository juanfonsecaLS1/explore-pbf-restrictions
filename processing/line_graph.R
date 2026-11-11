# Line (dual) graph of a directed road network with turn restrictions.
#
# Vertices of the line graph are directed edges of the road network and its
# edges are turns. Compared with igraph::make_line_graph() this
#   * keys everything on unique edge ids (parallel links stay separate),
#   * allows U-turns at dead ends, turning facilities and (optionally)
#     junctions, see uturn_rules(),
#   * deletes turns banned by node-via restrictions, and
#   * handles way-via restrictions by splitting the via vertices, so a banned
#     sequence from -> via ... -> to cannot be driven while each individual
#     turn still can, and
#   * adds OSRM car-profile turn penalties (turn angle, traffic signals,
#     U-turns) to the travel-time weight.

source("processing/osm_network.R")

# km/h, used where maxspeed is missing (urban Colombia: 50 km/h general limit)
default_speeds <- c(
  motorway = 80, motorway_link = 50, trunk = 60, trunk_link = 40,
  primary = 50, primary_link = 40, secondary = 50, secondary_link = 40,
  tertiary = 40, tertiary_link = 30, unclassified = 30, residential = 30,
  living_street = 10, service = 20, track = 15, road = 30
)

# Numeric km/h from maxspeed tags ("50", "30 mph", "50;60" -> first value)
parse_maxspeed <- function(x) {
  x <- sub(";.*", "", x)
  v <- suppressWarnings(as.numeric(sub("^([0-9.]+).*", "\\1", x)))
  mph <- grepl("mph", x, fixed = TRUE)
  v[mph] <- v[mph] * 1.609344
  v
}

# Travel time in seconds for each directed edge
edge_travel_time <- function(edges, speeds = default_speeds) {
  v <- parse_maxspeed(edges$maxspeed)
  fallback <- speeds[edges$highway]
  fallback[is.na(fallback)] <- min(speeds)
  v[is.na(v) | v <= 0] <- fallback[is.na(v) | v <= 0]
  3.6 * edges$length_m / v
}

# ---- Turn penalties (OSRM car profile) ---------------------------------------
# Defaults are those of OSRM's car.lua (seconds). The angle curve was fitted
# by OSRM developers to a small sample of measured turns (osrm-backend
# PR #2849); the U-turn and signal values are not documented as calibrated.
# Use `scale` for sensitivity runs (0 = no penalties).
osrm_turn_params <- function(max_pen = 7.5, bias = 1.075, signal = 2,
                             uturn = 20, lookahead_m = 10, scale = 1,
                             left_hand = FALSE) {
  list(max_pen = max_pen, bias = if (left_hand) 1 / bias else bias,
       signal = signal, uturn = uturn, lookahead_m = lookahead_m,
       scale = scale)
}

# Angle part of the penalty. angle: 0 = straight, + right, - left, 180 = U-turn.
# Only turns at junctions (> 2 roads) and U-turns get it, as in OSRM.
osrm_angle_penalty <- function(angle, n_roads, uturn, p = osrm_turn_params()) {
  pen <- numeric(length(angle))
  use <- (n_roads > 2 | uturn) & !is.na(angle)
  r <- use & angle >= 0
  l <- use & angle < 0
  b <- p$bias
  pen[r] <- p$max_pen / (1 + exp(-((13 / b) * angle[r] / 180 - 6.5 * b)))
  pen[l] <- p$max_pen / (1 + exp(-((13 * b) * -angle[l] / 180 - 6.5 / b)))
  pen
}

# Bearing (degrees clockwise from north) of each edge where it leaves its
# start node and where it arrives at its end node, measured over the first /
# last `d` metres (OSRM also looks a few metres past the junction so short
# digitising segments do not distort the angle)
edge_bearings <- function(edges, d = 10) {
  xy <- sf::st_coordinates(sf::st_geometry(edges))
  k <- xy[, "L1"]
  x <- xy[, "X"]
  y <- xy[, "Y"]
  if (isTRUE(sf::st_is_longlat(edges))) {
    m <- 6371008.8 * pi / 180
    x <- x * cos(y[match(k, k)] * pi / 180) * m
    y <- y * m
  }
  pts <- data.table::data.table(k, x, y)
  # point at distance d along the line (shorter lines: their far end)
  ahead <- function(x, y) {
    s <- c(0, cumsum(sqrt(diff(x)^2 + diff(y)^2)))
    if (s[length(s)] == 0) return(NA_real_)
    dd <- min(d, s[length(s)])
    j <- which(s >= dd)[1]
    f <- (dd - s[j - 1]) / (s[j] - s[j - 1])
    px <- x[j - 1] + f * (x[j] - x[j - 1])
    py <- y[j - 1] + f * (y[j] - y[j - 1])
    atan2(px - x[1], py - y[1]) * 180 / pi
  }
  b <- pts[, .(start = ahead(x, y),
               # leaving the end node backwards, then turned round
               end = (ahead(rev(x), rev(y)) + 180) %% 360), by = k]
  b[match(seq_len(nrow(edges)), b$k), .(start, end)]
}

# ---- U-turn rules ---------------------------------------------------------------
# A U-turn here means driving back along the same two-way link. Turning back
# on a divided road (separate one-way carriageways) is an ordinary sequence of
# turns through the median and is not affected by these rules.
#   * dead ends (no other exit from the edge): always allowed
#   * turning circles, turning loops, mini-roundabouts: always allowed
#   * junctions (>= 3 roads) if `junctions`, except on `exclude_highway` and,
#     unless `at_signals`, at signalised junctions
#   * never in the middle of a road (a node where only 2 links meet)
# no_u_turn relations are applied afterwards and remove U-turns as usual.
# uturn_rules(junctions = FALSE) reproduces the old dead-end-only rule.
turning_facilities <- c("turning_circle", "turning_loop", "mini_roundabout")

uturn_rules <- function(junctions = TRUE, at_signals = FALSE,
                        exclude_highway = c("motorway", "motorway_link",
                                            "trunk", "trunk_link")) {
  list(junctions = junctions, at_signals = at_signals,
       exclude_highway = exclude_highway)
}

# Restrictions that apply to cars, one row per relation with from/via/to
car_restrictions <- function(restrictions) {
  rel <- data.table::as.data.table(restrictions$relations)
  m <- data.table::as.data.table(restrictions$members)
  except <- if ("except" %in% names(rel)) rel$except else NA
  car_exempt <- grepl("(^|;)\\s*(motorcar|motor_vehicle|vehicle)\\s*(;|$)",
                      except)
  rel <- rel[!is.na(restriction) & !car_exempt]
  rel[, kind := sub("_.*", "", restriction)]
  rel <- rel[kind %in% c("no", "only")]

  m <- m[relation_id %in% rel$osm_id]
  data.table::setorder(m, relation_id, seq)
  out <- m[, .(
    from = list(member_id[role == "from" & member_type == "way"]),
    to = list(member_id[role == "to" & member_type == "way"]),
    via = list(member_id[role == "via"]),
    via_type = paste(unique(member_type[role == "via"]), collapse = ";")
  ), by = relation_id]
  out <- merge(out, rel[, .(relation_id = osm_id, restriction, kind)],
               by = "relation_id")
  # Well-formed only: one from, one to, via of a single type
  out[lengths(from) == 1 & lengths(to) == 1 & lengths(via) >= 1 &
        via_type %in% c("node", "way")]
}

#' @param edges        directed edges (from osm_driving_network()$edges)
#' @param restrictions output of oe_restriction_members()
#' @param weights      named list of numeric vectors (one value per edge) to
#'                     turn into line-graph weights; defaults to travel time
#' @param nodes        tagged network nodes (osm_driving_network()$nodes):
#'                     traffic signals and turning circles; NULL = none
#' @param uturns       uturn_rules(): where U-turns are allowed
#' @param turn_penalty osrm_turn_params(), or NULL for no turn penalties
#' @param penalise     names of the weights that get the turn penalty added
#'                     (seconds, so travel-time weights only)
#' @return list with
#'   vertices: vertex_id, edge_id (the road edge it represents; split copies
#'             share the edge_id of the original)
#'   turns:    from, to (vertex ids), turn_type, angle, n_roads, signal,
#'             penalty (s), and one column per weight: mean of the two edges'
#'             values (as in the MRN workflow), plus `penalty` for the
#'             weights in `penalise`
#'   report:   restriction-matching summary
build_line_graph <- function(edges, restrictions,
                             weights = list(tt = edge_travel_time(edges)),
                             nodes = NULL,
                             uturns = uturn_rules(),
                             turn_penalty = osrm_turn_params(),
                             penalise = "tt") {
  e <- data.table::as.data.table(sf::st_drop_geometry(edges))
  e[, signals := if ("signals" %in% names(e)) signals else 0L]
  e[, highway := if ("highway" %in% names(e)) highway else NA_character_]
  e <- e[, .(edge_id, link_id, way_id, from, to, signals, highway)]
  report <- list()

  # 1. All turns: e1 ends where e2 starts
  turns <- merge(
    e[, .(e1 = edge_id, node = to, l1 = link_id, n1 = from)],
    e[, .(e2 = edge_id, node = from, l2 = link_id, n2 = to)],
    by = "node", allow.cartesian = TRUE
  )
  turns <- turns[e1 != e2]

  # 2. U-turns (back along the same two-way link): see uturn_rules()
  nd <- if (is.null(nodes)) data.frame(node = numeric(), highway = character())
        else nodes
  signal_nodes <- nd$node[nd$highway == "traffic_signals"]
  facility_nodes <- nd$node[nd$highway %in% turning_facilities]
  # roads meeting at a node: links, counted once whatever their direction
  node_links <- unique(rbind(e[, .(node = from, link_id)],
                             e[, .(node = to, link_id)]))[, .N, by = node]
  turns[, uturn := l1 == l2]
  turns[, n_out := sum(!uturn), by = e1]
  turns[, n_roads := node_links$N[match(node, node_links$node)]]
  turns[, hw := e$highway[match(e1, e$edge_id)]]
  turns[, turn_type := data.table::fcase(
    !uturn, "turn",
    n_out == 0, "dead_end_uturn",
    node %in% facility_nodes, "facility_uturn",
    uturns$junctions & n_roads >= 3 & !hw %in% uturns$exclude_highway &
      (uturns$at_signals | !node %in% signal_nodes), "junction_uturn",
    default = NA_character_
  )]
  turns <- turns[!is.na(turn_type)]
  report$uturns <- table(turns$turn_type[turns$uturn])

  rx <- car_restrictions(restrictions)
  report$restrictions <- nrow(rx)

  # 3. Node-via restrictions: delete turns
  rn <- rx[via_type == "node" & lengths(via) == 1]
  # built from columns so an area with no via-node restrictions still gets
  # an empty table with the right columns
  rn <- data.table::data.table(
    relation_id = rn$relation_id, kind = rn$kind,
    from = as.numeric(unlist(rn$from)), via = as.numeric(unlist(rn$via)),
    to = as.numeric(unlist(rn$to))
  )
  rn <- merge(rn, e[, .(from = way_id, via = to, e1 = edge_id, l1 = link_id)],
              by = c("from", "via"))
  rn <- merge(rn, e[, .(to = way_id, via = from, e2 = edge_id, l2 = link_id)],
              by = c("to", "via"))
  # When from and to are the same way and it passes through the via node,
  # a U-turn restriction means going back along the same link, not on
  rn <- rn[!(restriction_is_uturn(relation_id, rx) & from == to & l1 != l2)]
  report$node_via <- sum(rx$via_type == "node")
  report$node_via_matched <- data.table::uniqueN(rn$relation_id)

  banned <- rbind(
    rn[kind == "no", .(e1, e2)],
    # only_*: every other turn out of the from edge at that node is banned
    merge(rn[kind == "only", .(e1, keep = e2, node = via)],
          turns[, .(e1, e2, node)], by = c("e1", "node"))[e2 != keep, .(e1, e2)]
  )
  n0 <- nrow(turns)
  turns <- turns[!banned, on = c("e1", "e2")]
  report$turns_removed_node_via <- n0 - nrow(turns)

  # Line graph as an edge list on vertex ids (copies are added below);
  # `fe` caches the road edge of `from` for the via-path searches
  lg <- turns[, .(from = e1, to = e2, turn_type, fe = e1)]
  vertices <- data.table::data.table(vertex_id = e$edge_id, edge_id = e$edge_id)

  # 4. Way-via restrictions: find the edge sequence and split the via vertices.
  #    Searches match vertices by road edge, so copies made by other
  #    restrictions are found too. Passes repeat until nothing changes, so the
  #    order in which overlapping restrictions are applied does not matter.
  rw <- rx[via_type == "way"]
  report$way_via <- nrow(rw)
  report$way_via_paths <- list()
  sets <- lapply(seq_len(nrow(rw)), function(i) list(
    from = e$edge_id[e$way_id %in% rw$from[[i]]],
    via = e$edge_id[e$way_id %in% rw$via[[i]]],
    to = e$edge_id[e$way_id %in% rw$to[[i]]]
  ))
  in_net <- vapply(seq_len(nrow(rw)), function(i) {
    all(rw$from[[i]] %in% e$way_id, rw$via[[i]] %in% e$way_id,
        rw$to[[i]] %in% e$way_id)
  }, logical(1))
  splits <- integer(nrow(rw))
  only_done <- vector("list", nrow(rw))  # from vertices already handled

  for (pass in 1:10) {
    changed <- FALSE
    for (i in which(in_net)) {
      s <- sets[[i]]
      id <- as.character(rw$relation_id[i])
      repeat {
        if (rw$kind[i] == "only") {
          # one split per from vertex (original or copy)
          verts <- unique(c(lg$from, lg$to))
          starts <- setdiff(verts[orig_id(verts) %in% s$from], only_done[[i]])
          if (!length(starts)) break
          only_done[[i]] <- c(only_done[[i]], starts[1])
          path <- find_via_path(lg, s$from, s$via, s$to, start = starts[1])
          if (is.null(path)) next
        } else {
          # keep splitting until no vertex can drive the banned sequence
          path <- find_via_path(lg, s$from, s$via, s$to)
          if (is.null(path)) break
        }
        lg <- split_via(lg, path, rw$kind[i], rw$relation_id[i])
        splits[i] <- splits[i] + 1L
        changed <- TRUE
        if (is.null(report$way_via_paths[[id]])) {
          report$way_via_paths[[id]] <- orig_id(path)
        }
        if (splits[i] >= 50L) {
          warning("Restriction ", id, " needed 50 splits; stopping it here")
          break
        }
      }
    }
    if (!changed) break
  }

  # Final check: no banned sequence may still be drivable
  violated <- vapply(seq_len(nrow(rw)), function(i) {
    if (!in_net[i] || rw$kind[i] != "no") return(FALSE)
    !is.null(find_via_path(lg, sets[[i]]$from, sets[[i]]$via, sets[[i]]$to))
  }, logical(1))
  report$way_via_status <- data.table::data.table(
    relation_id = rw$relation_id,
    restriction = rw$restriction,
    splits = splits,
    status = data.table::fcase(
      !in_net, "way not in network",
      violated, "VIOLATED",
      splits > 0L, "applied",
      default = "already impossible"
    )
  )
  if (any(violated)) {
    warning(sum(violated), " way-via restriction(s) still drivable: ",
            toString(rw$relation_id[violated]))
  }
  report$way_via_applied <- sum(splits > 0L)
  lg[, fe := NULL]
  copies <- unique(c(lg$from, lg$to))
  copies <- copies[!copies %in% vertices$vertex_id]
  vertices <- rbind(vertices, data.table::data.table(
    vertex_id = copies, edge_id = orig_id(copies)
  ))

  # 5. Turn attributes and penalties (split copies share their original's)
  idx_from <- match(orig_id(lg$from), e$edge_id)
  idx_to <- match(orig_id(lg$to), e$edge_id)
  via_node <- e$to[idx_from]
  uturn <- grepl("uturn$", lg$turn_type)
  lg[, n_roads := node_links$N[match(via_node, node_links$node)]]
  lg[, signal := via_node %in% signal_nodes]

  if (is.null(turn_penalty)) {
    lg[, `:=`(angle = NA_real_, penalty = 0)]
  } else {
    p <- turn_penalty
    b <- edge_bearings(edges, p$lookahead_m)
    a <- ((b$start[idx_to] - b$end[idx_from] + 180) %% 360) - 180
    a[uturn] <- 180
    facility <- via_node %in% facility_nodes
    lg[, angle := round(a, 1)]
    lg[, penalty := p$scale * (
      osrm_angle_penalty(a, n_roads, uturn, p) +
        p$signal * signal +                        # signal at the junction
        p$uturn * (uturn & !facility) +            # none at turning circles
        # signals inside the edges: half of each, like the edge weights
        0.5 * p$signal * (e$signals[idx_from] + e$signals[idx_to])
    )]
    report$penalty <- summary(lg$penalty)
  }

  # 6. Turn weights: mean of the two road edges (+ penalty)
  for (w in names(weights)) {
    lg[[w]] <- 0.5 * weights[[w]][idx_from] + 0.5 * weights[[w]][idx_to]
    if (w %in% penalise) lg[[w]] <- lg[[w]] + lg$penalty
  }

  list(vertices = vertices, turns = lg, report = report)
}

# igraph object of the line graph (vertex names = vertex ids)
line_graph_igraph <- function(lg, weight = "tt") {
  g <- igraph::graph_from_data_frame(
    data.frame(from = lg$turns$from, to = lg$turns$to,
               weight = lg$turns[[weight]]),
    directed = TRUE, vertices = data.frame(name = lg$vertices$vertex_id)
  )
  g
}

#' Fastest route between two road edges
#' @param g         output of line_graph_igraph()
#' @param from_edge,to_edge road edge ids (`edge_id`)
#' @return character vector of road edge ids (split copies collapsed), or
#'   NULL if unreachable. The cost is in attr(, "cost").
route_edges <- function(g, from_edge, to_edge) {
  vn <- igraph::V(g)$name
  targets <- which(sub("#.*", "", vn) == to_edge)  # original + copies
  d <- igraph::distances(g, v = from_edge, to = targets, mode = "out")
  if (all(is.infinite(d))) return(NULL)
  best <- targets[which.min(d)]
  p <- igraph::shortest_paths(g, from = from_edge, to = best, mode = "out",
                              output = "vpath")$vpath[[1]]
  out <- sub("#.*", "", names(p))
  attr(out, "vertices") <- names(p)
  attr(out, "cost") <- min(d)
  out
}

# TRUE for relations tagged no_u_turn / only_u_turn
restriction_is_uturn <- function(relation_id, rx) {
  grepl("_u_turn$", rx$restriction[match(relation_id, rx$relation_id)])
}

# Road edge a line-graph vertex stands for ("e#r123#r456" -> "e")
orig_id <- function(v) sub("#.*", "", v)

# A drivable vertex sequence from -> via edge(s) -> to, following the current
# turns. Vertices are matched by road edge, so split copies count.
# Returns the vertex ids c(v0, v1..vm, v_to) or NULL.
find_via_path <- function(lg, from_edges, via_edges, to_edges, start = NULL,
                          max_len = 12) {
  # only turns that can be part of the sequence
  rel <- lg[lg$fe %in% c(from_edges, via_edges)]
  out_by <- split(rel$to, rel$from)
  if (is.null(start)) {
    start <- unique(rel$from[orig_id(rel$from) %in% from_edges])
  }
  queue <- lapply(start, c)
  while (length(queue)) {
    p <- queue[[1]]
    queue <- queue[-1]
    nxt <- out_by[[p[length(p)]]]
    if (length(p) > 1) {
      hit <- nxt[orig_id(nxt) %in% to_edges]
      if (length(hit)) return(c(p, hit[1]))
    }
    if (length(p) >= max_len) next
    nxt <- nxt[orig_id(nxt) %in% via_edges & !orig_id(nxt) %in% orig_id(p)]
    queue <- c(queue, lapply(nxt, function(v) c(p, v)))
  }
  NULL
}

# Split the via vertices of one restriction. path = c(v0, v1..vm, v_to)
split_via <- function(lg, path, kind, relation_id) {
  n <- length(path)
  v0 <- path[1]
  via <- path[2:(n - 1)]
  vto <- path[n]
  copy <- paste0(via, "#r", relation_id)
  existing <- unique(lg$from)
  new <- list()

  for (k in seq_along(via)) {
    # already split for this restriction from another entry: reuse the copy
    if (copy[k] %in% existing) next
    nxt_orig <- if (k < length(via)) via[k + 1] else vto
    nxt_new <- if (k < length(via)) copy[k + 1] else vto
    out <- lg[lg$from == via[k]]  # vector scan: no index rebuilt after each split
    if (kind == "only") {
      # the copy may only continue along the sequence
      new[[k]] <- out[to == nxt_orig][, `:=`(from = copy[k], to = nxt_new)]
    } else {
      # the copy keeps every exit except the banned continuation;
      # inside the chain it continues to the next copy instead
      keep <- out[to != nxt_orig][, from := copy[k]]
      if (k < length(via)) {
        keep <- rbind(keep, out[to == nxt_orig][, `:=`(from = copy[k],
                                                        to = nxt_new)])
      }
      new[[k]] <- keep
    }
  }

  # Entering the chain from v0 now goes to the first copy
  if (kind == "only") {
    lg <- lg[lg$from != v0 | lg$to == via[1]]
  }
  lg[lg$from == v0 & lg$to == via[1], to := copy[1]]
  new <- data.table::rbindlist(new)
  if (nrow(new)) new[, fe := orig_id(from)]
  rbind(lg, new)
}
