# Build a directed, junction-to-junction driving network from an OSM .pbf
# file, keeping the OSM node ids at every junction so turn restrictions can be
# matched by id (as with the OS MRN).
#
# osmium does the heavy lifting (clipping the extract and attaching node
# coordinates to ways); R only parses the resulting OPL text. Ways are split
# at every node shared with another way, so no spatial subdivision or
# smoothing is needed afterwards.

source("processing/restriction_members.R")

# Value of one tag for each OPL tag string (NA if absent)
opl_tag <- function(tag_str, key) {
  pat <- paste0("(^|,)", key, "=([^,]*)")
  out <- rep(NA_character_, length(tag_str))
  hit <- grepl(pat, tag_str)
  out[hit] <- opl_unescape(sub(paste0(".*", pat, ".*"), "\\2", tag_str[hit]))
  out
}

# Same idea as oe_get_network(mode = "driving") in osmextract
driving_filter <- function(tg) {
  excluded <- c(
    "abandoned", "bus_guideway", "byway", "construction", "corridor",
    "elevator", "fixme", "escalator", "gallop", "historic", "no", "planned",
    "platform", "proposed", "cycleway", "pedestrian", "bridleway", "path",
    "footway", "steps", "raceway", "busway", "bus_stop"
  )
  !is.na(tg$highway) &
    !tg$highway %in% excluded &
    (is.na(tg$access) | !tg$access %in% c("private", "no") |
       tg$motor_vehicle %in% c("yes", "permissive", "designated")) &
    !tg$motor_vehicle %in% c("no", "private") &
    (is.na(tg$service) | !grepl("^private", tg$service))
}

# yes / -1 / no, applying the implied one-ways
normalise_oneway <- function(oneway, highway, junction) {
  out <- rep("no", length(oneway))
  out[oneway %in% c("yes", "true", "1")] <- "yes"
  out[oneway %in% c("-1", "reverse")] <- "-1"
  implied <- is.na(oneway) &
    (junction %in% c("roundabout", "circular") | highway %in% "motorway")
  out[implied] <- "yes"
  out
}

#' @param pbf      path to the .osm.pbf file
#' @param boundary sf/sfc polygon or lon/lat bbox (see oe_restriction_members)
#' @param osmium   command used to call osmium (character vector)
#' @param tags     way tags to keep as columns
#' @return list with
#'   edges:        sf, one row per directed edge (edge_id, link_id, way_id,
#'                 from, to, oneway, tags, length_m, signals = traffic
#'                 signals passed inside the edge in its direction)
#'   nodes:        network nodes tagged highway = traffic_signals,
#'                 turning_circle, turning_loop or mini_roundabout
#'                 (node, highway, direction)
#'   restrictions: output of oe_restriction_members() for the same area
#'   pbf:          path to the clipped extract (temporary file)
osm_driving_network <- function(
    pbf,
    boundary,
    osmium = "osmium",
    tags = c("highway", "name", "oneway", "junction", "access",
             "motor_vehicle", "service", "maxspeed", "lanes")) {
  stopifnot(file.exists(pbf))
  tags <- union(tags, c("highway", "oneway", "junction", "access",
                        "motor_vehicle", "service"))

  # 1. Clip once; complete ways that touch the boundary are kept whole
  area_pbf <- tempfile(fileext = ".osm.pbf")
  run_osmium(osmium, c("extract", osmium_area_args(boundary),
                       shQuote(normalizePath(pbf))), area_pbf)

  # 2. Highways with their node coordinates inlined:
  #    "w1 Thighway=primary,... Nn10x-74.1y4.7,n11x-74.1y4.7"
  hw_pbf <- tempfile(fileext = ".osm.pbf")
  opl <- tempfile(fileext = ".opl")
  on.exit(unlink(c(hw_pbf, opl)), add = TRUE)
  run_osmium(osmium, c("tags-filter", shQuote(area_pbf), "w/highway"), hw_pbf)
  run_osmium(osmium, c("add-locations-to-ways", shQuote(hw_pbf),
                       "-f", "opl,add_metadata=false"), opl)

  lines <- readLines(opl, encoding = "UTF-8")
  # Tagged nodes (traffic signals, turning circles) are kept by
  # add-locations-to-ways: "n123 Thighway=traffic_signals x-74.1 y4.7"
  nlines <- lines[startsWith(lines, "n")]
  lines <- lines[startsWith(lines, "w")]
  way_id <- as.numeric(sub("^w(\\d+) .*", "\\1", lines))
  tag_str <- ifelse(grepl(" T", lines), sub(".* T(\\S*).*", "\\1", lines), "")
  nd_str <- sub(".* N(\\S*).*", "\\1", lines)

  tg <- as.data.frame(sapply(tags, opl_tag, tag_str = tag_str,
                             simplify = FALSE), check.names = FALSE)
  keep <- driving_filter(tg)
  way_id <- way_id[keep]
  tg <- tg[keep, , drop = FALSE]
  nd <- strsplit(nd_str[keep], ",", fixed = TRUE)

  # 3. Long table of way nodes
  n_nd <- lengths(nd)
  nd <- unlist(nd)
  pts <- data.table::data.table(
    way = rep(seq_along(way_id), n_nd),
    pos = sequence(n_nd),
    n = n_nd[rep(seq_along(way_id), n_nd)],
    node = as.numeric(sub("^n(\\d+)x.*", "\\1", nd)),
    x = as.numeric(sub(".*x([-0-9.]+)y.*", "\\1", nd)),
    y = as.numeric(sub(".*y([-0-9.]+)$", "\\1", nd))
  )
  pts <- pts[!is.na(x)]  # nodes outside the extract (none with complete_ways)

  # 4. Split at way ends and at nodes used more than once
  #    (shared with another way, or a way that touches itself)
  pts[, uses := .N, by = node]
  pts[, split := pos == 1L | pos == n | uses > 1L]
  pts[, seg := cumsum(split), by = way]
  # A split node ends the previous segment and starts the next one
  ends <- pts[split & pos > 1L][, seg := seg - 1L]
  pts <- data.table::rbindlist(list(pts[pos < n], ends))
  data.table::setorder(pts, way, seg, pos)
  pts[, link := paste(way_id[way], seg, sep = "_")]

  # Node tags relevant to turn costs, for nodes on the network only
  ntag <- ifelse(grepl(" T", nlines), sub(".* T(\\S*).*", "\\1", nlines), "")
  nodes <- data.table::data.table(
    node = as.numeric(sub("^n(\\d+) .*", "\\1", nlines)),
    highway = opl_tag(ntag, "highway"),
    # direction of a signal relative to its way (forward / backward)
    direction = data.table::fcoalesce(opl_tag(ntag, "traffic_signals:direction"),
                                      opl_tag(ntag, "direction"))
  )
  nodes <- nodes[highway %in% c("traffic_signals", "turning_circle",
                                "turning_loop", "mini_roundabout") &
                   node %in% pts$node]
  nodes[!direction %in% c("forward", "backward"), direction := NA]

  # Signals in the middle of a link are passed by every vehicle on it, so
  # they belong to the edge; signals at link ends are handled per turn
  sig <- nodes[highway == "traffic_signals"]
  mid <- pts[!split & node %in% sig$node]
  mid[, direction := sig$direction[match(node, sig$node)]]
  link_sig <- mid[, .(sig_any = sum(is.na(direction)),
                      sig_fwd = sum(direction %in% "forward"),
                      sig_bwd = sum(direction %in% "backward")), by = link]

  # 5. One row per undirected link
  pts[, `:=`(dx = x - data.table::shift(x), dy = y - data.table::shift(y)),
      by = link]
  links <- pts[, .(
    way = way[1],
    from = node[1],
    to = node[.N],
    # equirectangular approximation, fine for short segments
    length_m = sum(sqrt((dx * cos(y * pi / 180))^2 + dy^2), na.rm = TRUE) *
      6371008.8 * pi / 180
  ), by = link]
  links <- cbind(links, tg[links$way, , drop = FALSE])
  links[, way_id := way_id[way]]
  links[, oneway := normalise_oneway(oneway, highway, junction)]
  links <- links[from != to | length_m > 0]
  links <- merge(links, link_sig, by = "link", all.x = TRUE, sort = FALSE)
  for (v in c("sig_any", "sig_fwd", "sig_bwd")) {
    data.table::set(links, which(is.na(links[[v]])), v, 0L)
  }

  geom <- sfheaders::sfc_linestring(as.data.frame(pts[, .(link, x, y)]),
                                    x = "x", y = "y", linestring_id = "link") |>
    sf::st_set_crs(4326)
  geom <- geom[match(links$link, unique(pts$link))]

  # 6. Directed edges; backward edges get reversed geometry
  fwd <- links$oneway %in% c("yes", "no")
  bwd <- links$oneway %in% c("-1", "no")
  edges <- data.table::rbindlist(list(
    links[fwd][, `:=`(edge_id = paste0(link, "_f"),
                      signals = sig_any + sig_fwd)],
    links[bwd][, `:=`(edge_id = paste0(link, "_b"), from = to, to = from,
                      signals = sig_any + sig_bwd)]
  ))
  edge_geom <- c(geom[fwd], sf::st_reverse(geom[bwd]))
  edges[, c("way", "sig_any", "sig_fwd", "sig_bwd") := NULL]
  data.table::setnames(edges, "link", "link_id")
  data.table::setcolorder(edges, c("edge_id", "link_id", "way_id", "from", "to"))
  edges <- sf::st_sf(as.data.frame(edges), geometry = edge_geom)

  list(
    edges = edges,
    nodes = as.data.frame(nodes),
    restrictions = oe_restriction_members(area_pbf, osmium = osmium),
    pbf = area_pbf
  )
}
