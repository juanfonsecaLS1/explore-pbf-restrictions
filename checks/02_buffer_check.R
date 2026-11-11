# Rebuild the network with a 1 km buffer around the study polygon and list the
# edges inside the polygon that are still outside the main strongly connected
# component (i.e. not cut off by the boundary clip), with the highway ways the
# driving filter removed next to them.
# Adds layer cut_off_with_buffer to outputs/routability_examples.gpkg
# (run checks/01_routability_checks.R first: it recreates that file).
source("processing/config.R")

poly <- study_polygon()
buf <- st_buffer(st_transform(poly, 3116), 1000) |> st_transform(4326)
net1 <- osm_driving_network(colombia_pbf(), buf, osmium = OSMIUM)
lg <- build_line_graph(net1$edges, net1$restrictions, nodes = net1$nodes)
g <- line_graph_igraph(lg)
cc <- components(g, "strong")
inbig <- cc$membership == which.max(cc$csize)
ok <- tapply(inbig, sub("#.*", "", V(g)$name), any)
inside <- as.vector(st_within(net1$edges, st_geometry(poly), sparse = FALSE))
still <- net1$edges$edge_id[inside & !ok[net1$edges$edge_id]]
cat("edges inside the polygon still cut off:", length(still), "\n")

# Highway ways removed by the driving filter that touch those edges
hw <- tempfile(fileext = ".osm.pbf"); opl <- tempfile(fileext = ".opl")
run_osmium(OSMIUM, c("tags-filter", shQuote(net1$pbf), "w/highway"), hw)
run_osmium(OSMIUM, c("cat", shQuote(hw), "-t", "way", "-f", "opl,add_metadata=false"), opl)
l <- readLines(opl, encoding = "UTF-8")
tag_str <- sub(".* T(\\S*) N.*", "\\1", l)
tg <- as.data.frame(sapply(c("highway", "access", "motor_vehicle", "service"),
                           opl_tag, tag_str = tag_str, simplify = FALSE))
excl <- !driving_filter(tg)
nodes_of <- strsplit(gsub("n", "", sub(".* N(\\S*)$", "\\1", l[excl])), ",")
e1 <- as.data.table(st_drop_geometry(net1$edges))
bad_nodes <- unique(c(e1[edge_id %in% still, from], e1[edge_id %in% still, to]))
touch <- vapply(nodes_of, function(n) any(as.numeric(n) %in% bad_nodes), logical(1))
cat("removed highway ways touching them:", sum(touch), "\n")
print(sort(table(paste0("highway=", tg$highway[excl][touch],
                        " access=", tg$access[excl][touch],
                        " motor_vehicle=", tg$motor_vehicle[excl][touch])),
           decreasing = TRUE))
st <- e1[edge_id %in% still]
print(table(st$highway, st$oneway))
cat("by component size:\n")
print(table(cc$csize[cc$membership[match(still, V(g)$name)]]))

st_write(net1$edges[net1$edges$edge_id %in% still, ], EXAMPLES_GPKG,
         layer = "cut_off_with_buffer", quiet = TRUE, append = FALSE)
cat("saved layer cut_off_with_buffer in", EXAMPLES_GPKG, "\n")
