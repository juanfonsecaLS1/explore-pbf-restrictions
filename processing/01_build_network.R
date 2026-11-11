# Build the directed driving network and its line graph (turn restrictions,
# U-turn rules and OSM turn penalties) for the study polygon.
# Output: outputs/bogota_network.rds = list(net, lg)
source("processing/config.R")

poly <- study_polygon()
t_net <- system.time(
  net <- osm_driving_network(colombia_pbf(), poly, osmium = OSMIUM)
)
cat("network:", nrow(net$edges), "directed edges in", round(t_net[3], 1), "s\n")
cat("tagged nodes:\n"); print(table(net$nodes$highway))

t_lg <- system.time(
  lg <- build_line_graph(net$edges, net$restrictions, nodes = net$nodes)
)
cat("line graph:", nrow(lg$vertices), "vertices,", nrow(lg$turns), "turns in",
    round(t_lg[3], 2), "s\n")
cat("U-turns:\n"); print(lg$report$uturns)
cat("node-via restrictions matched:", lg$report$node_via_matched, "of",
    lg$report$node_via, "\n")
cat("way-via restrictions:\n"); print(lg$report$way_via_status[, .N, by = status])

saveRDS(list(net = net, lg = lg), NETWORK_RDS)
cat("saved", NETWORK_RDS, "\n")
