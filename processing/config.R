# Shared settings for the processing and checks scripts. Run every script from
# the repository root (the folder with explore-pbf-restrictions.Rproj), e.g.
#   Rscript processing/01_build_network.R
#
# processing/  builds the data
#   restriction_members.R   turn-restriction relations from a .pbf (functions)
#   osm_network.R           directed driving network from a .pbf (functions)
#   line_graph.R            line graph, restrictions, U-turns, turn penalties
#                           (functions)
#   config.R                this file
#   01_build_network.R      network + line graph for the study polygon
# checks/      checks the data and produces results
#   01_routability_checks.R restriction, O-D and connectivity checks
#   02_buffer_check.R       edges still cut off when built with a buffer
#   03_plot_routes.R        static maps of the example routes
#   turn_restriction_schematics.qmd  toy-network schematics (Quarto page;
#                           schematic_helpers.R has its plotting functions)

suppressMessages({
  library(sf)
  library(data.table)
  library(igraph)
})
# also sources osm_network.R and restriction_members.R
source("processing/line_graph.R")

# osmium-tool command. Through pixi it needs no installation:
#   pixi exec -s osmium-tool osmium ...
# If osmium is on the PATH, use OSMIUM <- "osmium"
OSMIUM <- c("pixi", "exec", "-s", "osmium-tool", "osmium")

# Study area (Bogotá test polygon) and outputs
POLYGON <- "data/bogota_poly.geojson"
OUT_DIR <- "outputs"
NETWORK_RDS <- file.path(OUT_DIR, "bogota_network.rds")
EXAMPLES_GPKG <- file.path(OUT_DIR, "routability_examples.gpkg")
PLOT_DIR <- file.path(OUT_DIR, "route_plots")
dir.create(OUT_DIR, showWarnings = FALSE)

# Geofabrik extract for Colombia (downloaded once into osmextract's folder)
colombia_pbf <- function() {
  osmextract::oe_get("Colombia", provider = "geofabrik", download_only = TRUE,
                     skip_vectortranslate = TRUE, quiet = TRUE)
}

study_polygon <- function() st_read(POLYGON, quiet = TRUE)
