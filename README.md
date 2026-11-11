# explore-pbf-restrictions

Build a routable, directed road network with OSM turn restrictions from a
`.pbf` file (no Overpass API limits), and check that it routes correctly. The
test case is a polygon in Bogotá, cut from the Geofabrik extract for Colombia.

The network is turned into a **line graph** (each directed road edge is a
vertex, each allowed turn is an edge), in the same way as OSRM:

- via-node restrictions: the banned turn is left out (`only_*` keeps just the
  one allowed turn);
- via-way restrictions: the via edges are copied so that only the banned
  sequence is removed;
- U-turns: allowed at dead ends, turning circles and unsignalised junctions
  (`uturn_rules()`);
- turn costs: OSRM car-profile penalties for turn angle, traffic signals and
  U-turns (`osrm_turn_params()`).

## Structure

```
processing/   functions and the network build
  restriction_members.R   restriction relations from a .pbf
  osm_network.R           directed driving network from a .pbf
  line_graph.R            line graph, restrictions, U-turns, turn penalties
  config.R                settings: osmium command, paths, input .pbf
  01_build_network.R      -> outputs/bogota_network.rds
checks/       checks and results
  01_routability_checks.R -> outputs/routability_examples.gpkg
  02_buffer_check.R       -> layer cut_off_with_buffer
  03_plot_routes.R        -> outputs/route_plots/*.png
  04_network_size.R       primal vs line graph size -> outputs/network_size.*
  turn_restriction_schematics.qmd  how OTP, OSRM, GraphHopper and this
                          project apply restrictions, on toy networks
data/         study polygon (bogota_poly.geojson)
outputs/      generated files
```

## Requirements

- R with `sf`, `data.table`, `igraph`, `sfheaders`, `osmextract`; `tmap` and
  `maptiles` for the route maps; `sfnetworks`, `ggplot2` and `patchwork` for
  the schematics.
- [osmium-tool](https://osmcode.org/osmium-tool/). By default it is run
  through [pixi](https://pixi.sh) (`pixi exec -s osmium-tool osmium`), so it
  needs no installation; set `OSMIUM` in `processing/config.R` to use another
  one.
- [Quarto](https://quarto.org) to render the schematics page.

## Run

From the repository root:

```bash
Rscript processing/01_build_network.R
Rscript checks/01_routability_checks.R
Rscript checks/02_buffer_check.R
Rscript checks/03_plot_routes.R
Rscript checks/04_network_size.R
quarto render checks/turn_restriction_schematics.qmd
```

The first run downloads the Colombia extract with `osmextract`.

## License

AGPL-3.0, see [LICENSE](LICENSE).
