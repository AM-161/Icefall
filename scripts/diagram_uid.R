# scripts/00_build_plots_all.R

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

assign_path <- "data/AWS/icefalls_nearest_station.csv"
if (!file.exists(assign_path)) stop("Fehlt: ", assign_path)

assign <- read_csv(assign_path, show_col_types = FALSE)

if (!("uid" %in% names(assign))) stop("Spalte 'uid' fehlt in ", assign_path)

uids <- sort(unique(assign$uid))
uids <- uids[is.finite(uids)]

if (length(uids) == 0) stop("Keine UIDs gefunden.")

dir.create("site/plots", recursive = TRUE, showWarnings = FALSE)

# Diagramm-Skript Pfad
plot_script <- "scripts/diagram_uid.R"
if (!file.exists(plot_script)) stop("Fehlt Diagramm-Skript: ", plot_script)

message("Starte Plot-Build für ", length(uids), " UIDs in EINEM R-Prozess...")

args <- c("--vanilla", plot_script, as.character(uids))

out <- system2(
  "Rscript",
  args,
  stdout = TRUE,
  stderr = TRUE
)

status <- attr(out, "status")

if (!is.null(status) && status != 0) {
  cat(paste(out, collapse = "\n"), "\n")
  stop("❌ Plot-Build fehlgeschlagen (exit ", status, ").")
}

message("\n✅ Alle Plots gebaut: ", length(uids))
