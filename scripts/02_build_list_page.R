suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(lubridate)
  library(jsonlite)
})

TZ_LOCAL <- "Europe/Vienna"

# Inputs
PATH_ASSIGN   <- "data/AWS/icefalls_nearest_station.csv"
PATH_MODELS   <- "data/ModelRuns"
META_DIR      <- "data/Koordinaten_Wasserfaelle/oetztalclimbingice_clean"

# Outputs
OUT_JSON <- "site/icefalls_table.json"
OUT_HTML <- "site/list.html"

dir.create("site", recursive = TRUE, showWarnings = FALSE)

# -------------------- helpers --------------------
read_any_delim <- function(path) {
  x <- tryCatch(readr::read_delim(path, delim = "\t", show_col_types = FALSE, progress = FALSE), error = function(e) NULL)
  if (!is.null(x) && ncol(x) > 1) return(x)
  x <- tryCatch(readr::read_delim(path, delim = ";", show_col_types = FALSE, progress = FALSE), error = function(e) NULL)
  if (!is.null(x) && ncol(x) > 1) return(x)
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

pick_meta_file <- function(dir_path) {
  if (!dir.exists(dir_path)) stop("META_DIR existiert nicht: ", dir_path)
  files <- list.files(dir_path, full.names = TRUE)
  # bevorzugt .csv/.tsv/.txt
  cand <- files[grepl("\\.(csv|tsv|txt)$", files, ignore.case = TRUE)]
  if (length(cand) == 0) cand <- files
  if (length(cand) == 0) stop("Keine Dateien in META_DIR: ", dir_path)
  cand[1]
}

parse_time_local <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  if (all(is.na(x))) return(as.POSIXct(character(), tz = TZ_LOCAL))

  # wenn String Zeitzone/Offset/Z enthält -> als UTC parsen und nach Lokal konvertieren
  has_tz <- grepl("Z$|[+-]\\d\\d:?\\d\\d$", x)
  if (any(has_tz, na.rm = TRUE)) {
    t <- suppressWarnings(lubridate::parse_date_time(
      x, orders = c("Ymd HMS", "Ymd HM", "Y-m-d H:M:S", "Y-m-d H:M", "Y-m-d\\TH:M:S", "Y-m-d\\TH:M"),
      tz = "UTC"
    ))
    return(with_tz(as.POSIXct(t), TZ_LOCAL))
  }

  # sonst: als Lokalzeit interpretieren (deine model_uid*.csv sind typischerweise so)
  t <- suppressWarnings(lubridate::parse_date_time(
    x, orders = c("Ymd HMS", "Ymd HM", "Y-m-d H:M:S", "Y-m-d H:M", "Y-m-d\\TH:M:S", "Y-m-d\\TH:M"),
    tz = TZ_LOCAL
  ))
  as.POSIXct(t, tz = TZ_LOCAL)
}

pick_nearest_row <- function(df, target_time) {
  if (nrow(df) == 0) return(NULL)
  dt <- abs(as.numeric(difftime(df$time_local, target_time, units = "secs")))
  df[which.min(dt), , drop = FALSE]
}

uid_pad3 <- function(u) sprintf("%03d", as.integer(u))

# -------------------- read inputs --------------------
stopifnot(file.exists(PATH_ASSIGN))
assign <- readr::read_csv(PATH_ASSIGN, show_col_types = FALSE, progress = FALSE)

# wichtig: deine assign hat mind. uid + icefall_name + station stuff
assign <- assign %>%
  mutate(uid = as.integer(uid))

# Metadaten-Datei (deine Tabelle mit uid/name/lat/lon/hoehe_dgm5m/schwierigkeit/eisfallhhe/ausrichtung…)
meta_file <- pick_meta_file(META_DIR)
message("Meta file: ", meta_file)

meta <- read_any_delim(meta_file) %>%
  mutate(uid = as.integer(uid))

# normalize Spaltennamen minimal (falls Tippfehler/verschiedene Schreibweisen)
# wir gehen von deinen Spalten aus:
# uid name latitude longitude hoehe_dgm5m ... schwierigkeit eisfallhhe ... ausrichtung
if (!"icefall_name" %in% names(assign) && "name" %in% names(assign)) {
  assign$icefall_name <- assign$name
}

# joint: Stationzuordnung + Metadaten (Metadaten gewinnen)
meta_joined <- assign %>%
  select(any_of(c("uid","icefall_name","station_id","source","dist_km","elev_diff_m","ice_lat","ice_lon","topo_url"))) %>%
  left_join(meta, by = "uid") %>%
  mutate(
    # finaler Name: erst meta$name, sonst assign$icefall_name, sonst UID fallback
    name_clean = dplyr::coalesce(
      suppressWarnings(as.character(.data$name)),
      suppressWarnings(as.character(.data$icefall_name)),
      paste0("UID ", uid)
    )
  )

# uids: lieber aus Metadaten (vollständig), fallback assign
uids <- sort(unique(dplyr::coalesce(meta$uid, meta_joined$uid)))
uids <- uids[is.finite(uids)]

# -------------------- tomorrow window (Vienna) --------------------
now_local  <- with_tz(Sys.time(), TZ_LOCAL)
today      <- as.Date(now_local)
tomorrow   <- today + 1
t_target   <- as.POSIXct(paste(tomorrow, "07:00:00"), tz = TZ_LOCAL)

# -------------------- build rows --------------------
rows <- list()

for (u in uids) {
  model_path <- file.path(PATH_MODELS, sprintf("model_uid%s.csv", u))
  if (!file.exists(model_path)) next

  mod <- readr::read_csv(model_path, show_col_types = FALSE, progress = FALSE)
  if (!("time" %in% names(mod))) next

  mod <- mod %>%
    mutate(
      time_local  = parse_time_local(time),
      date_local  = as.Date(time_local),
      thickness_m = suppressWarnings(as.numeric(thickness_m)),
      climbability = suppressWarnings(as.numeric(climbability))
    ) %>%
    filter(!is.na(time_local))

  day_df <- mod %>% filter(date_local == tomorrow)

  # Eisdicke morgen ~07:00
  thick_row <- pick_nearest_row(day_df, t_target)
  thick_m <- if (!is.null(thick_row)) thick_row$thickness_m[1] else NA_real_

  # Max climbability morgen + Uhrzeit
  best_row <- day_df %>%
    filter(is.finite(climbability)) %>%
    slice_max(order_by = climbability, n = 1, with_ties = FALSE)

  ci_max  <- if (nrow(best_row) == 1) best_row$climbability[1] else NA_real_
  ci_time <- if (nrow(best_row) == 1) format(best_row$time_local[1], "%d.%m.%Y %H:%M") else NA_character_

  meta_u <- meta_joined %>% filter(uid == u) %>% slice(1)
  if (nrow(meta_u) == 0) {
    name_clean <- paste0("UID ", u)
    meta_u <- tibble(uid = u, name_clean = name_clean)
  }

  # Plot-Link (nur wenn vorhanden, sonst leer)
  plot_rel <- paste0("plots/uid_", uid_pad3(u), ".png")
  plot_exists <- file.exists(file.path("site", plot_rel))
  plot_link <- if (plot_exists) plot_rel else ""

  rows[[length(rows) + 1]] <- tibble(
    uid = u,
    name = meta_u$name_clean[1],

    eisdicke_morgen_m_0700 = thick_m,
    ci_max_morgen = ci_max,
    ci_time_morgen = ci_time,

    schwierigkeit = if ("schwierigkeit" %in% names(meta_u)) as.character(meta_u$schwierigkeit[1]) else NA_character_,
    laenge_m = NA_real_,  # falls du später Länge als Spalte hast -> hier mappen
    eisfallhoehe_m = if ("eisfallhhe" %in% names(meta_u)) suppressWarnings(as.numeric(meta_u$eisfallhhe[1])) else NA_real_,
    hoehe_ueNN_m = if ("hoehe_dgm5m" %in% names(meta_u)) suppressWarnings(as.numeric(meta_u$hoehe_dgm5m[1])) else NA_real_,
    ausrichtung = if ("ausrichtung" %in% names(meta_u)) as.character(meta_u$ausrichtung[1]) else NA_character_,

    station_id = if ("station_id" %in% names(meta_u)) as.character(meta_u$station_id[1]) else NA_character_,
    source     = if ("source" %in% names(meta_u)) as.character(meta_u$source[1]) else NA_character_,
    dist_km    = if ("dist_km" %in% names(meta_u)) suppressWarnings(as.numeric(meta_u$dist_km[1])) else NA_real_,
    dz_m       = if ("elev_diff_m" %in% names(meta_u)) suppressWarnings(as.numeric(meta_u$elev_diff_m[1])) else NA_real_,

    plot = plot_link,
    topo_url = if ("topo_url" %in% names(meta_u)) as.character(meta_u$topo_url[1]) else NA_character_
  )
}

df <- bind_rows(rows) %>%
  mutate(
    eisdicke_morgen_m_0700 = round(eisdicke_morgen_m_0700, 2),
    ci_max_morgen = round(ci_max_morgen, 2),
    eisfallhoehe_m = round(eisfallhoehe_m, 0),
    hoehe_ueNN_m = round(hoehe_ueNN_m, 0),
    dist_km = round(dist_km, 2),
    dz_m = round(dz_m, 0)
  )

write_json(df, OUT_JSON, pretty = TRUE, auto_unbox = TRUE)
message("✅ geschrieben: ", OUT_JSON, " (rows=", nrow(df), ")")

# -------------------- list.html (DataTables) --------------------
html <- paste0(
'<!doctype html>
<html lang="de">
<head>
  <meta charset="utf-8"/>
  <meta name="viewport" content="width=device-width,initial-scale=1"/>
  <title>Eisfälle – Übersicht</title>

  <link rel="stylesheet" href="https://cdn.datatables.net/1.13.8/css/jquery.dataTables.min.css"/>
  <style>
    body{font-family:system-ui, -apple-system, Segoe UI, Roboto, Helvetica, Arial; margin:18px;}
    h1{font-size:20px;margin:0 0 10px 0;}
    .hint{font-size:12px;color:#444;margin-bottom:12px;}
    table.dataTable thead th{white-space:nowrap;}
    a{color:#0b57d0;text-decoration:none;}
    a:hover{text-decoration:underline;}
    .toplinks{margin-bottom:12px;font-size:13px;}
  </style>
</head>
<body>
  <h1>Eisfälle – Übersicht (morgen)</h1>
  <div class="toplinks">
    <a href="index.html">← zurück zur Karte</a>
  </div>
  <div class="hint">
    Eisdicke ist Wert nahe <b>07:00</b> (Europe/Vienna). CI max = Maximum morgen + Uhrzeit. Sortieren: Spaltenkopf klicken.
  </div>

  <table id="tbl" class="display" style="width:100%">
    <thead>
      <tr>
        <th>UID</th>
        <th>Name</th>
        <th>Eisdicke morgen (m, 07:00)</th>
        <th>CI max morgen</th>
        <th>CI Zeit</th>
        <th>Schwierigkeit</th>
        <th>Eisfallhöhe (m)</th>
        <th>üNN (m)</th>
        <th>Ausrichtung</th>
        <th>Station</th>
        <th>Quelle</th>
        <th>Dist (km)</th>
        <th>dz (m)</th>
        <th>Plot</th>
        <th>Topo</th>
      </tr>
    </thead>
    <tbody></tbody>
  </table>

  <script src="https://code.jquery.com/jquery-3.7.1.min.js"></script>
  <script src="https://cdn.datatables.net/1.13.8/js/jquery.dataTables.min.js"></script>

  <script>
    async function main(){
      const res = await fetch("icefalls_table.json", {cache:"no-store"});
      const data = await res.json();

      $("#tbl").DataTable({
        data: data,
        pageLength: 25,
        order: [[3, "desc"]],  // default: CI max absteigend
        columns: [
          {data:"uid"},
          {data:"name"},
          {data:"eisdicke_morgen_m_0700"},
          {data:"ci_max_morgen"},
          {data:"ci_time_morgen"},
          {data:"schwierigkeit"},
          {data:"eisfallhoehe_m"},
          {data:"hoehe_ueNN_m"},
          {data:"ausrichtung"},
          {data:"station_id"},
          {data:"source"},
          {data:"dist_km"},
          {data:"dz_m"},
          {data:"plot", render:(d)=> d ? `<a href="${d}" target="_blank">Plot</a>` : "" },
          {data:"topo_url", render:(d)=> d ? `<a href="${d}" target="_blank">Topo</a>` : "" }
        ]
      });
    }
    main();
  </script>
</body>
</html>'
)

writeLines(html, OUT_HTML, useBytes = TRUE)
message("✅ geschrieben: ", OUT_HTML)
