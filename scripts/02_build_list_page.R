suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(lubridate)
  library(jsonlite)
})

TZ_LOCAL <- "Europe/Vienna"

PATH_ASSIGN <- "data/AWS/icefalls_nearest_station.csv"
PATH_MODELS <- "data/ModelRuns"
OUT_JSON    <- "site/icefalls_table.json"
OUT_HTML    <- "site/list.html"

dir.create("site", recursive = TRUE, showWarnings = FALSE)

# --- Helpers ----------------------------------------------------------
parse_time_local <- function(x) {
  # robust: versucht ISO / ymd_hms / ymd_hm; behandelt UTC/Offsets und bringt nach Europe/Vienna
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  t <- suppressWarnings(ymd_hms(x, tz = "UTC"))
  if (all(is.na(t))) t <- suppressWarnings(ymd_hm(x, tz = "UTC"))
  if (all(is.na(t))) t <- suppressWarnings(parse_date_time(
    x,
    orders = c("Ymd HMS", "Ymd HM", "Y-m-d H:M:S", "Y-m-d H:M", "Y-m-d\\TH:M:S", "Y-m-d\\TH:M"),
    tz = "UTC"
  ))
  with_tz(t, TZ_LOCAL)
}

pick_nearest <- function(df, target_time) {
  if (nrow(df) == 0) return(NULL)
  dt <- abs(as.numeric(difftime(df$time_local, target_time, units = "secs")))
  df[which.min(dt), , drop = FALSE]
}

# --- Inputs -----------------------------------------------------------
stopifnot(file.exists(PATH_ASSIGN))
assign <- read_csv(PATH_ASSIGN, show_col_types = FALSE)

# erwartete Metadaten-Spalten (wenn vorhanden, werden sie übernommen)
meta_cols <- c(
  "uid","name","icefall_name","ice_name",
  "icefall_elev_m","icefall_height_m","icefall_length_m","icefall_grade",
  "difficulty","length_m","height_m","elev_m",
  "topo_url","longitude","latitude",
  "station_id","source","dist_km","elev_diff_m"
)
assign_meta <- assign %>%
  select(any_of(meta_cols)) %>%
  mutate(uid = as.integer(uid)) %>%
  distinct(uid, .keep_all = TRUE)

# Name sauber setzen
assign_meta <- assign_meta %>%
  mutate(
    name_clean = dplyr::coalesce(icefall_name, paste0("UID ", uid))
  )

# --- Tomorrow window --------------------------------------------------
today_local   <- as.Date(with_tz(Sys.time(), TZ_LOCAL))
tomorrow      <- today_local + 1
t_target      <- as.POSIXct(paste(tomorrow, "07:00:00"), tz = TZ_LOCAL)  # <- Fixzeitpunkt
tomorrow_end  <- as.POSIXct(paste(tomorrow, "23:59:59"), tz = TZ_LOCAL)

# --- Build summary per UID -------------------------------------------
uids <- sort(unique(assign_meta$uid))
rows <- list()

for (u in uids) {
  f <- file.path(PATH_MODELS, sprintf("model_uid%s.csv", u))
  if (!file.exists(f)) next

  mod <- read_csv(f, show_col_types = FALSE)
  if (!("time" %in% names(mod))) next

  mod <- mod %>%
    mutate(
      time_local = parse_time_local(time),
      date_local = as.Date(time_local),
      thickness_m = as.numeric(thickness_m),
      climbability = as.numeric(climbability)
    ) %>%
    filter(!is.na(time_local))

  day_df <- mod %>% filter(date_local == tomorrow)

  # Eisdicke morgen: Wert nahe 07:00
  thick_row <- pick_nearest(day_df, t_target)
  thick_m <- if (!is.null(thick_row) && is.finite(thick_row$thickness_m[1])) thick_row$thickness_m[1] else NA_real_

  # Max Climbability morgen + Zeit
  best_row <- day_df %>%
    filter(is.finite(climbability)) %>%
    slice_max(order_by = climbability, n = 1, with_ties = FALSE)

  ci_max  <- if (nrow(best_row) == 1) best_row$climbability[1] else NA_real_
  ci_time <- if (nrow(best_row) == 1) format(best_row$time_local[1], "%d.%m.%Y %H:%M") else NA_character_

  meta <- assign_meta %>% filter(uid == u) %>% slice(1)

  uid_pad <- sprintf("%03d", u)
  plot_url <- paste0("plots/uid_", uid_pad, ".png")

  rows[[length(rows)+1]] <- tibble(
    uid = u,
    name = meta$name_clean,
    eisdicke_morgen_m_0700 = thick_m,
    ci_max_morgen = ci_max,
    ci_time_morgen = ci_time,
    # Metadaten (nur wenn vorhanden)
    grade = coalesce(meta$icefall_grade, meta$difficulty),
    length_m = coalesce(meta$icefall_length_m, meta$length_m),
    height_m = coalesce(meta$icefall_height_m, meta$height_m),
    elev_m = coalesce(meta$icefall_elev_m, meta$elev_m),
    station_id = meta$station_id,
    source = meta$source,
    dist_km = meta$dist_km,
    dz_m = coalesce(meta$elev_diff_m),
    topo_url = meta$topo_url,
    plot = plot_url
  )
}

df <- bind_rows(rows) %>%
  mutate(
    eisdicke_morgen_m_0700 = round(eisdicke_morgen_m_0700, 2),
    ci_max_morgen = round(ci_max_morgen, 2),
    length_m = suppressWarnings(as.numeric(length_m)),
    height_m = suppressWarnings(as.numeric(height_m)),
    elev_m = suppressWarnings(as.numeric(elev_m))
  )

# --- Write JSON for the table ----------------------------------------
write_json(df, OUT_JSON, pretty = TRUE, auto_unbox = TRUE)
message("✅ geschrieben: ", OUT_JSON, " (rows=", nrow(df), ")")

# --- list.html (DataTables) ------------------------------------------
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
    .small{font-size:12px;color:#333;}
    a{color:#0b57d0;text-decoration:none;}
    a:hover{text-decoration:underline;}
  </style>
</head>
<body>
  <h1>Eisfälle – Übersicht (morgen)</h1>
  <div class="hint">
    Eisdicke ist der Wert nahe <b>07:00</b> (Europe/Vienna). Climbability: <b>Maximum</b> am morgigen Tag + Uhrzeit.
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
        <th>Länge (m)</th>
        <th>Höhe (m)</th>
        <th>üNN (m)</th>
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
          {data:"grade"},
          {data:"length_m"},
          {data:"height_m"},
          {data:"elev_m"},
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
