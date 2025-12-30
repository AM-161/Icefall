# scripts/02_build_list_page.R
# ============================================================
# Build list page (summary table) for GitHub Pages
# - reads meta: data/Koordinaten_Wasserfaelle/oetztalclimbingice_clean.csv
# - reads assignments: data/AWS/icefalls_nearest_station.csv (optional extra fields)
# - reads topo urls: data/Koordinaten_Wasserfaelle/icefalls_sun_horizon.csv (optional)
# - reads model runs: data/ModelRuns/model_uid<uid>.csv
# - writes: site/icefalls_table.json + site/list.html
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(lubridate)
  library(jsonlite)
})

TZ_LOCAL <- "Europe/Vienna"

# ----------------------------
# Paths
# ----------------------------
PATH_ASSIGN <- "data/AWS/icefalls_nearest_station.csv"
PATH_META   <- "data/Koordinaten_Wasserfaelle/oetztalclimbingice_clean.csv"
PATH_SUN    <- "data/Koordinaten_Wasserfaelle/icefalls_sun_horizon.csv"
DIR_MODELS  <- "data/ModelRuns"

OUT_DIR  <- "site"
OUT_JSON <- file.path(OUT_DIR, "icefalls_table.json")
OUT_HTML <- file.path(OUT_DIR, "list.html")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# Helpers
# ----------------------------
to_num <- function(x) {
  if (is.null(x)) return(NA_real_)
  if (is.numeric(x)) return(x)
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  x <- gsub(",", ".", x, fixed = TRUE)
  suppressWarnings(as.numeric(x))
}

read_any_delim <- function(path) {
  # tries tab, semicolon, comma
  x <- tryCatch(readr::read_delim(path, delim = "\t", show_col_types = FALSE, progress = FALSE), error = function(e) NULL)
  if (!is.null(x) && ncol(x) > 1) return(x)
  x <- tryCatch(readr::read_delim(path, delim = ";", show_col_types = FALSE, progress = FALSE), error = function(e) NULL)
  if (!is.null(x) && ncol(x) > 1) return(x)
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

col_or_na <- function(df, col) {
  if (col %in% names(df)) df[[col]] else rep(NA, nrow(df))
}

coalesce_existing <- function(df, candidates, default = NA_character_) {
  # returns first non-NA among existing columns
  ok <- candidates[candidates %in% names(df)]
  if (length(ok) == 0) return(rep(default, nrow(df)))
  out <- df[[ok[1]]]
  if (length(ok) >= 2) {
    for (k in ok[-1]) out <- dplyr::coalesce(out, df[[k]])
  }
  out
}

parse_time_any <- function(x, tz = TZ_LOCAL) {
  if (inherits(x, "POSIXct")) return(with_tz(x, tz))
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  t <- suppressWarnings(lubridate::ymd_hms(x, tz = tz))
  if (all(is.na(t))) t <- suppressWarnings(lubridate::ymd_hm(x, tz = tz))
  if (all(is.na(t))) t <- suppressWarnings(lubridate::parse_date_time(
    x, orders = c("Ymd HMS", "Ymd HM", "Y-m-d H:M:S", "Y-m-d H:M"),
    tz = tz
  ))
  t
}

fmt_num <- function(x, digits = 2) {
  ifelse(is.finite(x), formatC(x, format = "f", digits = digits), NA_character_)
}

fmt_pct <- function(x, digits = 0) {
  ifelse(is.finite(x), paste0(round(x * 100, digits), "%"), NA_character_)
}

# ----------------------------
# 1) Load meta (CSV)
# ----------------------------
if (!file.exists(PATH_META)) {
  stop(
    "Fehlt: ", PATH_META,
    "\nVorhanden in data/Koordinaten_Wasserfaelle:\n",
    paste(list.files("data/Koordinaten_Wasserfaelle", recursive = TRUE), collapse = "\n")
  )
}

meta <- read_any_delim(PATH_META) %>%
  rename_with(~tolower(.x)) %>%
  mutate(
    uid = as.integer(.data$uid)
  )

# expected columns from your header (lowercased):
# uid, name, latitude, longitude, hoehe_dgm5m, erstbegehnung, schwierigkeit,
# eisfallhhe, zustieg, abstieg, beschreibung, ausrichtung

meta <- meta %>%
  transmute(
    uid = uid,
    name = as.character(coalesce_existing(cur_data(), c("name"), default = NA_character_)),
    latitude  = to_num(coalesce_existing(cur_data(), c("latitude","lat"), default = NA_character_)),
    longitude = to_num(coalesce_existing(cur_data(), c("longitude","lon"), default = NA_character_)),
    elev_m = to_num(coalesce_existing(cur_data(), c("hoehe_dgm5m","hoehe","höhe","elevation","elev_m"), default = NA_character_)),
    difficulty = as.character(coalesce_existing(cur_data(), c("schwierigkeit","difficulty","grad"), default = NA_character_)),
    icefall_height_m = to_num(coalesce_existing(cur_data(), c("eisfallhhe","eisfallhoehe","eisfallhöhe","height_m"), default = NA_character_)),
    aspect = as.character(coalesce_existing(cur_data(), c("ausrichtung","aspect"), default = NA_character_)),
    approach = as.character(coalesce_existing(cur_data(), c("zustieg","approach"), default = NA_character_)),
    descent  = as.character(coalesce_existing(cur_data(), c("abstieg","descent"), default = NA_character_)),
    first_ascent = as.character(coalesce_existing(cur_data(), c("erstbegehnung","first_ascent"), default = NA_character_)),
    description  = as.character(coalesce_existing(cur_data(), c("beschreibung","description"), default = NA_character_))
  ) %>%
  filter(!is.na(uid))

# ----------------------------
# 2) Load assign (optional) to add station/meta fields
# ----------------------------
assign <- NULL
if (file.exists(PATH_ASSIGN)) {
  assign <- read_any_delim(PATH_ASSIGN) %>%
    mutate(uid = as.integer(uid))
}

# ----------------------------
# 3) Load sun/topo urls (optional)
# ----------------------------
sun <- NULL
if (file.exists(PATH_SUN)) {
  sun <- read_any_delim(PATH_SUN) %>%
    mutate(uid = as.integer(uid)) %>%
    select(any_of(c("uid","topo_url","topo_slug","name")))
  # ensure at most 1 row per uid (take first non-empty topo_url)
  sun <- sun %>%
    group_by(uid) %>%
    summarise(
      topo_url  = dplyr::coalesce(first(topo_url[topo_url != ""]), first(topo_url)),
      topo_slug = dplyr::coalesce(first(topo_slug[topo_slug != ""]), first(topo_slug)),
      .groups = "drop"
    )
}

# ----------------------------
# 4) Read model run per uid and compute "tomorrow" metrics
# ----------------------------
tomorrow <- as.Date(with_tz(Sys.time(), TZ_LOCAL) + days(1))

summarise_uid_model <- function(uid) {
  f <- file.path(DIR_MODELS, sprintf("model_uid%s.csv", uid))
  if (!file.exists(f) || file.info(f)$size <= 0) {
    return(tibble(
      uid = uid,
      thickness_tomorrow_07_m = NA_real_,
      climb_max_tomorrow = NA_real_,
      climb_max_time_local = NA_character_,
      thickness_at_climb_max_m = NA_real_
    ))
  }

  df <- readr::read_csv(f, show_col_types = FALSE, progress = FALSE)

  if (!"time" %in% names(df)) {
    return(tibble(
      uid = uid,
      thickness_tomorrow_07_m = NA_real_,
      climb_max_tomorrow = NA_real_,
      climb_max_time_local = NA_character_,
      thickness_at_climb_max_m = NA_real_
    ))
  }

  df <- df %>%
    mutate(
      time = parse_time_any(.data$time, tz = TZ_LOCAL),
      date = as.Date(time),
      thickness_m = if ("thickness_m" %in% names(df)) to_num(thickness_m) else NA_real_,
      climbability = if ("climbability" %in% names(df)) to_num(climbability) else NA_real_
    ) %>%
    filter(!is.na(time))

  df_day <- df %>% filter(date == tomorrow)
  if (nrow(df_day) == 0) {
    return(tibble(
      uid = uid,
      thickness_tomorrow_07_m = NA_real_,
      climb_max_tomorrow = NA_real_,
      climb_max_time_local = NA_character_,
      thickness_at_climb_max_m = NA_real_
    ))
  }

  # thickness at ~07:00 local (closest)
  t07 <- as.POSIXct(paste0(format(tomorrow, "%Y-%m-%d"), " 07:00:00"), tz = TZ_LOCAL)
  i07 <- which.min(abs(as.numeric(difftime(df_day$time, t07, units = "mins"))))
  thickness_07 <- df_day$thickness_m[i07]

  # max climbability tomorrow
  if (all(!is.finite(df_day$climbability))) {
    climb_max <- NA_real_
    climb_time <- NA_character_
    thick_at_best <- NA_real_
  } else {
    imax <- which.max(df_day$climbability)
    climb_max <- df_day$climbability[imax]
    climb_time <- format(df_day$time[imax], "%d.%m.%Y %H:%M")
    thick_at_best <- df_day$thickness_m[imax]
  }

  tibble(
    uid = uid,
    thickness_tomorrow_07_m = thickness_07,
    climb_max_tomorrow = climb_max,
    climb_max_time_local = climb_time,
    thickness_at_climb_max_m = thick_at_best
  )
}

uids <- sort(unique(meta$uid))
model_sum <- bind_rows(lapply(uids, summarise_uid_model))

# ----------------------------
# 5) Merge everything
# ----------------------------
out <- meta %>%
  left_join(model_sum, by = "uid")

if (!is.null(assign)) {
  # bring some useful columns if present
  assign_slim <- assign %>%
    select(any_of(c(
      "uid","station_id","source","dist_km","elev_diff_m",
      "icefall_name","ice_lon","ice_lat","icefall_elev_m","icefall_height_m"
    ))) %>%
    mutate(uid = as.integer(uid))

  out <- out %>%
    left_join(assign_slim, by = "uid") %>%
    mutate(
      name = dplyr::coalesce(as.character(icefall_name), name),
      latitude  = dplyr::coalesce(to_num(ice_lat), latitude),
      longitude = dplyr::coalesce(to_num(ice_lon), longitude),
      elev_m = dplyr::coalesce(to_num(icefall_elev_m), elev_m),
      icefall_height_m = dplyr::coalesce(to_num(icefall_height_m), icefall_height_m)
    )
}

if (!is.null(sun)) {
  out <- out %>% left_join(sun, by = "uid")
} else {
  out$topo_url <- NA_character_
  out$topo_slug <- NA_character_
}

out <- out %>%
  mutate(
    # plot url (must exist under site/plots/)
    plot_url = sprintf("plots/uid_%03d.png", uid),

    thickness_tomorrow_07_txt = fmt_num(thickness_tomorrow_07_m, 2),
    climb_max_tomorrow_txt    = fmt_pct(climb_max_tomorrow, 0),
    thickness_at_best_txt     = fmt_num(thickness_at_climb_max_m, 2)
  ) %>%
  arrange(dplyr::desc(climb_max_tomorrow), dplyr::desc(thickness_tomorrow_07_m))

# ----------------------------
# 6) Write JSON (for JS table)
# ----------------------------
jsonlite::write_json(out, OUT_JSON, pretty = TRUE, auto_unbox = TRUE, na = "null")
message("✅ Wrote JSON: ", OUT_JSON)

# ----------------------------
# 7) Write list.html (sortable + search + fullscreen modal)
# ----------------------------
# Note: this page expects:
# - site/icefalls_table.json
# - site/plots/uid_###.png (generated by scripts/00_build_plots_all.R)
html <- sprintf(
'<!doctype html>
<html lang="de">
<head>
  <meta charset="utf-8"/>
  <meta name="viewport" content="width=device-width, initial-scale=1"/>
  <title>Icefalls – Übersicht</title>
  <style>
    body { font-family: system-ui, -apple-system, Segoe UI, Roboto, sans-serif; margin: 0; padding: 0; }
    header { padding: 10px 14px; border-bottom: 1px solid #ddd; display:flex; gap:12px; align-items:center; flex-wrap:wrap; }
    header a { text-decoration:none; padding:6px 10px; border:1px solid #ddd; border-radius:8px; color:#111; }
    header a:hover { background:#f4f4f4; }
    .wrap { padding: 12px 14px; }
    .controls { display:flex; gap:10px; flex-wrap:wrap; align-items:center; margin-bottom:10px; }
    input[type="search"] { padding:10px 12px; border:1px solid #ccc; border-radius:10px; min-width: 260px; font-size:16px; }
    table { width:100%%; border-collapse: collapse; }
    th, td { padding: 10px 8px; border-bottom: 1px solid #eee; vertical-align: top; }
    th { text-align:left; position: sticky; top: 0; background: #fff; z-index: 1; cursor:pointer; user-select:none; }
    tr:hover { background: #fafafa; }
    .muted { color:#666; font-size:12px; }
    .pill { display:inline-block; padding:2px 8px; border:1px solid #ddd; border-radius:999px; font-size:12px; }
    .btn { display:inline-flex; gap:6px; align-items:center; padding:6px 10px; border:1px solid #ddd; border-radius:10px; background:#fff; cursor:pointer; }
    .btn:hover { background:#f4f4f4; }
    .small { font-size: 12px; }

    /* Fullscreen modal */
    #modal { display:none; position:fixed; inset:0; background:rgba(0,0,0,0.8); z-index:9999; }
    #modal .inner { position:absolute; inset:0; display:flex; flex-direction:column; }
    #modal .bar { padding:10px; display:flex; gap:10px; align-items:center; justify-content:space-between; color:#fff; }
    #modal img { flex:1; width:100%%; height:100%%; object-fit: contain; }
    #modal .bar button, #modal .bar a {
      color:#fff; border:1px solid rgba(255,255,255,0.35);
      background: transparent; padding:8px 12px; border-radius:10px; cursor:pointer;
      text-decoration:none;
    }
    #modal .bar button:hover, #modal .bar a:hover { background: rgba(255,255,255,0.12); }

    @media (max-width: 720px) {
      th, td { padding: 12px 6px; }
      input[type="search"] { width: 100%%; min-width: 0; }
      header { gap:8px; }
    }
  </style>
</head>
<body>
  <header>
    <a href="index.html">🗺️ Karte</a>
    <a href="list.html"><b>📋 Übersicht</b></a>
    <span class="muted">Morgen: %s (TZ: Europe/Vienna)</span>
  </header>

  <div class="wrap">
    <div class="controls">
      <input id="q" type="search" placeholder="Suchen: Name, Schwierigkeit, Ausrichtung, Station …">
      <span class="muted">Klick auf Spaltenkopf = sortieren</span>
    </div>

    <div class="muted small" id="status">Lade Daten …</div>

    <table id="tbl">
      <thead>
        <tr>
          <th data-key="name">Eisfall</th>
          <th data-key="difficulty">Schwierigkeit</th>
          <th data-key="elev_m">Höhe (m)</th>
          <th data-key="icefall_height_m">Eisfallhöhe (m)</th>
          <th data-key="thickness_tomorrow_07_m">Eisdicke morgen ~07:00 (m)</th>
          <th data-key="climb_max_tomorrow">Max. Kletterbarkeit morgen</th>
          <th data-key="climb_max_time_local">Uhrzeit</th>
          <th>Diagramm</th>
          <th>Topo</th>
        </tr>
      </thead>
      <tbody></tbody>
    </table>
  </div>

  <!-- Fullscreen modal -->
  <div id="modal">
    <div class="inner">
      <div class="bar">
        <div id="modalTitle">Diagramm</div>
        <div style="display:flex; gap:10px; align-items:center;">
          <a id="openNewTab" href="#" target="_blank" rel="noopener">In neuem Tab</a>
          <button id="closeModal">Schließen</button>
        </div>
      </div>
      <img id="modalImg" src="" alt="Diagramm"/>
    </div>
  </div>

<script>
(function(){
  const status = document.getElementById("status");
  const q = document.getElementById("q");
  const tbody = document.querySelector("#tbl tbody");
  const ths = Array.from(document.querySelectorAll("th[data-key]"));

  const modal = document.getElementById("modal");
  const modalImg = document.getElementById("modalImg");
  const modalTitle = document.getElementById("modalTitle");
  const closeModal = document.getElementById("closeModal");
  const openNewTab = document.getElementById("openNewTab");

  let rows = [];
  let sortKey = "climb_max_tomorrow";
  let sortAsc = false;

  function num(x){
    if (x === null || x === undefined) return NaN;
    const n = Number(x);
    return isFinite(n) ? n : NaN;
  }

  function str(x){
    if (x === null || x === undefined) return "";
    return String(x);
  }

  function matches(r, query){
    if(!query) return true;
    const t = query.toLowerCase();
    const blob = [
      r.name, r.difficulty, r.aspect, r.station_id, r.source,
      r.approach, r.descent
    ].map(str).join(" | ").toLowerCase();
    return blob.includes(t);
  }

  function cmp(a,b){
    const va = a[sortKey];
    const vb = b[sortKey];

    // try numeric first
    const na = num(va), nb = num(vb);
    if (isFinite(na) && isFinite(nb)) return sortAsc ? (na-nb) : (nb-na);
    if (isFinite(na) && !isFinite(nb)) return sortAsc ? -1 : 1;
    if (!isFinite(na) && isFinite(nb)) return sortAsc ? 1 : -1;

    // fallback string
    const sa = str(va).toLowerCase();
    const sb = str(vb).toLowerCase();
    if (sa < sb) return sortAsc ? -1 : 1;
    if (sa > sb) return sortAsc ? 1 : -1;
    return 0;
  }

  function openFullscreen(plotUrl, title){
    modalImg.src = plotUrl;
    modalTitle.textContent = title || "Diagramm";
    openNewTab.href = plotUrl;
    modal.style.display = "block";
  }

  function closeFullscreen(){
    modal.style.display = "none";
    modalImg.src = "";
  }

  closeModal.addEventListener("click", closeFullscreen);
  modal.addEventListener("click", function(e){
    if (e.target === modal) closeFullscreen();
  });
  document.addEventListener("keydown", function(e){
    if (e.key === "Escape") closeFullscreen();
  });

  function render(){
    const query = q.value.trim();
    const view = rows.filter(r => matches(r, query)).sort(cmp);

    tbody.innerHTML = "";
    for(const r of view){
      const tr = document.createElement("tr");

      const topoLink = r.topo_url ? `<a href="${r.topo_url}" target="_blank" rel="noopener">Topo</a>` : "<span class=\\"muted\\">—</span>";
      const plotUrl = r.plot_url || "";
      const plotBtn = plotUrl
        ? `<button class="btn" data-plot="${plotUrl}" data-title="${str(r.name).replace(/"/g, "&quot;")}">🔍 Vollbild</button>`
        : `<span class="muted">—</span>`;

      tr.innerHTML = `
        <td>
          <div><b>${str(r.name) || ("UID " + r.uid)}</b></div>
          <div class="muted">
            ${r.station_id ? ("Station: " + str(r.station_id) + (r.source ? (" (" + str(r.source) + ")") : "")) : ""}
          </div>
        </td>
        <td>${str(r.difficulty) || "<span class=\\"muted\\">—</span>"}</td>
        <td>${isFinite(num(r.elev_m)) ? Math.round(num(r.elev_m)) : "<span class=\\"muted\\">—</span>"}</td>
        <td>${isFinite(num(r.icefall_height_m)) ? Math.round(num(r.icefall_height_m)) : "<span class=\\"muted\\">—</span>"}</td>
        <td>${r.thickness_tomorrow_07_txt || "<span class=\\"muted\\">—</span>"}</td>
        <td>${r.climb_max_tomorrow_txt || "<span class=\\"muted\\">—</span>"}</td>
        <td>${str(r.climb_max_time_local) || "<span class=\\"muted\\">—</span>"}</td>
        <td>${plotBtn} ${plotUrl ? `<a class="muted small" href="${plotUrl}" target="_blank" rel="noopener">neu Tab</a>` : ""}</td>
        <td>${topoLink}</td>
      `;
      tbody.appendChild(tr);
    }

    // bind plot buttons
    Array.from(document.querySelectorAll("button[data-plot]")).forEach(btn => {
      btn.addEventListener("click", () => {
        openFullscreen(btn.getAttribute("data-plot"), btn.getAttribute("data-title"));
      });
    });

    status.textContent = `Einträge: ${view.length} / ${rows.length}  | Sort: ${sortKey} ${sortAsc ? "↑" : "↓"}`;
  }

  ths.forEach(th => {
    th.addEventListener("click", () => {
      const key = th.getAttribute("data-key");
      if (key === sortKey) sortAsc = !sortAsc;
      else { sortKey = key; sortAsc = true; }
      render();
    });
  });

  q.addEventListener("input", render);

  fetch("icefalls_table.json", {cache: "no-store"})
    .then(r => r.json())
    .then(data => {
      rows = data || [];
      status.textContent = "Daten geladen.";
      render();
    })
    .catch(err => {
      status.textContent = "Fehler beim Laden: " + err;
    });
})();
</script>
</body>
</html>
', format(tomorrow, "%d.%m.%Y")
)

writeLines(html, OUT_HTML, useBytes = TRUE)
message("✅ Wrote HTML: ", OUT_HTML)

# Final info
message("Done. Outputs:")
message(" - ", normalizePath(OUT_JSON, winslash = "/", mustWork = FALSE))
message(" - ", normalizePath(OUT_HTML, winslash = "/", mustWork = FALSE))
