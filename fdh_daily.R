# =====================================================================
# FDH-basiertes Eiswachstumsmodell für Nordtirol
# - INCA 1h-Daten laden (ab 2025-11-01) in 2-Tages-Chunks
# - FDH & MDH berechnen (Freezing & Melting Degree Hours)
# - Wind- und Feuchteeinfluss (Konvektion / latente Wärme) berücksichtigen
# - zeitaufgelöste Sonnenhöhe (Datum + Tageszeit) + Exposition (DEM-basiert)
# - zeitlich gewichtete Temperaturgeschichte (jüngere Stunden stärker)
# - effektive FDH -> Eisdicke
# - interaktive Leaflet-Karte als eisdicke_nordtirol.html speichern
# =====================================================================

# Datenquellen:
#   INCA (GeoSphere Austria), DOI: https://doi.org/10.60669/6akt-5p05
#   DEM Tirol: https://www.data.gv.at/katalog/datasets/0454f5f3-1d8c-464e-847d-541901eb021a
# Autor: @antifascist_mountaineer (https://www.instagram.com/antifascist_mountaineer/)

library(httr)
library(raster)
library(leaflet)
library(htmlwidgets)
library(ncdf4)


# 0) Zeitraum & Gebiet -------------------------------------------------

start_all    <- as.Date("2025-11-01")        # Startdatum der Saison
end_all      <- Sys.Date()                    # bis heute
chunk_days   <- 2                             # 2-Tages-Chunks (API-Limit)
chunk_starts <- seq.Date(start_all, end_all, by = chunk_days)

# Nordtirol-BBox in WGS84
bbox <- c(
  46.7,  # lat_min
  10.1,  # lon_min
  47.7,  # lat_max
  12.2   # lon_max
)

# INCA-Parameter (erweitert für Wind, Feuchte etc.)
parameters <- c("RR", "T2M", "RH2M", "UU", "VV", "GL", "P0", "TD2M")
param_str  <- paste(parameters, collapse = ",")

base_url <- "https://dataset.api.hub.geosphere.at/v1/grid/historical/inca-v1-1h-1km"

out_dir <- "data/inca_nordtirol"  # etwas neutraler Ordnername
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

nc_files <- character(0)


# 1) INCA-Daten chunkweise laden (existierende Dateien überspringen) ---

for (cs in chunk_starts) {
  cs <- as.Date(cs)
  ce <- as.Date(min(cs + (chunk_days - 1), end_all))
  
  message("Chunk: ", as.character(cs), " bis ", as.character(ce))
  
  start_time <- paste0(format(cs, "%Y-%m-%d"), "T00:00")
  end_time   <- paste0(format(ce, "%Y-%m-%d"), "T23:00")
  
  outfile <- file.path(
    out_dir,
    sprintf("inca_nordtirol_%s_%s.nc",
            format(cs, "%Y%m%d"),
            format(ce, "%Y%m%d"))
  )
  
  # schon vorhanden -> Download überspringen
  if (file.exists(outfile) && file.info(outfile)$size > 0) {
    message("  -> übersprungen (existiert): ", outfile)
  } else {
    query <- list(
      parameters    = param_str,
      start         = start_time,
      end           = end_time,
      bbox          = paste(bbox, collapse = ","),
      output_format = "netcdf",
      filename      = "inca_nordtirol"
    )
    
    resp <- GET(
      url   = base_url,
      query = query,
      write_disk(outfile, overwrite = TRUE)
    )
    stop_for_status(resp)
    
    sz <- file.info(outfile)$size
    message("  -> gespeichert: ", outfile, " (", sz, " Bytes)")
  }
  
  nc_files <- c(nc_files, outfile)
}

nc_files <- sort(unique(nc_files))


# 2) Alle NetCDFs einlesen und in ein Objekt packen --------------------

# Zeit-Units in POSIXct umrechnen
convert_nc_time <- function(time_vals, time_units) {
  if (!grepl("since", time_units)) {
    origin <- as.POSIXct("1970-01-01 00:00:00", tz = "UTC")
    return(origin + time_vals * 3600)
  }
  unit_str   <- sub(" since.*", "", time_units)
  origin_str <- sub(".*since ", "", time_units)
  origin     <- as.POSIXct(origin_str, tz = "UTC")
  if (grepl("hour", unit_str, ignore.case = TRUE)) {
    origin + time_vals * 3600
  } else if (grepl("second", unit_str, ignore.case = TRUE)) {
    origin + time_vals
  } else if (grepl("day", unit_str, ignore.case = TRUE)) {
    origin + time_vals * 86400
  } else {
    origin + time_vals
  }
}

vars_wanted <- c("T2M", "RH2M", "RR", "GL", "UU", "VV", "TD2M", "P0")

# 2.1 erste Datei öffnen -> Dimensionen, lon/lat
nc0       <- nc_open(nc_files[1])
varnames0 <- names(nc0$var)

lon_name <- intersect(c("lon", "longitude", "x", "xc"), varnames0)[1]
lat_name <- intersect(c("lat", "latitude", "y", "yc"), varnames0)[1]

lon0 <- ncvar_get(nc0, lon_name)  # [nx, ny]
lat0 <- ncvar_get(nc0, lat_name)  # [nx, ny]

nx <- dim(lon0)[1]
ny <- dim(lon0)[2]

time_dim_name <- names(nc0$dim)[grep("time", tolower(names(nc0$dim)))[1]]

vars_found <- intersect(vars_wanted, varnames0)
nc_close(nc0)

# 2.2 Zeitlängen pro Datei
nt_per_file <- integer(length(nc_files))
for (i in seq_along(nc_files)) {
  nc <- nc_open(nc_files[i])
  time_vals <- ncvar_get(nc, time_dim_name)
  nt_per_file[i] <- length(time_vals)
  nc_close(nc)
}

nt_total <- sum(nt_per_file)

# 2.3 Arrays allozieren

time_all  <- as.POSIXct(rep(NA_real_, nt_total), origin = "1970-01-01", tz = "UTC")
inca_data <- lapply(vars_found, function(.) array(NA_real_, dim = c(nx, ny, nt_total)))
names(inca_data) <- vars_found

# 2.4 Daten einlesen

pos <- 1
for (i in seq_along(nc_files)) {
  nc <- nc_open(nc_files[i])
  
  time_vals  <- ncvar_get(nc, time_dim_name)
  time_units <- ncatt_get(nc, time_dim_name, "units")$value
  nt_i       <- length(time_vals)
  
  idx <- pos:(pos + nt_i - 1)
  time_all[idx] <- convert_nc_time(time_vals, time_units)
  
  for (v in vars_found) {
    arr_i <- ncvar_get(nc, v)    # [nx, ny, nt_i]
    inca_data[[v]][ , , idx] <- arr_i
  }
  
  nc_close(nc)
  pos <- pos + nt_i
}

inca_nordtirol_all <- list(
  lon  = lon0,
  lat  = lat0,
  time = time_all,
  data = inca_data
)


# 3) FDH/MDH (stundenweise) + Orientierung + Raster-Template ----------

# Basis-Arrays (C, %, m/s)
T2M_arr <- inca_nordtirol_all$data$T2M   # °C [nx, ny, nt]
RH_arr  <- inca_nordtirol_all$data$RH2M  # %  [nx, ny, nt]
UU_arr  <- inca_nordtirol_all$data$UU    # m/s
VV_arr  <- inca_nordtirol_all$data$VV    # m/s

lon <- inca_nordtirol_all$lon            # [nx, ny]
lat <- inca_nordtirol_all$lat            # [nx, ny]

nx <- dim(T2M_arr)[1]
ny <- dim(T2M_arr)[2]
nt <- dim(T2M_arr)[3]

# 3.1 Orientierung (S->N, W->E) für alle Felder

lat_mean_j <- colMeans(lat, na.rm = TRUE)
if (lat_mean_j[1] < lat_mean_j[length(lat_mean_j)]) {
  idx_j <- ny:1
  lon   <- lon[, idx_j]
  lat   <- lat[, idx_j]
  
  T2M_arr <- T2M_arr[, idx_j, , drop = FALSE]
  RH_arr  <- RH_arr[,  idx_j, , drop = FALSE]
  UU_arr  <- UU_arr[,  idx_j, , drop = FALSE]
  VV_arr  <- VV_arr[,  idx_j, , drop = FALSE]
}

lon_mean_i <- rowMeans(lon, na.rm = TRUE)
if (lon_mean_i[1] > lon_mean_i[length(lon_mean_i)]) {
  idx_i <- nx:1
  lon   <- lon[idx_i, ]
  lat   <- lat[idx_i, ]
  
  T2M_arr <- T2M_arr[idx_i, , , drop = FALSE]
  RH_arr  <- RH_arr[idx_i, , , drop = FALSE]
  UU_arr  <- UU_arr[idx_i, , , drop = FALSE]
  VV_arr  <- VV_arr[idx_i, , , drop = FALSE]
}

nx <- dim(T2M_arr)[1]
ny <- dim(T2M_arr)[2]
nt <- dim(T2M_arr)[3]

# 3.2 FDH & MDH (C·h) + Windgeschwindigkeit

FDH_hourly <- pmax(-T2M_arr, 0)   # T < 0
MDH_hourly <- pmax( T2M_arr, 0)   # T > 0

W_arr <- sqrt(UU_arr^2 + VV_arr^2)   # m/s, [nx, ny, nt]

# Referenz-Wind (für Normierung)
wind_ref <- mean(W_arr, na.rm = TRUE)
if (!is.finite(wind_ref) || wind_ref <= 0) wind_ref <- 5  # fallback

# 3.3 "Plain FDH" (ohne Korrekturen) + Raster-Template

FDH_sum_plain <- apply(FDH_hourly, c(1, 2), sum, na.rm = TRUE)  # [nx, ny]
FDH_mat_plain <- t(FDH_sum_plain)                               # [ny, nx]

r_template <- raster(FDH_mat_plain)
extent(r_template) <- c(min(lon, na.rm = TRUE),
                        max(lon, na.rm = TRUE),
                        min(lat, na.rm = TRUE),
                        max(lat, na.rm = TRUE))
crs(r_template) <- "EPSG:4326"


# 4) DEM auf INCA-Grid + Expositionsindex (solar_index_ij) -----------

# vorverarbeitetes DEM im INCA-Grid (EPSG:4326, ~1 km), vorher lokal aus 10-m-DEM abgeleitet
# Datei liegt im Repo unter data/DEM_Tirol_INCAgrid_1km_epsg4326.tif
dem_inca <- raster("data/DEM_Tirol_INCAgrid_1km_epsg4326.tif")
crs(dem_inca) <- "EPSG:4326"
names(dem_inca) <- "elev_m"

# Slope & Aspect im INCA-Grid
sl_as       <- terrain(dem_inca, opt = c("slope", "aspect"), unit = "degrees")
aspect_inca <- sl_as[["aspect"]]   # 0°=Ost, 90°=Nord, etc.

# Aspect so drehen, dass 0° = Nord
aspect_rad_from_north <- (aspect_inca - 90) * pi / 180

# +1 = Nordhang, -1 = Südhang
northness_r <- cos(aspect_rad_from_north)

# 0 = Nord (schattig), 1 = Süd (besonnt)
solar_index_r <- (1 - northness_r) / 2

# Raster -> Matrix [nx, ny] passend zu FDH_hourly
solar_mat      <- as.matrix(solar_index_r)  # [nrow = ny, ncol = nx]
solar_index_ij <- t(solar_mat)              # [nx, ny]


# 5) Zeitabhängige Sonnenhöhe + Zeit-Gewichtung ----------------------

# Zeitvektor (UTC) und lokale Zeit in Tirol
time_vec   <- inca_nordtirol_all$time          # Länge nt_total, aber
# wir gehen davon aus: Reihenfolge identisch zu FDH_hourly über die Saison
# -> hier kontrollieren wir nichts extra und nutzen alle Zeiten durch

# Beschränken: zur Sicherheit auf nt Einträge (falls minimaler Off-by-One)
if (length(time_vec) != nt) {
  # einfache Anpassung: kürzen / nehmen die letzten nt Werte
  if (length(time_vec) > nt) {
    time_vec <- tail(time_vec, nt)
  } else {
    stop("Länge von time_vec (", length(time_vec), ") != nt (", nt, ")")
  }
}

time_local <- as.POSIXlt(time_vec, tz = "Europe/Vienna")

doy  <- time_local$yday + 1
hour <- time_local$hour + time_local$min / 60 + time_local$sec / 3600

# Sonnen-Deklination (rad) pro Zeitpunkt
delta_t <- 23.44 * pi/180 * sin(2 * pi * (284 + doy) / 365)

# mittlere Breite des Gebiets
lat_center_deg <- (min(lat, na.rm = TRUE) + max(lat, na.rm = TRUE)) / 2
lat_center_rad <- lat_center_deg * pi / 180

# Stundenwinkel (rad): 0 bei 12:00, ± am Morgen/Abend
H_t <- (hour - 12) * 15 * pi / 180

# sinus der Sonnenhöhe (rad)
sin_alpha_t <- sin(lat_center_rad) * sin(delta_t) +
               cos(lat_center_rad) * cos(delta_t) * cos(H_t)

# Nacht -> 0
sin_alpha_t[sin_alpha_t < 0] <- 0

# Faktor 0–1 (0 = Nacht, 1 = max. Sonnenstand)
solar_height_factor_t <- sin_alpha_t

# Zeit-Gewichtung (neuere Stunden wichtiger)
max_time  <- max(time_vec, na.rm = TRUE)
age_days  <- as.numeric(difftime(max_time, time_vec, units = "days"))

tau_days  <- 10   # e-Falldauer ~10 Tage (ältere Beiträge werden stark abgewertet)
weight_time_t <- exp(-age_days / tau_days)  # 1 für jetzt, -> 0 für weit zurück


# 6) Effektive FDH (FDHm + Wind + Feuchte + Strahlung + Zeit) ---------

# Parameter für physikalische Korrekturen
k_expo  <- 0.5   # Stärke des Strahlungs-/Expositions-Effekts (0–1)

k_wind  <- 0.5   # Einfluss der Windabweichung auf den konvektiven Term
wind_min <- 0.5  # min. Faktor
wind_max <- 2.0  # max. Faktor

rh_opt  <- 0.7   # optimale rel. Feuchte (als Anteil, ~70 %)
rh_sig  <- 0.15  # Breite des Feuchte-Peaks

k_melt  <- 1.2   # Gewichtung der MDH (Schmelz-Effekt)

FDH_sum_eff <- matrix(0, nrow = nx, ncol = ny)

for (k in seq_len(nt)) {
  # Basisfelder für diesen Zeitpunkt
  FDH_k <- FDH_hourly[ , , k]
  MDH_k <- MDH_hourly[ , , k]
  W_k   <- W_arr[ , , k]
  RH_k  <- RH_arr[ , , k] / 100  # -> 0..1
  
  # Wind-Faktor (konvektive Kühlung)
  wind_norm <- (W_k - wind_ref) / wind_ref
  f_wind    <- 1 + k_wind * wind_norm
  f_wind[!is.finite(f_wind)] <- 1
  f_wind <- pmax(wind_min, pmin(wind_max, f_wind))
  
  # Feuchte-Faktor für Wachstum (Peak bei rh_opt, 0.5..1)
  f_rh_peak <- exp(- (RH_k - rh_opt)^2 / (2 * rh_sig^2))
  f_rh      <- 0.5 + 0.5 * f_rh_peak
  
  # Strahlungs-/Expositions-Faktor (reduziert Wachstum bei starker Sonne auf Südhängen)
  f_rad <- 1 - k_expo * solar_height_factor_t[k] * solar_index_ij
  f_rad[!is.finite(f_rad)] <- 1
  f_rad <- pmax(0, pmin(1, f_rad))
  
  # Zeit-Gewichtung
  f_time <- weight_time_t[k]
  if (!is.finite(f_time)) f_time <- 1
  
  # Effektive Freezing Degree Hours (mit Wind, Feuchte, Strahlung, Zeit)
  FDH_eff_k <- FDH_k * f_wind * f_rh * f_rad * f_time
  
  # Effektive Melting Degree Hours (Schmelz-Term, zeitgewichtet)
  MDH_eff_k <- MDH_k * k_melt * f_time
  
  # Netto-Beitrag (kann theoretisch auch negativ sein)
  FDH_sum_eff <- FDH_sum_eff + (FDH_eff_k - MDH_eff_k)
}

# Physikalisch: kein "negatives Eis" -> untere Grenze 0
FDH_sum_eff[FDH_sum_eff < 0] <- 0

# Effektive FDH -> Raster
FDH_mat_eff <- t(FDH_sum_eff)   # [ny, nx]

r_FDH_eff <- raster(FDH_mat_eff)
extent(r_FDH_eff) <- extent(r_template)
crs(r_FDH_eff)    <- crs(r_template)
names(r_FDH_eff)  <- "FDH_eff_C_h"



# 7) Effektive FDH -> Eisdicke + Climbability-Index ------------------

# Energie -> Eisdicke (einfacher lineare Ansatz)
h_c      <- 30       # W m^-2 K^-1
rho_i    <- 880      # kg m^-3
Lf       <- 334000   # J kg^-1
FDH_crit <- 50       # °C·h (Minimum für "relevante" Vereisung)

alpha <- h_c * 3600 / (rho_i * Lf)   # m / (°C·h)

ice_thick_expo <- r_FDH_eff * alpha
ice_thick_expo[r_FDH_eff < FDH_crit] <- NA
names(ice_thick_expo) <- "h_pot_expo_m"


# Matrix der Eisdicke (für Climbability)
h_mat_ij <- FDH_sum_eff * alpha   # [nx, ny]

# Aktuelle Bedingungen (letzter Zeitschritt)
T_last_ij <- T2M_arr[ , , nt]
RH_last_ij <- RH_arr[ , , nt] / 100

# 3-Tages-Mittel (Tr3) der Lufttemperatur
n_tr3 <- min(72, nt)   # 72 h = 3 Tage
idx_tr3 <- (nt - n_tr3 + 1):nt

Tr3_ij <- apply(T2M_arr[ , , idx_tr3, drop = FALSE], c(1, 2), mean, na.rm = TRUE)

# Scoring-Funktionen für Climbability-Index (alles sehr heuristisch)

# 7.1 Dicke-Score: 0 unter h_min, 1 ab h_opt
h_min <- 0.10   # 10 cm: untere Grenze für "überhaupt Eis da"
h_opt <- 0.50   # 50 cm: "gut kletterbar"

score_h <- (h_mat_ij - h_min) / (h_opt - h_min)
score_h[h_mat_ij <= h_min] <- 0
score_h[h_mat_ij >= h_opt] <- 1
score_h[!is.finite(score_h)] <- 0

# 7.2 Aktuelle Temperatur (T_last): optimum ~ -4 °C
T_opt  <- -4
T_min  <- -20
T_max  <- 0

range_T <- max(T_opt - T_min, T_max - T_opt)
score_T <- 1 - abs(T_last_ij - T_opt) / range_T
score_T[T_last_ij <= T_min | T_last_ij >= T_max] <- 0
score_T[score_T < 0] <- 0
score_T[!is.finite(score_T)] <- 0

# 7.3 Tr3 (3-Tages-Mittel): optimal etwas kälter (~ -6 °C)
T3_opt <- -6
T3_min <- -20
T3_max <- -1

range_T3 <- max(T3_opt - T3_min, T3_max - T3_opt)
score_T3 <- 1 - abs(Tr3_ij - T3_opt) / range_T3
score_T3[Tr3_ij <= T3_min | Tr3_ij >= T3_max] <- 0
score_T3[score_T3 < 0] <- 0
score_T3[!is.finite(score_T3)] <- 0

# 7.4 Luftfeuchte (mittlere Werte bevorzugt)
RH_opt_c <- 0.7
RH_sig_c <- 0.2

score_RH_peak <- exp(- (RH_last_ij - RH_opt_c)^2 / (2 * RH_sig_c^2))
score_RH      <- score_RH_peak
score_RH[!is.finite(score_RH)] <- 0

# 7.5 Gesamt-Climbability-Index (0..1)

climb_index_ij <- score_h * score_T * score_T3 * score_RH
climb_index_ij[h_mat_ij <= h_min] <- 0
climb_index_ij[!is.finite(climb_index_ij)] <- 0

# Raster bauen
climb_mat_yx <- t(climb_index_ij)   # [ny, nx]

climb_r <- raster(climb_mat_yx)
extent(climb_r) <- extent(r_template)
crs(climb_r)    <- crs(r_template)
names(climb_r)  <- "Climbability_0_1"


# (Optional) als GeoTIFF speichern
# writeRaster(ice_thick_expo, "Eisdicke_pot_expo_Tirol_1km.tif", format = "GTiff", overwrite = TRUE)
# writeRaster(climb_r,      "Climbability_Tirol_1km.tif",       format = "GTiff", overwrite = TRUE)


# 8) Interaktive Leaflet-Karte (Eisdicke & Climbability) -----------------

# Display-Version der Raster für eine weichere Darstellung (nur fürs Plotten)
# Original: ~1 km; hier z.B. auf ~250 m verfeinern (Faktor 4)
fact_disp <- 4

ice_thick_disp <- disaggregate(
  ice_thick_expo,
  fact   = fact_disp,
  method = "bilinear"  # weiche Interpolation nur für die Anzeige
)

climb_disp <- disaggregate(
  climb_r,
  fact   = fact_disp,
  method = "bilinear"
)

# Climbability = 0 unsichtbar machen, damit die Basiskarte gut sichtbar bleibt
climb_disp_masked <- climb_disp
climb_disp_masked[climb_disp_masked <= 0] <- NA

max_h   <- cellStats(ice_thick_expo, max, na.rm = TRUE)
col_fun <- colorRampPalette(c("white", "orange", "green", "darkgreen"))

pal_h <- colorNumeric(
  palette  = col_fun(100),
  domain   = c(0, max_h),
  na.color = "transparent"
)

# Palette für Climbability (0–1)
pal_ci <- colorNumeric(
  palette  = rev(terrain.colors(100)),
  domain   = c(0, 1),
  na.color = "transparent"
)

d_today <- Sys.Date()

m <- leaflet() |>
  addTiles(group = "OSM") |>
  addRasterImage(
    ice_thick_disp,
    colors  = pal_h,
    opacity = 0.8,
    project = TRUE,
    group   = "Eisdicke"
  ) |>
  addRasterImage(
    climb_disp_masked,
    colors  = pal_ci,
    opacity = 0.7,
    project = TRUE,
    group   = "Climbability"
  ) |>
  addLegend(
    pal       = pal_h,
    values    = c(0, max_h),
    title     = sprintf(
      "Eisdicke (m)
mit Exposition & Sonne (%s)",
      format(d_today, "%Y-%m-%d")
    ),
    labFormat = labelFormat(digits = 2),
    position  = "bottomright"
  ) |>
  addLegend(
    pal       = pal_ci,
    values    = c(0, 1),
    title     = "Climbability (0–1)",
    labFormat = labelFormat(digits = 1),
    position  = "bottomleft"
  ) |>
  addLayersControl(
    baseGroups    = c("OSM"),
    overlayGroups = c("Eisdicke", "Climbability"),
    options       = layersControlOptions(collapsed = FALSE)
  ) |>
  addControl(
    html = htmltools::HTML(
      "<div style='font-size: 10px; background: rgba(255,255,255,0.9); padding: 4px 6px; border-radius: 4px; max-width: 260px; line-height: 1.3; text-align: center; position: relative; left: 50%; transform: translateX(-50%);'><strong>Quellen:</strong> INCA (GeoSphere Austria, <a href='https://doi.org/10.60669/6akt-5p05' target='_blank'>doi:10.60669/6akt-5p05</a>); DEM Tirol (<a href='https://www.data.gv.at/katalog/datasets/0454f5f3-1d8c-464e-847d-541901eb021a' target='_blank'>data.gv.at</a>)<br/><strong>Autor:</strong> <a href='https://www.instagram.com/antifascist_mountaineer/' target='_blank'>@antifascist_mountaineer</a></div>"
    ),
    position = "bottomleft"
  )

# Popup beim Laden: kurze Erklärung zu Karte, Eisdicke & Climbability
m <- htmlwidgets::onRender(
  m,
  "
  function(el, x) {
    var map = this;
    var content =
      '<b>Eisfall-Karte Nordtirol – experimentelles Modell</b><br/>' +
      '<b>Eisdicke (m):</b> Abschätzung aus Freezing Degree Hours (Summe der Stunden mit Lufttemperatur unter 0&nbsp;°C), ' +
      'skaliert über einen einfachen Energiebilanz-Ansatz (konvektiver Wärmetransport) und korrigiert mit Hangexposition, Hangneigung ' +
      'und jahreszeitlicher Sonnenhöhe. Auflösung ~1&nbsp;km – beschreibt eher Gitterzellen/Regionen als konkrete einzelne Eislinien.<br/>' +
      '<b>Climbability (0–1):</b> heuristischer Index, der Eisdicke mit aktuellen Schmelzstunden (Temperaturen über 0&nbsp;°C), ' +
      'Strahlungsbelastung (Südexposition, flache Wintersonne) und groben Wind-/Feuchtebedingungen kombiniert. Werte nahe 1 = eher kalte, ' +
      'potenziell günstigere Phasen; nahe 0 = eher ungünstige, schmelzdominierte Phasen.<br/>' +
      '<b>Wichtiger Hinweis:</b> Kein Lawinen- oder Eisgutachten, nur grobe regionale Orientierung. Lokale Verhältnisse (Wasserführung, ' +
      'Lawinenkegel, Felsqualität, Eisstruktur etc.) sowie der aktuelle Lawinenlagebericht und deine eigene Erfahrung sind immer entscheidend.';
    L.popup({maxWidth: 340})({maxWidth: 320})
      .setLatLng(map.getCenter())
      .setContent(content)
      .openOn(map);
  }
  "
)

# für GitHub Pages / Web speichern
saveWidget(m, "eisdicke_nordtirol.html", selfcontained = TRUE)
