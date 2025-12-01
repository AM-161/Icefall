# =====================================================================
# Urheberrechtlicher Hinweis / Nutzungseinschränkung
# ---------------------------------------------------------------------
# © 2025 Carlos Wydra. Alle Rechte vorbehalten.
#
# Dieses Skript wurde im Rahmen der Masterarbeit zur Kletterbarkeit
# von Eisfällen in Nordtirol erstellt. Der Code darf ohne vorherige
# ausdrückliche schriftliche Zustimmung des Autors weder ganz noch
# teilweise kopiert, verändert, weiterverbreitet oder in eigene
# Projekte (insbesondere kommerzielle Anwendungen oder proprietäre
# Software) integriert werden.
#
# Eine Nutzung über das reine Lesen und Nachvollziehen im Kontext
# dieser Masterarbeit hinaus ist nicht gestattet. Bei Rückfragen
# oder einer gewünschten Nutzung bitte vorher Kontakt mit dem Autor
# aufnehmen.
# =====================================================================

# =====================================================================
# FDH-basiertes Eiswachstumsmodell für Nordtirol
# - INCA 1h-Daten laden (ab 2025-11-01) in 2-Tages-Chunks
# - FDH & MDH berechnen (Freezing & Melting Degree Hours)
# - Wind- und Feuchteeinfluss (Konvektion / latente Wärme) berücksichtigen
# - zeitaufgelöste Sonnenhöhe (Datum + Tageszeit) + Exposition (DEM-basiert)
# - zeitlich gewichtete Temperaturgeschichte (jüngere Stunden stärker)
# - effektive FDH -> Eisdicke
# - Climbability-Index (0–1) für aktuelle Bedingungen
# - NWP-Vorhersage (AROME, nwp-v1-1h-2500m) für ~60 h
# - Prognose-Eisdicke & -Climbability für mehrere Zeitschritte (0, +12, +24, +36, +48, +60 h)
# - interaktive Leaflet-Karte mit Slider für Prognose (eisdicke_nordtirol.html)
# =====================================================================

# Datenquellen:
#   INCA (GeoSphere Austria), DOI: https://doi.org/10.60669/6akt-5p05
#   DEM Tirol: https://www.data.gv.at/katalog/datasets/0454f5f3-1d8c-464e-847d-541901eb021a
#   NWP AROME (nwp-v1-1h-2500m): siehe GeoSphere API
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

base_url_inca <- "https://dataset.api.hub.geosphere.at/v1/grid/historical/inca-v1-1h-1km"

out_dir <- "data/inca_nordtirol"
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
      url   = base_url_inca,
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

# vorverarbeitetes DEM im INCA-Grid (EPSG:4326, ~1 km),
# vorher lokal aus 10-m-DEM abgeleitet
# Datei liegt im Repo im Hauptverzeichnis unter DEM_Tirol_INCAgrid_1km_epsg4326.tif

dem_inca <- raster("DEM_Tirol_INCAgrid_1km_epsg4326.tif")
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

# Zeitvektor (UTC)
time_vec <- inca_nordtirol_all$time

# Beschränken: zur Sicherheit auf nt Einträge
if (length(time_vec) != nt) {
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

# --- Zeitachsen-Infos für historischen Eisdicken-Zeitverlauf ---
t_start <- min(time_vec, na.rm = TRUE)
time_offset_hours <- as.numeric(difftime(time_vec, t_start, units = "hours"))
snap_hours_hist <- seq(
  from = 0,
  to   = floor(max(time_offset_hours, na.rm = TRUE)),
  by   = 24
)
FDH_hist_layers <- list()
FDH_hist_labels <- character(0)
FDH_hist_times  <- as.POSIXct(character(0), tz = "UTC")
snap_idx_hist <- 1


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
  
  # Historische FDH-Schnappschüsse für Zeit-Slider (täglich)
  t_k_hours <- time_offset_hours[k]
  while (snap_idx_hist <= length(snap_hours_hist) &&
         t_k_hours >= snap_hours_hist[snap_idx_hist] - 0.5) {
    
    FDH_k_snap <- FDH_sum_eff
    FDH_k_snap[FDH_k_snap < 0] <- 0
    
    FDH_mat_snap <- t(FDH_k_snap)
    r_FDH_snap <- raster(FDH_mat_snap)
    extent(r_FDH_snap) <- extent(r_template)
    crs(r_FDH_snap)    <- crs(r_template)
    
    label_time <- t_start + snap_hours_hist[snap_idx_hist] * 3600
    label_str  <- format(label_time, "%d.%m.%Y")
    
    FDH_hist_layers[[length(FDH_hist_layers) + 1]] <- r_FDH_snap
    FDH_hist_labels <- c(FDH_hist_labels, label_str)
    FDH_hist_times  <- c(FDH_hist_times, label_time)
    
    snap_idx_hist <- snap_idx_hist + 1
  }
}


# Physikalisch: kein "negatives Eis" -> untere Grenze 0
FDH_sum_eff[FDH_sum_eff < 0] <- 0

# Effektive FDH -> Raster
FDH_mat_eff <- t(FDH_sum_eff)   # [ny, nx]

r_FDH_eff <- raster(FDH_mat_eff)
extent(r_FDH_eff) <- extent(r_template)
crs(r_FDH_eff)    <- crs(r_template)
names(r_FDH_eff)  <- "FDH_eff_C_h"


# 7) Effektive FDH -> Eisdicke + Climbability-Index (aktuell) ---------

# Energie -> Eisdicke (einfacher linearer Ansatz)
h_c      <- 30       # W m^-2 K^-1
rho_i    <- 880      # kg m^-3
Lf       <- 334000   # J kg^-1
FDH_crit <- 30        # °C·h; Schwellwert hier deaktiviert (sonst im Frühwinter überall NA)

alpha <- h_c * 3600 / (rho_i * Lf)   # m / (°C·h)

# Historische Eisdicken-Layer aus FDH_hist_layers (falls vorhanden)
ice_hist_layers <- list()
ice_hist_labels <- character(0)

if (length(FDH_hist_layers) > 0) {
  ice_hist_layers <- vector("list", length(FDH_hist_layers))
  ice_hist_labels <- FDH_hist_labels
  for (i in seq_along(FDH_hist_layers)) {
    r_fd <- FDH_hist_layers[[i]]
    ice_i <- r_fd * alpha
    ice_i[r_fd < FDH_crit] <- NA
    names(ice_i) <- paste0("h_pot_expo_hist_", ice_hist_labels[i])
    ice_hist_layers[[i]] <- ice_i
  }
}


ice_thick_expo <- r_FDH_eff * alpha
ice_thick_expo[r_FDH_eff < FDH_crit] <- NA
names(ice_thick_expo) <- "h_pot_expo_m"

# Matrix der Eisdicke (für aktuellen Climbability)
h_mat_ij <- FDH_sum_eff * alpha   # [nx, ny]

# Aktuelle Bedingungen (letzter Zeitschritt)
T_last_ij  <- T2M_arr[ , , nt]
RH_last_ij <- RH_arr[ , , nt] / 100

# 3-Tages-Mittel (Tr3) der Lufttemperatur
n_tr3   <- min(72, nt)   # 72 h = 3 Tage
idx_tr3 <- (nt - n_tr3 + 1):nt

Tr3_ij <- apply(T2M_arr[ , , idx_tr3, drop = FALSE], c(1, 2), mean, na.rm = TRUE)

# Scoring-Parameter (wie bisher)

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

# Climbability-Historie (tägliche Snapshots entlang der Saison)
climb_hist_layers <- list()
climb_hist_labels <- character(0)

if (length(FDH_hist_layers) > 0 && length(FDH_hist_labels) == length(FDH_hist_layers)) {
  climb_hist_layers <- vector("list", length(FDH_hist_layers))
  climb_hist_labels <- FDH_hist_labels
  for (i in seq_along(FDH_hist_layers)) {
    r_fd <- FDH_hist_layers[[i]]
    # Eisdicke-Snapshot
    h_i <- r_fd * alpha
    h_i[r_fd < FDH_crit] <- NA
    
    # zugehöriger Zeitpunkt
    t_i <- FDH_hist_times[i]
    dt_vec <- abs(as.numeric(difftime(time_vec, t_i, units = "hours")))
    k_i <- which.min(dt_vec)
    if (!is.finite(k_i) || length(k_i) == 0) next
    
    # Temperatur / Feuchte als Raster
    T_i_mat  <- t(T2M_arr[ , , k_i])
    RH_i_mat <- t(RH_arr[ , , k_i] / 100)
    r_T_i    <- raster(T_i_mat);  extent(r_T_i)  <- extent(r_template); crs(r_T_i)  <- crs(r_template)
    r_RH_i   <- raster(RH_i_mat); extent(r_RH_i) <- extent(r_template); crs(r_RH_i) <- crs(r_template)
    
    # 3-Tages-Mittel bis zu diesem Zeitpunkt
    k_start <- max(1, k_i - (72 - 1))
    idx3    <- k_start:k_i
    Tr3_i_arr <- apply(T2M_arr[ , , idx3, drop = FALSE], c(1, 2), mean, na.rm = TRUE)
    Tr3_i_mat <- t(Tr3_i_arr)
    r_Tr3_i   <- raster(Tr3_i_mat); extent(r_Tr3_i) <- extent(r_template); crs(r_Tr3_i) <- crs(r_template)
    
    # Scores für diesen Tag
    score_h_i <- calc(h_i, fun = function(h) {
      s <- (h - h_min) / (h_opt - h_min)
      s[h <= h_min] <- 0
      s[h >= h_opt] <- 1
      s[!is.finite(s)] <- 0
      s
    })
    
    score_T_i <- calc(r_T_i, fun = function(Tv) {
      s <- 1 - abs(Tv - T_opt) / range_T
      s[Tv <= T_min | Tv >= T_max] <- 0
      s[s < 0] <- 0
      s[!is.finite(s)] <- 0
      s
    })
    
    score_T3_i <- calc(r_Tr3_i, fun = function(T3v) {
      s <- 1 - abs(T3v - T3_opt) / range_T3
      s[T3v <= T3_min | T3v >= T3_max] <- 0
      s[s < 0] <- 0
      s[!is.finite(s)] <- 0
      s
    })
    
    score_RH_i <- calc(r_RH_i, fun = function(rh) {
      peak <- exp(- (rh - RH_opt_c)^2 / (2 * RH_sig_c^2))
      s <- peak
      s[!is.finite(s)] <- 0
      s
    })
    
    climb_i <- overlay(score_h_i, score_T_i, score_T3_i, score_RH_i,
                       fun = function(a, b, c, d) a * b * c * d)
    
    climb_i[is.na(h_i) | h_i <= h_min] <- NA
    climb_i[!is.finite(climb_i)]       <- NA
    climb_i[climb_i <= 0]              <- NA
    names(climb_i) <- paste0("Climbability_hist_", climb_hist_labels[i])
    climb_hist_layers[[i]] <- climb_i
  }
}

# 7.5 Gesamt-Climbability-Index (0..1)

climb_index_ij <- score_h * score_T * score_T3 * score_RH
climb_index_ij[h_mat_ij <= h_min]        <- NA
climb_index_ij[!is.finite(climb_index_ij)] <- NA
climb_index_ij[climb_index_ij <= 0]        <- NA

# Raster für aktuellen Climbability
climb_mat_yx <- t(climb_index_ij)   # [ny, nx]

climb_r <- raster(climb_mat_yx)
extent(climb_r) <- extent(r_template)
crs(climb_r)    <- crs(r_template)
names(climb_r)  <- "Climbability_0_1"


# 8) Prognose-Eisdicke & -Climbability aus NWP (nwp-v1-1h-2500m) -----

base_url_nwp   <- "https://dataset.api.hub.geosphere.at/v1/grid/forecast/nwp-v1-1h-2500m"
parameters_nwp <- c("t2m", "rh2m", "u10m", "v10m")
param_str_nwp  <- paste(parameters_nwp, collapse = ",")

# Letzte INCA-Zeit als "jetzt" (UTC)
t_now    <- max(inca_nordtirol_all$time, na.rm = TRUE)
fc_h_max <- 60  # Stunden Vorhersage (max. ~60 h)

t_fc_end <- t_now + fc_h_max * 3600

start_fc <- format(t_now,    "%Y-%m-%dT%H:%M")
end_fc   <- format(t_fc_end, "%Y-%m-%dT%H:%M")

outfile_fc <- file.path(
  out_dir,
  sprintf("nwp_nordtirol_%s_%sh.nc",
          format(as.Date(t_now), "%Y%m%d"),
          fc_h_max)
)

if (!file.exists(outfile_fc) || file.info(outfile_fc)$size == 0) {
  message("NWP-Chunk: ", start_fc, " bis ", end_fc)
  
  query_fc <- list(
    parameters    = param_str_nwp,
    start         = start_fc,
    end           = end_fc,
    bbox          = paste(bbox, collapse = ","),
    output_format = "netcdf",
    filename      = "nwp_nordtirol"
  )
  
  resp_fc <- GET(
    url   = base_url_nwp,
    query = query_fc,
    write_disk(outfile_fc, overwrite = TRUE)
  )
  stop_for_status(resp_fc)
  message("  -> NWP-Vorhersage gespeichert: ", outfile_fc,
          " (", file.info(outfile_fc)$size, " Bytes)")
} else {
  message("  -> NWP-Vorhersage übersprungen (existiert): ", outfile_fc)
}

ice_fc_layers   <- NULL   # Liste von RasterLayern (Eisdicke-Prognose)
climb_fc_layers <- NULL   # Liste von RasterLayern (Climbability-Prognose)
fc_step_hours   <- NULL

if (file.exists(outfile_fc)) {
  
  # 8.1 Zeitvektor der Vorhersage aus NetCDF lesen
  nc_fc <- nc_open(outfile_fc)
  
  time_dim_name_fc <- names(nc_fc$dim)[grep("time", tolower(names(nc_fc$dim)))[1]]
  time_fc_vals     <- ncvar_get(nc_fc, time_dim_name_fc)
  time_fc_units    <- ncatt_get(nc_fc, time_dim_name_fc, "units")$value
  time_fc          <- convert_nc_time(time_fc_vals, time_fc_units)
  
  nc_close(nc_fc)
  
  # Führende Zeitdifferenz zur letzten INCA-Stunde (h)
  lead_hours <- as.numeric(difftime(time_fc, t_now, units = "hours"))
  
  keep_fc <- which(lead_hours >= 0 & lead_hours <= fc_h_max + 0.01)
  if (length(keep_fc) > 0) {
    
    time_fc    <- time_fc[keep_fc]
    lead_hours <- lead_hours[keep_fc]
    
    # 8.2 NWP-Variablen als RasterBrick lesen und auf INCA-Grid bringen
    T2M_fc_brick  <- brick(outfile_fc, varname = "t2m")[[keep_fc]]   # °C
    RH2M_fc_brick <- brick(outfile_fc, varname = "rh2m")[[keep_fc]]  # %
    U10_fc_brick  <- brick(outfile_fc, varname = "u10m")[[keep_fc]]  # m/s (Ost)
    V10_fc_brick  <- brick(outfile_fc, varname = "v10m")[[keep_fc]]  # m/s (Nord)
    
    crs(T2M_fc_brick)  <- crs(r_template)
    crs(RH2M_fc_brick) <- crs(r_template)
    crs(U10_fc_brick)  <- crs(r_template)
    crs(V10_fc_brick)  <- crs(r_template)
    
    # Auf INCA-Grid resamplen (1 km, EPSG:4326)
    T2M_fc  <- resample(T2M_fc_brick,  r_template, method = "bilinear")
    RH2M_fc <- resample(RH2M_fc_brick, r_template, method = "bilinear")
    U10_fc  <- resample(U10_fc_brick,  r_template, method = "bilinear")
    V10_fc  <- resample(V10_fc_brick,  r_template, method = "bilinear")
    
    # 8.3 Hilfsgrößen: Wind, Sonnenhöhe, Zeitgewicht (Zukunft = 1)
    W_fc <- overlay(U10_fc, V10_fc,
                    fun = function(u, v) sqrt(u^2 + v^2))
    
    time_local_fc <- as.POSIXlt(time_fc, tz = "Europe/Vienna")
    doy_fc  <- time_local_fc$yday + 1
    hour_fc <- time_local_fc$hour +
      time_local_fc$min / 60 + time_local_fc$sec / 3600
    
    delta_fc <- 23.44 * pi/180 * sin(2 * pi * (284 + doy_fc) / 365)
    H_fc     <- (hour_fc - 12) * 15 * pi/180
    sin_alpha_fc <- sin(lat_center_rad) * sin(delta_fc) +
      cos(lat_center_rad) * cos(delta_fc) * cos(H_fc)
    sin_alpha_fc[sin_alpha_fc < 0] <- 0
    solar_height_factor_fc <- sin_alpha_fc
    
    weight_time_fc <- rep(1, length(time_fc))   # Zukunft gleich gewichtet
    
    # 8.4 "3-Tage-Regime" für Prognose: Mittelwert über gesamten Vorhersagehorizont
    Tr3_fc <- calc(T2M_fc, fun = function(x) mean(x, na.rm = TRUE))
    
    score_T3_fc <- calc(Tr3_fc, fun = function(t3) {
      s <- 1 - abs(t3 - T3_opt) / range_T3
      s[t3 <= T3_min | t3 >= T3_max] <- 0
      s[s < 0] <- 0
      s[!is.finite(s)] <- 0
      s
    })
    
    # 8.5 Ziel-Zeitpunkte (Slider-Schritte, in Stunden)
    target_hours_all <- c(0, 12, 24, 36, 48, 60)
    target_hours_all <- target_hours_all[target_hours_all <= max(lead_hours) + 1e-6]
    
    target_times  <- t_now + target_hours_all * 3600
    target_labels <- ifelse(
      target_hours_all == 0,
      paste0(format(target_times, "%d.%m.%Y"), " (Jetzt)"),
      paste0(format(target_times, "%d.%m.%Y"), " (+", target_hours_all, " h)")
    )
    
    ice_fc_layers   <- vector("list", length(target_hours_all))
    climb_fc_layers <- vector("list", length(target_hours_all))
    names(ice_fc_layers)   <- target_labels
    names(climb_fc_layers) <- target_labels
    
    
    # Index 0 h: einfach aktueller Zustand
    if (length(target_hours_all) > 0 && target_hours_all[1] == 0) {
      ice_fc_layers[[1]]   <- ice_thick_expo
      climb_fc_layers[[1]] <- climb_r
    }
    
    # Ausgangs-FDH (bisher) als Raster
    FDH_base_r <- r_FDH_eff
    
    delta_cum_r <- raster(FDH_base_r)
    delta_cum_r[] <- 0
    
    next_target_idx <- which(target_hours_all > 0)[1]
    if (is.na(next_target_idx)) next_target_idx <- length(target_hours_all) + 1
    
    # 8.6 Schleife über alle Vorhersage-Stunden
    for (k in seq_along(time_fc)) {
      
      dt_k <- lead_hours[k]
      
      T_k  <- raster(T2M_fc,  layer = k)
      RH_k <- raster(RH2M_fc, layer = k) / 100
      W_k  <- raster(W_fc,    layer = k)
      
      # Basis-FDH / MDH (°C·h)
      FDH_k <- calc(T_k, fun = function(x) pmax(-x, 0))
      MDH_k <- calc(T_k, fun = function(x) pmax( x, 0))
      
      # Wind-Faktor
      f_wind <- calc(W_k, fun = function(w) {
        wn <- (w - wind_ref) / wind_ref
        f  <- 1 + k_wind * wn
        f[!is.finite(f)] <- 1
        pmax(wind_min, pmin(wind_max, f))
      })
      
      # Feuchte-Faktor (Peak bei rh_opt)
      f_rh <- calc(RH_k, fun = function(rh) {
        peak <- exp(- (rh - rh_opt)^2 / (2 * rh_sig^2))
        0.5 + 0.5 * peak
      })
      
      # Strahlungs-/Expositions-Faktor (Südhang + hohe Sonne -> weniger Wachstum)
      f_rad_k <- calc(solar_index_r, fun = function(s_index) {
        f <- 1 - k_expo * solar_height_factor_fc[k] * s_index
        f[!is.finite(f)] <- 1
        pmax(0, pmin(1, f))
      })
      
      # Zeitgewicht (Zukunft = 1)
      f_time_k <- weight_time_fc[k]
      
      FDH_eff_k <- FDH_k * f_wind * f_rh * f_rad_k * f_time_k
      MDH_eff_k <- MDH_k * k_melt * f_time_k
      
      # kumulative Änderung
      delta_cum_r <- delta_cum_r + (FDH_eff_k - MDH_eff_k)
      
      # Zielzeitpunkte setzen, sobald erreicht
      while (next_target_idx <= length(target_hours_all) &&
             dt_k >= target_hours_all[next_target_idx] - 0.5 &&
             is.null(ice_fc_layers[[next_target_idx]])) {
        
        FDH_fc_k <- FDH_base_r + delta_cum_r
        FDH_fc_k[FDH_fc_k < 0] <- 0
        
        ice_k <- FDH_fc_k * alpha
        ice_k[FDH_fc_k < FDH_crit] <- NA
        names(ice_k) <- paste0("h_pot_expo_m_fc_", target_hours_all[next_target_idx], "h")
        ice_fc_layers[[next_target_idx]] <- ice_k
        
        # Climbability für diesen Zeitschritt
        h_k <- ice_k
        
        score_h_k <- calc(h_k, fun = function(h) {
          s <- (h - h_min) / (h_opt - h_min)
          s[h <= h_min] <- 0
          s[h >= h_opt] <- 1
          s[!is.finite(s)] <- 0
          s
        })
        
        score_T_k <- calc(T_k, fun = function(Tv) {
          s <- 1 - abs(Tv - T_opt) / range_T
          s[Tv <= T_min | Tv >= T_max] <- 0
          s[s < 0] <- 0
          s[!is.finite(s)] <- 0
          s
        })
        
        score_RH_k <- calc(RH_k, fun = function(rh) {
          peak <- exp(- (rh - RH_opt_c)^2 / (2 * RH_sig_c^2))
          s <- peak
          s[!is.finite(s)] <- 0
          s
        })
        
        climb_k <- overlay(score_h_k, score_T_k, score_T3_fc, score_RH_k,
                           fun = function(a, b, c, d) a * b * c * d)
        
        # keine Climbability ohne Eis
        climb_k[is.na(h_k) | h_k <= h_min] <- NA
        climb_k[!is.finite(climb_k)]       <- NA
        climb_k[climb_k <= 0]              <- NA
        names(climb_k) <- paste0("Climbability_fc_", target_hours_all[next_target_idx], "h")
        climb_fc_layers[[next_target_idx]] <- climb_k
        
        next_target_idx <- next_target_idx + 1
      }
    }
    
    # am Ende: nicht gesetzte Targets ggf. mit letztem Zustand füllen
    if (any(vapply(ice_fc_layers, is.null, logical(1)))) {
      last_non_null <- max(which(!vapply(ice_fc_layers, is.null, logical(1))))
      for (i in seq_along(ice_fc_layers)) {
        if (is.null(ice_fc_layers[[i]])) {
          ice_fc_layers[[i]]   <- ice_fc_layers[[last_non_null]]
          climb_fc_layers[[i]] <- climb_fc_layers[[last_non_null]]
        }
      }
    }
    
    fc_step_hours <- target_hours_all
  }
}


# 9) Interaktive Leaflet-Karte (Eisdicke & Climbability – Zeitverlauf) -----

# 9.1 Zeitachsen-Stacks für Historie + Prognose ----------------------------

ice_time_layers   <- list()
climb_time_layers <- list()
time_labels       <- character(0)

# Historische Snapshots (täglich)
if (length(ice_hist_layers) > 0L) {
  ice_time_layers   <- c(ice_time_layers, ice_hist_layers)
  climb_time_layers <- c(climb_time_layers, climb_hist_layers)
  time_labels       <- c(time_labels, ice_hist_labels)
}

# Prognose-Snapshots (0, +12, +24, ...) falls vorhanden
if (!is.null(ice_fc_layers) && length(ice_fc_layers) > 0L) {
  ice_time_layers   <- c(ice_time_layers, ice_fc_layers)
  climb_time_layers <- c(climb_time_layers, climb_fc_layers)
  # Namen der Forecast-Layer wurden bereits als Labels gesetzt
  time_labels       <- c(time_labels, names(ice_fc_layers))
}

# Fallback: falls irgendetwas schiefging, mindestens aktueller Zustand
if (length(ice_time_layers) == 0L) {
  ice_time_layers   <- list(ice_thick_expo)
  climb_time_layers <- list(climb_r)
  time_labels       <- format(Sys.Date(), "%d.%m.%Y (Jetzt)")
}

# Längen konsistent halten
n_steps <- min(length(ice_time_layers), length(climb_time_layers), length(time_labels))
ice_time_layers   <- ice_time_layers[seq_len(n_steps)]
climb_time_layers <- climb_time_layers[seq_len(n_steps)]
time_labels       <- time_labels[seq_len(n_steps)]

# 9.2 Farbpaletten über alle Zeitschritte ----------------------------------

max_h_all <- vapply(
  ice_time_layers,
  function(r) {
    if (is.null(r)) return(NA_real_)
    suppressWarnings(cellStats(r, "max", na.rm = TRUE))
  },
  numeric(1)
)
max_h <- max(max_h_all, na.rm = TRUE)
if (!is.finite(max_h) || max_h <= 0) max_h <- 1

col_fun <- colorRampPalette(c("white", "orange", "green", "darkgreen"))

pal_h <- colorNumeric(
  palette  = col_fun(100),
  domain   = c(0, max_h),
  na.color = "transparent"
)

pal_ci <- colorNumeric(
  palette  = rev(terrain.colors(100)),
  domain   = c(0, 1),
  na.color = "transparent"
)

d_today     <- Sys.Date()
last_update <- format(Sys.time(), "%d.%m.%Y %H:%M UTC")

# 9.3 Leaflet-Karte mit allen Zeitschritten bauen -------------------------

ext <- extent(r_template)

m <- leaflet() |>
  addTiles(group = "OSM")

# Alle Zeitschritte als Layer in zwei Gruppen: "Eisdicke" & "Climbability"
if (length(ice_time_layers) > 0L) {
  visible_step_ice <- length(ice_time_layers)
  for (i in seq_along(ice_time_layers)) {
    r_i <- ice_time_layers[[i]]
    if (!is.null(r_i)) {
      m <- m |>
        addRasterImage(
          r_i,
          colors  = pal_h,
          opacity = if (i == visible_step_ice) 0.8 else 0.0,
          project = TRUE,
          group   = "Eisdicke",
          layerId = paste0("ice_step_", i)
        )
    }
  }
}

if (length(climb_time_layers) > 0L) {
  visible_step_climb <- length(climb_time_layers)
  for (i in seq_along(climb_time_layers)) {
    r_i <- climb_time_layers[[i]]
    if (!is.null(r_i)) {
      m <- m |>
        addRasterImage(
          r_i,
          colors  = pal_ci,
          opacity = if (i == visible_step_climb) 0.7 else 0.0,
          project = TRUE,
          group   = "Climbability",
          layerId = paste0("climb_step_", i)
        )
    }
  }
}

m <- m |>
  fitBounds(
    lng1 = ext@xmin,
    lat1 = ext@ymin,
    lng2 = ext@xmax,
    lat2 = ext@ymax
  ) |>
  addLegend(
    pal       = pal_h,
    values    = c(0, max_h),
    title     = sprintf(
      "Eisdicke (m) (%s)",
      format(d_today, "%d.%m.%Y")
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
      paste0(
        "<div style='font-size: 16px; background: rgba(255,255,255,0.9); padding: 4px 6px; border-radius: 4px; max-width: 300px; line-height: 1.4; margin-top: 4px;'>",
        "<strong>Quellen:</strong> INCA (GeoSphere Austria, ",
        "<a href='https://doi.org/10.60669/6akt-5p05' target='_blank'>doi:10.60669/6akt-5p05</a>); ",
        "DEM Tirol (<a href='https://www.data.gv.at/katalog/datasets/0454f5f3-1d8c-464e-847d-541901eb021a' target='_blank'>data.gv.at</a>); ",
        "NWP AROME (nwp-v1-1h-2500m)<br/>",
        "<strong>Autor:</strong> <a href='https://www.instagram.com/antifascist_mountaineer/' target='_blank'>@antifascist_mountaineer</a><br/>",
        "<em>Letztes Update: ", last_update, "</em>",
        "</div>"
      )
    ),
    position = "bottomleft"
  )

# Debug: Leaflet-Methoden
if (!is.null(m$x$calls)) {
  leaf_methods <- vapply(m$x$calls, `[[`, character(1), "method")
  msg <- paste(
    paste(names(table(leaf_methods)), table(leaf_methods), sep = "="),
    collapse = " | "
  )
  message("Leaflet-Methoden: ", msg)
}

# 9.4 Slider-Logik: Zeitverlauf (Historie + Prognose) ---------------------

if (length(time_labels) > 0L) {
  labels_js  <- paste0("['", paste(time_labels, collapse = "','"), "']")
  n_steps_js <- n_steps
  
  js_code <- sprintf(
    "
    function(el, x) {
      var map    = this;
      window._eisMap = map;
      var labels = %s;
      var nSteps = %d;
      if (!labels || labels.length === 0 || nSteps <= 0) {
        console.log('Zeit-Slider: keine Labels / nSteps', labels, nSteps);
        return;
      }
      if (labels.length !== nSteps) {
        console.log('Zeit-Slider Warnung: labels.length != nSteps', labels.length, nSteps);
      }
      var initial = nSteps - 1;
      if (initial < 0) initial = 0;

      function splitLayersByOrder() {
        var rasters = [];
        map.eachLayer(function(layer) {
          if (!layer || !layer.options) return;
          if (layer instanceof L.TileLayer) return;
          if (typeof layer.setOpacity !== 'function') return;
          rasters.push(layer);
        });
        console.log('splitLayersByOrder: rasters=', rasters.length, 'nSteps=', nSteps);

        var iceLayers   = [];
        var climbLayers = [];
        if (rasters.length >= 2 * nSteps) {
          for (var i = 0; i < nSteps; i++) iceLayers.push(rasters[i]);
          for (var j = 0; j < nSteps; j++) climbLayers.push(rasters[nSteps + j]);
        } else {
          var half = Math.floor(rasters.length / 2);
          for (var i = 0; i < half; i++) iceLayers.push(rasters[i]);
          for (var j = half; j < rasters.length; j++) climbLayers.push(rasters[j]);
          console.log('splitLayersByOrder Fallback: half=', half);
        }
        console.log('splitLayersByOrder: ice=', iceLayers.length, 'climb=', climbLayers.length);
        return { ice: iceLayers, climb: climbLayers };
      }

      var layerCache = splitLayersByOrder();

      function setTimeStep(step) {
        if (step < 0) step = 0;
        if (step >= nSteps) step = nSteps - 1;

        var iceLayers   = layerCache.ice || [];
        var climbLayers = layerCache.climb || [];

        iceLayers.forEach(function(layer, idx) {
          if (layer && typeof layer.setOpacity === 'function') {
            layer.setOpacity(idx === step ? 0.8 : 0.0);
          }
        });

        climbLayers.forEach(function(layer, idx) {
          if (layer && typeof layer.setOpacity === 'function') {
            layer.setOpacity(idx === step ? 0.7 : 0.0);
          }
        });

        var labelDiv = document.getElementById('time-label');
        if (labelDiv && step >= 0 && step < labels.length) {
          labelDiv.textContent = labels[step];
        }

        var slider = document.getElementById('time-slider');
        if (slider && slider.value !== String(step)) {
          slider.value = step;
        }
      }

      var sliderControl = L.control({position: 'topright'});
      sliderControl.onAdd = function() {
        var div = L.DomUtil.create('div', 'info leaflet-control');
        div.style.background   = 'rgba(255,255,255,0.9)';
        div.style.padding      = '8px 10px';
        div.style.borderRadius = '6px';
        div.style.minWidth     = '260px';

        var title = document.createElement('div');
        title.style.fontSize   = '16px';
        title.style.marginBottom = '4px';
        title.innerHTML = '<b>Eisdicke & Climbability – Zeitverlauf</b>';
        div.appendChild(title);

        var labelDiv = document.createElement('div');
        labelDiv.id = 'time-label';
        labelDiv.style.fontSize   = '14px';
        labelDiv.style.marginBottom = '4px';
        labelDiv.textContent = labels[initial] || labels[0];
        div.appendChild(labelDiv);

        var slider = document.createElement('input');
        slider.type  = 'range';
        slider.min   = 0;
        slider.max   = nSteps - 1;
        slider.step  = 1;
        slider.value = initial;
        slider.style.width = '240px';
        slider.id    = 'time-slider';
        div.appendChild(slider);

        slider.addEventListener('input', function(e) {
          var step = parseInt(e.target.value, 10);
          if (!isNaN(step)) {
            console.log('Slider input: step=', step);
            setTimeStep(step);
          }
        });

        slider.addEventListener('mousedown', function(e) {
          if (map && map.dragging) map.dragging.disable();
        });
        slider.addEventListener('mouseup', function(e) {
          if (map && map.dragging) map.dragging.enable();
        });

        L.DomEvent.disableClickPropagation(div);
        return div;
      };
      sliderControl.addTo(map);

      setTimeout(function() {
        setTimeStep(initial);
      }, 500);
    }
    ",
    labels_js,
    n_steps_js
  )
  
  m <- htmlwidgets::onRender(m, js_code)
}

# HTML speichern und Widget zurückgeben
saveWidget(m, "eisdicke_nordtirol.html", selfcontained = TRUE)
