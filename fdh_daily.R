# =====================================================================
# Urheberrechtlicher Hinweis / Nutzungseinsschränkung
# ---------------------------------------------------------------------
# © 2025 Alle Rechte vorbehalten.
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
# =====================================================================

library(httr)
library(raster)
library(leaflet)
library(htmlwidgets)
library(ncdf4)
library(readr)
library(dplyr)
library(tibble)
library(lubridate)

# 0) Zeitraum & Gebiet -------------------------------------------------

start_all    <- as.Date("2025-10-01") # Startdatum der Saison
end_all      <- Sys.Date()           # bis heute
chunk_days   <- 2                    # 2-Tages-Chunks (API-Limit)
chunk_starts <- seq.Date(start_all, end_all, by = chunk_days)

# Nordtirol-BBox in WGS84
bbox <- c(
  46.7,  # lat_min
  10.1,  # lon_min
  47.7,  # lat_max
  12.2   # lon_max
)

# INCA-Parameter (erweitert für Wind, Feuchte etc.)
parameters   <- c("RR", "T2M", "RH2M", "UU", "VV", "GL", "P0", "TD2M")
param_str    <- paste(parameters, collapse = ",")
base_url_inca <- "https://dataset.api.hub.geosphere.at/v1/grid/historical/inca-v1-1h-1km"

out_dir <- "data/inca_nordtirol"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

nc_files <- character(length(chunk_starts))

# 1) INCA-Daten chunkweise laden (existierende Dateien überspringen) ---

for (idx in seq_along(chunk_starts)) {
  cs <- as.Date(chunk_starts[idx])
  ce <- as.Date(min(cs + (chunk_days - 1), end_all))
  
  message("Chunk: ", cs, " bis ", ce)
  
  start_time <- sprintf("%sT00:00", format(cs, "%Y-%m-%d"))
  end_time   <- sprintf("%sT23:00", format(ce, "%Y-%m-%d"))
  
  outfile <- file.path(
    out_dir,
    sprintf("inca_nordtirol_%s_%s.nc",
            format(cs, "%Y%m%d"),
            format(ce, "%Y%m%d"))
  )
  
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
    
    message("  -> gespeichert: ", outfile, " (", file.info(outfile)$size, " Bytes)")
  }
  
  nc_files[idx] <- outfile
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
vars_found    <- intersect(vars_wanted, varnames0)
nc_close(nc0)

# 2.2 Zeitlängen pro Datei
nt_per_file <- integer(length(nc_files))
for (i in seq_along(nc_files)) {
  nc <- nc_open(nc_files[i])
  nt_per_file[i] <- length(ncvar_get(nc, time_dim_name))
  nc_close(nc)
}

nt_total <- sum(nt_per_file)

# 2.3 Arrays allozieren

time_all  <- as.POSIXct(rep(NA_real_, nt_total), origin = "1970-01-01", tz = "UTC")
inca_data <- lapply(vars_found, function(.) array(NA_real_, dim = c(nx, ny, nt_total)))
names(inca_data) <- vars_found

# 2.4 Daten einlesen

pos <- 1L
for (i in seq_along(nc_files)) {
  nc <- nc_open(nc_files[i])
  
  time_vals  <- ncvar_get(nc, time_dim_name)
  time_units <- ncatt_get(nc, time_dim_name, "units")$value
  nt_i       <- length(time_vals)
  
  idx <- pos:(pos + nt_i - 1L)
  time_all[idx] <- convert_nc_time(time_vals, time_units)
  
  for (v in vars_found) {
    inca_data[[v]][ , , idx] <- ncvar_get(nc, v)
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
if (lat_mean_j[1] < tail(lat_mean_j, 1)) {
  idx_j <- ny:1
  lon   <- lon[, idx_j]
  lat   <- lat[, idx_j]
  T2M_arr <- T2M_arr[, idx_j, , drop = FALSE]
  RH_arr  <- RH_arr[,  idx_j, , drop = FALSE]
  UU_arr  <- UU_arr[,  idx_j, , drop = FALSE]
  VV_arr  <- VV_arr[,  idx_j, , drop = FALSE]
}

lon_mean_i <- rowMeans(lon, na.rm = TRUE)
if (lon_mean_i[1] > tail(lon_mean_i, 1)) {
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

FDH_hourly <- pmax(-T2M_arr, 0)      # T < 0
MDH_hourly <- pmax( T2M_arr, 0)      # T > 0
W_arr      <- sqrt(UU_arr^2 + VV_arr^2)   # m/s, [nx, ny, nt]

wind_ref <- mean(W_arr, na.rm = TRUE)
if (!is.finite(wind_ref) || wind_ref <= 0) wind_ref <- 5  # fallback

# 3.3 "Plain FDH" (ohne Korrekturen) + Raster-Template

FDH_mat_plain <- t(apply(FDH_hourly, c(1, 2), sum, na.rm = TRUE)) # [ny, nx]

r_template <- raster(FDH_mat_plain)
extent(r_template) <- c(min(lon, na.rm = TRUE),
                        max(lon, na.rm = TRUE),
                        min(lat, na.rm = TRUE),
                        max(lat, na.rm = TRUE))
crs(r_template) <- "EPSG:4326"

# 4) DEM auf INCA-Grid + Expositionsindex (solar_index_ij) -----------

dem_inca <- raster("data/DEM/DEM_Tirol_INCAgrid_1km_epsg4326.tif")
crs(dem_inca) <- "EPSG:4326"
names(dem_inca) <- "elev_m"

sl_as       <- terrain(dem_inca, opt = c("slope", "aspect"), unit = "degrees")
aspect_inca <- sl_as[["aspect"]]

aspect_rad_from_north <- (aspect_inca - 90) * pi / 180
northness_r           <- cos(aspect_rad_from_north)

# 0 = Nord (schattig), 1 = Süd (besonnt)
solar_index_r  <- (1 - northness_r) / 2
solar_index_r[is.na(solar_index_r[])] <- 0.5
solar_index_ij <- t(as.matrix(solar_index_r))

# 5) Zeitabhängige Sonnenhöhe + Zeit-Gewichtung ----------------------

time_vec <- inca_nordtirol_all$time

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

delta_t <- 23.44 * pi/180 * sin(2 * pi * (284 + doy) / 365)

lat_center_deg <- (min(lat, na.rm = TRUE) + max(lat, na.rm = TRUE)) / 2
lat_center_rad <- lat_center_deg * pi / 180

H_t <- (hour - 12) * 15 * pi / 180

sin_alpha_t <- sin(lat_center_rad) * sin(delta_t) +
  cos(lat_center_rad) * cos(delta_t) * cos(H_t)
sin_alpha_t[sin_alpha_t < 0] <- 0

solar_height_factor_t <- sin_alpha_t

max_time <- max(time_vec, na.rm = TRUE)
age_days <- as.numeric(difftime(max_time, time_vec, units = "days"))

tau_days       <- 10
weight_time_t  <- exp(-age_days / tau_days)

t_start          <- min(time_vec, na.rm = TRUE)
time_offset_hours <- as.numeric(difftime(time_vec, t_start, units = "hours"))

# Lokalzeit (für Hist-Snapshots)
time_local_all <- as.POSIXct(format(time_vec, tz = "Europe/Vienna", usetz = TRUE),
                             tz = "Europe/Vienna")
t_min_local <- min(time_local_all, na.rm = TRUE)
t_max_local <- max(time_local_all, na.rm = TRUE)

t0_local <- as.POSIXct(
  paste0(format(t_min_local, "%Y-%m-%d"), " 07:00:00"),
  tz = "Europe/Vienna"
)
if (t0_local < t_min_local) t0_local <- t0_local + 24 * 3600

t_last_local <- as.POSIXct(
  paste0(format(t_max_local, "%Y-%m-%d"), " 07:00:00"),
  tz = "Europe/Vienna"
)
if (t_last_local > t_max_local) t_last_local <- t_last_local - 24 * 3600

if (t_last_local < t0_local) {
  snap_hours_hist <- numeric(0)
} else {
  snap_times_local <- seq(from = t0_local, to = t_last_local, by = "1 day")
  snap_times_utc   <- as.POSIXct(format(snap_times_local, tz = "UTC", usetz = TRUE),
                                 tz = "UTC")
  snap_hours_hist  <- as.numeric(difftime(snap_times_utc, t_start, units = "hours"))
}

FDH_hist_layers <- list()
FDH_hist_labels <- character(0)
FDH_hist_times  <- as.POSIXct(character(0), tz = "UTC")
snap_idx_hist   <- 1L

# 6) Effektive FDH (FDHm + Wind + Feuchte + Strahlung + Zeit) ---------

k_expo   <- 0.5
k_wind   <- 0.5
wind_min <- 0.5
wind_max <- 2.0
rh_opt   <- 0.7
rh_sig   <- 0.15
k_melt   <- 1.2

FDH_sum_eff <- matrix(0, nrow = nx, ncol = ny)

for (k in seq_len(nt)) {
  T_k  <- T2M_arr[ , , k]
  RH_k <- RH_arr[ , , k] / 100
  
  # 1) Ganzer Zeitschritt ohne Daten -> überspringen (kein Beitrag)
  if (all(is.na(T_k)) || all(is.na(RH_k))) {
    # Historische Snapshots dürfen trotzdem gesetzt werden
    t_k_hours <- time_offset_hours[k]
    while (snap_idx_hist <= length(snap_hours_hist) &&
           t_k_hours >= snap_hours_hist[snap_idx_hist] - 0.5) {
      
      FDH_k_snap <- FDH_sum_eff
      FDH_k_snap[FDH_k_snap < 0] <- 0
      
      r_FDH_snap <- raster(t(FDH_k_snap))
      extent(r_FDH_snap) <- extent(r_template)
      crs(r_FDH_snap)    <- crs(r_template)
      
      label_time <- t_start + snap_hours_hist[snap_idx_hist] * 3600
      label_str  <- format(label_time, "%d.%m.%Y")
      
      FDH_hist_layers[[length(FDH_hist_layers) + 1L]] <- r_FDH_snap
      FDH_hist_labels <- c(FDH_hist_labels, label_str)
      FDH_hist_times  <- c(FDH_hist_times, label_time)
      
      snap_idx_hist <- snap_idx_hist + 1L
    }
    next
  }
  
  # 2) Teilweise fehlende RH -> auf neutralen Optimalwert setzen
  RH_k_eff <- RH_k
  RH_k_eff[is.na(RH_k_eff)] <- rh_opt  # z.B. 0.7
  
  FDH_k <- pmax(-T_k, 0)
  MDH_k <- pmax( T_k, 0)
  W_k   <- W_arr[ , , k]
  
  wind_norm <- (W_k - wind_ref) / wind_ref
  f_wind    <- 1 + k_wind * wind_norm
  f_wind[!is.finite(f_wind)] <- 1
  f_wind <- pmax(wind_min, pmin(wind_max, f_wind))
  
  f_rh_peak <- exp(- (RH_k_eff - rh_opt)^2 / (2 * rh_sig^2))
  f_rh      <- 0.5 + 0.5 * f_rh_peak
  
  s_height <- solar_height_factor_t[k]
  if (s_height == 0) {
    f_rad <- 1
  } else {
    f_rad <- 1 - k_expo * s_height * solar_index_ij
  }
  f_rad[!is.finite(f_rad)] <- 1
  f_rad <- pmax(0, pmin(1, f_rad))
  
  f_time <- weight_time_t[k]
  if (!is.finite(f_time)) f_time <- 1
  
  FDH_eff_k <- FDH_k * f_wind * f_rh * f_rad * f_time
  MDH_eff_k <- MDH_k * k_melt * f_time
  
  FDH_sum_eff <- FDH_sum_eff + (FDH_eff_k - MDH_eff_k)
  
  # Historische FDH-Schnappschüsse (unverändert wie bei dir)
  t_k_hours <- time_offset_hours[k]
  while (snap_idx_hist <= length(snap_hours_hist) &&
         t_k_hours >= snap_hours_hist[snap_idx_hist] - 0.5) {
    
    FDH_k_snap <- FDH_sum_eff
    FDH_k_snap[FDH_k_snap < 0] <- 0
    
    r_FDH_snap <- raster(t(FDH_k_snap))
    extent(r_FDH_snap) <- extent(r_template)
    crs(r_FDH_snap)    <- crs(r_template)
    
    label_time <- t_start + snap_hours_hist[snap_idx_hist] * 3600
    label_str  <- format(label_time, "%d.%m.%Y")
    
    FDH_hist_layers[[length(FDH_hist_layers) + 1L]] <- r_FDH_snap
    FDH_hist_labels <- c(FDH_hist_labels, label_str)
    FDH_hist_times  <- c(FDH_hist_times, label_time)
    
    snap_idx_hist <- snap_idx_hist + 1L
  }
}

FDH_sum_eff[FDH_sum_eff < 0] <- 0

r_FDH_eff <- raster(t(FDH_sum_eff))
extent(r_FDH_eff) <- extent(r_template)
crs(r_FDH_eff)    <- crs(r_template)
names(r_FDH_eff)  <- "FDH_eff_C_h"

# 7) Effektive FDH -> Eisdicke + Climbability-Index (aktuell) ---------

h_c      <- 30
rho_i    <- 880
Lf       <- 334000
FDH_crit <- 1

alpha <- h_c * 3600 / (rho_i * Lf)

# Hilfsfunktion Scores (für mehrfachen Gebrauch)
score_h_fun <- function(h, h_min, h_opt) {
  s <- (h - h_min) / (h_opt - h_min)
  s[h <= h_min] <- 0
  s[h >= h_opt] <- 1
  s[!is.finite(s)] <- 0
  s
}

score_T_fun <- function(Tv, T_opt, T_min, T_max, range_T) {
  s <- 1 - abs(Tv - T_opt) / range_T
  s[Tv <= T_min | Tv >= T_max] <- 0
  s[s < 0] <- 0
  s[!is.finite(s)] <- 0
  s
}

score_RH_fun <- function(rh, RH_opt, RH_sig) {
  peak <- exp(- (rh - RH_opt)^2 / (2 * RH_sig^2))
  s <- peak
  s[!is.finite(s)] <- 0
  s
}

# Historische Eisdicken-Layer
ice_hist_layers <- list()
ice_hist_labels <- character(0)

if (length(FDH_hist_layers) > 0) {
  ice_hist_layers <- vector("list", length(FDH_hist_layers))
  ice_hist_labels <- FDH_hist_labels
  for (i in seq_along(FDH_hist_layers)) {
    r_fd  <- FDH_hist_layers[[i]]
    ice_i <- r_fd * alpha
    ice_i[r_fd < FDH_crit] <- NA
    names(ice_i) <- paste0("h_pot_expo_hist_", ice_hist_labels[i])
    ice_hist_layers[[i]] <- ice_i
  }
}

ice_thick_expo <- r_FDH_eff * alpha
ice_thick_expo[r_FDH_eff < FDH_crit] <- NA
names(ice_thick_expo) <- "h_pot_expo_m"

h_mat_ij <- FDH_sum_eff * alpha

T_last_ij  <- T2M_arr[ , , nt]
RH_last_ij <- RH_arr[ , , nt] / 100

n_tr3   <- min(72, nt)
idx_tr3 <- (nt - n_tr3 + 1):nt

Tr3_ij <- apply(T2M_arr[ , , idx_tr3, drop = FALSE], c(1, 2), mean, na.rm = TRUE)

# 7.1 Dicke-Score
h_min <- 0.10
h_opt <- 0.50
score_h <- score_h_fun(h_mat_ij, h_min, h_opt)

# 7.2 Aktuelle Temperatur (T_last)
T_opt  <- -4
T_min  <- -20
T_max  <- 0
range_T <- max(T_opt - T_min, T_max - T_opt)
score_T <- score_T_fun(T_last_ij, T_opt, T_min, T_max, range_T)

# 7.3 Tr3 (3-Tages-Mittel)
T3_opt  <- -6
T3_min  <- -20
T3_max  <- -1
range_T3 <- max(T3_opt - T3_min, T3_max - T3_opt)
score_T3 <- score_T_fun(Tr3_ij, T3_opt, T3_min, T3_max, range_T3)

# 7.4 Luftfeuchte
RH_opt_c <- 0.7
RH_sig_c <- 0.2
score_RH <- score_RH_fun(RH_last_ij, RH_opt_c, RH_sig_c)

# Climbability-Historie (tägliche Snapshots)
climb_hist_layers <- list()
climb_hist_labels <- character(0)

if (length(FDH_hist_layers) > 0 && length(FDH_hist_labels) == length(FDH_hist_layers)) {
  climb_hist_layers <- vector("list", length(FDH_hist_layers))
  climb_hist_labels <- FDH_hist_labels
  
  for (i in seq_along(FDH_hist_layers)) {
    r_fd <- FDH_hist_layers[[i]]
    h_i  <- r_fd * alpha
    h_i[r_fd < FDH_crit] <- NA
    
    t_i    <- FDH_hist_times[i]
    dt_vec <- abs(as.numeric(difftime(time_vec, t_i, units = "hours")))
    k_i    <- which.min(dt_vec)
    if (!length(k_i) || !is.finite(k_i)) next
    
    T_i_mat  <- t(T2M_arr[ , , k_i])
    RH_i_mat <- t(RH_arr[ , , k_i] / 100)
    
    r_T_i    <- raster(T_i_mat);  extent(r_T_i)  <- extent(r_template); crs(r_T_i)  <- crs(r_template)
    r_RH_i   <- raster(RH_i_mat); extent(r_RH_i) <- extent(r_template); crs(r_RH_i) <- crs(r_template)
    
    k_start   <- max(1, k_i - (72 - 1))
    idx3      <- k_start:k_i
    Tr3_i_arr <- apply(T2M_arr[ , , idx3, drop = FALSE], c(1, 2), mean, na.rm = TRUE)
    Tr3_i_mat <- t(Tr3_i_arr)
    r_Tr3_i   <- raster(Tr3_i_mat); extent(r_Tr3_i) <- extent(r_template); crs(r_Tr3_i) <- crs(r_template)
    
    score_h_i <- calc(h_i, fun = function(h) score_h_fun(h, h_min, h_opt))
    score_T_i <- calc(r_T_i, fun = function(Tv) score_T_fun(Tv, T_opt, T_min, T_max, range_T))
    score_T3_i <- calc(r_Tr3_i, fun = function(T3v) score_T_fun(T3v, T3_opt, T3_min, T3_max, range_T3))
    score_RH_i <- calc(r_RH_i, fun = function(rh) score_RH_fun(rh, RH_opt_c, RH_sig_c))
    
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
climb_index_ij[h_mat_ij <= h_min]          <- NA
climb_index_ij[!is.finite(climb_index_ij)] <- NA
climb_index_ij[climb_index_ij <= 0]        <- NA

climb_r <- raster(t(climb_index_ij))
extent(climb_r) <- extent(r_template)
crs(climb_r)    <- crs(r_template)
names(climb_r)  <- "Climbability_0_1"

# 8) Prognose-Eisdicke & -Climbability aus NWP ------------------------

base_url_nwp   <- "https://dataset.api.hub.geosphere.at/v1/grid/forecast/nwp-v1-1h-2500m"
parameters_nwp <- c("t2m", "rh2m", "u10m", "v10m")
param_str_nwp  <- paste(parameters_nwp, collapse = ",")

t_now    <- max(inca_nordtirol_all$time, na.rm = TRUE)
t0       <- t_now
fc_h_max <- 60

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

ice_fc_layers   <- NULL
climb_fc_layers <- NULL
fc_step_hours   <- NULL

if (file.exists(outfile_fc)) {
  nc_fc <- nc_open(outfile_fc)
  
  time_dim_name_fc <- names(nc_fc$dim)[grep("time", tolower(names(nc_fc$dim)))[1]]
  time_fc_vals     <- ncvar_get(nc_fc, time_dim_name_fc)
  time_fc_units    <- ncatt_get(nc_fc, time_dim_name_fc, "units")$value
  time_fc          <- convert_nc_time(time_fc_vals, time_fc_units)
  
  nc_close(nc_fc)
  
  lead_hours <- as.numeric(difftime(time_fc, t_now, units = "hours"))
  keep_fc    <- which(lead_hours >= 0 & lead_hours <= fc_h_max + 0.01)
  
  if (length(keep_fc) > 0) {
    time_fc    <- time_fc[keep_fc]
    lead_hours <- lead_hours[keep_fc]
    
    T2M_fc_brick  <- brick(outfile_fc, varname = "t2m")[[keep_fc]]
    RH2M_fc_brick <- brick(outfile_fc, varname = "rh2m")[[keep_fc]]
    U10_fc_brick  <- brick(outfile_fc, varname = "u10m")[[keep_fc]]
    V10_fc_brick  <- brick(outfile_fc, varname = "v10m")[[keep_fc]]
    
    crs(T2M_fc_brick)  <- crs(r_template)
    crs(RH2M_fc_brick) <- crs(r_template)
    crs(U10_fc_brick)  <- crs(r_template)
    crs(V10_fc_brick)  <- crs(r_template)
    
    T2M_fc  <- resample(T2M_fc_brick,  r_template, method = "bilinear")
    RH2M_fc <- resample(RH2M_fc_brick, r_template, method = "bilinear")
    U10_fc  <- resample(U10_fc_brick,  r_template, method = "bilinear")
    V10_fc  <- resample(V10_fc_brick,  r_template, method = "bilinear")
    
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
    
    weight_time_fc <- rep(1, length(time_fc))
    
    Tr3_fc <- calc(T2M_fc, fun = function(x) mean(x, na.rm = TRUE))
    
    score_T3_fc <- calc(Tr3_fc, fun = function(t3) {
      score_T_fun(t3, T3_opt, T3_min, T3_max, range_T3)
    })
    
    time_fc_local <- as.POSIXct(
      format(time_fc, tz = "Europe/Vienna", usetz = TRUE),
      tz = "Europe/Vienna"
    )
    
    t0_local_fc <- as.POSIXct(
      format(t0, tz = "Europe/Vienna", usetz = TRUE),
      tz = "Europe/Vienna"
    )
    
    # -------------------------------
    # Stundenweise Zeitstufen (Prognose)
    # 0 h  = Analysezeit (t0)
    # >0 h = jede NWP-Stunde
    # -------------------------------
    target_hours_all  <- c(0, lead_hours)
    target_labels_all <- c(
      paste0(format(t0_local_fc, "%d.%m.%Y, %H:%M"), " (Analyse)"),
      format(time_fc_local, "%d.%m.%Y, %H:%M")
    )
    
    
    ord <- order(target_hours_all)
    target_hours_all  <- target_hours_all[ord]
    target_labels_all <- target_labels_all[ord]
    
    max_lead <- max(lead_hours) + 1e-6
    keep_idx <- target_hours_all <= max_lead
    target_hours_all  <- target_hours_all[keep_idx]
    target_labels_all <- target_labels_all[keep_idx]
    
    target_times <- t0 + target_hours_all * 3600
    
    # Zuordnung jeder Ziel-Zeitstufe zum nächstgelegenen NWP-Index (k)
    step_to_k <- integer(length(target_hours_all))
    if (length(target_hours_all) > 0) {
      for (ii in seq_along(target_hours_all)) {
        diffs <- abs(lead_hours - target_hours_all[ii])
        if (all(!is.finite(diffs))) {
          step_to_k[ii] <- NA_integer_
        } else {
          step_to_k[ii] <- which.min(diffs)
        }
      }
    }
    
    ice_fc_layers   <- vector("list", length(target_hours_all))
    climb_fc_layers <- vector("list", length(target_hours_all))
    names(ice_fc_layers)   <- target_labels_all
    names(climb_fc_layers) <- target_labels_all
    
    if (length(target_hours_all) > 0 && target_hours_all[1] == 0) {
      ice_fc_layers[[1]]   <- ice_thick_expo
      climb_fc_layers[[1]] <- climb_r
    }
    
    FDH_base_r  <- r_FDH_eff
    delta_cum_r <- raster(FDH_base_r); delta_cum_r[] <- 0
    
    next_target_idx <- which(target_hours_all > 0)[1]
    if (is.na(next_target_idx)) next_target_idx <- length(target_hours_all) + 1
    
    for (k in seq_along(time_fc)) {
      dt_k <- lead_hours[k]
      
      T_k  <- raster(T2M_fc,  layer = k)
      RH_k <- raster(RH2M_fc, layer = k) / 100
      W_k  <- raster(W_fc,    layer = k)
      
      FDH_k <- calc(T_k, fun = function(x) pmax(-x, 0))
      MDH_k <- calc(T_k, fun = function(x) pmax( x, 0))
      
      f_wind <- calc(W_k, fun = function(w) {
        wn <- (w - wind_ref) / wind_ref
        f  <- 1 + k_wind * wn
        f[!is.finite(f)] <- 1
        pmax(wind_min, pmin(wind_max, f))
      })
      
      f_rh <- calc(RH_k, fun = function(rh) {
        0.5 + 0.5 * exp(- (rh - rh_opt)^2 / (2 * rh_sig^2))
      })
      
      s_height_k <- solar_height_factor_fc[k]
      f_rad_k <- if (s_height_k == 0) {
        solar_index_r*0 + 1
      } else {
        calc(solar_index_r, fun = function(s_index) {
          f <- 1 - k_expo * s_height_k * s_index
          f[!is.finite(f)] <- 1
          pmax(0, pmin(1, f))
        })
      }
      
      f_time_k <- weight_time_fc[k]
      
      FDH_eff_k <- FDH_k * f_wind * f_rh * f_rad_k * f_time_k
      MDH_eff_k <- MDH_k * k_melt * f_time_k
      
      delta_cum_r <- delta_cum_r + (FDH_eff_k - MDH_eff_k)
      
      while (next_target_idx <= length(target_hours_all) &&
             dt_k >= target_hours_all[next_target_idx] - 0.5 &&
             is.null(ice_fc_layers[[next_target_idx]])) {
        
        FDH_fc_k <- FDH_base_r + delta_cum_r
        FDH_fc_k[FDH_fc_k < 0] <- 0
        
        ice_k <- FDH_fc_k * alpha
        ice_k[FDH_fc_k < FDH_crit] <- NA
        names(ice_k) <- paste0("h_pot_expo_m_fc_", target_hours_all[next_target_idx], "h")
        ice_fc_layers[[next_target_idx]] <- ice_k
        
        h_k <- ice_k
        
        score_h_k <- calc(h_k, fun = function(h) score_h_fun(h, h_min, h_opt))
        score_T_k <- calc(T_k, fun = function(Tv) score_T_fun(Tv, T_opt, T_min, T_max, range_T))
        score_RH_k <- calc(RH_k, fun = function(rh) score_RH_fun(rh, RH_opt_c, RH_sig_c))
        
        climb_k <- overlay(score_h_k, score_T_k, score_T3_fc, score_RH_k,
                           fun = function(a, b, c, d) a * b * c * d)
        
        climb_k[is.na(h_k) | h_k <= h_min] <- NA
        climb_k[!is.finite(climb_k)]       <- NA
        climb_k[climb_k <= 0]              <- NA
        names(climb_k) <- paste0("Climbability_fc_", target_hours_all[next_target_idx], "h")
        climb_fc_layers[[next_target_idx]] <- climb_k
        
        next_target_idx <- next_target_idx + 1
      }
    }
    
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

# === Kompakte Tagesübersicht mit Piktogrammen + Kurz-Labels =========

best_html <- NULL

if (exists("climb_fc_layers") &&
    !is.null(climb_fc_layers) &&
    length(climb_fc_layers) > 0 &&
    exists("target_times") &&
    length(target_times) == length(climb_fc_layers)) {
  
  # 1) Alle Zeiten in lokale Zeit umrechnen
  times_local_all <- lubridate::with_tz(target_times, "Europe/Vienna")
  dates_all       <- as.Date(times_local_all)
  
  n_steps <- length(climb_fc_layers)
  
  # 2) Mittlere Kletterbarkeit pro Zeitschritt
  mean_all <- rep(NA_real_, n_steps)
  for (i in seq_len(n_steps)) {
    r_i <- climb_fc_layers[[i]]
    if (!is.null(r_i)) {
      mean_all[i] <- suppressWarnings(
        raster::cellStats(r_i, "mean", na.rm = TRUE)
      )
    }
  }
  
  # Dataframe für saubere Klassen
  df_fc <- data.frame(
    step       = seq_len(n_steps),
    time_local = times_local_all,
    date       = dates_all,
    mean_ci    = mean_all,
    stringsAsFactors = FALSE
  )
  
  # Nur gültige Zeilen
  df_fc <- df_fc[is.finite(df_fc$mean_ci), ]
  
  if (nrow(df_fc) > 0) {
    # Icon je nach Tageszeit
    hour_icon <- function(h) {
      if (h < 6 || h >= 22) {
        "🌙"
      } else if (h < 9) {
        "🌅"
      } else if (h < 18) {
        "☀️"
      } else {
        "🌇"
      }
    }
    
    # Temperatur-Tags (inkl. „evtl. zu kalt“)
    temp_tag <- function(T_mean, best = TRUE) {
      if (!is.finite(T_mean)) return(NA_character_)
      
      if (T_mean <= -18) {
        # sehr tiefe T -> sprödes Eis, eher heikel
        if (best) "evtl. zu kalt / spröde" else "sehr kalt, spröde"
      } else if (T_mean <= -10) {
        if (best) "sehr kalt" else "sehr kalt"
      } else if (T_mean <= -5) {
        if (best) "kalt, meist gut" else "kalt"
      } else if (T_mean <= -2) {
        # Optimum T_climbing ≈ -4…-3 °C
        if (best) "Optimalbereich (-4…-3 °C)" else "nahe Optimum"
      } else if (T_mean <= 0) {
        if (best) "mild, eher weich" else "mild"
      } else {
        "zu mild (>0 °C)"
      }
    }
    
    # Feuchte-Tags (40–50 % ideal)
    rh_tag <- function(RH_mean) {
      if (!is.finite(RH_mean)) return(NA_character_)
      
      if (RH_mean < 0.30) {
        "zu trocken"
      } else if (RH_mean < 0.40) {
        "eher trocken"
      } else if (RH_mean <= 0.50) {
        "opt. Feuchte (40–50 %)"
      } else if (RH_mean <= 0.70) {
        "leicht feucht"
      } else {
        "sehr feucht"
      }
    }
    
    build_short_tag <- function(T_mean, RH_mean, best = TRUE) {
      tt <- temp_tag(T_mean, best)
      rr <- rh_tag(RH_mean)
      tags <- na.omit(c(tt, rr))
      if (length(tags) == 0) return("")
      paste(tags, collapse = ", ")
    }
    
    pict_rows <- character(0)
    
    # Nach Tagen gruppieren
    split_by_date <- split(df_fc, df_fc$date)
    
    for (d_str in names(split_by_date)) {
      block <- split_by_date[[d_str]]
      d <- as.Date(d_str)  # explizit -> kein trim-Fehler
      
      # Beste / schlechteste Zeit an diesem Tag
      idx_best <- which.max(block$mean_ci)
      idx_min  <- which.min(block$mean_ci)
      
      row_best <- block[idx_best, ]
      row_min  <- block[idx_min, ]
      
      t_best <- row_best$time_local
      t_min  <- row_min$time_local
      
      h_best <- lubridate::hour(t_best) + lubridate::minute(t_best) / 60
      h_min  <- lubridate::hour(t_min)  + lubridate::minute(t_min)  / 60
      
      icon_best <- hour_icon(h_best)
      icon_min  <- hour_icon(h_min)
      
      # Kurz-Labels aus Prognose-T und RH holen (falls möglich)
      tag_best <- ""
      tag_min  <- ""
      
      if (exists("T2M_fc") && exists("RH2M_fc") && exists("step_to_k")) {
        # BESTE Zeit
        k_best <- NA_integer_
        if (!is.null(step_to_k) &&
            row_best$step <= length(step_to_k)) {
          k_best <- step_to_k[row_best$step]
        }
        if (is.finite(k_best) &&
            k_best >= 1 &&
            k_best <= raster::nlayers(T2M_fc)) {
          T_r_best  <- raster::raster(T2M_fc,  layer = k_best)
          RH_r_best <- raster::raster(RH2M_fc, layer = k_best) / 100
          
          T_mean_best  <- suppressWarnings(raster::cellStats(T_r_best,  "mean", na.rm = TRUE))
          RH_mean_best <- suppressWarnings(raster::cellStats(RH_r_best, "mean", na.rm = TRUE))
          
          tag_best <- build_short_tag(T_mean_best, RH_mean_best, best = TRUE)
        }
        
        # SCHLECHTESTE Zeit
        k_min <- NA_integer_
        if (!is.null(step_to_k) &&
            row_min$step <= length(step_to_k)) {
          k_min <- step_to_k[row_min$step]
        }
        if (is.finite(k_min) &&
            k_min >= 1 &&
            k_min <= raster::nlayers(T2M_fc)) {
          T_r_min  <- raster::raster(T2M_fc,  layer = k_min)
          RH_r_min <- raster::raster(RH2M_fc, layer = k_min) / 100
          
          T_mean_min  <- suppressWarnings(raster::cellStats(T_r_min,  "mean", na.rm = TRUE))
          RH_mean_min <- suppressWarnings(raster::cellStats(RH_r_min, "mean", na.rm = TRUE))
          
          tag_min <- build_short_tag(T_mean_min, RH_mean_min, best = FALSE)
        }
      }
      
      # Tabellenzeile:
      #  - erste Zeile: Uhrzeit + Icon
      #  - zweite Zeile klein: Kurz-Labels (inkl. "evtl. zu kalt" / "zu mild" etc.)
      pict_rows <- c(
        pict_rows,
        sprintf(
          "<tr>
             <td>%s</td>
             <td>
               ▲ %s&nbsp;%s<br/>
               <span style='font-size:10px;'>%s</span>
             </td>
             <td>
               ▼ %s&nbsp;%s<br/>
               <span style='font-size:10px;'>%s</span>
             </td>
           </tr>",
          format(d, "%d.%m."),
          format(t_best, "%H:%M"), icon_best, tag_best,
          format(t_min,  "%H:%M"), icon_min,  tag_min
        )
      )
    }
    
    if (length(pict_rows) > 0) {
      table_html <- paste0(
        "<table style='width:100%; border-collapse:collapse; font-size:12px;'>",
        "<tr>",
        "<th style='text-align:left;'>Tag</th>",
        "<th style='text-align:left;'>Beste Zeit</th>",
        "<th style='text-align:left;'>Schlechteste Zeit</th>",
        "</tr>",
        paste(pict_rows, collapse = ""),
        "</table>"
      )
      
      best_html <- htmltools::HTML(
        paste0(
          "<div id='climb-summary-box' ",
          "style='font-size: 13px; background: rgba(255,255,255,0.9);",
          "padding: 6px 8px; border-radius: 6px; max-width: 420px;",
          "line-height: 1.4;'>",
          
          # Header: anklickbar, Body zuerst versteckt
          "<div id='climb-summary-header' ",
          "style='font-weight:bold; cursor:pointer; margin-bottom:4px;'>",
          "Kletterbarkeit – Tagesübersicht (Prognose) ▾",
          "</div>",
          
          "<div id='climb-summary-body' style='display:none;'>",
          table_html,
          "<div style='font-size:11px; margin-top:4px;'>",
          "▲ höchster, ▼ niedrigster mittlerer Climbability-Wert (Gebietsmittel).",
          "<br/>",
          "<span style='font-size:10px;'>",
          "\"evtl. zu kalt\" bedeutet: sehr tiefe Lufttemperaturen (z.B. &lt; −18 °C), ",
          "Eis kann spröde sein.",
          "</span>",
          "</div>",
          "</div>",  # Ende body
          "</div>"   # Ende Box
        )
      )
    }
  }
}


# 9) Interaktive Leaflet-Karte (Eisdicke & Climbability – Zeitverlauf) -----

ice_time_layers   <- list()
climb_time_layers <- list()
time_labels       <- character(0)

if (length(ice_hist_layers) > 0L) {
  ice_time_layers   <- c(ice_time_layers, ice_hist_layers)
  climb_time_layers <- c(climb_time_layers, climb_hist_layers)
  time_labels       <- c(time_labels, ice_hist_labels)
}

if (!is.null(ice_fc_layers) && length(ice_fc_layers) > 0L) {
  ice_time_layers   <- c(ice_time_layers, ice_fc_layers)
  climb_time_layers <- c(climb_time_layers, climb_fc_layers)
  time_labels       <- c(time_labels, names(ice_fc_layers))
}

if (length(ice_time_layers) == 0L) {
  ice_time_layers   <- list(ice_thick_expo)
  climb_time_layers <- list(climb_r)
  time_labels       <- format(Sys.Date(), "%d.%m.%Y (Jetzt)")
}

n_steps <- min(length(ice_time_layers), length(climb_time_layers), length(time_labels))
ice_time_layers   <- ice_time_layers[seq_len(n_steps)]
climb_time_layers <- climb_time_layers[seq_len(n_steps)]
time_labels       <- time_labels[seq_len(n_steps)]

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

col_fun <- colorRampPalette(c("#ffffff", "#c6dbef", "#6baed6", "#08519c"))

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

last_update <- format(Sys.time(), "%d.%m.%Y %H:%M UTC")
ext         <- extent(r_template)

# 9a) Eisfall-Sonnendaten + vorhandene Topo-URLs ----------------------

sun_df <- readr::read_csv(
  "data/Koordinaten_Wasserfaelle/icefalls_sun_horizon_2024_2025.csv",
  show_col_types = FALSE
)

# Sicherstellen, dass sunrise_topo / sunset_topo als Datum-Zeit vorliegen
# (falls readr sie als Character eingelesen hat)
if (!inherits(sun_df$sunrise_topo, "POSIXt")) {
  sun_df <- sun_df %>%
    mutate(
      sunrise_topo = ymd_hms(sunrise_topo, tz = "UTC"),
      sunset_topo  = ymd_hms(sunset_topo,  tz = "UTC")
    )
}

sun_df <- sun_df %>%
  mutate(
    # CSV-Zeiten sind als UTC gespeichert -> zurück in lokale Zeit
    sunrise_topo = with_tz(sunrise_topo, "Europe/Vienna"),
    sunset_topo  = with_tz(sunset_topo,  "Europe/Vienna"),
    # Sonnendauer zur Sicherheit neu aus den korrigierten Zeiten berechnen
    sun_hours_topo = as.numeric(difftime(sunset_topo, sunrise_topo, units = "hours"))
  )

 sun_date <- Sys.Date()        

sun_today <- sun_df %>%
  dplyr::filter(date == sun_date)

# Fallback: falls kein Eintrag genau zu sun_date existiert -> letztes Datum
if (nrow(sun_today) == 0) {
  last_date <- max(sun_df$date, na.rm = TRUE)
  message("Keine Sonnendaten für ", sun_date,
          " gefunden. Verwende stattdessen ", last_date, ".")
  sun_today <- sun_df %>%
    dplyr::filter(date == last_date)
}

if (nrow(sun_today) == 0) {
  stop("sun_today ist leer – prüfe icefalls_sun_horizon_2024_2025.csv (Spaltennamen / date-Typ).")
}

sun_today <- sun_today %>%
  dplyr::mutate(
    # times so verwenden, wie sie in der CSV stehen
    sunrise_txt   = substr(as.character(sunrise_topo), 12, 16),  # "YYYY-MM-DD HH:MM:SS" -> HH:MM
    sunset_txt    = substr(as.character(sunset_topo),  12, 16),
    sun_hours_txt = sprintf("%.1f", sun_hours_topo),
    date_txt      = format(date, "%d.%m.%Y"),
    # Link-Text separat bauen
    link_txt = ifelse(
      !is.na(topo_url) & topo_url != "",
      paste0("<a href='", topo_url, "' target='_blank'>Topo öffnen</a>"),
      "(kein Topo-Link hinterlegt)"
    ),
    # Popup-Text: bei NAs -> "keine direkte Sonneneinstrahlung"
    popup = ifelse(
      is.na(sun_hours_topo) | is.na(sunrise_topo) | is.na(sunset_topo),
      sprintf(
        "<b>%s</b><br/>Sonne am %s: keine direkte Sonneneinstrahlung<br/>%s",
        name,
        date_txt,
        link_txt
      ),
      sprintf(
        "<b>%s</b><br/>Sonne am %s: %s – %s (%s h)<br/>%s",
        name,
        date_txt,
        sunrise_txt,
        sunset_txt,
        sun_hours_txt,
        link_txt
      )
    )
  )

m <- leaflet() |>
  addProviderTiles(providers$OpenStreetMap, group = "OSM") |>
  addProviderTiles(providers$OpenTopoMap,   group = "Gelände (Topo)")

if (!is.null(best_html)) {
  m <- m |>
    addControl(
      position = "bottomright", 
      html     = best_html
    )
}

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
    pal       = pal_ci,
    values    = c(0, 1),
    title     = "Climbability",
    labFormat = labelFormat(digits = 1),
    position  = "bottomleft"
  ) |>
  addLegend(
    pal       = pal_h,
    values    = c(0, max_h),
    title     = "Eisdicke (m)",
    labFormat = labelFormat(digits = 2),
    position  = "bottomleft"
  ) |>
  # === Eisfall-Markierungen ==========================================
addCircleMarkers(
  data        = sun_today,
  lng         = ~longitude,
  lat         = ~latitude,
  radius      = 5,
  color       = "black",
  weight      = 1,
  fillColor   = "orange",
  fillOpacity = 0.9,
  popup       = ~popup,
  group       = "Eisfälle"
) |>
  # ===================================================================
addLayersControl(
  baseGroups    = c("OSM", "Gelände (Topo)"),
  overlayGroups = c("Eisdicke", "Climbability", "Eisfälle"),
  options       = layersControlOptions(collapsed = FALSE)
) |>
  addControl(
    position = "bottomright",
    html = htmltools::HTML(
      paste0(
        "<div id='impressum-box' ",
        "style='font-size: 14px; background: rgba(255,255,255,0.9);",
        "padding: 6px 8px; border-radius: 6px; max-width: 340px;",
        "line-height: 1.4;'>",
        "<div id='impressum-header' ",
        "style='font-weight:bold; cursor:pointer; margin-bottom:4px;'>",
        "Impressum / Quellen ▾",
        "</div>",
        "<div id='impressum-body' style='display:none;'>",
        "<strong>Quellen</strong><br/>",
        "INCA (GeoSphere Austria, ",
        "<a href='https://doi.org/10.60669/6akt-5p05' target='_blank'>doi:10.60669/6akt-5p05</a>)<br/>",
        "NWP AROME (GeoSphere Austria, ",
        "<a href='https://doi.org/10.60669/9zm8-s664' target='_blank'>doi:10.60669/9zm8-s664</a>)<br/>",
        "DEM Tirol (",
        "<a href='https://www.data.gv.at/katalog/datasets/0454f5f3-1d8c-464e-847d-541901eb021a' target='_blank'>data.gv.at</a>)<br/><br/>",
        "<em>Letztes Update: ", last_update, "</em>",
        "</div>",
        "</div>"
      )
    )
  )

m <- htmlwidgets::onRender(
  m,
  "
  function(el, x) {
    var header = el.querySelector('#impressum-header');
    var body   = el.querySelector('#impressum-body');
    if (!header || !body) return;

    header.addEventListener('click', function() {
      var visible = body.style.display !== 'none';
      body.style.display = visible ? 'none' : 'block';
      header.innerHTML = visible
        ? 'Impressum / Quellen ▾'
        : 'Impressum / Quellen ▴';
    });
  }
  "
)

if (length(time_labels) > 0L) {
  labels_js  <- paste0("['", paste(time_labels, collapse = "','"), "']")
  n_steps_js <- n_steps
  
  js_code <- sprintf(
    "
  function(el, x) {
    var map    = this;
    var lc = el.getElementsByClassName('leaflet-control-layers-expanded')[0]
             || el.getElementsByClassName('leaflet-control-layers')[0];
    if (lc) {
       lc.style.marginTop   = '10px';
       lc.style.marginRight = '90px';
       lc.style.transform   = 'scale(1.5)';
       lc.style.transformOrigin = 'top left';
       lc.style.padding     = '12px 15px';
       lc.style.fontSize    = '16px';
    }

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

      var iceLayers   = [];
      var climbLayers = [];
      if (rasters.length >= 2 * nSteps) {
        for (var i = 0; i < nSteps; i++) iceLayers.push(rasters[i]);
        for (var j = 0; j < nSteps; j++) climbLayers.push(rasters[nSteps + j]);
      } else {
        var half = Math.floor(rasters.length / 2);
        for (var i = 0; i < half; i++) iceLayers.push(rasters[i]);
        for (var j = half; j < rasters.length; j++) climbLayers.push(rasters[j]);
      }
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
      div.style.marginTop    = '70px';
      div.style.marginRight  = '10px';

      var title = document.createElement('div');
      title.style.fontSize    = '16px';
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

m <- htmlwidgets::onRender(
  m,
  "
  function(el, x) {
    var header = el.querySelector('#climb-summary-header');
    var body   = el.querySelector('#climb-summary-body');
    if (!header || !body) return;

    // sicherheitshalber beim Laden zugeklappt
    body.style.display = 'none';

    header.addEventListener('click', function() {
      var visible = body.style.display !== 'none';
      body.style.display = visible ? 'none' : 'block';
      header.innerHTML = visible
        ? 'Kletterbarkeit – Tagesübersicht (Prognose) ▾'
        : 'Kletterbarkeit – Tagesübersicht (Prognose) ▴';
    });
  }
  "
)

saveWidget(m, "eisdicke_nordtirol.html", selfcontained = TRUE)
