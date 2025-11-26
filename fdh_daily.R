library(ncdf4)
library(httr)
library(raster)
library(leaflet)
library(htmlwidgets)

# Daten Herunterladen ####

# ---- 0) Zeitraum ##

start_all <- as.Date("2025-11-01")
end_all   <- Sys.Date()          

chunk_days <- 4

# Startdaten der Chunks explizit als Date-Sequenz
chunk_starts <- seq.Date(start_all, end_all, by = chunk_days)

# ---- 1) Nordtirol-BBox in WGS84 

bbox <- c(
  46.7,  # lat_min (south)
  10.1,  # lon_min (west)
  47.7,  # lat_max (north)
  12.2   # lon_max (east)
)

# ---- 2) Parameter & Endpoint #########-

parameters <- c("T2M", "TD2M", "RH2M", "RR", "GL")
param_str  <- paste(parameters, collapse = ",")

base_url <- "https://dataset.api.hub.geosphere.at/v1/grid/historical/inca-v1-1h-1km"

out_dir <- "data/inca_nordtirol_nov"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

nc_files <- character(0)

# ---- 3) Schleife über 4-Tages-Chunks ######--

for (cs in chunk_starts) {
  cs <- as.Date(cs)
  ce <- as.Date(min(cs + (chunk_days - 1), end_all))
  
  message("Chunk: ", as.character(cs), " bis ", as.character(ce))
  
  start_time <- paste0(base::format(cs, "%Y-%m-%d"), "T00:00")
  end_time   <- paste0(base::format(ce, "%Y-%m-%d"), "T23:00")  # 23:00, nicht 24:00
  
  # Dateiname pro Chunk
  outfile <- file.path(
    out_dir,
    sprintf("inca_nordtirol_%s_%s.nc",
            base::format(cs, "%Y%m%d"),
            base::format(ce, "%Y%m%d"))
  )
  
  query <- list(
    parameters    = param_str,
    start         = start_time,
    end           = end_time,
    bbox          = paste(bbox, collapse = ","),  # "46.7,10.1,47.7,12.2"
    output_format = "netcdf",
    filename      = sprintf("inca_nordtirol_%s_%s",
                            base::format(cs, "%Y%m%d"),
                            base::format(ce, "%Y%m%d"))
  )
  
  resp <- GET(
    url   = base_url,
    query = query,
    write_disk(outfile, overwrite = TRUE)
  )
  
  stop_for_status(resp)
  
  sz <- file.info(outfile)$size
  message("  -> gespeichert: ", outfile, " (", sz, " Bytes)")
  nc_files <- c(nc_files, outfile)
}




# ---- 4) Alle Chunks einlesen und zusammenfügen ###-
# ###
# 4) Alle Chunks mit ncdf4 einlesen & zusammenfügen
#    -> lon/lat + alle Variablen ("T2M", "TD2M", "RH2M", "RR", "GL")
# ###

# gewünschte Variablen (so wie du sie bei parameters übergibst)
vars_wanted <- c("T2M", "TD2M", "RH2M", "RR", "GL")

# kleine Hilfsfunktion: NetCDF-Zeit in POSIXct umrechnen
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

nc_files <- sort(nc_files)

# ### 4.1 Erster Durchlauf: Dimensionen & nt je Chunk bestimmen ###

# 1. Datei öffnen -> lon/lat, Dim
nc0 <- nc_open(nc_files[1])

varnames0 <- names(nc0$var)

# lon/lat-Variablen finden (hier minimal, ohne viel Kontrolle)
lon_name <- intersect(c("lon", "longitude", "x", "xc"), varnames0)[1]
lat_name <- intersect(c("lat", "latitude", "y", "yc"), varnames0)[1]

lon0 <- ncvar_get(nc0, lon_name)
lat0 <- ncvar_get(nc0, lat_name)

nx <- dim(lon0)[1]
ny <- dim(lon0)[2]

# Zeit-Dim
time_dim_name <- names(nc0$dim)[grep("time", tolower(names(nc0$dim)))[1]]
time_vals0    <- ncvar_get(nc0, time_dim_name)
time_units0   <- ncatt_get(nc0, time_dim_name, "units")$value

# Variablen, die wir wirklich haben
vars_found <- intersect(vars_wanted, varnames0)

nc_close(nc0)

# nt je Chunk holen
nt_per_file <- integer(length(nc_files))
for (i in seq_along(nc_files)) {
  nc <- nc_open(nc_files[i])
  time_vals <- ncvar_get(nc, time_dim_name)
  nt_per_file[i] <- length(time_vals)
  nc_close(nc)
}

nt_total <- sum(nt_per_file)

# ### 4.2 Arrays vorallozieren ---

# Zeitvektor
time_all <- as.POSIXct(rep(NA_real_, nt_total), origin = "1970-01-01", tz = "UTC")

# Daten-Arrays für jede Variable: [nx, ny, nt_total]
inca_data <- lapply(vars_found, function(.) array(NA_real_, dim = c(nx, ny, nt_total)))
names(inca_data) <- vars_found

# ### 4.3 Zweiter Durchlauf: Daten reinschreiben ---

pos <- 1
for (i in seq_along(nc_files)) {
  nc <- nc_open(nc_files[i])
  
  # Zeit
  time_vals  <- ncvar_get(nc, time_dim_name)
  time_units <- ncatt_get(nc, time_dim_name, "units")$value
  nt_i       <- length(time_vals)
  
  idx <- pos:(pos + nt_i - 1)
  time_all[idx] <- convert_nc_time(time_vals, time_units)
  
  # alle gewünschten Variablen einlesen
  for (v in vars_found) {
    arr_i <- ncvar_get(nc, v)    # [nx, ny, nt_i]
    inca_data[[v]][ , , idx] <- arr_i
  }
  
  nc_close(nc)
  pos <- pos + nt_i
}

# ### 4.4 Alles in ein Objekt packen ---

inca_nordtirol_all <- list(
  lon  = lon0,       # 2D: [nx, ny]
  lat  = lat0,       # 2D: [nx, ny]
  time = time_all,   # POSIXct: Länge nt_total
  data = inca_data   # Liste mit Arrays [nx, ny, nt_total] pro Variable
)

# Kurzer Check, aber ohne viel Output
str(inca_nordtirol_all, max.level = 1)




# Speichern
library(ncdf4)

# Ziel-Dateiname
out_nc <- "data/inca_nordtirol_all.nc"

lon  <- inca_nordtirol_all$lon
lat  <- inca_nordtirol_all$lat
time <- inca_nordtirol_all$time
data_list <- inca_nordtirol_all$data   # Liste mit RR, T2M, RH2M, UU, VV

nx <- dim(lon)[1]
ny <- dim(lon)[2]
nt <- length(time)

# ==== Dimensionen definieren (Index-Gitter + Zeit) ===

dim_i <- ncdim_def("i", "", 1:nx)   # x-Index
dim_j <- ncdim_def("j", "", 1:ny)   # y-Index

time_units <- "seconds since 1970-01-01 00:00:00"
time_vals  <- as.numeric(time - as.POSIXct("1970-01-01 00:00:00", tz = "UTC"))
dim_t <- ncdim_def("time", time_units, time_vals, unlim = TRUE)

# ==== Variablen-Definitionen: lon, lat 

var_lon <- ncvar_def(
  name   = "lon",
  units  = "degrees_east",
  dim    = list(dim_i, dim_j),
  missval = NA_real_,
  prec   = "double"
)

var_lat <- ncvar_def(
  name   = "lat",
  units  = "degrees_north",
  dim    = list(dim_i, dim_j),
  missval = NA_real_,
  prec   = "double"
)

# === Daten-Variablen (RR, T2M, RH2M, UU, VV)

var_defs <- list(var_lon, var_lat)

for (vn in names(data_list)) {
  vdef <- ncvar_def(
    name   = vn,
    units  = "",  # falls du willst: hier echte Units eintragen
    dim    = list(dim_i, dim_j, dim_t),
    missval = NA_real_,
    prec   = "float"
  )
  var_defs[[length(var_defs) + 1]] <- vdef
}

# ==== NetCDF erzeugen und füllen ===

nc <- nc_create(out_nc, var_defs)

# lon/lat schreiben
ncvar_put(nc, "lon", lon)
ncvar_put(nc, "lat", lat)

# alle Datenvariablen schreiben
for (vn in names(data_list)) {
  ncvar_put(nc, vn, data_list[[vn]])
}

nc_close(nc)


# Modell berechnen ###################################################################################################################
## 1) FDH berechnen ----------------------------------------

T2M_arr <- inca_nordtirol_all$data$T2M   # [nx, ny, nt]

FDH_hourly <- -T2M_arr
FDH_hourly[FDH_hourly < 0] <- 0

FDH_sum <- apply(FDH_hourly, c(1, 2), sum, na.rm = TRUE)  # [nx, ny]


## 2) Lon/Lat & Orientierung -------------------------------

lon <- inca_nordtirol_all$lon   # [nx, ny]
lat <- inca_nordtirol_all$lat   # [nx, ny]

nx <- dim(FDH_sum)[1]
ny <- dim(FDH_sum)[2]

lat_mean_j <- colMeans(lat, na.rm = TRUE)
if (lat_mean_j[1] < lat_mean_j[length(lat_mean_j)]) {
  idx_j   <- ny:1
  FDH_sum <- FDH_sum[, idx_j]
  lon     <- lon[, idx_j]
  lat     <- lat[, idx_j]
}

lon_mean_i <- rowMeans(lon, na.rm = TRUE)
if (lon_mean_i[1] > lon_mean_i[length(lon_mean_i)]) {
  idx_i   <- nx:1
  FDH_sum <- FDH_sum[idx_i, ]
  lon     <- lon[idx_i, ]
  lat     <- lat[idx_i, ]
}


## 3) FDH-Raster mit richtiger Extent ----------------------

FDH_mat <- t(FDH_sum)

r_FDH <- raster(FDH_mat)
extent(r_FDH) <- c(min(lon, na.rm = TRUE),
                   max(lon, na.rm = TRUE),
                   min(lat, na.rm = TRUE),
                   max(lat, na.rm = TRUE))
crs(r_FDH) <- "EPSG:4326"
names(r_FDH) <- "FDH_C_h"


## 4) Potentielle Eisdicke --------------------------------

h_c      <- 30       # W m^-2 K^-1
rho_i    <- 880      # kg m^-3
Lf       <- 334000   # J kg^-1
FDH_crit <- 50       # °C·h

alpha <- h_c * 3600 / (rho_i * Lf)   # m / (°C·h)

ice_thick_pot <- r_FDH * alpha
names(ice_thick_pot) <- "h_pot_m"

ice_thick_pot[r_FDH < FDH_crit] <- NA


## 5) Plot mit gewünschter Farbskala ------------------

max_h   <- cellStats(ice_thick_pot, max, na.rm = TRUE)
col_fun <- colorRampPalette(c("white", "orange", "green", "darkgreen"))

pal <- colorNumeric(
  palette = col_fun(100),
  domain  = c(0, max_h),
  na.color = "transparent"
)

# --- Leaflet-Karte bauen -------------------------
m <- leaflet() |>
  addTiles(group = "OSM") |>
  addRasterImage(
    ice_thick_pot,
    colors  = pal,
    opacity = 0.8,
    project = TRUE
  ) |>
  addLegend(
    pal    = pal,
    values = c(0, max_h),
    title  = "Eisdicke (m)",
    labFormat = labelFormat(digits = 2)
  )

# --- HTML speichern ------------------------------
out_file <- "eisdicke_nordtirol.html"
saveWidget(m, out_file, selfcontained = TRUE)

