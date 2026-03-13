##### Preprocessing VectAbundance + EID datasets v.1 10.03.26 #####

library(tidyverse)
library(sf)
library(fuzzySim)
library(purrr)
library(lubridate)

# Load datasets
data_VA <- read_csv("Vectabundace_v015.csv") # 42853 rows
data_nice <- data_all # 11667 rows

# Remove Cote Azur from VA dataset
VA_rem_nice <- data_VA |> 
  filter_out(Region == "Cote Azur") # 42540 rows

# Remove irrelevant columns to match Nice dataset
VA_clean <- VA_rem_nice |> 
  transmute(
    ID, 
    endDate = as.Date(date),
    startDate = as.Date(date - 7),
    year = isoyear(date),
    week = isoweek(date), 
    count = as.numeric(value),
    country = Country,
    region = Region,
    lat = latitude,
    long = longitude,
    EPSG
  ) # 42540 rows

# Convert lat long coords in Nice dataset
# helper: convert mixed DD / DMS column to numeric DD
toDD <- function(x, coord = c("lat","lon")) {
  coord <- match.arg(coord)
  if (is.numeric(x)) return(x)
  
  x <- as.character(x) |> str_trim()
  
  # identify DMS-like rows (degree symbol)
  is_dms <- str_detect(x, "°")
  is_dms[is.na(is_dms)] <- FALSE
  
  out <- rep(NA_real_, length(x))
  
  # Non-DMS: parse as decimal number
  if (any(!is_dms)) {
    out[!is_dms] <- parse_number(x[!is_dms])
  }
  
  # DMS: convert
  if (any(is_dms)) {
    x_dms <- x[is_dms] |>
      str_replace_all(",", ".") |>
      str_replace_all("''", "\"") |>
      str_replace_all("”|“", "\"") |>
      str_squish()
    
    has_hemi <- str_detect(x_dms, "[NSEW]")
    has_hemi[is.na(has_hemi)] <- FALSE
    
    # 1) With hemisphere letters: normal conversion
    if (any(has_hemi)) {
      out[is_dms][has_hemi] <- dms2dec(x_dms[has_hemi])
    }
    
    # 2) Without hemisphere letters: assume + (N for lat, E for lon)
    if (any(!has_hemi)) {
      x_nohemi <- x_dms[!has_hemi]
      # append assumed hemisphere so dms2dec behaves
      assumed <- if (coord == "lat") "N" else "E"
      x_nohemi2 <- paste(x_nohemi, assumed)
      vals <- dms2dec(x_nohemi2)
      out[is_dms][!has_hemi] <- abs(vals)  # enforce positive
    }
  }
  
  out
}

# Clean Nice dataset
Nice_clean <- data_nice  |> 
  transmute(
    ID = id_piege,
    endDate = as.Date(date_detection),
    startDate = as.Date(date_detection - trapping_days),
    year = isoyear(date_detection),
    week = isoweek(date_detection),
    count = as.numeric(eggs_per_day*7),
    country = "France",
    region = "Cote Azur",
    lat = round(toDD(lat), 5),
    long = round(toDD(lon), 5),
    EPSG = "EPSG"
  ) # 11667 rows

# Merge datasets
all_noeco <- rbind(VA_clean, Nice_clean) # 54207 rows
write_csv(all_noeco, file.path("VectAbundance_2024","VectAbundance_2024_raw.csv"))

# Remove NA lat long
missing_na <- all_noeco |>  filter_out(!is.na(lat), !is.na(long)) # 7757 rows
all_coords_na <- all_noeco |> filter(!is.na(lat), !is.na(long)) # 46450 rows
write_csv(all_coords_na, file.path("VectAbundance_2024", "VectAbundance_noNAcoords.csv"))

# Load bioecoregions
bioecoregions <- st_read("BiogeoRegions2016.sqlite") |> st_transform(4326)

# Create geometry points for the data
all_coords <- st_as_sf(all_coords_na, coords = c("long", "lat"), crs = 4326, remove = FALSE)

# Intersect points
all_intersect <- st_intersection(all_coords,bioecoregions)  |> st_drop_geometry() # 44135 rows

# Troubleshoot rows removed
nrow(all_coords) - nrow(all_intersect) # 2315 rows removed
missing_all_intersect <- anti_join(all_coords, all_intersect, by = c("ID", "endDate")) # 2315 rows
missing_all_intersect |> distinct(ID, region, lat, long, endDate) |> count(ID, region, lat, long, endDate) # check which rows are on the water using QGIS

# Mutate missing rows
missing_intersect_fix <- bind_rows(
  missing_all_intersect |>
    filter(ID %in% c("11867", "15856", "16090", "15353", "15588", "06N29")) |>
    mutate(short_name = "mediterranean", pre_2012 = "MED", code = "Mediterranean", name = "Mediterranean Bio-geographical Region"),
  missing_all_intersect |>
    filter(ID %in% c("8189", "9116")) |>
    mutate(short_name = "continental", pre_2012   = "CON", code = "Continental", name = "Continental Bio-geographical Region")
  ) |>
  st_drop_geometry()

# Join with others
all_intersect_final <- bind_rows(all_intersect, missing_intersect_fix)

# Save as csv
write_csv(all_intersect_final, file.path("VectAbundance_2024", "VectAbundance_2024_rawecoreg.csv"))

# Remove NA trap counts
count_NA <- all_intersect_final |> filter(!is.na(count)) # 30780 rows, 13355 rows with NA counts
write_csv(count_NA, file.path("VectAbundance_2024", "VectAbundance_2024_noNA.csv"))

# Keep only nonzero trap counts
count_nonzero <- count_NA |> filter_out(count <= 0) |>  mutate(count = round(count, 3)) # 14366 rows
write_csv(count_nonzero, file.path("VectAbundance_2024", "VectAbundance_2024_nonzero.csv"))

# Keep only counts above the threshold (POT from Da Re et al. 2025) which is 55
within_threshold <- count_nonzero |> filter_out(count < 55) # 8388 rows
write_csv(within_threshold, file.path("VectAbundance_2024", "VectAbundance_2024_POT.csv"))

#### simple visualization of trap data by region
data_aedes <- count_NA 
data_aedes |> 
  group_by(year, week, code)  |> 
  summarise(count = mean(count))|> 
  mutate(date = as.Date(paste(as.numeric(year), as.numeric(week), 1, sep="-"), "%Y-%U-%u")) |> 
  ggplot(aes(x=date, y = count, group = code, )) + 
  geom_line() + 
  facet_wrap(.~code, nrow = 3)

#### simple visualization of trap data above threshold (55 counts) by region
data_aedes2 <- within_threshold |> 
  group_by(year, week, code)  |> 
  summarise(count = mean(count))|> 
  mutate(date = as.Date(paste(as.numeric(year), as.numeric(week), 1, sep="-"), "%Y-%U-%u")) |> 
  ggplot(aes(x=date, y = count, group = code, )) + 
  geom_line() + 
  facet_wrap(.~code, nrow = 3)
