## Preprocessing to CCM pipeline v.130526
library(tidyverse)
library(sf)
library(fuzzySim)
library(purrr)
library(lubridate)
library(dplyr)
library(tidyr)
library(ISOweek)
install.packages("remotes")
library(remotes)
#remotes::install_github("wranglezone/tibblify", dependencies = TRUE, force = TRUE)
#remotes::install_github("tpisel/openmeteo", dependencies = TRUE, force = TRUE)
library(tibblify)
library(openmeteo) # openmeteo_0.2.4.tar.gz
library(scales)
library(patchwork)
library(cowplot)
library(grid)
library(ggh4x)
library(devtools)
#devtools::install_github("Nmoiroux/ecoXCorr")
library(ecoXCorr)
library(readr)

datasetName <- "VectAbundance_2024"

# 1. Preprocessing the VectAbundance dataset ----
## Load datasets ----
data_VA <- read_csv("Vectabundace_v015.csv") # 42853 rows

## Identify duplicates ----
# Define column groups
all_col <- colnames(data_VA)
ignore_col <- c("Institute", "contact_person", "contact_person_email", "DOI_website",
                "volume", "larvicide_presence", "larvicide_type")
main_col <- setdiff(all_col, c(ignore_col, "value"))

# Identify and flag duplicates
dup_flag <- data_VA |> 
  ### Dup 1: exact duplicates across columns ----
group_by(across(all_of(all_col))) |> 
  mutate(dup_1 = n() > 1) |> 
  relocate(dup_1, .before = ID) |> 
  ungroup() |> 
  
  ### Dup 2: duplicates across all except provider details ----
group_by(across(all_of(setdiff(all_col, ignore_col)))) |> 
  mutate(dup_2 = n() > 1) |> 
  relocate(dup_2, .after = dup_1) |> 
  ungroup() |> 
  
  ### Dup 3: duplicates across all but with a zero and nonzero value ---- 
group_by(across(all_of(main_col))) |>
  mutate(has_zero     = any(value == 0, na.rm = TRUE),
         has_nonzero  = any(value > 0, na.rm = TRUE),
         has_na       = any(is.na(value)),
         dup_3 = n() > 1 & ((has_zero & has_nonzero) |
                              (has_nonzero & has_na)   |
                              (has_zero & has_na))) |>
  relocate(dup_3, .after = dup_2) |> 
  ungroup() |> 
  
  ### Dup 4: duplicates across all but with differing nonzero values ----
group_by(across(all_of(main_col))) |>
  mutate(dup_4 = n() > 1 & n_distinct(value[value !=0]) > 1) |>
  relocate(dup_4, .after = dup_3) |>
  ungroup()

# From Daniele, remove Duplicate types 1-3
# For Duplicate type 4, remove trap ID 7258 (basically remove also)
VA_clean <- dup_flag |> # start with 42853 rows
  # Remove Duplicate type 4 (ID 7258 = 578 rows())
  filter(ID != 7258) |> # 42275 rows
  
  # Remove Duplicate type 3 (but they are from ID 7258 anyway)
  group_by(across(all_of(main_col))) |>
  filter(!(dup_3 & (value == 0 | # drop zero if mixed
                      (is.na(value) & has_nonzero) | # drop NA if nonzero exists
                      (is.na(value) & has_zero & !has_nonzero)))) |> # drop NA if only zero+NA
  ungroup() |> 
  
  # Remove Duplicates type 1 and 2
  distinct(across(all_of(setdiff(all_col, ignore_col))), .keep_all = TRUE) # 41596 rows

## Remove irrelevant columns ----
VA_clean <- VA_clean |> 
  transmute(ID,
            endDate = as.Date(date),
            startDate = as.Date(date - 7),
            year = isoyear(date),
            week = isoweek(date),
            count = as.numeric(value),
            country = Country,
            region = Region,
            lat = latitude,
            long = longitude,EPSG) |> 
  filter(!is.na(lat), !is.na(long)) #|> # 42741 rows

## Load bioecoregions ---
bioecoregions <- st_read("BiogeoRegions2016.sqlite") |> st_transform(4326)

# Create geometry points for the data
all_coords <- st_as_sf(VA_clean, coords = c("long", "lat"), crs = 4326, remove = FALSE)

## Intersect points ----
all_intersect <- st_intersection(all_coords,bioecoregions)  |> st_drop_geometry() # 40443 rows

# Troubleshoot rows removed
nrow(all_coords) - nrow(all_intersect) # 2298 rows removed
missing_all_intersect <- anti_join(all_coords, all_intersect, by = c("ID", "endDate")) # 2298 rows
missing_all_intersect |> distinct(ID, region, lat, long, endDate) |> count(ID, region, lat, long, endDate)
unique(missing_all_intersect$ID) # 7 traps missing
# check with QGIS after

## Mutate missing rows ----
missing_intersect_fix <- bind_rows(
  missing_all_intersect |>
    filter(ID %in% c("11867", "15353", "15588", "15856", "16090")) |>
    mutate(short_name = "mediterranean", pre_2012 = "MED", code = "Mediterranean", name = "Mediterranean Bio-geographical Region"),
  missing_all_intersect |>
    filter(ID %in% c("8189", "9116")) |>
    mutate(short_name = "continental", pre_2012 = "CON", code = "Continental", name = "Continental Bio-geographical Region")
) |>
  st_drop_geometry()

# Join with others
all_intersect_final <- bind_rows(all_intersect, missing_intersect_fix) # 42741 rows again

# Save as csv
# write_csv(all_intersect_final, file.path("VectAbundance_2024", "VectAbundance_2024_rawecoreg.csv"))

## Remove NA trap counts ----
count_noNA <- all_intersect_final |> filter(!is.na(count)) # 29159 rows, 12959 rows with NA counts
# write_csv(count_noNA, file.path("VectAbundance_2024", "VectAbundance_2024_noNA.csv"))

## Keep only nonzero trap counts ----
count_nonzero <- count_NA |> filter_out(count <= 0) |>  mutate(count = round(count, 3)) # 13302 rows
# write_csv(count_nonzero, file.path("VectAbundance_2024", "VectAbundance_2024_nonzero.csv"))

# 2. OpenMeteo retrieval ----
## Set up ----
# Prepare data
data_clim <- count_nonzero |> # do also with count_noNA
  transmute(ID, endDate, startDate, ecoreg = code, year, week, count, lat, long)
# data_clim <- read_csv(file.path(datasetName, "VectAbundance_2024_noNA.csv")) # for importing

## Set parameters ----
maxLag <- 12
daily <- "relative_humidity_2m_mean"

## Create a new output directory ----
# Keep in mind to also repeat for count_NA data, to plot full abundance timeseries
outdir_clim <- file.path(datasetName, "clim_daily_noNA", daily) 
dir.create(outdir_clim, recursive = TRUE, showWarnings = FALSE)

## Batch divide ----
# Identify start and end dates per ID
idDates <- data_clim |> 
  group_by(ID, lat, long, ecoreg) |> 
  summarise(
    idStart = min(startDate, na.rm = TRUE), # earliest sampling start
    idEnd   = max(endDate,   na.rm = TRUE), # latest trap date
    .groups = "drop")  |> 
  mutate(startLag   = idStart - weeks(maxLag),
         yearStart  = year(startLag),
         yearEnd    = year(idEnd))

sites <- idDates |>
  distinct(ID, lat, long, ecoreg)

# Create site batches based on years
idBatch <- idDates |> 
  rowwise() |> 
  mutate(years = list(seq(yearStart, yearEnd))) |> 
  unnest(years) |> 
  ungroup() |> 
  transmute(
    ID, lat, long, ecoreg,
    batchStart = as.character(as.Date(paste0(years, "-01-01"))),
    batchEnd = as.character(as.Date(paste0(years, "-12-31"))))

## Function to retrieve climate data ----
unlink("openmeteo_cache", recursive = TRUE)
dir.create("openmeteo_cache", showWarnings = FALSE)

cachePath <- function(ID, start, end) {
  file.path("openmeteo_cache", paste0("site_", ID, "__", start, "__", end, ".rds"))
}

retrieveBatchCached <- function(lat, long, start, end, ID, 
                                ecoreg = NULL, 
                                maxTries = 8,
                                baseWait = 0.8,
                                maxWait  = 60) {
  f <- cachePath(ID, start, end)
  if (file.exists(f)) return(readRDS(f)) # If already downloaded, load from disk
  wait <- baseWait
  for (k in seq_len(maxTries)) {
    Sys.sleep(wait)
    out <- tryCatch(
      weather_history( # function from OpenMeteo package
        location = c(lat, long),
        start = start,
        end   = end,
        daily = daily), # check this
      error = function(e) e)
    if (!inherits(out, "error")) {
      out$ID <- ID
      out$ecoreg <- ecoreg
      out$lat <- lat
      out$long <- long
      out$batchStart <- start
      out$batchEnd <- end
      saveRDS(out, f)
      return(out)
    }
    msg <- conditionMessage(out)
    message("FAILED ID=", ID, " ", start, " to ", end, " try=", k, " : ", msg)
    if (grepl("429", msg)) {
      wait <- min(max(wait * 2, 2), maxWait)
    } 
    else {wait <- min(wait * 1.5, maxWait)}
  }
  message("GAVE UP ID=", ID, " ", start, " to ", end)
  NULL
}

## Run function ----
climList <- map(seq_len(nrow(idBatch)), \(i) {
  retrieveBatchCached(ID = idBatch$ID[i],
                      ecoreg = idBatch$ecoreg[i],
                      lat   = idBatch$lat[i],
                      long  = idBatch$long[i],
                      start = idBatch$batchStart[i],
                      end   = idBatch$batchEnd[i])})

# Bind rows
climSites <- bind_rows(compact(climList))
#climSites <- climSites |> left_join(site_lookup, by = "ID")

# Save as .csv for backup
write.csv(climSites, file = file.path(outdir_clim, paste0(datasetName, "_", daily, "_clim.csv")), row.names = FALSE)

## Compile into table (lat/long already in climSites) ----
# Pick the date column
dateCol <- if ("date" %in% names(climSites)) {
  "date"
} else {
  intersect(names(climSites), c("time", "datetime", "date_raw"))[1]
}

if (is.na(dateCol)) {
  stop("Could not find a usable date column in climSites. Columns are: ",
       paste(names(climSites), collapse = ", "))
}

climSites <- read_csv(file.path(datasetName, "clim_daily_noNA", "temperature_2m_max", "VectAbundance_2024_temperature_2m_max_clim.csv"))

# Build daily table (keep lat/long)
climDaily <- climSites |>
  rename(dateRaw = all_of(dateCol)) |>
  transmute(
    ID,
    ecoreg,
    lat,
    long,
    date = as.Date(dateRaw),
    max_temp = daily_temperature_2m_max)

# De-duplicate overlapping year chunks
climDaily <- climDaily |>
  distinct(ID, date, .keep_all = TRUE) |>
  arrange(ID, date)

# Check if all dates have climate values
message("Missing climate values: ", sum(is.na(climDaily$max_temp)))

# Save as .csv (this is already "spatial")
write.csv(climDaily, 
          file = file.path(outdir_clim, paste0(datasetName, "_", daily, "_climSpatial.csv")),
          row.names = FALSE)

## Compile climate data ----
c1 <- read_csv(file.path(datasetName, "clim_daily_noNA", "temperature_2m_mean", "VectAbundance_2024_temperature_2m_mean_climSpatial.csv"))
c2 <- read_csv(file.path(datasetName, "clim_daily_noNA", "temperature_2m_max", "VectAbundance_2024_temperature_2m_max_climSpatial.csv"))
c3 <- read_csv(file.path(datasetName, "clim_daily_noNA", "temperature_2m_min", "VectAbundance_2024_temperature_2m_min_climSpatial.csv"))
c4 <- read_csv(file.path(datasetName, "clim_daily_noNA", "precipitation_sum", "VectAbundance_2024_precipitation_sum_climSpatial.csv"))
c5 <- read_csv(file.path(datasetName, "clim_daily_noNA", "relative_humidity_2m_mean", "VectAbundance_2024_relative_humidity_2m_mean_climSpatial.csv"))
data_clim <- c1 |> 
  left_join(c2, by = c("ID", "ecoreg", "lat", "long", "date")) |> 
  left_join(c3, by = c( "ID", "ecoreg", "lat", "long", "date")) |> 
  left_join(c4, by = c("ID", "ecoreg", "lat", "long", "date")) |> 
  left_join(c5, by = c("ID", "ecoreg", "lat", "long", "date")) |> 
  mutate(year = year(date),
         week = isoweek(date)) |> 
  relocate(year, .after = date) |> 
  relocate(week, .before = mean_temp)
#write_csv(data_clim, file.path("VectAbundance_2024", "VectAbundance_2024_allclimnoNA.csv"))

## Add GDD column ----
with_gdd <- data_clim |> 
  mutate(min_temp_adj = pmin(pmax(min_temp, 11), 30),
         max_temp_adj = pmin(pmax(max_temp, 11), 30),
         gdd = ((max_temp_adj + min_temp_adj) / 2) - 11) |>
  select(-min_temp_adj, -max_temp_adj) |> 
  relocate(gdd, .after = min_temp)
# Overwrite previous csv
write_csv(with_gdd, file.path("VectAbundance_2024", "VectAbundance_2024_allclimnoNA.csv"))

# Simple visualization of trap data by site
plot_aedes <- count_nonzero |> # also check count_noNA
  rename(ecoreg = code) |> 
  #filter(code == "Mediterranean") |> 
  group_by(ecoreg, year, week)  |> 
  summarise(count = mean(count))|> 
  mutate(date = as.Date(paste(as.numeric(year), as.numeric(week), 1, sep="-"), "%Y-%U-%u")) |> 
  ggplot(aes(x=date, y = count, group = ecoreg, )) + 
  geom_line() + 
  facet_wrap(.~ecoreg, nrow = 3)
plot_aedes

# At this point, we can create visualizations of the abundance and climate timeseries

# 3. Cross correlation maps ----
## Set up ----
### Set output directories ----
outdir_ccm <- file.path(datasetName, "CCM_outputs")
out_csv <- file.path(outdir_ccm, "fitmodels")
out_png <- file.path(outdir_ccm, "maps")
dir.create(out_csv, recursive = TRUE, showWarnings = FALSE)
dir.create(out_png, recursive = TRUE, showWarnings = FALSE)

### Load biological data ----
data_aedes <- count_nonzero |>
  transmute(ID = as.character(ID),
            date = as.Date(endDate),
            year, week,
            count,
            count_c = ceiling(count),
            ecoreg = code, lat, long)
# For imported data:
# data_aedes <- read_csv(file.path(datasetName, "VectAbundance_2024_nonzero.csv")) |>
#   transmute(ID = as.character(ID),
#             date = as.Date(endDate),
#             year, week,
#             count,
#             count_c = ceiling(count),
#             ecoreg = code, lat, long)


### Load climate data ----
data_clim <- with_gdd |> 
  mutate(ID = as.character(ID))
# For imported data:
# data_clim <- read_csv(file.path(datasetName, "VectAbundance_2024_allclimnonzero.csv")) |>
#   mutate(ID = as.character(ID))

# Make sure both IDs are the same class (character)

### Set parameters ----
interval <- 7
max_lag  <- 8
response <- "count_c"
random   <- "(1|ID)"
family   <- "truncated_nbinom2"
threshold_p <- 0.2 
clim_vars <- c("mean_temp", "max_temp", "min_temp", "gdd", "sum_precip", "mean_rh")
ecoregs   <- c("Global", "Alpine", "Continental", "Mediterranean")

## Aggregate climate data -----
# Filter by ecoregion - for Global model, no need to filter
eco <- "Continental" # change the ecoregion here, but also change OBJECT NAMES
aedes_eco <- filter(data_aedes, ecoreg == eco)
clim_eco <- filter(data_clim, ecoreg == eco)
sampling_dates_con <- unique(aedes_eco$date) # like this
# sampling_dates_glob <- unique(data_aedes$date) # for Global model

# split by ID and aggregate
clim_agg_list <- list()
for (id in unique(data_aedes$ID)) { # aedes_eco for Local model, data_aedes for Global model
  message("Processing ID:", id)
  
  clim_sub <- clim_eco |> filter(ID == id) # for Local model
  aedes_sub <- aedes_eco |> filter(ID == id)
  # clim_sub <- data_clim |> filter(ID == id) # for Global model
  # aedes_sub <- data_aedes |> filter(ID == id)
  
  # aggregate per ID using aggregate_lagged_intervals
  agg <- aggregate_lagged_intervals(
    data       = clim_sub,
    date_col   = "date",
    value_cols = clim_vars,
    ref_date   = unique(aedes_sub$date),
    interval   = interval,
    max_lag    = max_lag,
    id_col     = NULL,   # IMPORTANT: already split
    funs       = list(mean = mean, sum = sum),
    na.rm      = TRUE)
  
  agg$ID <- id
  clim_agg_list[[id]] <- agg
}
# Combine all IDs
clim_agg_con <- bind_rows(clim_agg_list) |> 
  mutate(ecoreg = eco) |> # for Local model
  # left_join(data_aedes |> distinct(ID, ecoreg), by = "ID") |> # for Global model
  relocate(ID, ecoreg, .before = date)

# Merge with count data
merged_data_con <- clim_agg_con |>
  left_join(aedes_eco, by = c("ID", "date"), # for Local model
            relationship = "many-to-many") |>
  # left_join(data_aedes, by = c("ID", "date", "ecoreg"), # for global
  #           relationship = "many-to-many") |>
  filter(!is.na(count_c)) |>
  select(-mean_temp_sum, -max_temp_sum, -min_temp_sum,
         -gdd_mean, -sum_precip_mean, -mean_rh_sum) |>
  mutate(ecoreg = eco) |> # for Local model
  relocate(ecoreg, year, week, count, count_c, lat, long, .after = date)

# Sanity check
all(paste(aedes_eco$ID, aedes_eco$date) %in% paste(merged_data_con$ID, merged_data_con$date)) # for Local model
# all(paste(data_aedes$ID, data_aedes$date) %in% paste(merged_data_glob$ID, merged_data_glob$date)) # for Global model

# save as .csv for later
write_csv(merged_data_con, file.path(out_csv, paste0(eco, "_aggclim.csv"))) # for Local model
# write_csv(merged_data_glob, file.path(out_csv, "Global_aggclim.csv")) # for Global model

## Fit models ----
predictors <- c("mean_temp_mean", "max_temp_mean", "min_temp_mean", "gdd_sum", "sum_precip_sum", "mean_rh_mean")
fits_tmp_con <- list() # store results here
for (pred in predictors) {
  message("Running predictor:", pred)
  # fit models
  fit_tbl <- fit_models_by_lag(
    data = merged_data_con,
    response = response,
    predictors = pred,
    random = random,
    family = family)
  # store in list
  fits_tmp_con[[pred]] <- fit_tbl
  #save as csv
  write_csv(fit_tbl, file.path(out_csv, paste0(eco, "_fitmodels_", pred, ".csv"))) # for Local model
  #write_csv(fit_tbl, file.path(out_csv, paste0("Global_fitmodels_", pred, ".csv"))) # for Global model
}

## 3: CREATE CCM PLOTS ##
plots_tmp_con <- list()
for (pred in names(fits_tmp_con)){
  fit <- fits_tmp_con[[pred]]
  # generate CCM
  p <- plotCCM(fit,
               model_outcome = "R2sign",
               threshold_p = 0.2) +
    labs(title = paste0(eco, " - ", pred), # for Local model
      #title = paste0("Global - ", pred), # for Global model
         x = "Lag start (weeks)", y = "Lag end (weeks)") +
    theme(plot.title = element_text(face = "bold"))

  #extract best lag directly from the plot
  best_lag <- p$layers[[2]]$data
  label <- paste0("r(", best_lag$lag_start, ",", best_lag$lag_end, ") = ",
                  formatC(best_lag$value, digits = 3, format = "f"))
  p <- p +
    annotate("text", x = 1.6, y = max_lag - 1.8,
             label = label, fontface = "bold",
             hjust = 1, size = 5)
  # store in the list
  plots_tmp_con[[pred]] <- p
  print(p)
  # save as png
  ggsave(filename = file.path(out_png, paste0(eco, "_", pred, ".png")), # for Local model
         #filename = file.path(out_png, paste0("Global - ", pred, ".png")), # for Global model
         plot = p, width = 8, height = 6, dpi = 300)
}

# combine all fits
fits_all <- list(Global = fits_tmp_glob, Alpine = fits_tmp_alp, Continental = fits_tmp_con, Mediterranean = fits_tmp_med)
ccm_data_list <- list()
ccm_best_list <- list()

for (eco in names(fits_all)) {
  fits_region <- fits_all[[eco]]
  
  for (pred in names(fits_region)) {
    df <- fits_region[[pred]]
    
    # mask non-significant
    df$R2sign[df$p_adj > threshold_p] <- NA
    
    # value used for plotting
    df$value <- df$R2sign
    
    # best lag (max absolute R2sign)
    abs_ind <- abs(df$value)
    idx <- which(abs_ind == max(abs_ind, na.rm = TRUE))
    best <- df[idx, ]
    
    # add identifiers
    df$ecoreg <- eco
    df$predictor <- pred
    
    best$ecoreg <- eco
    best$predictor <- pred
    
    ccm_data_list[[paste(eco, pred)]] <- df
    ccm_best_list[[paste(eco, pred)]] <- best
  }
}

# Combine
ccm_data <- bind_rows(ccm_data_list)
ccm_best <- bind_rows(ccm_best_list)
ccm_data$predictor <- factor(ccm_data$predictor, levels = predictors)
ccm_best$predictor <- factor(ccm_best$predictor, levels = predictors)
ccm_data$ecoreg <- factor(ccm_data$ecoreg, levels = c("Global", "Alpine", "Continental", "Mediterranean"))
ccm_best$ecoreg <- factor(ccm_best$ecoreg, levels = c("Global", "Alpine", "Continental", "Mediterranean"))

# Set global color limits
r2_limits <- c(min(ccm_data$value, na.rm = TRUE),
               max(ccm_data$value, na.rm = TRUE))

# CCM plot based on plotCCM function from ecoxCorr
p_all <- ggplot(ccm_data, aes(lag_start, lag_end, fill = value)) +
  geom_tile() +
  geom_tile(data = ccm_best,
            color = "deeppink3",
            linewidth = 1,
            fill = NA,
            inherit.aes = FALSE,
            aes(x = lag_start, y = lag_end)) +
  geom_text(data = ccm_best,
            aes(x = 1.6, y = max_lag - 1.8,
                label = paste0("r(", lag_start, ",", lag_end, ") = ",
                               formatC(value, digits = 3, format = "f"))),
            inherit.aes = FALSE,
            fontface = "bold",
            hjust = 1,
            size = 3.2) +
  scale_x_reverse(breaks = breaks_extended(4)) +
  scale_y_continuous(breaks = breaks_extended(4)) +
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red",
                       midpoint = 0,
                       limits = r2_limits,
                       name = "signed R²",
                       na.value = "grey") +
  facet_grid(rows = vars(ecoreg),
             cols = vars(predictor)) +
  labs(title = "CCMs across ecoregions",
       x = "Lag start (weeks)", y = "Lag end (weeks)") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey95", color = "grey40"),
        strip.text = element_text(face = "bold", size = 10),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 14))
# Print
p_all

# Save
ggsave(filename = file.path(out_png, ("Allecoreg_CCM.png")),
       plot = p_all, width = 18, height = 12, dpi = 300)
