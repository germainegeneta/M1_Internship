#devtools::install_github("Nmoiroux/ecoXCorr")
library(ecoXCorr)
library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)
library(viridis)
library(purrr)
library(tidyverse)

## SET UP ## 
# Set dataset name
datasetName <- "VectAbundance_2024"

# Set output directories
outdir_ccm <- file.path(datasetName, "CCM_outputs")
out_csv <- file.path(outdir_ccm, "fitmodels")
out_png <- file.path(outdir_ccm, "maps")
dir.create(out_csv, recursive = TRUE, showWarnings = FALSE)
dir.create(out_png, recursive = TRUE, showWarnings = FALSE)

# Load biological data
data_aedes <- read_csv(file.path(datasetName, "VectAbundance_2024_POT.csv")) |> 
  transmute(
    ID,
    date   = as.Date(endDate),
    count  = as.numeric(count),
    count_c = as.integer(ceiling(as.numeric(count))),
    lat, long,
    ecoreg = code) |>
  filter(!is.na(count_c), count_c > 0)

# Load climate data
data_clim <- read_csv(file.path(datasetName, "clim_daily_nonzero", "VectAbundance_2024_nonzero_allclim.csv")) 

# Set parameters
ecoregs   <- c("Alpine", "Continental", "Mediterranean")
clim_vars <- c("mean_temp", "max_temp", "min_temp", "sum_precip", "max_windspeed")
interval <- 7
max_lag  <- 8
response <- "count_c"
random   <- "(1|ID)"
family   <- "truncated_nbinom2"
threshold_p <- 0.2 

## 1: CREATE CCM PLOT FUNCTION ##
# Helper function to make CCMs
make_ccm_plot <- function(data,
                          mode = c("local", "global"),
                          r2_limits = NULL,
                          threshold_p = 0.2,
                          max_lag = 8,
                          show_title = TRUE,
                          title_text = NULL,
                          show_axes = TRUE) {
  mode <- match.arg(mode)
  indicator <- "R2sign"
  
  data[data$p_adj > threshold_p, indicator] <- NA
  
  df <- data |>
    filter(lag_start >= lag_end,
           lag_start >= 1, lag_start <= max_lag,
           lag_end   >= 1, lag_end   <= max_lag) |>
    mutate(value = .data[[indicator]])
  
  i_best <- which.max(replace(df$value, is.na(df$value), -Inf))
  best_s <- df$lag_start[i_best]
  best_e <- df$lag_end[i_best]
  best_r <- df$value[i_best]
  
  best_label <- paste0(
    "r(", best_s, ",", best_e, ") = ",
    formatC(best_r, format = "f", digits = 3)
  )
  
  p <- ggplot(df, aes(lag_start, lag_end, fill = value)) +
    geom_tile() +
    annotate("rect",
             xmin = best_s - 0.5, xmax = best_s + 0.5,
             ymin = best_e - 0.5, ymax = best_e + 0.5,
             fill = NA, color = "black", linewidth = 1.1) +
    annotate("text",
             x = 1.6, y = max_lag - 1.8, 
             label = best_label, size = 3.6, 
             hjust = 1, fontface = "bold") +
    scale_x_reverse(breaks = 1:max_lag) +
    scale_y_continuous(breaks = 1:max_lag) +
    labs(title = if (show_title) title_text else NULL,
         x = if (show_axes) "Lag start (weeks)" else NULL, 
         y = if (show_axes) "Lag end (weeks)" else NULL) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 11),
          plot.margin = margin(4, 4, 4, 4))
  
  if (mode == "local") {
    p <- p + scale_fill_viridis_c(option = "magma", name = expression("Signed " * R^2), na.value = "grey80")
  } else {
    p <- p + scale_fill_viridis_c(option = "magma", limits = r2_limits,
                                  oob = squish, name = expression("Signed " * R^2), na.value = "grey80")
  }
  
  if (!show_axes) {p <- p + theme(axis.title.x = element_blank(),
                                  axis.title.y = element_blank())
  }
  p <- p + theme(plot.margin = margin(1, 1, 1, 1),
                 axis.text.x = element_text(size = 7),
                 axis.text.y = element_text(size = 7))
  p
}

## 2: RUN MODELS FOR ALL TIME LAGS ##
fits <- list()
for (eco in ecoregs) {
  aedes_eco <- data_aedes |> filter(ecoreg == eco)
  clim_eco  <- data_clim  |> filter(ecoreg == eco)
  ref_date  <- sort(unique(aedes_eco$date))
  
  for (v in clim_vars) {
    clim_agg <- aggregate_lagged_intervals(data       = clim_eco,
                                           date_col   = "date",
                                           value_cols = v,
                                           ref_date   = ref_date,
                                           interval   = interval,
                                           max_lag    = max_lag,
                                           funs       = list(mean = mean),
                                           na.rm      = TRUE)
    merged_eco <- merge(clim_agg, aedes_eco, by = "date", all = TRUE)
    pred <- paste0(v, "_mean")
    if (!pred %in% names(merged_eco)) next
    fit_tbl <- fit_models_by_lag(data       = merged_eco,
                                 response   = response,
                                 predictors = pred,
                                 random     = random,
                                 family     = family)
    write_csv(fit_tbl, file.path(out_csv, paste0(eco, "_fitmodels_", v, "_POT.csv")))
    fits[[paste0(eco, "__", v)]] <- fit_tbl
  }
}

## 3: CREATE CCMS ## 
for (eco in ecoregs) {
  for (v in clim_vars) {
    
    key <- paste0(eco, "__", v)
    if (!key %in% names(fits)) next
    fit_tbl <- fits[[key]]
    
    p_ind <- make_ccm_plot(data = fit_tbl,
                           mode = "local",               # per-plot scale
                           threshold_p = threshold_p,
                           max_lag = max_lag,
                           show_title = TRUE,
                           title_text = paste0(eco, " (full timeseries) - ", v),
                           show_axes = TRUE)
    
    ggsave(filename = file.path(out_png, paste0(eco, "_", v, "_full_POT.png")),
           plot = p_ind, width = 8, height = 6, dpi = 300)
  }
}

## 4: COMPILE INDIVIDUAL PLOTS INTO ONE PNG ##
# build one combined DF
threshold_p <- 0.2
fits_long <- imap_dfr(fits, function(df, key) {
  parts <- strsplit(key, "__", fixed = TRUE)[[1]]
  eco <- parts[1]
  var <- parts[2]
  
  df |>
    transmute(ecoreg   = eco, clim_var = var,
              lag_start, lag_end,
              R2    = R2sign, p_adj = p_adj)
})

# mask by p_adj like plotCCM + keep Curriero triangle
fits_long <- fits_long |>
  mutate(R2_plot = ifelse(p_adj > threshold_p, NA_real_, R2)) |>
  filter(lag_start >= lag_end,
         lag_start >= 1, lag_start <= max_lag,
         lag_end   >= 1, lag_end   <= max_lag)

# set factors (force facet order, not alphabetical)
clim_order <- c("mean_temp", "max_temp", "min_temp", "sum_precip", "max_windspeed")
eco_order  <- c("Alpine", "Continental", "Mediterranean")

fits_long <- fits_long |>
  mutate(clim_var = factor(clim_var, levels = clim_order),
         ecoreg   = factor(ecoreg,   levels = eco_order))

# compute best lag per facet
best_pts <- fits_long |>
  group_by(ecoreg, clim_var) |>
  filter(!is.na(R2_plot)) |>
  slice_max(R2_plot, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(best_label = paste0("r(", lag_start, ",", lag_end, ") = ",
                             formatC(R2_plot, digits = 3, format = "f")))

# global limits for compiled plot
r2_limits <- c(min(fits_long$R2), max(fits_long$R2))
p_facet <- ggplot(fits_long, aes(lag_start, lag_end, fill = R2_plot)) +
  geom_tile() +
  geom_rect(data = best_pts, # best-lag light-blue box per facet
            aes(xmin = lag_start - 0.5, xmax = lag_start + 0.5,
                ymin = lag_end - 0.5,   ymax = lag_end + 0.5),
            inherit.aes = FALSE,
            fill = NA, color = "red", linewidth = 1.2) +
  geom_text(data = best_pts, # best-lag label per facet
            aes(x = 2.3, y = max_lag - 1.3, label = best_label),
            inherit.aes = FALSE,
            fontface = "bold",
            size = 3.2,
            hjust = 1) +
  scale_x_reverse(breaks = 1:max_lag) +
  scale_y_continuous(breaks = 1:max_lag) +
  scale_fill_viridis_c(option = "magma",
                       limits = r2_limits,
                       oob = squish,
                       name = expression("Signed " * R^2),
                       na.value = "grey80") +
  facet_grid(rows = vars(ecoreg), cols = vars(clim_var)) +
  labs(title = "Above-threshold-CCMs (n ≥ 55) across ecoregions (full timeseries)",
       x = "Lag start (weeks)", y = "Lag end (weeks)") +
  theme_bw() +
  theme(legend.position = "right",
        strip.background = element_rect(fill = "grey90", colour = "grey40", linewidth = 0.6),
        strip.text = element_text(face = "bold"),
        panel.spacing = unit(0.6, "lines"),
        plot.title = element_text(face = "bold"))

p_facet
ggsave(file.path(outdir_ccm, "Allecoreg_full_POT.png"),
       p_facet, width = 18, height = 10, dpi = 300)


## 4: FIT SPLIT TIMESERIES MODELS + SAVE INDIVIDUAL PNGs ##
# define splits
splits <- tibble(
  split_id = c("2008-2012", "2013-2017", "2018-2022"),
  start    = as.Date(c("2008-01-01", "2013-01-01", "2018-01-01")),
  end      = as.Date(c("2012-12-31", "2017-12-31", "2022-12-31"))
)

# object to store fit tables
fits_split <- list()
for (s in seq_len(nrow(splits))) {
  split_id <- splits$split_id[s]
  d1 <- splits$start[s]
  d2 <- splits$end[s]
  
  message("\n=== Running split: ", split_id, " ===")
  
  aedes_s <- data_aedes |> filter(date >= d1, date <= d2)
  clim_s <- data_clim |> filter(date >= d1, date <= d2)
  
  fits_s <- list()
  for (eco in ecoregs) {
    
    aedes_eco <- aedes_s |> filter(ecoreg == eco)
    clim_eco  <- clim_s  |> filter(ecoreg == eco)
    
    if (nrow(aedes_eco) == 0 || nrow(clim_eco) == 0) next
    
    ref_date <- sort(unique(aedes_eco$date))
    
    for (v in clim_vars) {
      clim_agg <- aggregate_lagged_intervals(data       = clim_eco,
                                             date_col   = "date",
                                             value_cols = v,
                                             ref_date   = ref_date,
                                             interval   = interval,
                                             max_lag    = max_lag,
                                             funs       = list(mean = mean),
                                             na.rm      = TRUE)
      merged_eco <- merge(clim_agg, aedes_eco, by = "date", all = TRUE)
      
      pred <- paste0(v, "_mean")
      if (!pred %in% names(merged_eco)) next
      
      message("Fitting: split = ", split_id, " | eco = ", eco, " | var = ", v)
      fit_tbl <- tryCatch(
        fit_models_by_lag(
          data       = merged_eco,
          response   = response,
          predictors = pred,
          random     = random,
          family     = family
        ),
        error = function(e) {
          message("FAILED: split = ", split_id,
                  " | eco = ", eco,
                  " | var = ", v,
                  " -> ", conditionMessage(e))
          return(NULL)
        }
      )
      
      if (is.null(fit_tbl)) next
      
      write_csv(
        fit_tbl,
        file.path(out_csv, paste0(eco, "_fitmodels_", v, "_", split_id, "_POT.csv"))
      )
      
      fits_s[[paste0(eco, "__", v)]] <- fit_tbl
    }
  }
  fits_split[[split_id]] <- fits_s
}

## 5: SAVE INDIVIDUAL PNGs FOR SPLIT MODELS ##
for (split_id in names(fits_split)) {
  fits_s <- fits_split[[split_id]]
  
  for (eco in ecoregs) {
    for (v in clim_vars) {
      key <- paste0(eco, "__", v)
      if (!key %in% names(fits_s)) next
      
      fit_tbl <- fits_s[[key]]
      p_ind <- make_ccm_plot(data        = fit_tbl,
                             mode        = "local",
                             threshold_p = threshold_p,
                             max_lag     = max_lag,
                             show_title  = TRUE,
                             title_text  = paste0(eco, " (", split_id, ") - ", v),
                             show_axes   = TRUE)
      ggsave(filename = file.path(out_png, paste0(eco, "_", v, "_", split_id, "_POT.png")),
             plot = p_ind, width = 8, height = 6, dpi = 300)
    }
  }
}

## 6: COMPILE SPLIT CCMs BY ECOREGION ##
for (eco in ecoregs) {
  fits_long_eco <- imap_dfr(fits_split, function(fits_s, split_id) {
    
    eco_keys <- names(fits_s)[grepl(paste0("^", eco, "__"), names(fits_s))]
    if (length(eco_keys) == 0) return(tibble())
    
    map_dfr(eco_keys, function(key) {
      var <- strsplit(key, "__", fixed = TRUE)[[1]][2]
      df  <- fits_s[[key]]
      
      df |> transmute(split_id = split_id,
                      ecoreg   = eco, clim_var = var,
                      lag_start, lag_end,
                      R2    = R2sign, p_adj = p_adj)
    })
  })
  
  if (nrow(fits_long_eco) == 0) {
    message("No fits found for ", eco, " — skipping compiled plot.")
    next
  }
  
  fits_long_eco <- fits_long_eco |>
    mutate(R2_plot = ifelse(p_adj > threshold_p, NA_real_, R2),
           clim_var = factor(clim_var, levels = clim_order),
           split_id = factor(split_id, levels = splits$split_id)) |>
    filter(lag_start >= lag_end,
           lag_start >= 1, lag_start <= max_lag,
           lag_end   >= 1, lag_end   <= max_lag)
  best_pts_eco <- fits_long_eco |>
    group_by(split_id, clim_var) |>
    filter(!is.na(R2_plot)) |>
    slice_max(R2_plot, n = 1, with_ties = FALSE) |>
    ungroup() |>
    mutate(best_label = paste0("r(", lag_start, ",", lag_end, ") = ",
                               formatC(R2_plot, digits = 3, format = "f")))
  r2_limits_eco <- c(min(fits_long_eco$R2, na.rm = TRUE), max(fits_long_eco$R2, na.rm = TRUE))
  p_facet_eco <- ggplot(fits_long_eco, aes(lag_start, lag_end, fill = R2_plot)) +
    geom_tile() +
    geom_rect(data = best_pts_eco,
              aes(xmin = lag_start - 0.5, xmax = lag_start + 0.5,
                  ymin = lag_end - 0.5,   ymax = lag_end + 0.5),
              inherit.aes = FALSE, fill = NA, color = "red", linewidth = 1.2) +
    geom_text(data = best_pts_eco, 
              aes(x = 2.3, y = max_lag - 1.3, label = best_label),
              inherit.aes = FALSE, fontface = "bold", size = 3.2, hjust = 1) +
    scale_x_reverse(breaks = 1:max_lag) +
    scale_y_continuous(breaks = 1:max_lag) +
    scale_fill_viridis_c(option = "magma", 
                         limits = r2_limits_eco,
                         oob = squish, 
                         name = expression("Signed " * R^2),
                         na.value = "grey80") +
    facet_grid(rows = vars(split_id), cols = vars(clim_var)) +
    labs(title = paste0("Above-threshold-CCMs (n ≥ 55) across timeseries splits in ", eco, " ecoregion"),
         x = "Lag start (weeks)", 
         y = "Lag end (weeks)") +
    theme_bw() +
    theme(legend.position = "right",
          strip.background = element_rect(fill = "grey90", colour = "grey40", linewidth = 0.6),
          strip.text = element_text(face = "bold"),
          panel.spacing = unit(0.6, "lines"),
          plot.title = element_text(face = "bold"))
  
  ggsave(filename = file.path(out_png, paste0(eco, "_allsplits_POT.png")),
         plot = p_facet_eco, width = 18, height = 10, dpi = 300)
}

## 7: PLOT THE BEST TIME LAG START AND END ##
extract_best_lags <- function(fits_s, split_id, threshold_p = 0.2, max_lag = 8) {
  
  imap_dfr(fits_s, function(df, key) {
    
    parts <- strsplit(key, "__", fixed = TRUE)[[1]]
    eco <- parts[1]
    var <- parts[2]
    
    df2 <- df |>
      mutate(R2_plot = if_else(p_adj > threshold_p, NA_real_, R2sign)) |>
      filter(lag_start >= lag_end,
             between(lag_start, 1, max_lag),
             between(lag_end,   1, max_lag))
    
    if (nrow(df2) == 0 || all(is.na(df2$R2_plot))) return(tibble())
    
    df2 |>
      filter(!is.na(R2_plot)) |>
      slice_max(R2_plot, n = 1, with_ties = FALSE) |>
      transmute(split_id  = split_id,
                ecoreg    = eco, clim_var  = var,
                lag_start = lag_start, lag_end   = lag_end)
  })
}

# extract best lags
best_lags_splits <- imap_dfr(
  fits_split,
  ~ extract_best_lags(.x, split_id = .y,
                      threshold_p = threshold_p,
                      max_lag = max_lag))

# prepare plotting dataframe
best_lags_plot <- best_lags_splits |>
  mutate(split_id = factor(split_id, levels = splits$split_id),
         clim_var = factor(clim_var, levels = clim_order),
         ecoreg   = factor(ecoreg,   levels = eco_order))

eco_cols <- c(Alpine        = "#f8766d", 
              Continental   = "#00ba39", 
              Mediterranean = "#609dff")

# convert to diverging format
best_lags_div <- best_lags_plot |>
  pivot_longer(cols = c(lag_start, lag_end),
               names_to = "metric",
               values_to = "lag_week") |>
  mutate(lag_plot = if_else(metric == "lag_start", lag_week, -lag_week))

# dodge spacing (3 bars per split)
pd <- position_dodge(width = 0.8)

# plot
p_lag_div <- ggplot(best_lags_div,
                    aes(x = split_id, y = lag_plot, fill = ecoreg)) +
  geom_col(position = pd, width = 0.22) +
  geom_hline(yintercept = 0, linewidth = 0.7, color = "black") +
  facet_wrap(~ clim_var, nrow = 1) +
  scale_fill_manual(values = eco_cols, drop = FALSE) +
  scale_y_continuous(limits = c(-max_lag, max_lag),
                     breaks = c(-max_lag:-1, 1:max_lag),
                     labels = c(max_lag:1, 1:max_lag)) +
  labs(title = "Best above-threshold-CCMs (n ≥ 55) lag start (top) and end (bottom) across timeseries splits",
       x = "Timeseries split", y = "Lag (weeks before sampling)", fill = NULL) +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "grey90", colour = "grey40", linewidth = 0.6),
        strip.text = element_text(face = "bold"),
        panel.spacing = unit(0.8, "lines"),
        plot.title = element_text(face = "bold"),
        panel.grid.major.y = element_line(color = "grey80", linewidth = 0.4),
        panel.grid.minor = element_blank())
p_lag_div

ggsave(file.path(outdir_ccm, "bestlag_diverging_POT.png"), p_lag_div,
       width = 16, height = 5, dpi = 300)


# 8: PLOT BEST LAG LENGTH  # 
best_laglen_plot <- best_lags_splits |>
  mutate(lag_len  = lag_start - lag_end + 1,
         split_id = factor(split_id, levels = splits$split_id),
         clim_var = factor(clim_var, levels = clim_order),
         ecoreg   = factor(ecoreg,   levels = eco_order))

# plot
p_laglen <- ggplot(best_laglen_plot,
                   aes(x = as.numeric(split_id),
                       y = lag_len,
                       group = ecoreg,
                       color = ecoreg)) +
  geom_line(linewidth = 1.3) +
  geom_point(size = 2.5) +
  scale_color_manual(values = eco_cols, drop = FALSE) +
  scale_x_continuous(breaks = seq_along(splits$split_id),
                     labels = splits$split_id) +
  scale_y_continuous(breaks = 1:max_lag,
                     limits = c(1, max_lag)) +
  facet_wrap(~ clim_var, nrow = 1) +
  labs(title = "Evolution of best lag length of above-threshold-CCMs (n ≥ 55) across timeseries splits",
       x = "Timeseries split", y = "Best lag length (weeks before sampling)", color = NULL) +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "grey90", colour = "grey40", linewidth = 0.6),
        strip.text = element_text(face = "bold"),
        panel.spacing = unit(1.5, "lines"),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 25, vjust = 0.7),
        plot.title = element_text(face = "bold"),
        plot.margin = margin(10, 15, 10, 15))

p_laglen
ggsave(file.path(outdir_ccm, "bestlag_length_POT.png"), p_laglen,
       width = 16, height = 5, dpi = 300)
