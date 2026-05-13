library(tidyverse)
library(lubridate)
library(scales)
library(patchwork)
library(cowplot)
library(grid)
library(ggh4x)

## SET UP ## 
# 1) Set dataset name
datasetName <- "VectAbundance_2024"

# 2) Set output directory
outdir_vis <- file.path(datasetName, "data_vis", "bioclim")
dir.create(outdir_vis, recursive = TRUE, showWarnings = FALSE)

# 3) Load bio data
data_bio <- read_csv(file.path(datasetName, "VectAbundance_2024_noNA.csv")) |>
  transmute(ID, lat, long,
            date  = as.Date(endDate),
            count = as.numeric(count),
            ecoreg = as.character(code)) |>
  mutate(week_date = floor_date(date, "week", week_start = 1)) |>
  filter(!is.na(ID), !is.na(date), !is.na(week_date), !is.na(count), !is.na(ecoreg))
# Load POT
data_bio_POT <- read_csv(file.path(datasetName, "VectAbundance_2024_POT.csv")) |>
  transmute(ID, lat, long,
            date  = as.Date(endDate),
            count = as.numeric(count),
            ecoreg = as.character(code)) |>
  mutate(week_date = floor_date(date, "week", week_start = 1)) |>
  filter(!is.na(ID), !is.na(date), !is.na(week_date), !is.na(count), !is.na(ecoreg))

# 4) Load clim data
data_clim <- read_csv(file.path(datasetName, "clim_daily_noNA", "VectAbundance_2024_allclim_noNA.csv")) |>
  filter(!is.na(ID), !is.na(date), !is.na(week_date), !is.na(ecoreg),
         !is.na(lat), !is.na(long))
  
# 5) Settings for clim
clim_vars <- c("mean_temp", "max_temp", "min_temp", "sum_precip", "max_windspeed")
eco_order <- c("Alpine", "Continental", "Mediterranean")
ecoreg_cols <- c("Alpine" = "firebrick", "Continental" = "forestgreen", "Mediterranean" = "steelblue", "Climate variable" = "orange")

# 6) Create keys for matching clim with POT timeseries
# exact event keys for matching climate to bio
bio_keys <- data_bio |>
  distinct(ID, lat, long, date, ecoreg)

bio_keys_POT <- data_bio_POT |>
  distinct(ID, lat, long, date, ecoreg)

# climate matched exactly to each bio version
data_clim_match <- data_clim |>
  inner_join(bio_keys, by = c("ID", "lat", "long", "date", "ecoreg"))

data_clim_match_POT <- data_clim |>
  inner_join(bio_keys_POT, by = c("ID", "lat", "long", "date", "ecoreg"))

## PLOT 1: WEEKLY ABUNDANCE AND CLIMATE BY ECOREGION ## 
# 1) Get weekly mean per ID, start with sum
bio_id_week <- data_bio |>
  group_by(ecoreg, ID, week_date) |>
  summarise(id_week_count = sum(count, na.rm = TRUE), .groups = "drop")
# Ecoreg-week mean (and sd if you want later)
bio_week <- bio_id_week |>
  group_by(ecoreg, week_date) |>
  summarise(mean_count = mean(id_week_count, na.rm = TRUE),
            sd_count = sd(id_week_count, na.rm = TRUE),
            .groups = "drop") |>
  arrange(ecoreg, week_date)

# Do the same for POT
# Get weekly mean per ID, start with sum
bio_id_week_POT <- data_bio_POT |>
  group_by(ecoreg, ID, week_date) |>
  summarise(id_week_count = sum(count, na.rm = TRUE), .groups = "drop")
# Ecoreg-week mean (and sd if you want later)
bio_week_POT <- bio_id_week_POT |>
  group_by(ecoreg, week_date) |>
  summarise(mean_count = mean(id_week_count, na.rm = TRUE),
            sd_count = sd(id_week_count, na.rm = TRUE),
            .groups = "drop") |>
  arrange(ecoreg, week_date)

# 2) Get the weekly aggregate of each clim variable
clim_week <- data_clim |>
  group_by(ecoreg, week_date) |>
  summarise(mean_temp     = mean(mean_temp, na.rm = TRUE),
            max_temp      = mean(max_temp, na.rm = TRUE),
            min_temp      = mean(min_temp, na.rm = TRUE),
            sum_precip    = sum(sum_precip, na.rm = TRUE),
            max_windspeed = mean(max_windspeed, na.rm = TRUE),
            .groups = "drop")

clim_long <- clim_week |>
  pivot_longer(all_of(clim_vars), names_to = "clim_var", values_to = "clim_value") |>
  mutate(clim_var = factor(clim_var, levels = clim_vars))

# 3) Join clim and bio (keep climate-only weeks; bio will be filtered per-layer) #
bioclim <- full_join(bio_week, clim_long, by = c("ecoreg", "week_date")) |>
  mutate(clim_var = factor(clim_var, levels = clim_vars),
         ecoreg   = factor(ecoreg, levels = eco_order)) |>
  filter(!is.na(clim_var), !is.na(ecoreg))
bioclim_POT <- full_join(bio_week_POT, clim_long, by = c("ecoreg", "week_date")) |>
  mutate(clim_var = factor(clim_var, levels = clim_vars),
         ecoreg   = factor(ecoreg, levels = eco_order)) |>
  filter(!is.na(clim_var), !is.na(ecoreg))

# 4) Helper: one panel = one ecoreg x one clim_var (TRUE per-panel sec axis)
split_lines <- as.Date(c("2013-01-01", "2018-01-01"))
plot_bioclim_facetgrid <- function(df) {
  
  df <- df |>
    filter(!is.na(clim_var), !is.na(ecoreg)) |>
    arrange(ecoreg, clim_var, week_date)
  
  # per-panel scaling (ecoreg x clim_var) so climate overlays trap counts nicely
  df_sc <- df |>
    group_by(ecoreg, clim_var) |>
    mutate(max_bio  = suppressWarnings(max(mean_count, na.rm = TRUE)),
           max_clim = suppressWarnings(max(clim_value, na.rm = TRUE)),
           scale_factor = if_else(is.finite(max_bio) & max_bio > 0 & is.finite(max_clim) & max_clim != 0,
                                  max_bio / max_clim,
                                  NA_real_),
           clim_scaled = if_else(is.finite(scale_factor), clim_value * scale_factor, NA_real_)) |>
    ungroup()
  
  xMin <- min(df_sc$week_date, na.rm = TRUE)
  xMax <- max(df_sc$week_date, na.rm = TRUE)
  
  ggplot(df_sc, aes(x = week_date)) +
    geom_vline(xintercept = split_lines, linetype = "dashed",
               color = "black", linewidth = 0.4) +
    geom_line(data = df_sc |> filter(!is.na(mean_count)), # abundance
              aes(y = mean_count, color = "Weekly mean trap counts", group = 1),
              linewidth = 0.50) +
    geom_line(data = df_sc |> filter(!is.na(clim_scaled)), # climate (scaled)
              aes(y = clim_scaled, color = "Weekly mean climate variable", group = 1),
              linewidth = 0.50) +
    facet_grid(rows = vars(clim_var), cols = vars(ecoreg), scales = "free_y") +
    scale_color_manual(values = c("Weekly mean trap counts" = "firebrick",
                                  "Weekly mean climate variable" = "steelblue"),
                       name = NULL) +
    scale_x_date(limits = c(xMin, xMax),
                 breaks = seq(floor_date(xMin, "year"), ceiling_date(xMax, "year"), by = "1 year"),
                 labels = date_format("%Y"),
                 expand = c(0, 0)) +
    scale_y_continuous(name = "Mean trap counts (weekly)",
                       sec.axis = sec_axis(~ ., name = "Climate (scaled to trap counts)")) +
    labs(title = "Relative abundance over climate, by ecoregion", # change accordingly
      x = "Year") +
    theme_minimal(base_size = 9) +
    theme(
      strip.text = element_text(size = 8, face = "bold"),
      strip.background = element_rect(fill = "grey88", color = "black"),
      strip.placement = "outside",
      strip.text.y.right = element_text(size = 8, angle = 270, face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.spacing = unit(0.25, "lines"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y.left = element_text(angle = 90),
      axis.text.y.right = element_text(angle = 270),
      axis.title = element_text(size = 9, face = "bold"),
      plot.title = element_text(size = 11, face = "bold"),
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      plot.margin = margin(5, 5, 5, 5)
    )
}

# 5) Plot and save
p1 <- plot_bioclim_facetgrid(bioclim)
p1
ggsave(filename = file.path(outdir_vis, "bioclim_full.png"), plot = p1, width = 12, height = 7, dpi = 300)
p1_POT <- plot_bioclim_facetgrid(bioclim_POT)
p1_POT
ggsave(filename = file.path(outdir_vis, "bioclim_full_POT.png"), plot = p1_POT, width = 12, height = 7, dpi = 300)

## PLOT 2: Individual plots per region ##
# 1) Plot helper
plot_one_ecoreg <- function(df, eco) {
  
  df_eco <- df |>
    filter(ecoreg == eco) |>
    arrange(clim_var, week_date)
  
  # scale climate onto trap-count scale within this ecoregion, per clim_var
  df_sc <- df_eco |>
    group_by(clim_var) |>
    mutate(max_bio  = suppressWarnings(max(mean_count, na.rm = TRUE)),
           max_clim = suppressWarnings(max(clim_value, na.rm = TRUE)),
           scale_factor = if_else(is.finite(max_bio) & max_bio > 0 & is.finite(max_clim) & max_clim != 0,
                                  max_bio / max_clim,
                                  NA_real_),
           clim_scaled = if_else(is.finite(scale_factor), clim_value * scale_factor, NA_real_)) |>
    ungroup()
  
  xMin <- min(df_sc$week_date, na.rm = TRUE)
  xMax <- max(df_sc$week_date, na.rm = TRUE)
  
  ggplot(df_sc, aes(x = week_date)) +
    geom_vline(xintercept = split_lines,
               linetype = "dashed", color = "black", linewidth = 0.4) +
    geom_line(data = df_sc |> filter(!is.na(mean_count)),
              aes(y = mean_count, color = "Weekly mean trap counts", group = 1),
              linewidth = 0.50) +
    geom_line(data = df_sc |> filter(!is.na(clim_scaled)),
              aes(y = clim_scaled, color = "Weekly mean climate variable", group = 1),
              linewidth = 0.50) +
    facet_wrap(~ clim_var, ncol = 1, scales = "free_y", strip.position = "right") +
    scale_color_manual(values = c("Weekly mean trap counts" = "firebrick",
                                  "Weekly mean climate variable" = "steelblue"),
                       name = NULL) +
    scale_x_date(limits = c(xMin, xMax),
                 breaks = seq(floor_date(xMin, "year"), ceiling_date(xMax, "year"), by = "1 year"),
                 labels = date_format("%Y"),
                 expand = c(0, 0)) +
    scale_y_continuous(name = "Mean trap counts (weekly)",
                       sec.axis = sec_axis(~ ., name = "Climate (scaled to trap counts)")) +
    labs(title = paste0("Relative abundance and climate timeseries for ", eco, " ecoregion"), 
         x = "Year") + # change accordingly
    theme_minimal(base_size = 10) +
    theme(strip.text = element_text(face = "bold"),
          strip.background = element_rect(fill = "grey88", color = "black"),
          strip.placement = "outside",
          strip.text.y.right = element_text(angle = 270, face = "bold"),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(0.4, "lines"),
          #panel.spacing = unit(0.25, "lines"),
          axis.text = element_text(size = 7),
          axis.text.x = element_text(hjust = 1),
          axis.text.y.left = element_text(angle = 90),
          axis.text.y.right = element_text(angle = 270),
          axis.title = element_text(face = "bold"),
          plot.title = element_text(size = 11, face = "bold"),
          legend.position = "bottom",
          plot.margin = margin(5, 5, 5, 5) 
    )
}

# 2) Plot and save
p2_alp <- plot_one_ecoreg(bioclim, "Alpine")
p2_alp
ggsave(filename = file.path(outdir_vis, "full_Alpine.png"), p2_alp, width = 12, height = 7, dpi = 300)
p2_con <- plot_one_ecoreg(bioclim, "Continental")
p2_con
ggsave(filename = file.path(outdir_vis, "full_Continental.png"), p2_con, width = 12, height = 7, dpi = 300)
p2_med <- plot_one_ecoreg(bioclim, "Mediterranean")
p2_med
ggsave(filename = file.path(outdir_vis, "full_Mediterranean.png"), p2_med, width = 12, height = 7, dpi = 300)
p2_alp_POT <- plot_one_ecoreg(bioclim_POT, "Alpine")
p2_alp_POT
ggsave(filename = file.path(outdir_vis, "full_Alpine_POT.png"), p2_alp_POT, width = 12, height = 7, dpi = 300)
p2_con_POT <- plot_one_ecoreg(bioclim_POT, "Continental")
p2_con_POT
ggsave(filename = file.path(outdir_vis, "full_Continental_POT.png"), p2_con_POT, width = 12, height = 7, dpi = 300)
p2_med_POT <- plot_one_ecoreg(bioclim_POT, "Mediterranean")
p2_med_POT
ggsave(filename = file.path(outdir_vis, "full_Mediterranean_POT.png"), p2_med_POT, width = 12, height = 7, dpi = 300)

## PLOT 3: Yearly mean ##
# 1) Get yearly abundance from weekly abundance
bio_year <- bio_week |>
  mutate(year_date = as.Date(paste0(year(week_date), "-01-01"))) |>
  group_by(ecoreg, year_date) |>
  summarise(mean_count = mean(mean_count, na.rm = TRUE),
            sd_count   = sd(mean_count, na.rm = TRUE),
            .groups = "drop") |>
  mutate(sd_count = replace_na(sd_count, 0)) |>
  arrange(ecoreg, year_date)

# Do the same for POT
bio_year_POT <- bio_week_POT |>
  mutate(year_date = as.Date(paste0(year(week_date), "-01-01"))) |>
  group_by(ecoreg, year_date) |>
  summarise(mean_count = mean(mean_count, na.rm = TRUE),
            sd_count   = sd(mean_count, na.rm = TRUE),
            .groups = "drop") |>
  mutate(sd_count = replace_na(sd_count, 0)) |>
  arrange(ecoreg, year_date)

# 2) Get yearly climate from weekly climate
clim_year <- clim_week |>
  mutate(year_date = as.Date(paste0(year(week_date), "-01-01"))) |>
  group_by(ecoreg, year_date) |>
  summarise(mean_temp     = mean(mean_temp, na.rm = TRUE),
            max_temp      = mean(max_temp, na.rm = TRUE),
            min_temp      = mean(min_temp, na.rm = TRUE),
            sum_precip    = sum(sum_precip, na.rm = TRUE),
            max_windspeed = mean(max_windspeed, na.rm = TRUE),
            .groups = "drop")

clim_long_year <- clim_year |>
  pivot_longer(all_of(clim_vars), names_to = "clim_var", values_to = "clim_value") |>
  mutate(clim_var = factor(clim_var, levels = clim_vars))

# 3) Join yearly abundance and yearly climate
bioclim_year <- full_join(bio_year, clim_long_year, by = c("ecoreg", "year_date")) |>
  mutate(clim_var = factor(clim_var, levels = clim_vars),
         ecoreg   = factor(ecoreg, levels = eco_order)) |>
  filter(!is.na(clim_var), !is.na(ecoreg))

bioclim_year_POT <- full_join(bio_year_POT, clim_long_year, by = c("ecoreg", "year_date")) |>
  mutate(clim_var = factor(clim_var, levels = clim_vars),
         ecoreg   = factor(ecoreg, levels = eco_order)) |>
  filter(!is.na(clim_var), !is.na(ecoreg))

# 4) Plot helper
plot_bioclim_facetgrid_year <- function(df) {
  
  df <- df |>
    filter(!is.na(clim_var), !is.na(ecoreg)) |>
    arrange(ecoreg, clim_var, year_date)
  
  # per-panel scaling
  df_sc <- df |>
    group_by(ecoreg, clim_var) |>
    mutate(
      max_bio  = suppressWarnings(max(mean_count, na.rm = TRUE)),
      max_clim = suppressWarnings(max(clim_value, na.rm = TRUE)),
      scale_factor = if_else(
        is.finite(max_bio) & max_bio > 0 & is.finite(max_clim) & max_clim != 0,
        max_bio / max_clim,
        NA_real_
      ),
      clim_scaled = if_else(is.finite(scale_factor), clim_value * scale_factor, NA_real_)
    ) |>
    ungroup()
  
  xMin <- min(df_sc$year_date, na.rm = TRUE)
  xMax <- max(df_sc$year_date, na.rm = TRUE)
  
  ggplot(df_sc, aes(x = year_date)) +
    geom_vline(
      xintercept = split_lines,
      linetype = "dashed",
      color = "black",
      linewidth = 0.4
    ) +
    geom_line(
      data = df_sc |> filter(!is.na(mean_count)),
      aes(y = mean_count, color = "Annual mean trap counts", group = 1),
      linewidth = 0.55
    ) +
    geom_line(
      data = df_sc |> filter(!is.na(clim_scaled)),
      aes(y = clim_scaled, color = "Annual mean climate variable", group = 1),
      linewidth = 0.55
    ) +
    facet_grid(rows = vars(clim_var), cols = vars(ecoreg), scales = "free_y") +
    scale_color_manual(
      values = c(
        "Annual mean trap counts" = "firebrick",
        "Annual mean climate variable" = "steelblue"
      ),
      name = NULL
    ) +
    scale_x_date(
      limits = c(xMin, xMax),
      breaks = seq(xMin, xMax, by = "1 year"),
      labels = date_format("%Y"),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      name = "Mean trap counts (annual)",
      sec.axis = sec_axis(~ ., name = "Climate (scaled to trap counts)")
    ) +
    labs(title = "Annual relative abundance above threshold (n ≥ 55) over climate, by ecoregion", # change accordingly
         x = "Year") +
    theme_minimal(base_size = 9) +
    theme(strip.text = element_text(size = 8, face = "bold"),
          strip.background = element_rect(fill = "grey88", color = "black"),
          strip.placement = "outside",
          strip.text.y.right = element_text(size = 8, angle = 270, face = "bold"),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          panel.spacing = unit(0.25, "lines"),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y.left = element_text(angle = 90),
          axis.text.y.right = element_text(angle = 270),
          axis.title = element_text(size = 9, face = "bold"),
          plot.title = element_text(size = 11, face = "bold"),
          legend.position = "bottom",
          legend.text = element_text(size = 8),
          plot.margin = margin(5, 5, 5, 5))
}

# 5) Plot and save
p3 <- plot_bioclim_facetgrid_year(bioclim_year)
p3
ggsave(filename = file.path(outdir_vis, "bioclim_yearly.png"), p3, width = 12, height = 7, dpi = 300)
p3_POT <- plot_bioclim_facetgrid_year(bioclim_year_POT)
p3_POT
ggsave(filename = file.path(outdir_vis, "bioclim_yearly_POT.png"), p3_POT, width = 12, height = 7, dpi = 300)

## PLOT 4: YEARLY MEAN BY TIMESERIES ##
# 1) Build yearly climate split tables directly from bioclim_year objects
clim_long_year_split <- bioclim_year |>
  select(ecoreg, year_date, clim_var, clim_value) |>
  distinct() |>
  mutate(year = year(year_date),
         period = case_when(year >= 2008 & year <= 2012 ~ "2008–2012",
                            year >= 2013 & year <= 2017 ~ "2013–2017",
                            year >= 2018 & year <= 2022 ~ "2018–2022",
                            TRUE ~ NA_character_)) |> 
  filter(!is.na(period))

clim_long_year_split_POT <- bioclim_year_POT |>
  select(ecoreg, year_date, clim_var, clim_value) |>
  distinct() |>
  mutate(year = year(year_date),
         period = case_when(year >= 2008 & year <= 2012 ~ "2008–2012",
                            year >= 2013 & year <= 2017 ~ "2013–2017",
                            year >= 2018 & year <= 2022 ~ "2018–2022",
                            TRUE ~ NA_character_)) |> 
  filter(!is.na(period))

# 2) Abundance summary per split
bio_split_mean <- bio_year |>
  mutate(year = year(year_date),
         period = case_when(year >= 2008 & year <= 2012 ~ "2008–2012",
                            year >= 2013 & year <= 2017 ~ "2013–2017",
                            year >= 2018 & year <= 2022 ~ "2018–2022",
                            TRUE ~ NA_character_)) |> 
  filter(!is.na(period)) |>
  group_by(ecoreg, period) |>
  summarise(mean_count = mean(mean_count, na.rm = TRUE),
            .groups = "drop") |> 
  mutate(period = factor(period, levels = c("2008–2012", "2013–2017", "2018–2022")))

bio_split_mean_POT <- bio_year_POT |>
  mutate(year = year(year_date),
         period = case_when(year >= 2008 & year <= 2012 ~ "2008–2012",
                            year >= 2013 & year <= 2017 ~ "2013–2017",
                            year >= 2018 & year <= 2022 ~ "2018–2022",
                            TRUE ~ NA_character_)) |> 
  filter(!is.na(period)) |>
  group_by(ecoreg, period) |>
  summarise(mean_count = mean(mean_count, na.rm = TRUE),
            .groups = "drop") |> 
  mutate(period = factor(period, levels = c("2008–2012", "2013–2017", "2018–2022")))

# 3) Climate summary per split
# one line per climate variable = mean across ecoregions within each split
clim_split_mean <- clim_long_year_split |>
  group_by(clim_var, period) |>
  summarise(mean_clim = mean(clim_value, na.rm = TRUE),
            .groups = "drop") |> 
  mutate(period = factor(period, levels = c("2008–2012", "2013–2017", "2018–2022")),
         clim_var = factor(clim_var, levels = clim_vars))
clim_split_mean_POT <- clim_long_year_split_POT |>
  group_by(clim_var, period) |>
  summarise(mean_clim = mean(clim_value, na.rm = TRUE),
            .groups = "drop") |> 
  mutate(period = factor(period, levels = c("2008–2012", "2013–2017", "2018–2022")),
         clim_var = factor(clim_var, levels = clim_vars))

# 4) Plot helper
plot_split_means_bioclim <- function(df_bio, df_clim) {
  
  # repeat abundance data across climate facets
  bio_plot <- crossing(df_bio, clim_var = factor(clim_vars, levels = clim_vars))
  
  # scale climate to abundance range within each climate-variable panel
  scale_tbl <- bio_plot |>
    group_by(clim_var) |>
    summarise(max_bio = max(mean_count, na.rm = TRUE), .groups = "drop") |>
    left_join(df_clim |>
                group_by(clim_var) |>
                summarise(max_clim = max(mean_clim, na.rm = TRUE), .groups = "drop"),
              by = "clim_var") |>
    mutate(scale_factor = if_else(is.finite(max_bio) & max_bio > 0 & is.finite(max_clim) & max_clim != 0,
                                  max_bio / max_clim,
                                  NA_real_))
  
  clim_plot <- df_clim |>
    left_join(scale_tbl, by = "clim_var") |>
    mutate(clim_scaled = if_else(is.finite(scale_factor), mean_clim * scale_factor, NA_real_))
  
  ggplot() +
    geom_line(data = bio_plot,
              aes(x = period, y = mean_count, color = ecoreg, group = ecoreg),
              linewidth = 0.9) +
    geom_point(data = bio_plot,
               aes(x = period, y = mean_count, color = ecoreg, group = ecoreg),
               size = 2) +
    geom_line(data = clim_plot,
              aes(x = period, y = clim_scaled, color = "Climate variable", group = 1),
              linewidth = 0.9) +
    geom_point(data = clim_plot,
               aes(x = period, y = clim_scaled, color = "Climate variable", group = 1),
               size = 2) +
    facet_wrap(~ clim_var, nrow = 1, scales = "free_y") +
    scale_color_manual(values = ecoreg_cols, name = NULL) +
    scale_y_continuous(name = "Mean trap counts",
                       sec.axis = sec_axis(~ ., name = "Climate (scaled to trap counts)")) +
    labs(title = "Mean abundance over threshold (n ≥ 55) and climate across timeseries splits", # change accordingly
      x = "Timeseries split") +
    theme_minimal(base_size = 10) +
    theme(strip.text = element_text(face = "bold"),
          strip.background = element_rect(fill = "grey88", color = "black"),
          strip.placement = "outside",
          strip.text.y.right = element_text(angle = 270, face = "bold"),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.spacing = unit(0.35, "lines"),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(color = "black"),
          axis.ticks = element_line(linewidth = 0.25, color = "black"),
          axis.ticks.length = unit(3, "pt"),
          axis.ticks.x.top = element_blank(),
          axis.ticks.y.right = element_blank(),
          plot.title = element_text(face = "bold"),
          legend.position = "bottom")
}

# 5) Plot and save
p4 <- plot_split_means_bioclim(bio_split_mean, clim_split_mean)
p4
ggsave(filename = file.path(outdir_vis, "bioclim_split_means.png"),
       plot = p4, width = 14, height = 5, dpi = 300)

p4_POT <- plot_split_means_bioclim(bio_split_mean_POT, clim_split_mean_POT)
p4_POT
ggsave(filename = file.path(outdir_vis, "bioclim_split_means_POT.png"),
       plot = p4_POT, width = 14, height = 5, dpi = 300)