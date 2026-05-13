library(tidyverse)
library(lubridate)
library(scales)

## SET UP ##

# Set dataset name
datasetName <- "VectAbundance_2024"

# Create output directory
outdir_vis <- file.path(datasetName, "data_vis", "bio_only")
dir.create(outdir_vis, recursive = TRUE, showWarnings = FALSE)

# Load abundance timeseries, keep original date, and add week_date (the ISO week of that date)
data_aedes <- read_csv(file.path(datasetName, "VectAbundance_2024_noNA.csv")) |>
  transmute(ID, date  = as.Date(endDate), count, ecoreg = code) |>
  mutate(week_date = floor_date(date, "week", week_start = 1)) |>    # Monday-anchored week
  filter(!is.na(ID), !is.na(date), !is.na(week_date), !is.na(count), !is.na(ecoreg))

# Also do the same for abundance data above the threshold (55)
data_aedes_POT <- read_csv(file.path(datasetName, "VectAbundance_2024_POT.csv")) |>
  transmute(ID, date  = as.Date(endDate), count, ecoreg = code) |>
  mutate(week_date = floor_date(date, "week", week_start = 1)) |>    # Monday-anchored week
  filter(!is.na(ID), !is.na(date), !is.na(week_date), !is.na(count), !is.na(ecoreg))

## PLOT 1: WEEKLY ABUNDANCE BY ECOREGION ##
# 1) Get weekly sum per ID (sampling unit)
id_week <- data_aedes |>
  group_by(ecoreg, ID, week_date) |>
  summarise(id_week_count = sum(count, na.rm = TRUE), .groups = "drop")
# Do the same for POT
id_week_POT <- data_aedes_POT |>
  group_by(ecoreg, ID, week_date) |>
  summarise(id_week_count = sum(count, na.rm = TRUE), .groups = "drop")

# 2) Get the ecoregion weekly mean ± SD across IDs
agg_aedes <- id_week |>
  group_by(ecoreg, week_date) |>
  summarise(mean_count = mean(id_week_count, na.rm = TRUE),
            sd_count   = sd(id_week_count,   na.rm = TRUE),
            n_id       = n_distinct(ID),
            .groups    = "drop") |>
  mutate(sd_count = if_else(n_id < 2 | is.na(sd_count), 0, sd_count)) |>
  arrange(ecoreg, week_date)

# Do the same for POT
agg_aedes_POT <- id_week_POT |>
  group_by(ecoreg, week_date) |>
  summarise(mean_count = mean(id_week_count, na.rm = TRUE),
            sd_count   = sd(id_week_count,   na.rm = TRUE),
            n_id       = n_distinct(ID),
            .groups    = "drop") |>
  mutate(sd_count = if_else(n_id < 2 | is.na(sd_count), 0, sd_count)) |>
  arrange(ecoreg, week_date)

# 3) Helper function year_vlines
year_vlines <- function(xMin, xMax) {
  v <- seq(floor_date(xMin, "year"), ceiling_date(xMax, "year"), by = "1 year")
  v[v >= xMin & v <= xMax]}

# 4) Plot (facet by ecoreg; legend = line vs ribbon; strip box filled)
y_max_cap <- 900
plot_full <- function(df) {
  
  xMin <- min(df$week_date, na.rm = TRUE)
  xMax <- max(df$week_date, na.rm = TRUE)
  vY   <- year_vlines(xMin, xMax)
  
  # global maximum across ALL ecoregs (including SD ribbon)
  global_max <- y_max_cap
  
  # distinct colors per ecoregion
  ecoreg_levels <- sort(unique(df$ecoreg))
  ecoreg_cols <- setNames(hue_pal()(length(ecoreg_levels)), ecoreg_levels)
  
  ggplot(df, aes(x = week_date, y = mean_count,
                 color = ecoreg, fill = ecoreg)) +
    geom_vline(xintercept = vY,
               linetype = "dashed",
               linewidth = 0.3,
               color = "grey50") +
    geom_ribbon(aes(ymin = pmax(0, mean_count - sd_count),
                    ymax = pmin(global_max, mean_count + sd_count)),
                alpha = 0.18,
                color = NA) +
    geom_line(linewidth = 0.5) +
    facet_wrap(~ ecoreg,
               nrow = 3,
               ncol = 1,
               scales = "fixed") +
    scale_color_manual(values = ecoreg_cols) +
    scale_fill_manual(values  = ecoreg_cols) +
    scale_y_continuous(limits = c(0, global_max),
                       expand = c(0, 0)) +
    scale_x_date(limits = c(xMin, xMax),
                 breaks = seq(floor_date(xMin, "year"), ceiling_date(xMax, "year"), by = "1 year"),
                 labels = date_format("%Y"),
                 expand = c(0, 0)) +
    labs(title = "Relative abundance (weekly mean ± SD across IDs) by ecoregion",
         x = "Year", y = "Mean trap counts") +
    theme_minimal(base_size = 11) +
    theme(strip.text = element_text(face = "bold"),
          strip.background = element_rect(fill = "grey85", color = "black"),
          legend.position = "none",
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(color = "black"))
}

# 5) Plot and save
p1 <- plot_full(agg_aedes)
p1
ggsave(filename = file.path(outdir_vis, paste0(datasetName, "_full_timeseries.png")), plot = p1, width = 12, height = 6, dpi = 300)
p1_POT <- plot_full(agg_aedes_POT)
p1_POT
ggsave(filename = file.path(outdir_vis, paste0(datasetName, "_full_timeseries_POT.png")), plot = p1_POT, width = 12, height = 6, dpi = 300)

## PLOT 2: WEEKLY ABUNDANCE BY TIMESERIES SPLIT ##
# 1) Period split for plotting
agg_split <- agg_aedes |>
  mutate(year = year(week_date),
         period = case_when(year >= 2008 & year <= 2012 ~ "2008–2012",
                            year >= 2013 & year <= 2017 ~ "2013–2017",
                            year >= 2018 & year <= 2022 ~ "2018–2022",
                            TRUE ~ NA_character_)) |>
  filter(!is.na(period))

# Do the same for POT
agg_split_POT <- agg_aedes_POT |>
  mutate(year = year(week_date),
         period = case_when(year >= 2008 & year <= 2012 ~ "2008–2012",
                            year >= 2013 & year <= 2017 ~ "2013–2017",
                            year >= 2018 & year <= 2022 ~ "2018–2022",
                            TRUE ~ NA_character_)) |>
  filter(!is.na(period))

# 2) Plot
plot_by_period <- function(df) {
  # same y-axis across all panels (optional but helpful)
  global_max <- y_max_cap
  
  ggplot(df, aes(x = week_date, y = mean_count, color = ecoreg, fill = ecoreg)) +
    geom_ribbon(aes(ymin = pmax(0, mean_count - sd_count),
                    ymax = pmin(global_max, mean_count + sd_count)),
                alpha = 0.18, color = NA) +
    geom_line(linewidth = 0.6) +
    facet_wrap(~ period, nrow = 3, ncol = 1, scales = "free_x") +
    scale_y_continuous(limits = c(0, global_max), expand = c(0, 0)) +
    scale_x_date(breaks = seq(as.Date("2008-01-01"), as.Date("2022-01-01"), by = "1 year"),
                 labels = date_format("%Y"),
                 expand = c(0, 0)) +
    labs(title = "Relative abundance above threshold (n ≥ 55) by timeseries split (weekly mean ± SD across IDs)",
         x = "Year", y = "Mean trap counts",
         color = "Ecoregion", fill = "Ecoregion") +
    theme_minimal(base_size = 11)+
    theme(strip.text = element_text(face = "bold"),
          strip.background = element_rect(fill = "grey85", color = "black"),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          panel.grid.major.x = element_line(linetype = "dashed", color = "grey50"),
          panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(color = "black"))
}

p2 <- plot_by_period(agg_split)
p2
ggsave(filename = file.path(outdir_vis, paste0(datasetName, "_split_timeseries.png")),
       plot = p2, width = 12, height = 7, dpi = 300)
p2_POT <- plot_by_period(agg_split_POT)
p2_POT
ggsave(filename = file.path(outdir_vis, paste0(datasetName, "_split_timeseries_POT.png")),
       plot = p2_POT, width = 12, height = 7, dpi = 300)


## PLOT 3: ON A YEARLY SCALE ##
# 1) Get the mean per year from the weekly series
agg_aedes_year <- agg_aedes |>
  mutate(year_date = as.Date(paste0(year(week_date), "-01-01"))) |>
  group_by(ecoreg, year_date) |>
  summarise(mean_count = mean(mean_count, na.rm = TRUE),
            sd_count   = sd(mean_count,   na.rm = TRUE),
            .groups = "drop") |>
  mutate(sd_count = replace_na(sd_count, 0)) |>
  arrange(ecoreg, year_date)


# Do the same for POT
agg_aedes_year_POT <- agg_aedes_POT |>
  mutate(year_date = as.Date(paste0(year(week_date), "-01-01"))) |>
  group_by(ecoreg, year_date) |>
  summarise(mean_count = mean(mean_count, na.rm = TRUE),
            sd_count   = sd(mean_count,   na.rm = TRUE),
            .groups = "drop") |>
  mutate(sd_count = replace_na(sd_count, 0)) |>
  arrange(ecoreg, year_date)


# 2) Plot helper
plot_yearly <- function(df) {
  
  xMin <- min(df$year_date, na.rm = TRUE)
  xMax <- max(df$year_date, na.rm = TRUE)
  
  split_lines <- as.Date(c("2013-01-01", "2018-01-01"))
  
  global_max <- max(df$mean_count + df$sd_count, na.rm = TRUE)
  
  ecoreg_levels <- sort(unique(df$ecoreg))
  ecoreg_cols <- setNames(hue_pal()(length(ecoreg_levels)), ecoreg_levels)
  
  ggplot(df, aes(x = year_date, y = mean_count, color = ecoreg, fill = ecoreg)) +
    geom_vline(xintercept = split_lines,
               linetype = "dashed",
               linewidth = 0.4,
               color = "black") +
    geom_ribbon(aes(ymin = pmax(0, mean_count - sd_count),
                    ymax = mean_count + sd_count),
                alpha = 0.18,
                color = NA) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~ ecoreg, nrow = 3, ncol = 1, scales = "fixed") +
    scale_color_manual(values = ecoreg_cols) +
    scale_fill_manual(values = ecoreg_cols) +
    scale_y_continuous(limits = c(0, global_max),
                       expand = c(0, 0)) +
    scale_x_date(limits = c(xMin, xMax),
                 breaks = seq(xMin, xMax, by = "1 year"),
                 labels = date_format("%Y"),
                 expand = c(0, 0)) +
    labs(title = "Annual mean relative abundance (of weekly mean trap counts) above threshold (n ≥ 55) by ecoregion",
         x = "Year", y = "Mean trap counts") +
    theme_minimal(base_size = 11) +
    theme(strip.text = element_text(face = "bold"),
          strip.background = element_rect(fill = "grey85", color = "black"),
          legend.position = "none",
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
          panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(color = "black"))
}

# 3) Plot and save
p3 <- plot_yearly(agg_aedes_year)
p3
ggsave(filename = file.path(outdir_vis, paste0(datasetName, "_yearly_scale.png")),
       plot = p3, width = 12, height = 6, dpi = 300)
p3_POT <- plot_yearly(agg_aedes_year_POT)
p3_POT
ggsave(filename = file.path(outdir_vis, paste0(datasetName, "_yearly_scale_POT.png")),
       plot = p3_POT, width = 12, height = 6, dpi = 300)

## 4: PLOT YEARLY ABUNDANCE WITH  TIMESERIES SPLITS ## 
# 1) Add period labels to yearly data
agg_year_split <- agg_aedes_year |>
  mutate(year = year(year_date),
         period = case_when(year >= 2008 & year <= 2012 ~ "2008–2012",
                            year >= 2013 & year <= 2017 ~ "2013–2017",
                            year >= 2018 & year <= 2022 ~ "2018–2022",
                            TRUE ~ NA_character_)) |>
  filter(!is.na(period))

# Do the same for POT
agg_year_split_POT <- agg_aedes_year_POT |>
  mutate(year = year(year_date),
         period = case_when(year >= 2008 & year <= 2012 ~ "2008–2012",
                            year >= 2013 & year <= 2017 ~ "2013–2017",
                            year >= 2018 & year <= 2022 ~ "2018–2022",
                            TRUE ~ NA_character_)) |>
  filter(!is.na(period))

# 2) Plot helper
plot_yearly_by_period <- function(df) {
  
  global_max <- max(df$mean_count + df$sd_count, na.rm = TRUE)
  
  ecoreg_levels <- sort(unique(df$ecoreg))
  ecoreg_cols <- setNames(hue_pal()(length(ecoreg_levels)), ecoreg_levels)
  
  ggplot(df, aes(x = year_date, y = mean_count, color = ecoreg, fill = ecoreg)) +
    geom_ribbon(aes(ymin = pmax(0, mean_count - sd_count),
                    ymax = mean_count + sd_count),
                alpha = 0.18,
                color = NA) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~ period, nrow = 3, ncol = 1, scales = "free_x") +
    scale_color_manual(values = ecoreg_cols) +
    scale_fill_manual(values = ecoreg_cols) +
    scale_y_continuous(limits = c(0, global_max),
                       expand = c(0, 0)) +
    scale_x_date(breaks = seq(as.Date("2008-01-01"), as.Date("2022-01-01"), by = "1 year"),
                 labels = date_format("%Y"),
                 expand = c(0, 0)) +
    labs(title = "Annual mean relative abundance (of weekly mean trap counts) by timeseries split",
         x = "Year", y = "Mean trap counts",
         color = "Ecoregion", fill = "Ecoregion") +
    theme_minimal(base_size = 11) +
    theme(strip.text = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey85", color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major.x = element_line(linetype = "dashed", color = "grey50"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"))
}

# 3) Plot and save
p4 <- plot_yearly_by_period(agg_year_split)
p4
ggsave(filename = file.path(outdir_vis, paste0(datasetName, "_yearly_split_timeseries.png")),
       plot = p4, width = 12, height = 7, dpi = 300)
p4_POT <- plot_yearly_by_period(agg_year_split_POT)
p4_POT
ggsave(filename = file.path(outdir_vis, paste0(datasetName, "_yearly_split_timeseries_POT.png")),
       plot = p4_POT, width = 12, height = 7, dpi = 300)

## PLOT 5: MEAN OF EACH TIMESERIES SPLIT ## 
# 1) Summarise yearly means within each split
agg_year_split_mean <- agg_year_split |>
  group_by(ecoreg, period) |>
  summarise(mean_count = mean(mean_count, na.rm = TRUE),
            sd_count   = sd(mean_count, na.rm = TRUE),
            .groups = "drop") |>
  mutate(period = factor(period, levels = c("2008–2012", "2013–2017", "2018–2022")),
         sd_count = replace_na(sd_count, 0))

# Do the same for POT
agg_year_split_mean_POT <- agg_year_split_POT |>
  group_by(ecoreg, period) |>
  summarise(mean_count = mean(mean_count, na.rm = TRUE),
            sd_count   = sd(mean_count, na.rm = TRUE),
            .groups = "drop") |>
  mutate(period = factor(period, levels = c("2008–2012", "2013–2017", "2018–2022")),
         sd_count = replace_na(sd_count, 0))

# 2) Plot helper
plot_split_means <- function(df) {
  
  ecoreg_cols <- c(
    "Alpine" = "firebrick",
    "Continental" = "forestgreen",
    "Mediterranean" = "steelblue"
  )
  
  ggplot(df, aes(x = period, y = mean_count, color = ecoreg, group = ecoreg)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    labs(title = "Relative mean abundance above threshold (n ≥ 55) across timeseries splits",
         x = "Timeseries split", y = "Mean trap counts",
         color = "Ecoregion") +
    scale_color_manual(values = ecoreg_cols) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(color = "black"),
          panel.grid.minor = element_blank(),
          legend.position = "right",
          legend.title = element_text(face = "bold"),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6))
}

# 3) Plot and save
p5 <- plot_split_means(agg_year_split_mean)
p5
ggsave(filename = file.path(outdir_vis, paste0(datasetName, "_split_means.png")),
       plot = p5, width = 8, height = 5, dpi = 300)
p5_POT <- plot_split_means(agg_year_split_mean_POT)
p5_POT
ggsave(filename = file.path(outdir_vis, paste0(datasetName, "_split_means_POT.png")),
       plot = p5_POT, width = 8, height = 5, dpi = 300)