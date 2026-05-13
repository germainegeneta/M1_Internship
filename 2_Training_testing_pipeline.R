## VectAbundance pipeline for modelling

# 0. Set up ----
## 0.1 Load libraries ----
library(tidyverse) ## Version ‘2.0.0’
library(caret) ## Version ‘6.0.94’
library(CAST) ## Version ‘1.0.2’
library(ranger) ## Version ‘0.16.0’
library(correlation) ## Version ‘0.8.5’
library(ggplot2)
library(mlr3measures)
library(ISOweek)
library(vip)
library(pdp)
library(ggtext)
library(ggnewscale)
library(lubridate)
library(scales)
library(forecast)
library(iml)
library(patchwork)
library(precrec)

## 0.2 Create output directory ----
datasetName <- "VectAbundance_2024"
outdir_pred <- file.path(datasetName, "prediction_res")
dir.create(outdir_pred, recursive = TRUE, showWarnings = FALSE)

## 0.3 Load climate data ----
# Load and merge aggregated climate data from CCM script
raw_glob <- read.csv(file.path(datasetName, "CCM_outputs", "fitmodels", "Global_aggclim.csv"))
raw_alp <- read.csv(file.path(datasetName, "CCM_outputs", "fitmodels", "Alpine_aggclim.csv"))
raw_con <- read.csv(file.path(datasetName, "CCM_outputs", "fitmodels", "Continental_aggclim.csv"))
raw_med <- read.csv(file.path(datasetName, "CCM_outputs", "fitmodels", "Mediterranean_aggclim.csv"))
# Put in a list
raw <- list(data_glob = raw_glob, data_alp = raw_alp, data_con = raw_con, data_med = raw_med)

# Reformat
data_list <- list()
for (nm in names(raw)) {
  cat("\nProcessing:", nm, "\n")
  df <- raw[[nm]] |> 
    mutate(lag_id = paste0(lag_start, "_", lag_end)) |> 
    arrange(lag_start, lag_end) |> 
    rename(AV_TEMP = mean_temp_mean, # rename clim vars
           MX_TEMP = max_temp_mean,
           MN_TEMP = min_temp_mean,
           GDD = gdd_sum,
           SM_RAIN = sum_precip_sum,
           AV_RH = mean_rh_mean) |> 
    pivot_wider(id_cols = c(ID, date, ecoreg, year, week, count, count_c, lat, long),
                names_from = lag_id,
                values_from = c(AV_TEMP, MX_TEMP, MN_TEMP,
                                GDD, SM_RAIN, AV_RH),
                names_glue = "{.value}_{lag_id}") |> 
    relocate(ID, date, ecoreg, year, week, count, count_c, lat, long,
             starts_with("AV_TEMP"),
             starts_with("MX_TEMP"),
             starts_with("MN_TEMP"),
             starts_with("GDD"),
             starts_with("SM_RAIN"),
             starts_with("AV_RH"))
  
  # Store
  data_list[[nm]] <- df
}
# list2env(data_list, envir = .GlobalEnv)

# Check climate data trends
stripchart(SM_RAIN_1_1 ~ ecoreg, data = data_alp, method = "jitter")

# 1. Calculate and add Fourier terms, per timeseries ----
data_model_list <- list()
for (nm in names(data_list)) {
  cat("\nAdding Fourier terms to:", nm, "\n")
  df <- data_list[[nm]] |> 
    select(-date, -count_c, -lat, -long) |> 
    arrange(year, week)
  
  # Compute Fourier; adopted from Da Re et al.'s code
  yt.ts <- ts(df$count, start = 1, frequency = 52)
  fourier_var <- as.data.frame(fourier(yt.ts, K = 2))
  fourier_col <- colnames(fourier_var)
  
  df <- cbind(df, fourier_var) |>
    relocate(all_of(fourier_col), .before = "AV_TEMP_1_1") |> 
    rename("S1" = "S1-52", 
           "S2" = "S2-52",
           "C1" = "C1-52",
           "C2" = "C2-52")
  # Store
  eco <- sub("^data_", "", nm)
  data_model_list[[paste0(eco, "_model")]] <- df
  # Save as .csv just in case
  write_csv(df, file.path(outdir_pred, paste0(eco, "_model_data.csv")))
  
}
list2env(data_model_list, envir = .GlobalEnv)

# 2. Abundance model variables ----
## 2.1 Select variables of the best time lags, including (1,1) time lag ----
# Step 2: to evaluate the correlation between these variables
# Step 3:to select the variables not correlated with the highest sense ecological. 
# The first selection is crossed with the other selection done with the VIF with the corSelect function of the fuzzySim package 
# to select variables with the lowest VIF. The final selection is a mixed of both methods.

# Global model
glob_preds <- c("S1", "C1", "S2", "C2",
                "AV_TEMP_3_1", "AV_TEMP_1_1",
                "MX_TEMP_4_1", "MX_TEMP_1_1",
                "MN_TEMP_3_1", "MN_TEMP_1_1",
                "GDD_3_1", "GDD_1_1",
                "SM_RAIN_8_1", "SM_RAIN_1_1",
                "AV_RH_6_1", "AV_RH_1_1")
glob_model <- glob_model |> select("ecoreg", "year", "week", "count", all_of(glob_preds))

# Alpine model
alp_preds <- c("S1", "C1", "S2", "C2",
               "AV_TEMP_6_1", "AV_TEMP_1_1",
               "MX_TEMP_6_1", "MX_TEMP_1_1",
               "MN_TEMP_5_1", "MN_TEMP_1_1",
               "GDD_6_1", "GDD_1_1",
               "SM_RAIN_8_4", "SM_RAIN_1_1",
               "AV_RH_6_4", "AV_RH_1_1")
alp_model <- alp_model |> select("ecoreg", "year", "week", "count", all_of(alp_preds))

# Continental model
con_preds <- c("S1", "C1", "S2", "C2",
               "AV_TEMP_4_1", "AV_TEMP_1_1",
               "MX_TEMP_4_1", "MX_TEMP_1_1",
               "MN_TEMP_3_1", "MN_TEMP_1_1",
               "GDD_4_1", "GDD_1_1",
               "SM_RAIN_6_1", "SM_RAIN_1_1",
               "AV_RH_5_1", "AV_RH_1_1")
con_model <- con_model |> select("ecoreg", "year", "week", "count", all_of(con_preds))

# Mediterranean model
med_preds <- c("S1", "C1", "S2", "C2",
               "AV_TEMP_5_1", "AV_TEMP_1_1",
               "MX_TEMP_5_1", "MX_TEMP_1_1",
               "MN_TEMP_4_1", "MN_TEMP_1_1",
               "GDD_4_1", "GDD_1_1",
               "SM_RAIN_8_1", "SM_RAIN_1_1",
               "AV_RH_8_1", "AV_RH_1_1")
med_model <- med_model |> select("ecoreg", "year", "week", "count", all_of(med_preds))

## 2.2 Plot the bivariate relationship between abundance and climate variables ----
# Compile variables
pred_list <- list(Global = glob_preds, Alpine = alp_preds, Mediterranean = med_preds, Continental = con_preds)

# Compile model data
model_list <- list(Global = glob_model, Alpine = alp_model, Mediterranean = med_model, Continental = con_model)

# Select colors (point + smooth)
color_list <- list(Global = c("gold", "goldenrod2"),
                   Alpine = c("lightcoral", "firebrick"),
                   Continental = c("palegreen3", "forestgreen"),
                   Mediterranean = c("skyblue", "steelblue"))

# Run plotting loop
plot_list <- list()
for (eco in names(pred_list)) {
  preds <- pred_list[[eco]]
  preds <- setdiff(preds, c("S1", "C1", "S2", "C2")) # don't plot Fourier
  cols  <- color_list[[eco]]
  df    <- model_list[[eco]]
  
  # long format (raw)
  df_long_raw <- df |>
    select(count, all_of(preds), year) |>
    pivot_longer(-c(count, year))
  # bin + count for bubble size
  df_long <- df_long_raw |>
    mutate(value_bin = round(value, 1),   # adjust if needed
           count_bin = round(count, 0)) |> 
    group_by(name, value_bin, count_bin) |>
    summarise(n = n(), .groups = "drop")
  # plot
  p <- ggplot() +
    geom_point(data = df_long, # Bubble points (frequency)
               aes(x = value_bin, y = count_bin, size = n),
               color = cols[1], alpha = 0.7) +
    geom_smooth(data = df_long_raw, # Smooth the raw data
                aes(x = value, y = count),
                color = cols[2], 
                se = FALSE) +
    facet_wrap(~name, scales = "free_x") +
    labs(title = paste("Bivariate relationship of abundance and selected climate variables in the", eco, "ecoregion"),
         x = "Value of climate variable", y = "Relative abundance (trap count)", size = "N") +
    ylim(c(0, 80)) +
    scale_size(range = c(0.5, 8)) +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "grey95", color = "grey40"),
      strip.text = element_text(face = "bold", size = 10),
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "bottom")
  
  plot_list[[eco]] <- p
  print(p)
  # save
  ggsave(filename = file.path(outdir_pred, paste0(eco, "_bivarplot.png")),
         plot = p, width = 12, height = 6, dpi = 300)
}

## 2.3 Identify correlated variables (greater than Pearson 0.8) ----
pear <- list()
for (eco in names(model_list)) {
  df    <- model_list[[eco]]
  preds <- pred_list[[eco]]
  
  m <- cor(df[, preds], method = "pearson", use = "pairwise.complete.obs") # correlation matrix
  
  index <- which(abs(m) > 0.8 & abs(m) < 1, arr.ind = TRUE) # find high correlations
  df_cor <- subset(as.data.frame(index), row < col)
  
  p <- cbind.data.frame(var1 = rownames(m)[df_cor$row], # build result dataframe
                        var2 = colnames(m)[df_cor$col],
                        corr = m[cbind(df_cor$row, df_cor$col)])
  
  pear[[eco]] <- p # store
  
  cat("\n---", eco, "---\n") # print
  print(p)
}

## 2.4 Select final variables based on signed R2 on CCM ----
# Don't forget (1,1) lags
# For explainability: did not use GDD_SUM
glob_preds_pear <- c("ecoreg", "S1", "C1", "S2", "C2", "MX_TEMP_4_1", "SM_RAIN_8_1", "SM_RAIN_1_1", "AV_RH_6_1", "AV_RH_1_1")
alp_preds_pear <- c("S1", "C1", "S2", "C2", "MN_TEMP_5_1", "SM_RAIN_8_4", "SM_RAIN_1_1", "AV_RH_6_4", "AV_RH_1_1")
con_preds_pear <- c("S1", "C1", "S2", "C2", "MN_TEMP_3_1", "SM_RAIN_6_1", "SM_RAIN_1_1", "AV_RH_5_1", "AV_RH_1_1")
med_preds_pear <- c("S1", "C1", "S2", "C2", "MX_TEMP_5_1", "SM_RAIN_8_1", "SM_RAIN_1_1", "AV_RH_8_1", "AV_RH_1_1")

# Create final data frame for the multivariate analysis
glob_model <- glob_model |> select("ecoreg", "year", "week", "count", all_of(glob_preds_pear))
alp_model <- alp_model |> select("ecoreg", "year", "week", "count", all_of(alp_preds_pear))
con_model <- con_model |> select("ecoreg", "year", "week", "count", all_of(con_preds_pear))
med_model <- med_model |> select("ecoreg", "year", "week", "count", all_of(med_preds_pear))

# 3. Multivariate analysis using leave-one-year-out CV ----
## 3.1 Compile model data ----
model_list_mv <- list(Global = glob_model,
                      Alpine = alp_model,
                      Continental = con_model,
                      Mediterranean = med_model)

## 3.2 Compile model variables ----
pred_list_mv <- list(Global = glob_preds_pear,
                     Alpine = alp_preds_pear,
                     Continental = con_preds_pear,
                     Mediterranean = med_preds_pear)

## 3.3 Define Spearman correlation function ----
spearmcor <- function(data, lev = NULL, model = NULL) {
  out <- cor(x = data$pred, y = data$obs, method = "spearman")
  names(out) <- "spearman"
  out
}

## 3.4 Define cross validation column (year) ----
cv_col <- "year"

## 3.5 Run training and testing loop ----
# To preserve computing power and troubleshoot, I will only run 2 ecoregions at a time
model_list_mv_1 <- list(Global = glob_model, Alpine = alp_model)
model_list_mv_2 <- list(Continental = con_model, Mediterranean = med_model)
pred_list_mv_1 <- list(Global = glob_preds_pear, Alpine = alp_preds_pear)
pred_list_mv_2 <- list(Continental = con_preds_pear, Mediterranean = med_preds_pear)

# Store results here
results_list <- list()
# if you also split into 2, do not run this line again or else it will clear the list

# Loop
for (name in names(model_list_mv_2)) {
  cat("\nRunning model:", name, "\n")
  
  df <- model_list_mv_2[[name]]
  preds <- pred_list_mv_2[[name]]
  
  # Log transform trap counts
  df$count <- log(df$count) 
  
  # Create CV folds
  indices_cv <- CAST::CreateSpacetimeFolds(df, 
                                           timevar = cv_col,
                                           k = length(unique(df[[cv_col]])))
  # trainControl
  tr <- trainControl(method = "cv",
                     index = indices_cv$index,
                     indexOut = indices_cv$indexOut,
                     savePredictions = "final",
                     summaryFunction = spearmcor)
  
  # Train model
  mod <- caret::train(x = df[, preds],
                      y = df$count,
                      method = "ranger",
                      tuneLength = 10,
                      trControl = tr,
                      metric = "spearman",
                      maximize = TRUE,
                      preProcess = c("center", "scale"),
                      importance = "permutation",
                      verbose = TRUE)
  
  # Attach predictions
  df$rowIndex <- seq_len(nrow(df))
  
  df_cv <- mod$pred |>
    left_join(df) |>
    select(pred, obs, ecoreg, week, year)
  
  # Store results
  res <- list(model = mod,
              df_cv = df_cv,
              df_mod = df)
  results_list[[name]] <- res
  
  # Save RDS
  saveRDS(res, file.path(outdir_pred, paste0(name, "_multiv_model_nowcasting2.rds")))
}

# OPTIONAL: assign all to environment
list2env(results_list, envir = .GlobalEnv)

# 4. Plot observed vs predicted ----
## 4.1. Year-faceted point-line ----
plot_model_list_year <- list()
for (name in names(results_list)) {
  cat("\nPlotting:", name, "\n")
  
  eco <- name
  # Skip Global plot
  if (eco == "Global") next
  
  # Aggregate observed abundance for plotting
  df_obs <- model_list[[eco]] |>
    dplyr::group_by(year, week) |>
    dplyr::summarise(obs_mean = mean(count, na.rm = TRUE),
                     obs_sd   = sd(count, na.rm = TRUE),
                     n        = dplyr::n(),
                     .groups  = "drop") |> 
    dplyr::mutate(se = obs_sd / sqrt(n)) # Compute SE
  
  # Back transform regional model prediction counts (ID-level)
  df_reg_raw <- results_list[[eco]]$df_cv |>
    dplyr::filter(ecoreg == eco) |>
    dplyr::mutate(pred_exp = exp(pred),
                  obs_exp  = exp(obs))
  
  # Aggregate regional model predictions for plotting
  df_reg <- df_reg_raw |>
    dplyr::group_by(year, week) |>
    dplyr::summarise(pred_reg_mean = mean(pred_exp, na.rm = TRUE),
                     .groups = "drop")
  
  # Merge for plotting
  df_plot <- df_obs |>
    dplyr::left_join(df_reg,  by = c("year", "week"))
  
  # Compute MAE and rho on back transformed, ID-level counts
  metrics <- df_reg_raw |>
    dplyr::group_by(year) |>
    dplyr::summarise(MAE_reg  = round(mae(obs_exp, pred_exp), 1),
                     RHO_reg  = round(cor(obs_exp, pred_exp, method = "spearman"), 2),
                     .groups = "drop") |> 
    dplyr::mutate(label = paste0("MAE = <span style='color:red;'>", MAE_reg, "</span>, ",
                                 "ρ = <span style='color:red;'>", RHO_reg, "</span>"))
  
  # Plot
  p <- ggplot() +
    # SE bars
    geom_errorbar(data = df_plot,
                  aes(x = week, ymin = obs_mean - se, ymax = obs_mean + se),
                  width = 0.25,
                  color = "grey80",
                  linewidth = 0.5) +
    # Observed line
    geom_line(data = df_plot,
              aes(x = week, y = obs_mean,
              color = "Observed abundance"),
              linewidth = 1) +
    # Observed points
    geom_point(data = df_plot,
               aes(x = week, y = obs_mean, size = n),
               color = "grey80",
               alpha = 0.8) +
    # Model lines
    geom_line(data = df_plot,
              aes(x = week, y = pred_reg_mean, color = "Predicted abundance (Local model)"),
              linewidth = 1) +
    facet_wrap(~year, ncol = 3, scales = "free_y") +
    # Metrics
    ggtext::geom_richtext(data = metrics,
                          aes(x = Inf, y = Inf, label = label),
                          inherit.aes = FALSE,
                          hjust = 1.1, vjust = 1.1,
                          size = 3.5, fill = NA, label.color = NA) +
    scale_color_manual(values = c("Observed" = "grey60", "Predicted abundance (Local model)" = "red"), name = NULL) +
    scale_size(range = c(0.15, 6), name = "No. of traps") +
    labs(title = paste("Observed vs. predicted abundance of Local models in the", eco, "ecoregion per year"),
         x = "Week", y = "Abundance") +
    theme(panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(color = "grey95"),
          panel.border = element_rect(color = "black", fill = NA),
          strip.background = element_rect(fill = "grey90", color = "black"),
          strip.text = element_text(face = "bold"),
          axis.text.y = element_text(angle = 90, hjust = 0),
          legend.position = "bottom",
          plot.title = element_text(face = "bold"))
  
  # Store + print + save
  plot_model_list_year[[name]] <- p
  print(p)
  ggsave(filename = file.path(outdir_pred, paste0(name, "_model_prediction_year.png")),
         plot = p, width = 13, height = 9, dpi = 300)
}

## 4.2 Full timeseries  ----
# Note: I made 4 versions of the plot
# Ver. 1: compressed year (x) axis, does not show no activity weeks, y scale not fixed, no SE
# Ver. 2: x axis accounts for missing weeks (shows seasonal activity), y scale not fixed, no SE
# Ver. 3: x axis accounts for missing weeks (shows seasonal activity), y scale is fixed and truncated/clipped to value of 300, no SE
# Ver. 4: x axis accounts for missing weeks (shows seasonal activity), y scale is fixed and truncated/clipped to value of 300, with SE
### Bar plot ----
plot_model_list_full <- list()
combined_list <- list()
for (name in names(results_list)) {

  eco <- name

  if (eco == "Global") next

  # # Observed abundance for ver. 1, 2, 3
  # df_obs <- model_list[[eco]] |>
  #   dplyr::group_by(year, week) |>
  #   dplyr::summarise(obs_mean = mean(count, na.rm = TRUE),
  #                    n = dplyr::n(),
  #                    .groups = "drop")
  # Observed abundance for ver. 4
  df_obs <- model_list[[eco]] |>
    dplyr::group_by(year, week) |>
    dplyr::summarise(
      obs_mean = mean(count, na.rm = TRUE),
      obs_sd = sd(count, na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop") |>
    dplyr::mutate(se = obs_sd / sqrt(n)) # calculate SE

  # Local model (ID-level)
  df_reg_raw <- results_list[[eco]]$df_cv |>
    dplyr::filter(ecoreg == eco) |>
    dplyr::mutate(pred_reg_exp = exp(pred),
                  obs_exp = exp(obs))

  df_reg <- df_reg_raw |>
    dplyr::group_by(year, week) |>
    dplyr::summarise(pred_reg_mean = mean(pred_reg_exp, na.rm = TRUE),
                     .groups = "drop")

  # # Merge (for ver. 1)
  # df_plot <- df_obs |>
  #   dplyr::left_join(df_reg, by = c("year", "week")) |>
  #   dplyr::arrange(year, week)

  # Merge + keep missing weeks (for ver. 2-4)
  df_plot <- df_obs |>
    dplyr::left_join(df_reg, by = c("year", "week")) |>
    tidyr::complete(year, week = 1:52) |>
    dplyr::arrange(year, week)

  df_plot$ecoreg <- eco

  combined_list[[eco]] <- df_plot
}

combined_df <- dplyr::bind_rows(combined_list)


# Continuous time index
combined_df <- combined_df |>
  dplyr::group_by(ecoreg) |>
  dplyr::arrange(year, week) |>
  dplyr::mutate(
    time_id = dplyr::row_number()
  ) |>
  dplyr::ungroup()


# Metrics
metrics <- dplyr::bind_rows(
  lapply(names(results_list), function(eco) {

    if (eco == "Global") return(NULL)

    df_reg_raw <- results_list[[eco]]$df_cv |>
      dplyr::filter(ecoreg == eco) |>
      dplyr::mutate(pred = exp(pred),
                    obs = exp(obs))

    data.frame(ecoreg = eco,
               MAE = round(mae(df_reg_raw$obs, df_reg_raw$pred), 1),
      RHO = round(cor(df_reg_raw$obs, df_reg_raw$pred, method = "spearman"), 2))
  })
)

metrics <- metrics |>
  dplyr::mutate(label = paste0("MAE = <span style='color:red;'>", MAE, "</span>",
                               ", ρ = <span style='color:red;'>", RHO, "</span>"))
# Year separators
year_breaks <- combined_df |>
  dplyr::group_by(ecoreg, year) |>
  dplyr::summarise(
    x = min(time_id),
    .groups = "drop")

# Plot
p <- ggplot() +
  # Observed bars
  geom_col(data = combined_df,
           aes(x = time_id, y = obs_mean, fill = n),
           width = 1, alpha = 0.8, na.rm = TRUE) +
  # SE error bars - only for ver. 4
  geom_errorbar(
    data = combined_df,
    aes(
      x = time_id,
      ymin = obs_mean - se,
      ymax = obs_mean + se,
      color = "Observed ± SE"
    ),
    width = 0.15,
    linewidth = 0.3,
    alpha = 0.8,
    na.rm = TRUE
  ) +

  # Local model line
  geom_line(data = combined_df,
            aes(x = time_id, y = pred_reg_mean, color = "Local model"),
            linewidth = 1, na.rm = TRUE) +
  # Year separators
  geom_vline(data = year_breaks,
             aes(xintercept = x),
             linetype = "dashed", color = "grey70") +
  # Year labels
  geom_text(data = year_breaks,
            aes(x = x, y = -Inf, label = year),
            inherit.aes = FALSE, vjust = 1.25,
            size = 3, color = "grey30") +
  # Metrics
  ggtext::geom_richtext(data = metrics,
                        aes(x = Inf, y = Inf, label = label),
                        inherit.aes = FALSE,
                        hjust = 1.1, vjust = 1.1,
                        size = 3.5, fill = NA, label.color = NA) +
  # Facet
  #facet_wrap(~ecoreg, ncol = 1, scales = "free") + # for ver. 1
  facet_wrap(~ecoreg, ncol = 1, scales = "fixed") + # for ver. 2-4
  # Colors
  scale_color_manual(values = c("Observed ± SE" = "black", # only for ver. 4
                                "Local model" = "red"),
                     name = NULL) +
  scale_fill_gradient(low = "grey90", high = "grey30", name = "No. of traps") +
  # X-axis
  scale_x_continuous(expand = c(0.005, 0), breaks = NULL) +
  # Labels
  labs(title = "Observed vs. predicted abundance of Local models across ecoregions",
       #subtitle = "Observed abundance shown as bars; shade intensity indicates sampling intensity",
       x = NULL, y = "Abundance") +
  # Coords
  #coord_cartesian(clip = "off") + # for ver. 1 and 2
  coord_cartesian(ylim = c(0, 300), clip = "on") + # for ver. 3-4
  # Theme
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "grey95"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black",fill = NA),
    panel.spacing.y = unit(1.2, "lines"),
    plot.margin = margin(10, 10, 20, 10),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(angle = 90),
    axis.ticks.x = element_blank())

print(p)

plot_model_list_full[["full_timeseries_v4"]] <- p # change version
ggsave(filename = file.path(outdir_pred, "full_timeseries_v4.png"), # change version
  plot = p, width = 12, height = 8, dpi = 300)

# 5. VIP ----
# Store results here
vip_list <- list()

## 5.1 Extract variable importance per ecoregion ----
for (eco in names(results_list)) {
  if (eco == "Global") next
  # Extract model
  model <- results_list[[eco]]$model
  # Extract variable importance
  imp <- model$finalModel$variable.importance
  imp <- as.data.frame(imp)
  imp$var <- rownames(imp)
  # Rename
  imp <- imp  |> 
    dplyr::rename(importance = imp) |> 
    dplyr::arrange(dplyr::desc(importance))
  # Add ecoregion
  imp$ecoreg <- eco
  # Store
  vip_list[[eco]] <- imp
}

#Combine all ecoregions
vip_df <- dplyr::bind_rows(vip_list)

# Reorder variables
vip_df$var <- forcats::fct_reorder(vip_df$var,
                                   vip_df$importance,
                                   .fun = max,
                                   .desc = TRUE)

## 5.2 Plot VIP ----
plot_imp_abundance <- ggplot(vip_df,
                             aes(x = var, y = importance, fill = ecoreg)) +
  # Bars
  geom_col(position = position_dodge(width = 0.8),
           width = 0.7) +
  # Colors by ecoregion
  scale_fill_manual(values = c("Alpine" = "firebrick",
                               "Continental" = "forestgreen",
                               "Mediterranean" = "steelblue")) +
  # Labels
  labs(title = "Variable importance of Local models",
       x = NULL, y = "Variable importance", fill = "Ecoregion") +
  # Theme similar to Taconet et al. figure
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45,
                               hjust = 1,
                               size = 9),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom")
print(plot_imp_abundance)

# Save
ggsave(filename = file.path(outdir_pred, "VIP_all_ecoregions_taconet_style.png"),
       plot = plot_imp_abundance, width = 10, height = 6, dpi = 300)

# 6. PDP ----
# Store results
pdp_list <- list()
rug_list <- list()

## Loop through ecoregions ----
for (eco in names(results_list)) {
  
  if (eco == "Global") next
  cat("\nGenerating PDPs for:", eco, "\n")
  # Extract model and data
  model  <- results_list[[eco]]$model
  df_mod <- results_list[[eco]]$df_mod
  # Set response variable
  response_var <- "count"
  # Create iML predictor
  predictor <- Predictor$new(
    model = model,
    data  = df_mod |> dplyr::select(-all_of(response_var)),
    y     = df_mod[[response_var]])
  # Get variable names
  vars <- names(model$finalModel$variable.importance)
  # Generate PDP
  pdp_all <- purrr::map_dfr(vars, function(v) {
    cat(" >", v, "\n")
    # Generate PDP
    pd <- FeatureEffect$new(
      predictor,
      feature = v,
      method = "pdp",
      grid.size = 50)
    # Extract PDP dataframe
    pd_df <- pd$results
    # Standardize names
    names(pd_df)[1] <- "x"
    names(pd_df)[2] <- "y"
    # Add metadata
    pd_df$variable <- v
    pd_df$ecoreg   <- eco
    # Create rug plot
    rug_df <- data.frame(
      rug_x   = df_mod[[v]],
      variable = v,
      ecoreg   = eco)
    rug_list[[paste(eco, v, sep = "_")]] <<- rug_df
    return(pd_df)
  })
  
  # Store PDPs
  pdp_list[[eco]] <- pdp_all
}

# Combine data
pdp_df <- bind_rows(pdp_list)
rug_df_all <- bind_rows(rug_list)

# Order variables
var_order <- unique(pdp_df$variable)
pdp_df$variable <- factor(
  pdp_df$variable,
  levels = var_order)

rug_df_all$variable <- factor(
  rug_df_all$variable,
  levels = var_order)

# Plot
p_pdp <- ggplot(pdp_df, 
                aes(x = x, y = y, color = ecoreg)) +
  # PDP lines
  geom_line(linewidth = 1) +
  # Rug plot
  geom_rug(data = rug_df_all, aes(x = rug_x, color = ecoreg),
           inherit.aes = FALSE,
           alpha = 0.2,
           sides = "b") +
  # Facets by predictor
  facet_wrap(~variable, scales = "free_x") +
  # Colors by ecoregion
  scale_color_manual(values = c("Alpine" = "firebrick",
                                "Continental" = "forestgreen",
                                "Mediterranean" = "steelblue")) +
  # Labels
  labs(title = "Partial dependence plots of Local models",
       x = NULL, y = "Predicted abundance", color = "Ecoregion") +
  # Theme
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey90", color = "black"),
        strip.text = element_text(face = "bold", size = 8),
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        legend.position = "bottom")
print(p_pdp)

# Save
ggsave(filename = file.path(outdir_pred, "PDP_all_ecoregions_iml.png"),
       plot = p_pdp, width = 14, height = 10, dpi = 300)

# # Supplementary figures ----
# # Output directory for supplementary figures
# outdir_supp <- file.path(outdir_pred, "supplementary_figs")
# dir.create(outdir_supp, recursive = TRUE, showWarnings = FALSE)
# 
# ##  Full timeseries  ----
# ### Bar plot ----
# plot_model_list_full <- list()
# combined_list <- list()
# 
# for (name in names(results_list)) {
#   
#   eco <- name
#   
#   if (eco == "Global") next
#   
#   # Observed abundance for ver. 1, 2, 3
#   df_obs <- model_list[[eco]] |>
#     dplyr::group_by(year, week) |>
#     dplyr::summarise(obs_mean = mean(count, na.rm = TRUE),
#                      n = dplyr::n(),
#                      .groups = "drop")
#   
#   # # Observed abundance for ver. 4
#   # df_obs <- model_list[[eco]] |>
#   #   dplyr::group_by(year, week) |>
#   #   dplyr::summarise(
#   #     obs_mean = mean(count, na.rm = TRUE),
#   #     obs_sd = sd(count, na.rm = TRUE),
#   #     n = dplyr::n(),
#   #     .groups = "drop") |>
#   #   dplyr::mutate(se = obs_sd / sqrt(n)) # calculate SE
#   
#   # Local model (ID-level)
#   df_reg_raw <- results_list[[eco]]$df_cv |>
#     dplyr::filter(ecoreg == eco) |>
#     dplyr::mutate(
#       pred_reg_exp = exp(pred),
#       obs_exp = exp(obs))
#   
#   df_reg <- df_reg_raw |>
#     dplyr::group_by(year, week) |>
#     dplyr::summarise(
#       pred_reg_mean = mean(pred_reg_exp, na.rm = TRUE),
#       .groups = "drop")
#   
#   # Global model (filtered by ecoregion)
#   df_glob_raw <- results_list[["Global"]]$df_cv |>
#     dplyr::filter(ecoreg == eco) |>
#     dplyr::mutate(
#       pred_glob_exp = exp(pred),
#       obs_exp = exp(obs))
#   
#   df_glob <- df_glob_raw |>
#     dplyr::group_by(year, week) |>
#     dplyr::summarise(
#       pred_glob_mean = mean(pred_glob_exp, na.rm = TRUE),
#       .groups = "drop")
#   
#   # Merge (for ver. 1)
#   df_plot <- df_obs |>
#     dplyr::left_join(df_reg, by = c("year", "week")) |>
#     dplyr::left_join(df_glob, by = c("year", "week")) |>
#     dplyr::arrange(year, week)
#   
#   # # Merge + keep missing weeks (for ver. 2-4)
#   # df_plot <- df_obs |>
#   #   dplyr::left_join(df_reg, by = c("year", "week")) |>
#   #   dplyr::left_join(df_glob, by = c("year", "week")) |>
#   #   tidyr::complete(year, week = 1:52) |>
#   #   dplyr::arrange(year, week)
#   
#   df_plot$ecoreg <- eco
#   
#   combined_list[[eco]] <- df_plot
# }
# 
# combined_df <- dplyr::bind_rows(combined_list)
# 
# 
# # Continuous time index
# combined_df <- combined_df |>
#   dplyr::group_by(ecoreg) |>
#   dplyr::arrange(year, week) |>
#   dplyr::mutate(
#     time_id = dplyr::row_number()
#   ) |>
#   dplyr::ungroup()
# 
# 
# # Metrics
# metrics <- dplyr::bind_rows(
#   lapply(names(results_list), function(eco) {
#     
#     if (eco == "Global") return(NULL)
#     
#     # Local model
#     df_reg_raw <- results_list[[eco]]$df_cv |>
#       dplyr::filter(ecoreg == eco) |>
#       dplyr::mutate(
#         pred = exp(pred),
#         obs = exp(obs))
#     
#     # Global model
#     df_glob_raw <- results_list[["Global"]]$df_cv |>
#       dplyr::filter(ecoreg == eco) |>
#       dplyr::mutate(
#         pred = exp(pred),
#         obs = exp(obs))
#     
#     data.frame(
#       ecoreg = eco,
#       MAE_reg = round(mae(df_reg_raw$obs, df_reg_raw$pred), 1),
#       MAE_glob = round(mae(df_glob_raw$obs, df_glob_raw$pred), 1),
#       RHO_reg = round(cor(df_reg_raw$obs, df_reg_raw$pred, method = "spearman"), 2),
#       RHO_glob = round(cor(df_glob_raw$obs, df_glob_raw$pred, method = "spearman"), 2))
#   })
# )
# 
# metrics <- metrics |>
#   dplyr::mutate(
#     label = paste0(
#       "MAE = <span style='color:red;'>", MAE_reg,
#       "</span>, <span style='color:blue;'>", MAE_glob,
#       "</span><br>",
#       "ρ = <span style='color:red;'>", RHO_reg,
#       "</span>, <span style='color:blue;'>", RHO_glob,
#       "</span>"))
# 
# # Year separators
# year_breaks <- combined_df |>
#   dplyr::group_by(ecoreg, year) |>
#   dplyr::summarise(
#     x = min(time_id),
#     .groups = "drop")
# 
# # Plot
# p <- ggplot() +
#   # Observed bars
#   geom_col(data = combined_df,
#            aes(x = time_id, y = obs_mean, fill = n),
#            width = 1, alpha = 0.8, na.rm = TRUE) +
#   # # SE error bars - only for ver. 4
#   # geom_errorbar(data = combined_df, 
#   #               aes(x = time_id,
#   #                   ymin = obs_mean - se,
#   #                   ymax = obs_mean + se,
#   #                   color = "Observed ± SE"),
#   #               width = 0.15, 
#   #               linewidth = 0.3,
#   #               alpha = 0.8,
#   #               na.rm = TRUE) +
#   # Local model line
#   geom_line(data = combined_df,
#     aes(x = time_id, y = pred_reg_mean, color = "Local model"),
#     linewidth = 1, na.rm = TRUE) +
#   # Global model line
#   geom_line(data = combined_df,
#             aes(x = time_id, y = pred_glob_mean, color = "Global model"),
#             linewidth = 1, na.rm = TRUE) +
#   # Year separators
#   geom_vline(data = year_breaks, 
#              aes(xintercept = x),
#              linetype = "dashed",
#              color = "grey70") +
#   # Year labels
#   geom_text(data = year_breaks,
#             aes(x = x, y = -Inf, label = year),
#             inherit.aes = FALSE,
#             vjust = 1.25,
#             size = 3,
#             color = "grey30") +
#   # Metrics
#   ggtext::geom_richtext(data = metrics,
#                         aes(x = Inf, y = Inf, label = label),
#                         inherit.aes = FALSE,
#                         hjust = 1.1, vjust = 1.1,
#                         size = 3.5, fill = NA, label.color = NA) +
#   # Facet
#   facet_wrap(~ecoreg, ncol = 1, scales = "free") + # for ver. 1
#   #facet_wrap(~ecoreg, ncol = 1, scales = "fixed") + # for ver. 2-4
#   # Colors
#   scale_color_manual(values = c(#"Observed ± SE" = "black", # only for ver. 4
#     "Local model" = "red", "Global model" = "blue"),
#     name = NULL) +
#   scale_fill_gradient(low = "grey90",
#                       high = "grey30",
#                       name = "No. of traps") +
#   # X-axis
#   scale_x_continuous(expand = c(0.005, 0), breaks = NULL) +
#   # Labels
#   labs(title = "Observed vs. predicted abundance of Global-trained and Local-trained models across ecoregions",
#        #subtitle = "Observed abundance shown as bars; shade intensity indicates sampling intensity",
#        x = NULL, y = "Abundance") +
#   # Coords
#   coord_cartesian(clip = "off") + # for ver. 1 and 2
#   #coord_cartesian(ylim = c(0, 300), clip = "on") + # for ver. 3-4
#   # Theme
#   theme(panel.background = element_rect(fill = "white"),
#         panel.grid.major = element_line(color = "grey95"),
#         panel.grid.minor = element_blank(),
#         panel.border = element_rect(color = "black", fill = NA),
#         panel.spacing.y = unit(1.2, "lines"),
#         plot.margin = margin(10, 10, 20, 10),
#         strip.background = element_rect(fill = "grey90", color = "black"),
#         strip.text = element_text(face = "bold", size = 10),
#         legend.position = "bottom",
#         plot.title = element_text(face = "bold"),
#         axis.title = element_text(face = "bold"),
#         axis.text.x = element_blank(),
#         axis.text.y = element_text(angle = 90),
#         axis.ticks.x = element_blank())
# print(p)
# plot_model_list_full[["full_timeseries_v1"]] <- p # change version
# 
# ggsave(filename = file.path(outdir_supp, "full_timeseries_v1.png"), # change version
#        plot = p, width = 12, height = 8, dpi = 300)
# 
# ## Year-faceted point-line ----
# plot_model_list_year <- list()
# for (name in names(results_list)) {
#   
#   cat("\nPlotting:", name, "\n")
#   eco <- name
#   
#   # Skip Global plot
#   if (eco == "Global") next
#   
#   # Aggregate observed abundance for plotting
#   df_obs <- model_list[[eco]] |>
#     dplyr::group_by(year, week) |>
#     dplyr::summarise(
#       obs_mean = mean(count, na.rm = TRUE),
#       obs_sd = sd(count, na.rm = TRUE),
#       n = dplyr::n(),
#       .groups = "drop") |>
#     dplyr::mutate(se = obs_sd / sqrt(n)) # Compute SE
#   
#   # Local model (ID-level)
#   df_reg_raw <- results_list[[eco]]$df_cv |>
#     dplyr::filter(ecoreg == eco) |>
#     dplyr::mutate(
#       pred_reg_exp = exp(pred),
#       obs_exp = exp(obs))
#   
#   # Aggregate local model predictions
#   df_reg <- df_reg_raw |>
#     dplyr::group_by(year, week) |>
#     dplyr::summarise(
#       pred_reg_mean = mean(pred_reg_exp, na.rm = TRUE),
#       .groups = "drop")
#   
#   # Global model (filtered by ecoregion)
#   df_glob_raw <- results_list[["Global"]]$df_cv |>
#     dplyr::filter(ecoreg == eco) |>
#     dplyr::mutate(
#       pred_glob_exp = exp(pred),
#       obs_exp = exp(obs))
#   
#   # Aggregate global model predictions
#   df_glob <- df_glob_raw |>
#     dplyr::group_by(year, week) |>
#     dplyr::summarise(
#       pred_glob_mean = mean(pred_glob_exp, na.rm = TRUE),
#       .groups = "drop")
#   
#   # Merge for plotting
#   df_plot <- df_obs |>
#     dplyr::left_join(df_reg, by = c("year", "week")) |>
#     dplyr::left_join(df_glob, by = c("year", "week"))
#   
#   # Metrics
#   metrics <- df_reg_raw |>
#     dplyr::group_by(year) |>
#     dplyr::summarise(
#       MAE_reg = round(mae(obs_exp, pred_reg_exp), 1),
#       RHO_reg = round(cor(obs_exp, pred_reg_exp, method = "spearman"), 2),
#       .groups = "drop") |>
#     dplyr::left_join(df_glob_raw |>
#         dplyr::group_by(year) |>
#         dplyr::summarise(MAE_glob = round(mae(obs_exp, pred_glob_exp), 1),
#                          RHO_glob = round(cor(obs_exp, pred_glob_exp, method = "spearman"), 2),
#                          .groups = "drop"),
#         by = "year") |>
#     dplyr::mutate(label = paste0("MAE = <span style='color:red;'>", MAE_reg,
#                                  "</span>, <span style='color:blue;'>", MAE_glob,
#                                  "</span><br>",
#                                  "ρ = <span style='color:red;'>", RHO_reg,
#                                  "</span>, <span style='color:blue;'>", RHO_glob,
#                                  "</span>"))
#   
#   # Plot
#   p <- ggplot() +
#     # SE bars
#     geom_errorbar(data = df_plot,
#       aes(x = week, ymin = obs_mean - se, ymax = obs_mean + se),
#       width = 0.25,
#       color = "grey80",
#       linewidth = 0.5) +
#     # Observed line
#     geom_line(data = df_plot,
#               aes(x = week,
#                   y = obs_mean,
#                   color = "Observed abundance"),
#               linewidth = 1) +
#     # Observed points
#     geom_point(data = df_plot,
#                aes(x = week, y = obs_mean, size = n),
#                color = "grey80", alpha = 0.8) +
#     # Local model line
#     geom_line(data = df_plot,
#               aes(x = week, y = pred_reg_mean, color = "Predicted abundance (Local model)"),
#               linewidth = 1) +
#     # Global model line
#     geom_line(data = df_plot, aes(x = week, y = pred_glob_mean, color = "Predicted abundance (Global model)"),
#               linewidth = 1) +
#     facet_wrap(~year, ncol = 3, scales = "free_y") +
#     # Metrics
#     ggtext::geom_richtext(data = metrics,
#                           aes(x = Inf, y = Inf, label = label),
#                           inherit.aes = FALSE,
#                           hjust = 1.1, vjust = 1.1,
#                           size = 3.5, fill = NA,
#                           label.color = NA) +
#     scale_color_manual(values = c("Observed abundance" = "grey60", "Predicted abundance (Local model)" = "red",
#                                   "Predicted abundance (Global model)" = "blue"),
#                        name = NULL) +
#     scale_size(range = c(0.15, 6), name = "No. of traps") +
#     labs(title = paste("Observed vs. predicted abundance of Global-trained and Local-trained models in the",
#                        eco, "ecoregion per year"),
#          x = "Week", y = "Abundance") +
#     theme(panel.background = element_rect(fill = "white"),
#           panel.grid.major = element_line(color = "grey95"),
#           panel.border = element_rect(color = "black", fill = NA),
#           strip.background = element_rect(fill = "grey90", color = "black"),
#           strip.text = element_text(face = "bold"),
#           axis.text.y = element_text(angle = 90, hjust = 0),
#           legend.position = "bottom",
#           plot.title = element_text(face = "bold"))
#   # Store + print + save
#   plot_model_list_year[[name]] <- p
#   print(p)
#   
#   ggsave(filename = file.path(outdir_supp, paste0(name, "_model_prediction_year.png")),
#          plot = p, width = 13, height = 9, dpi = 300)
# }
# 
# ## Year-faceted ver 1 bar plot ----
# plot_model_list_v1 <- list()
# for (name in names(results_list)) {
#   
#   cat("\nPlotting:", name, "\n")
#   eco <- name
#   
#   # Skip Global plot
#   if (eco == "Global") next
#   
#   # Aggregate observed abundance for plotting
#   df_obs <- model_list[[eco]] |>
#     dplyr::group_by(year, week) |>
#     dplyr::summarise(obs_mean = mean(count, na.rm = TRUE),
#                      obs_sd   = sd(count, na.rm = TRUE),
#                      n        = dplyr::n(),
#                      .groups  = "drop") |> 
#     dplyr::mutate(se = obs_sd / sqrt(n)) # Calculate SE
#   
#   # Back transform regional model
#   df_reg_raw <- results_list[[eco]]$df_cv |>
#     dplyr::filter(ecoreg == eco) |>
#     dplyr::mutate(pred_exp = exp(pred),
#                   obs_exp  = exp(obs))
#   
#   # Aggregate regional predictions
#   df_reg <- df_reg_raw |>
#     dplyr::group_by(year, week) |>
#     dplyr::summarise(pred_reg_mean = mean(pred_exp, na.rm = TRUE),
#                      .groups = "drop")
#   
#   # Back transform global model
#   df_glob_raw <- results_list[["Global"]]$df_cv |>
#     dplyr::filter(ecoreg == eco) |>
#     dplyr::mutate(pred_exp = exp(pred),
#                   obs_exp  = exp(obs))
#   
#   # Aggregate global predictions
#   df_glob <- df_glob_raw |>
#     dplyr::group_by(year, week) |>
#     dplyr::summarise(pred_glob_mean = mean(pred_exp, na.rm = TRUE),
#                      .groups = "drop")
#   
#   # Merge for plotting
#   df_plot <- df_obs |>
#     dplyr::left_join(df_reg,  by = c("year", "week")) |>
#     dplyr::left_join(df_glob, by = c("year", "week"))
#   
#   # Long format
#   df_long <- df_plot |>
#     dplyr::select(year, week, pred_reg_mean, pred_glob_mean) |>
#     tidyr::pivot_longer(cols = c("pred_reg_mean", "pred_glob_mean"),
#                         names_to = "type",
#                         values_to = "value")
#   
#   # Metrics
#   metrics <- df_reg_raw |>
#     dplyr::group_by(year) |>
#     dplyr::summarise(MAE_reg  = round(mae(obs_exp, pred_exp), 1),
#                      RHO_reg  = round(cor(obs_exp, pred_exp, method = "spearman"), 2),
#                      .groups = "drop") |> 
#     dplyr::left_join(df_glob_raw |>
#                        dplyr::group_by(year) |>
#                        dplyr::summarise(MAE_glob = round(mae(obs_exp, pred_exp), 1),
#                                         RHO_glob = round(cor(obs_exp, pred_exp, method = "spearman"), 2),
#                                         .groups = "drop"),
#                      by = "year") |> 
#     dplyr::mutate(label = paste0("MAE = <span style='color:red;'>", MAE_reg, "</span>, ",
#                                  "<span style='color:blue;'>", MAE_glob, "</span><br>",
#                                  "ρ = <span style='color:red;'>", RHO_reg, "</span>, ",
#                                  "<span style='color:blue;'>", RHO_glob, "</span>"))
#   
#   # Plot
#   p <- ggplot() +
#     # Observed bars (intensity = n)
#     geom_col(data = df_plot,
#              aes(x = week, y = obs_mean, fill = n),
#              color = "grey50",
#              width = 0.8) +
#     # SE bars
#     geom_errorbar(data = df_plot,
#                   aes(x = week, ymin = obs_mean - se, ymax = obs_mean + se),
#                   width = 0.25,
#                   color = "grey40",
#                   linewidth = 0.5) +
#     # Model lines
#     geom_line(data = df_long,
#               aes(x = week, y = value, color = type),
#               linewidth = 1) +
#     facet_wrap(~year, ncol = 3, scales = "free_y") +
#     # Metrics
#     ggtext::geom_richtext(data = metrics,
#                           aes(x = Inf, y = Inf, label = label),
#                           inherit.aes = FALSE,
#                           hjust = 1.1, vjust = 1.1,
#                           size = 3.5, fill = NA, label.color = NA) +
#     # Fill scale (intensity)
#     scale_fill_gradient(low = "grey30", high = "grey95", name = "No. of traps") +
#     scale_color_manual(values = c("pred_reg_mean" = "red",
#                                   "pred_glob_mean" = "blue"),
#                        labels = c("pred_reg_mean" = paste0("Local (", eco, ") Model"),
#                                   "pred_glob_mean" = "Global Model")) +
#     labs(title = paste("Observed vs. predicted abundance of Global-trained and Local-trained models in the", eco, "ecoregion per year"),
#          subtitle = "Model predictors selected using Pearson correlation-based Filtering Approach; Observed abundance in bars",
#          x = "Week", y = "Abundance", color = "Model Training Approach") +
#     theme(panel.background = element_rect(fill = "white"),
#           panel.grid.major = element_line(color = "grey95"),
#           panel.border = element_rect(color = "black", fill = NA),
#           strip.background = element_rect(fill = "grey90", color = "black"),
#           strip.text = element_text(face = "bold"),
#           axis.text.y = element_text(angle = 90, hjust = 0),
#           axis.title = element_text(face = "bold"),
#           legend.position = "bottom",
#           plot.title = element_text(face = "bold"))
#   
#   # Store + print + save
#   plot_model_list_v1[[name]] <- p
#   print(p)
#   ggsave(filename = file.path(outdir_supp, paste0(name, "_model_prediction_v1_pc.png")),
#          plot = p, width = 13, height = 9, dpi = 300)
# }