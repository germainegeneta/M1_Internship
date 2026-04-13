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
library(lubridate)
library(scales)

# Set up
datasetName <- "VectAbundance_2024"
outdir_pred <- file.path(datasetName, "prediction_res")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Load biological data
data_glob <- read.csv(file.path(datasetName, "VectAbundance_2024_nonzero.csv")) |>
  transmute(ID = as.character(ID),
            ecoreg,
            date = as.Date(endDate),
            year, week, count, preceding_week,
            lat, long)
data_aedes <- read.csv(file.path(datasetName, "VectAbundance_2024_nonzero.csv")) |> 
  transmute(ID = as.character(ID),
            ecoreg,
            date = as.Date(endDate),
            year, week, count, preceding_week,
            lat, long)
data_aedes <- bind_rows(data_aedes, data_glob)

# Combine all aggregated climate dataframes
data_clim_raw <- bind_rows(read.csv(file.path(datasetName, "CCM_outputs", "fitmodels", "global_aggclim.csv")) |> 
                             mutate(ecoreg = paste0(ecoreg, "_g")),
                           read.csv(file.path(datasetName, "CCM_outputs", "fitmodels", "Alpine_aggclim.csv")),
                           read.csv(file.path(datasetName, "CCM_outputs", "fitmodels", "Continental_aggclim.csv")),
                           read.csv(file.path(datasetName, "CCM_outputs", "fitmodels", "Mediterranean_aggclim.csv")))

data_clim <- data_clim_raw |> 
  mutate(lag_id = paste0(lag_start, "_", lag_end)) |> 
  arrange(lag_start, lag_end) |> 
  rename(AV_TEMP = mean_temp_mean, 
         MX_TEMP = max_temp_mean,
         MN_TEMP = min_temp_mean,
         GDD = gdd_mean,
         SM_RAIN = sum_precip_sum,
         AV_RH = mean_rh_mean,
         PREC_ABUN = preceding_week) |> 
  pivot_wider(id_cols = c(ID, date, ecoreg, year, week, count, count_c, lat, long),
              names_from = lag_id,
              values_from = c(AV_TEMP, MX_TEMP, MN_TEMP,
                              GDD, SM_RAIN, AV_RH, PREC_ABUN),
              names_glue = "{.value}_{lag_id}") |> 
  relocate(ID, date, ecoreg, year, week, count, count_c, lat, long,
           starts_with("AV_TEMP"),
           starts_with("MX_TEMP"),
           starts_with("MN_TEMP"),
           starts_with("GDD"),
           starts_with("SM_RAIN"),
           starts_with("AV_RH"),
           starts_with("PREC_ABUN"))
# check climate data
stripchart(SM_RAIN_1_1 ~ ecoreg, data = data_clim, method = "jitter")

# save as .csv just in case
write_csv(data_clim, file.path(outdir_pred, "bioclim_data_withglob.csv"))

# Group and summarize at site-level
data_model <- data_clim |>
  select(-ID, -date, -count_c, -lat, -long) |> 
  group_by(ecoreg, year, week) |>
  summarise(across(count:PREC_ABUN_8_8, ~ mean(.x, na.rm = TRUE))) |>
  ungroup()

## ABUNDANCE MODEL PREPARATION ##
# Step 1: to select for meteorological variables, for every type of variable, the time lag for which the r2 was the highest. 
# Step 2: to evaluate the correlation between these variables
# Step 3:to select the variables not correlated with the highest sense ecological. 
# The first selection is crossed with the other selection done with the VIF with the corSelect function of the fuzzySim package 
# to select variables with the lowest VIF. The final selection is a mixed of both methods.


## STEP 1: select variables for abundance models ##
# global
glob_preds <- c("AV_TEMP_3_1", "MX_TEMP_4_1", "MN_TEMP_3_1", "GDD_3_1", "SM_RAIN_8_1", "AV_RH_6_1")
glob_model <- data_model |> 
  filter(str_detect(ecoreg, "_g")) |> 
  select("ecoreg", "year", "week", "count", all_of(glob_preds))

# alpine
alp_preds <- c("AV_TEMP_6_1", "MX_TEMP_6_1", "MN_TEMP_5_1", "GDD_6_1", "SM_RAIN_8_4", "AV_RH_6_4")
alp_model <- data_model |> 
  filter(ecoreg == "Alpine") |> 
  select("ecoreg", "year", "week", "count", all_of(alp_preds))

# continental
con_preds <- c("AV_TEMP_4_1", "MX_TEMP_4_1", "MN_TEMP_3_1", "GDD_4_1", "SM_RAIN_6_1", "AV_RH_5_1")
con_model <- data_model |> 
  filter(ecoreg == "Continental") |> 
  select("ecoreg", "year", "week", "count", all_of(con_preds))

# mediterranean
med_preds <- c("AV_TEMP_5_1", "MX_TEMP_5_1", "MN_TEMP_4_1", "GDD_4_1", "SM_RAIN_8_1", "AV_RH_8_1")
med_model <- data_model |> 
  filter(ecoreg == "Mediterranean") |> 
  select("ecoreg", "year", "week", "count", all_of(med_preds))

# Plot the bivariate relationship between abundance and each selected predictor
# predictors
pred_list <- list(Global = glob_preds, Alpine = alp_preds, Mediterranean = med_preds, Continental = con_preds)

# model data
model_list <- list(Global = glob_model, Alpine = alp_model, Mediterranean = med_model, Continental = con_model)

# colors (point + smooth)
color_list <- list(Global = c("gold", "goldenrod2"),
                   Alpine = c("lightcoral", "firebrick"),
                   Continental = c("palegreen3", "forestgreen"),
                   Mediterranean = c("skyblue", "steelblue"))

# Run loop
plot_list <- list()
for (eco in names(pred_list)) {
  preds <- pred_list[[eco]]
  cols <- color_list[[eco]]
  df <- model_list[[eco]]
  # subset data
  # df <- data_model |>
  #   #filter(ecoreg == eco) |>
  #   select(ecoreg, year, week, count, all_of(preds))
  
  # plot
  p <- df |>
    select(count, all_of(preds), year) |>
    pivot_longer(-c(count, year)) |>
    ggplot(aes(y = count, x = value)) +
    geom_point(color = cols[1]) +
    geom_smooth(color = cols[2], se = FALSE) +
    labs(
      title = paste("Bivariate relationship of abundance and selected climate variables in the", eco, "region"),
      x = "Value of climate variable",
      y = "Relative abundance (trap count)"
    ) +
    ylim(c(0, 80)) +
    facet_wrap(~name, scales = "free_x") +
    theme_bw() +
    theme(strip.background = element_rect(fill = "grey95", color = "grey40"),
          strip.text = element_text(face = "bold", size = 10),
          plot.title = element_text(face = "bold", size = 14))
  
  # store plot
  plot_list[[eco]] <- p
  print(p)
  
  # save
  # ggsave(filename = file.path(outdir_pred, paste0(eco, "_bivarplot.png")),
  #        plot = p, width = 12, height = 6, dpi = 300)
}


## STEP 2: identify correlated variables that are greater than 0.8 (Pearson correlation coefficient) ##
pear <- list()
for (eco in names(model_list)) {
  df    <- model_list[[eco]]
  preds <- pred_list[[eco]]
  
  # correlation matrix
  m <- cor(df[, preds], method = "pearson", use = "pairwise.complete.obs")
  
  # find high correlations
  index <- which(abs(m) > 0.8 & abs(m) < 1, arr.ind = TRUE)
  df_cor <- subset(as.data.frame(index), row <= col)
  
  # build result dataframe
  p <- cbind.data.frame(var1 = rownames(m)[df_cor$row],
                        var2 = colnames(m)[df_cor$col],
                        corr = m[cbind(df_cor$row, df_cor$col)])
  # store
  pear[[eco]] <- p
  
  # print
  cat("\n---", eco, "---\n")
  print(p)
}

# select final variables based on Signed R2 on CCM
glob_preds_pear <- c("MX_TEMP_4_1", "SM_RAIN_8_1", "AV_RH_6_1")
alp_preds_pear <- c("AV_TEMP_6_1", "SM_RAIN_8_4", "AV_RH_6_4")
con_preds_pear <- c("MN_TEMP_3_1", "SM_RAIN_6_1", "AV_RH_5_1")
med_preds_pear <- c("MX_TEMP_5_1", "SM_RAIN_8_1")


# Final data frame for the multivariate analysis
glob_model_pear <- glob_model |> select("ecoreg", "year", "week", "count", all_of(glob_preds_pear))
alp_model_pear <- alp_model |> select("ecoreg", "year", "week", "count", all_of(alp_preds_pear))
con_model_pear <- con_model |> select("ecoreg", "year", "week", "count", all_of(con_preds_pear))
med_model_pear <- med_model |> select("ecoreg", "year", "week", "count", all_of(med_preds_pear))


### STAGE 2: multivariate anaylsis using a leave-one-site-out cross validation and a 
# STEP 0: organize all datasets
model_list_mv <- list(Global_pear = glob_model_pear,
                      Alpine_pear = alp_model_pear,
                      Continental_pear = con_model_pear,
                      Mediterranean_pear = med_model_pear)

# STEP 0b: predictor lists (IMPORTANT)
pred_list_mv <- list(Global_pear = glob_preds_pear,
                     Alpine_pear = alp_preds_pear,
                     Continental_pear = con_preds_pear,
                     Mediterranean_pear = med_preds_pear)

# STEP 1: Spearman function
spearmcor <- function(data, lev = NULL, model = NULL) {
  out <- cor(x = data$pred, y = data$obs, method = "spearman")
  names(out) <- "spearman"
  out
}

# STEP 2: CV column
cv_col <- "year"

# STEP 3: loop
results_list <- list()
for (name in names(model_list_mv)) {
  cat("\nRunning model:", name, "\n")
  
  df <- model_list_mv[[name]]
  preds <- pred_list_mv[[name]]
  
  # log transform trap counts
  df$count <- log(df$count)
  
  # create CV folds
  indices_cv <- CAST::CreateSpacetimeFolds(df,
                                           spacevar = cv_col,
                                           k = length(unique(df[[cv_col]])))
  
  # trainControl
  tr <- trainControl(
    method = "cv",
    index = indices_cv$index,
    indexOut = indices_cv$indexOut,
    savePredictions = "final",
    summaryFunction = spearmcor)
 
  # train model
  mod <- caret::train(
    x = df[, preds],
    y = df$count,
    method = "ranger",
    tuneLength = 10,
    trControl = tr,
    metric = "spearman",
    maximize = TRUE,
    preProcess = c("center", "scale"),
    importance = "permutation")

  # STEP 5: attach predictions
  df$rowIndex <- seq_len(nrow(df))
  
  df_cv <- mod$pred |>
    left_join(df) |>
    select(pred, obs, ecoreg, week, year)
  
  # STEP 6: store results
  res <- list(model = mod,
              df_cv = df_cv,
              df_mod = df)
  results_list[[name]] <- res
  
  # STEP 7: save RDS
  saveRDS(res, file.path(outdir_pred, paste0(name, "_multiv_model_abundance_nowcasting2.rds")))
}

# OPTIONAL: assign all to environment
list2env(results_list, envir = .GlobalEnv)

# Now test Global-trained model on each ecoregion
# extract Global-trained model
mod_global <- results_list$Global_pear$model

# define regions
reg <- c("Alpine", "Continental", "Mediterranean")

# loop: predict
results_glob <- list()
for (eco in reg){
  cat("\nRunning GLOBAL model on:", eco, "\n")
  
  df <- data_model |> 
    filter(str_detect(ecoreg, eco))
  
  # log transform trap counts
  df$count <- log(df$count)
  
  # create CV folds
  indices_cv <- CAST::CreateSpacetimeFolds(df, spacevar = cv_col, k = length(unique(df[[cv_col]])))
  
  preds_all <- list()
  
  # 🔁 manual CV loop (important!)
  for (i in seq_along(indices_cv$index)) {
    
    train_idx <- indices_cv$index[[i]]
    test_idx  <- indices_cv$indexOut[[i]]
    
    df_test <- df[test_idx, ]
    
    # ensure predictors exist
    missing <- setdiff(glob_preds_pear, colnames(df_test))
    if (length(missing) > 0) {
      stop(paste("Missing predictors in", eco, ":", paste(missing, collapse = ", ")))
    }
    
    # predict using GLOBAL model
    pred <- predict(mod_global, newdata = df_test[, glob_preds_pear])
    
    # store fold predictions
    preds_all[[i]] <- data.frame(
      pred = pred,
      obs = df_test$count,
      ecoreg = eco,
      week = df_test$week,
      year = df_test$year
    )
  }
  
  # combine all folds
  df_cv <- bind_rows(preds_all)
  
  # results dataframe
  results_glob[[eco]] <- list(model = mod_global,
                              df_cv = df_cv)
}

## Plot them
plot_model_list <- list()

for (name in names(results_list)) {
  
  cat("\nPlotting:", name, "\n")
  
  parts <- strsplit(name, "_")[[1]]
  eco <- parts[1]
  method <- parts[2]
  
  if (eco == "Global") next
  
  # -----------------------------
  # REGIONAL
  # -----------------------------
  df_reg <- results_list[[name]]$df_cv |>
    mutate(
      obs_exp = exp(obs),
      pred_reg_exp = exp(pred)
    )
  
  # -----------------------------
  # GLOBAL
  # -----------------------------
  df_glob <- results_glob[[eco]]$df |>
    mutate(pred_glob_exp = exp(pred))
  
  # -----------------------------
  # MERGE
  # -----------------------------
  df_plot <- df_reg |>
    left_join(
      df_glob |> select(year, week, pred_glob_exp),
      by = c("year", "week")
    )
  
  # -----------------------------
  # LONG
  # -----------------------------
  df_long <- df_plot |>
    select(year, week, obs_exp, pred_reg_exp, pred_glob_exp) |>
    pivot_longer(
      cols = c("obs_exp", "pred_reg_exp", "pred_glob_exp"),
      names_to = "type",
      values_to = "value"
    )
  
  # -----------------------------
  # METRICS (WITH HTML COLORS)
  # -----------------------------
  metrics <- df_plot |>
    group_by(year) |>
    summarise(
      MAE_reg = round(mae(obs_exp, pred_reg_exp), 1),
      MAE_glob = round(mae(obs_exp, pred_glob_exp), 1),
      RHO_reg = round(cor(obs, pred, method = "spearman"), 2),
      RHO_glob = round(cor(obs, log(pred_glob_exp), method = "spearman"), 2),
      .groups = "drop"
    ) |>
    mutate(
      label = paste0(
        "MAE = <span style='color:red;'>", MAE_reg, "</span>, ",
        "<span style='color:blue;'>", MAE_glob, "</span><br>",
        "ρ = <span style='color:red;'>", RHO_reg, "</span>, ",
        "<span style='color:blue;'>", RHO_glob, "</span>"
      )
    )
  
  # -----------------------------
  # PLOT
  # -----------------------------
  p <- ggplot(df_long, aes(x = week, y = value, color = type)) +
    geom_line() +
    facet_wrap(~year, ncol = 3, scales = "free_y") +
    
    # 🔥 CLEAN COLORED LABEL
    ggtext::geom_richtext(
      data = metrics,
      aes(x = Inf, y = Inf, label = label),
      inherit.aes = FALSE,
      hjust = 1.1,
      vjust = 1.1,
      size = 3.5,
      fill = NA,
      label.color = NA
    ) +
    
    labs(
      title = paste("Observed vs Regional vs Global (", eco, ")", sep = ""),
      x = "Week",
      y = "Abundance",
      color = "Type"
    ) +
    
    scale_color_manual(
      values = c(
        "obs_exp" = "black",
        "pred_reg_exp" = "red",
        "pred_glob_exp" = "blue"
      ),
      labels = c(
        "obs_exp" = "Observed",
        "pred_reg_exp" = "Regional Model",
        "pred_glob_exp" = "Global Model"
      )
    ) +
    
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
  
  plot_model_list[[name]] <- p
  print(p)
  
  # ggsave(
  #   filename = file.path(outdir_pred, paste0(name, "_model_comparison.png")),
  #   plot = p,
  #   width = 10,
  #   height = 6,
  #   dpi = 300
  # )
}

## Plot full timeseries
# -----------------------------
# STEP 1: combine REGIONAL results
# -----------------------------
combined_list <- list()

for (name in names(results_list)) {
  
  df <- results_list[[name]]$df_cv
  
  parts <- strsplit(name, "_")[[1]]
  eco <- parts[1]
  
  if (eco == "Global") next
  
  df2 <- df |>
    mutate(
      ecoreg = eco,
      obs_exp = exp(obs),
      pred_reg_exp = exp(pred),
      week_date = ISOweek2date(paste0(year, "-W", sprintf("%02d", week), "-1"))
    )
  
  combined_list[[eco]] <- df2
}

combined_df <- bind_rows(combined_list)

# -----------------------------
# STEP 2: add GLOBAL predictions
# -----------------------------
glob_list <- list()

for (eco in names(results_glob)) {
  
  df <- results_glob[[eco]]$df
  
  df2 <- df |>
    mutate(
      ecoreg = eco,
      pred_glob_exp = exp(pred),
      week_date = ISOweek2date(paste0(year, "-W", sprintf("%02d", week), "-1"))
    ) |>
    select(ecoreg, year, week, week_date, pred_glob_exp)
  
  glob_list[[eco]] <- df2
}

glob_df <- bind_rows(glob_list)

# -----------------------------
# STEP 3: merge
# -----------------------------
combined_df <- combined_df |>
  left_join(glob_df,
            by = c("ecoreg", "year", "week", "week_date")) |>
  arrange(ecoreg, year, week)

# -----------------------------
# STEP 4: CREATE COMPRESSED TIME INDEX 🔥
# -----------------------------
combined_df <- combined_df |>
  group_by(ecoreg) |>
  arrange(year, week) |>
  mutate(time_id = row_number()) |>
  ungroup()

# -----------------------------
# STEP 5: METRICS
# -----------------------------
metrics <- combined_df |>
  group_by(ecoreg) |>
  summarise(
    
    MAE_reg = round(
      mae(
        obs_exp[!is.na(obs_exp) & !is.na(pred_reg_exp)],
        pred_reg_exp[!is.na(obs_exp) & !is.na(pred_reg_exp)]
      ), 1),
    
    MAE_glob = round(
      mae(
        obs_exp[!is.na(obs_exp) & !is.na(pred_glob_exp)],
        pred_glob_exp[!is.na(obs_exp) & !is.na(pred_glob_exp)]
      ), 1),
    
    RHO_reg = round(
      cor(obs_exp, pred_reg_exp,
          method = "spearman",
          use = "complete.obs"),
      2),
    
    RHO_glob = round(
      cor(obs_exp, pred_glob_exp,
          method = "spearman",
          use = "complete.obs"),
      2),
    
    .groups = "drop"
  ) |>
  mutate(
    label = paste0(
      "MAE = <span style='color:red;'>", MAE_reg, "</span>, ",
      "<span style='color:blue;'>", MAE_glob, "</span><br>",
      "ρ = <span style='color:red;'>", RHO_reg, "</span>, ",
      "<span style='color:blue;'>", RHO_glob, "</span>"
    )
  )

# -----------------------------
# STEP 6: YEAR POSITIONS (for dashed lines)
# -----------------------------
year_positions <- combined_df |>
  group_by(ecoreg, year) |>
  summarise(x = min(time_id), .groups = "drop")

# -----------------------------
# STEP 7: reshape
# -----------------------------
df_long <- combined_df |>
  select(ecoreg, time_id, obs_exp, pred_reg_exp, pred_glob_exp) |>
  pivot_longer(
    cols = c("obs_exp", "pred_reg_exp", "pred_glob_exp"),
    names_to = "type",
    values_to = "value"
  )

# -----------------------------
# STEP 8: plot
# -----------------------------
p <- ggplot(df_long,
            aes(x = time_id,
                y = value,
                color = type,
                group = type)) +
  
  geom_line(linewidth = 0.7, na.rm = TRUE) +   # 🔥 breaks at NA automatically
  
  geom_vline(data = year_positions,
             aes(xintercept = x),
             linetype = "dashed",
             linewidth = 0.3,
             color = "grey50",
             inherit.aes = FALSE) +
  
  facet_wrap(~ecoreg,
             ncol = 1,
             scales = "free_y") +
  
  ggtext::geom_richtext(
    data = metrics,
    aes(x = Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = 1.1,
    vjust = 1.1,
    size = 3.5,
    fill = NA,
    label.color = NA
  ) +
  
  scale_color_manual(
    values = c(
      "obs_exp" = "black",
      "pred_reg_exp" = "red",
      "pred_glob_exp" = "blue"
    ),
    labels = c(
      "obs_exp" = "Observed",
      "pred_reg_exp" = "Regional Model",
      "pred_glob_exp" = "Global Model"
    )
  ) +
  
  labs(
    title = "Observed vs Regional vs Global (Compressed Time Series)",
    x = "Time (compressed weeks)",
    y = "Abundance",
    color = "Type"
  ) +
  
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey85", color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

# -----------------------------
# STEP 9: save
# -----------------------------
ggsave(
  filename = file.path(outdir_pred, "Full_TS_compressed.png"),
  plot = p,
  width = 10,
  height = 8,
  dpi = 300
)

p

## VIP
results_list_sg <- results_list[-1]
vi_list <- list()
vip_list <- list()

# -----------------------------
# LOOP: extract VI + individual VIP plots
# -----------------------------
for (name in names(results_list_sg)) {
  
  cat("\nProcessing:", name, "\n")
  
  mod <- results_list_sg[[name]]$model
  
  eco <- strsplit(name, "_")[[1]][1]
  
  vi_scores <- vi(mod) |>
    mutate(ecoreg = eco)
  
  vi_list[[name]] <- vi_scores
  
  p <- vip(mod) +
    ggtitle(paste("Variable importance plot for", eco, "model"))
  
  vip_list[[name]] <- p
  
  print(p)
  
  ggsave(
    filename = file.path(outdir_pred, paste0(name, "_vip.png")),
    plot = p,
    width = 6,
    height = 5,
    dpi = 300)
}

# -----------------------------
# COMBINE ALL VI TABLES
# -----------------------------
vi_all <- bind_rows(vi_list)

# -----------------------------
# CREATE ORDER BASED ON PREFIX
# -----------------------------
var_order <- c("AV_RH", "SM_RAIN", "GDD", "MN_TEMP", "MX_TEMP", "AV_TEMP")

vi_all <- vi_all |>
  mutate(
    var_group = case_when(
      str_detect(Variable, "^AV_TEMP") ~ "AV_TEMP",
      str_detect(Variable, "^MX_TEMP") ~ "MX_TEMP",
      str_detect(Variable, "^MN_TEMP") ~ "MN_TEMP",
      str_detect(Variable, "^GDD")     ~ "GDD",
      str_detect(Variable, "^SM_RAIN") ~ "SM_RAIN",
      str_detect(Variable, "^AV_RH")   ~ "AV_RH",
      TRUE ~ "OTHER"
    ),
    
    var_group = factor(var_group, levels = var_order)
  )

# -----------------------------
# ORDER VARIABLES WITHIN GROUP
# -----------------------------
vi_all <- vi_all |>
  arrange(ecoreg, var_group, desc(Importance)) |>
  group_by(ecoreg) |>
  mutate(Variable = factor(Variable, levels = unique(Variable))) |>
  ungroup()

# -----------------------------
# FACETED VIP PLOT (ALL GREY)
# -----------------------------
p_all <- ggplot(vi_all,
                aes(x = Importance,
                    y = Variable)) +
  
  geom_col(fill = "grey40", width = 0.7) +
  
  facet_wrap(~ecoreg,
             ncol = 1,
             scales = "free_y") +
  
  labs(
    title = "Variable importance across ecoregion-specific models",
    x = "Importance",
    y = NULL
  ) +
  
  theme_bw() +
  theme(
    # ALL BLACK TEXT
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    strip.text = element_text(face = "bold", color = "black"),
    plot.title = element_text(face = "bold", color = "black"),
    
    # NO LEGEND
    legend.position = "none"
  )

p_all

# -----------------------------
# SAVE FINAL PLOT
# -----------------------------
ggsave(
  filename = file.path(outdir_pred, "VIP_all_ecoregions.png"),
  plot = p_all,
  width = 10,
  height = 8,
  dpi = 300
)

### PDP
pred_wrapper_reg <- function(object, newdata) {
  p <- predict(object, newdata = newdata)
  #c("avg" = mean(p), "sd" = sd(p))
  c("avg" = mean(p))
}

d <- results_list_sg$Continental_pear$model
tr <- results_list_sg$Continental_pear$df_mod[, pred_list_mv$Continental_pear]
t <- partial(d,
             pred.var = pred_list$Continental_pear,
             pred.fun = pred_wrapper_reg,
             train = tr)

###############################
# Plot model results
plot_model_list <- list()
for (name in names(results_list)) {
  
  cat("\nPlotting:", name, "\n")
  
  df_cv <- results_list[[name]]$df_cv
  
  # split name
  parts <- strsplit(name, "_")[[1]]
  eco <- parts[1]
  method <- parts[2]
  
  # back-transform + prepare
  df_plot <- df_cv |>
    mutate(obs_exp = exp(obs),
           pred_exp = exp(pred))
  
  # compute metrics per year
  metrics_df <- df_plot |>
    group_by(year) |>
    summarise(MAE = mae(obs_exp, pred_exp),
              SRHO = srho(obs, pred),
              .groups = "drop") |> 
    mutate(label = paste0("MAE = ", round(MAE, 1), "\nρ = ", round(SRHO, 2)))
  
  # main plot
  p <- df_plot |>
    pivot_longer(c("pred", "obs")) |>
    mutate(value = exp(value)) |>
    ggplot(aes(x = week, y = value, group = name, color = name)) +
    geom_line() +. 
    facet_wrap(~year, ncol = 4, scales = "free_y") +
    geom_text(data = metrics_df,
              aes(x = Inf, y = Inf, label = label),
              inherit.aes = FALSE,
              hjust = 1.1, 
              vjust = 1.1,
              size = 3.5) +
    labs(title = paste("Observed vs Predicted Abundance per Year\n",
                       eco, " region (", toupper(method), ")", sep = ""),
         x = "Week", y = "Abundance", color = "Type") +
    scale_color_manual(values = c("obs" = "black", "pred" = "red"),
                       labels = c("obs" = "Observed", "pred" = "Predicted")) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")
  # store
  plot_model_list[[name]] <- p
  
  # save
  ggsave(filename = file.path(outdir_pred, paste0(name, "_model_vis.png")),
         plot = p, width = 10, height = 6, dpi = 300)
}


# Plot the full timeseries
# STEP 1: combine all model outputs
combined_list <- list()

for (name in names(results_list)) {
  df <- results_list[[name]]$df_cv
  
  parts <- strsplit(name, "_")[[1]]
  eco <- parts[1]
  method <- parts[2]
  
  df2 <- df |>
    mutate(ecoreg = eco,
           method = toupper(method),   # PEAR / FFS
           obs_exp = exp(obs),
           pred_exp = exp(pred),
           date = ISOweek2date(paste0(year, "-W", sprintf("%02d", week), "-1")))
  combined_list[[name]] <- df2
}

combined_df <- bind_rows(combined_list)

# FUNCTION to plot one method
plot_full_timeseries <- function(method_name) {
  df_sub <- combined_df |>
    filter(method == method_name) |>
    mutate(ecoreg = factor(ecoreg,
                           levels = c("Global", "Alpine", "Continental", "Mediterranean")))
  
  # metrics per region
  metrics <- df_sub |>
    group_by(ecoreg) |>
    summarise(MAE = mae(obs_exp, pred_exp),
              SRHO = srho(obs, pred),
              .groups = "drop") |> 
    mutate(label = paste0("MAE = ", round(MAE, 1), "\nρ = ", round(SRHO, 2)))
  
  # year separators
  year_lines <- df_sub |>
    distinct(year) |>
    mutate(date = ISOweek2date(paste0(year, "-W01-1")))
  
  # plot
  p <- df_sub |>
    pivot_longer(c("pred", "obs")) |>
    mutate(value = exp(value)) |>
    ggplot(aes(x = date, y = value, color = name, group = name)) +
    geom_line() +
    geom_vline(data = year_lines,
               aes(xintercept = date),
               linetype = "dashed",
               color = "grey40",
               alpha = 0.6,
               inherit.aes = FALSE) +
    facet_wrap(~ecoreg, ncol = 1, scales = "free_y") +
    geom_text(data = metrics,
              aes(x = Inf, y = Inf, label = label),
              inherit.aes = FALSE,
              hjust = 1.1,
              vjust = 1.1,
              size = 3.5) +
    scale_x_date(date_breaks = "1 year",
                 date_labels = "%Y") +
    scale_color_manual(values = c("obs" = "black", "pred" = "red"),
                       labels = c("obs" = "Observed", "pred" = "Predicted")) +
    labs(title = paste("Observed vs Predicted Abundance (Full Time Series -", method_name, ")"),
         x = "Year", y = "Abundance", color = "Type") +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(face = "bold"))
  return(p)
}


# STEP 2: create both plots
p_pear <- plot_full_timeseries("PEAR")
p_pear

p_ffs  <- plot_full_timeseries("FFS")
p_ffs

# STEP 3: save
ggsave(file.path(outdir_pred, "Full_timeseries_FFS.png"),
       p_ffs, width = 10, height = 12, dpi = 300)


###### TRIAL: MAE FOR PARAMETER OPTIMIZATION #####

### STAGE 2: multivariate anaylsis using a leave-one-site-out cross validation and a 

# STEP 1: Spearman function
maefunc <- function(data, lev = NULL, model = NULL) {
  out <- mean(abs(data$pred - data$obs))
  names(out) <- "MAE"
  return(out)
}

# STEP 2: CV column
cv_col <- "year"

# STEP 3: loop
results_list_m <- list()
for (name in names(model_data_list)) {
  cat("\nRunning model:", name, "\n")
  
  df <- model_data_list[[name]]
  preds <- pred_list[[name]]
  
  # STEP 1: log transform
  df$count <- log(df$count)
  
  # STEP 2: create CV folds
  indices_cv <- CAST::CreateSpacetimeFolds(df,
                                           spacevar = cv_col,
                                           k = length(unique(df[[cv_col]])))
  
  # STEP 3: trainControl
  tr <- trainControl(
    method = "cv",
    index = indices_cv$index,
    indexOut = indices_cv$indexOut,
    savePredictions = "final",
    summaryFunction = maefunc)
  
  # STEP 4: train model
  mod <- caret::train(
    x = df[, preds],
    y = df$count,
    method = "ranger",
    tuneLength = 10,
    trControl = tr,
    metric = "MAE",
    maximize = FALSE,
    preProcess = c("center", "scale"),
    importance = "permutation")
  
  # STEP 5: attach predictions
  df$rowIndex <- seq_len(nrow(df))
  
  df_cv <- mod$pred |>
    left_join(df) |>
    select(pred, obs, ecoreg, week, year)
  
  # STEP 6: store results
  res <- list(
    model = mod,
    df_cv = df_cv,
    df_mod = df)
  
  results_list_m[[name]] <- res
  
  # STEP 7: save RDS
  saveRDS(res, file.path(outdir_pred, paste0(name, "_multiv_model_abundance_nowcasting2_mae.rds")))
}

# OPTIONAL: assign all to environment
list2env(results_list, envir = .GlobalEnv)

# Plot model results
plot_model_list_m <- list()
for (name in names(results_list_m)) {
  
  cat("\nPlotting:", name, "\n")
  
  df_cv <- results_list_m[[name]]$df_cv
  
  # split name
  parts <- strsplit(name, "_")[[1]]
  eco <- parts[1]
  method <- parts[2]
  
  # back-transform + prepare
  df_plot <- df_cv |>
    mutate(obs_exp = exp(obs),
           pred_exp = exp(pred))
  
  # compute metrics per year
  metrics_df <- df_plot |>
    group_by(year) |>
    summarise(MAE = mae(obs_exp, pred_exp),
              SRHO = srho(obs, pred),
              .groups = "drop") |> 
    mutate(label = paste0("MAE = ", round(MAE, 1), "\nρ = ", round(SRHO, 2)))
  
  # main plot
  p <- df_plot |>
    pivot_longer(c("pred", "obs")) |>
    mutate(value = exp(value)) |>
    ggplot(aes(x = week, y = value, group = name, color = name)) +
    geom_line() +
    facet_wrap(~year, ncol = 4, scales = "free_y") +
    geom_text(data = metrics_df,
              aes(x = Inf, y = Inf, label = label),
              inherit.aes = FALSE,
              hjust = 1.1, 
              vjust = 1.1,
              size = 3.5) +
    labs(title = paste("Observed vs Predicted Abundance per Year using MAE\n",
                       eco, " region (", toupper(method), ")", sep = ""),
         x = "Week", y = "Abundance", color = "Type") +
    scale_color_manual(values = c("obs" = "black", "pred" = "red"),
                       labels = c("obs" = "Observed", "pred" = "Predicted")) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")
  # store
  plot_model_list_m[[name]] <- p
  
  # save
  ggsave(filename = file.path(outdir_pred, paste0(name, "_model_vis_mae.png")),
         plot = p, width = 10, height = 6, dpi = 300)
}

# Plot the full timeseries
# STEP 1: combine all model outputs
combined_list_m <- list()

for (name in names(results_list_m)) {
  df <- results_list_m[[name]]$df_cv
  
  parts <- strsplit(name, "_")[[1]]
  eco <- parts[1]
  method <- parts[2]
  
  df2 <- df |>
    mutate(ecoreg = eco,
           method = toupper(method),   # PEAR / FFS
           obs_exp = exp(obs),
           pred_exp = exp(pred),
           date = ISOweek2date(paste0(year, "-W", sprintf("%02d", week), "-1")))
  combined_list_m[[name]] <- df2
}

combined_df_m <- bind_rows(combined_list_m)

plot_full_timeseries_m <- function(method_name) {
  df_sub <- combined_df_m |>
    filter(method == method_name) |>
    mutate(ecoreg = factor(ecoreg,
                           levels = c("Global", "Alpine", "Continental", "Mediterranean")))
  
  # metrics per region
  metrics <- df_sub |>
    group_by(ecoreg) |>
    summarise(MAE = mae(obs_exp, pred_exp),
              SRHO = srho(obs, pred),
              .groups = "drop") |> 
    mutate(label = paste0("MAE = ", round(MAE, 1), "\nρ = ", round(SRHO, 2)))
  
  # year separators
  year_lines <- df_sub |>
    distinct(year) |>
    mutate(date = ISOweek2date(paste0(year, "-W01-1")))
  
  # plot
  p <- df_sub |>
    pivot_longer(c("pred", "obs")) |>
    mutate(value = exp(value)) |>
    ggplot(aes(x = date, y = value, color = name, group = name)) +
    geom_line() +
    geom_vline(data = year_lines,
               aes(xintercept = date),
               linetype = "dashed",
               color = "grey40",
               alpha = 0.6,
               inherit.aes = FALSE) +
    facet_wrap(~ecoreg, ncol = 1, scales = "free_y") +
    geom_text(data = metrics,
              aes(x = Inf, y = Inf, label = label),
              inherit.aes = FALSE,
              hjust = 1.1,
              vjust = 1.1,
              size = 3.5) +
    scale_x_date(date_breaks = "1 year",
                 date_labels = "%Y") +
    scale_color_manual(values = c("obs" = "black", "pred" = "red"),
                       labels = c("obs" = "Observed", "pred" = "Predicted")) +
    labs(title = paste("Observed vs Predicted Abundance MAE (Full Time Series -", method_name, ")"),
         x = "Year", y = "Abundance", color = "Type") +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(face = "bold"))
  return(p)
}


# STEP 2: create both plots
p_pear_m <- plot_full_timeseries_m("PEAR")
p_pear_m

p_ffs_m  <- plot_full_timeseries_m("FFS")
p_ffs_m

# STEP 3: save
ggsave(file.path(outdir_pred, "Full_timeseries_FFS_mae.png"),
       p_ffs_m, width = 10, height = 12, dpi = 300)


########################### VIP ############################

df <- results_list$Global_ffs$model
vi_scores <- vi(df)
print(vi_scores)
vip(df)

vi_list_m <- list()
vip_list_m <- list()
for (name in names(results_list_m)) {
  
  cat("\nProcessing:", name, "\n")
  
  mod <- results_list_m[[name]]$model
  
  # get VI table
  vi_scores <- vi(mod)
  print(vi_scores)
  
  # store VI table
  vi_list_m[[name]] <- vi_scores
  
  # create plot
  p <- vip(mod) + ggtitle(name)
  
  # store plot
  vip_list_m[[name]] <- p
  
  # optionally print
  print(p)
  
  # save plot
  # ggsave(
  #   filename = file.path(outdir_pred, paste0(name, "_vip.png")),
  #   plot = p,
  #   width = 6,
  #   height = 5,
  #   dpi = 300
  # )
}
