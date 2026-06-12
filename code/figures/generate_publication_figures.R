# ==============================================================================
# GENERATE PUBLICATION FIGURES, TABLES, AND LEGEND DOCUMENT
# For: ASF ENM Research Paper
# Output: outputs/publication_figures/
# ==============================================================================

library(officer)
library(terra)
library(geodata)
library(sf)
library(ggplot2)
library(tidyterra)
library(scales)
library(cowplot)
library(dplyr)
library(viridis)
library(caret)
library(ggspatial)   # scale bar + north arrow
library(xgboost)     # for XGBoost variable importance regeneration
library(tidyr)       # for pivot_longer in Figure S4

setwd("C:/Users/s9240/OneDrive/Desktop/ENM project/ASF_ENM_Asia")

pub_dir <- "outputs/publication_figures"
dir.create(pub_dir, showWarnings = FALSE, recursive = TRUE)

# Load common data
env_stack  <- rast("data_processed/env_stack_asia_FINAL.tif")
e_asia     <- ext(90, 150, -10, 50)
land_asia  <- world(resolution = 3, path = "data_raw") |> crop(e_asia)
asf_clean  <- read.csv("outputs/ASF_cleaned_points.csv")
perf_table <- read.csv("outputs/Model_Performance_Summary.csv")

# Helper: clean display names (replace underscores with spaces)
clean_name <- function(x) gsub("_", " ", x)


# ==============================================================================
# FIGURE 1: Study Area & ASF Outbreak Distribution
# ==============================================================================
cat("Generating Figure 1: Study Area...\n")

# Country label positions for key nations
country_labels <- data.frame(
  label = c("China", "Vietnam", "Thailand", "Philippines",
            "Indonesia", "Myanmar", "Cambodia", "Laos",
            "Malaysia", "South Korea", "Japan"),
  lon   = c(105, 106, 101, 122,
            115, 97, 105, 103,
            109, 128, 138),
  lat   = c(35, 16, 15, 12,
            -2, 20, 12.5, 19,
            4, 36, 37)
)

p_fig1 <- ggplot() +
  geom_spatvector(data = land_asia, fill = "#f5f5f5", color = "grey40",
                  linewidth = 0.3) +
  geom_point(data = asf_clean, aes(x = Longitude, y = Latitude),
             color = "#bd0026", fill = "#ef4444", shape = 21,
             size = 1.0, stroke = 0.3, alpha = 0.5) +
  geom_text(data = country_labels, aes(x = lon, y = lat, label = label),
            size = 2.8, color = "grey25", fontface = "italic") +
  annotation_scale(location = "bl", width_hint = 0.2,
                   pad_x = unit(0.5, "cm"), pad_y = unit(0.5, "cm"),
                   style = "ticks") +
  annotation_north_arrow(location = "tr", which_north = "true",
                         style = north_arrow_fancy_orienteering(),
                         height = unit(1.2, "cm"), width = unit(1.2, "cm"),
                         pad_x = unit(0.3, "cm"), pad_y = unit(0.3, "cm")) +
  coord_sf(xlim = c(90, 150), ylim = c(-10, 50), expand = FALSE) +
  labs(
    title = "Distribution of Confirmed ASF Outbreaks in East & Southeast Asia",
    subtitle = paste0(nrow(asf_clean),
                      " unique outbreak locations (WOAH EMPRES-i, 2015\u20132025)"),
    x = "Longitude (\u00b0E)", y = "Latitude (\u00b0N)",
    caption = "Source: FAO EMPRES-i | Deduplicated to 2.5 arc-min grid"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, color = "grey30"),
    plot.caption  = element_text(size = 7, color = "grey50", hjust = 0),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2)
  )

ggsave(file.path(pub_dir, "Figure_1_Study_Area_Outbreaks.png"),
       p_fig1, width = 12, height = 8, dpi = 300)


# ==============================================================================
# FIGURE 2A: Correlation Matrix (copy existing — already ggcorrplot)
# ==============================================================================
cat("Copying Figure 2A: Correlation Matrix...\n")
file.copy("outputs/Correlation_Matrix.png",
          file.path(pub_dir, "Figure_2A_Correlation_Matrix.png"), overwrite = TRUE)


# ==============================================================================
# FIGURE 2B: Boruta Feature Selection — REGENERATED in ggplot2
# ==============================================================================
cat("Generating Figure 2B: Boruta Selection (ggplot2)...\n")

# Load the final data to reconstruct Boruta-style importance
data_final <- readRDS("data_processed/asf_final_data.rds")

# Run a quick Boruta to extract importance stats, or reconstruct from the
# selected variables list (the original Boruta object was not saved to disk).
# Strategy: train a quick RF and compute Z-scores to replicate the Boruta display.
if (requireNamespace("Boruta", quietly = TRUE)) {
  library(Boruta)
  set.seed(124)
  data_boruta <- data_final
  data_boruta$Class <- factor(make.names(data_boruta$Class))
  boruta_rerun <- Boruta(Class ~ ., data = data_boruta,
                          doTrace = 0, maxRuns = 100)

  # Extract importance statistics
  boruta_imp <- as.data.frame(boruta_rerun$ImpHistory)

  # Reshape for ggplot
  boruta_long <- boruta_imp |>
    pivot_longer(cols = everything(), names_to = "Variable", values_to = "Importance") |>
    filter(!is.na(Importance))

  # Get decision for each variable
  decisions <- data.frame(
    Variable = names(boruta_rerun$finalDecision),
    Decision = as.character(boruta_rerun$finalDecision),
    stringsAsFactors = FALSE
  )
  # Add shadow variables
  shadow_vars <- grep("^shadow", names(boruta_imp), value = TRUE)
  shadow_decisions <- data.frame(
    Variable = shadow_vars,
    Decision = "Shadow",
    stringsAsFactors = FALSE
  )
  decisions <- rbind(decisions, shadow_decisions)

  boruta_long <- boruta_long |>
    left_join(decisions, by = "Variable") |>
    mutate(Decision = ifelse(is.na(Decision), "Shadow", Decision))

  # Order by median importance
  var_order <- boruta_long |>
    group_by(Variable) |>
    summarise(med = median(Importance, na.rm = TRUE)) |>
    arrange(med) |>
    pull(Variable)
  boruta_long$Variable <- factor(boruta_long$Variable, levels = var_order)

  # Clean display names
  boruta_long$DisplayName <- factor(clean_name(as.character(boruta_long$Variable)),
                                     levels = clean_name(var_order))

  p_fig2b <- ggplot(boruta_long, aes(x = DisplayName, y = Importance, fill = Decision)) +
    geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.5, linewidth = 0.3) +
    scale_fill_manual(values = c("Confirmed" = "#2ca02c",
                                  "Rejected"  = "#d62728",
                                  "Tentative" = "#ff7f0e",
                                  "Shadow"    = "#1f77b4"),
                      name = "Boruta Decision") +
    labs(
      title = "Boruta Feature Importance \u2014 ASF Asia ENM",
      subtitle = "maxRuns = 100 | Green = Confirmed, Red = Rejected, Blue = Shadow reference",
      x = NULL, y = "Z-Score Importance"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, color = "grey30"),
      axis.text.x   = element_text(angle = 45, hjust = 1, size = 9),
      legend.position = "top"
    )

  ggsave(file.path(pub_dir, "Figure_2B_Boruta_Selection.png"),
         p_fig2b, width = 12, height = 7, dpi = 300)
} else {
  cat("  Boruta package not available. Copying original plot instead.\n")
  file.copy("outputs/Boruta_Feature_Selection.png",
            file.path(pub_dir, "Figure_2B_Boruta_Selection.png"), overwrite = TRUE)
}


# ==============================================================================
# FIGURE 3: Wild Boar Density Map (copy existing)
# ==============================================================================
cat("Copying Figure 3: Wild Boar Density...\n")
file.copy("outputs/Wild_Boar_Density_Map.png",
          file.path(pub_dir, "Figure_3_Wild_Boar_Density.png"), overwrite = TRUE)


# ==============================================================================
# FIGURE 4: Model Performance — Random CV vs Spatial Block CV (FIXED)
# ==============================================================================
cat("Generating Figure 4: Spatial CV Comparison...\n")

# Load caret models to extract standard CV SD from resampling results
rf_fit  <- readRDS("data_processed/model_RF.rds")
svm_fit <- readRDS("data_processed/model_SVM.rds")
glm_fit <- readRDS("data_processed/model_GLM.rds")

# Extract standard CV SD from caret resampling results
# For caret models, use the SD of ROC across resampling folds
cv_results <- resamples(list(
  "Random Forest"    = rf_fit,
  "SVM"              = svm_fit,
  "GLM"              = glm_fit
))
cv_summary <- summary(cv_results)
# SD of ROC across 30 resamples (10-fold x 3 repeats)
std_cv_sds <- cv_summary$statistics$ROC[, "SD"]

# Build long-format data with SDs for BOTH CV methods
perf_long <- data.frame(
  Model  = rep(perf_table$Model, 2),
  Metric = rep(c("Standard Random CV", "Spatial Block CV"),
               each = nrow(perf_table)),
  AUC    = c(perf_table$CV_AUC, perf_table$Spatial_CV_AUC),
  SD     = c(
    # Standard CV SDs: match from caret resamples; XGBoost gets NA (native API)
    std_cv_sds["Random Forest"],
    NA,  # XGBoost — native API, no caret resamples available
    std_cv_sds["SVM"],
    std_cv_sds["GLM"],
    # Spatial CV SDs
    perf_table$Spatial_CV_SD
  )
)
perf_long$Model  <- factor(perf_long$Model,
  levels = c("Random Forest", "XGBoost (native)", "SVM", "GLM"))
perf_long$Metric <- factor(perf_long$Metric,
  levels = c("Standard Random CV", "Spatial Block CV"))

p_fig4 <- ggplot(perf_long, aes(x = Model, y = AUC, fill = Metric)) +
  geom_col(position = position_dodge(0.7), width = 0.6,
           color = "grey30", linewidth = 0.3) +
  geom_errorbar(aes(ymin = AUC - ifelse(is.na(SD), 0, SD),
                     ymax = AUC + ifelse(is.na(SD), 0, SD)),
                position = position_dodge(0.7), width = 0.15, linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.3f", AUC)),
            position = position_dodge(0.7), vjust = -0.8,
            size = 3, fontface = "bold") +
  scale_fill_manual(values = c("Standard Random CV" = "#4393c3",
                                "Spatial Block CV"   = "#d6604d")) +
  scale_y_continuous(limits = c(0.5, 1.05), breaks = seq(0.5, 1.0, 0.1),
                     expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(clip = "off") +
  labs(
    title    = "Model Performance: Standard CV vs. Spatial Block CV",
    subtitle = "Spatial block CV (500 km hexagonal blocks, k = 5) tests generalization to unseen geographic regions",
    y = "Area Under the ROC Curve (AUC)", x = NULL,
    fill    = "Validation Method",
    caption = "Error bars = \u00b1 1 SD across folds | XGBoost uses native API (no caret resamples for std CV) | Ref: Valavi et al. (2019)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, color = "grey30"),
    plot.caption  = element_text(size = 7, color = "grey50", hjust = 0),
    legend.position = "top",
    axis.text.x   = element_text(size = 10, face = "bold"),
    plot.margin   = margin(10, 15, 10, 10)
  )

ggsave(file.path(pub_dir, "Figure_4_Spatial_CV_Comparison.png"),
       p_fig4, width = 10, height = 7, dpi = 300)


# ==============================================================================
# FIGURE 5–8: Copy existing high-quality figures
# ==============================================================================
cat("Copying Figures 5-8...\n")
file.copy("outputs/ASF_RiskMaps_All4Models.png",
          file.path(pub_dir, "Figure_5_Individual_Model_Maps.png"), overwrite = TRUE)
file.copy("outputs/ASF_RiskMap_ENSEMBLE_FINAL.png",
          file.path(pub_dir, "Figure_6_Ensemble_Risk_Map.png"), overwrite = TRUE)
file.copy("outputs/PDP_Top6_Variables.png",
          file.path(pub_dir, "Figure_7_PDP_Top6_Variables.png"), overwrite = TRUE)
file.copy("outputs/ICE_Pig_Density.png",
          file.path(pub_dir, "Figure_8_ICE_Pig_Density.png"), overwrite = TRUE)


# ==============================================================================
# SUPPLEMENTARY FIGURES
# ==============================================================================
cat("Generating Supplementary Figures...\n")

# --- Figure S1: Soft Tick Occurrences (copy existing) ---
file.copy("outputs/Soft_Tick_Occurrences_Map.png",
          file.path(pub_dir, "Figure_S1_Soft_Tick_Occurrences.png"), overwrite = TRUE)


# --- Figure S2: RF Variable Importance — REGENERATED in ggplot2 ---
cat("  Regenerating Figure S2: RF Variable Importance (ggplot2)...\n")

rf_imp   <- varImp(rf_fit, scale = TRUE)
# Safely extract importance column (works for both caret naming conventions)
imp_vals <- if ("Positive" %in% names(rf_imp$importance)) {
  rf_imp$importance$Positive
} else if ("Overall" %in% names(rf_imp$importance)) {
  rf_imp$importance$Overall
} else {
  rf_imp$importance[, 1]
}

rf_imp_df <- data.frame(
  Variable   = rownames(rf_imp$importance),
  Importance = imp_vals,
  stringsAsFactors = FALSE
) |>
  mutate(DisplayName = clean_name(Variable)) |>
  arrange(Importance) |>
  mutate(DisplayName = factor(DisplayName, levels = DisplayName))

p_s2 <- ggplot(rf_imp_df, aes(x = DisplayName, y = Importance)) +
  geom_segment(aes(xend = DisplayName, y = 0, yend = Importance),
               color = "#4393c3", linewidth = 0.8) +
  geom_point(color = "#2166ac", size = 3) +
  coord_flip() +
  labs(
    title = "Random Forest Variable Importance",
    subtitle = "Scaled importance (0\u2013100) | ASF Asia ENM",
    x = NULL, y = "Scaled Importance"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, color = "grey30"),
    panel.grid.major.y = element_blank()
  )

ggsave(file.path(pub_dir, "Figure_S2_RF_Variable_Importance.png"),
       p_s2, width = 8, height = 6, dpi = 300)


# --- Figure S3: XGBoost Variable Importance — REGENERATED in ggplot2 ---
cat("  Regenerating Figure S3: XGBoost Variable Importance (ggplot2)...\n")

xgb_fit  <- readRDS("data_processed/model_XGBoost.rds")
xgb_imp_mat <- xgb.importance(model = xgb_fit)

xgb_imp_df <- data.frame(
  Variable   = xgb_imp_mat$Feature,
  Gain       = xgb_imp_mat$Gain,
  stringsAsFactors = FALSE
) |>
  mutate(DisplayName = clean_name(Variable)) |>
  arrange(Gain) |>
  mutate(DisplayName = factor(DisplayName, levels = DisplayName))

p_s3 <- ggplot(xgb_imp_df, aes(x = DisplayName, y = Gain)) +
  geom_col(fill = "#d6604d", color = "grey30", linewidth = 0.3, width = 0.7) +
  coord_flip() +
  labs(
    title = "XGBoost Variable Importance (Gain)",
    subtitle = "Total gain contribution | ASF Asia ENM",
    x = NULL, y = "Gain"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, color = "grey30"),
    panel.grid.major.y = element_blank()
  )

ggsave(file.path(pub_dir, "Figure_S3_XGBoost_Variable_Importance.png"),
       p_s3, width = 8, height = 6, dpi = 300)


# --- Figure S4: CV AUC Distribution — REGENERATED in ggplot2 ---
cat("  Regenerating Figure S4: CV AUC Boxplot (ggplot2)...\n")

# Extract per-fold AUC values from caret resamples
cv_vals <- cv_results$values
# Columns are named "ModelName~ROC", "ModelName~Sens", "ModelName~Spec"
roc_cols <- grep("~ROC$", names(cv_vals), value = TRUE)
cv_auc_df <- cv_vals[, roc_cols, drop = FALSE]
names(cv_auc_df) <- gsub("~ROC$", "", names(cv_auc_df))

# Reshape to long format
cv_auc_long <- cv_auc_df |>
  pivot_longer(cols = everything(),
               names_to = "Model", values_to = "AUC")

# Ensure consistent naming and ordering
cv_auc_long$Model <- factor(cv_auc_long$Model,
  levels = c("Random Forest", "SVM", "GLM"))

p_s4 <- ggplot(cv_auc_long, aes(x = Model, y = AUC, fill = Model)) +
  geom_boxplot(width = 0.5, outlier.shape = 21, outlier.size = 2,
               color = "grey30", linewidth = 0.4) +
  geom_jitter(width = 0.12, alpha = 0.3, size = 1.5, color = "grey40") +
  scale_fill_manual(values = c("Random Forest" = "#4393c3",
                                "SVM"           = "#f4a582",
                                "GLM"           = "#d6604d")) +
  scale_y_continuous(limits = c(min(cv_auc_long$AUC, na.rm = TRUE) - 0.02, 1.0),
                     breaks = seq(0.8, 1.0, 0.02)) +
  labs(
    title = "Cross-Validation AUC Distribution (10-Fold \u00d7 3 Repeats)",
    subtitle = paste0("RF, SVM, and GLM via caret (n = 30 resamples each) | ",
                      "XGBoost omitted (native API with separate 10-fold CV)"),
    x = NULL, y = "Area Under the ROC Curve (AUC)",
    caption = "Individual resample values shown as jittered points"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, color = "grey30"),
    plot.caption  = element_text(size = 7, color = "grey50", hjust = 0),
    legend.position = "none",
    axis.text.x   = element_text(size = 11, face = "bold")
  )

ggsave(file.path(pub_dir, "Figure_S4_CV_AUC_Boxplot.png"),
       p_s4, width = 8, height = 6, dpi = 300)


# ==============================================================================
# TABLE 1: Environmental Covariates
# ==============================================================================
cat("Generating Table 1: Covariates...\n")

boruta_selected <- c("BIO13", "BIO14", "Elevation", "Wind_Speed",
                      "Population_Density", "Pig_Density", "Chicken_Density",
                      "Soil_Clay_Content", "Tree_Cover_Fraction",
                      "Cropland_Fraction", "Grassland_Fraction",
                      "Wild_Boar_Density")

table1 <- data.frame(
  Variable = c(
    "BIO1","BIO2","BIO3","BIO4","BIO5","BIO6","BIO7",
    "BIO10","BIO11","BIO12","BIO13","BIO14","BIO15","BIO16","BIO17",
    "Elevation","Wind_Speed","Population_Density",
    "Pig_Density","Chicken_Density",
    "Soil_pH","Soil_Clay_Content",
    "Tree_Cover_Fraction","Cropland_Fraction","Grassland_Fraction",
    "Wild_Boar_Density","Soft_Tick_Density"),
  Description = c(
    "Annual Mean Temperature","Mean Diurnal Range","Isothermality",
    "Temperature Seasonality","Max Temp of Warmest Month",
    "Min Temp of Coldest Month","Temperature Annual Range",
    "Mean Temp of Warmest Quarter","Mean Temp of Coldest Quarter",
    "Annual Precipitation","Precipitation of Wettest Month",
    "Precipitation of Driest Month","Precipitation Seasonality",
    "Precipitation of Wettest Quarter","Precipitation of Driest Quarter",
    "Elevation (m a.s.l.)","Mean Wind Speed (m/s)",
    "Human Population Density (persons/km2)",
    "Domestic Pig Density (heads/km2)","Chicken Density (heads/km2)",
    "Soil pH (H2O, 0-5 cm)","Soil Clay Content (%, 0-5 cm)",
    "Tree Cover Fraction (%)","Cropland Fraction (%)","Grassland Fraction (%)",
    "Wild Boar KDE Density (GBIF)","Soft Tick KDE Density (GBIF)"),
  Source = c(rep("WorldClim 2.1", 15),
             "SRTM via WorldClim","WorldClim 2.1","WorldPop 2020",
             "FAO GLW4","FAO GLW4","SoilGrids","SoilGrids",
             "ESA WorldCover","ESA WorldCover","ESA WorldCover",
             "GBIF + KDE","GBIF + KDE"),
  Resolution = c(rep("2.5 arc-min", 18),
                 "10 km (resampled)","10 km (resampled)",
                 "250 m (resampled)","250 m (resampled)",
                 "100 m (resampled)","100 m (resampled)","100 m (resampled)",
                 "2.5 arc-min (KDE)","2.5 arc-min (KDE)"),
  Final_Status = ifelse(
    c("BIO1","BIO2","BIO3","BIO4","BIO5","BIO6","BIO7",
      "BIO10","BIO11","BIO12","BIO13","BIO14","BIO15","BIO16","BIO17",
      "Elevation","Wind_Speed","Population_Density",
      "Pig_Density","Chicken_Density",
      "Soil_pH","Soil_Clay_Content",
      "Tree_Cover_Fraction","Cropland_Fraction","Grassland_Fraction",
      "Wild_Boar_Density","Soft_Tick_Density") %in% boruta_selected,
    "Selected", "Excluded"),
  stringsAsFactors = FALSE
)

write.csv(table1, file.path(pub_dir, "Table_1_Covariates.csv"), row.names = FALSE)


# ==============================================================================
# TABLE 2: Model Performance (copy existing, already has spatial CV)
# ==============================================================================
cat("Copying Table 2: Model Performance...\n")
file.copy("outputs/Model_Performance_Summary.csv",
          file.path(pub_dir, "Table_2_Model_Performance.csv"), overwrite = TRUE)


# ==============================================================================
# TABLE 3: Variable Importance Rankings (FIXED column access)
# ==============================================================================
cat("Generating Table 3: Variable Importance...\n")

# Safely extract importance — handle both "Positive" and "Overall" column names
rf_imp_t3 <- varImp(rf_fit, scale = TRUE)
imp_col <- if ("Positive" %in% names(rf_imp_t3$importance)) {
  rf_imp_t3$importance$Positive
} else if ("Overall" %in% names(rf_imp_t3$importance)) {
  rf_imp_t3$importance$Overall
} else {
  rf_imp_t3$importance[, 1]
}

imp_ord <- order(imp_col, decreasing = TRUE)
imp_df  <- data.frame(
  Rank       = seq_along(imp_ord),
  Variable   = rownames(rf_imp_t3$importance)[imp_ord],
  RF_Importance_Scaled = round(imp_col[imp_ord], 2)
)
write.csv(imp_df, file.path(pub_dir, "Table_3_Variable_Importance.csv"),
          row.names = FALSE)


# ==============================================================================
# WORD DOCUMENT: Figure & Table Legend Document
# ==============================================================================
cat("Generating Word document with all legends...\n")

doc <- read_docx()

# ------ TITLE ------
doc <- body_add_par(doc,
  "Publication Figures & Tables \u2014 Legend Document", style = "heading 1")
doc <- body_add_par(doc,
  "African Swine Fever Ecological Niche Modeling in East & Southeast Asia",
  style = "heading 2")
doc <- body_add_par(doc, paste("Generated:", Sys.Date()), style = "Normal")
doc <- body_add_par(doc, "", style = "Normal")


# ------ MAIN TEXT FIGURES ------
doc <- body_add_par(doc, "Main Text Figures", style = "heading 1")

# Helper: safely add image
add_fig <- function(doc, path, w = 6, h = 4) {
  if (file.exists(path)) doc <- body_add_img(doc, src = path, width = w, height = h)
  return(doc)
}

# Figure 1
doc <- body_add_par(doc,
  "Figure 1. Study Area and ASF Outbreak Distribution", style = "heading 2")
doc <- add_fig(doc, file.path(pub_dir, "Figure_1_Study_Area_Outbreaks.png"), 6, 4)
doc <- body_add_par(doc, paste0(
  "Figure 1. Spatial distribution of ", nrow(asf_clean),
  " confirmed African Swine Fever (ASF) outbreak locations across East and ",
  "Southeast Asia (90-150 E, 10 S - 50 N), sourced from the FAO EMPRES-i/WOAH ",
  "database (2015-2025). Points were spatially deduplicated to one unique ",
  "location per 2.5 arc-minute (~4.5 km) grid cell to reduce spatial ",
  "pseudoreplication. The study extent encompasses major swine-producing ",
  "nations spanning tropical, subtropical, and temperate climate zones. ",
  "Country labels indicate key nations within the study region."),
  style = "Normal")
doc <- body_add_par(doc, "", style = "Normal")

# Figure 2
doc <- body_add_par(doc,
  "Figure 2. Feature Selection: Correlation Filtering and Boruta Analysis",
  style = "heading 2")
doc <- add_fig(doc,
  file.path(pub_dir, "Figure_2A_Correlation_Matrix.png"), 5.5, 4.5)
doc <- body_add_par(doc,
  "(A) Pearson correlation matrix of the initial 27 environmental covariates. Variables with |r| > 0.7 were removed to reduce multicollinearity prior to modeling.",
  style = "Normal")
doc <- add_fig(doc,
  file.path(pub_dir, "Figure_2B_Boruta_Selection.png"), 6, 3.5)
doc <- body_add_par(doc,
  "(B) Boruta feature importance analysis (maxRuns = 500). Green = confirmed important; red = rejected; yellow = tentative. Twelve variables were confirmed as statistically significant predictors and retained for modeling. Notably, Wild Boar Density was confirmed, while Soft Tick Density was rejected due to extreme geographic reporting bias (see Figure S1).",
  style = "Normal")
doc <- body_add_par(doc, "", style = "Normal")

# Figure 3
doc <- body_add_par(doc,
  "Figure 3. Wild Boar (Sus scrofa) Density Surface", style = "heading 2")
doc <- add_fig(doc,
  file.path(pub_dir, "Figure_3_Wild_Boar_Density.png"), 6, 4)
doc <- body_add_par(doc,
  "Figure 3. Continuous density surface of wild boar (Sus scrofa) across the study region, derived from 5,000 GBIF occurrence records using 2D Gaussian Kernel Density Estimation (bandwidth = 2.0 degrees). Wild boar density emerged as the second most important predictor of ASF suitability in the Random Forest model (Table 3), supporting the hypothesis that wild suids act as environmental reservoirs maintaining ASFV at the domestic-wildlife interface.",
  style = "Normal")
doc <- body_add_par(doc, "", style = "Normal")

# Figure 4
doc <- body_add_par(doc,
  "Figure 4. Model Validation: Standard CV vs. Spatial Block CV",
  style = "heading 2")
doc <- add_fig(doc,
  file.path(pub_dir, "Figure_4_Spatial_CV_Comparison.png"), 6, 4.2)
doc <- body_add_par(doc, paste0(
  "Figure 4. Comparison of model discrimination under standard 10-fold ",
  "repeated cross-validation (blue) and 5-fold spatial block cross-validation ",
  "using 500 km hexagonal blocks (red). Error bars represent +/- 1 standard ",
  "deviation across folds (XGBoost standard CV error bar omitted as it uses ",
  "the native API rather than caret resampling). Random Forest and XGBoost ",
  "maintained AUC values > 0.93 under spatial CV, with a modest decline of ",
  "~0.04 from standard CV, confirming that these models capture genuine ",
  "environmental suitability signals rather than artifacts of spatial ",
  "autocorrelation. SVM and GLM showed larger declines (~0.11), indicating ",
  "lower spatial transferability. Spatial block CV follows Valavi et al. (2019)."),
  style = "Normal")
doc <- body_add_par(doc, "", style = "Normal")

# Figure 5
doc <- body_add_par(doc,
  "Figure 5. Individual Model Risk Suitability Maps (4-Panel)",
  style = "heading 2")
doc <- add_fig(doc,
  file.path(pub_dir, "Figure_5_Individual_Model_Maps.png"), 6.5, 4.3)
doc <- body_add_par(doc,
  "Figure 5. Predicted ASF ecological niche suitability (0 = low, 1 = high) generated by each classifier: (A) Random Forest, (B) XGBoost, (C) Support Vector Machine, and (D) Logistic Regression. Tree ensemble models (A, B) produce sharper risk boundaries consistent with non-linear ecological thresholds, while the GLM baseline (D) produces a diffuse, oversmoothed surface, confirming that linear models are insufficient to capture the complex environmental niche of ASF.",
  style = "Normal")
doc <- body_add_par(doc, "", style = "Normal")

# Figure 6
doc <- body_add_par(doc,
  "Figure 6. Ensemble ASF Suitability Risk Map", style = "heading 2")
doc <- add_fig(doc,
  file.path(pub_dir, "Figure_6_Ensemble_Risk_Map.png"), 6.5, 4.6)
doc <- body_add_par(doc, paste0(
  "Figure 6. Consensus ensemble suitability map for ASF in East and Southeast ",
  "Asia, computed as the unweighted mean of predictions from all four ",
  "classifiers. White dots represent ", nrow(asf_clean), " confirmed outbreak ",
  "locations (WOAH EMPRES-i, 2015-2025). High-suitability zones (red) are ",
  "concentrated in eastern China, the Mekong Delta (Vietnam), the Philippine ",
  "archipelago, and Java (Indonesia), corresponding to areas of high domestic ",
  "pig density, elevated precipitation, and significant wild boar habitat. ",
  "The ensemble approach reduces individual model bias and provides a robust ",
  "spatial risk surface for targeted surveillance planning."),
  style = "Normal")
doc <- body_add_par(doc, "", style = "Normal")

# Figure 7
doc <- body_add_par(doc,
  "Figure 7. Partial Dependence Plots \u2014 Top 6 Predictors",
  style = "heading 2")
doc <- add_fig(doc,
  file.path(pub_dir, "Figure_7_PDP_Top6_Variables.png"), 6.5, 3.8)
doc <- body_add_par(doc,
  "Figure 7. Partial Dependence Plots (PDP) showing the marginal effect of each of the six most important predictors on the probability of ASF presence, derived from the Random Forest model. Rug plots along the x-axis show the training data distribution. (A) BIO13 shows a steep suitability increase above ~200 mm, identifying a moisture threshold. (B) Wild Boar Density shows rapid increase above moderate densities, confirming the wildlife reservoir role. (C) Pig Density demonstrates a clear step-threshold effect. (D-F) BIO14, Wind Speed, and Population Density show more gradual, non-linear responses.",
  style = "Normal")
doc <- body_add_par(doc, "", style = "Normal")

# Figure 8
doc <- body_add_par(doc,
  "Figure 8. Individual Conditional Expectation (ICE) \u2014 Pig Density",
  style = "heading 2")
doc <- add_fig(doc,
  file.path(pub_dir, "Figure_8_ICE_Pig_Density.png"), 6, 3.6)
doc <- body_add_par(doc,
  "Figure 8. ICE plot for domestic pig density from the Random Forest model. Thin blue lines represent individual response curves for 100 randomly sampled locations; the thick red line is the overall average (PDP). High variability at intermediate pig densities indicates that the pig density effect is modulated by local environmental context (precipitation, wild boar presence), supporting the existence of significant variable interactions.",
  style = "Normal")
doc <- body_add_par(doc, "", style = "Normal")


# ------ SUPPLEMENTARY FIGURES ------
doc <- body_add_par(doc, "Supplementary Figures", style = "heading 1")

doc <- body_add_par(doc,
  "Figure S1. Soft Tick (Argasidae/Ornithodoros) GBIF Occurrences",
  style = "heading 2")
doc <- add_fig(doc,
  file.path(pub_dir, "Figure_S1_Soft_Tick_Occurrences.png"), 6, 4)
doc <- body_add_par(doc,
  "Figure S1. Geographic distribution of soft tick occurrences from GBIF. Records are overwhelmingly concentrated in South Korea, reflecting a severe geographic reporting bias rather than true distribution. This extreme sampling bias renders the KDE density surface uninformative for continental-scale modeling, and Boruta correctly rejected Soft Tick Density as a non-significant predictor.",
  style = "Normal")
doc <- body_add_par(doc, "", style = "Normal")

doc <- body_add_par(doc,
  "Figure S2. Random Forest Variable Importance", style = "heading 2")
doc <- add_fig(doc,
  file.path(pub_dir, "Figure_S2_RF_Variable_Importance.png"), 5, 3.5)
doc <- body_add_par(doc,
  "Figure S2. Scaled variable importance from the Random Forest classifier, ranked by contribution to classification accuracy. Lollipop chart generated with ggplot2 for publication consistency.",
  style = "Normal")
doc <- body_add_par(doc, "", style = "Normal")

doc <- body_add_par(doc,
  "Figure S3. XGBoost Variable Importance (Gain)", style = "heading 2")
doc <- add_fig(doc,
  file.path(pub_dir, "Figure_S3_XGBoost_Variable_Importance.png"), 5, 3.5)
doc <- body_add_par(doc,
  "Figure S3. Feature importance from the native XGBoost model measured by total gain. The ranking is broadly consistent with Random Forest (Figure S2), with BIO13 and Wild Boar Density dominating.",
  style = "Normal")
doc <- body_add_par(doc, "", style = "Normal")

doc <- body_add_par(doc,
  "Figure S4. Cross-Validation AUC Distribution", style = "heading 2")
doc <- add_fig(doc,
  file.path(pub_dir, "Figure_S4_CV_AUC_Boxplot.png"), 5, 3.5)
doc <- body_add_par(doc,
  "Figure S4. Distribution of AUC values across 30 cross-validation resamples (10-fold x 3 repeats) for Random Forest, SVM, and GLM. Individual resample values are overlaid as jittered points. XGBoost is omitted as it uses the native API with a separate 10-fold CV procedure.",
  style = "Normal")
doc <- body_add_par(doc, "", style = "Normal")


# ------ TABLES ------
doc <- body_add_par(doc, "Tables", style = "heading 1")

# Table 1
doc <- body_add_par(doc, "Table 1. Environmental Covariates", style = "heading 2")
doc <- body_add_par(doc,
  "Table 1. Summary of the 27 environmental and anthropogenic covariates assembled for the ASF ecological niche model. All variables were resampled to a common 2.5 arc-minute grid (~4.5 km). After correlation filtering (|r| > 0.7) and Boruta feature selection (maxRuns = 500), 12 variables were retained as statistically significant predictors.",
  style = "Normal")
doc <- body_add_table(doc, value = table1, style = "table_template")
doc <- body_add_par(doc, "", style = "Normal")

# Table 2
doc <- body_add_par(doc,
  "Table 2. Model Performance Summary", style = "heading 2")
doc <- body_add_par(doc, paste0(
  "Table 2. Classification metrics for the four ML models on the held-out ",
  "test set (20%). Spatial block CV AUC ",
  "(5-fold, 500 km hexagonal blocks; Valavi et al., 2019) is included as a ",
  "robustness check for spatial transferability. RF and XGBoost maintained ",
  "spatial CV AUC > 0.93, confirming strong generalization."),
  style = "Normal")
doc <- body_add_table(doc, value = perf_table, style = "table_template")
doc <- body_add_par(doc, "", style = "Normal")

# Table 3
doc <- body_add_par(doc,
  "Table 3. Variable Importance Rankings (Random Forest)", style = "heading 2")
doc <- body_add_par(doc,
  "Table 3. Scaled variable importance from the Random Forest model. BIO13 (Precipitation of Wettest Month) is the dominant predictor, followed by Wild Boar Density, confirming the importance of the wildlife reservoir in shaping the ASF ecological niche.",
  style = "Normal")
doc <- body_add_table(doc, value = imp_df, style = "table_template")


# ------ SAVE ------
out_docx <- file.path(pub_dir, "Figure_Table_Legends.docx")
print(doc, target = out_docx)

cat("\n====================================================\n")
cat("  Publication Figures Generated Successfully!\n")
cat("====================================================\n")
cat("  Directory:", file.path(getwd(), pub_dir), "\n")
cat("  Main Figures:  Figure_1 through Figure_8\n")
cat("  Supplementary: Figure_S1 through Figure_S4\n")
cat("  Tables:        Table_1, Table_2, Table_3 (.csv)\n")
cat("  Legend Doc:    Figure_Table_Legends.docx\n")
cat("====================================================\n")
cat("\n  FIXES APPLIED in this version:\n")
cat("  [+] Figure 1:  Scale bar, north arrow, country labels added\n")
cat("  [+] Figure 2B: Boruta plot regenerated in ggplot2\n")
cat("  [+] Figure 4:  Error bars added for Standard CV (RF, SVM, GLM)\n")
cat("  [+] Figure 4:  Caption updated to note XGBoost std CV omission\n")
cat("  [+] Figure S2: RF importance regenerated as ggplot2 lollipop chart\n")
cat("  [+] Figure S3: XGBoost importance regenerated as ggplot2 bar chart\n")
cat("  [+] Figure S4: CV AUC boxplot regenerated in ggplot2 with jitter\n")
cat("  [+] Table 3:   Safe column access (handles Positive/Overall)\n")
cat("  [+] Removed unused yor_colors palette\n")
cat("  [+] Added clean_name() helper for display names\n")
cat("====================================================\n")
