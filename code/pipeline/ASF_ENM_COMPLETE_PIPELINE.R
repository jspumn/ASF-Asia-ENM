# ==============================================================================
# PROJECT : Ecological Niche Modeling (ENM) for African Swine Fever (ASF)
# REGION  : East & Southeast Asia (90-150 E, -10-50 N)
# DATA    : EMPRES-i WOAH confirmed outbreaks + WorldClim + GLW4
# MODELS  : Random Forest | XGBoost (native) | SVM | Logistic Regression
# OUTPUT  : Risk maps (GeoTIFF + PNG), performance metrics, PDP/ICE plots
# DATE    : June 2026
# ==============================================================================
# HOW TO USE:
#   1. Open ASF_ENM_Asia.Rproj in RStudio (sets working directory automatically)
#   2. Run STEP 1 ONLY ONCE to download climate data
#   3. From second run onwards, start from STEP 2
# ==============================================================================


# ==============================================================================
# STEP 0: ENVIRONMENT SETUP
# ==============================================================================
rm(list = ls())

if (!require("pacman")) install.packages("pacman")
library(pacman)

p_load(
  terra, geodata, sf,
  caret, randomForest, xgboost, kernlab,
  Boruta,
  iml, pdp,
  ggplot2, tidyterra, scales, gridExtra, viridis, ggcorrplot, cowplot,
  dplyr, tidyverse,
  pROC, mltools,
  doParallel,
  rgbif, MASS,
  blockCV
)

setwd("C:/Users/s9240/OneDrive/Desktop/ENM project/ASF_ENM_Asia")

raw_path  <- "data_raw"
proc_path <- "data_processed"
out_path  <- "outputs"
walk(c(raw_path, proc_path, out_path), ~dir.create(.x, showWarnings = FALSE))

registerDoParallel(cores = max(1, parallel::detectCores() - 2))
cat("Cores available:", parallel::detectCores(), "\n")

e_asia <- ext(90, 150, -10, 50)

cat("=== STEP 0 COMPLETE ===\n")


# ==============================================================================
# STEP 1: DOWNLOAD & PREPARE ENVIRONMENTAL COVARIATES
# ==============================================================================
# RUN ONLY ONCE. Skip to STEP 2 on subsequent runs.
# WHY each variable: see COVARIATE_JUSTIFICATION.md

if (!file.exists(file.path(proc_path, "env_stack_asia_FINAL.tif"))) {
  cat("\n--- STEP 1: Downloading & preparing all covariates (including advisor requests) ---\n")

  # 1. WorldClim Climatic Variables
  cat("Downloading WorldClim data...\n")
  clim_raw  <- worldclim_global(var = "bio", res = 2.5, path = raw_path)
  names(clim_raw) <- paste0("BIO", 1:19)
  clim_asia <- crop(clim_raw[[c(1:7, 10:17)]], e_asia)  # exclude BIO8,9,18,19

  # Helper function to align and resample any raster to the template (clim_asia[[1]])
  template <- clim_asia[[1]]
  align_raster <- function(x, name_val, method_val = "bilinear") {
    cat(paste("Aligning & resampling:", name_val, "...\n"))
    x_cropped <- crop(x, e_asia)
    x_resampled <- resample(x_cropped, template, method = method_val)
    names(x_resampled) <- name_val
    return(x_resampled)
  }

  # 2. Elevation, Wind, and Population
  elev_raw  <- elevation_global(res = 2.5, path = raw_path)
  elev_asia <- align_raster(elev_raw, "Elevation")

  wind_raw  <- worldclim_global(var = "wind", res = 2.5, path = raw_path)
  wind_asia <- align_raster(mean(wind_raw), "Wind_Speed")

  pop_raw  <- population(year = 2020, res = 2.5, path = raw_path)
  pop_asia <- align_raster(pop_raw, "Population_Density")

  # 3. Livestock Densities (GLW4)
  pig_raw <- align_raster(rast("pig density/GLW4-2020.D-DA.PGS.tif"), "Pig_Density")
  chk_raw <- align_raster(rast("chicken density/GLW4-2020.D-DA.CHK.tif"), "Chicken_Density")

  # 4. Soil Grids (Soil Type proxies: pH and Clay fraction)
  cat("Downloading SoilGrids pH and Clay content...\n")
  soil_ph_raw   <- soil_world(var = "phh2o", depth = 5, path = raw_path)
  soil_clay_raw <- soil_world(var = "clay", depth = 5, path = raw_path)
  
  soil_ph   <- align_raster(soil_ph_raw, "Soil_pH")
  soil_clay <- align_raster(soil_clay_raw, "Soil_Clay_Content")

  # 5. Land Cover Fractions (Vegetation Density / NDVI proxies)
  cat("Downloading ESA WorldCover land use layers...\n")
  lc_trees_raw    <- landcover(var = "trees", path = raw_path)
  lc_crop_raw     <- landcover(var = "cropland", path = raw_path)
  lc_grass_raw    <- landcover(var = "grassland", path = raw_path)
  
  lc_trees    <- align_raster(lc_trees_raw, "Tree_Cover_Fraction")
  lc_cropland <- align_raster(lc_crop_raw, "Cropland_Fraction")
  lc_grass    <- align_raster(lc_grass_raw, "Grassland_Fraction")

  # 6. Wild Boar Abundance (GBIF Point Occurrences -> Gaussian KDE)
  cat("Querying GBIF for Wild Boar (Sus scrofa) occurrences in Asia...\n")
  boar_key <- name_backbone(name = "Sus scrofa")$usageKey
  res_boar <- occ_search(
    taxonKey = boar_key,
    decimalLongitude = "90,150",
    decimalLatitude = "-10,50",
    limit = 5000,
    hasCoordinate = TRUE
  )
  
  boar_coords <- na.omit(res_boar$data[, c("decimalLongitude", "decimalLatitude")])
  colnames(boar_coords) <- c("lon", "lat")
  boar_vec <- vect(boar_coords, geom = c("lon", "lat"), crs = "+proj=longlat +datum=WGS84")
  
  # Rasterize counts & apply Gaussian KDE (2.0 degree bandwidth)
  boar_counts <- rasterize(boar_vec, template, fun = "length", background = 0)
  w_mat <- focalMat(template, d = 2.0, type = "Gauss")
  wild_boar_kde <- focal(boar_counts, w = w_mat, fun = sum, na.rm = TRUE)
  names(wild_boar_kde) <- "Wild_Boar_Density"

  # 7. Soft Tick Abundance (GBIF Point Occurrences -> Gaussian KDE)
  cat("Querying GBIF for Soft Tick (Ornithodoros) occurrences in Asia...\n")
  # Use taxonKeys for both Ornithodoros genus (2184972) and Argasidae family (9168)
  res_tick_genus <- occ_search(
    taxonKey = 2184972,
    decimalLongitude = "90,150",
    decimalLatitude = "-10,50",
    limit = 2000,
    hasCoordinate = TRUE
  )
  res_tick_fam <- occ_search(
    taxonKey = 9168,
    decimalLongitude = "90,150",
    decimalLatitude = "-10,50",
    limit = 2000,
    hasCoordinate = TRUE
  )
  
  tick_data <- rbind(
    if (!is.null(res_tick_genus$data) && nrow(res_tick_genus$data) > 0) res_tick_genus$data[, c("decimalLongitude", "decimalLatitude")] else NULL,
    if (!is.null(res_tick_fam$data) && nrow(res_tick_fam$data) > 0) res_tick_fam$data[, c("decimalLongitude", "decimalLatitude")] else NULL
  )
  
  if (!is.null(tick_data) && nrow(tick_data) > 0) {
    tick_coords <- na.omit(tick_data)
    colnames(tick_coords) <- c("lon", "lat")
    tick_vec <- vect(tick_coords, geom = c("lon", "lat"), crs = "+proj=longlat +datum=WGS84")
    tick_counts <- rasterize(tick_vec, template, fun = "length", background = 0)
    tick_kde <- focal(tick_counts, w = w_mat, fun = sum, na.rm = TRUE)
  } else {
    # Fallback to zero grid if no tick records are returned
    cat("Warning: No soft tick occurrences found on GBIF. Creating blank density grid.\n")
    tick_kde <- template
    values(tick_kde) <- 0
  }
  names(tick_kde) <- "Soft_Tick_Density"

  # 8. Combine everything into one master SpatRaster stack
  env_stack <- c(
    clim_asia,
    elev_asia,
    wind_asia,
    pop_asia,
    pig_raw,
    chk_raw,
    soil_ph,
    soil_clay,
    lc_trees,
    lc_cropland,
    lc_grass,
    wild_boar_kde,
    tick_kde
  )

  writeRaster(env_stack,
              file.path(proc_path, "env_stack_asia_FINAL.tif"),
              overwrite = TRUE)
  cat("=== STEP 1 COMPLETE: Expanded stack saved ===\n")
} else {
  cat("\n--- STEP 1: Expanded covariate stack already exists. Skipping download. ---\n")
}


# ==============================================================================
# STEP 2: LOAD PRE-PROCESSED DATA
# ==============================================================================
# Start here on all runs after the first.

cat("\n--- STEP 2: Loading pre-processed data ---\n")

# Load raster stack (prefer FINAL, fall back to v3)
stack_file <- if (file.exists(file.path(proc_path, "env_stack_asia_FINAL.tif"))) {
  file.path(proc_path, "env_stack_asia_FINAL.tif")
} else {
  file.path(proc_path, "env_stack_asia_v3.tif")
}
env_stack <- rast(stack_file)

land_asia <- world(resolution = 3, path = raw_path) |> crop(e_asia)

cat("Stack loaded:", nlyr(env_stack), "layers\n")
cat("Names:", paste(names(env_stack), collapse = ", "), "\n")
cat("=== STEP 2 COMPLETE ===\n")


# ==============================================================================
# STEP 3: CLEAN ASF OCCURRENCE DATA
# ==============================================================================
# WHY: Raw EMPRES-i file has 11 metadata header lines to skip.
#      We deduplicate by spatial coordinates to avoid spatial pseudoreplication.

cat("\n--- STEP 3: Loading ASF occurrence data ---\n")

asf_raw <- read.csv(
  "EMPRESi_Raw data/overview-raw-data_202601062012.csv",
  skip = 11, stringsAsFactors = FALSE
)
cat("Raw records:", nrow(asf_raw), "\n")

asf_clean <- asf_raw |>
  filter(
    Diagnosis.status == "Confirmed",
    !is.na(Latitude), !is.na(Longitude),
    Longitude >= 90, Longitude <= 150,
    Latitude  >= -10, Latitude  <= 50
  ) |>
  dplyr::select(Latitude, Longitude, Country,
                Date = Observation.date..dd.mm.yyyy.) |>
  distinct(Latitude, Longitude, .keep_all = TRUE)

cat("After cleaning:", nrow(asf_clean), "unique outbreak locations\n")
write.csv(asf_clean, file.path(out_path, "ASF_cleaned_points.csv"),
          row.names = FALSE)
cat("=== STEP 3 COMPLETE ===\n")


# ==============================================================================
# STEP 4: BUILD THE MODELING DATASET
# ==============================================================================
# WHY: Presence-Background approach.
#   - Presence = confirmed ASF outbreak locations
#   - Background = random landscape samples (pseudo-absences)
# The model learns what environmental conditions distinguish outbreak sites
# from random background conditions.

cat("\n--- STEP 4: Building modeling dataset ---\n")

coords_presence <- asf_clean[, c("Longitude", "Latitude")]
presence_vals   <- terra::extract(env_stack, coords_presence) |>
  as.data.frame() |>
  mutate(Class = "Positive") |>
  dplyr::select(-ID)

n_presence <- nrow(presence_vals)
cat("Presence records:", n_presence, "\n")

set.seed(123)
background_sample <- spatSample(env_stack, size = n_presence,
                                na.rm = TRUE, as.df = TRUE, xy = TRUE)
bg_coords <- background_sample[, c("x", "y")]
colnames(bg_coords) <- c("Longitude", "Latitude")
background_vals <- background_sample |>
  dplyr::select(-x, -y) |>
  mutate(Class = "Negative")
cat("Background points:", nrow(background_vals), "\n")

# Combine data AND preserve coordinates for spatial block CV (blockCV)
all_coords <- rbind(
  data.frame(Longitude = coords_presence$Longitude,
             Latitude  = coords_presence$Latitude),
  bg_coords
)

full_data <- rbind(presence_vals, background_vals)
complete_rows <- complete.cases(full_data)
full_data   <- full_data[complete_rows, ] |>
  mutate(Class = factor(Class, levels = c("Negative", "Positive")))
all_coords  <- all_coords[complete_rows, ]
rownames(all_coords) <- seq_len(nrow(all_coords))

cat("Total dataset:", nrow(full_data), "rows\n")
cat("Coordinates preserved for spatial CV:", nrow(all_coords), "rows\n")
print(table(full_data$Class))
cat("=== STEP 4 COMPLETE ===\n")


# ==============================================================================
# STEP 5: FEATURE SELECTION
# ==============================================================================
# WHY:
#   5a. Correlation filter (|r| > 0.7): removes redundant variables
#   5b. Boruta: removes variables that add no predictive signal beyond chance

cat("\n--- STEP 5: Feature selection ---\n")

# 5a. Correlation filter
numeric_cols  <- full_data |> dplyr::select(-Class)
cor_matrix    <- cor(numeric_cols, use = "complete.obs")
high_cor_vars <- findCorrelation(cor_matrix, cutoff = 0.7, verbose = TRUE)

data_filtered <- if (length(high_cor_vars) > 0) {
  cat("Removing", length(high_cor_vars), "correlated variables\n")
  full_data[, -high_cor_vars]
} else {
  cat("No highly correlated variables found\n")
  full_data
}

# Save correlation plot
png(file.path(out_path, "Correlation_Matrix.png"), width = 1200, height = 1000, res = 120)
ggcorrplot(cor_matrix, hc.order = TRUE, type = "lower",
           lab = TRUE, lab_size = 2.5,
           colors = c("#6D9EC1", "white", "#E46726"),
           title = "Predictor Correlation Matrix - ASF Asia ENM") |> print()
dev.off()

cat("Variables before filter:", ncol(numeric_cols), "\n")
cat("Variables after  filter:", ncol(data_filtered) - 1, "\n")

# 5b. Boruta selection
cat("Running Boruta feature selection (3-10 min)...\n")
set.seed(124)
boruta_output <- Boruta(Class ~ ., data = data_filtered,
                        doTrace = 1, maxRuns = 500)
print(boruta_output)

png(file.path(out_path, "Boruta_Feature_Selection.png"),
    width = 1400, height = 700, res = 120)
par(mar = c(10, 4, 4, 2))
plot(boruta_output, las = 2, cex.axis = 0.8,
     main = "Boruta Feature Importance - ASF Asia ENM")
dev.off()

final_vars <- getSelectedAttributes(boruta_output, withTentative = FALSE)
data_final <- data_filtered[, c("Class", final_vars)]

cat("Boruta kept", length(final_vars), "variables:\n")
cat(paste(final_vars, collapse = ", "), "\n")

saveRDS(data_final, file.path(proc_path, "asf_final_data.rds"))
saveRDS(all_coords, file.path(proc_path, "asf_all_coords.rds"))
cat("Saved data_final and coordinate index for spatial CV.\n")
cat("=== STEP 5 COMPLETE ===\n")


# ==============================================================================
# STEP 6: TRAIN/TEST SPLIT & CROSS-VALIDATION SETUP
# ==============================================================================
# WHY: 80/20 split with 10-fold repeated CV.
#   - Training set (80%): used to build models
#   - Test set (20%):     held-out for unbiased final evaluation
#   - CV inside training: tunes hyperparameters without touching test set

cat("\n--- STEP 6: Data splitting and CV setup ---\n")

data_final <- readRDS(file.path(proc_path, "asf_final_data.rds"))
all_coords <- readRDS(file.path(proc_path, "asf_all_coords.rds"))
cc_mask    <- complete.cases(data_final)
data_final <- data_final[cc_mask, ]
all_coords <- all_coords[cc_mask, ]
data_final$Class <- factor(make.names(data_final$Class))

cat("Class levels:", paste(levels(data_final$Class), collapse = ", "), "\n")

set.seed(456)
train_idx  <- createDataPartition(data_final$Class, p = 0.8, list = FALSE)
train_data <- data_final[train_idx, ]
test_data  <- data_final[-train_idx, ]

cat("Training samples:", nrow(train_data), "\n")
cat("Testing  samples:", nrow(test_data),  "\n")

X <- data_final |> dplyr::select(-Class) |> as.data.frame()
Y <- as.character(data_final$Class)

# Cross-validation control for caret models (RF, SVM, GLM)
ctrl <- trainControl(
  method          = "repeatedcv",
  number          = 10,
  repeats         = 3,
  classProbs      = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final",
  verboseIter     = FALSE
)

cat("=== STEP 6 COMPLETE ===\n")


# ==============================================================================
# STEP 7: RANDOM FOREST
# ==============================================================================
# WHY: Ensemble of decision trees. Robust to outliers and non-linearity.
#      Tune: mtry = number of variables randomly tried at each node.

cat("\n=== STEP 7: Training Random Forest ===\n")
rf_file <- file.path(proc_path, "model_RF.rds")
if (file.exists(rf_file)) {
  cat("Loading existing Random Forest model from disk...\n")
  rf_fit <- readRDS(rf_file)
} else {
  set.seed(111)
  rf_grid <- expand.grid(mtry = 2:min(6, length(final_vars)))
  rf_fit <- train(
    Class ~ .,
    data      = train_data,
    method    = "rf",
    metric    = "ROC",
    trControl = ctrl,
    tuneGrid  = rf_grid,
    ntree     = 500,
    importance = TRUE
  )
  saveRDS(rf_fit, rf_file)
}

print(rf_fit)
cat("Best mtry:", rf_fit$bestTune$mtry, "\n")
cat("Best CV AUC:", round(max(rf_fit$results$ROC), 4), "\n")

rf_pred <- predict(rf_fit, test_data)
rf_prob <- predict(rf_fit, test_data, type = "prob")
rf_cm   <- confusionMatrix(rf_pred, test_data$Class, positive = "Positive")
rf_auc  <- roc(as.numeric(test_data$Class == "Positive"),
               rf_prob$Positive)$auc
rf_mcc  <- mcc(preds   = as.numeric(rf_pred == "Positive"),
               actuals = as.numeric(test_data$Class == "Positive"))

cat("\n--- RF Test Performance ---\n")
print(rf_cm)
cat("AUC:", round(rf_auc, 4), "| MCC:", round(rf_mcc, 4), "\n")

rf_imp   <- varImp(rf_fit, scale = TRUE)
top_vars <- rownames(rf_imp$importance)[
              order(rf_imp$importance$Positive, decreasing = TRUE)]

png(file.path(out_path, "RF_Variable_Importance.png"),
    width = 900, height = 700, res = 120)
plot(rf_imp, main = "Random Forest Variable Importance - ASF Asia")
dev.off()
# RF model already saved in Step 7
cat("=== STEP 7 COMPLETE ===\n")


# ==============================================================================
# STEP 8: XGBOOST (NATIVE API)
# ==============================================================================
# WHY native XGBoost instead of caret's xgbTree:
#   xgboost v2.0+ changed its internal C++ object type (XGBAltrepPointerClass).
#   caret's xgbTree tries to assign: modelFit$xNames <- colnames(x)
#   This assignment fails on the new object => all CV fold predictions fail
#   => all ROC metrics = NA => caret stops with "Error: Stopping".
#
# SOLUTION: Use xgboost's own native functions:
#   xgb.cv()    -> finds best nrounds via early stopping (replaces grid search)
#   xgb.train() -> trains final model on all training data
#   This works with xgboost v1.x AND v2.x.

cat("\n=== STEP 8: Training XGBoost (native API) ===\n")

# 8a. Convert data to XGBoost DMatrix format
# WHY: XGBoost requires its own data structure. Labels must be 0/1 numeric.
feature_cols <- setdiff(names(train_data), "Class")

train_mat   <- as.matrix(train_data[, feature_cols])
test_mat    <- as.matrix(test_data[, feature_cols])
train_label <- as.numeric(train_data$Class == "Positive")
test_label  <- as.numeric(test_data$Class  == "Positive")

dtrain <- xgb.DMatrix(data = train_mat, label = train_label)
dtest  <- xgb.DMatrix(data = test_mat,  label = test_label)

cat("DMatrix ready. Train:", nrow(train_mat), "rows |",
    "Test:", nrow(test_mat), "rows\n")

# 8b. XGBoost hyperparameters
# WHY each parameter:
#   objective = "binary:logistic" -> outputs probabilities for binary outcome
#   eval_metric = "auc"           -> optimize for AUC during CV
#   eta = 0.1                     -> learning rate (slower = more robust)
#   max_depth = 6                 -> tree complexity (6 is a safe default)
#   colsample_bytree = 0.8        -> 80% of features per tree (reduces overfitting)
#   subsample = 0.8               -> 80% of rows per tree (reduces overfitting)
#   nthread = 1                   -> single thread (avoids Windows OpenMP issues)
xgb_params <- list(
  objective        = "binary:logistic",
  eval_metric      = "auc",
  eta              = 0.1,
  max_depth        = 6,
  colsample_bytree = 0.8,
  subsample        = 0.8,
  min_child_weight = 1,
  gamma            = 0,
  nthread          = 1
)

# 8c. 10-fold CV with early stopping to find optimal nrounds
# WHY: Instead of a slow grid search, early stopping automatically finds the
#      number of trees where performance peaks and then starts to overfit.
cat("Running 10-fold CV with early stopping...\n")
set.seed(222)
xgb_cv <- xgb.cv(
  params                = xgb_params,
  data                  = dtrain,
  nrounds               = 500,
  nfold                 = 10,
  verbose               = 0,
  early_stopping_rounds = 20,
  prediction            = FALSE
)

# Extract best nrounds robustly across different xgboost package versions
best_nrounds <- xgb_cv$best_iteration
if (is.null(best_nrounds) || length(best_nrounds) == 0) {
  best_nrounds <- xgb_cv$early_stop$best_iteration
}
if (is.null(best_nrounds) || length(best_nrounds) == 0) {
  best_nrounds <- xgb_cv$niter
}

best_cv_auc  <- round(max(xgb_cv$evaluation_log$test_auc_mean, na.rm = TRUE), 4)
cat("Best nrounds:", best_nrounds, "| Best 10-fold CV AUC:", best_cv_auc, "\n")

# 8d. Train final model on all training data
xgb_file <- file.path(proc_path, "model_XGBoost.rds")
if (file.exists(xgb_file)) {
  cat("Loading existing XGBoost model from disk...\n")
  xgb_final <- readRDS(xgb_file)
} else {
  cat("Training final XGBoost model...\n")
  set.seed(222)
  xgb_final <- xgb.train(
    params  = xgb_params,
    data    = dtrain,
    nrounds = best_nrounds,
    verbose = 0
  )
  saveRDS(xgb_final, xgb_file)
}
cat("XGBoost loaded/trained with", best_nrounds, "trees.\n")

# 8e. Variable importance
xgb_imp_mat <- xgb.importance(model = xgb_final)
cat("\nTop XGBoost Variable Importance:\n")
print(head(xgb_imp_mat, 10))

png(file.path(out_path, "XGBoost_Variable_Importance.png"),
    width = 1000, height = 700, res = 120)
xgb.plot.importance(xgb_imp_mat,
                    top_n = min(15, nrow(xgb_imp_mat)),
                    main  = "XGBoost Variable Importance - ASF Asia ENM")
dev.off()

# 8f. Test set performance
xgb_prob_vec <- predict(xgb_final, dtest)
xgb_pred     <- factor(ifelse(xgb_prob_vec >= 0.5, "Positive", "Negative"),
                        levels = c("Negative", "Positive"))
xgb_prob     <- data.frame(Negative = 1 - xgb_prob_vec,
                            Positive = xgb_prob_vec)

xgb_cm  <- confusionMatrix(xgb_pred, test_data$Class, positive = "Positive")
xgb_auc <- roc(test_label, xgb_prob_vec)$auc
xgb_mcc <- mcc(preds   = as.numeric(xgb_pred == "Positive"),
               actuals = test_label)

cat("\n--- XGBoost Test Performance ---\n")
print(xgb_cm)
cat("Test AUC:", round(xgb_auc, 4),
    "| 10-fold CV AUC:", best_cv_auc,
    "| MCC:", round(xgb_mcc, 4), "\n")

# 8g. Prediction wrapper for spatial mapping (used in STEP 13)
pred_wrapper_xgb <- function(model, newdata) {
  newdata <- newdata[, feature_cols, drop = FALSE]
  dmat    <- xgb.DMatrix(data = as.matrix(newdata))
  predict(model, dmat)
}

# XGBoost model already saved in Step 8d
cat("=== STEP 8 COMPLETE ===\n")


# ==============================================================================
# STEP 9: SUPPORT VECTOR MACHINE (SVM)
# ==============================================================================
# WHY: Finds the optimal boundary between classes in high-dimensional space.
#      Radial Basis Function (RBF) kernel handles non-linear relationships.
#      Tune: C (penalty) and sigma (kernel bandwidth).

cat("\n=== STEP 9: Training SVM ===\n")
svm_file <- file.path(proc_path, "model_SVM.rds")
if (file.exists(svm_file)) {
  cat("Loading existing SVM model from disk...\n")
  svm_fit <- readRDS(svm_file)
} else {
  set.seed(333)
  svm_fit <- train(
    Class ~ .,
    data       = train_data,
    method     = "svmRadial",
    metric     = "ROC",
    trControl  = ctrl,
    tuneLength = 8
  )
  saveRDS(svm_fit, svm_file)
}

print(svm_fit)

svm_pred <- predict(svm_fit, test_data)
svm_prob <- predict(svm_fit, test_data, type = "prob")
svm_cm   <- confusionMatrix(svm_pred, test_data$Class, positive = "Positive")
svm_auc  <- roc(as.numeric(test_data$Class == "Positive"),
               svm_prob$Positive)$auc
svm_mcc  <- mcc(preds   = as.numeric(svm_pred == "Positive"),
               actuals = as.numeric(test_data$Class == "Positive"))

cat("\n--- SVM Test Performance ---\n")
print(svm_cm)
cat("AUC:", round(svm_auc, 4), "| MCC:", round(svm_mcc, 4), "\n")

# SVM model already saved in Step 9
cat("=== STEP 9 COMPLETE ===\n")


# ==============================================================================
# STEP 10: LOGISTIC REGRESSION (GLM)
# ==============================================================================
# WHY: Simple, interpretable baseline. Assumes linear relationship between
#      log-odds of ASF and each predictor. Useful for comparison to show
#      how much non-linear models (RF, XGBoost) improve over a linear model.

cat("\n=== STEP 10: Training Logistic Regression (GLM) ===\n")
glm_file <- file.path(proc_path, "model_GLM.rds")
if (file.exists(glm_file)) {
  cat("Loading existing GLM model from disk...\n")
  glm_fit <- readRDS(glm_file)
} else {
  set.seed(444)
  glm_fit <- train(
    Class ~ .,
    data      = train_data,
    method    = "glm",
    family    = "binomial",
    metric    = "ROC",
    trControl = ctrl
  )
  saveRDS(glm_fit, glm_file)
}

print(glm_fit)

glm_pred <- predict(glm_fit, test_data)
glm_prob <- predict(glm_fit, test_data, type = "prob")
glm_cm   <- confusionMatrix(glm_pred, test_data$Class, positive = "Positive")
glm_auc  <- roc(as.numeric(test_data$Class == "Positive"),
               glm_prob$Positive)$auc
glm_mcc  <- mcc(preds   = as.numeric(glm_pred == "Positive"),
               actuals = as.numeric(test_data$Class == "Positive"))

cat("\n--- GLM Test Performance ---\n")
print(glm_cm)
cat("AUC:", round(glm_auc, 4), "| MCC:", round(glm_mcc, 4), "\n")

# GLM model already saved in Step 10
cat("=== STEP 10 COMPLETE ===\n")


# ==============================================================================
# STEP 11: MODEL COMPARISON
# ==============================================================================
# WHY: Side-by-side performance comparison helps choose the best model
#      and understand the benefit of complexity (GLM -> RF -> XGBoost).
# NOTE: XGBoost uses native API, so it cannot go into resamples().
#       Its CV AUC comes from xgb.cv() in STEP 8 (best_cv_auc).

cat("\n--- STEP 11: Model Comparison ---\n")

# Compare the 3 caret models
caret_models <- list(RandomForest = rf_fit, SVM = svm_fit, GLM = glm_fit)
cv_results   <- resamples(caret_models)
summary(cv_results)

png(file.path(out_path, "Model_Comparison_AUC.png"),
    width = 900, height = 600, res = 120)
print(bwplot(cv_results, metric = "ROC",
             main = "CV AUC Comparison: RF vs SVM vs GLM"))
dev.off()

# Full performance table including XGBoost
perf_table <- data.frame(
  Model = c("Random Forest", "XGBoost (native)", "SVM", "GLM"),
  CV_AUC = round(c(max(rf_fit$results$ROC),
                   best_cv_auc,
                   max(svm_fit$results$ROC),
                   max(glm_fit$results$ROC)), 4),
  Test_AUC  = round(c(rf_auc,  xgb_auc, svm_auc, glm_auc), 4),
  Test_MCC  = round(c(rf_mcc,  xgb_mcc, svm_mcc, glm_mcc), 4),
  Test_Acc  = round(c(
    rf_cm$overall["Accuracy"],  xgb_cm$overall["Accuracy"],
    svm_cm$overall["Accuracy"], glm_cm$overall["Accuracy"]), 4),
  Test_Sens = round(c(
    rf_cm$byClass["Sensitivity"],  xgb_cm$byClass["Sensitivity"],
    svm_cm$byClass["Sensitivity"], glm_cm$byClass["Sensitivity"]), 4),
  Test_Spec = round(c(
    rf_cm$byClass["Specificity"],  xgb_cm$byClass["Specificity"],
    svm_cm$byClass["Specificity"], glm_cm$byClass["Specificity"]), 4)
)

print(perf_table)
write.csv(perf_table, file.path(out_path, "Model_Performance_Summary.csv"),
          row.names = FALSE)
cat("=== STEP 11 COMPLETE ===\n")


# ==============================================================================
# STEP 11b: SPATIAL BLOCK CROSS-VALIDATION (Robustness Check)
# ==============================================================================
# WHY: Random CV can overestimate AUC when outbreaks are spatially clustered,
#      because nearby training and test points share similar environments
#      purely due to geographic proximity (spatial autocorrelation).
#      Spatial block CV partitions data into geographically separated folds,
#      so the model is always tested on regions it has never seen during
#      training. This provides an honest estimate of spatial transferability.
#
# REFERENCE: Valavi R, Elith J, Lahoz-Monfort JJ, Guillera-Arroita G (2019).
#      "blockCV: An R package for generating spatially or environmentally
#       separated folds for k-fold cross-validation of species distribution
#       models." Methods in Ecology and Evolution, 10(2), 225-232.
# ==============================================================================

cat("\n--- STEP 11b: Spatial Block Cross-Validation ---\n")
cat("This provides an honest estimate of model generalization to new regions.\n\n")

# Convert data to sf object for blockCV
data_sf <- sf::st_as_sf(
  cbind(all_coords, data_final),
  coords = c("Longitude", "Latitude"),
  crs = 4326
)

# Create spatial blocks (5-fold)
# Block size = 500 km. At the continental scale of this study (60 x 60 degrees),
# 500 km blocks are large enough to exceed the spatial autocorrelation range
# of climatic variables (~100-300 km), ensuring genuine geographic independence.
cat("Creating spatial blocks (500 km hexagonal grid)...\n")
set.seed(42)
sb <- cv_spatial(
  x         = data_sf,
  column    = "Class",
  size      = 500000,   # 500 km blocks
  k         = 5,
  hexagon   = TRUE,
  selection = "random",
  iteration = 100,
  progress  = FALSE
)

cat("Spatial block folds created successfully.\n")
for (fi in seq_along(sb$folds_list)) {
  cat(sprintf("  Fold %d: train=%d, test=%d\n", fi,
              length(sb$folds_list[[fi]][[1]]),
              length(sb$folds_list[[fi]][[2]])))
}

# Function to run spatial CV: retrain on each fold, evaluate on held-out block
run_spatial_cv <- function(model_name, train_fn) {
  cat(sprintf("  Spatial CV: %s ...\n", model_name))
  fold_aucs <- sapply(seq_along(sb$folds_list), function(i) {
    idx_train <- sb$folds_list[[i]][[1]]
    idx_test  <- sb$folds_list[[i]][[2]]
    fold_train <- data_final[idx_train, ]
    fold_test  <- data_final[idx_test, ]
    
    tryCatch({
      probs  <- train_fn(fold_train, fold_test)
      actual <- as.numeric(fold_test$Class == "Positive")
      as.numeric(pROC::auc(pROC::roc(actual, probs, quiet = TRUE)))
    }, error = function(e) {
      cat(sprintf("    Fold %d error: %s\n", i, conditionMessage(e)))
      NA_real_
    })
  })
  
  m <- round(mean(fold_aucs, na.rm = TRUE), 4)
  s <- round(sd(fold_aucs, na.rm = TRUE), 4)
  cat(sprintf("    -> %s Spatial CV AUC: %.4f +/- %.4f\n", model_name, m, s))
  c(mean = m, sd = s)
}

# Minimal internal CV for hyperparameter tuning within each spatial fold
ctrl_inner <- trainControl(
  method = "cv", number = 5,
  classProbs = TRUE, summaryFunction = twoClassSummary, verboseIter = FALSE
)

# Random Forest
rf_sp <- run_spatial_cv("Random Forest", function(tr, te) {
  set.seed(111)
  fit <- train(Class ~ ., data = tr, method = "rf", metric = "ROC",
               trControl = ctrl_inner,
               tuneGrid = expand.grid(mtry = rf_fit$bestTune$mtry),
               ntree = 500)
  predict(fit, te, type = "prob")$Positive
})

# XGBoost (native API)
xgb_sp <- run_spatial_cv("XGBoost", function(tr, te) {
  fc <- setdiff(names(tr), "Class")
  dtr <- xgb.DMatrix(data = as.matrix(tr[, fc]),
                      label = as.numeric(tr$Class == "Positive"))
  dte <- xgb.DMatrix(data = as.matrix(te[, fc]))
  set.seed(222)
  fit <- xgb.train(params = xgb_params, data = dtr,
                    nrounds = best_nrounds, verbose = 0)
  predict(fit, dte)
})

# SVM
svm_sp <- run_spatial_cv("SVM", function(tr, te) {
  set.seed(333)
  fit <- train(Class ~ ., data = tr, method = "svmRadial", metric = "ROC",
               trControl = ctrl_inner, tuneLength = 5)
  predict(fit, te, type = "prob")$Positive
})

# GLM
glm_sp <- run_spatial_cv("GLM", function(tr, te) {
  set.seed(444)
  fit <- train(Class ~ ., data = tr, method = "glm", family = "binomial",
               metric = "ROC", trControl = ctrl_inner)
  predict(fit, te, type = "prob")$Positive
})

# Add spatial CV columns to performance table
perf_table$Spatial_CV_AUC <- c(rf_sp["mean"], xgb_sp["mean"],
                                svm_sp["mean"], glm_sp["mean"])
perf_table$Spatial_CV_SD  <- c(rf_sp["sd"], xgb_sp["sd"],
                                svm_sp["sd"], glm_sp["sd"])

cat("\n--- Updated Performance Table (with Spatial Block CV) ---\n")
print(perf_table)
write.csv(perf_table, file.path(out_path, "Model_Performance_Summary.csv"),
          row.names = FALSE)
cat("=== STEP 11b COMPLETE ===\n")


# ==============================================================================
# STEP 12: VARIABLE IMPORTANCE & PARTIAL DEPENDENCE PLOTS
# ==============================================================================
# WHY: "Black box" models need explanation.
#   - PDP (Partial Dependence Plot): average effect of one variable on risk,
#     holding all others at their mean. Shows the DIRECTION of effect.
#   - ICE (Individual Conditional Expectation): PDP per observation.
#     Shows whether the effect is consistent or varies by location.

cat("\n--- STEP 12: Interpretation (PDP + ICE) ---\n")

top6 <- top_vars[1:min(6, length(top_vars))]
cat("Top 6 variables for PDP:", paste(top6, collapse = ", "), "\n")

# Partial Dependence Plots for top 6 variables (RF model)
pdp_plots <- lapply(top6, function(v) {
  pd <- pdp::partial(rf_fit, pred.var = v,
                which.class = "Positive", prob = TRUE,
                train = train_data)
  ggplot(pd, aes_string(x = v, y = "yhat")) +
    geom_line(color = "#bd0026", linewidth = 1.2) +
    geom_rug(data = train_data, aes_string(x = v), sides = "b",
             alpha = 0.15, color = "grey40", inherit.aes = FALSE) +
    theme_minimal(base_size = 10) +
    labs(title = v, y = "P(ASF Positive)") +
    theme(plot.title = element_text(face = "bold", size = 9))
})

p_pdp <- plot_grid(plotlist = pdp_plots, ncol = 3, labels = "AUTO")
ggsave(file.path(out_path, "PDP_Top6_Variables.png"),
       p_pdp, width = 14, height = 8, dpi = 200)
cat("PDP plots saved.\n")

# ICE plot for Pig_Density (most biologically important variable)
if ("Pig_Density" %in% names(train_data)) {
  cat("Computing ICE plot for Pig_Density...\n")
  ice_pig <- pdp::partial(rf_fit, pred.var = "Pig_Density",
                     ice = TRUE, center = FALSE, prob = TRUE,
                     which.class = "Positive",
                     train = train_data, nsamples = 100)

  p_ice <- ggplot(ice_pig, aes(x = Pig_Density, y = yhat, group = yhat.id)) +
    geom_line(alpha = 0.2, color = "steelblue", linewidth = 0.4) +
    stat_summary(aes(group = 1), fun = mean, geom = "line",
                 color = "#bd0026", linewidth = 1.8) +
    theme_minimal(base_size = 12) +
    labs(
      title    = "ICE Plot: ASF Risk vs. Pig Density",
      subtitle = "Thin blue = individual location | Thick red = average trend (PDP)",
      x = "Pig Density (heads/km2)", y = "P(ASF Positive)"
    ) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(file.path(out_path, "ICE_Pig_Density.png"),
         p_ice, width = 10, height = 6, dpi = 200)
  cat("ICE plot saved.\n")
}

cat("=== STEP 12 COMPLETE ===\n")


# ==============================================================================
# STEP 13: GENERATE SPATIAL RISK MAPS (GeoTIFF)
# ==============================================================================
# WHY: Apply each trained model to every pixel in the raster stack to produce
#      a continuous risk surface (0 = low risk, 1 = high risk).
# NOTE: XGBoost uses pred_wrapper_xgb (native API). RF/SVM/GLM use pred_wrapper_caret.

cat("\n--- STEP 13: Generating spatial risk maps ---\n")

# Generic wrapper for caret models
pred_wrapper_caret <- function(model, newdata) {
  predict(model, newdata, type = "prob")[, "Positive"]
}

# XGBoost-specific wrapper (native API requires DMatrix conversion)
pred_wrapper_xgb <- function(model, newdata) {
  newdata <- newdata[, feature_cols, drop = FALSE]
  dmat    <- xgb.DMatrix(data = as.matrix(newdata))
  predict(model, dmat)
}

# Subset raster to Boruta-selected variables
final_predictors <- env_stack[[names(env_stack) %in% final_vars]]
cat("Mapping with layers:", paste(names(final_predictors), collapse = ", "), "\n")

risk_rasters <- list()

# RF
cat("Predicting: Random Forest...\n")
risk_rasters[["RF"]] <- terra::predict(final_predictors, rf_fit,
                                        fun = pred_wrapper_caret, na.rm = TRUE)
writeRaster(risk_rasters[["RF"]],
            file.path(out_path, "ASF_Risk_RF.tif"), overwrite = TRUE)
cat("  Saved ASF_Risk_RF.tif\n")

# XGBoost (native)
cat("Predicting: XGBoost...\n")
risk_rasters[["XGBoost"]] <- terra::predict(final_predictors, xgb_final,
                                             fun = pred_wrapper_xgb, na.rm = TRUE)
writeRaster(risk_rasters[["XGBoost"]],
            file.path(out_path, "ASF_Risk_XGBoost.tif"), overwrite = TRUE)
cat("  Saved ASF_Risk_XGBoost.tif\n")

# SVM
cat("Predicting: SVM...\n")
risk_rasters[["SVM"]] <- terra::predict(final_predictors, svm_fit,
                                         fun = pred_wrapper_caret, na.rm = TRUE)
writeRaster(risk_rasters[["SVM"]],
            file.path(out_path, "ASF_Risk_SVM.tif"), overwrite = TRUE)
cat("  Saved ASF_Risk_SVM.tif\n")

# GLM
cat("Predicting: GLM...\n")
risk_rasters[["GLM"]] <- terra::predict(final_predictors, glm_fit,
                                         fun = pred_wrapper_caret, na.rm = TRUE)
writeRaster(risk_rasters[["GLM"]],
            file.path(out_path, "ASF_Risk_GLM.tif"), overwrite = TRUE)
cat("  Saved ASF_Risk_GLM.tif\n")

cat("=== STEP 13 COMPLETE ===\n")


# ==============================================================================
# STEP 14: PUBLICATION-QUALITY RISK MAP FIGURES
# ==============================================================================
# WHY: GeoTIFFs are data files. For papers/reports we need annotated figures
#      with colour scale, country borders, titles, and captions.
# COLOUR SCALE: Yellow -> Orange -> Dark Red (standard epidemiological palette)

cat("\n--- STEP 14: Creating publication figures ---\n")

yor_colors <- c("#ffffcc", "#ffeda0", "#fed976", "#feb24c",
                "#fd8d3c", "#f03b20", "#bd0026")

auc_vals <- c(RF = rf_auc, XGBoost = xgb_auc, SVM = svm_auc, GLM = glm_auc)

make_risk_map <- function(risk_rast, model_name, auc_val) {
  ggplot() +
    geom_spatraster(data = risk_rast) +
    scale_fill_gradientn(
      colors   = yor_colors,
      values   = rescale(c(0, 0.15, 0.30, 0.45, 0.60, 0.80, 1.0)),
      limits   = c(0, 1),
      name     = "Suitability\nIndex",
      breaks   = c(0, 0.25, 0.50, 0.75, 1.0),
      labels   = c("0.00 (Very Low)", "0.25 (Low)", "0.50 (Moderate)",
                   "0.75 (High)", "1.00 (Very High)"),
      na.value = "transparent",
      guide    = guide_colorbar(barwidth = 1.2, barheight = 8)
    ) +
    geom_spatvector(data = land_asia, fill = NA,
                    color = "grey20", linewidth = 0.35) +
    coord_sf(xlim = c(90, 150), ylim = c(-10, 50), expand = FALSE) +
    labs(
      title    = paste("ASF Ecological Niche -", model_name),
      subtitle = paste0("Test AUC: ", round(auc_val, 3),
                        " | Confirmed outbreaks 2015-2025"),
      x = "Longitude (E)", y = "Latitude (N)",
      caption  = paste0("Resolution: 2.5 arc-min (~4.5 km) | ",
                        "Covariates: WorldClim 2.1, GLW4, WorldPop 2020")
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title       = element_text(face = "bold", size = 14),
      plot.subtitle    = element_text(size = 10, color = "grey30"),
      plot.caption     = element_text(size = 7,  color = "grey50", hjust = 0),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
      legend.position  = "right",
      legend.title     = element_text(face = "bold", size = 10)
    )
}

map_list <- list()
for (m in names(risk_rasters)) {
  p_map <- make_risk_map(risk_rasters[[m]], m, auc_vals[m])
  map_list[[m]] <- p_map
  ggsave(file.path(out_path, paste0("ASF_RiskMap_", m, ".png")),
         p_map, width = 12, height = 8, dpi = 300)
  cat("Map saved:", m, "\n")
}

# 4-panel combined figure
p_4panel <- plot_grid(
  map_list$RF      + theme(legend.position = "none"),
  map_list$XGBoost + theme(legend.position = "none"),
  map_list$SVM     + theme(legend.position = "none"),
  map_list$GLM     + theme(legend.position = "none"),
  ncol = 2,
  labels     = c("(A) Random Forest", "(B) XGBoost",
                 "(C) SVM",           "(D) Logistic Regression"),
  label_size = 10
)
ggsave(file.path(out_path, "ASF_RiskMaps_All4Models.png"),
       p_4panel, width = 18, height = 12, dpi = 300)
cat("4-panel combined map saved.\n")
cat("=== STEP 14 COMPLETE ===\n")


# ==============================================================================
# STEP 15: ENSEMBLE RISK MAP
# ==============================================================================
# WHY: Averaging predictions from all 4 models reduces individual model bias.
#      The ensemble is more stable and generally more accurate than any single model.
#      Outbreak points are overlaid to validate spatial prediction.

cat("\n--- STEP 15: Ensemble risk map ---\n")

risk_stack    <- rast(list(risk_rasters$RF, risk_rasters$XGBoost,
                           risk_rasters$SVM, risk_rasters$GLM))
ensemble_risk <- mean(risk_stack)
names(ensemble_risk) <- "Ensemble_Risk"

writeRaster(ensemble_risk,
            file.path(out_path, "ASF_Risk_ENSEMBLE.tif"), overwrite = TRUE)

p_ensemble <- ggplot() +
  geom_spatraster(data = ensemble_risk) +
  scale_fill_gradientn(
    colors   = yor_colors,
    values   = rescale(c(0, 0.15, 0.30, 0.45, 0.60, 0.80, 1.0)),
    limits   = c(0, 1),
    name     = "Ensemble\nSuitability",
    breaks   = c(0, 0.25, 0.50, 0.75, 1.0),
    labels   = c("0.00", "0.25", "0.50", "0.75", "1.00"),
    na.value = "transparent"
  ) +
  geom_point(data  = asf_clean,
             aes(x = Longitude, y = Latitude),
             color = "black", fill = "white",
             shape = 21, size = 1.2, stroke = 0.4, alpha = 0.7) +
  geom_spatvector(data = land_asia, fill = NA,
                  color = "grey20", linewidth = 0.35) +
  coord_sf(xlim = c(90, 150), ylim = c(-10, 50), expand = FALSE) +
  labs(
    title    = "Ensemble ASF Ecological Niche - East & Southeast Asia",
    subtitle = paste0("Average of RF + XGBoost + SVM + GLM | ",
                      nrow(asf_clean), " confirmed outbreaks (WOAH EMPRES-i, Jan 2026)"),
    x = "Longitude (E)", y = "Latitude (N)",
    caption = "White dots = confirmed ASF outbreaks | Color = ensemble suitability (0-1)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title        = element_text(face = "bold", size = 15),
    plot.subtitle     = element_text(size = 10, color = "grey30"),
    plot.caption      = element_text(size = 8,  color = "grey50", hjust = 0),
    panel.grid.major  = element_line(color = "grey90", linewidth = 0.2),
    legend.position   = c(0.06, 0.75),
    legend.background = element_rect(fill = "white", color = "grey50",
                                     linewidth = 0.4),
    legend.title      = element_text(face = "bold", size = 10)
  )

ggsave(file.path(out_path, "ASF_RiskMap_ENSEMBLE_FINAL.png"),
       p_ensemble, width = 14, height = 10, dpi = 300)
cat("Ensemble map saved.\n")
cat("=== STEP 15 COMPLETE ===\n")


# ==============================================================================
# FINAL SUMMARY
# ==============================================================================
cat("\n")
cat("====================================================\n")
cat("       ASF ENM PIPELINE - COMPLETE SUMMARY          \n")
cat("====================================================\n")
cat(sprintf("  Outbreak records used     : %d\n", nrow(asf_clean)))
cat(sprintf("  Predictor variables kept  : %d\n", length(final_vars)))
cat(sprintf("  Variables               : %s\n", paste(final_vars, collapse = ", ")))
cat("----------------------------------------------------\n")
cat("  MODEL PERFORMANCE (Test Set)\n")
cat("----------------------------------------------------\n")
for (i in seq_len(nrow(perf_table))) {
  cat(sprintf("  %-18s  AUC: %.4f  MCC: %.4f  Acc: %.4f\n",
              perf_table$Model[i],
              perf_table$Test_AUC[i],
              perf_table$Test_MCC[i],
              perf_table$Test_Acc[i]))
}
cat("----------------------------------------------------\n")
if ("Spatial_CV_AUC" %in% names(perf_table)) {
  cat("  SPATIAL BLOCK CV (Robustness Check)\n")
  cat("----------------------------------------------------\n")
  for (i in seq_len(nrow(perf_table))) {
    cat(sprintf("  %-18s  Spatial AUC: %.4f +/- %.4f\n",
                perf_table$Model[i],
                perf_table$Spatial_CV_AUC[i],
                perf_table$Spatial_CV_SD[i]))
  }
  cat("----------------------------------------------------\n")
}
cat("  OUTPUT FILES (outputs/)\n")
cat("    ASF_RiskMap_RF.png\n")
cat("    ASF_RiskMap_XGBoost.png\n")
cat("    ASF_RiskMap_SVM.png\n")
cat("    ASF_RiskMap_GLM.png\n")
cat("    ASF_RiskMaps_All4Models.png\n")
cat("    ASF_RiskMap_ENSEMBLE_FINAL.png\n")
cat("    Model_Performance_Summary.csv\n")
cat("    PDP_Top6_Variables.png\n")
cat("    ICE_Pig_Density.png\n")
cat("====================================================\n")
cat("  Pipeline complete.\n")
cat("====================================================\n")

# ==============================================================================
# END OF SCRIPT
# ==============================================================================
