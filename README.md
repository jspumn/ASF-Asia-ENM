# African Swine Fever (ASF) in Asia – Ecological Niche Modelling (ENM)

This repository contains the R pipeline and code infrastructure for the ASF Ecological Niche Modelling project.

The project is structured to be **completely portable**—you can clone or copy this folder to any computer, and the R scripts will automatically detect their location without requiring manual path changes.

## 📂 Project Structure

For full functionality (data processing and modelling), your local folder should follow this hierarchy (large data files are excluded from Git):

\\\	ext
ENM project/
├─ raw_data/                # Original source files (WorldClim, EMPRESi, GBIF, etc.)
├─ data/
│   ├─ raster/              # Processed raster stacks (*.tif, *.rds)
│   └─ tables/              # CSV tables (model performance, variable importance)
├─ code/                    # 👈 TRACKED IN GITHUB
│   ├─ pipeline/            # Main modelling pipeline
│   ├─ figures/             # Figure-generation script
│   ├─ utils/               # Helper R functions
│   └─ config.R             # Centralised path definitions (portable)
├─ output/                  # Outputs from the model
│   ├─ publication_figures/ # PNG/TIFF figures for the manuscript
│   ├─ tables/              # CSV tables exported for the manuscript
│   └─ logs/                # Run-time logs, Session_Summary_Log.md
├─ report/                  # Word report, legend document, manuscript drafts
└─ misc/                    # Rproj, .Rhistory, .RData, LICENSE, etc.
\\\

## 🚀 Quick-Start Guide (Running on a new computer)

1. **Install Required R Packages** (Run once in R console):
   \\\
   install.packages(c(
     "terra", "geodata", "sf", "ggplot2", "tidyterra", "scales", "cowplot",
     "dplyr", "viridis", "caret", "ggspatial", "xgboost", "tidyr", "blockCV",
     "Boruta", "pROC", "mltools", "doParallel", "rgbif", "MASS"
   ))
   \\\

2. **Open the Project:**
   Double-click the \.Rproj\ file inside the \misc/\ or root folder, or open the folder directly in RStudio.

3. **Run the Full Pipeline:**
   This script will read raw data, process rasters, train the XGBoost model, and generate prediction outputs.
   \\\
   source(file.path("code", "pipeline", "ASF_ENM_COMPLETE_PIPELINE.R"))
   \\\

4. **Generate Publication-Ready Figures:**
   Once the pipeline is complete, run the figure script to generate all manuscript visuals.
   \\\
   source(file.path("code", "figures", "generate_publication_figures.R"))
   \\\

*All figures, tables, and logs will be automatically saved into the \output/\ directory.*
