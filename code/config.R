# -------------------------------------------------------------
# ENM Project – Central configuration (portable paths)
# -------------------------------------------------------------
PROJECT_ROOT <- normalizePath(dirname(sys.frame(1)$ofile))

RAW_DIR      <- file.path(PROJECT_ROOT, "raw_data")
DATA_RASTER  <- file.path(PROJECT_ROOT, "data", "raster")
DATA_TABLES  <- file.path(PROJECT_ROOT, "data", "tables")
CODE_PIPE    <- file.path(PROJECT_ROOT, "code", "pipeline")
CODE_FIG     <- file.path(PROJECT_ROOT, "code", "figures")
UTILS_DIR    <- file.path(PROJECT_ROOT, "code", "utils")
OUTPUT_FIG   <- file.path(PROJECT_ROOT, "output", "publication_figures")
OUTPUT_TABLE <- file.path(PROJECT_ROOT, "output", "tables")
OUTPUT_LOG   <- file.path(PROJECT_ROOT, "output", "logs")
REPORT_DIR   <- file.path(PROJECT_ROOT, "report")
MISC_DIR     <- file.path(PROJECT_ROOT, "misc")
