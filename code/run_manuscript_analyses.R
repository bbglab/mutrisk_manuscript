# Run manuscript analyses:
options(warn = -1)  # suppress all warnings
library(tictoc) # for speed testing

# -- utils
run_script <- function(script, stage) {
  print(script)
  tic()
  status <- system2("Rscript", args = c("--slave", shQuote(script)))
  toc()
  gc()

  if (!identical(status, 0L)) {
    message(sprintf("[FAIL] %s (stage: %s, exit: %s)", script, stage, status))
    quit(save = "no", status = 1)
  }
}

# 0 (download data sources) automatically:
if(!all(dir.exists("raw_data/blood/mutational_signatures_analysis/burden_all.txt/"),
        dir.exists("raw_data/colon//"))) {
source("code/01_download_data/Download_data.R")
}

# additionally, download the data from BoostDM, BoostDM-CH, GENIE and COSMIC
# 1. Run Mutrisk mutation rate estimates for all tissues (expected time ~1h)
scripts = list.files("code/1_tissues/", recursive = TRUE, pattern = "pre", full.names = TRUE)
for (script in scripts) {
  run_script(script, "1_tissues")
}

# 2. Run resources scripts
resource_scripts = list.files("code/2_resources//",  full.names = TRUE)
for (script in resource_scripts) {
  run_script(script, "2_resources")
}

# 3. Run Figure scripts (including Supplementary, and Supplementary Note figures)
figure_scripts = list.files("code/3_figures/",  full.names = TRUE, recursive = TRUE)
for (script in figure_scripts) {
  run_script(script, "3_figures")
}

message("\nPipeline completed successfully.")
