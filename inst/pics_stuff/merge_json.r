# Load helper
library(jsonlite)

# Where your JPEG + JSON files live
#setwd("D:/GPhotos_Clean")       # <-- change path
dir <- here::here("inst","pics_stuff","exiftool.exe","experiment")
#setwd(dir)

# Get ExifTool path
et <- here::here("inst","pics_stuff","exiftool.exe")

# Loop through every JPEG
jpgs <- list.files(path = dir, pattern = "\\.jpg$", ignore.case = TRUE)

for (jpg in jpgs) {

  # Grab one or more matching JSONs
  json_candidates <- list.files(pattern = paste0("^", jpg, "\\.supplemental-metadata\\.json$"),
                                ignore.case = TRUE, full.names = TRUE)

  if (length(json_candidates) == 0) next  # nothing to do

  # --- Pick the JSON with most face regions (ties ⇒ newest file) ------------
  meta_list <- lapply(json_candidates, fromJSON)
  scores    <- sapply(meta_list, function(x)
    length(x$imageViews %||% list()))  # null-safe
  best_meta <- meta_list[[which.max(scores)]]

  # --- 1. Timestamps ---------------------------------------------------------
  ts <- as.numeric(best_meta$photoTakenTime$timestamp)
  dt <- format(as.POSIXct(ts, origin = "1970-01-01", tz = "UTC"),
               "%Y:%m:%d %H:%M:%S")

  # --- 2. GPS (if any) -------------------------------------------------------
  geo <- best_meta$geoData
  gps_args <- character()
  if (!is.null(geo$latitude) && !is.na(geo$latitude)) {
    latRef <- ifelse(geo$latitude  >= 0, "N", "S")
    lonRef <- ifelse(geo$longitude >= 0, "E", "W")
    gps_args <- c(
      sprintf("-GPSLatitude=%s",  geo$latitude),
      sprintf("-GPSLongitude=%s", geo$longitude),
      paste0("-GPSLatitudeRef=",  latRef),
      paste0("-GPSLongitudeRef=", lonRef)
    )
  }

  # --- 3. Face regions → MWG-RS XMP -----------------------------------------
  fr_args <- character()
  if (length(best_meta$imageViews) > 0) {
    iv <- best_meta$imageViews
    for (i in seq_along(iv)) {
      p <- iv[[i]]$person$displayName %||% "Unknown"
      r <- iv[[i]]$region
      fr_args <- c(fr_args,
                   sprintf('-RegionName[%d]=%s', i-1, shQuote(p, type="sh")),
                   sprintf('-RegionAreaX[%d]=%.6f', i-1, r$left),
                   sprintf('-RegionAreaY[%d]=%.6f', i-1, r$top),
                   sprintf('-RegionAreaW[%d]=%.6f', i-1, r$width),
                   sprintf('-RegionAreaH[%d]=%.6f', i-1, r$height)
      )
    }
  }

  # --- 4. Build & run ExifTool command --------------------------------------
  cmd <- c(
    et, "-overwrite_original",
    sprintf('-DateTimeOriginal=%s', dt),
    sprintf('-CreateDate=%s',      dt),
    sprintf('-ModifyDate=%s',      dt),
    gps_args, fr_args,
    shQuote(jpg, type = "cmd")          # the target file
  )

  system(paste(cmd, collapse = " "), show.output.on.console = FALSE)
  message("✓ patched ", jpg)
}
