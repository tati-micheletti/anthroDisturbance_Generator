saveDisturbances <- function(disturbanceList,
                             currentTime,
                             overwrite,
                             runName) {
  # Helper to write one spatial object
  write_one <- function(obj, sector, layer = NULL) {
    # decide extension and write func
    if (is.null(obj)) {
      warning(sprintf("The layer for %s%s is NULL. Not saving.", sector,
                      if (!is.null(layer)) paste0(" -- ", layer) else ""),
              immediate. = TRUE)
      return()
    }
    # Coerce and write vector
    if (any(class(obj) %in% c("SpatVector", "sf")) || is(obj, "Spatial")) {
      if (!inherits(obj, "SpatVector")) obj <- terra::vect(obj)
      fname <- paste0("disturbances_", sector,
                      if (!is.null(layer)) paste0("_", layer) else "",
                      "_", currentTime, "_", runName, ".shp")
      out_path <- file.path(Paths$outputPath, fname)
      message(sprintf("Writing vector disturbance to %s", out_path))
      tryCatch(
        terra::writeVector(obj, out_path, filetype = "ESRI Shapefile", overwrite = overwrite),
        error = function(e) stop(sprintf("writeVector failed for %s (%s): %s", out_path, class(obj)[1], conditionMessage(e)))
      )
    }
    # Raster branch
    else if (any(class(obj) %in% c("RasterLayer", "SpatRaster"))) {
      if (!inherits(obj, "SpatRaster")) obj <- terra::rast(obj)
      lyr_count <- tryCatch(terra::nlyr(obj), error = function(...) NA_integer_)
      fname <- paste0("disturbances_", sector,
                      if (!is.null(layer)) paste0("_", layer) else "",
                      "_", currentTime, "_", runName, ".tif")
      out_path <- file.path(Paths$outputPath, fname)
      if (length(out_path) != 1L) {
        lay_lbl <- if (is.null(layer)) "NA" else paste(layer, collapse = ";")
        stop(sprintf("saveDisturbances: expected single filename but got %d for sector=%s layer=%s (run=%s, time=%s)",
                     length(out_path), sector, lay_lbl, runName, currentTime))
      }
      message(sprintf("Writing raster disturbance to %s (layers: %s)",
                      out_path, lyr_count))
      tryCatch(
        terra::writeRaster(obj, out_path, filetype = "GTiff", overwrite = overwrite),
        error = function(e) stop(sprintf("writeRaster failed for %s (%s, %d layers): %s", out_path, class(obj)[1], lyr_count, conditionMessage(e)))
      )
    } else {
      stop(sprintf("Objects of class %s can't be used. Please use raster, sp, sf, or terra formats.",
                   paste(class(obj), collapse=",")))
    }
  }
  
  for (sector in names(disturbanceList)) {
    obj <- disturbanceList[[sector]]
    # flat list: single object
    if (any(class(obj) %in% c("SpatVector", "sf", "RasterLayer", "SpatRaster")) ||
        is(obj, "Spatial")) {
      message(sprintf("Layer: %s was likely not generated. Saving current disturbance.", sector))
      write_one(obj, sector)
    } else if (is.list(obj)) {
      # multi-layer
      layers <- setdiff(names(obj), grep("potential", names(obj), value=TRUE))
      for (lay in layers) {
        message(sprintf("Saving layer: %s -- %s", sector, lay))
        layer_obj <- obj[[lay]]
        # preserve or assign Class
        if (any(class(layer_obj) %in% c("SpatVector", "sf")) || is(layer_obj, "Spatial")) {
          if (!"Class" %in% names(layer_obj) || all(is.na(layer_obj$Class)))
            layer_obj$Class <- lay
        }
        write_one(layer_obj, sector, lay)
      }
    }
  }
  message(crayon::green(sprintf("All disturbances saved for %s", currentTime)))
}
