#===============================================================================
# TRAINING: Geospatial Data Analysis for Decision-Making
# PRACTICAL SCRIPT -- RWANDA LAND USE / LAND COVER (LULC) MAPPING
#-------------------------------------------------------------------------------
# Purpose : Starting from a satellite-derived land cover raster that has
#           ALREADY been downloaded and shared with you, this script:
#             1) sets up the working folders automatically
#             2) loads the raster
#             3) calculates how many km2 (and %) of Rwanda each land cover
#                class covers
#             4) produces a clean, labelled map with a colour legend
#             5) saves the map as a PNG and the statistics as a CSV
#
# Dataset : ESA WorldCover 10 m, already clipped to Rwanda
#           (file: landcover_rwanda_2021.tif -- shared by the trainer)
#
# No internet connection or Google Earth Engine account is required to run
# this script -- everything below runs on the file already on your computer.
#===============================================================================

## STEP 1 -- Install / load the R packages we need and set working directory ---
# pacman checks whether a package is already installed before installing it,
# so this line is safe to run every time.

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(here, terra, ggplot2, dplyr, scales)

setwd(here())

## STEP 2 -- Set up the working folders (runs automatically) ------------------
# Keeping raw data, and outputs in separate folders is good GIS practice:
# it makes the project easy to navigate and safe to re-run.
#
#   data/raw       -> place the shared .tif file here (see STEP 3)
#   outputs/maps   -> the final PNG map will be saved here
#   outputs/tables -> the statistics CSV will be saved here
#   script

project_folders <- c( "scripts","data/raw", "outputs/maps", "outputs/tables")
for (f in project_folders) {
  if (!dir.exists(f)) dir.create(f, recursive = TRUE)
}
cat("Project folders ready:", paste(project_folders, collapse = ", "), "\n\n")


## STEP 3 -- Load the pre-downloaded land cover raster -------------------------
# Copy the file the trainer shared with you ('landcover_rwanda_2021.tif')
# into the data/raw folder created above, then run this section.

raster_path <- "data/raw/landcover_rwanda_2025.tif"

if (!file.exists(raster_path)) {
  stop("File not found: ", raster_path,
       "\nPlease copy 'landcover_rwanda_2025.tif' (shared by the trainer) ",
       "into the data/raw folder, then re-run this script.")
}

land_cover_raster <- rast(raster_path)
cat("Raster loaded successfully.\n")
print(land_cover_raster)

## Data Exploration and Understanding

#1. Plot the data
plot(land_cover_raster , main = "Land Cover - Rwanda 2025")
# Get unique values and their frequencies
freq_table <- freq(land_cover_raster)
print(freq_table)

# Check the color table
coltab <- coltab(land_cover_raster)
print(coltab)
# Extract data for a specific location
point_data <- extract(land_cover_raster, 
                     cbind(30.0, -2.0))  # longitude, latitude

plot(point_data)

## STEP 4 -- Reproject to an accurate metric coordinate system -----------------
# The raster arrives in geographic coordinates (degrees), where a "pixel" is
# not a fixed size in metres. To calculate real areas (km2) we reproject to
# UTM Zone 35S -- the correct metric coordinate system for Rwanda -- using
# nearest-neighbour resampling ("near"), which is mandatory for categorical
# data such as land cover classes (it never invents new class values).

land_cover_proj <- project(land_cover_raster, "EPSG:32735", method = "near")

res_proj <- res(land_cover_proj)                 # pixel width/height, in metres
pixel_area_km2 <- (res_proj[1] * res_proj[2]) / 1e6   # m2 -> km2

cat("\nRaster resolution (projected):", round(res_proj[1], 1), "x",
    round(res_proj[2], 1), "metres\n")
cat("Pixel area:", round(pixel_area_km2, 6), "km2\n\n")


## STEP 5 -- Land cover class legend (ESA WorldCover official codes/colours) ---
# Every value in the raster (10, 20, 30 ...) corresponds to one of these
# classes. Official ESA colours are used so the map matches published
# WorldCover maps trainees may have seen before.

esa_classes <- data.frame(
  value = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 100),
  class = c("Tree cover", "Shrubland", "Grassland", "Cropland", "Built-up",
            "Bare/sparse vegetation", "Snow and ice", "Permanent water bodies",
            "Herbaceous wetland", "Mangroves", "Moss and lichen"),
  color = c("#006400", "#FFBB22", "#FFFF4C", "#F096FF", "#FA0000", "#B4B4B4",
            "#F0F0F0", "#0064C8", "#0096A0", "#00CF75", "#FAE6A0"),
  stringsAsFactors = FALSE
)


## STEP 6 -- Calculate area and percentage per class ---------------------------
# freq() counts how many pixels belong to each class -- this is the raster
# equivalent of a "zonal statistic" at national scale. Multiplying pixel
# count by pixel_area_km2 (from STEP 4) converts that into real km2.

freq_df <- freq(land_cover_proj)
freq_df <- merge(freq_df, esa_classes, by = "value", all.x = TRUE)
freq_df <- freq_df[!is.na(freq_df$class), ]              # drop no-data values

freq_df$area_km2 <- freq_df$count * pixel_area_km2
total_area_km2 <- sum(freq_df$area_km2)
freq_df$percentage <- (freq_df$area_km2 / total_area_km2) * 100
freq_df <- freq_df[order(-freq_df$area_km2), ]           # largest class first

# Print a clean summary table to the console
cat(strrep("=", 60), "\n")
cat("RWANDA LAND COVER STATISTICS (ESA WorldCover, 10 m)\n")
cat(strrep("=", 60), "\n")
cat(sprintf("Total mapped area: %s km2\n", format(round(total_area_km2), big.mark = ",")))
cat(strrep("-", 60), "\n")
for (i in seq_len(nrow(freq_df))) {
  cat(sprintf("%-25s %10s km2  (%5.1f%%)\n",
              freq_df$class[i],
              format(round(freq_df$area_km2[i]), big.mark = ","),
              freq_df$percentage[i]))
}
cat(strrep("=", 60), "\n\n")

# Save the statistics table
write.csv(freq_df[, c("value", "class", "area_km2", "percentage")],
          "outputs/tables/rwanda_landcover_statistics.csv", row.names = FALSE)
cat("Statistics saved to outputs/tables/rwanda_landcover_statistics.csv\n\n")


## STEP 7 -- Prepare the raster for mapping -------------------------------------
# Large rasters can be slow to plot. If needed, we downsample (aggregate)
# purely for the VISUAL -- the statistics above were already calculated on
# the full-resolution raster, so accuracy is not affected.

max_pixels <- 1e6
plot_raster <- land_cover_proj

if (ncell(plot_raster) > max_pixels) {
  agg_factor <- ceiling(sqrt(ncell(plot_raster) / max_pixels))
  plot_raster <- aggregate(plot_raster, fact = agg_factor, fun = "modal", na.rm = TRUE)
  cat("Raster downsampled for display only (factor", agg_factor, ")\n\n")
}

plot_df <- as.data.frame(plot_raster, xy = TRUE, na.rm = TRUE)
colnames(plot_df)[3] <- "value"

# Give the legend labels the area/percentage so the map is self-explanatory
freq_df$label <- paste0(freq_df$class, " (", sprintf("%.1f", freq_df$percentage), "%)")
plot_df <- merge(plot_df, freq_df[, c("value", "label")], by = "value")

class_colors <- setNames(freq_df$color, freq_df$label)


## STEP 8 -- Build the final map -------------------------------------------------

lulc_map <- ggplot(plot_df, aes(x = x, y = y, fill = label)) +
  geom_raster() +
  scale_fill_manual(values = class_colors, na.value = "transparent") +
  coord_equal() +
  theme_void(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10),
    legend.key = element_rect(color = "black", linewidth = 0.2),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, margin = margin(b = 10)),
    plot.caption = element_text(size = 8, hjust = 0.5, margin = margin(t = 10)),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(
    title = "Rwanda Land Use / Land Cover",
    subtitle = paste0("Total area: ", format(round(total_area_km2), big.mark = ","),
                       " km2  |  ESA WorldCover, 10 m resolution"),
    fill = "Land Cover Class",
    caption = "Source: ESA WorldCover  |  Training exercise: Geospatial Data Analysis for Decision-Making | Prepared by Murera Gisa"
  )

print(lulc_map)


## STEP 9 -- Save the map as PNG --------------------------------------------------

ggsave("outputs/maps/rwanda_landcover_map.png", lulc_map,
       width = 12, height = 8, dpi = 300, bg = "white")

cat("Map saved to outputs/maps/rwanda_landcover_map.png\n")
cat("\nPractical complete -- statistics table and map are in the outputs/ folder.\n")
