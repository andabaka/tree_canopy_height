# ZAGREB TREE CANOPY HEIGHT DERIVED FROM Meta/WRI data (https://www.sciencedirect.com/science/article/pii/S003442572300439X)
# This script analyzes tree canopy data for Zagreb and creates an interactive
# map with multiple layers showing different analyses


# 1. PACKAGE INSTALLATION AND LOADING
# -------------------------------------
# Install required packages (uncomment if needed)
# install.packages("devtools")
# install.packages("pacman")
# install.packages("tidyverse")
# install.packages("osmdata")
# install.packages("viridis")
# install.packages("future")
# install.packages("future.apply")
# install.packages("exactextractr")
# install.packages("RColorBrewer")
# install.packages("leaflet.extras")
# devtools::install_github("rstudio/leaflet")
# devtools::install_github("TESS-Laboratory/chmloader")


# Load all required packages with pacman
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
    chmloader,     # For canopy height model data
    terra,         # For raster data processing
    sf,            # For spatial features
    osmdata,       # For OpenStreetMap data
    leaflet,       # For interactive maps
    leaflet.extras,# For extra leaflet functionality
    htmlwidgets,   # For saving web widgets
    tidyverse,     # For data manipulation and visualization
    viridis,       # For color palettes
    exactextractr, # For extracting raster data by polygon
    RColorBrewer   # For color palettes
)

# -----------------------------------------------------------------------------
# 2. CONFIGURATION AND OPTIMIZATION
# -----------------------------------------------------------------------------
# Enable disk compression for terra to reduce disk space requirements
options(terra.disk.compression = TRUE)

# Set up parallel processing to speed up computations (adjusted for MacBook Pro M3, 18GB)
# We leave 2 cores free for system operations to maintain responsiveness
future::plan(future::multisession, workers = future::availableCores() - 2)

# Set memory limit for parallel processing (adjust based on available RAM)
options(future.globals.maxSize = 10 * 1024^3) # 10GB limit for parallel processing


# 3. DEFINE ZAGREB BOUNDARIES
# ---------------------------
# Define a bounding box that covers the Zagreb area
zagreb_bbox <- c(15.8, 45.7, 16.15, 45.9)

# Get the official City of Zagreb boundary (admin_level=6)
zagreb_query <- opq(bbox = zagreb_bbox) %>%
    add_osm_feature(key = "admin_level", value = "6") %>%
    add_osm_feature(key = "name", value = "Grad Zagreb")

zagreb_data <- osmdata_sf(zagreb_query)

# Extract the city boundary
zagreb_city <- zagreb_data$osm_multipolygons %>%
    filter(name == "Grad Zagreb") %>%
    select(name, osm_id) %>%
    st_transform(4326)  # Ensure WGS84 projection


# Download all districts (admin_level=9) within the bounding box
districts_query <- opq(bbox = zagreb_bbox) %>%
    add_osm_feature(key = "admin_level", value = "9") %>%
    add_osm_feature(key = "name")

districts_data <- osmdata_sf(districts_query)

# Extract districts with "Gradska četvrt" name pattern
all_districts <- districts_data$osm_multipolygons %>%
    filter(grepl("Gradska četvrt", name)) %>%
    select(name, osm_id) %>%
    st_transform(4326)

# Keep only districts that intersect with the official city boundary
zagreb_districts <- st_intersection(all_districts, zagreb_city)

# Clean up district names for better display
zagreb_districts <- zagreb_districts %>%
    mutate(district_name = str_replace(name, "Gradska četvrt ", ""))

# Print the districts to verify
print(zagreb_districts$district_name)


# Create a visualization of the administrative boundaries
# This plots the city boundary and district boundaries
city_plot <- ggplot() +
    geom_sf(data = zagreb_city, fill = "lightgrey", color = "black") +
    geom_sf(data = zagreb_districts, fill = NA, color = "red", size = 0.5) +
    theme_void() +
    theme(plot.title = element_text(size = 14, hjust = 0.5)) +
    labs(title = "Official Administrative Boundaries of Zagreb")

# Save the plot as a PNG image
ggsave(
    "zagreb-boundaries.png",
    city_plot,
    width = 10,
    height = 8,
    units = "in",
    bg = "white",
    dpi = 300
)

# 4. DOWNLOAD TREE CANOPY HEIGHT MODEL
# -------------------------------
# Download canopy height model for the entire city using the Meta/WRI dataset

zagreb_chm <- chmloader::download_chm(
    zagreb_city,
    filename = "zagreb-full-chm.tif"
)

# Replace zero values with NA to exclude areas without vegetation
zagreb_chm_clean <- terra::ifel(
    zagreb_chm == 0,
    NA,
    zagreb_chm
)

# Transform city boundary to match the CRS of the raster
zagreb_city_transformed <- st_transform(zagreb_city, crs(zagreb_chm_clean))

# Mask the raster to only include areas within Zagreb
zagreb_chm_masked <- terra::mask(
    zagreb_chm_clean,
    terra::vect(zagreb_city_transformed)
)

# Create an aggregated version (10m resolution) for faster processing
zagreb_chm_agg <- terra::aggregate(zagreb_chm_masked, fact=10, fun="mean")

chm_vis <- terra::plot(
    zagreb_chm_masked,
    col = viridis::viridis(64),
    main = "Tree Canopy Height in Zagreb (1m Resolution)"
)

# 5. DISTRICT-LEVEL CANOPY METRICS
# -------------------------------
# Transform districts to match the CRS of the raster
zagreb_districts_transformed <- st_transform(zagreb_districts, crs(zagreb_chm_agg))

# Function to extract key metrics from the raster for each district
calculate_district_metrics <- function(districts, chm_raster) {

    # Make sure the CRS matches
    districts <- st_transform(districts, crs(chm_raster))

    # Create a list to store results
    district_results <- list()

    for (i in 1:nrow(districts)) {
        district <- districts[i,]
        district_name <- district$district_name

        message(paste("Processing district:", district_name))

        # Extract all values within the district boundary
        # This gives us the raw values without calculating statistics
        tryCatch({
            # Extract raw values
            extracted_values <- exactextractr::exact_extract(
                chm_raster,
                district,
                include_xy = TRUE
            )

            if (length(extracted_values) > 0 && nrow(extracted_values[[1]]) > 0) {
                values_df <- extracted_values[[1]]

                # Filter to cells with significant coverage
                values_df <- values_df[values_df$coverage_fraction > 0.5, ]

                # Calculate statistics manually
                non_na_values <- values_df$value[!is.na(values_df$value)]

                if (length(non_na_values) > 0) {
                    mean_height <- mean(non_na_values)
                    max_height <- max(non_na_values)

                    # Calculate canopy cover (proportion of non-NA cells)
                    total_cells <- nrow(values_df)
                    canopy_cells <- sum(!is.na(values_df$value))
                    canopy_cover_pct <- (canopy_cells / total_cells) * 100
                } else {
                    mean_height <- NA
                    max_height <- NA
                    canopy_cover_pct <- 0
                }

                # Calculate area of district in km²
                area_km2 <- as.numeric(st_area(district) / 1e6)

                # Store results
                district_results[[i]] <- data.frame(
                    district_name = district_name,
                    mean_height = mean_height,
                    max_height = max_height,
                    canopy_cover_pct = canopy_cover_pct,
                    area_km2 = area_km2
                )

                message(paste("  - Successfully processed district:", district_name))
            } else {
                message(paste("  - No data found for district:", district_name))
                district_results[[i]] <- data.frame(
                    district_name = district_name,
                    mean_height = NA,
                    max_height = NA,
                    canopy_cover_pct = 0,
                    area_km2 = as.numeric(st_area(district) / 1e6)
                )
            }
        }, error = function(e) {
            message(paste("  - Error processing district:", district_name, "-", e$message))
            district_results[[i]] <- data.frame(
                district_name = district_name,
                mean_height = NA,
                max_height = NA,
                canopy_cover_pct = NA,
                area_km2 = as.numeric(st_area(district) / 1e6)
            )
        })
    }

    # Combine results
    district_metrics <- do.call(rbind, district_results)
    rownames(district_metrics) <- NULL

    # Join metrics back to the spatial data
    districts_with_metrics <- districts %>%
        left_join(district_metrics, by = "district_name")

    return(districts_with_metrics)
}

# Calculate metrics for each district
zagreb_districts_metrics <- calculate_district_metrics(
    zagreb_districts %>%
        mutate(district_name = str_replace(name, "Gradska četvrt ", "")),
    zagreb_chm_agg
)

# 6. URBAN HEAT ISLAND MITIGATION POTENTIAL
# ----------------------------------------
# Calculate cooling potential for each district

# Define cooling effect factors (simplified for demonstration)
cooling_factor_base <- 0.1  # Base cooling effect in °C
cooling_factor_height <- 0.05  # Additional cooling per meter of height

# Calculate cooling potential
zagreb_districts_cooling <- zagreb_districts_metrics %>%
    mutate(
        # Calculate cooling effect based on canopy cover and mean height
        cooling_potential = cooling_factor_base * canopy_cover_pct/100 +
            cooling_factor_height * mean_height * canopy_cover_pct/100,

        # Set cooling potential to 0 for NA values
        cooling_potential = ifelse(is.na(cooling_potential), 0, cooling_potential),

        # Classify cooling potential into categories
        cooling_category = case_when(
            cooling_potential < 0.2 ~ "Low",
            cooling_potential < 0.4 ~ "Moderate",
            cooling_potential < 0.6 ~ "High",
            TRUE ~ "Very High"
        )
    )


# 7. CALCULATE TREE HEIGHT DISTRIBUTION BY DISTRICT
# --------------------------------------------------

# Define height categories
height_categories <- c(0, 5, 15, 25, Inf)
category_names <- c("Small (1-5m)", "Medium (5-15m)", "Large (15-25m)", "Very Large (>25m)")

# Function to calculate area covered by each height category in each district
calculate_height_distribution <- function(districts, chm_raster) {
    message("Calculating tree height distribution by district...")

    # Ensure consistent CRS
    districts_proj <- st_transform(districts, crs(chm_raster))

    # Create separate rasters for each height category
    height_rasters <- list()
    for(i in 1:(length(height_categories)-1)) {
        height_rasters[[i]] <- terra::ifel(
            chm_raster >= height_categories[i] & chm_raster < height_categories[i+1],
            1, NA
        )
    }

    # Create result dataframe
    result_list <- list()

    # Process each district
    for(i in 1:nrow(districts_proj)) {
        district <- districts_proj[i,]
        district_name <- district$district_name
        message(paste("Processing district:", district_name))

        # Calculate total area with trees
        district_vect <- terra::vect(district)
        total_cells <- terra::extract(terra::ifel(!is.na(chm_raster), 1, NA),
                                      district_vect, fun='sum', na.rm=TRUE)[,2]

        # Calculate area for each height category
        category_counts <- numeric(length(category_names))
        category_percentages <- numeric(length(category_names))

        for(j in 1:length(height_rasters)) {
            category_cells <- terra::extract(height_rasters[[j]],
                                             district_vect, fun='sum', na.rm=TRUE)[,2]
            category_counts[j] <- category_cells
            category_percentages[j] <- if(!is.na(total_cells) && total_cells > 0) {
                (category_cells / total_cells) * 100
            } else {
                0
            }
        }

        # Store results
        result_list[[i]] <- data.frame(
            district_name = district_name,
            height_category = category_names,
            tree_count = category_counts,
            percentage = category_percentages
        )
    }

    # Combine all results
    all_results <- do.call(rbind, result_list)

    return(all_results)
}

# Calculate height distribution
zagreb_height_distribution <- calculate_height_distribution(
    zagreb_districts_cooling,
    zagreb_chm_agg
)

# Create a stacked bar chart visualization
height_distribution_plot <- ggplot(zagreb_height_distribution,
                                   aes(x = reorder(district_name, -percentage), y = percentage, fill = height_category)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_viridis_d() +
    theme_minimal() +
    labs(
        title = "Tree Height Distribution by District",
        subtitle = "Percentage of trees in each height category",
        x = "District",
        y = "Percentage (%)",
        fill = "Height Category"
    )

# Save the plot
ggsave(
    "zagreb-height-distribution.png",
    height_distribution_plot,
    width = 10,
    height = 8,
    units = "in",
    bg = "white",
    dpi = 300
)


# Reshape data for easy lookup and join to district data
height_dist_wide <- zagreb_height_distribution %>%
    # First create separate columns for both count and percentage
    pivot_wider(
        id_cols = district_name,
        names_from = height_category,
        values_from = c(percentage, tree_count),
        names_sep = "_"
    )

# Join to district data
zagreb_districts_with_heights <- zagreb_districts_cooling %>%
    left_join(height_dist_wide, by = "district_name")


# 8. PREPARE DATA FOR INTERACTIVE MAP
# ----------------------------------
# Create a more heavily aggregated version specifically for Leaflet visualization
zagreb_chm_leaflet <- terra::aggregate(zagreb_chm_agg, fact=3, fun="mean")

# Convert to raster format for Leaflet compatibility
zagreb_chm_raster_leaflet <- raster::raster(zagreb_chm_leaflet)

# Define color palettes for the map
height_pal <- colorNumeric(
    palette = viridis(64),
    domain = values(zagreb_chm_raster_leaflet),
    na.color = "transparent"
)

cover_pal <- colorNumeric(
    palette = "YlGn",
    domain = zagreb_districts_with_heights$canopy_cover_pct,
    na.color = "transparent"
)

cooling_pal <- colorNumeric(
    palette = "RdYlBu",
    domain = zagreb_districts_with_heights$cooling_potential,
    reverse = FALSE,
    na.color = "transparent"
)

# 9. CREATE THE INTERACTIVE LEAFLET MAP
# -----------------------------------
# Ensure all spatial data is in WGS84 (EPSG:4326) for Leaflet

zagreb_city <- st_transform(zagreb_city, 4326)
zagreb_districts_with_heights <- st_transform(zagreb_districts_with_heights, 4326)


# First, create a named list of popup content for each district
# This avoids having to extract values inside the popup function
district_popups <- lapply(1:nrow(zagreb_districts_with_heights), function(i) {
    district <- zagreb_districts_with_heights[i,]

    # Get all values we need directly from the data frame
    name <- as.character(district$district_name)
    cover <- round(as.numeric(district$canopy_cover_pct), 1)
    mean_h <- round(as.numeric(district$mean_height), 1)
    max_h <- round(as.numeric(district$max_height), 1)
    cooling <- round(as.numeric(district$cooling_potential), 2)
    area <- round(as.numeric(district$area_km2), 1)

    # Get height distribution data - handle potential NA values with correct column names
    small_pct <- ifelse(is.na(district$`percentage_Small (1-5m)`), 0,
                        round(as.numeric(district$`percentage_Small (1-5m)`), 1))
    medium_pct <- ifelse(is.na(district$`percentage_Medium (5-15m)`), 0,
                         round(as.numeric(district$`percentage_Medium (5-15m)`), 1))
    large_pct <- ifelse(is.na(district$`percentage_Large (15-25m)`), 0,
                        round(as.numeric(district$`percentage_Large (15-25m)`), 1))
    very_large_pct <- ifelse(is.na(district$`percentage_Very Large (>25m)`), 0,
                             round(as.numeric(district$`percentage_Very Large (>25m)`), 1))

    # Get tree counts - handle potential NA values with correct column names
    small_count <- ifelse(is.na(district$`tree_count_Small (1-5m)`), 0,
                          format(round(as.numeric(district$`tree_count_Small (1-5m)`)), big.mark=","))
    medium_count <- ifelse(is.na(district$`tree_count_Medium (5-15m)`), 0,
                           format(round(as.numeric(district$`tree_count_Medium (5-15m)`)), big.mark=","))
    large_count <- ifelse(is.na(district$`tree_count_Large (15-25m)`), 0,
                          format(round(as.numeric(district$`tree_count_Large (15-25m)`)), big.mark=","))
    very_large_count <- ifelse(is.na(district$`tree_count_Very Large (>25m)`), 0,
                               format(round(as.numeric(district$`tree_count_Very Large (>25m)`)), big.mark=","))

    # Create an inline SVG bar chart
    chart_width <- 200
    chart_height <- 120
    bar_width <- 35
    spacing <- 10
    y_axis_width <- 25

    # Define a visually appealing color palette (viridis-inspired)
    colors <- c("#440154", "#3b528b", "#21908c", "#5dc963")

    # Find the maximum percentage to scale the chart
    max_pct <- max(small_pct, medium_pct, large_pct, very_large_pct, 25)  # Minimum scale to 25% for visibility

    # Standardize y-axis scale to improve consistency across charts
    # Round up to nearest 20% for cleaner tick values
    max_pct <- ceiling(max_pct / 20) * 20

    # Use fewer ticks with whole numbers for better readability
    num_ticks <- 4  # 0%, 20%, 40%, 60%, etc.
    tick_interval <- max_pct / num_ticks

    # Calculate bar heights based on percentages
    # Reserve 35px for the bottom area (x-axis labels) and 20px for the top (title)
    chart_area_height <- chart_height - 35 - 20
    small_height <- small_pct * chart_area_height / max_pct
    medium_height <- medium_pct * chart_area_height / max_pct
    large_height <- large_pct * chart_area_height / max_pct
    very_large_height <- very_large_pct * chart_area_height / max_pct

    # Create SVG element with a clean, modern style
    svg <- paste0(
        '<svg width="', chart_width, '" height="', chart_height, '" style="background-color: #f8f9fa; border-radius: 5px; padding: 5px;">\n',

        # Add title higher up in the chart to avoid overlapping with bars
        '<text x="', chart_width/2, '" y="15" text-anchor="middle" font-size="11" font-weight="bold" fill="black">Tree Height Distribution</text>\n',

        # Add y-axis line
        '<line x1="', y_axis_width, '" y1="20" x2="', y_axis_width, '" y2="', chart_height - 35,
        '" stroke="black" stroke-width="1" />\n',

        # Add x-axis line
        '<line x1="', y_axis_width, '" y1="', chart_height - 35, '" x2="', chart_width - 5, '" y2="', chart_height - 35,
        '" stroke="black" stroke-width="1" />\n',

        # Add y-axis tick marks and labels - starting with 0%
        '<line x1="', y_axis_width - 3, '" y1="', chart_height - 35, '" x2="', y_axis_width, '" y2="', chart_height - 35,
        '" stroke="black" stroke-width="1" />\n',
        '<text x="', y_axis_width - 5, '" y="', chart_height - 32,
        '" text-anchor="end" font-size="9" fill="black">0%</text>\n'
    )

    # Add additional y-axis ticks - using whole number percentages at even intervals
    for (tick in 1:num_ticks) {
        tick_value <- tick * tick_interval
        tick_position <- chart_height - 35 - (tick_value * chart_area_height / max_pct)

        # Only show the tick if it's within the chart area
        if (tick_position >= 20) {
            svg <- paste0(svg,
                          '<line x1="', y_axis_width - 3, '" y1="', tick_position, '" x2="', y_axis_width, '" y2="', tick_position,
                          '" stroke="black" stroke-width="1" />\n',
                          '<text x="', y_axis_width - 5, '" y="', tick_position + 3,
                          '" text-anchor="end" font-size="9" fill="black">', round(tick_value), '%</text>\n'
            )

            # Add light grid line
            svg <- paste0(svg,
                          '<line x1="', y_axis_width, '" y1="', tick_position, '" x2="', chart_width - 5, '" y2="', tick_position,
                          '" stroke="#ccc" stroke-width="0.5" stroke-dasharray="3,3" />\n'
            )
        }
    }

    # Draw the bars - adjusted positions to account for y-axis and title
    svg <- paste0(svg,
                  # Small trees bar
                  '<rect x="', y_axis_width + 5, '" y="', chart_height - 35 - small_height,
                  '" width="', bar_width, '" height="', small_height,
                  '" fill="', colors[1], '" rx="2" ry="2" />\n',

                  # Medium trees bar
                  '<rect x="', y_axis_width + 5 + bar_width + spacing, '" y="', chart_height - 35 - medium_height,
                  '" width="', bar_width, '" height="', medium_height,
                  '" fill="', colors[2], '" rx="2" ry="2" />\n',

                  # Large trees bar
                  '<rect x="', y_axis_width + 5 + 2 * (bar_width + spacing), '" y="', chart_height - 35 - large_height,
                  '" width="', bar_width, '" height="', large_height,
                  '" fill="', colors[3], '" rx="2" ry="2" />\n',

                  # Very large trees bar
                  '<rect x="', y_axis_width + 5 + 3 * (bar_width + spacing), '" y="', chart_height - 35 - very_large_height,
                  '" width="', bar_width, '" height="', very_large_height,
                  '" fill="', colors[4], '" rx="2" ry="2" />\n'
    )

    # Add category labels below bars - adjusted positions for y-axis
    svg <- paste0(svg,
                  '<text x="', y_axis_width + 5 + bar_width/2, '" y="', chart_height - 20,
                  '" text-anchor="middle" font-size="9" fill="black">1-5m</text>\n',

                  '<text x="', y_axis_width + 5 + bar_width + spacing + bar_width/2, '" y="', chart_height - 20,
                  '" text-anchor="middle" font-size="9" fill="black">5-15m</text>\n',

                  '<text x="', y_axis_width + 5 + 2 * (bar_width + spacing) + bar_width/2, '" y="', chart_height - 20,
                  '" text-anchor="middle" font-size="9" fill="black">15-25m</text>\n',

                  '<text x="', y_axis_width + 5 + 3 * (bar_width + spacing) + bar_width/2, '" y="', chart_height - 20,
                  '" text-anchor="middle" font-size="9" fill="black">>25m</text>\n',

                  '</svg>'
    )

    # Create HTML popup content with the bar chart
    paste0(
        '<div style="min-width:220px; max-width:250px;">',
        "<strong style='font-size:16px;'>", name, "</strong><br>",
        "Canopy Cover: ", cover, "%<br>",
        "Average Height: ", mean_h, " m<br>",
        "Max Height: ", max_h, " m<br>",
        "Cooling Potential: ", cooling, " °C<br>",
        "Area: ", area, " km²<br>",
        "<hr style='margin:8px 0;'>",
        svg,  # Include the SVG chart
        "<div style='font-size:11px; margin-top:5px;'>",
        "Small (1-5m): ", small_count, " trees<br>",
        "Medium (5-15m): ", medium_count, " trees<br>",
        "Large (15-25m): ", large_count, " trees<br>",
        "Very Large (>25m): ", very_large_count, " trees",
        "</div>",
        "</div>"
    )
})

# Add the popups to the spatial object
zagreb_districts_with_heights$popup_content <- unlist(district_popups)


# Create the leaflet map
zagreb_map <- leaflet() %>%
    # Add base map layers
    addProviderTiles("CartoDB.Positron", group = "Carto Light") %>%
    addProviderTiles("Esri.WorldImagery", group = "Satellite") %>%
    addProviderTiles("OpenStreetMap", group = "OpenStreetMap") %>%

    # Add city boundary
    addPolygons(
        data = zagreb_city,
        fillColor = "transparent",
        color = "black",
        weight = 3,
        group = "City Boundary",
        label = "City of Zagreb"
    ) %>%

    # Add district canopy cover layer
    addPolygons(
        data = zagreb_districts_with_heights,
        fillColor = ~cover_pal(canopy_cover_pct),
        weight = 1,
        opacity = 1,
        color = "white",
        dashArray = "3",
        fillOpacity = 0.7,
        label = ~district_name,
        popup = ~popup_content,
        highlight = highlightOptions(
            weight = 3,
            color = "#666",
            dashArray = "",
            fillOpacity = 0.7,
            bringToFront = TRUE
        ),
        group = "Canopy Cover by District"
    ) %>%

    # Add cooling potential layer
    addPolygons(
        data = zagreb_districts_with_heights,
        fillColor = ~cooling_pal(cooling_potential),
        weight = 1,
        opacity = 1,
        color = "white",
        dashArray = "3",
        fillOpacity = 0.7,
        label = ~district_name,
        popup = ~popup_content,
        highlight = highlightOptions(
            weight = 3,
            color = "#666",
            dashArray = "",
            fillOpacity = 0.7,
            bringToFront = TRUE
        ),
        group = "Cooling Potential"
    ) %>%

    # Add tree height layer using the more aggregated raster optimized for Leaflet
    addRasterImage(
        zagreb_chm_raster_leaflet,
        colors = height_pal,
        opacity = 0.7,
        group = "Tree Heights"
    ) %>%

    # Add layer control
    addLayersControl(
        baseGroups = c("Carto Light", "OpenStreetMap", "Satellite"),
        overlayGroups = c(
            "City Boundary",
            "Canopy Cover by District",
            "Cooling Potential",
            "Tree Heights"
        ),
        options = layersControlOptions(collapsed = FALSE)
    ) %>%

    # Add legends
    addLegend(
        position = "bottomright",
        pal = cover_pal,
        values = zagreb_districts_with_heights$canopy_cover_pct,
        title = "Canopy Cover (%)",
        opacity = 0.7,
        group = "Canopy Cover by District"
    ) %>%
    addLegend(
        position = "bottomright",
        pal = cooling_pal,
        values = zagreb_districts_with_heights$cooling_potential,
        title = "Cooling Potential (°C)",
        opacity = 0.7,
        group = "Cooling Potential"
    ) %>%
    addLegend(
        position = "bottomright",
        pal = height_pal,
        values = values(zagreb_chm_raster_leaflet),
        title = "Tree Height (m)",
        opacity = 0.7,
        group = "Tree Heights"
    ) %>%

    # Add scale bar and minimap
    addScaleBar(position = "bottomleft") %>%
    addMiniMap(
        tiles = providers$CartoDB.Positron,
        toggleDisplay = TRUE,
        position = "bottomleft"
    ) %>%

    addResetMapButton() %>%


    # Hide some layers by default
    hideGroup("Cooling Potential") %>%

    # Set the initial view
    setView(
        lng = mean(st_bbox(zagreb_city)[c(1,3)]),
        lat = mean(st_bbox(zagreb_city)[c(2,4)]),
        zoom = 11
    )

# Add title and description
zagreb_map <- zagreb_map %>%
    htmlwidgets::onRender("
    function(el, x) {
      // Create a custom title div
      var titleDiv = document.createElement('div');
      titleDiv.innerHTML = '<div style=\"background-color: rgba(255, 255, 255, 0.8); padding: 6px 8px; border-radius: 4px; box-shadow: 0 1px 5px rgba(0,0,0,0.2);\"><h3 style=\"text-align:center; margin: 0; font-size: 16px;\">Zagreb Tree Height Canopy Derived from Meta/WRI data</h3><p style=\"text-align:center; margin: 5px 0 0 0; font-size: 12px;\">Use the layers control to toggle between different analyses.</p></div>';
      titleDiv.style.position = 'absolute';
      titleDiv.style.top = '10px';
      titleDiv.style.left = '60px';
      titleDiv.style.zIndex = '999';

      // Add it to the map container
      document.querySelector('.leaflet').appendChild(titleDiv);
    }
  ")

# 10. SAVE THE INTERACTIVE MAP
# -----------------------------
# Save the map as an HTML file
saveWidget(
    zagreb_map,
    file = "zagreb-interactive-tree-canopy.html",
    selfcontained = TRUE
)

