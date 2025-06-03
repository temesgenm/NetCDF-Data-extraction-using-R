#Final code for extracting the historical precipitation data from cordex

# Load libraries
library(ncdf4)
library(sp)
library(geosphere)
library(xts)
library(ggplot2)
library(dplyr)
library(lubridate)
library(zoo)
library(dygraphs)
library(RColorBrewer)
library(hydroTSM)
library(hydrostats)
library(trend)
library(zyp)

##process and etractnetCDF files
# Function to process netCDF files
process_netCDF <- function(file_path, lat_city, lon_city) {
  # Read file
  ncin <- nc_open(file_path, write = FALSE, readunlim = TRUE, verbose = TRUE, auto_GMT = TRUE, suppress_dimvals = FALSE)
  
  # Get longitude and latitude matrix
  lat <- ncvar_get(ncin, "lat")
  lon <- ncvar_get(ncin, "lon")
  
  # Get time
  t <- ncvar_get(ncin, "time")
  obsdatadates <- as.Date(t, origin = '1949-12-01')
  
  # Create mesh points coordinates
  coords_matrix <- cbind(as.vector(lon), as.vector(lat))
  mesh_points <- SpatialPoints(coords_matrix, proj4string = CRS("+proj=longlat +datum=WGS84"), bbox = NULL)
  city_point <- SpatialPoints(matrix(c(lon_city, lat_city), nrow = 1, ncol = 2), proj4string = CRS("+proj=longlat +datum=WGS84"), bbox = NULL)
  
  # ID of the closest point to the city in the mesh
  closest_point <- which.min(distGeo(mesh_points, city_point))
  
  # Create point_ID_matrix
  point_ID_matrix <- matrix(1:length(lat), nrow = nrow(lat), ncol = ncol(lat))
  
  # Get index i, j of the closest point
  ij <- which(point_ID_matrix == closest_point, arr.ind = TRUE)
  
  # Get values at location i, j
  timeseries <- ncvar_get(ncin, varid = 'pr', start = c(ij[1], ij[2], 1), count = c(1, 1, -1))
  
  # Construct the time series
  timeseries_ts <- xts(timeseries, obsdatadates)
  
  # Convert time series to data frame
  data_frame <- fortify.zoo(timeseries_ts)
  colnames(data_frame) <- c("Date", "pr")
  
  # Return the data frame
  return(data_frame)
}

# Set the folder containing netCDF files
data_folder <- "C:/Users/Temesgen/Desktop/Hydrological_Drought/Precipitation_historical/"


# Set the city coordinates
lat_city = 11.8666
lon_city = 37.9954

# Initialize an empty list to store data frames
data_frames <- list()

# Loop through each netCDF file in the folder
nc_files <- list.files(data_folder, pattern = ".nc", full.names = TRUE)
for (file_path in nc_files) {
  data <- process_netCDF(file_path, lat_city, lon_city)
  data_frames[[file_path]] <- data
}

# Print the list of data frames
print(data_frames)

# Merge all data frames into one
Prec_hist_data <- do.call(rbind, data_frames)
View(Prec_hist_data)

# Converting unit
Date=Prec_hist_data$Date
Prec=(Prec_hist_data$pr*86400)

# Create a data frame
timeseriesxy <- data.frame(Date, Prec)

# Convert Date column to POSIXct format
timeseriesxy$Date <- as.POSIXct(timeseriesxy$Date)

# Add additional columns for Year, Month, and Day
timeseriesxy$Year <- lubridate::year(timeseriesxy$Date)
timeseriesxy$Month_Names <- lubridate::month(timeseriesxy$Date)
timeseriesxy$Day <- lubridate::day(timeseriesxy$Date)

# Display the data frame with additional columns
new_columns_order <- c("Date", "Year", "Month_Names", "Day", "Prec")
timeseriesxy <- timeseriesxy[new_columns_order]

# Print the resulting data frame
print(timeseriesxy)

summary(timeseriesxy$Prec)

# Create xts object
timeseries_xts <- xts(timeseriesxy$Prec, order.by = timeseriesxy$Date)

# Convert to monthly averages
monthly_averages <- apply.monthly(timeseries_xts, sum)
monthly_data_frame <- fortify.zoo(monthly_averages)
colnames(monthly_data_frame) <- c("Date", "Monthly_Prec")

# Display summary statistics for monthly averages
summary(monthly_data_frame$Monthly_Prec)
# Print the resulting data frames
print(monthly_data_frame)

# Convert to yearly averages
yearly_averages <- apply.yearly(timeseries_xts, sum)
yearly_data_frame <- fortify.zoo(yearly_averages)
colnames(yearly_data_frame) <- c("Date", "Yearly_Prec")

# Display summary statistics for yearly averages
summary(yearly_data_frame$Yearly_Prec)
print(yearly_data_frame)


##time series data plot for different time step
# Plot daily time series with trend line and legend above
ggplot(timeseriesxy, aes(x = Date, y = Prec)) +
  geom_line(aes(color = "Daily Precipitation"), size = 1) +
  geom_smooth(aes(y = Prec, color = "Trend Line"), method = "lm", se = FALSE, linetype = "solid") +  # Add trend line
  labs(x = "Time", y = "Daily Precipitation (mm)",
       title = "Time Series Line Plot for Daily Precipitation (mm)",
       color = "Legend") +  # Set legend title
  scale_color_manual(values = c("Daily Precipitation" = "blue", "Trend Line" = "red")) +  # Set colors for the legend
  theme(legend.position = "top",  # Move legend to the top
        panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add panel border
        text = element_text(face = "bold", family = "Times New Roman"),  # Set text in bold and Times New Roman font
        plot.title = element_text(hjust = 0.5))  # Center the title


# Monthly Plot with Correct Colors
sum_monthly <- timeseriesxy %>%
  group_by(Month_Names) %>%
  summarise(SumMonthly = sum(Prec, na.rm = TRUE))

# Convert month names to characters
sum_monthly$Month_Names <- as.character(sum_monthly$Month_Names)
month_info$Month_Names<- as.character(month_info$Month_Names)


# Define month names and corresponding colors
month_info <- data.frame(
  Month_Names = factor(c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"), levels = month.name),
  Color = c("orange", "blue", "green", "purple", "red", "cyan", "pink", "brown", "yellow", "gray", "lightblue", "darkgreen")
)

ggplot(sum_monthly, aes(x = Month_Names, y = SumMonthly)) +
  geom_col(fill = month_info$Color, color = "black") +
  scale_fill_manual(values = month_info$Color) +
  scale_x_discrete(labels = month_info$Month_Names) +
  labs(x = "Month", y = "Monthly Precipitation (mm)",
       title = "Monthly Precipitation (mm)") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add panel border
        text = element_text(face = "bold", family = "Times New Roman"),  # Set text in bold and Times New Roman font
        plot.title = element_text(hjust = 0.5))  # Center the title

# Plot yearly Precipitation (mm) 
yearly_sum <- timeseriesxy %>%
  group_by(Year) %>%
  summarise(YearlySum = sum(Prec, na.rm = TRUE))

# Plot yearly with trend line and legend at the top
ggplot(yearly_sum, aes(x = as.Date(paste(Year, 1, 1, sep = "-")), y = YearlySum)) +
  geom_line(aes(color = "Yearly Precipitation (mm)"), size = 1) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", aes(group = 1, color = "Trend Line")) +
  labs(x = "Time in year", y = "Yearly Precipitation (mm)",
       title = "Yearly Precipitation (mm)") +
  scale_color_manual(values = c("Yearly Precipitation (mm)" = "blue", "Trend Line" = "red"), name = "Legend") +
  theme_minimal() +
  theme(legend.position = c(0.5, 0.9), legend.direction = "horizontal",
        plot.title = element_text(hjust = 0.5, vjust = 2, size = 14),
        text = element_text(face = "bold", family = "Times New Roman"),  # Set text in bold and Times New Roman font
        panel.border = element_rect(color = "black", fill = NA, size = 1))

###HydroTSM tool 

#generatea time series prec by converting the prec time series into a "zoo" class which is ordered observations which includes irrigular time series
timeseriesxy.ts <- zoo(timeseriesxy$Prec, order.by = timeseriesxy$Date)

# Set up a larger plotting device
options(repr.plot.width = 12, repr.plot.height = 8)  # Adjust width and height as needed

# Plotting daily, monthly, and annual time series precipitation data using hydroplot
hydroplot(timeseriesxy.ts,var.type = "Precipitation",var.unit = "mm", xlab = "Time",ylab = "Precipitation in mm", panel.border = element_rect(color = "black", fill = NA, size = 1)) + geom_line(aes(color = "Yearly Precipitation (mm) "), size = 1)

# Custom season names with xlab and ylab
hydroplot(timeseriesxy.ts, pfreq = "seasonal", FUN = sum, 
          stype = "default", season.names = c("Summer", "Autumn", "Winter", "Spring"),
          xlab = "year", ylab = "Precipitation (mm/season)")

# Perform Mann-Kendall trend test
mk.test(timeseriesxy$Prec)
result <- mk.test(timeseriesxy$Prec)

# Print the test result
print(result)


##export time series data
# Specify the file path where we want to save the CSV file in home directory
csv_file_path <- "C:/Users/Temesgen/Desktop/Hydrological_Drought/Precipitation_historical/Historical_precipitation.csv"

# Write the data frame to a CSV file
write.csv(timeseriesxy, file = csv_file_path, row.names = FALSE)

# Print a message indicating that the file has been written
cat(sprintf("CSV file '%s' has been successfully created.\n", csv_file_path))