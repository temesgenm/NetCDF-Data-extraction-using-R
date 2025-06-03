
#Final code for extracting the Average Tempreature for RCP8.5  data from cordex

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

# Function to process netCDF files
process_netCDF <- function(file_path, lat_city, lon_city) {
  # Read file
  ncin <- nc_open(file_path, write = FALSE, readunlim = TRUE, verbose = TRUE, auto_GMT = TRUE, suppress_dimvals = FALSE)
  
  # Get longitude and latitude matrix
  lat <- ncvar_get(ncin, "lat")
  lon <- ncvar_get(ncin, "lon")
  
  # Get time
  Time <- ncvar_get(ncin, "time")
  Datadates <- as.Date(Time, origin = '1949-12-01')
  
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
  timeseries <- ncvar_get(ncin, varid = 'tas', start = c(ij[1], ij[2], 1), count = c(1, 1, -1))
  
  # Construct the time series
  timeseries_ts <- xts(timeseries, Datadates)
  
  # Convert time series to data frame
  data_frame <- fortify.zoo(timeseries_ts)
  colnames(data_frame) <- c("Date", "tas")
  
  # Return the data frame
  return(data_frame)
}

# Set the folder containing netCDF files
data_folder <- "C:/Users/Temesgen/Desktop/Hydrological_Drought/Temperature_historical/"

# Set the city coordinates
lat_city = 10.741
lon_city = 37.608


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
Tasavg_data <- do.call(rbind, data_frames)

# Converting unit
Date <- Tasavg_data$Date
Tasavg <- (Tasavg_data$tas - 273.15)

# Create a data frame
timeseriesxy <- data.frame(Date, Tasavg)

# Convert Date column to POSIXct format
timeseriesxy$Date <- as.POSIXct(timeseriesxy$Date)

# Add additional columns for Year, Month, and Day
timeseriesxy$Year <- lubridate::year(timeseriesxy$Date)
timeseriesxy$Month_Names <- lubridate::month(timeseriesxy$Date)
timeseriesxy$Day <- lubridate::day(timeseriesxy$Date)

# Display the data frame with additional columns
new_columns_order <- c("Date", "Year", "Month_Names", "Day", "Tasavg")
timeseriesxy <- timeseriesxy[new_columns_order]

# Print the resulting data frame
print(timeseriesxy)

summary(timeseriesxy$Tasavg)

# Create xts object
timeseries_xts <- xts(timeseriesxy$Tasavg, order.by = timeseriesxy$Date)

# Convert to monthly averages
monthly_averages <- apply.monthly(timeseries_xts, mean)
monthly_data_frame <- fortify.zoo(monthly_averages)
colnames(monthly_data_frame) <- c("Date", "Monthly_Tasavg")

# Display summary statistics for monthly averages
summary(monthly_data_frame$Monthly_Tasavg)
# Print the resulting data frames
print(monthly_data_frame)

# Convert to yearly averages
yearly_averages <- apply.yearly(timeseries_xts, mean)
yearly_data_frame <- fortify.zoo(yearly_averages)
colnames(yearly_data_frame) <- c("Date", "Yearly_Tasavg")

# Display summary statistics for yearly averages
summary(yearly_data_frame$Yearly_Tasavg)
print(yearly_data_frame)


###Time series plot

# Plot daily time series with trend line and legend above
ggplot(timeseriesxy, aes(x = Date, y = Tasavg)) +
  geom_line(aes(color = "Daily average Temperature"), size = 1) +
  geom_smooth(aes(y = Tasavg, color = "Trend Line"), method = "lm", se = FALSE, linetype = "solid") +  # Add trend line
  labs(x = "Time", y = "Daily average Temperature (°C)",
       title = "Time Series Line Plot for Daily average Temperature (°C)",
       color = "Legend") +  # Set legend title
  scale_color_manual(values = c("Daily average Temperature" = "blue", "Trend Line" = "red")) +  # Set colors for the legend
  theme(legend.position = "top",  # Move legend to the top
        panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add panel border
        text = element_text(face = "bold", family = "Times New Roman"),  # Set text in bold and Times New Roman font
        plot.title = element_text(hjust = 0.5))  # Center the title

# Monthly Plot with Correct Colors
avg_monthly <- timeseriesxy %>%
  group_by(Month_Names) %>%
  summarise(AvgMonthly = mean(Tasavg, na.rm = TRUE))

# Convert month names to characters
avg_monthly$Month_Names <- as.character(avg_monthly$Month_Names)


# Define month names and corresponding colors
month_info <- data.frame(
  Month_Names = factor(c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"), levels = month.name),
  Color = c("orange", "blue", "green", "purple", "red", "cyan", "pink", "brown", "yellow", "gray", "lightblue", "darkgreen")
)

ggplot(avg_monthly, aes(x = Month_Names, y = AvgMonthly)) +
  geom_col(fill = month_info$Color, color = "black") +
  scale_fill_manual(values = month_info$Color) +
  scale_x_discrete(labels = month_info$Month_Names) +
  labs(x = "Month", y = "Monthly average Temperature (°C)",
       title = "Monthly average Temperature (°C)") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add panel border
        text = element_text(face = "bold", family = "Times New Roman"),  # Set text in bold and Times New Roman font
        plot.title = element_text(hjust = 0.5))  # Center the title

# Plot yearly Average Temperature (°C)
yearly_avg <- timeseriesxy %>%
  group_by(Year) %>%
  summarise(YearlyAvg = mean(Tasavg, na.rm = TRUE))

# Plot yearly with trend line and legend at the top
ggplot(yearly_avg, aes(x = as.Date(paste(Year, 1, 1, sep = "-")), y = YearlyAvg)) +
  geom_line(aes(color = "Yearly average Temperature (°C)"), size = 1) +
  geom_smooth(method = "lm", se = FALSE,linetype = "dashed", aes(group = 1, color = "Trend Line")) +
  labs(x = "Time in year", y = "Yearly average Temperature (°C)",
       title = "Yearly average Temperature (°C)") +
  scale_color_manual(values = c("Yearly average Temperature (°C)" = "blue", "Trend Line" = "red"), name = "Legend") +
  theme_minimal() +
  theme(legend.position = c(0.5, 0.9), legend.direction = "horizontal",
        plot.title = element_text(hjust = 0.5, vjust = 2, size = 14),
        text = element_text(face = "bold", family = "Times New Roman"),  # Set text in bold and Times New Roman font
        panel.border = element_rect(color = "black", fill = NA, size = 1))


###HydroTSM tool 

#generatea time series prec by converting the Tasavg time series into a "zoo" class which is ordered observations which includes irrigular time series
timeseriesxy.ts <- zoo(timeseriesxy$Tasavg, order.by = timeseriesxy$Date)

# Plotting daily, monthly, and annual time series precipitation data using hydroplot
hydroplot(timeseriesxy.ts,var.type = "Temperature",var.unit = "°C", xlab = "Time",ylab = "Temperature in °C", panel.border = element_rect(color = "black", fill = NA, size = 1)) + geom_line(aes(color = "Yearly Temperature (°C) "), size = 1)

# Custom season names with xlab and ylab
hydroplot(timeseriesxy.ts, pfreq = "seasonal", FUN = mean, 
          stype = "default", season.names = c("Summer", "Autumn", "Winter", "Spring"),
          xlab = "year", ylab = "seasonal average Temperature (°C)",
          oma = c(0, 0, 2, 0))  # Adjust the outer margin at the bottom (oma[3])

# Seasonal analysis
# Average seasonal values of precipitation
seasonal_avg <- seasonalfunction(timeseriesxy.ts, FUN = mean, na.rm = TRUE)

# Extracting the seasonal values for each year
DJF <- dm2seasonal(timeseriesxy.ts, season = "DJF", FUN = mean)
MAM <- dm2seasonal(timeseriesxy.ts, season = "MAM", FUN = mean)
JJA <- dm2seasonal(timeseriesxy.ts, season = "JJA", FUN = mean)
SON <- dm2seasonal(timeseriesxy.ts, season = "SON", FUN = mean)

# Plotting the time evolution of the seasonal precipitation values
hydroplot(timeseriesxy.ts, pfreq = "seasonal", FUN = mean, 
          stype = "default", season.names = c("Summer", "Autumn", "Winter", "Spring"),
          xlab = "year", ylab = "Temperature (°C)",
          oma = c(0, 0, 2, 0))  # Adjust the outer margin at the bottom (oma[3])


# Perform Mann-Kendall trend test
mk.test(timeseriesxy$Tasavg)                    #daily
Daily_trend <- mk.test(timeseriesxy$Tasavg)     #monthly

# Monthly trend
mk.test(avg_monthly$AvgMonthly)               #monthly
mk_monthly <- mk.test(avg_monthly$AvgMonthly)
print(mk_monthly)

# Yearly trend plot
mk.test(yearly_avg$YearlyAvg)                 #yearly
mk_Yearly<-mk.test(yearly_avg$YearlyAvg) 
print(mk_Yearly)


###Export timeseries data

# Specify the file path where we want to save the CSV file in the home directory
csv_file_path <- "C:/Users/Temesgen/Desktop/Hydrological_Drought/Temperature_historical/Dingayr_Tas.csv"

# Write the data frame to a CSV file
write.csv(timeseriesxy, file = csv_file_path, row.names = FALSE)

# Print a message indicating that the file has been written
cat(sprintf("CSV file '%s' has been successfully created.\n", csv_file_path))




