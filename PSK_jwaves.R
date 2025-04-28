install.packages("signal")
library(readxl)
library(dplyr)
library(signal)
library(tidyverse)
library(gridExtra)

# -------------------------------------
# 1. Load dataset

PSKj <- read_excel("20241119 Prabav/Prabav [i] 1.xlsx", sheet = 2, range = "A1:W241")

# -------------------------------------
# 2. Define key parameters
# -------------------------------------

# Time window for peak detection (65-85ms)
start_row <- (58 * 2) + 41  # rows start at -20, and each millisecond goes over 2 rows. Not sure why +1.
end_row <- (85 * 2) + 41    # 85ms = row 211

# -------------------------------------
# 3. Smoothing function
# -------------------------------------
smooth_ecg <- function(data, method = "savgol", window_size = 7) {
  if (method == "savgol") {
    # Savitzky-Golay filter
    sgolayfilt(data, p = 3, n = window_size)
  } else {
    # Moving average fallback
    stats::filter(data, rep(1/window_size, window_size), sides = 2)
  }
}

# -------------------------------------
# 4. Enhanced peak detection
# -------------------------------------
extract_peak_info <- function(data_vector, prominence_threshold = 7000) {
  # Smooth the data first
  smoothed <- smooth_ecg(data_vector)
  
  # Find peaks
  peaks <- which(diff(sign(diff(smoothed))) == -2) + 1
  if (length(peaks) == 0) return(list(peak_present = 0, peak_height = NA, peak_idx = NA, trough_idx = NA))
  
  # Find troughs
  troughs <- which(diff(sign(diff(smoothed))) == 2) + 1
  if (length(troughs) == 0) return(list(peak_present = 0, peak_height = NA, peak_idx = NA, trough_idx = NA))
  
  # Calculate prominence for each peak
  prominences <- numeric(length(peaks))
  peak_indices <- numeric(length(peaks))
  trough_indices <- numeric(length(peaks))
  
  for (i in 1:length(peaks)) {
    peak_idx <- peaks[i]
    peak_height <- smoothed[peak_idx]
    
    # Find preceding trough
    preceding_troughs <- troughs[troughs < peak_idx]
    preceding_trough_idx <- if (length(preceding_troughs) > 0) {
      max(preceding_troughs)  # Get the closest preceding trough
    } else {
      1  # Use first point if no preceding trough
    }
    
    # Get the height of the preceding trough
    preceding_trough_height <- smoothed[preceding_trough_idx]
    
    # Calculate prominence directly from preceding trough
    prominences[i] <- peak_height - preceding_trough_height
    peak_indices[i] <- peak_idx
    trough_indices[i] <- preceding_trough_idx
  }
  
  # Determine the highest peak prominence
  max_index <- which.max(prominences)
  max_prominence <- prominences[max_index]
  peak_detected <- as.integer(max_prominence > prominence_threshold)
  
  return(list(peak_present = peak_detected, peak_height = max_prominence, 
              peak_idx = peak_indices[max_index], trough_idx = trough_indices[max_index]))
}

# -------------------------------------
# 5. Process all patient columns in twins
# -------------------------------------
# Create empty results dataframe
j_resultsPSK <- tibble(patient_id = character(), peak_R = integer(), peak_height = numeric(), peak_idx = integer(), trough_idx = integer())

# Iterate through patient columns (assuming time is first column)
for (col in colnames(PSKj)[-1]) {
  # Extract the relevant time window
  window_data <- PSKj[start_row:end_row, col] %>% 
    pull() %>%
    na.omit()  # Remove NA values
  
  # Skip invalid traces
  if (length(window_data) < 5) next 
  
  # Detect peak and extract height
  peak_info <- extract_peak_info(window_data)
  
  # Store result
  j_resultsPSK <- j_resultsPSK %>%
    add_row(
      patient_id = col,
      peak_R = peak_info$peak_present,
      peak_height = peak_info$peak_height,
      peak_idx = peak_info$peak_idx,
      trough_idx = peak_info$trough_idx
    )
}

# Print first few results
print(head(j_resultsPSK))
# prevalence?
mean(j_resultsPSK$peak_R)

# Visualisation checker ----------------------------------------------------

plot_multiple_ecg_windowsPSK <- function(patient_ids) {
  plots <- lapply(patient_ids, function(patient_id) {
    patient_col <- as.character(patient_id)
    
    # Extract full ERG data
    patient_data <- PSKj %>% select(time, all_of(patient_col)) %>% na.omit()
    
    # Get peak and trough information from results
    peak_trough_info <- j_resultsPSK %>% dplyr::filter(as.character(patient_id) == patient_col)
    
    if (nrow(peak_trough_info) == 0) {
      return(NULL)
    }
    
    # Convert relative peak/trough index to absolute index in full data
    peak_idx_absolute <- start_row + peak_trough_info$peak_idx - 1
    trough_idx_absolute <- start_row + peak_trough_info$trough_idx - 1
    
    ggplot(patient_data, aes(x = time, y = .data[[patient_col]])) +
      geom_line(color = "steelblue") +
      
      # Overlay detected peak
      geom_point(aes(x = time[peak_idx_absolute], 
                     y = .data[[patient_col]][peak_idx_absolute]), 
                 color = "red", size = 3) +
      
      # Overlay preceding trough
      geom_point(aes(x = time[trough_idx_absolute], 
                     y = .data[[patient_col]][trough_idx_absolute]), 
                 color = "green", size = 3) +
      
      labs(
        title = paste("Full ERG for Patient", patient_id),
        x = "Time (ms)",
        y = "Amplitude"
      ) +
      theme_minimal()
  })
  
  plots <- Filter(Negate(is.null), plots) # Remove NULL entries
  
  if (length(plots) > 0) {
    gridExtra::grid.arrange(grobs = plots, ncol = 2)
  } else {
    message("No valid plots to display.")
  }
}


# Refinement --------------------------------------------------------------


plot_multiple_ecg_windowsPSK(c("R1", "R2", "R3", "R4"))
subset(j_resultsR, patient_id %in% c("90421", "90422"))



write_excel_csv(j_resultsPSK, "j_PSK5.csv")


#this positive result 

# Finding a threshold -----------------------------------------------------



#Save to CSV, when happy with parameters

write.csv(j_resultsR, "peak_detection_results.csv", row.names = FALSE)



# Repeat all of that for the Left -----------------------------------------

#Should/Can I create a machine learning model to do this??
#training and validation data set is the Right eye, the left eye is test set!!

twinsL <- read_excel("Twin study data.xlsx", sheet = 4, )

