install.packages("signal")
library(readxl)
library(dplyr)
library(signal) 
library(tidyverse)
library(gridExtra)

# 1. Load dataset

twinsR <- read_excel("Twin study data.xlsx", sheet = 3)


# 2. Define key parameters

# Time window for peak detection (65-85ms)
start_row <- (58 * 2) + 41  # rows start at -20, and each millisecond goes over 2 rows.
end_row <- (85 * 2) + 41    # 85ms = row 211


# 3. Smoothing function

smooth_ecg <- function(data, method = "savgol", window_size = 7) {
  if (method == "savgol") {
    # Savitzky-Golay filter
    sgolayfilt(data, p = 3, n = window_size)
  } else {
    # Moving average fallback
    stats::filter(data, rep(1/window_size, window_size), sides = 2)
  }
}


# 4. Enhanced peak detection

extract_peak_info <- function(data_vector, prominence_threshold = 25000) {
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


# 5. Process all patient columns in twins


# Create empty results dataframe
j_resultsR <- tibble(patient_id = character(), peak_R = integer(), peak_height = numeric(), peak_idx = integer(), trough_idx = integer())

# Iterate through patient columns (assuming time is first column)
for (col in colnames(twinsR)[-1]) {
  # Extract the relevant time window
  window_data <- twinsR[start_row:end_row, col] %>% 
    pull() %>%
    na.omit()  # Remove NA values
  
  # Skip invalid traces
  if (length(window_data) < 5) next 
  
  # Detect peak and extract height
  peak_info <- extract_peak_info(window_data)
  
  # Store result
  j_resultsR <- j_resultsR %>%
    add_row(
      patient_id = col,
      peak_R = peak_info$peak_present,
      peak_height = peak_info$peak_height,
      peak_idx = peak_info$peak_idx,
      trough_idx = peak_info$trough_idx
    )
}

# Print first few results
print(head(j_resultsR))
# prevalence?
mean(j_resultsR$peak_R)

# Visualisation checker ----------------------------------------------------

plot_multiple_ecg_windows3 <- function(patient_ids) {
  plots <- lapply(patient_ids, function(patient_id) {
    patient_col <- as.character(patient_id)
    
    # Extract full ERG data
    patient_data <- twinsR %>% select(time, all_of(patient_col)) %>% na.omit()
    
    # Get peak and trough information from results
    peak_trough_info <- j_resultsR %>% dplyr::filter(as.character(patient_id) == patient_col)
    
    if (nrow(peak_trough_info) == 0) {
      return(NULL)
    }
    
    # Convert relative peak/trough index to absolute index in full data
    peak_idx_absolute <- start_row + peak_trough_info$peak_idx - 1
    trough_idx_absolute <- start_row + peak_trough_info$trough_idx - 1
    
    ggplot(patient_data, aes(x = time, y = .data[[patient_col]])) +
      geom_line(color = "steelblue") +
      
      # Overlay detected peak
      #geom_point(aes(x = time[peak_idx_absolute], 
                    # y = .data[[patient_col]][peak_idx_absolute]), 
                # color = "red", size = 3) +
      
      # Overlay preceding trough
      #geom_point(aes(x = time[trough_idx_absolute], 
                 #    y = .data[[patient_col]][trough_idx_absolute]), 
                # color = "green", size = 3) +
      
      labs(
        title = paste("Patient", patient_id),
        x = "Time /ms",
        y = "Amplitude /nV"
      ) +
      theme_minimal()+
      theme(
        plot.title = element_text(size = 26, face = "bold"),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 14))
  })
  
  plots <- Filter(Negate(is.null), plots) # Remove NULL entries
  
  if (length(plots) > 0) {
    gridExtra::grid.arrange(grobs = plots, ncol = 2)
  } else {
    message("No valid plots to display.")
  }
}


# Refinement --------------------------------------------------------------


# Plot multiple patients and check for errors
plot_multiple_ecg_windows(c(371, 372))
subset(j_resultsR, patient_id %in% c("371", "372"))
# these both show a strong phnr2, unlike what is seen in harriet's recordings. + a return that a is less sharp than prabav's.
plot_multiple_ecg_windows(c(431, 432))
subset(j_resultsR, patient_id %in% c("431", "432"))

plot_multiple_ecg_windows(c(851, 852)) #852 is a good example of a noist dataset.
subset(j_resultsR, patient_id %in% c("851", "852"))


plot_multiple_ecg_windows(c(32151, 32152)) #both of these should be positive.
subset(j_resultsR, patient_id %in% c("32151", "32152"))

plot_multiple_ecg_windows(c(64161, 64162))
subset(j_resultsR, patient_id %in% c("64161", "64162"))

plot_multiple_ecg_windows(c(70901, 70902))
subset(j_resultsR, patient_id %in% c("70901", "70902"))

plot_multiple_ecg_windows(c(72941, 72942))
subset(j_resultsR, patient_id %in% c("72941", "72942"))

plot_multiple_ecg_windows(c(81061, 81062))
subset(j_resultsR, patient_id %in% c("81061", "81062"))

plot_multiple_ecg_windows(c(90421, 90422))
subset(j_resultsR, patient_id %in% c("90421", "90422"))

plot_multiple_ecg_windows(c(96221, 96222))
subset(j_resultsR, patient_id %in% c("96221", "96222"))

plot_multiple_ecg_windows(c(99601, 99602))
subset(j_resultsR, patient_id %in% c("99601", "99602"))



# Finding a threshold -----------------------------------------------------
#Random number generator was used to sample from 20 of the 200 patients.
#These were sorted into true negatives and true positives. 
#Method for Rejection of false positives: When the threshold is super high (25000) are definitely noise. Can subtract that from our actual prevalance.


#Sample that are not noisy, and close to the border, to see finalise the selection of threshold.
#True positives
plot_multiple_ecg_windows(c(4421, 32151, 32152, 371))
plot_multiple_ecg_windows(c(432, 81061, 90422, 96221))
plot_multiple_ecg_windows(c(96222, 64162))
subset(j_resultsR, patient_id %in% c("4421", "32151", "32152", "371", "432", "81061", "90422", "96221", "96222", "64162"))
#True negatives
plot_multiple_ecg_windows(c(64161, 4422, 70902))
plot_multiple_ecg_windows(c(72941, 72942, 81062, 90421))
plot_multiple_ecg_windows(c(99601, 99602, 852, 70901))
subset(j_resultsR, patient_id %in% c("64161", "4422", "70902", "72941", "72942", "81062", "90421", "99601", "99602", "852", "70901"))
#Difficult/Problematic traces
plot_multiple_ecg_windows3(c(852, 72942))
plot_multiple_ecg_windows3(c(372, 66761))
subset(j_resultsR, patient_id %in% c("70901", "70902", "852"))
#Prevalence at certain thresholds.
#4000nV = 0.6794258.
#5000nV = 0.5645933.
#6000nV = 0.4354067.
#7000nV = 0.3062201.
#8000nV = 0.215311
#9000nV = 0.1674641
#10000nV = 0.1244019 

#Adjusting the SAVGOL FITLER: CASE STUDY. 70901 WAS 11847nV. Does it become NA? It does.

# Checks under 20,000 nV prominance threshold. 5% have a j-wave.
plot_ecg_window(372) #true negative
plot_ecg_window(432) #false negative
plot_ecg_window(3541) #false negative #the messy one
plot_ecg_window(3572) #false positive
plot_ecg_window(4421) #true positive
plot_ecg_window(60731) #true postive

#checking 15,000nV threshold. 33% have a j-wave.

#Save to CSV, when happy with parameters

write.csv(j_resultsR, "peak_detection_results.csv", row.names = FALSE)

# Repeat all of that for the Left -----------------------------------------

#Should/Can I create a machine learning model to do this??
#training and validation data set is the Right eye, the left eye is 

twinsL <- read_excel("Twin study data.xlsx", sheet = 4)

j_resultsL <- tibble(patient_id = character(), peak_L = integer(), peak_height = numeric(), peak_idx = integer(), trough_idx = integer())

for (col in colnames(twinsL)[-1]) {
  # Extract the relevant time window
  window_data <- twinsL[start_row:end_row, col] %>% 
    pull() %>%
    na.omit()  # Remove NA values
  
  # Skip invalid traces
  if (length(window_data) < 5) next 
  
  # Detect peak and extract height
  peak_info <- extract_peak_info(window_data)
  
  # Store result
  j_resultsL <- j_resultsL %>%
    add_row(
      patient_id = col,
      peak_L = peak_info$peak_present,
      peak_height = peak_info$peak_height,
      peak_idx = peak_info$peak_idx,
      trough_idx = peak_info$trough_idx
    )
}

plot_multiple_ecg_windowsL <- function(patient_ids) {
  plots <- lapply(patient_ids, function(patient_id) {
    patient_col <- as.character(patient_id)
    
    # Extract full ERG data
    patient_data <- twinsL %>% select(time, all_of(patient_col)) %>% na.omit()
    
    # Get peak and trough information from results
    peak_trough_info <- j_resultsL %>% dplyr::filter(as.character(patient_id) == patient_col)
    
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

mean(j_resultsL$peak_L)

# Refinement --------------------------------------------------------------
#VALIDATION SET !

sample(1:210, 20) #158  54  60  73 194 151 185 105  57 169 113 156 193  25   2  47 114 138 181  50



#TP
plot_multiple_ecg_windowsL(75181)
plot_multiple_ecg_windowsL(60731)
plot_multiple_ecg_windowsL(22902)
plot_multiple_ecg_windowsL(17662)



#FP
plot_multiple_ecg_windowsL(95612) # data invalid. does that count?
plot_multiple_ecg_windowsL(82701) # na fair. Sadge.
plot_multiple_ecg_windowsL(75552) # unfairly noisy.
plot_multiple_ecg_windowsL(4031)
plot_multiple_ecg_windowsL(20751)


#FN
plot_multiple_ecg_windowsL(94531)
plot_multiple_ecg_windowsL(95611)
plot_multiple_ecg_windowsL(90671)
plot_multiple_ecg_windowsL(50662)



#True Negative
plot_multiple_ecg_windowsL(76312)
plot_multiple_ecg_windowsL(21431)
plot_multiple_ecg_windowsL(30301)
plot_multiple_ecg_windowsL(65251)
plot_multiple_ecg_windowsL(132)
plot_multiple_ecg_windowsL(70902)
plot_multiple_ecg_windowsL( 99612)



# Correlating Left and right eyes -----------------------------------------


both_eyes <- data.frame(peak_sizeL= j_resultsL$peak_height , peak_sizeR = j_resultsR$peak_height, peak_L = j_resultsL$peak_L, peak_R = j_resultsR$peak_R)

ggplot(both_eyes, aes(x = peak_sizeL, y = peak_sizeR)) +
  geom_point()+
  labs(
    title = "Relationship Between Left and Right Eye Peak Amplitide",
    x = "Left Eye Peak Amplitude / nV",
    y = "Right Eye Peak Amplitude /nV"
  )+
  geom_abline(slope = 1, intercept = 0, color = "red", size = 0.4, linetype = "solid")+
  #geom_smooth(method = "lm", se = FALSE, color = "red")+
  theme_minimal()+
  coord_cartesian(xlim = c(0, 20000), ylim = c(0, 20000)) +
  geom_vline(xintercept = 7000, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = 7000, linetype = "dashed", color = "blue") +
  annotate("text", x = 7100, y = 24000, label = "Threshold = 7000", color = "blue", hjust = 0)

model <- lm(peak_sizeL ~ peak_sizeR, data = both_eyes)
summary(model) 

#Quick t-test in excel to see whether there is significant difference.
write.csv(j_resultsR, "twinsRj.csv", row.names = FALSE)
write.csv(j_resultsL, "twinsLj.csv", row.names = FALSE)

