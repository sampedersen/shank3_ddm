library(writexl)
library(fs)
library(readr)  
library(openxlsx)


##################

workspace <- file.path("C:\\Users\\Feede\\Documents\\Sam_Temp\\Shank3\\Outputs\\May-06\\VALUE_HE_May-06_13-27.Rdata")
output <- file.path("C:\\Users\\Feede\\Documents\\Sam_Temp\\Shank3\\Outputs\\May-06\\")

# Function to load the workspace, convert Summary to Excel and save
convert_to_excel <- function(workspace_path, output_directory) {
  # Load the RData workspace
  load(workspace_path)
  
  # Ensure the 'Summary' object exists in the environment
  if (!exists("Summary")) {
    stop("Summary variable not found in the workspace.")
  }
  
  # Convert Summary to a data frame for Excel export
  summary_df <- as.data.frame(Summary)
  
  # Create the output Excel file path
  base_name <- tools::file_path_sans_ext(basename(workspace_path))
  output_filename <- paste0(base_name, "_Summary.xlsx")
  output_path <- file.path(output_directory, output_filename)
  # Write the data frame to an Excel file
  write.xlsx(summary_df, output_path, rowNames = TRUE, colNames = TRUE)
  
  # Provide feedback
  message("Summary successfully saved as Excel at: ", output_path)
}

# Call the function to convert and save
convert_to_excel(workspace, output)


#############





base_dir <- path("C:\\Users\\Feede\\Documents\\Sam_Temp\\Shank3\\Outputs\\")
run_date <- "May-06"
workspace_name <- "VALUE_WT_May-06_13-27"

load_in <- paste0(base_dir,run_date,"//",workspace_name)
print(load_in)
target_Variable <- as.data.frame(balanced_subset)
target_dir <- path("C:\\Users\\Feede\\Documents\\Sam_Temp\\Shank3\\Outputs\\forBrian\\")
filename <- "subsetted_data.xlsx"
write_xlsx(target_Variable,paste0(path(target_dir,"\\",filename)))

load("C:\\Users\\Feede\\Documents\\Sam_Temp\\Shank3\\Outputs\\May-06\\VALUE_WT_May-06_13-27.RData")
value_HE_data <- Data
load("C:\\Users\\Feede\\Documents\\Sam_Temp\\Shank3\\Outputs\\April-25\\VALUE_WT_April-25_15-53.RData")
value_WT_data <- Data
load("C:\\Users\\Feede\\Documents\\Sam_Temp\\Shank3\\Outputs\\April-24\\HE_April-24_09-26.RData")
offer_HE_data <- Data
load("C:\\Users\\Feede\\Documents\\Sam_Temp\\Shank3\\Outputs\\April-24\\WT_April-24_09-26.RData")
offer_WT_data <- Data
disp("Subset of mice for value-based model (HE): ", unique(value_HE_data$subj_idx))
disp("Subset of mice for value-based model (WT): ", unique(value_WT_data$subj_idx))
disp("Subset of mice for offer-based model (HE): ", unique(offer_HE_data$subj_idx))
disp("Subset of mice for offer-based model (WT): ", unique(offer_WT_data$subj_idx))