library(ggplot2)
library(dplyr)
library(readr)
library(stringr)

# The path to the directory containing the TSV files
tsv_dir <- ""

tsv_files <- list.files(path = tsv_dir, pattern = "\\.tsv$", full.names = TRUE)

all_data_list <- list()

# Read and process each file
for (i in seq_along(tsv_files)) {
  file_path <- tsv_files[i]
  
  df <- read_tsv(file_path, show_col_types = FALSE)
  
  if ("Location" %in% colnames(df)) {
    # Classify the variants based on the Location column
    df <- df %>%
      mutate(Classification = case_when(
        str_detect(Location, "exon") ~ "Exonic",
        str_detect(Location, "intron") ~ "Intronic",
        #str_detect(Location, "txStart|txEnd") ~ "Downstream/Upstream",
        TRUE ~ "Intergenic"  # Assumes that if none of the above match, it's Intergenic
      ))
    
    df_summary <- df %>%
      count(Classification) %>%
      mutate(Percentage = n / sum(n) * 100,
             File = i)  # Use the loop index i as the File number
    
    all_data_list[[i]] <- df_summary
  } else {
    stop("The column 'Location' was not found in the file: ", file_path)
  }
}

# Combine all the summaries into one data frame
all_data <- bind_rows(all_data_list)

# Calculate overall percentages
overall_summary <- all_data %>%
  group_by(Classification) %>%
  summarise(Total_Count = sum(n)) %>%
  mutate(Overall_Percentage = Total_Count / sum(Total_Count) * 100)

# Plot the data per file
ggplot(all_data, aes(x = File, y = Percentage, color = Classification)) +
  geom_point() +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 100)) +  # Correct labels
  labs(x = "Individuals", y = "Genome Section in %", title = "Variant Classification Across Files") +
  theme_minimal() +
  theme(legend.title = element_blank(), legend.position = "bottom") +
  scale_color_manual(values = c("Exonic" = "green", "Intronic" = "red", "Downstream/Upstream" = "blue", "Intergenic" = "purple"))

print(overall_summary)