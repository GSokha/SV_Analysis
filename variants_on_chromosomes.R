library(plotly)
library(dplyr)

# Read the VCF file
vcf_data <- read.table("", #specify file path"
                       comment.char = "#", header = FALSE, stringsAsFactors = FALSE)

variants <- data.frame(chromosome = vcf_data$V1, position = vcf_data$V2)

chromosome_order <- c(paste0("chr", 1:22), "chrX", "chrY")

variants$chromosome <- factor(variants$chromosome, levels = chromosome_order)

# Count the number of variants per chromosome
variant_counts <- variants %>%
  group_by(chromosome) %>%
  summarise(count = n())

variants <- merge(variants, variant_counts, by = "chromosome")

# Ensure all chromosomes are present in the plot, even if some have no variants
all_chromosomes <- data.frame(chromosome = chromosome_order)
variant_counts <- merge(all_chromosomes, variant_counts, by = "chromosome", all = TRUE)
variant_counts$count[is.na(variant_counts$count)] <- 0  # replace NA with 0

# Create plot with plotly
p <- plot_ly(data = variants, x = ~position, y = ~chromosome, type = 'scatter', mode = 'markers',
             marker = list(opacity = 0.6),
             text = ~paste("Position:", position, "<br>Count:", count),  # hover text
             hoverinfo = "text") %>%
  layout(title = "<i>De novo</i> genomic variants",
         xaxis = list(title = "Position",
                      tickmode = "linear",
                      tick0 = 0,
                      dtick = 25e6,  # Adjusted to 25 million for granularity
                      tickformat = "s",
                      tickangle = 0),
         yaxis = list(title = "Chromosome", 
                      tickmode = "array", 
                      tickvals = 1:length(chromosome_order), 
                      ticktext = chromosome_order),
         margin = list(l = 50, r = 50, t = 40, b = 40)) %>%
  add_annotations(xref = "paper", yref = "y", 
                  x = 1.05, y = ~chromosome, text = ~paste(count), showarrow = FALSE, font = list(color = "blue"))

# Show plot
p
