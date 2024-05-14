# Structural variants analysis

This repository contains scripts for analyzing structural variants (SVs) using various tools and comparing SVs across different groups. Below are the details for each script, including their functionality, input, and output.

## Scripts Overview

### 1. `variants_on_chromosomes.R`
#### Description
Analyzes and visualizes variants on different chromosomes using a VCF file.

#### Input
- **VCF File**: The script requires a VCF file containing variant data.

#### Output
- **Interactive Plot**: An interactive Plotly scatter plot showing the positions of variants on different chromosomes.

### 2. `de_novo_detection.py`
#### Description
Detects de novo variants in a child's genome by comparing it with the genomes of two parents.

#### Input
- **Child VCF File**: VCF file with the child's variants.
- **Parent 1 VCF File**: VCF file with the first parent's variants.
- **Parent 2 VCF File**: VCF file with the second parent's variants.

#### Output
- **De Novo Variants VCF File**: VCF file containing only the de novo variants identified in the child's genome.

### 3. `genome_location.R`
#### Description
Classifies variants based on their genomic location (e.g., exonic, intronic) and visualize the classification.

#### Input
- **TSV Files Directory**: Directory containing TSV files with variant data.

#### Output
- **Summary Data Frame**: Summary of the classifications with counts and percentages.
- **Plot**: ggplot showing the distribution of variant classifications.

### 4. `gnomad.py`
#### Description
Checks for overlaps between de novo variants and variants in the gnomAD database.

#### Input
- **De Novo VCF Directory**: Directory with de novo VCF files.
- **gnomAD VCF File**: Reference VCF file from the gnomAD database.

#### Output
- **Statistics**: Percentage of de novo variants found in gnomAD.
- **Pie Chart**: Visualization of the percentage of de novo variants found and not found in gnomAD.
- **Found/Not Found VCF Files**: VCF files for variants found and not found in gnomAD.

### 5. `statistics_lumpy.py`
#### Description
Calculates statistics for structural variants detected by Lumpy.

#### Input
- **Directory**: Directory containing VCF files generated by Lumpy.

#### Output
- **Statistics DataFrame**: Statistics for each type of structural variant.
- **Table Image**: Visual summary of the statistics.

### 6. `statistics_delly.py`
#### Description
Calculates statistics for structural variants detected by Delly.

#### Input
- **Directory**: Directory containing VCF files generated by Delly.

#### Output
- **Statistics DataFrame**: Statistics for each type of structural variant.
- **Table Image**: Visual summary of the statistics.

### 7. `compare_SV_in_groups.py`
#### Description
Compares structural variants between two groups and performs statistical tests.

#### Input
- **Parent Directory**: Directory with VCF files for the parent group.
- **Children Directory**: Directory with VCF files for the children group.

#### Output
- **Statistics DataFrames**: Statistics for each type of structural variant for both groups.
- **Plots**: Visualizations of the statistics for each group.
- **T-Test Results**: Results of T-tests comparing SV counts between the groups.

## Usage
1. Ensure all required libraries are installed.
2. Update the file paths and parameters in each script as needed.
3. Run the scripts in the order that fits your analysis workflow.

