import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

# Specify the directory containing files
directory = ""  # Replace with folder path

vcf_paths = [os.path.join(directory, file) for file in os.listdir(directory) if file.endswith(".vcf")]

# Function to parse a single VCF file
def parse_vcf(file_path, quality_threshold):
    mutation_types = ["DUP", "DEL", "INS", "INV", "BND"]
    mutation_counts = {mut_type: 0 for mut_type in mutation_types}
    deletion_lengths = []

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue

            columns = line.strip().split('\t')
            if len(columns) < 8:
                continue

            if columns[6] != "PASS":
                continue

            qual = float(columns[5])
            info = columns[7]
            if qual < quality_threshold or "IMPRECISE" in info:
                continue

            info_dict = {entry.split('=')[0]: entry.split('=')[1] if '=' in entry else None for entry in info.split(';')}
            sv_type = info_dict.get('SVTYPE', None)
            pos = int(columns[1])

            if sv_type == "DEL":
                end = int(info_dict.get('END', pos))
                deletion_length = end - pos
                deletion_lengths.append(deletion_length)

            if sv_type in mutation_counts:
                mutation_counts[sv_type] += 1

    return mutation_counts, deletion_lengths

# Aggregate counts and deletion lengths from all files
aggregate_counts = {mut_type: [] for mut_type in ["DUP", "DEL", "INS", "INV", "BND"]}
aggregate_deletion_lengths = []

for vcf_path in vcf_paths:
    counts, deletion_lengths = parse_vcf(vcf_path, quality_threshold=100)
    for mut_type, count in counts.items():
        aggregate_counts[mut_type].append(count)
    aggregate_deletion_lengths.extend(deletion_lengths)

# Calculate statistics for each mutation type
statistics = {}
average_deletion_length = np.mean(aggregate_deletion_lengths) if aggregate_deletion_lengths else 0

for mut_type, counts in aggregate_counts.items():
    if len(counts) > 0:
        statistics[mut_type] = {
            "MEAN": np.mean(counts),
            "SE": stats.sem(counts),
            "SD": np.std(counts),
            "MEDIAN": np.median(counts),
            "COUNTS": np.sum(counts)
        }

# Add average deletion length to statistics
statistics['DEL']['AVG_DELETION_LENGTH'] = average_deletion_length

# Convert the statistics to a DataFrame
df_stats = pd.DataFrame(statistics).T

# Check if 'DEL' type exists and if 'AVG_DELETION_LENGTH' is calculated
if 'DEL' in df_stats.index and 'AVG_DELETION_LENGTH' in statistics['DEL']:
    # Add 'AVG_DELETION_LENGTH' to 'DEL' row in DataFrame and rename to 'AVG_SVLEN'
    df_stats.loc['DEL', 'AVG_SVLEN'] = statistics['DEL']['AVG_DELETION_LENGTH']

# Reorder and include all columns
columns = ['MEAN', 'SE', 'SD', 'COUNTS'] #, 'AVG_SVLEN'] if 'AVG_SVLEN' in df_stats.columns else ['MEAN', 'SE', 'SD', 'COUNTS']
df_stats = df_stats[columns]

# Create a table for visualization
fig, ax = plt.subplots(figsize=(12, 4))
ax.axis('tight')
ax.axis('off')
table = ax.table(cellText=df_stats.values, colLabels=df_stats.columns, rowLabels=df_stats.index, cellLoc='center', loc='center', fontsize=10)
table.auto_set_font_size(False)
table.set_fontsize(10)
plt.title("Statistics of Structural Variants VCF Files (Delly)")

# Save the table as an image
table_image_path = ""
plt.savefig(table_image_path)
