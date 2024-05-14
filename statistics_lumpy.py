import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

# Specify the directory containing files
directory = ""

vcf_paths = [os.path.join(directory, file) for file in os.listdir(directory) if file.endswith(".vcf")]

# Function to parse a single VCF file
def parse_vcf(file_path, quality_threshold, mapq_threshold=100):
    mutation_types = ["DUP", "DEL", "INV", "BND"]
    mutation_counts = {mut_type: 0 for mut_type in mutation_types}
    mutation_lengths = {mut_type: [] for mut_type in mutation_types}

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue

            columns = line.strip().split('\t')
            if len(columns) < 8:
                continue

            qual = float(columns[5])
            info = columns[7]
            if qual < quality_threshold or "IMPRECISE" in info:
                continue

            info_dict = {entry.split('=')[0]: entry.split('=')[1] for entry in info.split(';') if '=' in entry}
            sv_type = info_dict.get('SVTYPE')
            svlen = int(info_dict['SVLEN']) if 'SVLEN' in info_dict else 0

            if sv_type in mutation_counts:
                mutation_counts[sv_type] += 1
                mutation_lengths[sv_type].append(svlen)

    return mutation_counts, mutation_lengths

# Aggregate counts and lengths from all files
aggregate_counts = {mut_type: [] for mut_type in ["DUP", "DEL", "INV", "BND"]}
aggregate_lengths = {mut_type: [] for mut_type in ["DUP", "DEL", "INV", "BND"]}

for vcf_path in vcf_paths:
    counts, lengths = parse_vcf(vcf_path, quality_threshold=100)
    for mut_type in counts:
        aggregate_counts[mut_type].append(counts[mut_type])
        aggregate_lengths[mut_type].extend(lengths[mut_type])

# Calculate statistics including average SVLEN
statistics = {}
for mut_type in aggregate_counts:
    if len(aggregate_counts[mut_type]) > 0:
        avg_length = np.mean(aggregate_lengths[mut_type]) if len(aggregate_lengths[mut_type]) > 0 else 0
        statistics[mut_type] = {
            "MEAN": np.mean(aggregate_counts[mut_type]),
            "SE": stats.sem(aggregate_counts[mut_type]),
            "SD": np.std(aggregate_counts[mut_type]),
            "MEDIAN": np.median(aggregate_counts[mut_type]),
            "COUNTS": np.sum(aggregate_counts[mut_type]),
            "AVG_SVLEN": avg_length
        }

# Converting the statistics to a DataFrame
df_stats = pd.DataFrame(statistics).T
df_stats = df_stats[['MEAN', 'SE', 'SD', 'MEDIAN', 'COUNTS']]

# Creating a table for visualization
fig, ax = plt.subplots(figsize=(10, 3))
ax.axis('tight')
ax.axis('off')
ax.table(cellText=df_stats.values, colLabels=df_stats.columns, rowLabels=df_stats.index, cellLoc = 'center', loc='center', fontsize=24)

plt.title("Statistics of Structural Variants Across VCF Files (Smoove)")

# Saving the table as an image
table_image_path = ""
plt.savefig(table_image_path)
