import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats


def parse_vcf(file_path, quality_threshold): # Function to parse VCF file
    mutation_types = ["DUP", "DEL", "INS", "INV", "BND"]
    mutation_counts = {mut_type: 0 for mut_type in mutation_types}

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

            info_dict = {entry.split('=')[0]: entry.split('=')[1] if '=' in entry else None for entry in
                         info.split(';')}
            sv_type = info_dict.get('SVTYPE', None)

            if sv_type in mutation_counts:
                mutation_counts[sv_type] += 1

    return mutation_counts


def process_directory(directory, quality_threshold): # Function to process a directory of VCF files
    vcf_paths = [os.path.join(directory, file) for file in os.listdir(directory) if file.endswith(".vcf")]
    aggregate_counts = {mut_type: [] for mut_type in ["DUP", "DEL", "INS", "INV", "BND"]}

    for vcf_path in vcf_paths:
        counts = parse_vcf(vcf_path, quality_threshold)
        for mut_type, count in counts.items():
            aggregate_counts[mut_type].append(count)

    return aggregate_counts


def calculate_statistics(aggregate_counts): # Function to calculate statistics
    statistics = {}
    for mut_type, counts in aggregate_counts.items():
        if len(counts) > 0:
            statistics[mut_type] = {
                "MEAN": np.mean(counts),
                "SE": stats.sem(counts),
                "SD": np.std(counts),
                "MEDIAN": np.median(counts),
                "COUNTS": np.sum(counts)
            }
    return statistics


def plot_statistics(df_stats, title, save_path): # Function to plot statistics
    fig, ax = plt.subplots(figsize=(10, 3))
    ax.axis('tight')
    ax.axis('off')
    ax.table(cellText=df_stats.values, colLabels=df_stats.columns, rowLabels=df_stats.index, cellLoc='center',
             loc='center', fontsize=24)
    plt.title(title)
    plt.savefig(save_path)


# Replace these with folder paths
parent_directory = ""
children_directory = ""

# Process each directory
parent_counts = process_directory(parent_directory, quality_threshold=100)
children_counts = process_directory(children_directory, quality_threshold=100)

# Calculate statistics
parent_stats = calculate_statistics(parent_counts)
children_stats = calculate_statistics(children_counts)

# Convert to DataFrame and plot
df_parent_stats = pd.DataFrame(parent_stats).T
plot_statistics(df_parent_stats, "Parent Structural Variants", "~/parent.png")

df_children_stats = pd.DataFrame(children_stats).T
plot_statistics(df_children_stats, "Children Structural Variants", "~/hildren.png")

# Perform T-test for each mutation type
t_test_results = {}
for mut_type in parent_counts:
    parent_values = parent_counts[mut_type]
    children_values = children_counts[mut_type]
    if len(parent_values) > 0 and len(children_values) > 0:
        t_stat, p_value = stats.ttest_ind(parent_values, children_values, equal_var=False)
        t_test_results[mut_type] = {"T-Statistic": t_stat, "P-Value": p_value}

df_t_test_results = pd.DataFrame(t_test_results).T
print(df_t_test_results)
