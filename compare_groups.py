import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats


def process_directory(directory, quality_threshold): # function to analyse all vcf files in a specific folder
    vcf_paths = [os.path.join(directory, file) for file in os.listdir(directory) if file.endswith(".vcf")]
    aggregate_counts = {mut_type: [] for mut_type in ["DUP", "DEL", "INS", "INV", "BND"]}

    for vcf_path in vcf_paths:
        counts = parse_vcf(vcf_path, quality_threshold)
        for mut_type, count in counts.items():
            aggregate_counts[mut_type].append(count)

    return aggregate_counts

def parse_vcf(file_path, quality_threshold): #function to parse vcf and get relevant fields
    mutation_types = ["DUP", "DEL", "INS", "INV", "BND"]
    mutation_counts = {mut_type: 0 for mut_type in mutation_types}

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('#'): # skip header lines
                continue  

            columns = line.strip().split('\t')
            if len(columns) < 8: # ensuring there are enough columns
                continue  

            if columns[6] != "PASS": # checking for 'PASS' in the FILTER column
                continue

            # extracting relevant information
            qual = float(columns[5])
            info = columns[7]

            if qual < quality_threshold or "IMPRECISE" in info: # skipping variants with QUAL < quality_threshold or labeled as IMPRECISE

                continue

            # extracting SVTYPE from the INFO column
            info_dict = {entry.split('=')[0]: entry.split('=')[1] if '=' in entry else None for entry in
                         info.split(';')}
            sv_type = info_dict.get('SVTYPE', None)

            # counting the number of different svtypes
            if sv_type in mutation_counts:
                mutation_counts[sv_type] += 1

    return mutation_counts


def calculate_statistics(aggregate_counts): # function to calculate general statistics
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


def plot_statistics(df_stats, title, save_path): # function to plot statistics
    fig, ax = plt.subplots(figsize=(10, 3))
    ax.axis('tight')
    ax.axis('off')
    ax.table(cellText=df_stats.values, colLabels=df_stats.columns, rowLabels=df_stats.index, cellLoc='center',
             loc='center', fontsize=24)
    plt.title(title)
    plt.savefig(save_path)


# path to specific folders
male_directory = ""
female_directory = ""

# calculating svtypes in each directory
male_counts = process_directory(male_directory, quality_threshold=40)
female_counts = process_directory(female_directory, quality_threshold=40)

# calculating statistics
male_stats = calculate_statistics(male_counts)
female_stats = calculate_statistics(female_counts)

# converting results in dataframe and plot
df_male_stats = pd.DataFrame(male_stats).T
plot_statistics(df_male_stats, "Male Structural Variants", "")

df_female_stats = pd.DataFrame(female_stats).T
plot_statistics(df_female_stats, "Female Structural Variants", "")

# performing T test between the groups
t_test_results = {}
for mut_type in male_counts:
    male_values = male_counts[mut_type]
    female_values = female_counts[mut_type]
    if len(male_values) > 0 and len(female_values) > 0:
        t_stat, p_value = stats.ttest_ind(male_values, female_values, equal_var=False)
        t_test_results[mut_type] = {"T-Test": t_stat, "P-Value": p_value}

df_t_test_results = pd.DataFrame(t_test_results).T
print(df_t_test_results)
