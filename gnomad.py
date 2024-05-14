import os
import pysam
import matplotlib.pyplot as plt
from concurrent.futures import ThreadPoolExecutor

def read_variants(vcf_file): "Reads VCF file and extracts variant information"
    variants = []
    with pysam.VariantFile(vcf_file) as vcf:
        for record in vcf:
            variant_data = {
                'chrom': record.chrom,
                'start': record.start,
                'end': record.stop,
                'svtype': record.info.get('SVTYPE', None),
                'qual': record.qual,
                'record': record
            }
            variants.append(variant_data)
    return variants

"Calculates the overlap percentage between two genomic intervals in de novo VCF and reference"
def find_overlaps_in_reference(de_novo_variants, reference_variants, overlap_threshold=70):
    found_in_reference = []
    not_found_in_reference = []
    for de_novo_variant in de_novo_variants:
        found = False
        for reference_variant in reference_variants:
            overlap_start = max(de_novo_variant['start'], reference_variant['start'])
            overlap_end = min(de_novo_variant['end'], reference_variant['end'])
            if overlap_start <= overlap_end:
                overlap_length = overlap_end - overlap_start
                de_novo_length = de_novo_variant['end'] - de_novo_variant['start']
                reference_length = reference_variant['end'] - reference_variant['start']
                max_length = max(de_novo_length, reference_length)
                overlap_percentage = (overlap_length / max_length) * 100

                if overlap_percentage >= overlap_threshold:
                    found = True
                    break
        if found:
            found_in_reference.append(de_novo_variant['record'])
        else:
            not_found_in_reference.append(de_novo_variant['record'])
    return found_in_reference, not_found_in_reference

# Directory with de novo VCF files and reference VCF file
denovo_directory = "/home/quantori/vcf_denovo/denovos/both_delly_lumpy_denovo/"
gnomad_vcf_file = "/home/quantori/Downloads/gnomad_v2.1_sv.sites.vcf"
gnomad_variants = read_variants(gnomad_vcf_file)

def process_file(filename):
    if filename.endswith(".vcf"):
        file_path = os.path.join(denovo_directory, filename)
        denovo_variants = read_variants(file_path)
        return find_overlaps_in_reference(denovo_variants, gnomad_variants)
    return [], []

# Use ThreadPoolExecutor to process each file in parallel
all_found_in_gnomad = []
all_not_found_in_gnomad = []
total_de_novo_count = 0

with ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
    futures = [executor.submit(process_file, filename) for filename in os.listdir(denovo_directory)]
    for future in futures:
        found_in_gnomad, not_found_in_gnomad = future.result()
        all_found_in_gnomad.extend(found_in_gnomad)
        all_not_found_in_gnomad.extend(not_found_in_gnomad)
        total_de_novo_count += len(found_in_gnomad) + len(not_found_in_gnomad)

# Calculate statistics and plot
found_count = len(all_found_in_gnomad)
not_found_count = len(all_not_found_in_gnomad)
found_percentage = found_count / total_de_novo_count * 100 if total_de_novo_count else 0
not_found_percentage = not_found_count / total_de_novo_count * 100 if total_de_novo_count else 0

print(f"Total percentage of de novo variants found in gnomAD: {found_percentage:.2f}%")
print(f"Total percentage of de novo variants not found in gnomAD: {not_found_percentage:.2f}%")

# Plotting the results
labels = ['De novo Variants found in gnomAD', 'De novo Variants not found in gnomAD']
sizes = [found_percentage, not_found_percentage]
colors = ['lightgreen', 'lightcoral']
explode = (0.1, 0)

plt.figure(figsize=(7, 7))
plt.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=140, explode=explode, shadow=True)
plt.title('Percentage of de novo variants found in gnomAD database (all)')
plt.axis('equal')
plt.savefig('~/all_de_novo_gnomad_variants.png', bbox_inches='tight')
plt.close()

# Output to aggregate VCF files
found_vcf_file = ""
not_found_vcf_file = ""

# Output two VCFs - one containing variants found in gnomad, the other - not found in gnomad
template_file_path = os.path.join(denovo_directory, next(f for f in os.listdir(denovo_directory) if f.endswith(".vcf")))
with pysam.VariantFile(template_file_path) as template_vcf:
    with pysam.VariantFile(found_vcf_file, 'w', header=template_vcf.header) as found_vcf:
        for record in all_found_in_gnomad:
            found_vcf.write(record)
    with pysam.VariantFile(not_found_vcf_file, 'w', header=template_vcf.header) as not_found_vcf:
        for record in all_not_found_in_gnomad:
            not_found_vcf.write(record)