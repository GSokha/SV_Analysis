import pysam

def read_variants(vcf_file): "Read VCF file and extract variant information"
    variants = []
    with pysam.VariantFile(vcf_file) as vcf:
        for record in vcf:
            variant_data = {
                'chrom': record.chrom,
                'start': record.start,
                'end': record.stop,
                'svtype': record.info['SVTYPE'] if 'SVTYPE' in record.info else None,
                'qual': record.qual,
                'record': record
            }
            variants.append(variant_data)
    return variants

def overlap_percentage(interval1, interval2):"Calculate the overlap percentage between two genomic intervals"
    start1, end1 = interval1
    start2, end2 = interval2
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    if overlap_end < overlap_start:
        return 0.0
    overlap_length = overlap_end - overlap_start
    interval1_length = max(end1 - start1, 1)
    interval2_length = max(end2 - start2, 1)
    return (overlap_length / min(interval1_length, interval2_length)) * 100

"Identify de novo variants in a child's genome and compare with parent genomes"
def find_denovo_variants_in_trio(child_variants, parent1_variants, parent2_variants, quality_threshold=100, overlap_threshold=70):
    denovo_variant_records = []
    for child_variant in child_variants:
        if child_variant['qual'] < quality_threshold:
            continue
        is_denovo = True
        for parent_variant in parent1_variants + parent2_variants:
            if (child_variant['chrom'] == parent_variant['chrom'] and
                child_variant['svtype'] == parent_variant['svtype']):
                overlap = overlap_percentage((child_variant['start'], child_variant['end']),
                                             (parent_variant['start'], parent_variant['end']))
                if overlap >= overlap_threshold:
                    is_denovo = False
                    break
        if is_denovo:
            denovo_variant_records.append(child_variant['record'])
    return denovo_variant_records

# Replace with file paths
child_vcf_file = ""
parent1_vcf_file = ""
parent2_vcf_file = ""
output_vcf_file = ""

child_variants = read_variants(child_vcf_file)
parent1_variants = read_variants(parent1_vcf_file)
parent2_variants = read_variants(parent2_vcf_file)

denovo_variants = find_denovo_variants_in_trio(child_variants, parent1_variants, parent2_variants)

# Output results in VCF format
with pysam.VariantFile(child_vcf_file) as template_vcf, pysam.VariantFile(output_vcf_file, 'w', header=template_vcf.header) as output_vcf:
    for record in denovo_variants:
        if 'IMPRECISE' in record.info:
            continue  # Skip writing this variant if it's marked as IMPRECISE
        output_vcf.write(record)
